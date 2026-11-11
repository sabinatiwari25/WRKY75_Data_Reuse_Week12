from pathlib import Path
import gzip, hashlib, json, csv, re, sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'.validation_deps'))
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '02_Validation'
OUT.mkdir(exist_ok=True)
report = {'files': [], 'workbooks': [], 'gene_lists': []}
raw=ROOT/'01_Raw_Data'
original_names=['GSE102302_WRKY75RNAi_vs_Col_down-regulatedgenes.xls.gz','GSE102302_WRKY75RNAi_vs_Col_up-regulatedgenes.xls.gz','GSE184528_1_raw_reads.xlsx','GSE184528_2_raw_reads.xlsx']
for p in (raw/name for name in original_names):
    report['files'].append({'name': p.name, 'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
    if p.suffix == '.gz':
        target = p.with_suffix('')
        content = gzip.decompress(p.read_bytes())
        if target.exists() and target.read_bytes() != content:
            raise RuntimeError(f'Existing extraction differs: {target}')
        if not target.exists():
            target.write_bytes(content)
        report['files'].append({'name':target.name,'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest(),'extracted_from':p.name})
        sheets=[]
        if content[:2] == b'PK':
            with target.open('rb') as handle:
                book=load_workbook(handle,read_only=True,data_only=True)
                sheets=[(s.title,list(s.values)) for s in book]
                book.close()
        else:
            import xlrd
            book=xlrd.open_workbook(target)
            sheets=[(s.name,[s.row_values(i) for i in range(s.nrows)]) for s in book.sheets()]
        info=[]
        for name, rows in sheets:
            ids=[str(r[0]) for r in rows[1:] if r and r[0] is not None]
            if rows and rows[0][0]=='GeneID':
                direction='down' if 'down-regulated' in p.name else 'up'
                with (OUT/f'senescence_{direction}_published.csv').open('w',newline='',encoding='utf-8') as f:
                    writer=csv.writer(f); writer.writerow(['original_id','gene_id','valid_agi','published_log2FC','published_pvalue','published_FDR'])
                    for r in rows[1:]:
                        original=str(r[0]); gene=re.sub(r'\.\d+$','',original.upper())
                        writer.writerow([original,gene,bool(re.fullmatch(r'AT[1-5CM]G\d{5}',gene)),r[6],r[8],r[9]])
            info.append({'name':name,'rows':len(rows),'columns':max(map(len,rows),default=0), 'preview':rows[:4], 'first_column_unique':len(set(ids)), 'first_column_duplicates':len(ids)-len(set(ids)), 'agi_first_column':sum(bool(re.fullmatch(r'AT[1-5CM]G\d{5}(?:\.\d+)?',x,re.I)) for x in ids)})
        report['gene_lists'].append({'name':target.name,'format':'OOXML' if content[:2]==b'PK' else 'BIFF/OLE','sheets':info})
    if p.suffix == '.xlsx':
        wb = load_workbook(p,read_only=True,data_only=True)
        entry = {'file':p.name,'sheets':[]}
        for ws in wb:
            rows = list(ws.values)
            ids = [str(r[0]) for r in rows[1:] if r and r[0] is not None]
            nonagi=[x for x in ids if not re.fullmatch(r'AT[1-5CM]G\d{5}(?:\.\d+)?',x,re.I)]
            if p.name=='GSE184528_2_raw_reads.xlsx':
                selected=[(i,x) for i,x in enumerate(rows[0]) if x and re.fullmatch(r'(Col-0|wrky75)_(cntr|HN)_exp2_rep[1-3]',x)]
                assert len(selected)==12
                with (OUT/'counts_exp2.csv').open('w',newline='') as f:
                    writer=csv.writer(f); writer.writerow(['gene_id']+[x for i,x in selected])
                    writer.writerows([[r[0]]+[r[i] for i,x in selected] for r in rows[1:]])
                accessions={f'{g}_{n}_exp2_rep{rep}':f'GSM{base+rep-1}' for g,n,base in [('Col-0','cntr',5591950),('Col-0','HN',5591953),('wrky75','cntr',5591961),('wrky75','HN',5591964)] for rep in range(1,4)}
                assert set(x for _,x in selected)==set(accessions)
                with (OUT/'samples_exp2.csv').open('w',newline='') as f:
                    writer=csv.writer(f); writer.writerow(['sample','GEO_accession','genotype','nitrogen','replicate','excel_column','source'])
                    for i,x in selected:
                        genotype,nitrogen,experiment,rep=x.split('_')
                        writer.writerow([x,accessions[x],genotype,nitrogen,rep,i+1,'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184528'])
            entry.setdefault('noncanonical_identifier_examples',nonagi[:20])
            numeric = [v for r in rows[1:] for v in r[1:] if isinstance(v,(int,float))]
            assert len(ids)==len(rows)-1==len(set(ids)), f'Invalid feature identifiers in {p.name} / {ws.title}'
            assert len(numeric)==(len(rows)-1)*(ws.max_column-1) and all(v>=0 and v==int(v) for v in numeric), f'Invalid raw counts in {p.name} / {ws.title}'
            entry['sheets'].append({'name':ws.title,'rows':len(rows),'columns':ws.max_column,'preview':rows[:5],'first_column_unique':len(set(ids)), 'first_column_duplicates':len(ids)-len(set(ids)), 'agi_first_column':sum(bool(re.fullmatch(r'AT[1-5CM]G\d{5}(?:\.\d+)?',x,re.I)) for x in ids), 'numeric_cells':len(numeric),'negative_numeric':sum(v<0 for v in numeric),'noninteger_numeric':sum(v!=int(v) for v in numeric)})
        report['workbooks'].append(entry)
        wb.close()
(OUT/'input_audit.json').write_text(json.dumps(report,indent=2,default=str),encoding='utf-8')
print(json.dumps(report,indent=2,default=str))
