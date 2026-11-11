"""Render the verified results as one offline HTML file with embedded figures."""
from base64 import b64encode
from csv import DictReader
from html import escape
from pathlib import Path
import json
import shutil

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'03_Results'
FIGURES=ROOT/'04_Figures'

def read_json(name):
    return json.loads((RESULTS/name).read_text(encoding='utf-8'))

def read_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as file:
        return list(DictReader(file))

def result_csv(name):
    return read_csv(RESULTS/name)

def fmt(value):
    return f'{int(value):,}'

def num(value,digits=2):
    return f'{float(value):.{digits}f}' if value not in (None,'','NA') else '—'

def sci(value):
    return f'{float(value):.2e}' if value not in (None,'','NA') else '—'

def table(headers,rows):
    h=''.join(f'<th scope="col">{escape(str(x))}</th>' for x in headers)
    b=''.join('<tr>'+''.join(f'<td>{escape(str(x))}</td>' for x in row)+'</tr>' for row in rows)
    return f'<div class="table-wrap" tabindex="0"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'

def figure(name,title,caption,alt):
    data=b64encode((FIGURES/f'{name}.png').read_bytes()).decode('ascii')
    return f'<figure><img alt="{escape(alt,quote=True)}" src="data:image/png;base64,{data}"><figcaption><b>{escape(title)}</b> {escape(caption)}</figcaption></figure>'

def download(name,label):
    data=b64encode((RESULTS/name).read_bytes()).decode('ascii')
    return f'<a class="download" download="{name}" href="data:text/csv;base64,{data}">{escape(label)}</a>'

s=read_json('analysis_summary.json')
v=read_json('verification.json')
aud=json.loads((ROOT/'02_Validation/input_audit.json').read_text(encoding='utf-8'))
overlaps=result_csv('descriptive_overlaps.csv')
o=overlaps[0]
primary=s['contrasts'][0]
shared=result_csv('shared_interaction_senescence.csv')
pc1,pc2=s['PCA_variance_percent']
samples=read_csv(ROOT/'02_Validation/samples_exp2.csv')
qc={r['sample']:r for r in result_csv('sample_qc.csv')}
annot={r['TAIR']:r for r in result_csv('candidate_annotations.csv')}
candidate={r['gene_id']:r for r in result_csv('candidate_statistics.csv')}
go=result_csv('GO_BP_interaction.csv')
top_go=go[0]
go_summary=s['GO']
sr=s['senescence']
labels={'interaction':'Genotype × nitrogen interaction','nitrogen_Col0':'High N vs control in Col-0',
        'nitrogen_wrky75':'High N vs control in wrky75','genotype_control':'wrky75 vs Col-0 under control',
        'genotype_HN':'wrky75 vs Col-0 under high N'}
top_shared=shared[:8]
examples=[]
for r in top_shared:
    gene=r['gene_id']
    a=annot.get(gene,{})
    c=candidate[gene]
    symbol={'AT3G22840':'ELIP1','AT2G16060':'PGB1','AT3G54140':'PTR1 / NPF8.1'}.get(gene,a.get('SYMBOL',''))
    examples.append((gene,symbol if symbol not in ('','NA',gene) else '—',num(c['nitrogen_Col0_log2FC']),
                     num(c['nitrogen_wrky75_log2FC']),num(r['log2FoldChange']),sci(r['padj']),r['senescence_direction']))

sample_rows=[(r['GEO_accession'],r['sample'],r['genotype'],'High N' if r['nitrogen']=='HN' else 'Control',r['replicate'][-1],fmt(qc[r['sample']]['library_reads']),num(qc[r['sample']]['size_factor'],3)) for r in samples]
contrast_rows=[(labels[r['contrast']],fmt(r['tested_adjusted']),fmt(r['significant']),fmt(r['positive']),fmt(r['negative'])) for r in s['contrasts']]
overlap_rows=[(labels[r['contrast']],fmt(r['significant']),fmt(r['shared_up']),fmt(r['shared_down']),fmt(r['shared_total'])) for r in overlaps]
go_rows=[(r['GO'],r['description'],fmt(r['query_hits']),num(r['fold_enrichment']),sci(r['padj'])) for r in go[:10]]
hash_rows=[(r['name'],fmt(r['bytes']),r['sha256']) for r in aud['files'] if r['name'].endswith(('.xlsx','.gz'))]
conflicts=result_csv('senescence_direction_conflicts.csv')
conflict_ids=', '.join(sorted({r['gene_id'] for r in conflicts}))
library_min,library_max=(f'{x/1e6:.2f}' for x in s['library_size_range'])
pgb=candidate['AT2G16060']
ptr=candidate['AT3G54140']

CSS='''
:root{--ink:#182d30;--muted:#43585d;--paper:#f7f8f4;--line:#d8e0dc;--teal:#1e6569;--red:#9b3d42;--card:#fff;--serif:Georgia,Cambria,serif;--sans:system-ui,-apple-system,"Segoe UI",sans-serif}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.65 var(--sans)}a{color:#155c60;text-underline-offset:.18em}a:focus-visible,[tabindex]:focus-visible{outline:3px solid #b44748;outline-offset:4px}.wrap{width:min(1120px,calc(100% - 40px));margin-inline:auto}.hero{background:#153f43;color:#fff;padding:62px 0 52px;border-bottom:6px solid #80aa9d}.hero .eyebrow{color:#b6d6ce}.eyebrow{font-size:.75rem;letter-spacing:.14em;text-transform:uppercase;font-weight:750;margin:0 0 14px}.hero h1{font:500 clamp(2.2rem,5vw,4.1rem)/1.08 var(--serif);max-width:1000px;margin:0}.hero .deck{max-width:850px;color:#e2eeea;font-size:1.08rem;margin-bottom:0}.meta{display:flex;gap:.6rem 1.5rem;flex-wrap:wrap;font-size:.85rem;margin-top:25px;color:#c9e1da}.nav{background:#fff;border-bottom:1px solid var(--line)}.nav .wrap{display:flex;flex-wrap:wrap;gap:.45rem 1.4rem;padding:12px 0;font-size:.87rem}.nav a{font-weight:650;text-decoration:none}.nav a:hover{text-decoration:underline}main{padding:25px 0 60px}section{padding:28px 0 10px;scroll-margin-top:16px}section+section{border-top:1px solid var(--line);margin-top:24px}h2{font:500 clamp(1.8rem,3.3vw,2.45rem)/1.16 var(--serif);margin:0 0 15px}h3{font-size:1.15rem;line-height:1.35;margin:25px 0 8px}p{margin:.65rem 0 1rem}p,li{max-width:91ch}.lead{font-size:1.12rem;max-width:87ch}.metrics{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:20px 0}.metric{background:#fff;border:1px solid var(--line);border-top:4px solid var(--teal);padding:14px 18px}.metric strong{display:block;font:500 2.3rem/1.16 var(--serif);color:var(--teal)}.metric span{font-size:.9rem;color:var(--muted)}.callout{padding:14px 18px;background:#eef4ef;border-left:4px solid var(--teal);margin:18px 0}.callout p{margin:0}.warning{background:#fbf1e9;border-color:#ad6b4e}.table-wrap{overflow:auto;margin:13px 0 17px;border:1px solid var(--line);border-radius:4px;background:#fff}table{border-collapse:collapse;width:100%;font-size:.87rem;font-variant-numeric:tabular-nums}th,td{text-align:left;padding:9px 12px;border-bottom:1px solid #e1e8e3;vertical-align:top}th{background:#e8f0ec;color:#153f43;font-weight:700}tr:last-child td{border-bottom:0}tr:nth-child(even){background:#f9fbf9}td:last-child{overflow-wrap:anywhere}figure{margin:20px 0;background:var(--card);border:1px solid var(--line);padding:12px;border-radius:4px;break-inside:avoid}figure img{display:block;width:100%;max-width:880px;height:auto;margin:0 auto}figcaption{font-size:.88rem;color:var(--muted);padding:12px 8px 4px;border-top:1px solid var(--line)}figcaption b{color:var(--ink)}.grid{display:grid;grid-template-columns:1fr 1fr;gap:14px}.grid figure{margin:8px 0}details{margin:14px 0;background:#fff;border:1px solid var(--line);padding:11px 14px}summary{cursor:pointer;font-weight:650;color:var(--teal)}code{font-size:.88em;background:#e8efeb;padding:.1em .22em;border-radius:3px;overflow-wrap:anywhere}.download{display:inline-block;padding:7px 11px;border:1px solid #a0bab2;border-radius:4px;background:white;margin:.25em .4em .25em 0;font-size:.85rem;font-weight:650}.refs{font-size:.9rem}footer{color:#dceae6;background:#153f43;padding:22px 0;font-size:.88rem}.muted{color:var(--muted)}.formula{font-family:var(--serif);font-size:1.06rem}.small{font-size:.85rem}
@media(max-width:760px){.grid,.metrics{grid-template-columns:1fr}.hero{padding:37px 0}.wrap{width:calc(100% - 28px)}.metric strong{font-size:1.85rem}th,td{padding:7px 9px}}
@media print{body{background:white;font-size:10.5pt}.hero{color:#182d30;background:white;border-color:#80aa9d;padding:15px 0}.hero .eyebrow,.hero .deck,.hero .meta{color:#182d30}.nav,.download,footer{display:none}main{padding:0}.wrap{width:100%}section{padding-top:12px}figure{break-inside:avoid}.grid{display:block}}
'''

html=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="A reproducible reanalysis of WRKY75-associated responses in Arabidopsis roots and senescing leaves."><title>WRKY75 and the nitrogen response | Week 12 Data Reuse</title><style>{CSS}</style></head>
<body><header class="hero"><div class="wrap"><p class="eyebrow">Week 12 Data Reuse Project · Arabidopsis transcriptomics</p><h1>Shared WRKY75-associated genes in nitrogen response and leaf senescence</h1><p class="deck">Reanalysis of GEO GSE184528 and GSE102302 asks which genes respond differently to high nitrogen in <i>wrky75</i> roots and also appear in a published WRKY75 RNAi leaf comparison.</p><div class="meta"><span>12 root RNA-seq libraries</span><span>2 independent studies</span><span>Factorial DESeq2 reanalysis</span></div></div></header>
<nav class="nav" aria-label="Report contents"><div class="wrap"><a href="#summary">Summary</a><a href="#design">Data and methods</a><a href="#results">Nitrogen response</a><a href="#overlap">Cross-study overlap</a><a href="#go">Biological processes</a><a href="#interpretation">Interpretation</a><a href="#reproducibility">Reproducibility</a></div></nav>
<main class="wrap">
<section id="summary"><p class="eyebrow">01 · Abstract</p><h2>A shared set of candidate genes</h2>
<p class="lead">The nitrogen response depended on WRKY75 genotype for <b>{fmt(primary['significant'])} genes</b> at a 5% false discovery rate. Of those, <b>{fmt(s['shared_interaction'])}</b> also appeared in published WRKY75 RNAi leaf lists. The overlap identifies candidates for follow-up across tissues and conditions. A cross-study enrichment P-value could not be calculated because the complete tested-gene universe was unavailable for the leaf study.</p>
<div class="metrics"><div class="metric"><strong>{fmt(primary['significant'])}</strong><span>Significant genotype × nitrogen interactions</span></div><div class="metric"><strong>{fmt(s['shared_interaction'])}</strong><span>Also listed in the leaf study</span></div><div class="metric"><strong>{fmt(go_summary['significant'])}</strong><span>Enriched GO Biological Process terms in the nitrogen interaction set</span></div></div>
</section>
<section id="design"><p class="eyebrow">02 · Data and methods</p><h2>Two studies, two distinct comparisons</h2>
<p><b>Nitrogen response.</b> GSE184528 experiment 2 profiled roots of seven-day-old <i>Arabidopsis thaliana</i> seedlings after two hours with or without 20 mM NH<sub>4</sub>NO<sub>3</sub> at pH 5.8. The twelve matched GEO accessions make a balanced Col-0 / <i>wrky75</i> × control / high N design with three biological replicates per combination. The five <i>wrky40</i> libraries and experiment 1 (Col-0 / <i>agb1-2</i>) were excluded by design. Control is the study baseline without added high nitrogen. <a href="#ref1">[1]</a></p>
<p><b>Leaf senescence.</b> GSE102302 profiled the third and fourth leaves of 28-day-old Col-0 and WRKY75 RNAi plants, with one deposited library per condition. This project reused the authors’ upregulated and downregulated lists and did not test differential expression from these two libraries. The lists capture associations with RNAi in a senescence study; they do not isolate a senescence time effect or prove direct transcriptional regulation. <a href="#ref2">[2]</a></p>
<details><summary>Sample manifest and count checks</summary>{table(['GEO','Count column','Genotype','Treatment','Rep.','Assigned counts','Size factor'],sample_rows)}<p>The two count workbooks each have one worksheet and {fmt(s['input_genes'])} unique original feature identifiers. All count cells are nonnegative integers. Experiment 2 contains 17 sample columns and experiment 1 contains 12. The experiment 2 IDs include 32,833 AGI-pattern IDs and 1,385 other identifiers, including ENSRNA IDs. All original count IDs were preserved for modeling; exact AGI locus matching was used for cross-study analysis and annotation.</p></details>
<h3>Factorial model and contrasts</h3>
<p>The count matrix was filtered to features with at least 10 reads in at least three selected samples, retaining {fmt(s['retained_genes'])} of {fmt(s['input_genes'])} features. DESeq2 fitted a negative binomial model with median-ratio size factors and the design <code>~ genotype * nitrogen</code>. Col-0 and control were reference levels. The primary Wald test was the interaction; secondary contrasts estimated treatment within genotype and genotype within treatment. P-values were adjusted using Benjamini–Hochberg separately for each contrast, with DESeq2 independent filtering at α = 0.05. Reported log<sub>2</sub> fold changes are unshrunken. No additional fold-change cutoff was imposed. <a href="#ref3">[3]</a></p>
<div class="callout"><p class="formula"><b>Interaction</b> = (<i>wrky75</i> high N − <i>wrky75</i> control) − (Col-0 high N − Col-0 control), on the log<sub>2</sub> scale. A positive value means a more positive high-N response in the mutant; the gene can still decrease in both genotypes.</p></div>
</section>
<section id="results"><p class="eyebrow">03 · Root RNA-seq</p><h2>High nitrogen separates the samples, with a replicate-order pattern</h2>
<p>Assigned counts ranged from {library_min} to {library_max} million per selected library. On the 500 most variable variance-stabilized features, PC1 ({pc1:.1f}%) shifted in the same direction under high N within both genotypes. PC2 ({pc2:.1f}%) tracked replicate number within all {v['PC2_groups']} groups. Because GEO does not establish a matched block from replicate labels, the primary fit kept the factorial design and retained all 12 samples. A sensitivity fit included replicate index as an exploratory covariate.</p>
{figure('PCA','Figure 1. Sample structure.','Variance-stabilized expression, top 500 variable features. Point labels show replicate number.','PCA colored by genotype and shaped by nitrogen condition.')}
<details><summary>Additional quality control figures</summary><div class="grid">{figure('library_sizes','Library sizes.','Raw assigned counts in the 12 selected libraries.','Bar chart of raw library sizes.')}{figure('sample_distance','Sample distances.','Euclidean distance on variance-stabilized counts.','Heatmap of sample-to-sample distances.')}</div>{figure('dispersion','Dispersion estimates.','Gene-wise estimates, fitted trend and final estimates.','DESeq2 mean versus dispersion plot.')}</details>
<h3>Differential expression</h3>
{table(['Contrast','Adjusted P available','FDR < 0.05','Positive LFC','Negative LFC'],contrast_rows)}
<p>The wild type had more significant high-N responses than the mutant ({fmt(s['contrasts'][1]['significant'])} versus {fmt(s['contrasts'][2]['significant'])}). Comparing two significance counts alone does not establish a genotype difference. The interaction directly tests that question and identified {fmt(primary['positive'])} more positive and {fmt(primary['negative'])} less positive high-N responses in <i>wrky75</i> relative to Col-0.</p>
{figure('volcano_interaction','Figure 2. Interaction effect sizes and significance.','Only features with an adjusted P-value are shown. The dashed line marks BH FDR 0.05.','Volcano plot of genotype-by-nitrogen interaction log2 fold change versus adjusted P-value.')}
<div class="callout"><p><b>Replicate-index sensitivity.</b> The exploratory fit found {fmt(v['sensitivity_significant'])} significant interactions, with {fmt(v['significant_in_both'])} significant in both fits and {fmt(v['sign_agreement_in_both'])} of those agreeing in sign. Interaction estimates correlated at r = {v['interaction_log2FC_pearson']:.3f}. Of the {fmt(s['shared_interaction'])} primary shared candidates, {fmt(v['shared_senescence_in_both'])} also met the interaction threshold in the sensitivity fit. The choice of model changes which genes pass the cutoff; the overlap is not confined to a single fit.</p></div>
</section>
<section id="overlap"><p class="eyebrow">04 · Cross-study comparison</p><h2>{fmt(s['shared_interaction'])} interaction genes are also in the leaf lists</h2>
<p>After uppercasing identifiers and stripping terminal transcript suffixes such as <code>.1</code>, {fmt(sr['up_rows'])} upregulated rows mapped to {fmt(sr['up_unique'])} loci, and {fmt(sr['down_rows'])} downregulated rows mapped to {fmt(sr['down_unique'])} loci after excluding the unresolvable numeric ID <code>8</code>. Three loci appeared in both lists after collapse ({escape(conflict_ids)}); each was counted once in the union of {fmt(sr['unique_total'])} genes. Published rows, including transcript suffixes and conflicting directions, remain in the audit tables. An extra down-workbook worksheet contains previously computed GO output; two extra worksheets in the up workbook are empty. Neither was treated as a gene list.</p>
<p>The primary intersection contained {fmt(o['shared_up'])} loci listed as upregulated in RNAi leaves and {fmt(o['shared_down'])} listed as downregulated ({100*s['shared_interaction']/primary['significant']:.1f}% of significant interactions). The published leaf log fold change and the root interaction are different quantities, so their signs cannot be interpreted as replication of the same effect.</p>
{figure('overlap','Figure 3. Descriptive overlap.','The union of mapped leaf loci is compared with the significant nitrogen interaction set. Circle sizes are schematic.','Venn diagram of significant interaction genes and published WRKY75 RNAi leaf loci.')}
<div class="callout warning"><p><b>No cross-study enrichment test.</b> The GEO supplements provide selected leaf hits, but no complete list of genes that were eligible for detection and testing in that study. A genome-wide background, the published hit union, or the nitrogen-only tested set would not be a defensible common universe. The overlap has no enrichment P-value or claim of excess over chance.</p></div>
<h3>Examples among the strongest shared interactions</h3>
{table(['AGI locus','Symbol','Col-0 HN LFC','wrky75 HN LFC','Interaction LFC','Interaction FDR','Leaf RNAi list'],examples)}
<p>For example, {escape(top_shared[0]['gene_id'])} had a high-N log<sub>2</sub> fold change of {num(candidate[top_shared[0]['gene_id']]['nitrogen_Col0_log2FC'])} in Col-0 and {num(candidate[top_shared[0]['gene_id']]['nitrogen_wrky75_log2FC'])} in <i>wrky75</i>, and was listed as {escape(top_shared[0]['senescence_direction'])}. These numbers describe a change in root response and a separate association in leaves. The heatmap shows patterns among selected genes and is not independent confirmation of the interaction.</p>
{figure('expression_heatmap','Figure 4. Expression of the 30 strongest shared genes.','Rows are z-scored variance-stabilized counts and clustered; columns are in study group order. The genes were selected by interaction FDR.','Heatmap of variance-stabilized, standardized expression of 30 shared genes across 12 samples.')}
<details><summary>Descriptive overlaps for secondary contrasts</summary>{table(['Contrast','FDR < 0.05','Leaf up','Leaf down','Union overlap'],overlap_rows)}<p>Direction columns count membership in published leaf lists. The union counts a locus once even if it occurs in both lists.</p></details>
</section>
<section id="go"><p class="eyebrow">05 · Functional enrichment</p><h2>Nutrient response and stress-related processes are over-represented</h2>
<p>GO Biological Process analysis tested the nitrogen interaction set on its own. {fmt(go_summary['query'])} significant AGI loci had BP annotations among {fmt(go_summary['universe'])} annotated genes with an adjusted interaction P-value. We used the Arabidopsis annotation database <code>org.At.tair.db</code> and propagated GO BP annotations (<code>GOALL</code>), retaining terms with 10–500 background genes. A one-sided hypergeometric test and BH adjustment across all {fmt(go_summary['terms_tested'])} eligible terms identified {fmt(go_summary['significant'])} terms at FDR &lt; 0.05. No GO test was reported for the cross-study intersection because its common sampling universe was unavailable. <a href="#ref4">[4]</a></p>
<p>The top term was <b>{escape(top_go['description'])}</b> ({fmt(top_go['query_hits'])} genes, {num(top_go['fold_enrichment'])}× enrichment, FDR {sci(top_go['padj'])}). Several related oxygen-response terms were also prominent, alongside carbohydrate biosynthesis and secondary metabolism. Related GO terms share annotations; oxygen-response enrichment does not show that these roots experienced hypoxia.</p>
{figure('GO_BP','Figure 5. Top 15 enriched GO Biological Process terms.','Fold enrichment compares the term frequency among annotated significant interaction genes with its frequency in the annotated testable background.','Dot plot of GO Biological Process over-representation results.')}
{table(['GO ID','Biological process','Hits','Fold enrichment','BH FDR'],go_rows)}
</section>
<section id="interpretation"><p class="eyebrow">06 · Discussion</p><h2>The overlap narrows the search for a shared program</h2>
<p>The new question here is whether a WRKY75-dependent response to an acute nitrogen perturbation in roots shares genes with a WRKY75-associated expression signature in mature leaves. The {fmt(s['shared_interaction'])} candidates show that these two sets intersect. The GO results suggest that WRKY75 genotype alters expression of genes involved in nutrient response and broader stress-related processes in the root experiment. The overlap and GO results address different questions, and neither establishes a causal link between high nitrogen and leaf senescence.</p>
<p>Two shared loci illustrate distinct response patterns. The class 1 phytoglobin locus <b>AT2G16060</b> (<i>PGB1</i>) decreased under high N in both root genotypes, but less strongly in <i>wrky75</i> (log<sub>2</sub> fold change {num(pgb['nitrogen_Col0_log2FC'])} in Col-0 versus {num(pgb['nitrogen_wrky75_log2FC'])} in the mutant). The peptide transporter locus <b>AT3G54140</b> (<i>PTR1</i>/<i>NPF8.1</i>) increased in both, but less strongly in <i>wrky75</i> ({num(ptr['nitrogen_Col0_log2FC'])} versus {num(ptr['nitrogen_wrky75_log2FC'])}). Both appeared in the leaf RNAi lists. These loci are candidates for asking how nitrogen and oxygen-related responses intersect with WRKY75-associated leaf expression. <a href="#ref4">[4]</a></p>
<p><b>Limitations.</b> The studies differ in tissue, plant age, perturbation (mutant versus RNAi), and sequencing and analysis pipelines. The senescence data have no replicate-based estimate of variability in the deposited contrast. The selected leaf lists omit non-hits and therefore cannot support an overlap enrichment test. In the root study, three replicates per group limit precision, and the PC2 pattern raises the possibility of a technical or sampling gradient not specified in GEO. GO annotation coverage and gene detectability affect over-representation estimates; GO terms are correlated and should be read as groups of related functions.</p>
<p><b>Next steps.</b> Retrieve the complete original leaf differential-expression table and identify loci eligible for testing in both experiments before testing overlap. Prioritize shared genes that retain their interaction in the sensitivity fit. Replicated experiments measuring root and leaf expression under matched nitrogen treatments and developmental stages could distinguish tissue specificity from a shared response. A binding assay and complementation would help separate direct WRKY75 regulation from downstream effects.</p>
</section>
<section id="reproducibility"><p class="eyebrow">07 · Reproducibility</p><h2>Trace every result to its input</h2>
<p>The workflow uses R 4.4.1, DESeq2 1.44.0 and Arabidopsis/GO annotations from Bioconductor 3.19. The project includes the validation script, R analysis and sensitivity scripts, raw file hashes, result tables, and session information. The scripts verify all 12 sample assignments and nonnegative integer counts, check the design rank, reconstruct interaction estimates from both pairs of simple contrasts, independently recalculate BH and GO results, and confirm the gene intersection. The two contrast identity errors were below 10<sup>−8</sup>.</p>
<p>{download('contrast_summary.csv','Download contrast summary')}{download('shared_interaction_senescence.csv','Download shared genes')}{download('GO_BP_interaction.csv','Download GO results')}</p>
<details><summary>Original GEO files and SHA-256 hashes</summary>{table(['GEO supplementary file','Bytes','SHA-256'],hash_rows)}</details>
<h3>References</h3><ol class="refs">
<li id="ref1">Boonyaves K, Wu T-Y, Dong Y, Urano D. Interplay between ARABIDOPSIS Gβ and WRKY transcription factors differentiates environmental stress responses. <i>Plant Physiology</i> 190, 813–827 (2022). <a href="https://doi.org/10.1093/plphys/kiac305">DOI</a> · <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184528">GEO GSE184528</a>.</li>
<li id="ref2">Guo P et al. A tripartite amplification loop involving the transcription factor WRKY75, salicylic acid, and reactive oxygen species accelerates leaf senescence. <i>The Plant Cell</i> 29, 2854–2870 (2017). <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5728132/">Article</a> · <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102302">GEO GSE102302</a>.</li>
<li id="ref3">Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. <i>Genome Biology</i> 15, 550 (2014). <a href="https://doi.org/10.1186/s13059-014-0550-8">DOI</a>.</li>
<li id="ref4"><a href="https://bioconductor.org/packages/3.19/data/annotation/html/org.At.tair.db.html">org.At.tair.db 3.19.1</a> and <a href="https://bioconductor.org/packages/3.19/data/annotation/html/GO.db.html">GO.db 3.19.1</a>, Bioconductor annotation packages.</li>
</ol></section></main><footer><div class="wrap">Week 12 Data Reuse Project · WRKY75, nitrogen response and leaf senescence</div></footer></body></html>'''

target=ROOT/'WRKY75_Data_Reuse.html'
target.write_text(html,encoding='utf-8')
submission=ROOT/'Week12_Box_Submission'
submission.mkdir(exist_ok=True)
for file in (target,ROOT/'README.md',ROOT/'Github_repo_link.txt'):
    if file.exists():
        shutil.copy2(file,submission/file.name)
print(f'Rendered {target.name} ({target.stat().st_size:,} bytes), with offline figures and tables.')
