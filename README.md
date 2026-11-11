# WRKY75, nitrogen response, and leaf senescence

This Week 12 project asks whether genes whose response to high nitrogen depends on WRKY75 genotype in Arabidopsis roots also appear in a published WRKY75 RNAi leaf senescence comparison.

Open [WRKY75_Data_Reuse.html](WRKY75_Data_Reuse.html) in a browser. The report contains its figures, tables, methods, and selected CSV downloads in one offline file. On GitHub, download the HTML file to view it locally.

## Findings

- In the balanced experiment 2 design (12 libraries), DESeq2 detected **1,709 genotype × nitrogen interactions** at BH FDR < 0.05, from 19,871 features with an adjusted P-value.
- **194** interaction loci occurred among 2,159 unique mapped loci in the published leaf lists. The complete tested-gene universe for the leaf experiment was unavailable, so this is a descriptive intersection.
- A sensitivity fit including replicate index retained **183** of the 194 shared candidates. The index tracked PC2 within all four groups, but is not documented as a matched block.
- GO Biological Process analysis of the nitrogen interaction set found **163** significant terms against the testable, annotated nitrogen background. Response to nutrient levels ranked first (58 hits, 2.57× enrichment, BH FDR 2.55 × 10⁻⁸).

## Data and analysis

- [GSE184528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184528) supplies raw feature counts. Experiment 2 uses Col-0 and `wrky75` under control or high nitrogen (20 mM added NH₄NO₃ for two hours), with three replicates per group. See the [verified sample manifest](02_Validation/samples_exp2.csv). Experiment 1 and `wrky40` samples were excluded by design.
- [GSE102302](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102302) supplies published RNAi-versus-Col gene lists. The archives are preserved and the worksheets are read according to their actual formats. Transcript suffixes are removed for locus matching and one unresolvable numeric entry is excluded. No new differential-expression test was run on the two deposited leaf libraries.
- The primary fit uses `~ genotype * nitrogen`, with Col-0 and control as reference levels. The interaction is the high-N effect in `wrky75` minus the high-N effect in Col-0. Features with ≥10 counts in ≥3 samples are retained. Results use Wald tests and BH adjustment with DESeq2 independent filtering.
- GO BP uses `org.At.tair.db` and `GO.db`, a one-sided hypergeometric test and BH adjustment for annotated terms of size 10–500. The background is restricted to AGI loci with a nonmissing adjusted P-value and BP annotation in the primary nitrogen interaction analysis.

The studies differ in tissue, age, and perturbation. Overlap does not establish enrichment, direct WRKY75 targets, or a causal effect of high nitrogen on leaf senescence.

## Reproduce the analysis

The analysis ran with R 4.4.1, Bioconductor 3.19, DESeq2 1.44.0, `org.At.tair.db` 3.19.1 and `GO.db` 3.19.1. See [session information](06_Documentation/sessionInfo.txt) and the [software audit](02_Validation/software.csv). Use Python 3 and R with Cairo support for PDF figures. On a new machine:

```sh
python -m pip install -r requirements.txt
Rscript scripts/install_packages.R
python scripts/run_all.py --fetch
```

If `Rscript` is not on `PATH`, invoke the installer by its full path and pass `--rscript "C:/path/to/Rscript.exe"` to `run_all.py`. Omit `--fetch` if the four original GEO files are already in `01_Raw_Data`. The fetch script verifies the raw file hashes in [input_audit.json](02_Validation/input_audit.json). The analysis requires no credentials.

| Path | Contents |
| --- | --- |
| `scripts/` | GEO retrieval, workbook validation, factorial model, sensitivity checks, HTML rendering |
| `02_Validation/` | Input inventory and hashes, sample assignments, software and ID audit |
| `03_Results/` | Five DESeq2 contrasts, QC, overlap, GO results, sensitivity fit and checks |
| `04_Figures/` | PCA, volcano, overlap, expression heatmap, GO, and QC figures |
| `06_Documentation/` | [Project plan](06_Documentation/Week12_Project_Plan.md), R session and package versions |

The raw downloads and large regenerable R objects are kept locally. Run the workflow to regenerate them and the report.
