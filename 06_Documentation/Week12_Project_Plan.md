# Week 12 Data Reuse Project Plan

## Project title
Uncovering Hidden Links Between Nitrogen Signaling and Leaf Senescence Through WRKY75: Reanalysis of Public Arabidopsis Transcriptomic Datasets

## Research objective
Reuse two independently published Arabidopsis transcriptomic datasets to investigate whether WRKY75-associated genes involved in nitrogen response overlap with genes associated with WRKY75-mediated leaf senescence.

## Datasets
1. GSE184528 — Nitrogen-response RNA-seq dataset containing wild-type Col-0 and wrky75 mutant samples under control and high-nitrogen conditions.
2. GSE102302 — Published upregulated and downregulated gene lists comparing WRKY75 RNAi with Col-0 in a leaf-senescence study.

All four downloaded supplementary files are located in the project's 01_Raw_Data folder.

## Analysis plan
Phase 1: Inspect the downloaded datasets, verify sample metadata, identify gene identifiers, and check data quality.

Phase 2: Analyze GSE184528 experiment 2 using DESeq2, including genotype, nitrogen treatment, and their interaction. Select the 12 relevant samples comprising three replicates per genotype-treatment combination.

Phase 3: Extract the published WRKY75-associated gene lists from GSE102302. Do not perform new differential-expression testing on the unreplicated samples.

Phase 4: Harmonize gene identifiers and identify genes shared between the two datasets. Calculate gene overlap and perform an enrichment test only if an appropriate common gene universe is available.

Phase 5: Perform GO Biological Process enrichment when statistically appropriate, using Arabidopsis gene annotations.

Phase 6: Generate scientific figures, including PCA, volcano plot, gene-overlap diagram, expression heatmap, and GO plot where supported.

Phase 7: Interpret the findings, highlighting the new cross-study research question, unexpected biological relationships, limitations, and future research opportunities.

## Final deliverables
- Reproducible R analysis scripts
- Real statistical results and figures
- A professional self-contained HTML report named WRKY75_Data_Reuse.html
- A GitHub README.md
- A GitHub repository named WRKY75_Data_Reuse_Week12
- A Github_repo_link.txt containing the actual published repository URL
- A Week12_Box_Submission folder containing the three files above

## Scientific requirements
Use only real data. Verify all sample assignments and statistical analyses. Do not invent findings, references, or figures. Clearly distinguish completed results from proposed work.

First prioritize a valid pilot analysis, then generate the final HTML report and GitHub files. Ask for approval before publishing the GitHub repository.
