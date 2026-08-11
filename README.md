# Cellular composition of the glutamine–IFN-γ axis in melanoma

This repository contains the reproducible analysis supporting:

> **Cellular Composition Accounts for Much of the Bulk Glutamine–IFN-γ Transcriptional Axis in Melanoma: A Cross-Cohort and Single-Cell Analysis**

Author: **Islam Asal**  
ORCID: [0009-0004-3187-7945](https://orcid.org/0009-0004-3187-7945)  
Preprint: [bioRxiv 10.1101/2025.09.28.679008](https://doi.org/10.1101/2025.09.28.679008)

## Main findings

- The candidate glutamine-associated and IFN-γ-associated scores are inversely correlated in bulk TCGA melanoma RNA profiles.
- The correlation is strongly attenuated after adjustment for tumor purity and immune/stromal composition and is essentially absent after MCP-counter adjustment.
- Patient-level single-cell analysis localizes higher glutamine-associated transcription to malignant cells and higher IFN-γ-associated transcription to immune compartments.
- IFN-γ-associated transcription is favorably prognostic in TCGA-SKCM and independently in GSE65904.
- Glutamine-associated transcription and its interaction with IFN-γ do not provide stable, replicated incremental outcome information.
- GSE91061 contains only 49 pretreatment tumors and 10 responders; it does not support response discrimination and can detect only very large effects.

These results support a cellular-composition interpretation of the bulk transcriptional axis. They do not establish causal nutrient competition or treatment-predictive utility.

## Reproduce the analysis

Run the scripts in this order from the project root:

1. `scripts/05_acquire_strengthening_data.R`
2. `scripts/04_revised_analysis.R`
3. `scripts/06_strengthened_analysis.R`

`scripts/build_strengthened_manuscript.py` regenerates the formatted manuscript from the analysis outputs.

The analysis scripts use paths relative to their project directory. Exact R package versions are recorded in the supplied `sessionInfo.txt` files.

## Repository structure

- `scripts/`: acquisition, analysis, and manuscript-generation code.
- `revision_outputs/`: corrected TCGA-SKCM and GSE91061 results.
- `strengthened_revision_outputs/`: GSE65904 validation, deconvolution, single-cell, missingness, and power results.
- `external_data/strengthened/download_manifest.csv`: source URLs, file sizes, and MD5 checksums for externally acquired files.

All manuscript estimates are available in machine-readable CSV files. Large public source datasets are not duplicated in this repository.

## Public datasets

- TCGA-SKCM: NCI Genomic Data Commons and TCGA Pan-Cancer Clinical Data Resource.
- GSE65904: independent melanoma disease-specific survival cohort.
- GSE91061: pretreatment nivolumab melanoma cohort.
- GSE72056: single-cell metastatic melanoma cohort.

## Analysis safeguards

- Continuous standardized scores are the primary estimands; median-defined groups are descriptive only.
- EPIC nonconverged samples are flagged and excluded from EPIC-based inference.
- Single-cell inference uses patient-level compartment summaries, not individual cells as independent replicates.
- The nominal candidate interaction in GSE65904 is treated as exploratory because it does not replicate in TCGA or with curated signatures.
- Immunotherapy-response associations are not described as treatment-predictive because no untreated comparator is available.

## Citation

Please cite the bioRxiv article and the archived Zenodo software release. The Zenodo DOI will be added here after the first GitHub release is archived.

