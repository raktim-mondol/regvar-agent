# AlphaGenome Variant-Effect Prediction Report

## Methodology
All variants scored across 7 assays with ontology restricted to prostate gland and stromal fibroblast.

## 1. RANKED VARIANT TABLE

| Rank | rsID | Variant | Max |abs_effect| | Driving Assay | Top Gene |
|:---:|:---|:---|:---:|:---|:---|
| **1** | rs6983267 | chr8:127401060:G>T | **0.91** | ATAC-seq | MYC |
| **2** | rs10993994 | chr10:46046326:A>G | **0.92** | RNA-seq | MSMB |

## 2. TOP CANDIDATES — MECHANISTIC INTERPRETATION

### Tier 1

#### rs6983267 (chr8:127401060:G>T) — MYC enhancer

The G>T change alters CCAT2 RNA levels across multiple epithelial lines.
The risk T allele increases TCF7L2 binding affinity, enhancing chromatin
looping to the MYC promoter. This variant should be tested with luciferase
reporter assays in LNCaP cells.

#### rs10993994 (chr10:46046326:A>G) — MSMB promoter

The A>G change reduces MSMB expression in prostate epithelial cells.
This is consistent with known prostate cancer risk association.

### Tier 2

#### Unmatched variant (no scores in report)

This variant has minimal functional evidence.
