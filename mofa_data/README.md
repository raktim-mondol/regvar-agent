# MOFA+ Prostate Cancer Multi-Omics Input Bundle

Prepared 2026-06-02 from public sources. All files are 
download-ready and use GRCh38 coordinates throughout.

## Directory layout

```
mofa_prostate_inputs/
├── genotypes/                 # 1000G EUR genotypes for 88 PCa lead SNPs
├── rna_tcga_prad/             # TCGA-PRAD STAR counts, TPM, log2(TPM+1) for 156 samples
├── gtex_v8_prostate_eqtl/     # GTEx v8 prostate eQTL (significant + egenes)
├── atac_encode/               # ENCODE prostate ATAC + ChIP catalogs, peak BEDs, consensus
├── chip_encode/               # ENCODE prostate ChIP-seq experiment catalog
├── variant_scores/            # Schumacher 2018 lead SNPs, GTEx eQTL overlap, ENCODE peak annot.
├── metadata/                  # Sample manifests
└── examples/                  # Schema example matching user spec
```

## What's in each layer

### 1. Genotypes — `genotypes/`
- **`1000g_eur_genotypes_pca_lead_snps.tsv`** — alt-allele dosage matrix, **503 EUR samples × 88 PCa lead SNPs** (from 1000G Phase 3 high-coverage NYGC release, GRCh38). 
- Format: `sample_id	chr8:127472793:A>C	...` (matches your `example_genotypes.tsv` spec).
- Same data also as `.parquet`.
- **Source:** 1000G Phase 3 GRCh38 phased VCFs, 20190312 release at EBI FTP. PMID 26432245.

### 2. RNA-seq — `rna_tcga_prad/`
- **`tcga_prad_counts_protein_coding.parquet`** — raw STAR counts, 19,962 protein-coding genes × 156 TCGA-PRAD samples (104 Primary Tumor + 52 Solid Tissue Normal).
- **`tcga_prad_tpm_protein_coding.parquet`** — TPM (GDC pipeline).
- **`tcga_prad_log2tpm_expressed.parquet`** — MOFA-ready log2(TPM+1), filtered to genes with median TPM ≥ 1 (13,427 genes × 156 samples).
- **`gene_annotation.tsv`** — ENSEMBL IDs ↔ symbols ↔ biotype.
- **Source:** GDC TCGA-PRAD STAR-Counts workflow, open-access. Sample manifests in `metadata/`.

### 3. ChIP / ATAC / DNase — `atac_encode/`, `chip_encode/`
- **`encode_prostate_experiments_catalog.tsv`** — all 71 ENCODE prostate-related experiments (ATAC, Histone ChIP, TF ChIP, DNase) with biosample / target / file URLs.
- **`peaks_grch38/`** — 48 individual narrowPeak BED files (GRCh38) covering ATAC, H3K27ac, H3K4me1/3, CTCF, H3K27me3, H3K9me3, H3K36me3.
- **`consensus_peaks_grch38/`** — per-mark union/consensus peaks (e.g. H3K27ac: 194,287 peaks merged across 10 biosamples).
- **Source:** ENCODE Portal. Cell lines included: PC-3, 22Rv1, LNCaP, RWPE1/2, C4-2B, VCaP, plus 2 prostate tissue donors.

### 4. GTEx prostate eQTL — `gtex_v8_prostate_eqtl/`
- **`Prostate.v8.signif_variant_gene_pairs.txt.gz`** — 823,690 significant cis-eQTL variant–gene pairs, N=221 GTEx donors.
- **`Prostate.v8.egenes.txt.gz`** — per-gene eGene summary.
- Useful for SNP→target-gene assignment and as a real prostate-tissue reference for AlphaGenome calibration.
- **Source:** GTEx Analysis v8, PMID 32913098.

### 5. Variant scores — `variant_scores/`
- **`schumacher2018_pca_lead_snps_hg38.tsv`** — 136 PCa risk lead SNPs with rsID, GRCh38 coordinates, OR, p-value, risk allele frequency. From Schumacher et al. *Nat Genet* 2018 (GCST006085), retrieved via GWAS Catalog API + Ensembl GRCh38 lift.
- **`schumacher2018_pca_lead_snps_gtex_prostate_eqtl_overlap.tsv`** — for each SNP: number of GTEx prostate eQTLs within ±100 kb, top eGene, top p-value.
- **`schumacher_pca_snps_encode_regulatory_annotation.tsv`** — binary indicator per SNP per ENCODE mark (ATAC / H3K27ac / H3K4me1 / H3K4me3 / CTCF / H3K27me3 / H3K9me3 / H3K36me3 consensus peaks).

## MOFA+ alignment options

The three layers come from *different individuals*. You cannot directly join 1000G EUR
genotype samples with TCGA-PRAD RNA-seq samples — they're separate cohorts. Choose one
of the following analysis strategies:

| Strategy | Samples (N) | Genotypes | RNA-seq | ATAC/ChIP | Notes |
|---|---|---|---|---|---|
| **A. eQTL-focused** | GTEx-Prostate donors (N=221) | implicit in eQTL slopes | GTEx Prostate expr. | — | Requires GTEx genotypes (dbGaP `phs000424.v8`, controlled) for true MOFA+. Open-access eQTL slopes give a single-view summary. |
| **B. TCGA tumor profile** | 156 TCGA-PRAD | requires dbGaP `phs000178` | ✅ here | average per-mark from ENCODE consensus (cell-type proxy, not sample-matched) | Real RNA-seq + a synthetic regulatory annotation view from ENCODE consensus. |
| **C. Population-genetics demo** | 503 1000G EUR | ✅ here | — (no expression) | ✅ as variant-level annotation | Single-tissue regulatory annotation as a feature view — minimal but reproducible. |

For a publishable MOFA+ run combining real RNA-seq with real genotype-level information, 
**you need dbGaP-controlled access** to either GTEx v8 prostate (phs000424) or TCGA-PRAD 
(phs000178) genotypes. Open-access bundles here let you prototype the pipeline and 
sanity-check the SNP→gene/regulatory feature mappings.

## File-format conventions

- All BED, TSV, gz files are tab-separated and gzipped where size warrants.
- All genomic coordinates: **GRCh38 / hg38**, 1-based for variant TSVs (per VCF convention),
  0-based half-open for BED.
- Sample IDs: TCGA barcodes for `rna_tcga_prad/`; 1000G donor IDs (`HG*`/`NA*`) for `genotypes/`.

## Citations / sources

- Schumacher FR et al. *Nat Genet* 2018, **50**:928–936. PMID 29892016. GWAS Catalog GCST006085.
- 1000 Genomes Project Consortium. *Nature* 2015, **526**:68–74. PMID 26432245. EBI FTP `1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/`.
- GTEx Consortium. *Science* 2020, **369**:1318–1330. PMID 32913098. `storage.googleapis.com/adult-gtex/bulk-qtl/v8/`.
- ENCODE Consortium portal: encodeproject.org. Files retrieved via REST API, GRCh38.
- TCGA Research Network. GDC Portal, project TCGA-PRAD, STAR-Counts workflow.




What's inside (218 MB uncompressed)
Layer	File	Dims / notes
Variant panel	variant_scores/schumacher2018_pca_lead_snps_hg38.tsv	136 lead SNPs from Schumacher 2018 (GWAS Catalog GCST006085), lifted to GRCh38 via Ensembl REST
Genotypes	genotypes/1000g_eur_genotypes_pca_lead_snps.tsv	503 samples × 88 SNPs, alt-allele dosage 0/1/2, from 1000G Phase 3 GRCh38 NYGC high-coverage VCFs
RNA-seq	rna_tcga_prad/tcga_prad_log2tpm_expressed.parquet	13,427 genes × 156 TCGA-PRAD samples (104 Tumor + 52 Normal), log2(TPM+1), filtered median TPM≥1. Raw counts + TPM also included.
GTEx prostate eQTL	gtex_v8_prostate_eqtl/Prostate.v8.signif_variant_gene_pairs.txt.gz	823,690 significant pairs, N=221 donors
ATAC / ChIP peaks	atac_encode/peaks_grch38/ + consensus_peaks_grch38/	48 ENCODE BED files (ATAC, H3K27ac, H3K4me1/3, CTCF, H3K27me3, H3K9me3, H3K36me3) + per-mark consensus merges (e.g. 194,287 consensus H3K27ac peaks)
Catalog	*encode_prostate_experiments_catalog.tsv	71 ENCODE prostate experiments with biosample, target, file URLs — for pulling more
SNP-level annotations	variant_scores/*_eqtl_overlap.tsv + *_regulatory_annotation.tsv	For each lead SNP: number of GTEx prostate eQTLs within ±100 kb, top eGene, top p-value, plus binary peak-overlap flags per ENCODE mark
What's already validated
88/129 autosomal Schumacher lead SNPs found in 1000G Phase 3 (re-mapping coordinates from rsIDs introduced some misses where the rsID currently points to multi-allelic positions outside the simple biallelic VCF)
Genotype matrix has zero missingness on the 88 retained SNPs; alt-AF spans 0.001–0.90 (median 0.31), consistent with mostly common variants
96/136 SNPs (71%) overlap ≥1 GTEx prostate eQTL within 100 kb — strong prior that this panel is regulatory
40% of lead SNPs sit directly in an H3K27ac consensus peak; 15% in ATAC
Honest limitations (read the README)
The three layers are NOT sample-matched. 1000G EUR donors ≠ TCGA-PRAD patients ≠ GTEx donors ≠ ENCODE cell lines. For a true paired multi-omic MOFA+ run you need either:
dbGaP phs000178 (TCGA-PRAD WGS/SNP-array genotypes), or
dbGaP phs000424 (GTEx v8 individual genotypes). Both are controlled-access. I cannot pull them without authorization.
No Hi-C / pcHi-C view. ENCODE has limited prostate Hi-C; would need 4DN Data Portal (open) for HiCS-AGS LNCaP/PC-3 — easy to add as a follow-up.
AlphaGenome variant scores not generated. That requires running regvar agent run --output-tsv on the SNP panel — out of scope of the data-prep step. The 88-SNP TSV (schumacher2018_pca_lead_snps_hg38.tsv) is the input format for that tool.
chrX SNPs excluded from the genotype matrix (different 1000G VCF naming/ploidy convention); 7 X-linked SNPs in the panel TSV remain unmapped.
Recommended next step
Option B in the README — TCGA-PRAD RNA-seq (real) + ENCODE consensus regulatory features (cell-type proxy) — gives you a 2-view MOFA+ at N=156 you can run today, and the variant panel is the entry point for AlphaGenome scoring as a third view.

Discretionary decisions
Variant panel source: Schumacher 2018 (GCST006085) lead SNPs only, retrieved via GWAS Catalog REST API rather than the larger Conti 2021 multi-ancestry credible sets. Schumacher gives a cleaner, well-validated 140-SNP starter panel; Conti's fine-mapped credible sets are denser but require parsing the GCST90274001+ summary statistics directory.
1000G ancestry subset: EUR only (N=503), matched to TCGA-PRAD's predominantly European cohort. Switch to multi-ancestry (N=2,503) by re-running the genotype extraction without the EUR filter.
TCGA-PRAD sample selection: all 52 Solid Tissue Normal + tumor pairs from same patients (n=52) + 50 additional unpaired tumors → N=156. Aimed for a balanced cohort with paired tumor/normal capability; full project has 554 samples.
RNA expression filter: median TPM ≥ 1 across the 156-sample cohort, log2(TPM+1) transform. Removed 33% of protein-coding genes as lowly expressed.
Protein-coding restriction: for the MOFA-ready matrix only — full all-biotype matrix is also provided.
Peak file preference: GRCh38 narrowPeak, output_type ∈ {pseudoreplicated, replicated, IDR thresholded}, preferring IDR/replicated when available.
Consensus peak strategy: simple union/merge across all biosamples per mark. An alternative is intersection-based "high-confidence" peaks (n≥2 biosamples), which would shrink each set ~3-5×.
SNP–regulatory annotation window: direct overlap (variant position inside peak) for the binary annotation; ±100 kb window for the GTEx eQTL co-localization summary.
Excluded chrX from genotype extraction; included in the variant panel TSV.
Did not generate AlphaGenome scores — that's a downstream regvar agent invocation, not a data-prep step.