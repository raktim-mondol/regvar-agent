# regvar-agent — Feature Roadmap

> Prioritised list of additional features to extend the current system.

---

## High Value / Low Effort

### 1. VCF Input Support
Parse standard `.vcf` / `.vcf.gz` files directly (via `cyvcf2` or `pysam`) instead of requiring the custom TSV format. Eliminates a manual conversion step for most real-world GWAS / fine-mapping outputs.

### 2. Output Persistence
Save the agent's final report to a timestamped `.md` file and a structured `.tsv` of ranked scores after every run. Currently everything is lost when the process exits.

### 3. Expose `score_variants()` as a Tool
The batch method already exists in `AlphaGenomeClient` but is wired to nothing. Exposing it as a tool lets the agent score multiple variants in one call, reducing turn count and total latency.

### 4. Progress Bar / Logging
Use `tqdm` or `rich` to display per-variant scoring progress to stderr. Particularly useful for inputs with many variants where the agent loop can run silently for several minutes.

---

## Medium Effort / High Impact

### 5. Multi-Tissue / Cell-Type Comparison
Accept multiple UBERON/CL ontology terms and score each variant across all of them, surfacing tissue-specific regulatory effects. Currently hardcoded to prostate gland + fibroblast.

### 6. Gene Annotation Overlay
Annotate variants with nearby gene names, distance to TSS, and overlap with known regulatory elements (promoter, enhancer, UTR) using a BED/GTF file — before or after AlphaGenome scoring.

### 7. LD-Aware Variant Grouping
Group input variants by linkage disequilibrium (LD) blocks using `ldstore`, `plink`, or a precomputed LD matrix so the agent reasons about haplotypes rather than isolated SNPs.

### 8. Structured JSON Report Output
Emit a machine-readable JSON report (`variant → ranked effects → hypotheses → experiments`) in addition to the markdown narrative, enabling downstream programmatic use and database ingestion.

---

## Larger Features

### 9. Full Assay Coverage by Default
`DNase-seq`, `ChIP-seq TF`, and `splicing` are defined in `ASSAY_TO_OUTPUT` but excluded from the default assay list. Enable them via a `--assays all` CLI flag.

### 10. Parallel Variant Scoring
Use `asyncio` or `concurrent.futures` with a semaphore to respect the rate limit while scoring multiple variants concurrently, instead of the current strict 1-per-second sequential loop.

### 11. GWAS Integration
Accept a GWAS summary statistics file and auto-extract fine-mapped credible set variants (e.g. from `susieR` outputs), removing the manual TSV curation step entirely.

### 12. Integration Test Suite
A lightweight integration test using a well-known published variant against the real API, gated behind an environment flag (`RUN_INTEGRATION_TESTS=1`), to catch API/schema regressions early.

### 13. Web UI / Dashboard
A `streamlit` or `gradio` front-end that lets non-programmers upload a TSV, watch scoring progress in real time, and download the final report — without touching the CLI.

---

## Quick Wins (Code Cleanup)

| Item | Detail |
|------|--------|
| Remove dead code | `iter_chunks()` in `variants.py` is defined but never called — use it for batching or delete it |
| Export `score_variants()` | Add it to `regvar/__init__.py` if it is intended to be part of the public API |
| CLI flags | Add `--assays`, `--tissue`, `--top-n`, `--output` arguments to `agent.py` instead of hardcoded defaults |
| Dry-run mode | `--dry-run` flag that validates the TSV and prints what would be scored without making any API calls |

---

## Priority Summary

| Priority | Feature |
|----------|---------|
| P1 | VCF Input Support |
| P1 | Output Persistence (markdown + TSV) |
| P1 | Expose `score_variants()` as a Tool |
| P1 | Progress Bar / Logging |
| P2 | Multi-Tissue Comparison |
| P2 | Gene Annotation Overlay |
| P2 | LD-Aware Variant Grouping |
| P2 | Structured JSON Report Output |
| P3 | Full Assay Coverage (`--assays all`) |
| P3 | Parallel Variant Scoring |
| P3 | GWAS Integration |
| P3 | Integration Test Suite |
| P3 | Web UI / Dashboard |
| P4 | Code cleanup (dead code, CLI flags, dry-run) |
