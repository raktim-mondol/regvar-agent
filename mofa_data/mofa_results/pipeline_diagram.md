# MOFA+ Multi-Omics Integration Pipeline — Prostate Cancer Risk Variants

```mermaid
flowchart TB
    subgraph DATA["📦 Data Sources (Public)"]
        direction TB
        G1000["🧬 1000 Genomes Phase 3<br/>503 EUR individuals<br/>88 PCa lead SNPs"]
        TCGA["🧫 TCGA-PRAD<br/>156 samples<br/>13,427 expressed genes"]
        ENCODE["📊 ENCODE Portal<br/>48 peak BED files<br/>8 chromatin marks"]
        GTEx["🔗 GTEx v8 Prostate<br/>823K cis-eQTL pairs<br/>N=221 donors"]
        GWAS["📋 GWAS Catalog<br/>Schumacher 2018<br/>136 lead SNPs"]
    end

    subgraph PREP["🗂️ Data Preparation"]
        direction TB
        CLEAN["🧹 Clean & Organize<br/>385MB → 48MB<br/>70 files → 19 files"]
        ALIGN["🎯 SNP Alignment<br/>136 → 88 common SNPs<br/>across all 4 datasets"]
        MATRICES["📐 Build Input Matrices<br/>• Genotype dosage (503×88)<br/>• Regulatory burden (503×8)<br/>• eQTL burden (503×1)"]
        
        CLEAN --> ALIGN --> MATRICES
    end

    subgraph MOFA["🔬 MOFA+ Factor Analysis"]
        direction TB
        
        subgraph SAMPLE["Sample-Level MOFA"]
            S_IN["503 individuals × 2 views<br/>• Genotypes (88 SNPs)<br/>• Regulatory burden (8 marks)"]
            S_TRAIN["Train MOFA+<br/>8→6 factors retained<br/>ARD + Spike-and-Slab"]
            S_OUT["Results:<br/>• Genotypes R² = 9.4%<br/>• Regulatory R² = 89.3%<br/>• F1 ρ=0.933 regulatory"]
            S_IN --> S_TRAIN --> S_OUT
        end
        
        subgraph SNP["SNP-Level MOFA"]
            SN_IN["88 SNPs × 4 views<br/>• GWAS stats (3)<br/>• Regulatory binary (8)<br/>• eQTL stats (2)<br/>• Population genetics (3)"]
            SN_TRAIN["Train MOFA+<br/>5→3 factors retained<br/>Bernoulli + Gaussian"]
            SN_OUT["Results:<br/>• eQTL R² = 76.8%<br/>• PopGen R² = 71.5%<br/>• 3 SNP clusters found"]
            SN_IN --> SN_TRAIN --> SN_OUT
        end
        
        S_OUT -.-> SELECT["📊 Select Top 25 SNPs<br/>by factor loading magnitude"]
        SN_OUT -.-> SELECT
    end

    subgraph AGENT["🤖 AlphaGenome + LLM Agent Pipeline"]
        direction TB
        CANDIDATES["📝 Create Candidates TSV<br/>25 SNPs with rsID, OR,<br/>MOFA factor loadings"]
        SCORE["⚡ AlphaGenome Scoring<br/>regvar agent run<br/>4 assays per SNP"]
        LLM["🧠 DeepSeek v4 Pro Agent<br/>• Rank by functional impact<br/>• Predict mechanisms<br/>• Propose validation plan"]
        
        CANDIDATES --> SCORE --> LLM
    end

    subgraph SCORES_MOFA["📊 AlphaGenome MOFA Views"]
        direction TB
        BUILD["regvar mofa build<br/>Scores × Genotypes"]
        VIEWS["4 MOFA Views (503×33)<br/>• DNase-seq (2)<br/>• RNA-seq (25)<br/>• ChIP-seq histone (2)<br/>• ChIP-seq TF (4)"]
        BUILD --> VIEWS
    end

    subgraph OUTPUTS["📁 Outputs"]
        direction TB
        O1["📄 agent_analysis.md<br/>17KB ranked analysis<br/>+ 6 validation experiments"]
        O2["📊 alphagenome_scores.tsv<br/>375 effect predictions"]
        O3["📈 variance_explained.png<br/>R² bar charts per view"]
        O4["🔥 factor_heatmaps.png<br/>Latent factor visualization"]
        O5["🧬 snp_clustered_heatmap.png<br/>88 SNPs by factor profile"]
        O6["📋 snp_factor_loadings.tsv<br/>88 SNPs annotated"]
        O7["💾 *.npz model files<br/>Trained MOFA+ models"]
    end

    DATA --> PREP
    MATRICES --> MOFA
    SELECT --> AGENT
    SCORE --> SCORES_MOFA
    MOFA --> OUTPUTS
    AGENT --> OUTPUTS
    SCORES_MOFA --> OUTPUTS

    style DATA fill:#1a1a2e,stroke:#16213e,color:#e0e0e0
    style PREP fill:#0f3460,stroke:#16213e,color:#e0e0e0
    style MOFA fill:#533483,stroke:#3b2667,color:#e0e0e0
    style AGENT fill:#c0392b,stroke:#96281b,color:#fff
    style SCORES_MOFA fill:#2471a3,stroke:#1a5276,color:#fff
    style OUTPUTS fill:#1e8449,stroke:#145a32,color:#fff
```

## Pipeline Summary

| Phase | What | Input | Output |
|-------|------|-------|--------|
| **Data Prep** | Clean, align, organize 4 public datasets | 70 files / 385MB | 19 files / 48MB |
| **MOFA Sample** | Factor analysis on 503 individuals | Genotypes + Regulatory burden | 6 latent factors |
| **MOFA SNP** | Factor analysis on 88 SNPs | GWAS + Reg + eQTL + PopGen | 3 SNP clusters |
| **AlphaGenome** | Score 25 top SNPs across assays | SNP positions + tissue context | 375 effect predictions |
| **LLM Agent** | DeepSeek ranks + interprets + plans | AlphaGenome scores | Tiered analysis + 6 experiments |
| **MOFA Views** | Build views from real scores | Scores TSV + Genotype TSV | 4 assay views (503×33) |

## Key Biological Findings

1. **rs10993994** (chr10, MSMB) — strongest eQTL signal (888 GTEx pairs), internal validation of pipeline
2. **rs183373024** (chr8q24, OR=2.91) — **prostate-specific ChIP-seq effects** in 4 prostate cancer lines
3. **rs138213197** (chr17, OR=3.85) — HOXB13 pioneer factor disruption, highest-risk variant
4. **MOFA Factor 1** captures a genetic→regulatory burden axis (ρ=0.933) — the primary latent dimension
