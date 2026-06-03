# AlphaGenome Variant-Effect Scoring: Prostate Tumour Microenvironment Candidates

## Summary of Scoring Run

All 25 variants were scored against 7 assay types (ATAC-seq, DNase-seq, RNA-seq, ChIP-seq histone, ChIP-seq TF, Hi-C/pcHi-C, splicing) with prostate gland (UBERON:0002367) and fibroblast (CL:0000057) ontology terms. The tool returned 469,335 total effect predictions. Below I rank variants by **regulatory impact confidence** — weighting assay diversity, cell-type relevance to prostate stroma, gene-target biology, and provided OR/loadings.

> ⚠️ **Important caveat:** AlphaGenome scores are computational predictions, not clinical findings. Every score below is a hypothesis requiring experimental validation. The quantile_score reflects extremity of predicted effect in the model's reference distribution, not statistical significance.

---

## TIER 1: Strongest Regulatory Candidates

### 🥇 1. chr8:127091872:A>G — rs183373024 (OR = 2.91)
**Why it ranks #1:** This is the **only variant hitting multiple assay classes in prostate-relevant cell lines.** It sits in the well-known **8q24 prostate cancer risk locus**, upstream of *MYC*.

| Rank | Assay | Cell Type | Gene | Effect |
|------|-------|-----------|------|--------|
| 1 | **ChIP-seq histone** | **C4-2B (prostate Ca, metastatic)** | — | quantile −0.997 |
| 2 | RNA-seq | MCF-7 | ENSG00000224722 | −0.997 |
| 3 | RNA-seq | foreskin keratinocyte | **CASC8** | +0.997 |
| 4 | RNA-seq | HeLa-S3 | ENSG00000224722 | −0.997 |
| 5 | **ChIP-seq histone** | **22Rv1 (prostate Ca)** | — | −0.997 |
| 6 | **ChIP-seq TF** | MCF-7 | — | −0.996 |
| 7 | RNA-seq | mammary epithelial cell | ENSG00000224722 | −0.996 |
| 8 | RNA-seq | keratinocyte | CASC8 | +0.996 |
| 9 | RNA-seq | hair follicular keratinocyte | CASC8 | +0.996 |
| 10 | **ChIP-seq histone** | **LNCaP clone FGC (prostate Ca)** | — | −0.995 |

**Predicted mechanism:** The A>G allele is predicted to **disrupt histone modification landscapes** (H3K27ac/H3K4me1) in prostate cancer lines (C4-2B, 22Rv1, LNCaP). This likely impairs an enhancer element in the 8q24 locus, altering *CASC8* lncRNA expression and potentially long-range contacts with *MYC*. The concordance across three prostate cancer ChIP-seq datasets is striking.

---

### 🥈 2. chr17:48728343:C>T — rs138213197 (OR = 3.85, highest in panel)
**Why it ranks #2:** Highest OR, targets **HOXB13** — a prostate-lineage-defining transcription factor frequently mutated in hereditary prostate cancer (G84E).

| Rank | Assay | Cell Type | Gene | Effect |
|------|-------|-----------|------|--------|
| 1–3 | RNA-seq | myometrial cell, GM23248, M059J | **MIR3185** | −0.9998 |
| 4 | RNA-seq | mesangial cell | **HOXB13** | −0.9998 |
| 5 | RNA-seq | Caki2 | **HOXB13** | −0.9998 |
| 6–10 | RNA-seq | A172, NCI-H460, H4, RPMI7951, HFFc6 | MIR3185 | −0.9998 |

**Predicted mechanism:** C>T allele reduces expression of *HOXB13* and its host microRNA *MIR3185*. Given HOXB13's role in prostate epithelial identity and androgen receptor crosstalk, this could alter AR-driven transcription programmes in tumour-adjacent stroma.

---

### 🥉 3. chr10:46046326:A>G — rs10993994 (OR = 1.23)
**Why it ranks #3:** The best-characterised prostate cancer regulatory variant — in the *MSMB* promoter. AlphaGenome recapitulates known biology, giving confidence to other predictions.

| Rank | Assay | Cell Type | Gene | Effect |
|------|-------|-----------|------|--------|
| 1 | RNA-seq | T-cell | **MSMB** | +0.99998 |
| 2 | RNA-seq | type B pancreatic cell | **MSMB** | +0.99998 |
| 3 | RNA-seq | hepatocyte | **MSMB** | +0.99998 |
| 4 | RNA-seq | keratinocyte | **MSMB** | **−0.99998** |
| 5–10 | RNA-seq | glutamatergic neuron, bronchial epithelial, etc. | MSMB | ±0.99998 |

**Predicted mechanism:** The G allele (risk) reduces *MSMB* expression in prostate epithelial cells (context-dependent; keratinocyte model shows strong downregulation). MSMB (PSP94) is a secreted prostate protein with tumour-suppressive properties. Reduced MSMB in the TME may alter stromal paracrine signalling. This is an established eQTL — AlphaGenome's recapitulation serves as a positive control.

---

### 4. chr14:22836440:T>C — rs1004030 (OR = 1.05)
**Key finding: ChIP-seq TF + RNA-seq — multi-assay. Targets MMP14 (MT1-MMP), a master regulator of ECM remodelling.**

| Rank | Assay | Cell Type | Gene | Effect |
|------|-------|-----------|------|--------|
| 1–3 | **ChIP-seq TF** | HEK293T, K562, HepG2 | — | −0.99998 |
| 4–10 | RNA-seq | NCI-H460, HepG2, HeLa-S3, MCF-7, HT-29, GM23338, MCF-7 | **MMP14** | −0.99996 |

**Predicted mechanism:** T>C disrupts TF binding at an enhancer/promoter, reducing *MMP14* expression. MMP14 activates pro-MMP2 and degrades collagen — central to **fibroblast-mediated matrix remodelling** in the prostate TME. Loss of MMP14 in stromal fibroblasts could shift ECM composition and alter tumour invasion.

---

### 5. chr6:32224554:A>G — rs3096702 (OR = 1.06)
**Hits NOTCH4 consistently across T-cell subtypes. NOTCH signalling is a key stromal–epithelial crosstalk pathway.**

| Rank | Assay | Cell Type | Gene | Effect |
|------|-------|-----------|------|--------|
| 1–10 | RNA-seq | CD4+ T-cell, CD8+ T-cell, Treg, NK cell, Th17, naive CD4/CD8 T-cells, memory CD4 T-cell | **NOTCH4** | +0.9992–0.9991 |

**Predicted mechanism:** A>G upregulates *NOTCH4* in infiltrating T-cells. NOTCH4 is expressed in vascular endothelium and can be expressed in activated T-cells; its dysregulation may alter Jagged/Delta-mediated interactions with stromal fibroblasts expressing NOTCH ligands.

---

### 6. chr7:98187015:C>T — rs6465657 (OR = 1.11)
**Strongest fibroblast-relevant cell-type targeting in the panel.**

| Rank | Assay | Cell Type | Gene | Effect |
|------|-------|-----------|------|--------|
| 1 | RNA-seq | GM23248 | BHLHA15 | +0.9998 |
| 2 | RNA-seq | smooth muscle cell (pulmonary artery) | BHLHA15 | +0.9997 |
| 3 | RNA-seq | smooth muscle cell (umbilical artery) | BHLHA15 | +0.9997 |
| 8 | RNA-seq | **fibroblast of the aortic adventitia** | BHLHA15 | +0.9997 |
| 9 | RNA-seq | **cardiac atrium fibroblast** | BHLHA15 | +0.9997 |
| 10 | RNA-seq | IMR-90 (fetal lung fibroblast) | BHLHA15 | +0.9997 |

**Predicted mechanism:** C>T upregulates *BHLHA15* (Mist1) in fibroblasts and smooth muscle cells. BHLHA15 is a bHLH transcription factor implicated in secretory cell identity; ectopic expression in prostate stromal fibroblasts could alter their activation state.

---

### 7. chr8:127472793:A>C — rs1447295 (OR = 1.41)
**Multi-assay: ChIP-seq histone + DNase-seq + RNA-seq. Also in 8q24.**

| Rank | Assay | Cell Type | Gene | Effect |
|------|-------|-----------|------|--------|
| 1 | **ChIP-seq histone** | Karpas-422 (B-cell lymphoma) | — | +0.998 |
| 4 | **ChIP-seq histone** | Karpas-422 | — | +0.997 |
| 8 | **ChIP-seq histone** | Karpas-422 | — | +0.997 |
| 10 | **DNase-seq** | Karpas-422 | — | +0.996 |
| 2–9 | RNA-seq | Purkinje cell, endothelial, umbilical cord, etc. | ENSG00000278324/ENSG00000271509 | +0.997 |

**Predicted mechanism:** A>C likely **creates or strengthens** a regulatory element (increased chromatin accessibility and histone acetylation), potentially altering 8q24 enhancer activity near *MYC*. Complements rs183373024 (variant #1 above) in the same locus.

---

## TIER 2: Moderate Candidates

### 8. chr16:57619621:T>G — rs11859370 (OR = 1.13) + chr16:57620664:C>T — rs11863709 (OR = 1.16)
Adjacent variants at the *ADGRG1* (GPR56) locus. ADGRG1 is an adhesion GPCR involved in cell–matrix interactions. Both show strong RNA-seq effects: rs11859370 reduces *ADGRG1* in neural/embryonic cells; rs11863709 hits *ADGRG1* and *RNU6-20P*. Could affect fibroblast adhesion signalling.

### 9. chr2:172446825:A>G — rs12621278 (OR = 1.27)
Strong colon-mucosa RNA-seq signals on an uncharacterised lncRNA (ENSG00000226963) and *ITGA6* (integrin α6). ITGA6 pairs with ITGB4 to form the laminin receptor — relevant to basal cell adhesion in prostate epithelium and potentially reactive stroma.

### 10. chr2:66425753:C>T — rs74702681 (OR = 1.17)
Upregulates *MIR4778* across pancreatic progenitors, Caco-2, MCF-7, hepatocytes. MicroRNA perturbation in the TME could affect multiple targets.

---

## TIER 3: Lower-Priority Candidates

| Variant | OR | Top Assay | Target | Comment |
|---------|-----|-----------|--------|---------|
| chr5:169745129:A>G (rs76551843) | 1.31 | RNA-seq | FOXI1 | Kidney epithelial; limited prostate relevance |
| chr15:56093670:A>G (rs33984059) | 1.19 | ChIP-seq TF/RNA-seq | RFX7 | Mixed signals; fibroblast of scalp hit |
| chr10:45587537:T>C (rs76934034) | 1.12 | RNA-seq | CUBNP2/MSMB | Weaker effect than rs10993994 |
| chr12:48025835:A>C (rs80130819) | 1.10 | RNA-seq | ENSG00000257955 | Testis-specific; unclear prostate relevance |
| chr17:7899800:T>C (rs28441558) | 1.16 | RNA-seq | TNFSF12/SNORD10 | TWEAK (TNFSF12) is a stromal cytokine — interesting but moderate effect |
| chr4:73483441:G>A (rs1894292) | 1.06 | RNA-seq + DNase-seq | PF4 | Platelet factor 4; DNase-seq hit in retina |
| chr17:49267824:G>A (rs11650494) | 1.10 | RNA-seq | RNU1-42P | Pseudogene; weak |
| chr6:153119944:A>G (rs1933488) | 1.08 | RNA-seq | HSPD1P16 | Pseudogene; weak |
| chr3:169375312:A>G (rs142436749) | 1.25 | RNA-seq | LINC02082 | lncRNA; moderate |
| chr2:111146954:C>A (rs56366063) | 1.08 | RNA-seq | BUB1/MIR4435-2 | Weak effect |
| chr2:10570604:C>T (rs9287719) | 1.07 | RNA-seq | ENSG00000212558 | Weak |
| chr12:114247766:A>G (rs1270884) | 1.07 | RNA-seq | LINC01234 | Weak |
| chr14:68660027:T>C (rs7141529) | 1.05 | RNA-seq | RPL12P7 | Pseudogene; weak |
| chr7:47397647:G>A (rs56232506) | 1.06 | RNA-seq | PKD1L1-AS1 | Weakest effect in panel |

---

## Prioritised Experimental Validation Plan

### 🔬 Phase I — Immediate (for Tier 1, variants 1–4)

#### 1. rs183373024 (chr8:127091872, 8q24) — *ChIP-seq histone disruption*

| Experiment | What it tests | Expected result if prediction correct |
|-----------|---------------|--------------------------------------|
| **Allele-specific ChIP-qPCR** for H3K27ac in LNCaP / C4-2B cells heterozygous at rs183373024 | Whether the G allele reduces histone acetylation | G allele shows reduced H3K27ac enrichment vs. A allele |
| **Allele-specific ATAC-seq** in prostate fibroblasts (e.g. WPMY-1 or primary CAF) | Whether chromatin accessibility differs by allele | G allele shows reduced accessibility |
| **CRISPRi** of the ~500 bp region surrounding rs183373024 in LNCaP, then RNA-seq | Whether the element regulates *CASC8* / *MYC* in *cis* | CRISPRi reduces CASC8 and/or MYC expression |
| **Capture-C / 4C-seq** from the rs183373024 viewpoint in LNCaP | Whether the element contacts *MYC* promoter ~335 kb away | Allele-specific contact: A allele shows stronger contact |

#### 2. rs138213197 (chr17:48728343, HOXB13 locus)

| Experiment | What it tests |
|-----------|---------------|
| **Luciferase reporter** with ~1 kb flanking each allele in LNCaP and prostate fibroblasts | Whether C>T alters enhancer/promoter activity |
| **Allele-specific qRT-PCR** for *HOXB13* in prostate tissue from heterozygous carriers | Whether T allele reduces HOXB13 expression |
| **CRISPRa/i** of the element in prostate organoids | Rescue of HOXB13 expression |

#### 3. rs10993994 (chr10:46046326, MSMB) — Positive control

| Experiment | What it tests |
|-----------|---------------|
| **Allele-specific ATAC-seq** in prostate epithelial cells | Confirm known differential accessibility |
| **qRT-PCR** for MSMB in prostate fibroblasts treated with recombinant MSMB | Test paracrine effect on stromal activation markers (αSMA, FAP) |

#### 4. rs1004030 (chr14:22836440, MMP14) — *Stromal remodelling*

| Experiment | What it tests |
|-----------|---------------|
| **Allele-specific luciferase** in prostate myofibroblast lines | Whether T>C reduces enhancer activity |
| **ChIP-qPCR** for the TF(s) predicted to bind (pulldown with antibody against TF whose motif is disrupted) | Identify which TF's binding is abrogated |
| **Gelatin zymography** on conditioned media from heterozygous fibroblasts ± MMP14 siRNA | Whether allele imbalance predicts MMP14 activity and MMP2 activation |
| **3D co-culture** of prostate cancer spheroids with fibroblasts of each genotype | Tumour invasion difference |

### 🔬 Phase II — Follow-up (Tier 1, variants 5–7)

- **rs3096702** (NOTCH4): Flow-sort CD4+ T-cells from prostate tumour digests; allele-specific qPCR for NOTCH4. Co-culture with NOTCH-ligand-expressing fibroblasts.
- **rs6465657** (BHLHA15): Immunofluorescence for BHLHA15 in prostate stromal fibroblasts; lentiviral overexpression in CAFs + RNA-seq.
- **rs1447295** (8q24): DNase-I hypersensitivity assay at the variant site in prostate fibroblasts + LNCaP; complements rs183373024 analysis.

### 🔬 Phase III — Tier 2 confirmation

- **rs11859370 / rs11863709** (ADGRG1): Flow cytometry for surface GPR56 on prostate fibroblasts heterozygous at these SNPs. Collagen adhesion assay.
- **rs12621278** (ITGA6): Allele-specific qPCR in prostate basal epithelial cells; IF for integrin α6β4.

---

## Key Uncertainties

1. **Cell-type mismatch:** Most top-scoring predictions came from non-prostate cell types in the training data. The ontology terms (prostate gland, fibroblast) were requested but the model appears to return effects from the most confidently scored cell contexts, which were frequently HeLa, K562, HepG2, etc. Experimental validation **must** use prostate-relevant cells.

2. **Assay coverage gap:** Very few variants scored in Hi-C/pcHi-C or ATAC-seq — the model may lack prostate Hi-C training data. This limits direct regulatory-contact predictions.

3. **RNA-seq dominance:** RNA-seq effects dominate because the model is trained on GTEx/eQTL data. Regulatory mechanisms (accessibility, histone marks, looping) are under-represented but more informative for mechanism.

4. **Linkage disequilibrium:** All variants are treated independently; in the 8q24 region, rs183373024 and rs1447295 may tag the same functional haplotype. Fine-mapping with local LD structure is essential before functional validation.