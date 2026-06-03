#!/usr/bin/env python3
"""
MOFA+ Multi-Omics Integration Experiment
=========================================
Prostate cancer risk SNP analysis integrating:
  - 1000G EUR genotypes (503 individuals × 88 PCa lead SNPs)
  - ENCODE prostate regulatory peaks (8 chromatin marks)
  - GTEx v8 prostate eQTL statistics
  - Schumacher 2018 GWAS summary statistics

Two experimental designs:
  1. SAMPLE-level MOFA: 503 individuals, Views = genotypes + regulatory burden scores
  2. SNP-level MOFA: 88 SNPs, Views = GWAS + regulatory + eQTL + population genetics
"""

import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = Path("/mnt/c/Users/rakti/Downloads/alphaGenome/mofa_data")
GENOTYPES_TSV    = BASE / "genotypes/1000g_eur_genotypes_pca_lead_snps.tsv"
REGULATORY_TSV   = BASE / "variant_scores/schumacher_pca_snps_encode_regulatory_annotation.tsv"
EQTL_TSV         = BASE / "variant_scores/schumacher2018_pca_lead_snps_gtex_prostate_eqtl_overlap.tsv"
GWAS_TSV         = BASE / "variant_scores/schumacher2018_pca_lead_snps_hg38.tsv"
OUT_DIR          = BASE / "mofa_results"
OUT_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("MOFA+ PROSTATE CANCER MULTI-OMICS INTEGRATION")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════════════
# 1. LOAD AND ALIGN DATA
# ═══════════════════════════════════════════════════════════════════════════
print("\n[1] Loading data layers...")

# 1a. Genotype matrix (503 EUR × 88 SNPs)
geno_df = pd.read_csv(GENOTYPES_TSV, sep='\t', index_col=0)
print(f"    Genotypes:       {geno_df.shape[0]} samples × {geno_df.shape[1]} SNPs")

# 1b. Regulatory annotation (136 SNPs × 8 marks)
reg_df = pd.read_csv(REGULATORY_TSV, sep='\t', index_col=0)
reg_marks = [c for c in reg_df.columns if c.startswith('in_')]
print(f"    Regulatory:      {reg_df.shape[0]} SNPs × {len(reg_marks)} marks")

# 1c. eQTL overlap (136 SNPs)
eqtl_df = pd.read_csv(EQTL_TSV, sep='\t', index_col=0)
print(f"    eQTL overlap:    {eqtl_df.shape[0]} SNPs")

# 1d. GWAS summary stats (136 SNPs) — keep both rsid and variant_id as columns
gwas_df = pd.read_csv(GWAS_TSV, sep='\t')
print(f"    GWAS stats:      {gwas_df.shape[0]} SNPs")

# ── Align SNPs across all datasets ─────────────────────────────────────────
# Genotypes columns = variant_id (e.g. chr8:127472793:A>C)
# Regulatory/eQTL index = variant_id
# GWAS has variant_id column (not index — index is rsid)
geno_snps = set(geno_df.columns)
reg_snps  = set(reg_df.index)
eqtl_snps = set(eqtl_df.index)
gwas_snps = set(gwas_df['variant_id'].values)  # use variant_id column

common_snps = sorted(geno_snps & reg_snps & eqtl_snps & gwas_snps)
print(f"\n    SNPs in all 4 datasets: {len(common_snps)}")

# Filter to common SNPs
geno_aligned = geno_df[common_snps]  # 503 × N_common
reg_aligned  = reg_df.loc[common_snps]
eqtl_aligned = eqtl_df.loc[common_snps]
gwas_aligned = gwas_df.set_index('variant_id').loc[common_snps]

# ═══════════════════════════════════════════════════════════════════════════
# 2. EXPERIMENT 1 — SAMPLE-LEVEL MOFA (503 1000G EUR individuals)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("[2] EXPERIMENT 1: Sample-Level MOFA (503 individuals)")
print("=" * 70)

# View 1: Genotype dosages (503 samples × K SNPs)
view1_geno = geno_aligned.values.astype(np.float32)
print(f"    View 1 (Genotypes): {view1_geno.shape}")

# View 2: Polygenic regulatory burden scores per ENCODE mark
# For each individual i and mark m:
#   burden[i,m] = Σ_snp dosage[i,snp] × in_peak[snp,m]
reg_binary = reg_aligned[reg_marks].values.astype(np.float32)
view2_reg_burden = view1_geno @ reg_binary  # 503 × 8
print(f"    View 2 (Regulatory burden): {view2_reg_burden.shape}")

# View 3: eQTL burden score per individual
#   burden[i] = Σ_snp dosage[i,snp] × n_eqtls[snp]
eqtl_counts = eqtl_aligned['n_gtex_prostate_eqtls_within_100kb'].values.astype(np.float32)
# Log-transform to reduce skew (add 1 pseudo-count)
eqtl_counts_log = np.log1p(eqtl_counts)
view3_eqtl_burden = (view1_geno @ eqtl_counts_log.reshape(-1, 1)).astype(np.float32)
print(f"    View 3 (eQTL burden): {view3_eqtl_burden.shape}")

# Sample-level view names
sample_views = {
    "Genotypes":      view1_geno,
    "Regulatory":     view2_reg_burden,
    "eQTL_burden":    view3_eqtl_burden,
}
sample_ids = list(geno_aligned.index)

# ═══════════════════════════════════════════════════════════════════════════
# 3. EXPERIMENT 2 — SNP-LEVEL MOFA (K SNPs)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("[3] EXPERIMENT 2: SNP-Level MOFA (88 SNPs)")
print("=" * 70)

# View 1: GWAS features
gwas_view1 = pd.DataFrame({
    'log10_pvalue': -np.log10(gwas_aligned['p_value'].values.astype(float)),
    'odds_ratio':   gwas_aligned['or'].values.astype(float),
    'risk_af':      gwas_aligned['risk_freq'].values.astype(float),
}, index=common_snps).values.astype(np.float32)
print(f"    View 1 (GWAS): {gwas_view1.shape}")

# View 2: Regulatory binary annotations (8 marks)
snp_view2 = reg_aligned[reg_marks].values.astype(np.float32)
# Clean column names
reg_mark_names = [c.replace('in_', '').replace('_peak', '') for c in reg_marks]
print(f"    View 2 (Regulatory): {snp_view2.shape} — {reg_mark_names}")

# View 3: eQTL features
eqtl_view3 = pd.DataFrame({
    'n_eqtls_log1p': np.log1p(eqtl_aligned['n_gtex_prostate_eqtls_within_100kb'].values.astype(float)),
}, index=common_snps)
# Add top eQTL p-value if available
top_pval = eqtl_aligned['top_eqtl_pval'].replace('', np.nan).values.astype(float)
eqtl_view3['top_eqtl_log10p'] = np.where(np.isfinite(-np.log10(top_pval)), -np.log10(top_pval), 0)
eqtl_view3 = eqtl_view3.values.astype(np.float32)
print(f"    View 3 (eQTL): {eqtl_view3.shape}")

# View 4: Population genetics from 1000G
# Allele frequency and variance per SNP in the 1000G EUR cohort
af_eur = view1_geno.mean(axis=0) / 2.0  # dosage 0/1/2 → AF
var_eur = view1_geno.var(axis=0)
het_eur = (view1_geno == 1).mean(axis=0)  # heterozygosity rate
snp_view4 = np.column_stack([af_eur, var_eur, het_eur]).astype(np.float32)
print(f"    View 4 (Population genetics): {snp_view4.shape} — AF, variance, het_rate")

snp_views = {
    "GWAS":                gwas_view1,
    "Regulatory":          snp_view2,
    "eQTL":                eqtl_view3,
    "Population_genetics": snp_view4,
}

# ═══════════════════════════════════════════════════════════════════════════
# 4. SAVE PREPARED MATRICES
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("[4] Saving prepared matrices...")

# Sample-level
np.savez(OUT_DIR / "sample_level_views.npz",
         genotypes=view1_geno,
         regulatory_burden=view2_reg_burden,
         eqtl_burden=view3_eqtl_burden,
         sample_ids=np.array(sample_ids, dtype=str))

# SNP-level
np.savez(OUT_DIR / "snp_level_views.npz",
         gwas=gwas_view1,
         regulatory=snp_view2,
         eqtl=eqtl_view3,
         popgen=snp_view4,
         snp_ids=np.array(common_snps, dtype=str),
         reg_mark_names=np.array(reg_mark_names, dtype=str))

# Metadata
pd.DataFrame({'sample_id': sample_ids}).to_csv(OUT_DIR / "sample_ids.tsv", sep='\t', index=False)
pd.DataFrame({'variant_id': common_snps}).to_csv(OUT_DIR / "snp_ids.tsv", sep='\t', index=False)

print("    ✓ sample_level_views.npz")
print("    ✓ snp_level_views.npz")
print("    ✓ sample_ids.tsv, snp_ids.tsv")

# ═══════════════════════════════════════════════════════════════════════════
# 5. RUN MOFA+ — SAMPLE LEVEL
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("[5] RUNNING MOFA+: Sample-Level Model (503 individuals)")
print("=" * 70)

from mofapy2.run.entry_point import entry_point

# Build sample-level MOFA input: nested list [view][group]
# Using a single group (all 503 samples together)
sample_data = [
    [view1_geno],           # View 0: Genotypes (503 × 88)
    [view2_reg_burden],     # View 1: Regulatory burden (503 × 8)
]

# Initialize and train
ent = entry_point()
ent.set_data_options(scale_views=False, center_groups=True)
ent.set_data_matrix(
    sample_data,
    likelihoods=["gaussian", "gaussian"],
    views_names=["Genotypes", "Regulatory_burden"],
)
ent.set_model_options(
    factors=6,
    spikeslab_weights=True,
    ard_weights=True,
)
ent.set_train_options(
    iter=2000,
    convergence_mode="fast",
    startELBO=1,
    freqELBO=200,
    dropR2=0.01,
    startDrop=50,
    freqDrop=10,
    verbose=False,
    quiet=True,
    seed=42,
    outfile=str(OUT_DIR / "sample_mofa.hdf5"),
)
ent.build()
ent.run()

# Extract results: Z is (N×K), W is list of (D_m×K) arrays
Z_sample = ent.model.nodes['Z'].getExpectation()
W_sample = ent.model.nodes['W'].getExpectation()
K_sample = Z_sample.shape[1]

# Compute R² using ORIGINAL data (sample_data modified in-place by MOFA)
sample_raw = [view1_geno, view2_reg_burden]  # raw matrices, not list-of-lists
sample_intercepts = [Y.mean(axis=0) for Y in sample_raw]
sample_r2_total = np.array([
    max(0, 1 - np.sum((sample_raw[m] - (sample_intercepts[m] + Z_sample @ W_sample[m].T))**2) /
         np.sum((sample_raw[m] - sample_intercepts[m])**2))
    for m in range(2)
])

print(f"\n    ✓ Sample-level MOFA trained: {K_sample} factors")
print(f"    Factors Z: {Z_sample.shape}")
view_labels_sample = ["Genotypes", "Regulatory_burden"]
for v, name in enumerate(view_labels_sample):
    print(f"    Weights {name}: {W_sample[v].shape}  R² = {sample_r2_total[v]:.3f}")

# Save
np.savez(OUT_DIR / "sample_mofa_factors.npz", Z=Z_sample, r2_per_view=sample_r2_total)
for v, name in enumerate(view_labels_sample):
    np.savetxt(OUT_DIR / f"sample_weights_{name}.tsv", W_sample[v], delimiter='\t')
np.savetxt(OUT_DIR / "sample_factor_matrix.tsv", Z_sample, delimiter='\t')

# ═══════════════════════════════════════════════════════════════════════════
# 6. RUN MOFA+ — SNP LEVEL
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("[6] RUNNING MOFA+: SNP-Level Model (88 SNPs)")
print("=" * 70)

# SNP-level: 88 SNPs as samples, 4 views
snp_data = [
    [gwas_view1],       # View 0: GWAS (88 × 3)
    [snp_view2],        # View 1: Regulatory binary (88 × 8)
    [eqtl_view3],       # View 2: eQTL (88 × 2)
    [snp_view4],        # View 3: Population genetics (88 × 3)
]

ent2 = entry_point()
ent2.set_data_options(scale_views=False, center_groups=True)
ent2.set_data_matrix(
    snp_data,
    likelihoods=["gaussian", "bernoulli", "gaussian", "gaussian"],
    views_names=["GWAS", "Regulatory", "eQTL", "Population_genetics"],
)
ent2.set_model_options(
    factors=5,
    spikeslab_weights=True,
    ard_weights=True,
)
ent2.set_train_options(
    iter=2000,
    convergence_mode="fast",
    startELBO=1,
    freqELBO=200,
    dropR2=0.01,
    startDrop=50,
    freqDrop=10,
    verbose=False,
    quiet=True,
    seed=42,
    outfile=str(OUT_DIR / "snp_mofa.hdf5"),
)
ent2.build()
ent2.run()

# Extract results: Z is (N×K), W is list of (D_m×K) arrays
Z_snp = ent2.model.nodes['Z'].getExpectation()
W_snp = ent2.model.nodes['W'].getExpectation()
K_snp = Z_snp.shape[1]

# Compute R² using ORIGINAL data (snp_data modified in-place by MOFA)
snp_raw = [gwas_view1, snp_view2, eqtl_view3, snp_view4]
snp_intercepts_raw = [Y.mean(axis=0) for Y in snp_raw]
snp_r2_total = np.array([
    max(0, 1 - np.sum((snp_raw[m] - (snp_intercepts_raw[m] + Z_snp @ W_snp[m].T))**2) /
         np.sum((snp_raw[m] - snp_intercepts_raw[m])**2))
    for m in range(4)
])

print(f"\n    ✓ SNP-level MOFA trained: {K_snp} factors")
print(f"    Factors Z: {Z_snp.shape}")
snp_view_names = ["GWAS", "Regulatory", "eQTL", "Population_genetics"]
for v, name in enumerate(snp_view_names):
    print(f"    Weights {name}: {W_snp[v].shape}  R² = {snp_r2_total[v]:.3f}")

# Save
np.savez(OUT_DIR / "snp_mofa_factors.npz", Z=Z_snp, r2_per_view=snp_r2_total)
for v, name in enumerate(snp_view_names):
    np.savetxt(OUT_DIR / f"W_snp_{name}.tsv", W_snp[v], delimiter='\t')
np.savetxt(OUT_DIR / "snp_factor_matrix.tsv", Z_snp, delimiter='\t')

# ═══════════════════════════════════════════════════════════════════════════
# 7. VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("[7] VISUALIZATION")
print("=" * 70)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import spearmanr

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})

# ── 7a. Variance Explained ──────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Sample-level total R² per view
colors_sample = ['#2196F3', '#4CAF50']
axes[0].bar(view_labels_sample, sample_r2_total, color=colors_sample, alpha=0.85, edgecolor='white')
for i, v in enumerate(sample_r2_total):
    axes[0].text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=10, fontweight='bold')
axes[0].set_ylabel('Variance Explained (R²)')
axes[0].set_title(f'Sample-Level MOFA\nN=503, K={K_sample} factors')
axes[0].set_ylim(0, max(0.5, sample_r2_total.max() * 1.2))

# SNP-level total R² per view
colors_snp = ['#E91E63', '#009688', '#FF5722', '#673AB7']
axes[1].bar(snp_view_names, snp_r2_total, color=colors_snp, alpha=0.85, edgecolor='white')
for i, v in enumerate(snp_r2_total):
    axes[1].text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=10, fontweight='bold')
axes[1].set_ylabel('Variance Explained (R²)')
axes[1].set_title(f'SNP-Level MOFA\n{len(common_snps)} SNPs, K={K_snp} factors')
axes[1].set_ylim(0, max(0.5, snp_r2_total.max() * 1.2))
axes[1].tick_params(axis='x', rotation=30)

plt.tight_layout()
fig.savefig(OUT_DIR / 'variance_explained.png', bbox_inches='tight')
plt.close()
print("    ✓ variance_explained.png")

# ── 7b. Factor Heatmaps ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Sample-level factor heatmap
im0 = axes[0].imshow(Z_sample.T, aspect='auto', cmap='RdBu_r',
                      norm=Normalize(vmin=-3, vmax=3))
axes[0].set_xlabel('Sample (1000G EUR individual)')
axes[0].set_ylabel('Factor')
axes[0].set_title(f'Sample-Level Factors ({K_sample} factors × {Z_sample.shape[0]} samples)')
plt.colorbar(im0, ax=axes[0], shrink=0.8, label='Factor Value')

# SNP-level factor heatmap
im1 = axes[1].imshow(Z_snp.T, aspect='auto', cmap='RdBu_r',
                      norm=Normalize(vmin=-3, vmax=3))
axes[1].set_xlabel('SNP index')
axes[1].set_ylabel('Factor')
axes[1].set_title(f'SNP-Level Factors ({K_snp} factors × {len(common_snps)} SNPs)')
plt.colorbar(im1, ax=axes[1], shrink=0.8, label='Factor Value')

plt.tight_layout()
fig.savefig(OUT_DIR / 'factor_heatmaps.png', bbox_inches='tight')
plt.close()
print("    ✓ factor_heatmaps.png")

# ── 7c. Sample-level: Factor 1 vs Factor 2 scatter ──────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc0 = axes[0].scatter(Z_sample[:, 0], Z_sample[:, 1],
                       c=view3_eqtl_burden.flatten(), cmap='viridis',
                       alpha=0.5, s=8, edgecolors='none')
axes[0].set_xlabel('Factor 1')
axes[0].set_ylabel('Factor 2')
axes[0].set_title('Sample-Level: Factor 1 vs Factor 2')
plt.colorbar(sc0, ax=axes[0], label='eQTL burden score')

sc1 = axes[1].scatter(Z_snp[:, 0], Z_snp[:, 1],
                       c=gwas_view1[:, 0], cmap='plasma',
                       alpha=0.8, s=40, edgecolors='black', linewidth=0.3)
axes[1].set_xlabel('Factor 1')
axes[1].set_ylabel('Factor 2')
axes[1].set_title('SNP-Level: Factor 1 vs Factor 2')
plt.colorbar(sc1, ax=axes[1], label='-log10(GWAS p-value)')

plt.tight_layout()
fig.savefig(OUT_DIR / 'factor_scatters.png', bbox_inches='tight')
plt.close()
print("    ✓ factor_scatters.png")

# ── 7d. SNP-level: Weight heatmap across views ──────────────────────────
fig, axes = plt.subplots(1, len(snp_view_names), figsize=(4 * len(snp_view_names), 6))

for v, (name, weights) in enumerate(zip(snp_view_names, W_snp)):
    w = weights  # features × factors
    im = axes[v].imshow(w.T, aspect='auto', cmap='RdBu_r', norm=Normalize(vmin=-2, vmax=2))
    axes[v].set_title(f'{name}\nweights')
    axes[v].set_xlabel('Feature')
    axes[v].set_ylabel('Factor')

    # Add feature labels for regulatory view
    if name == 'Regulatory':
        axes[v].set_xticks(range(len(reg_mark_names)))
        axes[v].set_xticklabels(reg_mark_names, rotation=45, ha='right', fontsize=7)

plt.tight_layout()
fig.savefig(OUT_DIR / 'snp_weight_heatmaps.png', bbox_inches='tight')
plt.close()
print("    ✓ snp_weight_heatmaps.png")

# ── 7e. SNP clustering by factor profile ────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 8))
row_linkage = linkage(Z_snp, method='ward')
col_linkage = linkage(Z_snp.T, method='ward')
row_order = dendrogram(row_linkage, no_plot=True)['leaves']
col_order = dendrogram(col_linkage, no_plot=True)['leaves']

im = ax.imshow(Z_snp[np.ix_(row_order, col_order)],
               aspect='auto', cmap='RdBu_r', norm=Normalize(vmin=-2.5, vmax=2.5))
ax.set_xlabel('Factor')
ax.set_ylabel(f'SNP (n={len(common_snps)})')
ax.set_title('SNP-Level Factors — Clustered Heatmap\n(Ward linkage on factor profiles)')
ax.set_xticks(range(len(col_order)))
ax.set_xticklabels([f'F{c+1}' for c in col_order])
plt.colorbar(im, ax=ax, shrink=0.8, label='Factor Value')
plt.tight_layout()
fig.savefig(OUT_DIR / 'snp_clustered_heatmap.png', bbox_inches='tight')
plt.close()
print("    ✓ snp_clustered_heatmap.png")

# ── 7f. Top SNPs per factor ─────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes = axes.flatten()

for k in range(min(K_snp, 6)):
    factor_loadings = Z_snp[:, k]
    top_idx = np.argsort(np.abs(factor_loadings))[-15:][::-1]
    top_snps = [common_snps[i] for i in top_idx]
    top_vals = factor_loadings[top_idx]

    colors = ['#E53935' if v > 0 else '#1E88E5' for v in top_vals]
    axes[k].barh(range(len(top_snps)), top_vals, color=colors, alpha=0.85, edgecolor='white')
    axes[k].set_yticks(range(len(top_snps)))
    axes[k].set_yticklabels([s.split(':')[0] for s in top_snps], fontsize=6)
    axes[k].set_xlabel('Factor Loading')
    axes[k].set_title(f'Factor {k+1} — Top 15 SNPs')
    axes[k].axvline(x=0, color='black', linewidth=0.5)
    axes[k].invert_yaxis()

plt.tight_layout()
fig.savefig(OUT_DIR / 'snp_top_loadings_per_factor.png', bbox_inches='tight')
plt.close()
print("    ✓ snp_top_loadings_per_factor.png")

# ═══════════════════════════════════════════════════════════════════════════
# 8. INTERPRETATION SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("[8] INTERPRETATION SUMMARY")
print("=" * 70)

# Sample-level: variance explained per view
print("\n─── Sample-Level MOFA (503 individuals) ───")
if sample_r2_total is not None:
    for v, name in enumerate(view_labels_sample):
        total_v = sample_r2_total[v]
        print(f"  {name:20s}: total R² = {total_v:.4f}")

# SNP-level: variance explained per view
print("\n─── SNP-Level MOFA (88 SNPs) ───")
if snp_r2_total is not None:
    for v, name in enumerate(snp_view_names):
        total_v = snp_r2_total[v]
        print(f"  {name:25s}: total R² = {total_v:.4f}")

# Top SNPs per factor (SNP-level)
print("\n─── SNP-Level: Top contributing SNPs per factor ───")
for k in range(K_snp):
    factor_loadings = Z_snp[:, k]
    top_idx = np.argsort(np.abs(factor_loadings))[-5:][::-1]
    print(f"\n  Factor {k+1}:")
    for i in top_idx:
        snp_id = common_snps[i]
        rsid = gwas_aligned.loc[snp_id, 'rsid'] if snp_id in gwas_aligned.index else '?'
        direction = '↑' if factor_loadings[i] > 0 else '↓'
        print(f"    {direction} {rsid:20s}  ({snp_id})  loading={factor_loadings[i]:+.3f}")

# Correlate sample factors with regulatory burden
print("\n─── Sample-Level: Factor–Regulatory correlations (Spearman ρ) ───")
for k in range(K_sample):
    for v, vname in enumerate(['Regulatory_burden', 'eQTL_burden']):
        if v == 0:
            burden = view2_reg_burden.mean(axis=1)  # mean across marks
        else:
            burden = view3_eqtl_burden.flatten()
        rho, pval = spearmanr(Z_sample[:, k], burden)
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
        if abs(rho) > 0.1:
            print(f"  Factor {k+1} vs {vname}: ρ = {rho:+.3f} (p={pval:.2e}) {sig}")

# ═══════════════════════════════════════════════════════════════════════════
# 9. EXPORT INTERPRETABLE TABLES
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("[9] Exporting interpretable tables...")

# SNP factor loadings with annotations
snp_loadings_df = pd.DataFrame(
    Z_snp,
    index=common_snps,
    columns=[f'Factor_{i+1}' for i in range(K_snp)]
)
# Add metadata
snp_loadings_df['rsid'] = [gwas_aligned.loc[s, 'rsid'] if s in gwas_aligned.index else '' for s in common_snps]
snp_loadings_df['chr'] = [gwas_aligned.loc[s, 'chrom'] if s in gwas_aligned.index else '' for s in common_snps]
snp_loadings_df['pos_hg38'] = [gwas_aligned.loc[s, 'pos_hg38'] if s in gwas_aligned.index else np.nan for s in common_snps]
snp_loadings_df['-log10p'] = gwas_view1[:, 0]
snp_loadings_df['OR'] = gwas_view1[:, 1]
snp_loadings_df['n_regulatory_marks'] = snp_view2.sum(axis=1)
snp_loadings_df['n_eqtls'] = eqtl_aligned['n_gtex_prostate_eqtls_within_100kb'].values

# Sort by absolute loading on Factor 1
snp_loadings_df = snp_loadings_df.sort_values('Factor_1', key=abs, ascending=False)
snp_loadings_df.to_csv(OUT_DIR / 'snp_factor_loadings_annotated.tsv', sep='\t')
print(f"    ✓ snp_factor_loadings_annotated.tsv ({len(snp_loadings_df)} SNPs)")

# Sample factor loadings
sample_loadings_df = pd.DataFrame(
    Z_sample,
    index=sample_ids,
    columns=[f'Factor_{i+1}' for i in range(K_sample)]
)
sample_loadings_df['regulatory_burden_mean'] = view2_reg_burden.mean(axis=1)
sample_loadings_df['eqtl_burden'] = view3_eqtl_burden.flatten()
sample_loadings_df.to_csv(OUT_DIR / 'sample_factor_loadings_annotated.tsv', sep='\t')
print(f"    ✓ sample_factor_loadings_annotated.tsv ({len(sample_loadings_df)} samples)")

# Summary statistics
summary = {
    'experiment': ['Sample-level', 'Sample-level', 'SNP-level', 'SNP-level'],
    'metric': ['n_samples', 'n_factors', 'n_snps', 'n_factors'],
    'value': [Z_sample.shape[0], K_sample, len(common_snps), K_snp],
}
summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUT_DIR / 'experiment_summary.tsv', sep='\t', index=False)
print("    ✓ experiment_summary.tsv")

print("\n" + "=" * 70)
print("MOFA+ EXPERIMENT COMPLETE")
print(f"Results saved to: {OUT_DIR}")
print("=" * 70)
