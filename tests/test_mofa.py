"""Tests for the MOFA+ view builder and integration layer.

Pure-Python unit tests — no API calls, no mofapy2, no AlphaGenome.
Follows the same pattern as tests/test_variants.py.

    pip install pytest pandas numpy
    pytest -q tests/test_mofa.py
"""

import numpy as np
import pandas as pd
import pytest

from regvar.mofa_view import (
    build_mofa_view,
    read_genotypes_tsv,
    scored_variants_from_tsv,
    _sanitise_view_name,
    _feature_name,
)


# ---------------------------------------------------------------------------
# Helpers — build small, realistic fixture DataFrames
# ---------------------------------------------------------------------------

def _make_score_df(
    variant_id: str = "chr8:127401060:G>T",
    assays: list[str] | None = None,
    include_quantile: bool = True,
) -> pd.DataFrame:
    """Build a minimal score DataFrame matching AlphaGenomeClient output."""
    if assays is None:
        assays = ["ATAC-seq", "RNA-seq", "ChIP-seq histone", "Hi-C / pcHi-C"]

    rows = []
    for assay in assays:
        rows.append({
            "query_variant": variant_id,
            "variant_id": variant_id,
            "gene_name": f"GENE_{assay[:4]}",
            "output_type": f"OUT_{assay[:6]}",
            "biosample_name": "prostate",
            "raw_score": 0.5,
            "quantile_score": 0.75 if include_quantile else np.nan,
            "assay": assay,
            "abs_effect": 0.75,
        })
    return pd.DataFrame(rows)


def _make_multi_track_score_df(variant_id: str = "chr8:127401060:G>T") -> pd.DataFrame:
    """Score DF with multiple tracks per assay (e.g. multiple genes for RNA-seq)."""
    rows = [
        # Three RNA-seq tracks with different genes/scores
        {"query_variant": variant_id, "gene_name": "MYC",  "output_type": "RNA_SEQ", "raw_score": 0.1, "quantile_score": 0.10, "assay": "RNA-seq"},
        {"query_variant": variant_id, "gene_name": "MYC",  "output_type": "RNA_SEQ", "raw_score": -0.5, "quantile_score": -0.50, "assay": "RNA-seq"},
        {"query_variant": variant_id, "gene_name": "BRCA1","output_type": "RNA_SEQ", "raw_score": 0.9, "quantile_score": 0.95, "assay": "RNA-seq"},
        # One ATAC track
        {"query_variant": variant_id, "gene_name": None,   "output_type": "ATAC",    "raw_score": 0.3, "quantile_score": 0.30, "assay": "ATAC-seq"},
    ]
    return pd.DataFrame(rows)


def _make_genotype_df(
    variant_ids: list[str] | None = None,
    n_samples: int = 5,
    with_missing: bool = True,
) -> pd.DataFrame:
    """Build a genotype DataFrame with realistic 0/1/2/NaN values."""
    if variant_ids is None:
        variant_ids = ["chr8:127401060:G>T"]

    rng = np.random.default_rng(42)
    data: dict[str, list[float]] = {}
    for vid in variant_ids:
        dosages = rng.choice([0.0, 1.0, 2.0], size=n_samples).tolist()
        if with_missing:
            dosages[0] = np.nan   # first sample always missing
        data[vid] = dosages

    sample_ids = [f"SAMPLE_{i:02d}" for i in range(1, n_samples + 1)]
    return pd.DataFrame(data, index=pd.Index(sample_ids, name="sample_id"))


# ---------------------------------------------------------------------------
# Sanitisation & naming
# ---------------------------------------------------------------------------

def test_sanitise_view_name():
    assert _sanitise_view_name("ATAC-seq") == "ATAC-seq"
    assert _sanitise_view_name("Hi-C / pcHi-C") == "Hi-C_pcHi-C"
    assert _sanitise_view_name("ChIP-seq histone") == "ChIP-seq_histone"


def test_feature_naming():
    name = _feature_name("chr8:127401060:G>T", "MYC", "ATAC")
    assert name == "chr8:127401060:G>T_MYC_ATAC"
    # Missing gene name → "NA"
    name2 = _feature_name("chr8:127401060:G>T", None, "RNA_SEQ")
    assert name2 == "chr8:127401060:G>T_NA_RNA_SEQ"


# ---------------------------------------------------------------------------
# Genotype I/O
# ---------------------------------------------------------------------------

def test_read_genotypes_tsv(tmp_path):
    p = tmp_path / "geno.tsv"
    p.write_text(
        "sample_id\tchr8:127401060:G>T\tchr10:46046326:A>G\n"
        "SAMP_01\t0\t1\n"
        "SAMP_02\t1\t2\n"
        "SAMP_03\t\t0\n"
    )
    df = read_genotypes_tsv(p)
    assert df.shape == (3, 2)
    assert df.index.name == "sample_id"
    assert list(df.index) == ["SAMP_01", "SAMP_02", "SAMP_03"]
    assert df.loc["SAMP_01", "chr8:127401060:G>T"] == 0.0
    assert df.loc["SAMP_02", "chr10:46046326:A>G"] == 2.0
    assert np.isnan(df.loc["SAMP_03", "chr8:127401060:G>T"])


def test_read_genotypes_tsv_missing_sample_id(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("variant1\tvariant2\n0\t1\n")
    with pytest.raises(ValueError, match="sample_id"):
        read_genotypes_tsv(p)


# ---------------------------------------------------------------------------
# Scores TSV I/O
# ---------------------------------------------------------------------------

def test_scored_variants_from_tsv(tmp_path):
    p = tmp_path / "scores.tsv"
    p.write_text(
        "query_variant\tgene_name\tassay\traw_score\tquantile_score\n"
        "chr8:127401060:G>T\tMYC\tATAC-seq\t0.5\t0.75\n"
        "chr8:127401060:G>T\tMYC\tRNA-seq\t0.3\t0.45\n"
        "chr10:46046326:A>G\tMSMB\tRNA-seq\t-0.8\t-0.90\n"
    )
    result = scored_variants_from_tsv(p)
    assert len(result) == 2  # two unique variants
    # First variant should be chr8
    assert result[0].iloc[0]["query_variant"] == "chr8:127401060:G>T"
    assert len(result[0]) == 2  # two rows
    # Second variant
    assert result[1].iloc[0]["query_variant"] == "chr10:46046326:A>G"
    assert len(result[1]) == 1


def test_scored_variants_from_tsv_missing_column(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("gene_name\tassay\nMYC\tATAC-seq\n")
    with pytest.raises(ValueError, match="query_variant"):
        scored_variants_from_tsv(p)


# ---------------------------------------------------------------------------
# Core: build_mofa_view — shapes & structure
# ---------------------------------------------------------------------------

def test_build_view_shapes():
    """Correct (N, D) per view with known input."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_score_df(variant_id)]
    geno = _make_genotype_df([variant_id], n_samples=5, with_missing=False)

    views = build_mofa_view(scored, geno, by="max_abs")
    assert len(views) == 4  # ATAC-seq, RNA-seq, ChIP-seq histone, Hi-C / pcHi-C
    for vname, vdf in views.items():
        assert vdf.shape == (5, 1)  # 5 samples, 1 feature
        assert list(vdf.index) == [f"SAMPLE_{i:02d}" for i in range(1, 6)]
        assert vdf.index.name == "sample_id"
        # Feature name format
        col = vdf.columns[0]
        assert col.startswith(variant_id)


def test_build_view_with_assay_subset():
    """Only requested assays are included."""
    scored = [_make_score_df("chr8:127401060:G>T")]
    geno = _make_genotype_df(["chr8:127401060:G>T"], n_samples=3, with_missing=False)

    views = build_mofa_view(scored, geno, assays=["ATAC-seq", "RNA-seq"])
    assert set(views.keys()) == {"ATAC-seq", "RNA-seq"}


def test_build_view_unknown_assay_raises():
    """Passing an unrecognized assay raises ValueError."""
    scored = [_make_score_df("chr8:127401060:G>T")]
    geno = _make_genotype_df(["chr8:127401060:G>T"], n_samples=3)
    with pytest.raises(ValueError, match="Unknown assay"):
        build_mofa_view(scored, geno, assays=["NotAnAssay"])


# ---------------------------------------------------------------------------
# Dosage model
# ---------------------------------------------------------------------------

def test_dosage_zero():
    """0 alt copies → 0.0 feature value."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_score_df(variant_id)]
    geno = pd.DataFrame(
        {variant_id: [0, 0]},
        index=pd.Index(["S_01", "S_02"], name="sample_id"),
    )
    views = build_mofa_view(scored, geno)
    for vdf in views.values():
        assert (vdf.values == 0.0).all()


def test_dosage_one():
    """1 alt copy → quantile_score value."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_score_df(variant_id)]  # quantile_score = 0.75
    geno = pd.DataFrame(
        {variant_id: [1, 1]},
        index=pd.Index(["S_01", "S_02"], name="sample_id"),
    )
    views = build_mofa_view(scored, geno)
    for vdf in views.values():
        assert np.allclose(vdf.values, 0.75)


def test_dosage_two():
    """2 alt copies → 2 × quantile_score."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_score_df(variant_id)]  # quantile_score = 0.75
    geno = pd.DataFrame(
        {variant_id: [2, 2]},
        index=pd.Index(["S_01", "S_02"], name="sample_id"),
    )
    views = build_mofa_view(scored, geno)
    for vdf in views.values():
        assert np.allclose(vdf.values, 1.5)  # 2 × 0.75


def test_dosage_fallback_to_raw():
    """Falls back to raw_score when quantile_score absent."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_score_df(variant_id, include_quantile=False)]
    geno = pd.DataFrame(
        {variant_id: [1, 1]},
        index=pd.Index(["S_01", "S_02"], name="sample_id"),
    )
    views = build_mofa_view(scored, geno)
    for vdf in views.values():
        assert np.allclose(vdf.values, 0.5)  # raw_score = 0.5


def test_missing_genotype_is_nan():
    """NaN genotype → NaN in output matrix."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_score_df(variant_id)]
    geno = pd.DataFrame(
        {variant_id: [np.nan, 0]},
        index=pd.Index(["S_miss", "S_has"], name="sample_id"),
    )
    views = build_mofa_view(scored, geno)
    for vdf in views.values():
        assert np.isnan(vdf.loc["S_miss"].values[0])
        assert vdf.loc["S_has"].values[0] == 0.0


# ---------------------------------------------------------------------------
# Aggregation strategies
# ---------------------------------------------------------------------------

def test_aggregation_max_abs():
    """Multiple tracks → single strongest by |quantile_score|."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_multi_track_score_df(variant_id)]
    geno = pd.DataFrame(
        {variant_id: [1, 1]},
        index=pd.Index(["S_01", "S_02"], name="sample_id"),
    )
    # RNA-seq has 3 tracks: 0.1, -0.5, 0.95 → max_abs picks 0.95
    # ATAC-seq has 1 track: 0.3 → picks 0.3
    views = build_mofa_view(scored, geno, by="max_abs")

    rna_view = views["RNA-seq"]
    assert rna_view.shape[1] == 1  # single feature
    assert np.allclose(rna_view.values, 0.95)

    atac_view = views["ATAC-seq"]
    assert np.allclose(atac_view.values, 0.3)


def test_aggregation_all():
    """Multiple tracks → expanded features (no collapsing)."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_multi_track_score_df(variant_id)]
    geno = pd.DataFrame(
        {variant_id: [1, 1]},
        index=pd.Index(["S_01", "S_02"], name="sample_id"),
    )
    views = build_mofa_view(scored, geno, by="all")

    # RNA-seq should have 3 features
    rna_view = views["RNA-seq"]
    assert rna_view.shape[1] == 3  # three RNA-seq tracks
    # ATAC-seq has 1 feature
    atac_view = views["ATAC-seq"]
    assert atac_view.shape[1] == 1


def test_aggregation_top_gene_naming():
    """Feature name uses actual top-scoring gene name."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_multi_track_score_df(variant_id)]
    geno = pd.DataFrame(
        {variant_id: [1, 1]},
        index=pd.Index(["S_01", "S_02"], name="sample_id"),
    )
    views = build_mofa_view(scored, geno, by="top_gene")
    # RNA-seq top track is BRCA1 (quantile_score=0.95)
    rna_col = views["RNA-seq"].columns[0]
    assert "BRCA1" in rna_col


# ---------------------------------------------------------------------------
# Error conditions
# ---------------------------------------------------------------------------

def test_variant_genotype_mismatch_raises():
    """Scores without matching genotype column → ValueError."""
    scored = [_make_score_df("chr8:127401060:G>T")]
    geno = _make_genotype_df(["chr10:99999999:A>C"], n_samples=3)
    with pytest.raises(ValueError, match="missing from the genotype table"):
        build_mofa_view(scored, geno)


def test_single_sample_raises():
    """N=1 → error, MOFA+ needs ≥2."""
    scored = [_make_score_df("chr8:127401060:G>T")]
    geno = pd.DataFrame(
        {"chr8:127401060:G>T": [1]},
        index=pd.Index(["ONLY_ONE"], name="sample_id"),
    )
    with pytest.raises(ValueError, match="at least 2 samples"):
        build_mofa_view(scored, geno)


def test_no_scored_variants_raises():
    """Empty scored_variants list → ValueError."""
    geno = _make_genotype_df(["chr8:127401060:G>T"], n_samples=3)
    with pytest.raises(ValueError, match="No scored variants"):
        build_mofa_view([], geno)


def test_invalid_aggregation_strategy():
    """Bad 'by' value → ValueError."""
    scored = [_make_score_df("chr8:127401060:G>T")]
    geno = _make_genotype_df(["chr8:127401060:G>T"], n_samples=3)
    with pytest.raises(ValueError, match="Unknown aggregation strategy"):
        build_mofa_view(scored, geno, by="nonsense")


def test_all_nan_genotype_drops_feature():
    """When ALL samples have NaN for all variants, no views can be built."""
    variant_id = "chr8:127401060:G>T"
    scored = [_make_score_df(variant_id)]
    geno = pd.DataFrame(
        {variant_id: [np.nan, np.nan]},
        index=pd.Index(["S_01", "S_02"], name="sample_id"),
    )
    with pytest.raises(ValueError, match="No views could be built"):
        build_mofa_view(scored, geno)


# ---------------------------------------------------------------------------
# Feature naming uniqueness
# ---------------------------------------------------------------------------

def test_feature_naming_unique():
    """All feature names are unique across views."""
    scored = [
        _make_score_df("chr8:127401060:G>T"),
        _make_score_df("chr10:46046326:A>G"),
    ]
    geno = _make_genotype_df(
        ["chr8:127401060:G>T", "chr10:46046326:A>G"],
        n_samples=5, with_missing=False,
    )
    views = build_mofa_view(scored, geno, by="max_abs")
    all_features: list[str] = []
    for vdf in views.values():
        all_features.extend(vdf.columns.tolist())
    assert len(all_features) == len(set(all_features))


# ---------------------------------------------------------------------------
# Multiple variants, multiple assays integration
# ---------------------------------------------------------------------------

def test_build_mofa_view_multi_variant():
    """End-to-end: 2 variants × 4 assays = 8 features across 4 views."""
    v1 = "chr8:127401060:G>T"
    v2 = "chr10:46046326:A>G"
    scored = [_make_score_df(v1), _make_score_df(v2)]
    geno = _make_genotype_df([v1, v2], n_samples=5, with_missing=False)

    views = build_mofa_view(scored, geno, by="max_abs")
    assert len(views) == 4
    for vname, vdf in views.items():
        assert vdf.shape == (5, 2)  # 5 samples, 2 features (one per variant)
        assert vdf.index.name == "sample_id"
        # Both variant IDs should appear in column names
        cols_str = " ".join(vdf.columns)
        assert v1 in cols_str
        assert v2 in cols_str
