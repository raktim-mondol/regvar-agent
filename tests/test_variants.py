"""Tests for the pure-Python parts: parsing, coordinate conventions, ranking.

These run without an API key or network -- they exercise exactly the logic that
silently corrupts variant-effect analyses when it is wrong.

    pip install pytest pandas
    pytest -q
"""

import pandas as pd
import pytest

from regvar.variants import (
    CandidateVariant,
    bed_to_one_based,
    rank_by_effect,
    read_candidates_tsv,
)


def test_bed_to_one_based():
    assert bed_to_one_based(0) == 1
    assert bed_to_one_based(127401059) == 127401060


def test_candidate_vcf_id():
    v = CandidateVariant("chr8", 127401060, "G", "T")
    assert v.vcf_id == "chr8:127401060:G>T"
    assert v.as_dict() == {
        "chromosome": "chr8", "position": 127401060, "ref": "G", "alt": "T",
    }


def test_read_candidates_tsv(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text(
        "chromosome\tposition\tref\talt\tregion_id\tnote\n"
        "chr8\t127401060\tg\tt\tenh1\thello\n"
    )
    variants = read_candidates_tsv(p)
    assert len(variants) == 1
    # alleles upper-cased, position kept 1-based and integer
    assert variants[0].ref == "G" and variants[0].alt == "T"
    assert variants[0].position == 127401060
    assert variants[0].region_id == "enh1"


def test_read_candidates_missing_column(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("chromosome\tposition\tref\n chr8\t1\tG\n")
    with pytest.raises(ValueError):
        read_candidates_tsv(p)


def test_rank_by_effect_orders_by_abs_quantile():
    df = pd.DataFrame({
        "gene_name": ["A", "B", "C"],
        "raw_score": [0.1, -0.9, 0.2],
        "quantile_score": [0.10, -0.95, 0.30],
    })
    ranked = rank_by_effect(df, top_n=2)
    assert list(ranked["gene_name"]) == ["B", "C"]
    assert ranked.iloc[0]["abs_effect"] == pytest.approx(0.95)


def test_rank_by_effect_falls_back_to_raw():
    df = pd.DataFrame({"gene_name": ["A", "B"], "raw_score": [0.2, -0.5]})
    ranked = rank_by_effect(df)
    assert list(ranked["gene_name"]) == ["B", "A"]
