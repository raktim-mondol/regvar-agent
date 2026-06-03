"""Tests for the variant annotation module (regvar/annotation.py)."""

import pytest

from regvar.annotation import (
    VariantAnnotation,
    _ANNOTATIONS,
    annotate_nearest_gene,
    annotate_regulatory_overlap,
    annotate_variants,
    get_annotation,
)
from regvar.variants import CandidateVariant


@pytest.fixture(autouse=True)
def _clear_registry():
    """Reset the module-level annotation registry between tests."""
    _ANNOTATIONS.clear()
    yield
    _ANNOTATIONS.clear()


@pytest.fixture
def variants():
    return [
        CandidateVariant("chr8", 127401060, "G", "T"),
        CandidateVariant("chr8", 127410500, "A", "C"),
    ]


@pytest.fixture
def gtf_file(tmp_path):
    gtf = (
        "##gff-version 3\n"
        'chr8\tensembl\tgene\t127400000\t127410000\t.\t+\t.\tgene_id "ENSG1"; gene_name "MYC";\n'
        'chr8\tensembl\tgene\t127500000\t127510000\t.\t-\t.\tgene_id "ENSG2"; gene_name "TP53";\n'
    )
    p = tmp_path / "genes.gtf"
    p.write_text(gtf, encoding="utf-8")
    return p


@pytest.fixture
def bed_file(tmp_path):
    bed = (
        "chr8\t127401000\t127402000\tenh1\n"
        "chr8\t127410000\t127411000\tenh2\n"
    )
    p = tmp_path / "peaks.bed"
    p.write_text(bed, encoding="utf-8")
    return p


def test_annotate_nearest_gene(variants, gtf_file):
    result = annotate_nearest_gene(variants, gtf_file)
    ann = result["chr8:127401060:G>T"]
    assert ann.nearest_gene == "MYC"
    assert ann.tss_distance == 1060   # 127401060 - 127400000 (TSS on + strand)


def test_annotate_nearest_gene_minus_strand(variants, gtf_file):
    # Move the second variant close to TP53 (minus strand)
    vs = [CandidateVariant("chr8", 127510100, "A", "C")]
    result = annotate_nearest_gene(vs, gtf_file)
    ann = result["chr8:127510100:A>C"]
    assert ann.nearest_gene == "TP53"
    # TSS of TP53 on - strand is at end=127510000
    assert ann.tss_distance is not None


def test_annotate_regulatory_overlap(variants, bed_file):
    result = annotate_regulatory_overlap(variants, bed_file)
    ann1 = result["chr8:127401060:G>T"]
    assert "enh1" in ann1.regulatory_overlaps
    ann2 = result["chr8:127410500:A>C"]
    assert "enh2" in ann2.regulatory_overlaps


def test_annotate_regulatory_overlap_no_match(tmp_path):
    bed = tmp_path / "no.bed"
    bed.write_text("chr1\t0\t100\tpeak1\n", encoding="utf-8")
    result = annotate_regulatory_overlap(
        [CandidateVariant("chr8", 127401060, "G", "T")], bed,
    )
    assert result["chr8:127401060:G>T"].regulatory_overlaps == []


def test_annotate_variants_merges_sources(variants, gtf_file, bed_file):
    result = annotate_variants(variants, gtf_path=gtf_file, bed_path=bed_file)
    ann = result["chr8:127401060:G>T"]
    assert ann.nearest_gene == "MYC"
    assert "enh1" in ann.regulatory_overlaps
    assert get_annotation("chr8:127401060:G>T") is ann


def test_annotate_variants_empty_call_is_noop(variants):
    result = annotate_variants(variants)
    assert all(not ann.to_note() for ann in result.values())
    # Registry still populated (with empty annotations)
    assert get_annotation("chr8:127401060:G>T") is not None


def test_variant_annotation_to_note():
    ann = VariantAnnotation(
        variant=CandidateVariant("chr8", 100, "A", "G"),
        nearest_gene="MYC",
        tss_distance=-1500,
        regulatory_overlaps=["enh1", "enh2"],
        rsid="rs12345",
    )
    note = ann.to_note()
    assert "MYC" in note
    assert "-1500bp" in note
    assert "enh1" in note
    assert "rsid=rs12345" in note


def test_variant_annotation_to_note_empty():
    ann = VariantAnnotation(variant=CandidateVariant("chr1", 1, "A", "G"))
    assert ann.to_note() == ""


def test_gtf_gz_supported(tmp_path):
    import gzip
    gz = tmp_path / "genes.gtf.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as fh:
        fh.write(
            'chr1\tensembl\tgene\t100\t200\t.\t+\t.\tgene_name "FOO";\n'
        )
    vs = [CandidateVariant("chr1", 150, "A", "G")]
    result = annotate_nearest_gene(vs, gz)
    assert result["chr1:150:A>G"].nearest_gene == "FOO"
