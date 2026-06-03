"""Tests for VCF input support.

These tests use inline VCF content written to a file, so they run without
any external fixtures or API keys.  cyvcf2 is a required dependency for
these tests (``pip install regvar[vcf]``).
"""

import pytest

cyvcf2 = pytest.importorskip("cyvcf2")

from regvar.variants import read_candidates_vcf

MINIMAL_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=chr8,length=145138631>
##contig=<ID=chr10,length=133797422>
##contig=<ID=17,length=83257441>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr8\t127401060\trs6983267\tG\tT\t.\tPASS\t.
chr10\t46046326\trs10993994\tA\tG\t.\tPASS\t.
chr8\t127472793\t.\tA\tC\t.\tPASS\t.
"""

VCF_NO_CHR_PREFIX = """\
##fileformat=VCFv4.2
##contig=<ID=8,length=145138631>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
8\t127401060\trs6983267\tG\tT\t.\tPASS\t.
"""

VCF_MULTI_ALLELIC = """\
##fileformat=VCFv4.2
##contig=<ID=chr17,length=83257441>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr17\t715725\trs684232\tT\tC,A\t.\tPASS\t.
"""


def test_read_vcf_basic(tmp_path):
    p = tmp_path / "test.vcf"
    p.write_text(MINIMAL_VCF)
    variants = read_candidates_vcf(p)

    assert len(variants) == 3
    assert variants[0].chromosome == "chr8"
    assert variants[0].position == 127401060
    assert variants[0].ref == "G"
    assert variants[0].alt == "T"
    assert variants[0].region_id == "rs6983267"


def test_read_vcf_missing_id_is_empty(tmp_path):
    p = tmp_path / "test.vcf"
    p.write_text(MINIMAL_VCF)
    variants = read_candidates_vcf(p)

    assert variants[2].region_id == ""


def test_read_vcf_auto_chr_prefix(tmp_path):
    p = tmp_path / "no_chr.vcf"
    p.write_text(VCF_NO_CHR_PREFIX)
    variants = read_candidates_vcf(p)

    assert len(variants) == 1
    assert variants[0].chromosome == "chr8"


def test_read_vcf_multi_allelic(tmp_path):
    p = tmp_path / "multi.vcf"
    p.write_text(VCF_MULTI_ALLELIC)
    variants = read_candidates_vcf(p)

    assert len(variants) == 2
    assert variants[0].ref == "T"
    assert variants[0].alt == "C"
    assert variants[1].ref == "T"
    assert variants[1].alt == "A"
    assert variants[0].region_id == "rs684232"
    assert variants[1].region_id == "rs684232"


def test_read_vcf_positions_are_1_based(tmp_path):
    p = tmp_path / "test.vcf"
    p.write_text(MINIMAL_VCF)
    variants = read_candidates_vcf(p)

    for v in variants:
        assert v.position > 0
