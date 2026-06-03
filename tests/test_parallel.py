"""Tests for parallel AlphaGenomeClient.score_variants()."""

import time
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from regvar.alphagenome_client import AlphaGenomeClient, ClientConfig


def _make_tidy(variant_label: str) -> pd.DataFrame:
    return pd.DataFrame({
        "gene_name": ["GENE1"],
        "raw_score": [0.5],
        "assay": ["ATAC-seq"],
        "query_variant": [variant_label],
    })


def test_score_variants_serial_default():
    """max_concurrent=1 (default) runs calls sequentially and returns a
    concatenated DataFrame."""
    cfg = ClientConfig(min_seconds_between_calls=0)
    cfg.api_key = "dummy"
    client = AlphaGenomeClient(cfg)

    client._model = MagicMock()
    client._genome = MagicMock()
    client._variant_scorers = MagicMock()
    client._dna_client = MagicMock()

    call_order: list[int] = []

    def fake_score_variant(chromosome, position, ref, alt, **kw):
        call_order.append(position)
        label = f"{chromosome}:{position}:{ref}>{alt}"
        return _make_tidy(label)

    client.score_variant = fake_score_variant  # type: ignore[method-assign]

    variants = [
        {"chromosome": "chr1", "position": 100, "ref": "A", "alt": "G"},
        {"chromosome": "chr1", "position": 200, "ref": "C", "alt": "T"},
    ]
    df = client.score_variants(variants)

    assert len(df) == 2
    assert call_order == [100, 200]   # serial order


@patch("regvar.alphagenome_client.AlphaGenomeClient.score_variant")
def test_score_variants_parallel_dispatches(mock_score):
    """max_concurrent>1 dispatches through a ThreadPoolExecutor; results match."""
    mock_score.side_effect = lambda chromosome, position, ref, alt, **kw: _make_tidy(
        f"{chromosome}:{position}:{ref}>{alt}"
    )
    cfg = ClientConfig(min_seconds_between_calls=0, max_concurrent=4)
    cfg.api_key = "dummy"
    client = AlphaGenomeClient(cfg)

    variants = [
        {"chromosome": "chr1", "position": i, "ref": "A", "alt": "G"}
        for i in range(10)
    ]
    df = client.score_variants(variants)

    assert len(df) == 10
    assert mock_score.call_count == 10


@patch("regvar.alphagenome_client.AlphaGenomeClient.score_variant")
def test_score_variants_parallel_respects_semaphore(mock_score):
    """At most max_concurrent score_variant calls run simultaneously."""
    cfg = ClientConfig(min_seconds_between_calls=0, max_concurrent=2)
    cfg.api_key = "dummy"
    client = AlphaGenomeClient(cfg)

    active = {"n": 0, "max": 0}
    import threading
    lock = threading.Lock()

    def slow_score(chromosome, position, ref, alt, **kw):
        with lock:
            active["n"] += 1
            active["max"] = max(active["max"], active["n"])
        time.sleep(0.05)
        with lock:
            active["n"] -= 1
        return _make_tidy(f"{chromosome}:{position}:{ref}>{alt}")

    mock_score.side_effect = slow_score

    variants = [
        {"chromosome": "chr1", "position": i, "ref": "A", "alt": "G"}
        for i in range(6)
    ]
    client.score_variants(variants)

    assert active["max"] <= 2, f"semaphore allowed {active['max']} concurrent calls"


def test_client_config_max_concurrent_floor():
    """max_concurrent < 1 gets clamped to 1 in __post_init__."""
    cfg = ClientConfig(max_concurrent=0)
    assert cfg.max_concurrent == 1
    cfg2 = ClientConfig(max_concurrent=-5)
    assert cfg2.max_concurrent == 1


def test_score_variants_empty_returns_empty_df():
    cfg = ClientConfig(min_seconds_between_calls=0)
    cfg.api_key = "dummy"
    client = AlphaGenomeClient(cfg)
    df = client.score_variants([])
    assert df.empty
