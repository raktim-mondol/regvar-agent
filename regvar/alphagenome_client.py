"""Hardened client wrapper around the AlphaGenome API.

AlphaGenome is a *hosted* model: you call a remote service, you don't run the
weights locally. The public API is free for non-commercial use but rate-limited
and shared, and DeepMind state it is suited to thousands of predictions rather
than very large batches. That reality dictates the design here:

    * on-disk caching keyed by the exact request, so reruns and shared work cost
      nothing and are reproducible;
    * exponential backoff with jitter on transient/rate-limit errors;
    * a simple client-side rate limiter so an agent looping over a VCF cannot
      stampede the service;
    * one place that owns coordinate conventions, so the rest of the codebase
      (and any LLM agent calling it) cannot get REF/ALT or 0-/1-based wrong.

The genomic logic lives behind this boundary on purpose: an agent decides *what*
to score, this module guarantees the calls are correct and economical.

Requires:  pip install alphagenome
API key:   https://deepmind.google.com/science/alphagenome  (set ALPHAGENOME_API_KEY)
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import pickle
import random
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from dataclasses import dataclass, field
from threading import Semaphore
from typing import Any, Iterable, Sequence

logger = logging.getLogger("regvar.alphagenome")

# AlphaGenome output types mapped to the assays the Clark/Garvan multi-omic
# dataset actually generates. This is the bridge between "model capability" and
# "our data" -- a regulatory variant can be scored on exactly the readouts the
# wet lab measures, which is what makes a prediction interpretable to them.
ASSAY_TO_OUTPUT = {
    "ATAC-seq": "ATAC",            # chromatin accessibility
    "DNase-seq": "DNASE",          # chromatin accessibility (alt assay)
    "RNA-seq": "RNA_SEQ",          # gene expression
    "ChIP-seq histone": "CHIP_HISTONE",  # H3K4me3/me1, H3K27ac, H3K27me3
    "ChIP-seq TF": "CHIP_TF",      # transcription-factor binding
    "Hi-C / pcHi-C": "CONTACT_MAPS",     # enhancer-promoter looping
    "splicing": "SPLICE_SITES",    # splice-site usage
}

# Prostate microenvironment ontology terms (UBERON / Cell Ontology).
PROSTATE_GLAND = "UBERON:0002367"
FIBROBLAST = "CL:0000057"


@dataclass
class ClientConfig:
    api_key: str | None = None
    cache_dir: Path = field(default_factory=lambda: Path.home() / ".cache" / "regvar")
    min_seconds_between_calls: float = 1.0   # client-side throttle
    max_retries: int = 5
    base_backoff: float = 2.0                # seconds; doubles each retry
    sequence_length: str = "SEQUENCE_LENGTH_1MB"
    max_concurrent: int = 1                  # >1 enables parallel score_variants()

    def __post_init__(self) -> None:
        self.api_key = self.api_key or os.environ.get("ALPHAGENOME_API_KEY")
        self.cache_dir = Path(self.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        if self.max_concurrent < 1:
            self.max_concurrent = 1


class AlphaGenomeClient:
    """Thin, defensive wrapper. Lazily imports `alphagenome` so the rest of the
    package (tests, tool schemas, MCP wiring) can be imported without the heavy
    dependency or a network connection present."""

    def __init__(self, config: ClientConfig | None = None) -> None:
        self.config = config or ClientConfig()
        self._model = None
        self._dna_client = None
        self._variant_scorers = None
        self._genome = None
        self._last_call_ts = 0.0
        self._throttle_lock = __import__("threading").Lock()

    # -- lazy backend ------------------------------------------------------
    def _ensure_model(self) -> None:
        if self._model is not None:
            return
        if not self.config.api_key:
            raise RuntimeError(
                "No AlphaGenome API key. Get one (free, non-commercial) at "
                "https://deepmind.google.com/science/alphagenome and set "
                "ALPHAGENOME_API_KEY."
            )
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers

        self._genome = genome
        self._dna_client = dna_client
        self._variant_scorers = variant_scorers
        self._model = dna_client.create(self.config.api_key)
        logger.info("AlphaGenome model handle created.")

    # -- caching -----------------------------------------------------------
    def _cache_path(self, key: dict[str, Any]) -> Path:
        blob = json.dumps(key, sort_keys=True).encode()
        digest = hashlib.sha256(blob).hexdigest()[:24]
        return self.config.cache_dir / f"{digest}.parquet"

    def _cache_path_legacy(self, key: dict[str, Any]) -> Path:
        blob = json.dumps(key, sort_keys=True).encode()
        digest = hashlib.sha256(blob).hexdigest()[:24]
        return self.config.cache_dir / f"{digest}.pkl"

    def _load_cache(self, key: dict[str, Any]):
        path = self._cache_path(key)
        if path.exists():
            logger.debug("cache hit (parquet) %s", path.name)
            import pandas as pd
            return pd.read_parquet(path)
        legacy = self._cache_path_legacy(key)
        if legacy.exists():
            logger.debug("cache hit (pickle) %s", legacy.name)
            with legacy.open("rb") as fh:
                return pickle.load(fh)
        return None

    def _save_cache(self, key: dict[str, Any], value: Any) -> None:
        path = self._cache_path(key)
        try:
            value.to_parquet(path)
        except Exception:
            with path.with_suffix(".pkl").open("wb") as fh:
                pickle.dump(value, fh)

    # -- throttle + retry --------------------------------------------------
    def _throttle(self) -> None:
        with self._throttle_lock:
            wait = self.config.min_seconds_between_calls - (time.time() - self._last_call_ts)
            if wait > 0:
                time.sleep(wait)
            self._last_call_ts = time.time()

    def _call_with_retry(self, fn, *args, **kwargs):
        for attempt in range(self.config.max_retries):
            self._throttle()
            try:
                return fn(*args, **kwargs)
            except Exception as exc:  # noqa: BLE001 - surface, classify, back off
                transient = any(
                    tok in str(exc).lower()
                    for tok in ("rate", "quota", "timeout", "unavailable", "503", "429")
                )
                if not transient or attempt == self.config.max_retries - 1:
                    raise
                sleep = self.config.base_backoff * (2 ** attempt) + random.uniform(0, 1)
                logger.warning("transient error (%s); retry %d in %.1fs",
                               exc, attempt + 1, sleep)
                time.sleep(sleep)

    # -- public API --------------------------------------------------------
    def score_variant(
        self,
        chromosome: str,
        position: int,
        ref: str,
        alt: str,
        assays: Sequence[str] = ("ATAC-seq", "RNA-seq", "ChIP-seq histone"),
        ontology_terms: Sequence[str] = (PROSTATE_GLAND,),
    ):
        """Score one variant across the requested assays and return a tidy
        pandas DataFrame (one row per variant x gene x track x scorer).

        `position` is 1-based (VCF convention). Coordinate handling is owned
        here so callers never have to think about it.
        """
        cache_key = {
            "op": "score_variant", "chrom": chromosome, "pos": position,
            "ref": ref, "alt": alt, "assays": sorted(assays),
            "ontology": sorted(ontology_terms),
        }
        cached = self._load_cache(cache_key)
        if cached is not None:
            return cached

        self._ensure_model()
        genome = self._genome
        vs = self._variant_scorers
        seq_len = getattr(self._dna_client, self.config.sequence_length)

        variant = genome.Variant(
            chromosome=chromosome,
            position=int(position),
            reference_bases=ref.upper(),
            alternate_bases=alt.upper(),
        )
        interval = variant.reference_interval.resize(seq_len)

        scorers = []
        for assay in assays:
            output_name = ASSAY_TO_OUTPUT.get(assay)
            if output_name is None:
                raise ValueError(f"Unknown assay {assay!r}; choose from {list(ASSAY_TO_OUTPUT)}")
            scorers.append(vs.RECOMMENDED_VARIANT_SCORERS[output_name])

        raw = self._call_with_retry(
            self._model.score_variant,
            interval=interval,
            variant=variant,
            variant_scorers=scorers,
        )
        tidy = vs.tidy_scores(raw, match_gene_strand=True)
        # Attach the human-readable assay label for downstream interpretation.
        out_to_assay = {v: k for k, v in ASSAY_TO_OUTPUT.items()}
        if "output_type" in tidy.columns:
            tidy["assay"] = tidy["output_type"].map(out_to_assay).fillna(tidy["output_type"])
        self._save_cache(cache_key, tidy)
        return tidy

    def score_variants(
        self,
        variants: Iterable[dict[str, Any]],
        max_concurrent: int | None = None,
        **kwargs,
    ):
        """Batch helper. `variants` is an iterable of dicts with keys
        chromosome, position, ref, alt. Returns a single concatenated DataFrame.
        Caching means re-running a partially-completed batch is cheap.

        Pass ``max_concurrent`` to override ``ClientConfig.max_concurrent`` for
        this call. When the value is ``1`` (default), calls run serially and
        behaviour is identical to the pre-parallel implementation. Higher
        values dispatch calls through a ``ThreadPoolExecutor`` bounded by a
        ``Semaphore``; the per-call throttle still acts as the floor.
        """
        import pandas as pd

        variant_list = list(variants)
        concurrency = max_concurrent if max_concurrent is not None else self.config.max_concurrent
        if concurrency < 1:
            concurrency = 1

        def _score_one(v: dict[str, Any]) -> "pd.DataFrame":
            df = self.score_variant(
                chromosome=v["chromosome"], position=v["position"],
                ref=v["ref"], alt=v["alt"], **kwargs,
            )
            df = df.copy()
            df["query_variant"] = f"{v['chromosome']}:{v['position']}:{v['ref']}>{v['alt']}"
            return df

        if concurrency == 1:
            frames = [_score_one(v) for v in variant_list]
        else:
            semaphore = Semaphore(concurrency)

            def _bounded(v: dict[str, Any]) -> "pd.DataFrame":
                with semaphore:
                    return _score_one(v)

            with ThreadPoolExecutor(max_workers=concurrency) as pool:
                frames = list(pool.map(_bounded, variant_list))

        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    @staticmethod
    def supported_assays() -> dict[str, str]:
        return dict(ASSAY_TO_OUTPUT)
