"""MOFA+ view builder — bridge AlphaGenome variant scores to multi-omic
factor analysis.

MOFA+ (Multi-Omics Factor Analysis) expects one matrix per "view", each
shaped ``(N_samples, D_features)``. AlphaGenome produces per-variant
predicted effect scores across assays. This module owns the convention
that bridges the two:

    1. Group per-variant scores by assay → one feature column per variant ×
       assay per view.
    2. Apply an additive dosage model (0/1/2 alt copies → 0/score/2×score)
       to produce per-sample values.
    3. Output DataFrames ready for ``mofapy2.run_mofa.MOFA()``.

The genotype bridge is the standard additive genetic model used across
eQTL, TWAS, and multi-omic variant-integration studies.

Dependencies: pandas, numpy (already in the project). No mofapy2 required.
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .alphagenome_client import ASSAY_TO_OUTPUT

logger = logging.getLogger("regvar.mofa_view")

# ---------------------------------------------------------------------------
# Assay label → MOFA+-safe view name
# ---------------------------------------------------------------------------

def _sanitise_view_name(assay_label: str) -> str:
    """Convert an assay label into a MOFA+-safe view name."""
    return assay_label.replace(" / ", "_").replace("/", "_").replace(" ", "_")


# ---------------------------------------------------------------------------
# Genotype I/O
# ---------------------------------------------------------------------------

def read_genotypes_tsv(path: str | Path) -> pd.DataFrame:
    """Read a tab-separated genotype file.

    Required format::

        sample_id    chr8:127401060:G>T   chr8:127472793:A>C   ...
        SAMP_01      0                    1
        SAMP_02      1                    2
        SAMP_03      NaN                  0
        ...

    Values must be 0, 1, 2, or missing (empty / NaN).

    Returns
    -------
    pd.DataFrame
        ``sample_id`` column set as the index, one column per variant.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "sample_id" not in df.columns:
        raise ValueError(
            f"{path}: missing required column 'sample_id'. "
            "Genotype file must have a sample_id column as the first column."
        )
    df = df.set_index("sample_id")
    # Coerce to float (0/1/2/NaN)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# ---------------------------------------------------------------------------
# Score I/O
# ---------------------------------------------------------------------------

def scored_variants_from_tsv(
    scores_tsv: str | Path,
) -> list[pd.DataFrame]:
    """Reconstruct per-variant DataFrames from a combined scores TSV.

    The TSV is expected to have been produced by ``regvar agent run
    --output-tsv`` (or ``regvar score --output-tsv``), which includes a
    ``query_variant`` column linking each row back to its originating
    variant.

    Parameters
    ----------
    scores_tsv:
        Path to a TSV file with a ``query_variant`` column.

    Returns
    -------
    list[pd.DataFrame]
        One DataFrame per unique variant, in the order they first appear.
    """
    path = Path(scores_tsv)
    df = pd.read_csv(path, sep="\t", dtype=str)

    if "query_variant" not in df.columns:
        raise ValueError(
            f"{path}: missing required column 'query_variant'. "
            "Scores file must have been produced by 'regvar agent run --output-tsv' "
            "or 'regvar score --output-tsv'."
        )

    # Convert numeric columns
    for col in ("raw_score", "quantile_score", "abs_effect", "position"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    variant_ids = df["query_variant"].unique()
    return [df[df["query_variant"] == vid].copy() for vid in variant_ids]


# ---------------------------------------------------------------------------
# Core: build the view matrices
# ---------------------------------------------------------------------------

def _variant_id_from_record(rec: dict | pd.Series) -> str:
    """Build ``chr:pos:ref>alt`` from a score row.  Prefers the
    ``query_variant`` column; falls back to ``variant_id``."""
    if "query_variant" in rec and pd.notna(rec["query_variant"]):
        return str(rec["query_variant"])
    if "variant_id" in rec and pd.notna(rec["variant_id"]):
        return str(rec["variant_id"])
    return "unknown"


def _feature_name(
    variant_id: str,
    gene_name: str | None,
    output_type: str | None,
) -> str:
    """Build a MOFA+ feature name: ``chr:pos:ref>alt_gene_output``."""
    gene = str(gene_name) if gene_name and pd.notna(gene_name) else "NA"
    out = str(output_type) if output_type and pd.notna(output_type) else "NA"
    return f"{variant_id}_{gene}_{out}"


def _aggregate_scores(
    variant_df: pd.DataFrame,
    assay_label: str,
    by: str,
) -> list[dict[str, Any]]:
    """Given a single variant's score DataFrame, return the feature row(s)
    for *one* assay after applying the aggregation strategy."""
    assay_rows = variant_df[variant_df["assay"] == assay_label]
    if assay_rows.empty:
        return []

    # Pick the best score column: prefer quantile_score if it has any
    # non-NaN values; otherwise fall back to raw_score.
    score_col: str | None = None
    for candidate in ("quantile_score", "raw_score"):
        if candidate in assay_rows.columns and assay_rows[candidate].notna().any():
            score_col = candidate
            break

    if score_col is None:
        return []

    # Drop rows where the chosen score is NaN
    valid = assay_rows.dropna(subset=[score_col])
    if valid.empty:
        return []

    variant_id = _variant_id_from_record(valid.iloc[0])

    if by == "all":
        features: list[dict[str, Any]] = []
        seen_names: dict[str, int] = {}
        for _, row in valid.iterrows():
            base_name = _feature_name(
                variant_id,
                row.get("gene_name"),
                row.get("output_type"),
            )
            # Deduplicate feature names by appending a counter suffix
            if base_name in seen_names:
                seen_names[base_name] += 1
                fname = f"{base_name}_{seen_names[base_name]}"
            else:
                seen_names[base_name] = 0
                fname = base_name
            features.append({
                "feature_name": fname,
                "score": float(row[score_col]),
            })
        return features

    # max_abs or top_gene: pick the single strongest row
    best_idx = valid[score_col].abs().idxmax()
    best = valid.loc[best_idx]

    if by == "top_gene":
        gene = str(best.get("gene_name")) if pd.notna(best.get("gene_name")) else "top"
    else:
        gene = "combined"

    return [{
        "feature_name": _feature_name(
            variant_id,
            gene if by == "top_gene" else "combined",
            best.get("output_type"),
        ),
        "score": float(best[score_col]),
    }]


def build_mofa_view(
    scored_variants: list[pd.DataFrame],
    genotypes: pd.DataFrame,
    assays: list[str] | None = None,
    by: str = "max_abs",
    output_path: str | Path | None = None,
) -> dict[str, pd.DataFrame]:
    """Convert per-variant AlphaGenome scores into MOFA+-ready view matrices.

    This is the primary entry point. It takes a list of per-variant score
    DataFrames (from ``AlphaGenomeClient.score_variant``) and a genotype
    table, and returns one ``(N_samples, D_features)`` DataFrame per assay
    view, with the additive dosage model applied.

    Parameters
    ----------
    scored_variants:
        One DataFrame per variant, each returned by
        ``AlphaGenomeClient.score_variant()``. The variants should be in
        the same order as the columns of ``genotypes`` (after
        ``sample_id``).

    genotypes:
        DataFrame with ``sample_id`` as the index and one column per
        variant.  Column names must match the variant identifiers in the
        score DataFrames (i.e. ``chr:pos:ref>alt``). Values: 0, 1, 2, or
        NaN for missing genotypes.

    assays:
        Subset of assay labels to include.  ``None`` means all assays
        present in the scored variants.  Valid labels are the keys of
        ``ASSAY_TO_OUTPUT`` (e.g. ``"ATAC-seq"``, ``"RNA-seq"``).

    by:
        Aggregation strategy when a variant has multiple score rows for a
        single assay:

        - ``"max_abs"`` (default): pick the single row with the highest
          ``|quantile_score|``; label the feature with gene="combined".
        - ``"top_gene"``: same as max_abs, but the feature name uses the
          actual gene name from the top-scoring track.
        - ``"all"``: expand every track into its own feature column.

    output_path:
        If given, also save the views as an AnnData ``.h5ad`` file
        (requires ``anndata`` to be installed).

    Returns
    -------
    dict[str, pd.DataFrame]
        View name → ``(N_samples × D_features)`` DataFrame.
        Index = sample_id (from the genotype table).
        Columns = feature names, parseable back to variant/gene/assay.
        Suitable for direct input to ``mofapy2.run_mofa.MOFA()``.

    Raises
    ------
    ValueError
        If a variant in the scores is missing from the genotype table,
        if fewer than 2 samples are present, or if no valid views can be
        built.
    """
    # --- validation ----------------------------------------------------------
    if by not in ("max_abs", "top_gene", "all"):
        raise ValueError(
            f"Unknown aggregation strategy {by!r}. "
            "Choose 'max_abs', 'top_gene', or 'all'."
        )

    if len(genotypes) < 2:
        raise ValueError(
            f"MOFA+ requires at least 2 samples. Only {len(genotypes)} found."
        )

    if not scored_variants:
        raise ValueError("No scored variants provided.")

    # --- resolve assay list --------------------------------------------------
    all_assay_labels = list(ASSAY_TO_OUTPUT.keys())
    if assays is not None:
        unknown = set(assays) - set(all_assay_labels)
        if unknown:
            raise ValueError(
                f"Unknown assay(s): {sorted(unknown)}. "
                f"Choose from: {all_assay_labels}"
            )
        assay_labels = assays
    else:
        # Detect which assays are actually present in the score data
        present: set[str] = set()
        for vdf in scored_variants:
            if "assay" in vdf.columns:
                present.update(vdf["assay"].dropna().unique())
        assay_labels = [a for a in all_assay_labels if a in present]

    if not assay_labels:
        raise ValueError(
            "No valid assays found in the scored variant data. "
            "Check that the score DataFrames have an 'assay' column."
        )

    # --- build variant → genotype column map --------------------------------  -
    genotype_variant_ids = set(genotypes.columns)

    # Collect all variant IDs from the scored DataFrames and validate coverage
    variant_scores_map: dict[str, pd.DataFrame] = {}
    missing_geno: list[str] = []
    for vdf in scored_variants:
        vid = _variant_id_from_record(vdf.iloc[0]) if len(vdf) > 0 else "unknown"
        if vid == "unknown":
            logger.warning("Skipping a variant DataFrame with no identifiable variant ID.")
            continue
        if vid not in genotype_variant_ids:
            missing_geno.append(vid)
        else:
            variant_scores_map[vid] = vdf

    if missing_geno:
        raise ValueError(
            f"{len(missing_geno)} variant(s) found in scores but missing "
            f"from the genotype table:\n  {', '.join(missing_geno[:10])}"
            f"{'...' if len(missing_geno) > 10 else ''}\n"
            "Add these variants to the genotype file, or remove them from "
            "the scores."
        )

    if not variant_scores_map:
        raise ValueError(
            "No variants could be matched between scores and genotype table. "
            "Check that genotype column names match variant identifiers "
            "(e.g. 'chr8:127401060:G>T')."
        )

    # Warn about genotype columns with no corresponding scores
    extra_geno = genotype_variant_ids - set(variant_scores_map.keys())
    if extra_geno:
        logger.warning(
            "%d variant(s) in genotype table have no score data and will be "
            "skipped: %s",
            len(extra_geno),
            ", ".join(sorted(extra_geno)[:5])
            + ("..." if len(extra_geno) > 5 else ""),
        )

    # --- build views ----------------------------------------------------------
    sample_ids = genotypes.index.tolist()

    # Determine all assays actually present
    present_assays: set[str] = set()
    for vdf in variant_scores_map.values():
        if "assay" in vdf.columns:
            present_assays.update(vdf["assay"].dropna().unique())

    working_assays = [a for a in assay_labels if a in present_assays]
    skipped = set(assay_labels) - set(working_assays)
    if skipped:
        logger.warning(
            "No score data for assay(s): %s. These views will be omitted.",
            ", ".join(sorted(skipped)),
        )

    if not working_assays:
        raise ValueError(
            "No assays have score data after filtering. "
            "Check that the scored variants include the requested assays."
        )

    views: dict[str, pd.DataFrame] = {}

    for assay_label in working_assays:
        view_name = _sanitise_view_name(assay_label)
        columns: dict[str, list[float]] = {}
        dropped_features: list[str] = []

        # Collect feature scores per variant for this assay
        for vid, vdf in variant_scores_map.items():
            genotype_col = genotypes[vid]
            features = _aggregate_scores(vdf, assay_label, by)

            for feat in features:
                fname = feat["feature_name"]
                score = feat["score"]
                # Apply additive dosage: 0→0, 1→score, 2→2×score, NaN→NaN
                values = []
                for dosage in genotype_col:
                    if pd.isna(dosage):
                        values.append(np.nan)
                    else:
                        values.append(float(dosage) * score)
                columns[fname] = values

        if not columns:
            logger.warning(
                "Assay %r has no valid features. View %r will be empty and is "
                "omitted.",
                assay_label, view_name,
            )
            continue

        view_df = pd.DataFrame(columns, index=sample_ids)
        view_df.index.name = "sample_id"

        # Drop features that are all-NaN
        all_nan = view_df.columns[view_df.isna().all(axis=0)]
        if len(all_nan) > 0:
            logger.warning(
                "Dropping %d feature(s) with all-NaN values in view %r: %s",
                len(all_nan), view_name,
                ", ".join(all_nan[:5])
                + ("..." if len(all_nan) > 5 else ""),
            )
            view_df = view_df.drop(columns=all_nan)

        if view_df.shape[1] == 0:
            logger.warning(
                "View %r has no features after dropping all-NaN columns. Omitted.",
                view_name,
            )
            continue

        views[view_name] = view_df

        # Count missing genotypes and warn if substantial
        n_missing = genotypes[sorted(variant_scores_map.keys())].isna().any(axis=1).sum()
        if n_missing > 0:
            logger.info(
                "%d/%d sample(s) have at least one missing genotype in view %r. "
                "MOFA+ will natively impute these values.",
                n_missing, len(sample_ids), view_name,
            )

    if not views:
        raise ValueError(
            "No views could be built. All assays produced empty or all-NaN views. "
            "Check that scores exist for the requested variants and assays."
        )

    # --- optional AnnData output ----------------------------------------------
    if output_path is not None:
        _save_views_as_anndata(views, output_path)

    return views


# ---------------------------------------------------------------------------
# AnnData I/O helper
# ---------------------------------------------------------------------------

def _save_views_as_anndata(
    views: dict[str, pd.DataFrame],
    output_path: str | Path,
) -> None:
    """Save view matrices as an AnnData .h5ad file.

    Uses the ``anndata`` package (must be installed separately).  Views
    are stored in ``adata.obsm`` so they retain their named identity;
    the ``.X`` slot holds a concatenated matrix for convenience.
    """
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError(
            "Saving to .h5ad requires anndata. Install with: pip install anndata"
        ) from exc

    sample_ids = list(next(iter(views.values())).index)
    adata = ad.AnnData(X=np.empty((len(sample_ids), 0)))
    adata.obs_names = [str(s) for s in sample_ids]

    for view_name, view_df in views.items():
        # Ensure sample order matches
        aligned = view_df.loc[sample_ids]
        adata.obsm[view_name] = aligned.values
        # Store feature names
        adata.uns[f"{view_name}_feature_names"] = list(aligned.columns)

    adata.uns["mofa_view_names"] = list(views.keys())
    adata.uns["mofa_view_shapes"] = {k: v.shape for k, v in views.items()}

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path)
    logger.info("Views saved to %s", output_path)


def load_views_from_anndata(
    path: str | Path,
) -> dict[str, pd.DataFrame]:
    """Load MOFA+ views previously saved with ``build_mofa_view(output_path=...)``.

    Parameters
    ----------
    path:
        Path to an ``.h5ad`` file saved by ``_save_views_as_anndata``.

    Returns
    -------
    dict[str, pd.DataFrame]
        View name → (N × D) DataFrame, reconstructed from ``obsm``.
    """
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError(
            "Loading from .h5ad requires anndata. Install with: pip install anndata"
        ) from exc

    adata = ad.read_h5ad(path)
    views: dict[str, pd.DataFrame] = {}

    view_names = adata.uns.get("mofa_view_names", list(adata.obsm.keys()))
    for vn in view_names:
        if vn not in adata.obsm:
            continue
        feature_names = adata.uns.get(
            f"{vn}_feature_names",
            [f"feat_{i}" for i in range(adata.obsm[vn].shape[1])],
        )
        views[vn] = pd.DataFrame(
            adata.obsm[vn],
            index=adata.obs_names,
            columns=feature_names,
        )
        views[vn].index.name = "sample_id"

    return views
