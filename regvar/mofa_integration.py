"""Optional MOFA+ integration runner.

This module wraps ``mofapy2`` to train and compare factor models on the
views produced by ``regvar.mofa_view.build_mofa_view()``.  All ``mofapy2``
imports are lazy — you can safely import this module without having MOFA+
installed; only the train/compare functions will raise if it's missing.

Usage::

    from regvar.mofa_view import build_mofa_view
    from regvar.mofa_integration import train_mofa, compare_variant_contribution

    views = build_mofa_view(scores, genotypes)
    model = train_mofa(views, n_factors=10)
    comparison = compare_variant_contribution(base_views, variant_views)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger("regvar.mofa_integration")


# ---------------------------------------------------------------------------
# Internal — lazy mofapy2 helpers
# ---------------------------------------------------------------------------

def _get_mofa():
    """Lazy-import and return the mofapy2 entry point."""
    try:
        from mofapy2.run_mofa import MOFA
        return MOFA
    except ImportError as exc:
        raise ImportError(
            "MOFA+ integration requires mofapy2. Install with:\n"
            "  pip install mofapy2\n\n"
            "If you only need to format variant scores (not train models), "
            "use regvar.mofa_view.build_mofa_view() directly — it has no "
            "heavy dependencies."
        ) from exc


def _views_to_mofa_data(views: dict[str, pd.DataFrame]):
    """Convert pandas DataFrames to the nested list-of-arrays format
    mofapy2 expects for ``set_data_from_df()``.  We use a two-column
    (sample, feature, value) long-form approach per view, which avoids
    the need to transpose huge matrices."""
    import numpy as np

    data: dict[str, np.ndarray] = {}
    for view_name, df in views.items():
        data[view_name] = df.values.astype(np.float64)
    return data


def _validate_views(views: dict[str, pd.DataFrame]) -> None:
    """Check that views are compatible for joint MOFA+ training."""
    if len(views) < 1:
        raise ValueError("At least one view is required to train MOFA+.")

    sample_sets = [set(df.index) for df in views.values()]
    common = sample_sets[0]
    for i, s in enumerate(sample_sets[1:], 1):
        common = common & s

    if len(common) < 2:
        raise ValueError(
            f"MOFA+ requires at least 2 shared samples across all views. "
            f"Only {len(common)} found."
        )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# Standard likelihoods for common omic data types
LIKELIHOODS = {
    "gaussian":   "gaussian",   # continuous / normalised counts
    "poisson":    "poisson",    # raw counts
    "bernoulli":  "bernoulli",  # binary / methylation
}


def train_mofa(
    views: dict[str, pd.DataFrame],
    n_factors: int = 10,
    likelihoods: dict[str, str] | None = None,
    convergence_mode: str = "fast",
    seed: int = 42,
    output_path: str | Path | None = None,
    **kwargs,
):
    """Train a MOFA+ model on the given views.

    Parameters
    ----------
    views:
        View name → (N × D) DataFrame, as returned by
        ``build_mofa_view()``. The index must be sample IDs and must
        match across views.

    n_factors:
        Number of latent factors to learn.

    likelihoods:
        View-name → likelihood string mapping.  Defaults to
        ``"gaussian"`` for every view (appropriate for continuous
        variant-effect scores).

    convergence_mode:
        MOFA+ convergence mode: ``"fast"``, ``"medium"``, or ``"slow"``.

    seed:
        Random seed for reproducibility.

    output_path:
        If given, save the trained model as an ``.hdf5`` file.

    **kwargs:
        Passed through to ``MOFA.train_model()``.

    Returns
    -------
    MOFA model object
        The trained model. Use ``model.get_factors()``,
        ``model.get_weights()``, and ``model.get_variance_explained()``
        to inspect results.
    """
    import numpy as np

    _validate_views(views)

    # Align all views to shared sample set
    sample_sets = [set(df.index) for df in views.values()]
    common_samples = sorted(sample_sets[0].intersection(*sample_sets[1:]))
    if len(common_samples) < len(sample_sets[0]):
        logger.warning(
            "Views have different sample sets. Using %d shared samples "
            "(dropped %d sample(s) not present in all views).",
            len(common_samples),
            len(sample_sets[0]) - len(common_samples),
        )

    # Prepare data in mofapy2 format
    data_matrices: dict[str, np.ndarray] = {}
    for view_name, df in views.items():
        aligned = df.loc[common_samples]
        data_matrices[view_name] = aligned.values.astype(np.float64)

    if likelihoods is None:
        likelihoods = {vn: "gaussian" for vn in views}

    MOFA = _get_mofa()
    model = MOFA()
    model.set_data_matrix(data_matrices, likelihoods=likelihoods)
    model.set_model_options(
        factors=n_factors,
        convergence_mode=convergence_mode,
        seed=seed,
        **kwargs,
    )
    model.build()
    model.train()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        model.save(str(output_path))
        logger.info("Trained model saved to %s", output_path)

    return model


def compare_variant_contribution(
    base_views: dict[str, pd.DataFrame],
    variant_views: dict[str, pd.DataFrame],
    n_factors: int = 10,
    convergence_mode: str = "fast",
    seed: int = 42,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Train two MOFA+ models — with and without variant-effect views — and
    compare their factor decompositions.

    This is the key analysis function: it quantifies how much additional
    structure the regulatory variant scores contribute to the multi-omic
    factor model.

    Parameters
    ----------
    base_views:
        Views from real multi-omic data (e.g. RNA-seq, ATAC-seq, ChIP-seq
        matrices).  View name → (N × D) DataFrame.

    variant_views:
        Views built from variant-effect scores via ``build_mofa_view()``.
        Must have the same sample IDs as ``base_views``.

    n_factors:
        Number of factors for both models.

    convergence_mode, seed:
        Passed to both MOFA+ runs.

    output_path:
        If given, save a markdown comparison report to this path.

    Returns
    -------
    dict
        Keys: ``report`` (markdown string), ``model_with`` (trained MOFA
        model including variant views), ``model_without`` (trained MOFA
        model excluding variant views), ``variance_explained_with``,
        ``variance_explained_without`` (DataFrames).
    """
    import numpy as np

    # --- train without variant views -----------------------------------------
    logger.info("Training MOFA+ WITHOUT variant views (%d factors)...", n_factors)
    model_without = train_mofa(
        base_views,
        n_factors=n_factors,
        convergence_mode=convergence_mode,
        seed=seed,
    )

    # --- train with variant views --------------------------------------------
    combined_views = {**base_views, **variant_views}
    logger.info("Training MOFA+ WITH variant views (%d factors)...", n_factors)
    model_with = train_mofa(
        combined_views,
        n_factors=n_factors,
        convergence_mode=convergence_mode,
        seed=seed,
    )

    # --- extract variance explained ------------------------------------------
    def _var_explained_df(model) -> pd.DataFrame:
        """Extract variance explained per factor per view as a tidy DataFrame."""
        try:
            r2 = model.get_variance_explained()
            if r2 is None:
                return pd.DataFrame()
            rows = []
            for view_name, per_factor in r2.items():
                for factor_idx, val in enumerate(per_factor):
                    rows.append({
                        "view": view_name,
                        "factor": f"Factor {factor_idx + 1}",
                        "r2": float(val),
                    })
            return pd.DataFrame(rows)
        except Exception:
            return pd.DataFrame()

    var_with = _var_explained_df(model_with)
    var_without = _var_explained_df(model_without)

    # --- build comparison report ---------------------------------------------
    report_lines = [
        "# MOFA+ Variant-View Contribution Analysis",
        "",
        f"**Factors:** {n_factors}",
        f"**Convergence mode:** {convergence_mode}",
        f"**Seed:** {seed}",
        "",
        "## Views used",
        "",
        "### Base views (real multi-omic data)",
        *(f"- **{k}**: {v.shape[0]} samples × {v.shape[1]} features"
          for k, v in base_views.items()),
        "",
        "### Variant views (regulatory variant-effect scores)",
        *(f"- **{k}**: {v.shape[0]} samples × {v.shape[1]} features"
          for k, v in variant_views.items()),
        "",
        "## Variance Explained Comparison",
        "",
    ]

    if not var_with.empty and not var_without.empty:
        # Total variance explained per model
        total_with = var_with.groupby("view")["r2"].sum()
        total_without = var_without.groupby("view")["r2"].sum()

        report_lines.append("### Total variance explained by view")
        report_lines.append("")
        report_lines.append("| View | Without variants | With variants | Δ |")
        report_lines.append("|------|-----------------|---------------|----|")

        all_views = sorted(set(total_with.index) | set(total_without.index))
        for vn in all_views:
            r2_with = total_with.get(vn, 0.0)
            r2_without = total_without.get(vn, 0.0)
            delta = r2_with - r2_without
            sign = "+" if delta >= 0 else ""
            report_lines.append(
                f"| {vn} | {r2_without:.4f} | {r2_with:.4f} | {sign}{delta:.4f} |"
            )

        report_lines.append("")
        report_lines.append(
            "*A positive Δ indicates the variant view(s) helped explain "
            "additional variance in that view's features.*"
        )
    else:
        report_lines.append(
            "*(Variance explained could not be extracted from one or both models.)*"
        )

    report = "\n".join(report_lines) + "\n"

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(report, encoding="utf-8")
        logger.info("Comparison report saved to %s", output_path)

    return {
        "report": report,
        "model_with": model_with,
        "model_without": model_without,
        "variance_explained_with": var_with,
        "variance_explained_without": var_without,
    }
