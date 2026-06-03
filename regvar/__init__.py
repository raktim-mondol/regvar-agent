"""regvar: AlphaGenome-powered regulatory-variant triage with LLM-agent orchestration."""

# Load .env from the project root (one level up from this package) if present.
# Runs once on package import so agent.py and alphagenome_client.py don't need to.
try:
    from dotenv import load_dotenv
    from pathlib import Path
    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
except ImportError:
    pass

from .alphagenome_client import AlphaGenomeClient, ClientConfig, ASSAY_TO_OUTPUT
from .config import get_config, get_default
from .variants import CandidateVariant, read_candidates_tsv, read_candidates_vcf, rank_by_effect, validate_chromosome
from .annotation import VariantAnnotation, annotate_variants, get_annotation
from .mofa_view import (
    build_mofa_view,
    read_genotypes_tsv,
    scored_variants_from_tsv,
    load_views_from_anndata,
)

__all__ = [
    "AlphaGenomeClient", "ClientConfig", "ASSAY_TO_OUTPUT",
    "CandidateVariant", "read_candidates_tsv", "read_candidates_vcf",
    "rank_by_effect", "validate_chromosome",
    "get_config", "get_default",
    "VariantAnnotation", "annotate_variants", "get_annotation",
    # MOFA+ integration
    "build_mofa_view", "read_genotypes_tsv", "scored_variants_from_tsv",
    "load_views_from_anndata",
]
__version__ = "0.1.0"
