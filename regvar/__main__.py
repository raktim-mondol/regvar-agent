"""Allow ``python -m regvar`` to invoke the CLI.

This replaces the old ``python -m regvar.agent <file>`` entry point with the
richer click-based CLI, while keeping backward-compatible behaviour:

    python -m regvar run examples/candidate_variants.tsv   # new
    python -m regvar assays                                 # new
    python -m regvar score chr8 127401060 G T               # new

The old entry point (``python -m regvar.agent <file>``) still works unchanged.
"""

from regvar.cli import cli

cli()
