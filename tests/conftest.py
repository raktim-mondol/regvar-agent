"""Test configuration — force non-interactive matplotlib backend."""

import os

os.environ.setdefault("MPLBACKEND", "Agg")
