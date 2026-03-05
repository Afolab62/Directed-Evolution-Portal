"""
conftest.py — pytest session-wide configuration for the backend test suite.

Ensures the ``backend/`` package root is on sys.path so all ``services.*``
and ``routes.*`` imports resolve correctly regardless of how pytest is invoked.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Insert the backend/ directory at the start of sys.path.
# pytest.ini already declares ``pythonpath = .`` (requires pytest ≥ 7),
# this guard works for any pytest version.
_backend_root = Path(__file__).parent.resolve()
if str(_backend_root) not in sys.path:
    sys.path.insert(0, str(_backend_root))
