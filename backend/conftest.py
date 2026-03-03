"""
conftest.py — adds backend/ to sys.path so that
  `from services.X import Y` works in all unit tests,
  regardless of which directory pytest is invoked from.
"""
import sys
import os

# Insert the backend/ directory at the front of sys.path
_backend_dir = os.path.dirname(os.path.abspath(__file__))
if _backend_dir not in sys.path:
    sys.path.insert(0, _backend_dir)
