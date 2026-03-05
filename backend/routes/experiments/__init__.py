"""
experiments package — all route sub-modules are imported here so their
@experiments_bp.route decorators fire and register with the shared Blueprint.

Sub-module layout:
  _base.py       shared Blueprint instance + tiny helpers
  core.py        CRUD  (POST/GET/PATCH/DELETE /api/experiments/...)
  upload.py      preview-mapping + upload-data
  variants.py    GET variants + top-performers
  fingerprint.py fingerprint / fingerprint3d / fingerprint_linear
  analysis.py    analyze-sequences + background worker
  export.py      mutations/export + plots/activity-distribution
"""

from ._base import experiments_bp  # noqa: F401 — must come first

# Import sub-modules to trigger route registration (order doesn't matter)
from . import core        # noqa: F401
from . import upload      # noqa: F401
from . import variants    # noqa: F401
from . import fingerprint # noqa: F401
from . import analysis    # noqa: F401
from . import export      # noqa: F401

__all__ = ['experiments_bp']
