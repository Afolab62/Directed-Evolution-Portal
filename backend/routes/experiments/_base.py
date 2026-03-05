"""
Shared Blueprint instance and tiny helpers used by every sub-module.
"""
from flask import Blueprint, session
import math

experiments_bp = Blueprint('experiments', __name__, url_prefix='/api/experiments')


def require_auth():
    """Return user_id if authenticated, else None."""
    return session.get('user_id') or None


def clean_dict_for_json(obj):
    """Recursively convert NaN/Inf float values to None for JSON serialisation."""
    if isinstance(obj, dict):
        return {k: clean_dict_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_dict_for_json(item) for item in obj]
    elif isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj
