"""
Auth routes: login/register placeholders (to be implemented).
"""
from flask import Blueprint

auth_bp = Blueprint("auth", __name__, url_prefix="/auth")


@auth_bp.get("/login")
def login():
    """Placeholder login route."""
    return "Login page placeholder"
