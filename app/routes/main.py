"""
Main routes: lightweight landing page.
"""
from flask import Blueprint, render_template

main_bp = Blueprint("main", __name__)


@main_bp.get("/")
def home():
    """Render homepage."""
    return render_template("home.html")
