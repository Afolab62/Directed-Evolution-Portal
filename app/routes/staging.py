from flask import Blueprint, render_template, request

from app.services.staging import stage_experiment_validate_plasmid

staging_bp = Blueprint("staging", __name__, url_prefix="/staging")


@staging_bp.get("/")
def staging_form():
    return render_template("staging.html", result=None, accession="", plasmid_fasta="", fetch_features=True)


@staging_bp.post("/")
def staging_submit():
    accession = request.form.get("accession", "").strip()
    fetch_features = bool(request.form.get("fetch_features"))

    plasmid_fasta = request.form.get("plasmid_fasta", "").strip()
    uploaded = request.files.get("plasmid_file")
    if uploaded and uploaded.filename:
        plasmid_fasta = uploaded.read().decode("utf-8")

    if not accession or not plasmid_fasta:
        result = {"error": "Please provide a UniProt accession and plasmid FASTA (paste or upload)."}
        return render_template(
            "staging.html",
            result=result,
            accession=accession,
            plasmid_fasta=plasmid_fasta,
            fetch_features=fetch_features,
        )

    result = stage_experiment_validate_plasmid(
        accession=accession,
        plasmid_fasta_text=plasmid_fasta,
        fetch_features=fetch_features,
    )

    return render_template(
        "staging.html",
        result=result,
        accession=accession,
        plasmid_fasta=plasmid_fasta,
        fetch_features=fetch_features,
    )
