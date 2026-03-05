"""
landscape_service.py
--------------------
Activity landscape computation matching the logic in 3d_landscape.py.

Key design decisions (mirroring the original):
  - Amino-acid trigram CountVectorizer encoding captures local sequence
    patterns (e.g. active-site motifs) without requiring structural data.
  - TruncatedSVD operates directly on the sparse trigram matrix avoiding
    the memory cost of converting to a dense array for standard PCA.
  - t-SNE embeddings are cached to disk keyed by a SHA-1 of the sequences
    and parameters; recomputation on every request would be unusable.
  - One surface frame is generated per generation using cumulative data
    (gen <= x), so the slider shows the landscape growing over rounds.
  - Three z-mode transforms share a single global colour scale so the
    animation is directly comparable across generations.
  - Returns a complete go.Figure (via to_json()) so the frontend only
    needs to pass it straight to react-plotly.js, identical to the
    fingerprint endpoint pattern.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.manifold import TSNE

# ── embedding cache directory ──────────────────────────────────────────────────
_CACHE_ROOT = Path(__file__).resolve().parents[1] / "instance" / "embedding_cache"
_CACHE_ROOT.mkdir(parents=True, exist_ok=True)

_K = 3          # trigram length
# Grid size and initial z-mode respect the same env vars as the original 3d_landscape.py
_GRID_SIZE = int(os.getenv("ACTIVITY_GRID_SIZE", "120"))  # 120×120 matches original
_INITIAL_Z_MODE = os.getenv("ACTIVITY_Z_MODE", "robust").strip().lower()
if _INITIAL_Z_MODE not in {"raw", "normalized", "robust"}:
    _INITIAL_Z_MODE = "robust"


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------

def _trigram_encode(sequences: list[str]):
    """Encode protein sequences as character-trigram count vectors (sparse)."""
    vec = CountVectorizer(analyzer="char", ngram_range=(_K, _K), lowercase=False)
    return vec.fit_transform(sequences)


# ---------------------------------------------------------------------------
# Dimensionality reduction with disk cache
# ---------------------------------------------------------------------------

def _cache_path(sequences: list[str], method: str, **params) -> Path:
    sig = hashlib.sha1()
    sig.update(str(_K).encode())
    sig.update(method.encode())
    for s in sequences:
        sig.update(s.encode("utf-8", errors="ignore"))
        sig.update(b"|")
    for k, v in sorted(params.items()):
        sig.update(f"{k}:{v}".encode())
    return _CACHE_ROOT / f"{method}_{sig.hexdigest()}.npy"


def _reduce_to_2d(sequences: list[str], method: str) -> np.ndarray:
    """
    Reduce sequences to 2-D coordinates via trigram encoding + chosen method.
    Cached on disk to avoid recomputing expensive t-SNE runs.
    """
    X = _trigram_encode(sequences)

    if method == "tsne":
        perplexity = float(os.getenv("ACTIVITY_TSNE_PERPLEXITY", "15"))
        max_iter   = int(os.getenv("ACTIVITY_TSNE_MAX_ITER", "500"))
        n_samples  = X.shape[0]
        # sklearn requires perplexity < n_samples; clamp so small datasets work.
        perplexity = min(perplexity, max(1.0, n_samples - 1))
        cache = _cache_path(sequences, "tsne", perplexity=perplexity, max_iter=max_iter)
        if cache.exists():
            return np.load(cache)
        try:
            XY = TSNE(n_components=2, perplexity=perplexity,
                      random_state=0, init="pca",
                      learning_rate="auto", max_iter=max_iter).fit_transform(X)
        except TypeError:
            # sklearn < 1.4 used n_iter
            XY = TSNE(n_components=2, perplexity=perplexity,
                      random_state=0, init="pca",
                      learning_rate="auto", n_iter=max_iter).fit_transform(X)
        np.save(cache, XY)
        return XY

    if method == "umap":
        try:
            import umap as _umap
            cache = _cache_path(sequences, "umap")
            if cache.exists():
                return np.load(cache)
            XY = _umap.UMAP(n_components=2, random_state=0).fit_transform(X.toarray())
            np.save(cache, XY)
            return XY
        except ImportError:
            pass  # fall through to PCA

    # PCA — TruncatedSVD works directly on sparse matrix
    cache = _cache_path(sequences, "pca")
    if cache.exists():
        return np.load(cache)
    XY = TruncatedSVD(n_components=2, random_state=0).fit_transform(X)
    np.save(cache, XY)
    return XY


# ---------------------------------------------------------------------------
# Surface interpolation
# ---------------------------------------------------------------------------

def _surface(df: pd.DataFrame, x_col: str, y_col: str, z_col: str,
             grid_x: np.ndarray, grid_y: np.ndarray) -> np.ndarray:
    """
    Linear interpolation onto regular grid, NaN edges filled with nearest-
    neighbour, then Gaussian smoothing (sigma=1.5).
    """
    points = np.column_stack([df[x_col].to_numpy(), df[y_col].to_numpy()])
    values = df[z_col].to_numpy()
    z = griddata(points, values, (grid_x, grid_y), method="linear")
    if np.isnan(z).any():
        z_near = griddata(points, values, (grid_x, grid_y), method="nearest")
        z = np.where(np.isnan(z), z_near, z)
    return gaussian_filter(z, sigma=1.5)


# ---------------------------------------------------------------------------
# Z-mode transforms (shared global colour scale)
# ---------------------------------------------------------------------------

def _apply_mode(values: np.ndarray, mode: str,
                raw_min: float, raw_max: float,
                robust_lo: float, robust_hi: float) -> np.ndarray:
    if mode == "raw":
        return values
    if mode == "robust":
        return np.clip(values, robust_lo, robust_hi)
    out = (values - raw_min) / (raw_max - raw_min + 1e-9)
    return np.clip(out, 0, 1) ** 1.3


def _mode_range(mode: str, raw_min: float, raw_max: float,
                robust_lo: float, robust_hi: float) -> tuple[float, float]:
    if mode == "raw":
        return raw_min, raw_max
    if mode == "robust":
        return robust_lo, robust_hi
    return 0.0, 1.0


def _mode_axis_title(mode: str) -> str:
    titles = {
        "raw":        "Activity Score",
        "normalized": "Activity Score (Normalized)",
        "robust":     "Activity Score (Robust-clipped)",
    }
    return titles.get(mode, "Activity Score")


def _fname(mode: str, generation: Any) -> str:
    return f"{mode}:{generation}"


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def build_landscape_figure(
    sequences: list[str],
    activity_scores: list[float],
    generations: list[int],
    variant_indices: list[Any],
    method: str = "pca",
) -> go.Figure:
    """
    Build a complete animated Plotly Figure matching 3d_landscape.py output.

    Returns a go.Figure with per-generation animation frames, three z-mode
    transforms, play buttons, and a generation slider.
    """
    if len(sequences) < 3:
        raise ValueError("Need at least 3 variants to compute the landscape.")

    XY = _reduce_to_2d(sequences, method)

    data = pd.DataFrame({
        "x":              XY[:, 0],
        "y":              XY[:, 1],
        "activity_score": activity_scores,
        "generation":     pd.to_numeric(generations, errors="coerce"),
        "variant_index":  variant_indices,
    }).dropna(subset=["activity_score", "generation"])

    g_x = np.linspace(data["x"].min(), data["x"].max(), _GRID_SIZE)
    g_y = np.linspace(data["y"].min(), data["y"].max(), _GRID_SIZE)
    grid_x, grid_y = np.meshgrid(g_x, g_y)

    # global colour bounds shared by all frames
    all_act  = data["activity_score"].to_numpy(dtype=float)
    raw_min  = float(np.nanmin(all_act))
    raw_max  = float(np.nanmax(all_act))
    robust_lo = float(np.nanpercentile(all_act, 1))
    robust_hi = float(np.nanpercentile(all_act, 99))

    gens    = sorted(data["generation"].unique())
    z_modes = ["robust", "raw", "normalized"]
    frames: list[go.Frame] = []

    for mode in z_modes:
        m_min, m_max = _mode_range(mode, raw_min, raw_max, robust_lo, robust_hi)
        for gen in gens:
            df_f = data[data["generation"] <= gen]

            z_surf = _apply_mode(
                _surface(df_f, "x", "y", "activity_score", grid_x, grid_y),
                mode, raw_min, raw_max, robust_lo, robust_hi,
            )
            z_scatter = _apply_mode(
                df_f["activity_score"].to_numpy(dtype=float),
                mode, raw_min, raw_max, robust_lo, robust_hi,
            )

            frames.append(go.Frame(
                name=_fname(mode, gen),
                data=[
                    go.Surface(
                        x=g_x.tolist(), y=g_y.tolist(), z=z_surf.tolist(),
                        opacity=0.65, colorscale="Hot", showscale=True,
                        cmin=m_min, cmax=m_max,
                        colorbar=dict(title=_mode_axis_title(mode)),
                    ),
                    go.Scatter3d(
                        x=df_f["x"].tolist(),
                        y=df_f["y"].tolist(),
                        z=z_scatter.tolist(),
                        customdata=np.column_stack([
                            df_f["variant_index"].to_numpy(),
                            df_f["generation"].to_numpy(),
                            df_f["activity_score"].to_numpy(),
                        ]),
                        hovertemplate=(
                            "Variant: %{customdata[0]}<br>"
                            "Generation: %{customdata[1]}<br>"
                            "Activity Score: %{customdata[2]:.3f}<br>"
                            "<extra></extra>"
                        ),
                        mode="markers",
                        marker=dict(size=4, opacity=0.5, color="mediumpurple"),
                    ),
                ],
                layout=go.Layout(scene=dict(
                    zaxis_title=_mode_axis_title(mode),
                    zaxis=dict(range=[m_min, m_max]),
                )),
            ))

    initial_frame = next(f for f in frames if f.name == _fname(_INITIAL_Z_MODE, gens[0]))

    # slider steps: all (mode × generation) combinations
    slider_steps: list[dict] = []
    for label, key in [("Robust", "robust"), ("Raw", "raw"), ("Norm", "normalized")]:
        for gen in gens:
            g_label = int(gen) if float(gen).is_integer() else gen
            slider_steps.append(dict(
                label=f"{label} G{g_label}",
                method="animate",
                args=[[_fname(key, gen)],
                      {"mode": "immediate", "frame": {"duration": 0, "redraw": True}}],
            ))

    robust_fnames = [_fname("robust", g) for g in gens]
    raw_fnames    = [_fname("raw", g)    for g in gens]
    norm_fnames   = [_fname("normalized", g) for g in gens]

    axis_label = {"tsne": "t-SNE", "umap": "UMAP"}.get(method, "PC")

    fig = go.Figure(data=initial_frame.data, frames=frames)
    fig.update_layout(
        title="3D Activity Landscape",
        margin=dict(l=0, r=20, b=0, t=45),
        scene=dict(
            xaxis_title=f"{axis_label} Dim 1",
            yaxis_title=f"{axis_label} Dim 2",
            zaxis_title=_mode_axis_title(_INITIAL_Z_MODE),
            zaxis=dict(range=list(_mode_range(
                _INITIAL_Z_MODE, raw_min, raw_max, robust_lo, robust_hi))),
        ),
        scene_camera=dict(eye=dict(x=-1.8, y=-1.8, z=1.2), up=dict(x=0, y=0, z=1)),
        scene_aspectmode="manual",
        scene_aspectratio=dict(x=1, y=1, z=0.7),
        uirevision="keep",
        height=660,
        updatemenus=[dict(
            type="buttons", showactive=False,
            x=0.0, xanchor="left", y=0.02, yanchor="bottom",
            buttons=[
                dict(label="▶ Robust", method="animate",
                     args=[robust_fnames, {"frame": {"duration": 500, "redraw": True},
                                           "transition": {"duration": 200}}]),
                dict(label="▶ Raw", method="animate",
                     args=[raw_fnames, {"frame": {"duration": 500, "redraw": True},
                                        "transition": {"duration": 200}}]),
                dict(label="▶ Norm", method="animate",
                     args=[norm_fnames, {"frame": {"duration": 500, "redraw": True},
                                         "transition": {"duration": 200}}]),
                dict(label="Robust", method="animate",
                     args=[[robust_fnames[0]], {"mode": "immediate",
                                                 "frame": {"duration": 0, "redraw": True}}]),
                dict(label="Raw", method="animate",
                     args=[[raw_fnames[0]], {"mode": "immediate",
                                              "frame": {"duration": 0, "redraw": True}}]),
                dict(label="Norm", method="animate",
                     args=[[norm_fnames[0]], {"mode": "immediate",
                                               "frame": {"duration": 0, "redraw": True}}]),
                dict(label="■ Stop", method="animate",
                     args=[[None], {"frame": {"duration": 0, "redraw": False},
                                    "mode": "immediate"}]),
            ],
        )],
        sliders=[dict(
            currentvalue=dict(prefix="Frame: "),
            steps=slider_steps,
        )],
    )

    return fig


# (old one_hot_encode / get_esm_embeddings / compute_fitness_landscape removed;
#  use build_landscape_figure() instead)

