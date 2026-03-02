import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def one_hot_encode(sequence: str) -> np.ndarray:
    """Encode a protein sequence as a one-hot vector."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
    encoding = np.zeros((len(sequence), len(amino_acids)))
    for i, aa in enumerate(sequence):
        if aa in aa_to_idx:
            encoding[i, aa_to_idx[aa]] = 1
    return encoding.flatten()


def get_esm_embeddings(sequences: list) -> np.ndarray:
    """
    Get ESM-2 embeddings for protein sequences.
    Falls back to one-hot encoding if ESM is not available.
    """
    try:
        import esm
        import torch

        model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()  # smallest ESM-2 model
        model.eval()
        batch_converter = alphabet.get_batch_converter()

        data = [(f"seq_{i}", seq) for i, seq in enumerate(sequences)]
        _, _, batch_tokens = batch_converter(data)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[6])

        # Mean pool over sequence length to get fixed-size embedding
        embeddings = results["representations"][6].mean(dim=1).numpy()
        return embeddings

    except ImportError:
        # Fallback: one-hot encoding padded to the same length
        print("ESM not available, falling back to one-hot encoding")
        max_len = max(len(s) for s in sequences)
        padded = [s.ljust(max_len, 'X') for s in sequences]
        return np.array([one_hot_encode(seq) for seq in padded])


def compute_fitness_landscape(
    sequences: list,
    activity_scores: list,
    method: str = "pca",   # "pca" | "tsne" | "umap"
    grid_resolution: int = 50
) -> dict:
    """
    Compute a 3D fitness landscape from protein sequences and activity scores.
    Returns grid data suitable for Plotly go.Surface.
    """
    if len(sequences) < 3:
        raise ValueError("Need at least 3 variants to compute landscape")

    # Step 1: Encode sequences
    embeddings = get_esm_embeddings(sequences)

    # Step 2: Dimensionality reduction to 2D
    if method == "tsne":
        perplexity = min(30, len(sequences) - 1)
        reducer = TSNE(n_components=2, perplexity=perplexity, random_state=42)
        coords_2d = reducer.fit_transform(embeddings)

    elif method == "umap":
        try:
            import umap
            reducer = umap.UMAP(n_components=2, random_state=42)
            coords_2d = reducer.fit_transform(embeddings)
        except ImportError:
            print("UMAP not available, falling back to PCA")
            reducer = PCA(n_components=2)
            coords_2d = reducer.fit_transform(embeddings)

    else:  # default: PCA
        reducer = PCA(n_components=2)
        coords_2d = reducer.fit_transform(embeddings)

    x_coords = coords_2d[:, 0]
    y_coords = coords_2d[:, 1]
    z_coords = np.array(activity_scores)

    # Step 3: Interpolate onto a regular grid for surface rendering
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()

    grid_x, grid_y = np.mgrid[
        x_min:x_max:complex(grid_resolution),
        y_min:y_max:complex(grid_resolution)
    ]

    grid_z = griddata(
        points=(x_coords, y_coords),
        values=z_coords,
        xi=(grid_x, grid_y),
        method="linear"
    )

    # Fill NaN edges with nearest-neighbour interpolation
    grid_z_filled = griddata(
        points=(x_coords, y_coords),
        values=z_coords,
        xi=(grid_x, grid_y),
        method="nearest"
    )
    nan_mask = np.isnan(grid_z)
    grid_z[nan_mask] = grid_z_filled[nan_mask]

    # Gaussian smoothing to produce a smooth rolling landscape surface
    grid_z = gaussian_filter(grid_z, sigma=1.5)

    return {
        "x": grid_x.tolist(),
        "y": grid_y.tolist(),
        "z": grid_z.tolist(),
        "scatter_points": {
            "x": x_coords.tolist(),
            "y": y_coords.tolist(),
            "z": z_coords.tolist(),
        },
        "method": method,
        "variant_count": len(sequences),
    }
