"Shared pieces of the ML networks in mlgrow and mlside."

import numpy as np

# N, CA and C defining the frame of the ML networks (CA at the origin)
REFERENCE = np.array(
    [
        [+0.8431, +0.0000, -1.2040],
        [+0.0000, +0.0000, +0.0000],
        [+0.8775, +0.0000, +1.2530],
    ]
)


def reference_transforms(residues: np.ndarray):
    """Rotations (b, 3, 3) and translations (b, 3) taking REFERENCE onto residues.

    The superpositions are least-squares onto N, CA and C positions (b, 3, 3).
    """
    centre = REFERENCE.mean(axis=0)
    centres = residues.mean(axis=1)
    h = np.einsum("pi,bpj->bij", REFERENCE - centre, residues - centres[:, None])
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(np.einsum("bji,bkj->bik", vt, u)))
    vt[:, 2] *= d[:, None]
    rotations = np.einsum("bji,bkj->bik", vt, u)
    return rotations, centres - np.einsum("bij,j->bi", rotations, centre)


def place_reference(residues: np.ndarray) -> np.ndarray:
    "REFERENCE superposed onto each of the residues (b, 3, 3)"
    rotations, translations = reference_transforms(residues)
    return np.matmul(REFERENCE, rotations.transpose(0, 2, 1)) + translations[:, None]


def sample_probes(grid, residues: np.ndarray, probes) -> np.ndarray:
    """Standardised density (b, probes) float32 at probe points defined about
    REFERENCE, moved onto each of the residues (b, 3, 3) by reference_transforms;
    values are interpolated linearly."""
    rotations, translations = reference_transforms(residues)
    points = np.matmul(probes, rotations.transpose(0, 2, 1)) + translations[:, None]
    values = grid.interpolate_position_array(points.reshape(-1, 3))
    values = values.reshape(len(residues), -1).astype(np.float32)
    std = values.std(axis=1, keepdims=True)
    return (values - values.mean(axis=1, keepdims=True)) / np.where(std > 0, std, 1)
