"Shared pieces of the ML networks in mlgrow and mlside."

import numpy as np

from .geometry import superpose
from .stdres import STD_C, STD_CA, STD_N

# N, CA and C defining the frame of the ML networks
REFERENCE = np.array([STD_N, STD_CA, STD_C])


def place_reference(residues: np.ndarray) -> np.ndarray:
    "REFERENCE superposed onto each of the residues (b, 3, 3)"
    rotations, translations = superpose(REFERENCE, residues)
    return np.matmul(REFERENCE, rotations.transpose(0, 2, 1)) + translations[:, None]


def sample_probes(grid, residues: np.ndarray, probes) -> np.ndarray:
    """Standardised density (b, probes) float32 at probe points defined about
    REFERENCE, superposed onto each of the residues (b, 3, 3), interpolated
    linearly."""
    rotations, translations = superpose(REFERENCE, residues)
    points = np.matmul(probes, rotations.transpose(0, 2, 1)) + translations[:, None]
    values = grid.interpolate_position_array(points.reshape(-1, 3))
    values = values.reshape(len(residues), -1).astype(np.float32)
    std = values.std(axis=1, keepdims=True)
    return (values - values.mean(axis=1, keepdims=True)) / np.where(std > 0, std, 1)
