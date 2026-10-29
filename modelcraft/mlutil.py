"Shared pieces of the ML networks in mlgrow and mltype."

import gemmi
import numpy as np

# N, CA and C defining the frame of the ML networks (CA at the origin)
REFERENCE = [
    gemmi.Position(+0.8431, +0.0000, -1.2040),
    gemmi.Position(+0.0000, +0.0000, +0.0000),
    gemmi.Position(+0.8775, +0.0000, +1.2530),
]


def reference_transform(ncac):
    """Rotation and translation (numpy) taking REFERENCE onto N, CA and C.

    The superposition is least-squares, with positions as a gemmi.Position list.
    """
    transform = gemmi.superpose_positions(ncac, REFERENCE).transform
    return np.array(transform.mat.tolist()), np.array(transform.vec.tolist())


def sample_probes(grid, ncac, probes) -> np.ndarray:
    """Standardised density (float32) at probe points defined about REFERENCE.

    The probes are moved by reference_transform onto the residue's N, CA and C
    (gemmi.Position list); values are interpolated linearly.
    """
    rotation, translation = reference_transform(ncac)
    points = probes @ rotation.T + translation
    values = grid.interpolate_position_array(points).astype(np.float32)
    if values.std() > 0:
        values = (values - values.mean()) / values.std()
    return values
