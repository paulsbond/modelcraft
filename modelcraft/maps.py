import gemmi
import numpy as np


def read_map(path: str) -> gemmi.Ccp4Map:
    ccp4 = gemmi.read_ccp4_map(path, setup=True)
    if any(abs(ccp4.header_float(w)) > 0 for w in (50, 51, 52)):
        print(f"Non-zero origin in {path} is not supported and will be ignored.")
    for word_number in (50, 51, 52):
        ccp4.set_header_float(word_number, 0.0)
    return ccp4


def sample_probes(grid, ncac, reference, probes) -> np.ndarray:
    """Standardised density (float32) at probe points defined about reference N/CA/C.

    The reference is least-squares superposed onto the residue's N, CA and C
    (gemmi.Position lists) and the probes moved with it; values are interpolated
    linearly.
    """
    transform = gemmi.superpose_positions(ncac, reference).transform
    rotation = np.array(transform.mat.tolist())
    points = probes @ rotation.T + np.array(transform.vec.tolist())
    values = grid.interpolate_position_array(points).astype(np.float32)
    if values.std() > 0:
        values = (values - values.mean()) / values.std()
    return values
