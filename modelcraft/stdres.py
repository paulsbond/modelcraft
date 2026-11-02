"Standard residue frame."

import numpy as np

STD_N = np.array([0.8431, 0.0, -1.204])
STD_CA = np.array([0.0, 0.0, 0.0])
STD_C = np.array([0.8775, 0.0, 1.253])
STD_CB = np.array([-1.03, -1.11, 0.0])


def _unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def frame(n, ca, c):
    "Rotation and translation (numpy) taking the standard frame onto a residue"
    n, ca, c = (np.asarray(p, dtype=float) for p in (n, ca, c))
    c1, c2 = c - ca, n - ca
    v1 = _unit(_unit(c1) + _unit(c2))
    v2 = _unit(np.cross(c1, c2))
    v3 = _unit(np.cross(v1, v2))
    return np.stack([v1, v2, v3], axis=-1), ca
