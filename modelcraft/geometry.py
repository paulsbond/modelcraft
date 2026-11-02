"Coordinate geometry on numpy arrays."

import numpy as np


def _unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def torsion_position(abc, length, angle, torsion):
    """Position d with |cd| = length, angle bcd and torsion abcd (radians).

    The positions in abc may be (3,) or batches (..., 3), and the angle and torsion
    may be arrays broadcasting with them.
    """
    a, b, c = (np.asarray(x, dtype=float) for x in abc)
    bc = _unit(c - b)
    n = _unit(np.cross(b - a, bc))
    m = np.cross(n, bc)
    angle = np.asarray(angle, dtype=float)[..., None]
    torsion = np.asarray(torsion, dtype=float)[..., None]
    return c + length * (
        -np.cos(angle) * bc
        + np.sin(angle) * np.cos(torsion) * m
        + np.sin(angle) * np.sin(torsion) * n
    )


def superpose(source, target):
    """Rotations (..., 3, 3) and translations (..., 3) superposing source points
    (..., n, 3) onto target points (..., n, 3) by least squares."""
    source, target = np.asarray(source, dtype=float), np.asarray(target, dtype=float)
    sc, tc = source.mean(axis=-2), target.mean(axis=-2)
    h = np.swapaxes(source - sc[..., None, :], -1, -2) @ (target - tc[..., None, :])
    u, _, vt = np.linalg.svd(h)
    v, ut = np.swapaxes(vt, -1, -2), np.swapaxes(u, -1, -2)
    v[..., 2] *= np.sign(np.linalg.det(v @ ut))[..., None]
    rotation = v @ ut
    return rotation, tc - (rotation @ sc[..., None])[..., 0]
