"The next or previous residue's N, CA and C from the density around a residue."

import functools
import math
from pathlib import Path

import gemmi
import numpy as np

from .geometry import pos4
from .maps import sample_probes

_DATA = Path(__file__).parent / "data"

_REFERENCE = [
    gemmi.Position(+0.8431, +0.0000, -1.2040),  # N
    gemmi.Position(+0.0000, +0.0000, +0.0000),  # CA
    gemmi.Position(+0.8775, +0.0000, +1.2530),  # C
]


@functools.cache
def _network(direction: str):
    "Probes and layers (float32) of one direction; loaded on first use."
    w = np.load(_DATA / f"mlgrow-{direction}.npz")

    def array(key):
        return w[key].astype(np.float32)

    probes = w["probes"].astype(np.float64)
    layers = [
        (array(f"{layer}_kernel"), array(f"{layer}_bias"))
        for layer in ("dense", "dense_1", "dense_2", "dense_3")
    ]
    heads = tuple(
        array(key) for key in ("sin_kernel", "sin_bias", "cos_kernel", "cos_bias")
    )
    return probes, layers, heads


def _predict(grid: gemmi.FloatGrid, residue: gemmi.Residue, direction: str):
    "N, CA, C as numpy arrays, (psi, omega, phi) and the confidence, else None."
    atoms = [residue.find_atom(name, "*") for name in ("N", "CA", "C")]
    if any(atom is None for atom in atoms):
        return None
    probes, layers, (sk, sb, ck, cb) = _network(direction)
    x = sample_probes(grid, [a.pos for a in atoms], _REFERENCE, probes)
    for w, b in layers:
        x = np.maximum(x @ w + b, 0)
    sin, cos = x @ sk + sb, x @ ck + cb
    angles = np.arctan2(sin, cos).astype(np.float64)
    ncac = [np.array(a.pos.tolist()) for a in atoms]
    return ncac, angles, math.hypot(sin[0], cos[0])


def next_residue(grid: gemmi.FloatGrid, residue: gemmi.Residue):
    """((N, CA, C), confidence) of the next residue, else None.

    Positions are numpy arrays. The confidence is the length of the (sin, cos) output
    for psi: a regression onto unit vectors averages over the modes it cannot tell
    apart, so a short vector means the network is unsure.
    """
    predicted = _predict(grid, residue, "next")
    if predicted is None:
        return None
    (n, ca, c), (psi, omega, phi), confidence = predicted
    n2 = pos4(n, ca, c, length=1.32, angle=1.99, torsion=psi)
    ca2 = pos4(ca, c, n2, length=1.47, angle=2.15, torsion=omega)
    c2 = pos4(c, n2, ca2, length=1.53, angle=1.92, torsion=phi)
    return (n2, ca2, c2), confidence


def prev_residue(grid: gemmi.FloatGrid, residue: gemmi.Residue):
    "((N, CA, C), confidence) of the previous residue, else None; see next_residue."
    predicted = _predict(grid, residue, "prev")
    if predicted is None:
        return None
    (n, ca, c), (psi, omega, phi), confidence = predicted
    c0 = pos4(c, ca, n, length=1.32, angle=2.15, torsion=phi)
    ca0 = pos4(ca, n, c0, length=1.53, angle=1.99, torsion=omega)
    n0 = pos4(n, c0, ca0, length=1.47, angle=1.92, torsion=psi)
    return (n0, ca0, c0), confidence
