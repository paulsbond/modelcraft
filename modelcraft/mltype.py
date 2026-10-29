"Amino-acid type probabilities for a residue from the density around it."

import functools
from pathlib import Path

import gemmi
import numpy as np

from .mlutil import sample_probes

_WEIGHTS = Path(__file__).parent / "data" / "mltype.npz"

# The network's output order (findMySequence's)
_LABELS = [
    "CYS", "GLN", "ILE", "SER", "VAL", "LYS", "TRP", "PRO", "GLY", "THR",
    "ALA", "PHE", "HIS", "ASN", "LEU", "ARG", "ASP", "GLU", "TYR", "MET",
]  # fmt: skip


@functools.cache
def _network():
    "Probes and layers (float32); loaded on first use."
    w = np.load(_WEIGHTS)

    def layer(name):
        kernel, bias = w[f"{name}_kernel"], w[f"{name}_bias"]
        return kernel.astype(np.float32), bias.astype(np.float32)

    probes = w["probes"].astype(np.float64)
    names = ("dense", "dense_1", "dense_2", "dense_3", "residue")
    return probes, [layer(name) for name in names]


def type_probabilities(grid: gemmi.FloatGrid, residue: gemmi.Residue):
    "{type: probability} for a residue with N, CA and C, else None."
    atoms = [residue.find_atom(name, "*") for name in ("N", "CA", "C")]
    if any(atom is None for atom in atoms):
        return None
    probes, layers = _network()
    x = sample_probes(grid, [a.pos for a in atoms], probes)
    for w, b in layers[:-1]:
        x = np.maximum(x @ w + b, 0)
    logits = (x @ layers[-1][0] + layers[-1][1]).astype(np.float64)
    exp = np.exp(logits - logits.max())
    return dict(zip(_LABELS, (exp / exp.sum()).tolist()))
