"""
Amino-acid type probabilities for a residue from the density around it.
Weights from findMySequence (EMBL/Grzegorz Chojnowski, BSD-3).
"""

import functools
from pathlib import Path

import gemmi
import numpy as np

from .maps import sample_probes

_WEIGHTS = Path(__file__).parent / "data" / "findmysequence.npz"

_REFERENCE = [
    gemmi.Position(11.596, 14.410, 9.372),  # N
    gemmi.Position(13.064, 14.381, 9.353),  # CA
    gemmi.Position(13.625, 13.209, 8.604),  # C
]


@functools.cache
def _network(em: bool):
    "Labels, probe cloud and layers; loaded on first use."
    w = np.load(_WEIGHTS)
    prefix = "Em" if em else "Xray"
    labels = list(w["Labels"])
    hidden = len(w[prefix + "Bias1"])
    cloud = w["Cloud"].astype(np.float64)
    layers = (
        w[prefix + "Kernel1"].reshape(324, hidden),
        w[prefix + "Bias1"],
        w[prefix + "Kernel2"].reshape(hidden, len(labels)),
        w[prefix + "Bias2"],
    )
    return labels, cloud, layers


def type_probabilities(grid: gemmi.FloatGrid, residue: gemmi.Residue, em: bool = False):
    "{type: probability} for a residue with N, CA and C, else None."
    atoms = [residue.find_atom(name, "*") for name in ("N", "CA", "C")]
    if any(atom is None for atom in atoms):
        return None
    labels, cloud, (kernel1, bias1, kernel2, bias2) = _network(em)
    values = sample_probes(grid, [a.pos for a in atoms], _REFERENCE, cloud)
    hidden = np.maximum(values @ kernel1 + bias1, 0)
    logits = hidden @ kernel2 + bias2
    exp = np.exp(logits - logits.max())
    return dict(zip(labels, (exp / exp.sum()).tolist()))
