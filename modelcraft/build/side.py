"""Residue types and side-chain rotamers from the density around a residue.

The network gives the probabilities of (type, rotamer) classes for a residue from
the density about its N, CA and C. The classes are the types in LABELS order,
each with its rotamers in the order of data/rotamers.json (GLY and ALA have one).
"""

import functools
import json
from pathlib import Path

import numpy as np

from .network import sample_probes

_DATA = Path(__file__).parent / "data"

LABELS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]  # fmt: skip
_SAME_AS = {"MSE": "MET"}


@functools.cache
def _classes():
    "First class (20,) and number of rotamers (20,) of each type in LABELS"
    with open(_DATA / "rotamers.json") as stream:
        library = json.load(stream)
    counts = np.array([len(library.get(name, [None])) for name in LABELS])
    return np.cumsum(counts) - counts, counts


@functools.cache
def _network():
    "Probes and layers (float32); loaded on first use."
    w = np.load(_DATA / "side.npz")

    def layer(name):
        kernel, bias = w[f"{name}_kernel"], w[f"{name}_bias"]
        return kernel.astype(np.float32), bias.astype(np.float32)

    probes = w["probes"].astype(np.float64)
    names = ("dense", "dense_1", "dense_2", "dense_3", "rotamer")
    layers = [layer(name) for name in names]
    if layers[-1][1].shape[0] != _classes()[1].sum():
        raise ValueError("side.npz doesn't match the classes of rotamers.json")
    return probes, layers


def probabilities(grid, residues: np.ndarray) -> np.ndarray:
    "Probabilities (b, classes) of the (type, rotamer) classes for residues (b, 3, 3)"
    probes, layers = _network()
    if len(residues) == 0:
        return np.zeros((0, len(layers[-1][1])))
    x = sample_probes(grid, residues, probes)
    for w, b in layers[:-1]:
        x = np.maximum(x @ w + b, 0)
    logits = (x @ layers[-1][0] + layers[-1][1]).astype(np.float64)
    exp = np.exp(logits - logits.max(axis=1, keepdims=True))
    return exp / exp.sum(axis=1, keepdims=True)


def type_probabilities(joint: np.ndarray) -> np.ndarray:
    "Probabilities (b, 20) of the types in LABELS from class probabilities (b, classes)"
    starts, _ = _classes()
    return np.add.reduceat(joint, starts, axis=1)


def rotamer_probabilities(joint: np.ndarray, name: str):
    """Probabilities (b, rotamers) of the rotamers of a type, given the type, in the
    order of data/rotamers.json, from class probabilities (b, classes). None for a
    type the network doesn't predict."""
    name = _SAME_AS.get(name, name)
    if name not in LABELS:
        return None
    starts, counts = _classes()
    k = LABELS.index(name)
    rows = joint[:, starts[k] : starts[k] + counts[k]]
    return rows / np.maximum(rows.sum(axis=1, keepdims=True), 1e-300)
