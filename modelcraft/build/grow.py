"The next or previous residues' N, CA and C from the density around residues."

import functools
from pathlib import Path

import gemmi
import numpy as np

from ..geometry import superpose
from .network import REFERENCE, place_reference, sample_probes

_DATA = Path(__file__).parent / "data"


@functools.cache
def _network(direction: str):
    "Probes, trunk layers and heads (float32) of one direction; loaded on first use."
    w = np.load(_DATA / f"grow-{direction}.npz")

    def layer(name):
        kernel, bias = w[f"{name}_kernel"], w[f"{name}_bias"]
        return kernel.astype(np.float32), bias.astype(np.float32)

    probes = w["probes"].astype(np.float64)
    layers = [layer(name) for name in ("dense", "dense_1", "dense_2", "dense_3")]
    return probes, layers, layer("frame"), layer("error")


def _grow(grid: gemmi.FloatGrid, residues: np.ndarray, direction: str):
    if len(residues) == 0:
        return np.zeros((0, 3, 3)), np.zeros(0)
    probes, layers, (fk, fb), (ek, eb) = _network(direction)
    x = sample_probes(grid, residues, probes)
    for w, b in layers:
        x = np.maximum(x @ w + b, 0)
    rotations, translations = superpose(REFERENCE, residues)
    predicted = (x @ fk + fb).astype(np.float64).reshape(-1, 3, 3)
    predicted = (
        np.matmul(predicted, rotations.transpose(0, 2, 1)) + translations[:, None]
    )
    return place_reference(predicted), (x @ ek + eb)[:, 0].astype(np.float64)


def next_residues(grid: gemmi.FloatGrid, residues: np.ndarray):
    """The next residues (b, 3, 3) and predicted errors (b,).

    The residues are N, CA and C positions (b, 3, 3). The network predicts the
    three atoms in the frame of each residue, and REFERENCE is superposed onto
    them, so the peptide bond is left for refinement. The error is the network's
    estimate of the RMSD of its prediction (Å).
    """
    return _grow(grid, residues, "next")


def prev_residues(grid: gemmi.FloatGrid, residues: np.ndarray):
    "The previous residues (b, 3, 3) and predicted errors (b,); see next_residues."
    return _grow(grid, residues, "prev")
