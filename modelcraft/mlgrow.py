"The next or previous residue's N, CA and C from the density around a residue."

import functools
from pathlib import Path

import gemmi
import numpy as np

from .mlutil import REFERENCE, reference_transform, sample_probes

_DATA = Path(__file__).parent / "data"


@functools.cache
def _network(direction: str):
    "Probes, trunk layers and heads (float32) of one direction; loaded on first use."
    w = np.load(_DATA / f"mlgrow-{direction}.npz")

    def layer(name):
        kernel, bias = w[f"{name}_kernel"], w[f"{name}_bias"]
        return kernel.astype(np.float32), bias.astype(np.float32)

    probes = w["probes"].astype(np.float64)
    layers = [layer(name) for name in ("dense", "dense_1", "dense_2", "dense_3")]
    return probes, layers, layer("frame"), layer("error")


def _grow(grid: gemmi.FloatGrid, residue: gemmi.Residue, direction: str):
    atoms = [residue.find_atom(name, "*") for name in ("N", "CA", "C")]
    if any(atom is None for atom in atoms):
        return None
    probes, layers, (fk, fb), (ek, eb) = _network(direction)
    positions = [a.pos for a in atoms]
    x = sample_probes(grid, positions, probes)
    for w, b in layers:
        x = np.maximum(x @ w + b, 0)
    rotation, translation = reference_transform(positions)
    predicted = (x @ fk + fb).astype(np.float64).reshape(3, 3)
    predicted = predicted @ rotation.T + translation
    rotation, translation = reference_transform([gemmi.Position(*p) for p in predicted])
    reference = np.array([p.tolist() for p in REFERENCE])
    placed = reference @ rotation.T + translation
    return tuple(placed), float((x @ ek + eb)[0])


def next_residue(grid: gemmi.FloatGrid, residue: gemmi.Residue):
    """((N, CA, C), error) of the next residue, else None.

    Positions are numpy arrays. The network predicts the three atoms in the frame of
    the residue, and REFERENCE is superposed onto them, so the peptide bond is left
    for refinement. The error is the network's estimate of the RMSD of its
    prediction (Å).
    """
    return _grow(grid, residue, "next")


def prev_residue(grid: gemmi.FloatGrid, residue: gemmi.Residue):
    "((N, CA, C), error) of the previous residue, else None; see next_residue."
    return _grow(grid, residue, "prev")
