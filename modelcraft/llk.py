"""Buccaneer's log-likelihood (LLK) map targets.

A target is a sphere of density statistics in the frame of a residue, accumulated
from a reference map at every residue of a reference structure (LLK_map_target and
Ca_prep in Buccaneer). It is scored in a work map with either the log-likelihood or,
as Buccaneer always does now, the negative weighted correlation (lower is better).
"""

import os

import gemmi
import numpy as np

from .mapsim import simulate
from .reflections import DataItem

# N, CA, C and CB in Buccaneer's standard residue frame (Ca_group)
STD_N = np.array([0.8431, 0.0000, -1.2040])
STD_CA = np.array([0.0000, 0.0000, 0.0000])
STD_C = np.array([0.8775, 0.0000, 1.2530])
STD_CB = np.array([-1.03, -1.11, 0.00])

# Buccaneer's order of the residue types for the side-chain targets
TYPES = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]  # fmt: skip


def frame(n, ca, c):
    """Rotation and translation (numpy) taking the standard frame onto a residue.

    The frame is built from the geometry (Ca_group::rtop_from_std_ori): x bisects
    CA->C and CA->N, y is their cross product and the origin is CA.
    """
    n, ca, c = (np.asarray(p, dtype=float) for p in (n, ca, c))
    c1, c2 = c - ca, n - ca
    v1 = c1 / np.linalg.norm(c1) + c2 / np.linalg.norm(c2)
    v1 /= np.linalg.norm(v1)
    v2 = np.cross(c1, c2)
    v2 /= np.linalg.norm(v2)
    v3 = np.cross(v1, v2)
    v3 /= np.linalg.norm(v3)
    return np.column_stack([v1, v2, v3]), ca


def _cubic(array: np.ndarray, point: np.ndarray) -> float:
    "Catmull-Rom interpolation of a non-periodic array at a grid point"
    base = np.floor(point).astype(int)
    t = point - base
    weights = [
        np.array([
            -0.5 * x * (1 - x) ** 2,
            (1 - x) * (-1.5 * x * x + x + 1),
            x * (-1.5 * (1 - x) ** 2 + (1 - x) + 1),
            -0.5 * x * x * (1 - x),
        ])
        for x in t
    ]  # fmt: skip
    block = array[tuple(slice(b - 1, b + 3) for b in base)]
    return float(np.einsum("ijk,i,j,k", block, *weights))


class LLKTarget:
    "Density statistics within a radius of a frame, sampled on a cubic grid"

    def __init__(self, radius: float, sampling: float = 0.5):
        self.radius = radius
        self.sampling = sampling
        self.ng = int(np.rint((radius + 2 * sampling) / sampling))
        steps = (np.arange(2 * self.ng + 1) - self.ng) * sampling
        self.points = np.stack(np.meshgrid(steps, steps, steps, indexing="ij"), -1)
        self.target = np.zeros(self.points.shape[:3])
        self.weight = np.zeros(self.points.shape[:3])
        self.count = 0
        self.fast = self.slow = None

    def accumulate(self, grid: gemmi.FloatGrid, frames) -> None:
        """Add density (cubic interpolation) about each (rotation, translation).

        Before prepare(), target and weight hold sums of the density and its square.
        """
        extent = self.ng * self.sampling
        inside = (self.points**2).sum(-1) <= extent**2
        local = self.points[inside]
        for rotation, translation in frames:
            xyz = local @ rotation.T + translation
            rho = grid.interpolate_position_array(xyz, order=3).astype(np.float64)
            self.target[inside] += rho
            self.weight[inside] += rho * rho
            self.count += 1

    def prepare(self) -> None:
        "Turn the sums into a target and weight, then make the fast and slow lists"
        if self.count > 0:
            used = self.weight > 0
            sn = self.count * used.sum()
            sr = self.target[used].sum()
            sr2 = self.weight[used].sum()
            rmap = sr / sn
            smap = np.sqrt(sn * sr2 - sr * sr) / sn
            rho = self.target[used] / self.count
            std = np.sqrt(
                np.maximum(self.weight[used] / self.count - rho * rho, 0.001 * smap**2)
            )
            w1 = np.maximum(smap**2 / std**2 - 1.0, 0.001)
            w2 = np.minimum(1.0 / w1, 2.0)
            self.target[used] = rho + w2 * (rho - rmap)
            self.weight[used] = w1 * 0.5 / smap**2
            self.count = 0
        outside = (self.points**2).sum(-1) > self.radius**2
        self.target[outside] = self.weight[outside] = 0.0

        # Fast list: 13 points about the centre (cubic interpolation)
        r0 = self.radius * 3.0 / 8.0
        offsets = [
            (u, v, w)
            for u in (-1, 0, 1)
            for v in (-1, 0, 1)
            for w in (-1, 0, 1)
            if (u + v + w) % 2 == 0
        ]
        xyz = r0 * np.array(offsets, dtype=float)
        index = xyz / self.sampling + self.ng
        self.fast = (
            xyz,
            np.array([_cubic(self.target, i) for i in index]),
            np.array([_cubic(self.weight, i) for i in index]),
        )

        # Slow list: alternate grid points within the radius
        irad = int(np.rint(self.radius / self.sampling))
        index = np.indices(self.target.shape)
        mask = (index.sum(0) % 2 == 0) & (((index - self.ng) ** 2).sum(0) <= irad**2)
        self.slow = (self.points[mask], self.target[mask], self.weight[mask])

    def _values(self, grid, rotations, translations, fast):
        "Points, target, weight and (frames, points) linear map values"
        xyz, target, weight = self.fast if fast else self.slow
        rotations = np.asarray(rotations).reshape(-1, 3, 3)
        translations = np.asarray(translations).reshape(-1, 3)
        positions = np.einsum("fij,pj->fpi", rotations, xyz) + translations[:, None]
        values = grid.interpolate_position_array(positions.reshape(-1, 3), order=1)
        return target, weight, values.reshape(len(rotations), -1).astype(np.float64)

    def llk(self, grid, rotations, translations, fast=False) -> np.ndarray:
        "Weighted mean squared difference for a batch of frames (lower is better)"
        target, weight, y = self._values(grid, rotations, translations, fast)
        return (weight * (y - target) ** 2).sum(1) / weight.sum()

    def correl(self, grid, rotations, translations, fast=False) -> np.ndarray:
        "Negative weighted correlation for a batch of frames (lower is better)"
        x, w, y = self._values(grid, rotations, translations, fast)
        sw, swx, swxx = w.sum(), (w * x).sum(), (w * x * x).sum()
        swy, swyy, swxy = y @ w, (y * y) @ w, y @ (w * x)
        denominator = (sw * swxx - swx * swx) * (sw * swyy - swy * swy)
        return -(sw * swxy - swx * swy) / np.sqrt(np.maximum(denominator, 1e-20))


def _peptide_bond(r1: gemmi.Residue, r2: gemmi.Residue) -> bool:
    c, n = r1.find_atom("C", "*"), r2.find_atom("N", "*")
    return c is not None and n is not None and c.pos.dist(n.pos) < 1.5


def prepare_targets(
    structure: gemmi.Structure,
    grid: gemmi.FloatGrid,
    main_radius: float = 4.0,
    side_radius: float = 5.5,
):
    """The main-chain target and {type: side-chain target} (Ca_prep).

    Residues must have a peptide bond on both sides. The main-chain target is about
    CA in all non-GLY/PRO residues, and the side-chain targets are about the
    standard CB position of each type (in the same orientation).
    """
    main_frames = []
    side_frames = {name: [] for name in TYPES}
    for chain in structure[0]:
        residues = list(chain.first_conformer())
        for prev, res, next_ in zip(residues, residues[1:], residues[2:]):
            name = "MET" if res.name == "MSE" else res.name
            atoms = [res.find_atom(a, "*") for a in ("N", "CA", "C")]
            if name not in side_frames or any(a is None for a in atoms):
                continue
            if not (_peptide_bond(prev, res) and _peptide_bond(res, next_)):
                continue
            rotation, translation = frame(*(a.pos.tolist() for a in atoms))
            if name not in ("GLY", "PRO"):
                main_frames.append((rotation, translation))
            side_frames[name].append((rotation, translation + rotation @ STD_CB))
    main = LLKTarget(main_radius)
    main.accumulate(grid, main_frames)
    main.prepare()
    sides = {}
    for name, frames in side_frames.items():
        sides[name] = LLKTarget(side_radius)
        sides[name].accumulate(grid, frames)
        sides[name].prepare()
    return main, sides


def reference_targets(
    fsigf: DataItem, phases: DataItem, resolution: float = 2.0, seed: int = 0
):
    """Targets from Buccaneer's reference structure (1tqw) degraded like the work data.

    As in Buccaneer, the resolution is limited by both the reference and work data.
    """
    directory = os.path.join(os.environ["CLIBD"], "reference_structures")
    mtz = gemmi.read_mtz_file(os.path.join(directory, "reference-1tqw.mtz"))
    resolution = max(resolution, mtz.resolution_high(), fsigf.resolution_high())
    grid = simulate(
        DataItem(mtz, "FP.F_sigF.F,FP.F_sigF.sigF"),
        DataItem(mtz, "FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D"),
        fsigf,
        phases,
        resolution,
        seed=seed,
    )
    structure = gemmi.read_structure(os.path.join(directory, "reference-1tqw.pdb"))
    return prepare_targets(structure, grid)
