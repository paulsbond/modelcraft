"Distances and neighbours of points under crystal symmetry."

import gemmi
import numpy as np
from scipy.spatial import cKDTree


class Symmetry:
    "Find symmetry copies of points near a query, with the operator to apply"

    def __init__(self, cell: gemmi.UnitCell, spacegroup: gemmi.SpaceGroup):
        self.frac = np.array(cell.frac.mat.tolist())
        self.orth = np.array(cell.orth.mat.tolist())
        self.ops = [
            (np.array(op.rot) / op.DEN, np.array(op.tran) / op.DEN)
            for op in spacegroup.operations()
        ]

    def index(self, points: np.ndarray, radius: float):
        "Prepare to search points (n, 3) within a radius"
        frac = points @ self.frac.T
        margin = radius * np.linalg.norm(self.frac, axis=1)
        coords, owners, ops, offsets = [], [], [], []
        for k, (rot, tran) in enumerate(self.ops):
            copy = frac @ rot.T + tran
            wrap = np.floor(copy)
            copy -= wrap
            for cell_shift in np.ndindex(3, 3, 3):
                shift = np.array(cell_shift) - 1
                moved = copy + shift
                keep = np.all((moved > -margin) & (moved < 1 + margin), axis=1)
                coords.append(moved[keep])
                owners.append(np.flatnonzero(keep))
                ops.append(np.full(keep.sum(), k))
                offsets.append(tran - wrap[keep] + shift)
        self.coords = np.concatenate(coords)
        self.owners = np.concatenate(owners)
        self.op_index = np.concatenate(ops)
        self.offsets = np.concatenate(offsets)
        self.tree = cKDTree(self.coords @ self.orth.T)

    def near(self, query: np.ndarray, radius: float):
        "[(point index, transform)] for copies within the radius of a query point"
        fq = self.frac @ query
        wq = np.floor(fq)
        hits = self.tree.query_ball_point(self.orth @ (fq - wq), radius)
        return [
            (self.owners[h], (self.op_index[h], self.offsets[h] + wq)) for h in hits
        ]

    def any_near(self, points: np.ndarray, radius: float) -> np.ndarray:
        "Whether each of the points (n, 3) has an indexed point within the radius"
        fq = points @ self.frac.T
        query = (fq - np.floor(fq)) @ self.orth.T
        return self.tree.query_ball_point(query, radius, return_length=True) > 0

    def within(self, points: np.ndarray, radius: float):
        """Query indices, indexed point indices and distances (arrays) of the copies
        within the radius of each of the points (n, 3)"""
        fq = points @ self.frac.T
        query = (fq - np.floor(fq)) @ self.orth.T
        hits = self.tree.query_ball_point(query, radius)
        counts = np.array([len(h) for h in hits], dtype=int)
        queries = np.repeat(np.arange(len(points)), counts)
        found = np.array([h for row in hits for h in row], dtype=int)
        distances = np.linalg.norm(self.tree.data[found] - query[queries], axis=1)
        return queries, self.owners[found], distances

    def apply(self, transform, positions: np.ndarray) -> np.ndarray:
        "Positions moved by a (operator index, fractional offset) transform"
        k, offset = transform
        rot = self.ops[k][0]
        return ((positions @ self.frac.T) @ rot.T + offset) @ self.orth.T

    def distance(self, a: np.ndarray, b: np.ndarray) -> float:
        "Shortest distance between points a and b under symmetry"
        fa, fb = self.frac @ a, self.frac @ b
        best = np.inf
        for rot, tran in self.ops:
            diff = rot @ fb + tran - fa
            diff -= np.rint(diff)
            best = min(best, float(np.linalg.norm(self.orth @ diff)))
        return best

    def nearest_copy(self, positions: np.ndarray, atom, reference: np.ndarray):
        "Symmetry copy of positions (..., 3) with positions[atom] nearest reference"
        frac = positions @ self.frac.T
        fref = self.frac @ reference
        best, best_d2 = None, np.inf
        for rot, tran in self.ops:
            copy = frac @ rot.T + tran
            copy += np.rint(fref - copy[atom])
            d2 = np.sum(((copy[atom] - fref) @ self.orth.T) ** 2)
            if d2 < best_d2:
                best, best_d2 = copy, d2
        return best @ self.orth.T
