"""Find residues in a map with Buccaneer's fast search (Ca_find with SSfind).

Three-residue helix and strand templates are matched to the map in every rotation
on a 24 degree grid, using the gap between the lowest density at points expected to
be high and the highest density at points expected to be low. The best orientation
at each grid point is rescored with the fast LLK target, a prior from an existing
model is added, clashes are removed and the best matches are refined with the slow
target.
"""

import dataclasses

import gemmi
import numpy as np

from .llk import STD_C, STD_CA, STD_N, LLKTarget

# (high, low) density points about a residue for each SSfind template
_ALPHA3 = [
    ((0.00, 0.00, 0.00), (3.00, 0.00, 0.00)),
    ((0.87, 0.00, 1.23), (-1.00, 2.00, 0.25)),
]
_BETA3 = [
    ((0.00, 0.00, 0.00), (2.00, 2.00, -0.25)),
    ((0.87, 0.00, 1.23), (3.00, 0.50, -0.25)),
]


@dataclasses.dataclass
class Found:
    score: float
    rotation: np.ndarray
    translation: np.ndarray

    def positions(self):
        "N, CA and C (numpy)"
        return [self.rotation @ p + self.translation for p in (STD_N, STD_CA, STD_C)]


def _place(abc, length, angle, torsion):
    "Position d with |cd| = length, angle bcd and torsion abcd (radians)"
    a, b, c = abc
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array(
        [
            -length * np.cos(angle),
            length * np.sin(angle) * np.cos(torsion),
            length * np.sin(angle) * np.sin(torsion),
        ]
    )
    return c + d[0] * bc + d[1] * m + d[2] * n


def _superpose(source, target):
    "Rotation and translation taking source (n, 3) onto target (n, 3)"
    sc, tc = source.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((source - sc).T @ (target - tc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1, 1, d]) @ u.T
    return rotation, tc - rotation @ sc


def _template(pairs, phi, psi, length=3):
    "(high, low) point arrays for a secondary structure template (SSfind::Target)"
    residue = [np.array(p) for p in ((0.83, 0.00, -1.18), (0, 0, 0), (0.87, 0, 1.23))]
    residues = []
    for _ in range(length):
        n, ca, c = residue
        n = _place((n, ca, c), 1.32, 1.99, psi)
        ca = _place((ca, c, n), 1.47, 2.15, np.pi)
        c = _place((c, n, ca), 1.53, 1.92, phi)
        residue = [n, ca, c]
        residues.append(np.array(residue))
    middle = residues[(length - 1) // 2]
    ops = [_superpose(middle, r) for r in residues]
    high = np.array([r @ np.array(h) + t for h, _ in pairs for r, t in ops])
    low = np.array([r @ np.array(lo) + t for _, lo in pairs for r, t in ops])
    return high, low


def _rotation(axis: int, angle: float) -> np.ndarray:
    "Rotation by an angle (radians) about the x (0), y (1) or z (2) axis"
    i, j = [k for k in range(3) if k != axis]
    matrix = np.eye(3)
    matrix[i, i] = matrix[j, j] = np.cos(angle)
    matrix[i, j], matrix[j, i] = -np.sin(angle), np.sin(angle)
    if axis == 1:
        matrix = matrix.T
    return matrix


def _euler_ccp4(alpha, beta, gamma):
    "Rotation from CCP4 (z-y-z) Euler angles in degrees"
    a, b, c = np.radians([alpha, beta, gamma]).astype(float)
    return _rotation(2, a) @ _rotation(1, b) @ _rotation(2, c)


def rotations(step: float = 24.0) -> np.ndarray:
    "Uniformly sampled rotations (LLK_map_target::rtop_list in P1)"
    f = np.float32
    step = f(step)
    result = []
    bdeg = step / f(2)
    while bdeg < f(180):
        beta = np.radians(bdeg)
        spl = f(360) / np.floor(np.cos(f(0.5) * beta) * f(360) / step + f(1))
        smi = f(360) / np.floor(np.sin(f(0.5) * beta) * f(360) / step + f(1))
        thpl = spl / f(2)
        while thpl < f(720):
            thmi = smi / f(2)
            while thmi < f(360):
                adeg = np.mod(f(0.5) * (thpl + thmi), f(360))
                gdeg = np.mod(f(0.5) * (thpl - thmi), f(360))
                if adeg <= 360 and bdeg <= 180 and gdeg <= 360:
                    result.append(_euler_ccp4(adeg, bdeg, gdeg).T)
                thmi += smi
            thpl += spl
        bdeg += step
    return np.array(result)


def _asu_points(grid: gemmi.FloatGrid) -> np.ndarray:
    "Grid indices (n, 3) of one copy of each symmetry-equivalent point"
    mask = np.array(grid.masked_asu().mask_array, copy=False)
    return np.argwhere(mask == 0)


def _ssfind(grid: gemmi.FloatGrid, points: np.ndarray, ops: np.ndarray):
    "Best (score, rotation index) at each point over both templates"
    array = np.array(grid, copy=False)
    shape = np.array(array.shape)
    to_grid = np.array(grid.unit_cell.frac.mat.tolist()) * shape[:, None]
    best_score = np.zeros(len(points))
    best_rot = np.full(len(points), -1)
    for pairs, phi, psi in ((_ALPHA3, -58.0, -47.0), (_BETA3, -120.0, 120.0)):
        high, low = _template(pairs, np.radians(phi), np.radians(psi))
        score = np.zeros(len(points))
        rot = np.full(len(points), -1)
        for r, op in enumerate(ops):
            hi_off = np.rint((high @ op.T) @ to_grid.T).astype(int)
            lo_off = np.rint((low @ op.T) @ to_grid.T).astype(int)
            hi = np.full(len(points), np.inf)
            lo = np.full(len(points), -np.inf)
            for off in hi_off:
                idx = (points + off) % shape
                hi = np.minimum(hi, array[idx[:, 0], idx[:, 1], idx[:, 2]])
            for off in lo_off:
                idx = (points + off) % shape
                lo = np.maximum(lo, array[idx[:, 0], idx[:, 1], idx[:, 2]])
            better = hi - lo > score
            score[better] = (hi - lo)[better]
            rot[better] = r
        # Buccaneer keeps the beta result unless the alpha score is higher
        use = (score >= best_score) if pairs is _BETA3 else (score > best_score)
        best_score[use] = score[use]
        best_rot[use] = rot[use]
    return best_score, best_rot


def _prob_dist(x):
    return 0.999 * np.exp(-75.0 * (x - 3.5) ** 2 * x**-2.5) + 0.001


def _prior(structure: gemmi.Structure, positions: np.ndarray, radius: float):
    "-log of the prior for finding a residue near an existing model"
    if structure is None or structure[0].count_atom_sites() == 0:
        return np.zeros(len(positions))
    structure.setup_cell_images()
    ns = gemmi.NeighborSearch(structure[0], structure.cell, radius).populate()
    distances = np.full(len(positions), radius)
    for i, p in enumerate(positions):
        position = gemmi.Position(*p)
        mark = ns.find_nearest_atom(position, radius)
        if mark is not None:
            image = structure.cell.find_nearest_pbc_image(position, mark.pos, 0)
            distances[i] = min(radius, image.dist())
    return -np.log(_prob_dist(distances))


class _Clashes:
    "Symmetry-aware check for points within a distance of earlier accepted points"

    def __init__(self, cell: gemmi.UnitCell, spacegroup: gemmi.SpaceGroup, dist):
        self.frac = np.array(cell.frac.mat.tolist())
        self.orth = np.array(cell.orth.mat.tolist())
        ops = spacegroup.operations()
        self.ops = [
            (np.array(op.rot) / op.DEN, np.array(op.tran) / op.DEN) for op in ops
        ]
        self.dist = dist
        self.copies = np.zeros((0, 3))

    def clashes(self, position) -> bool:
        if len(self.copies) == 0:
            return False
        diff = self.copies - self.frac @ position
        diff -= np.rint(diff)
        return bool(((diff @ self.orth.T) ** 2).sum(1).min() < self.dist**2)

    def add(self, position) -> None:
        frac = self.frac @ position
        copies = [rot @ frac + tran for rot, tran in self.ops]
        self.copies = np.vstack([self.copies, copies])


def _euler_xyz(a, b, c):
    "Rotation by a about x, then b about y, then c about z (static axes)"
    return _rotation(2, c) @ _rotation(1, b) @ _rotation(0, a)


def _simplex_step(function, centre, shift, worst, best):
    "Reflected, extended or contracted point (and value) better than the worst"
    t1 = centre + shift
    f1 = function(t1)
    if f1 < worst:
        if f1 < best:
            t2 = centre + 2 * shift
            f2 = function(t2)
            if f2 < f1:
                return t2, f2
        return t1, f1
    t0 = centre - 0.5 * shift
    f0 = function(t0)
    if f0 < worst:
        return t0, f0
    return None, None


def _simplex(function, simplex, tolerance, max_cycles=50):
    "Buccaneer's gradient-weighted simplex (Optimiser_simplex::GRADIENT)"
    params = [np.array(p, dtype=float) for p in simplex]
    values = [function(p) for p in params]
    for _ in range(max_cycles):
        iw = int(np.argmax(values))
        ib = int(np.argmin(values))
        if values[iw] - values[ib] < tolerance:
            break
        others = [i for i in range(len(params)) if i != iw]
        weights = [
            (values[iw] - values[i]) / np.linalg.norm(params[i] - params[iw])
            for i in others
        ]
        centre = sum(w * params[i] for w, i in zip(weights, others)) / sum(weights)
        shift = centre - params[iw]
        point, value = _simplex_step(function, centre, shift, values[iw], values[ib])
        if point is not None:
            params[iw], values[iw] = point, value
        else:
            for i in others:
                if i != ib:
                    params[i] = 0.5 * (params[i] + params[ib])
                    values[i] = function(params[i])
    return params[int(np.argmin(values))]


def refine(grid, target: LLKTarget, rotation, translation, step=0.2):
    "Rotation and translation refined against the slow target (simplex)"

    def frame(args):
        return _euler_xyz(*args[:3]) @ rotation, translation + args[3:]

    def function(args):
        return float(target.correl(grid, *frame(args))[0])

    simplex = [np.zeros(6)]
    for i in range(6):
        args = np.zeros(6)
        args[i] = 0.5 * step
        simplex.append(args)
    best = _simplex(function, simplex, tolerance=0.005 * function(simplex[0]))
    return frame(best)


def default_num_find(grid: gemmi.FloatGrid, max_find=500, per_100_residues=20):
    "Buccaneer's number of residues to find (assuming 320 A^3 per residue)"
    volume = grid.unit_cell.volume / len(grid.spacegroup.operations())
    return min(max_find, per_100_residues * int(volume / 320) // 100)


def find(
    grid: gemmi.FloatGrid,
    target: LLKTarget,
    num_find: int = None,
    structure: gemmi.Structure = None,
    ops: np.ndarray = None,
) -> list:
    """Residues found in the map, best first, as a list of Found.

    An existing model, if given, makes positions about 3.5 A from it more likely.
    """
    if num_find is None:
        num_find = default_num_find(grid)
    if ops is None:
        ops = rotations(24.0)
    array = np.array(grid, copy=False)
    to_orth = np.array(grid.unit_cell.orth.mat.tolist()) / np.array(array.shape)

    # Search at points above 1 sigma, then rescore with the fast target
    asu = _asu_points(grid)
    values = array[asu[:, 0], asu[:, 1], asu[:, 2]]
    cut = values.mean() + values.std()
    above = asu[values > cut]
    _, rot = _ssfind(grid, above, ops)
    found = rot >= 0
    points, rot = above[found], rot[found]
    positions = points @ to_orth.T
    llk = target.correl(grid, ops[rot], positions, fast=True)
    zstep = max(len(asu) // 1000, 1)
    sample = asu[::zstep] @ to_orth.T
    identity = np.broadcast_to(np.eye(3), (len(sample), 3, 3))
    zvalues = target.correl(grid, identity, sample, fast=True)
    zscore = (llk - zvalues.mean()) / zvalues.std()

    # Add the prior and keep the best non-clashing matches
    score = 2.0 * zscore + _prior(structure, positions, 6.0)
    order = np.argsort(score, kind="stable")[: 20 * num_find]
    clashes = _Clashes(grid.unit_cell, grid.spacegroup, 5.0)
    kept = []
    for i in order:
        if len(kept) == num_find:
            break
        if not clashes.clashes(positions[i]):
            clashes.add(positions[i])
            kept.append(i)

    # Refine against the slow target
    results = []
    for i in kept:
        rotation, translation = refine(grid, target, ops[rot[i]], positions[i])
        final = float(target.correl(grid, rotation, translation)[0])
        results.append(Found(final, rotation, translation))
    results.sort(key=lambda r: r.score)
    return results


def found_structure(results, cell, spacegroup) -> gemmi.Structure:
    "A structure with one UNK residue (N, CA, C) per result"
    structure = gemmi.Structure()
    structure.cell = cell
    structure.spacegroup_hm = spacegroup.hm
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for i, result in enumerate(results):
        residue = gemmi.Residue()
        residue.name = "UNK"
        residue.seqid = gemmi.SeqId(2 * i + 2, " ")
        u_iso = np.exp(-10.0 * (result.score - results[0].score))
        for name, position in zip(("N", "CA", "C"), result.positions()):
            atom = gemmi.Atom()
            atom.name = name
            atom.element = gemmi.Element(name)
            atom.pos = gemmi.Position(*position)
            atom.b_iso = 8 * np.pi**2 * u_iso
            residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    return structure
