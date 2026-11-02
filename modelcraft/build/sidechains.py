"""Side chains from the rotamers predicted by the side network.

Rotamers are chosen to balance the network's probabilities against clashes with
other atoms, including symmetry copies. Side chains that still clash badly are cut
back to CB, and carbonyl oxygens are added.

The rotamer library (data/rotamers.json) is in the order of the side network's
classes, so it can't be reordered without retraining. The rotamers are the Top8000
central values (Richardson Lab, CC BY 4.0), with ASN from the Penultimate library.
"""

import functools
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial.distance import cdist

from ..geometry import torsion_position
from ..stdres import frame
from ..symmetry import Symmetry
from . import side

_LIBRARY = Path(__file__).parent / "data" / "rotamers.json"
# Atoms closer than these distances (A) overlap: between two N or O atoms, and others
_POLAR_CONTACT = 2.5
_CONTACT = 3.4
# Penalty per A of overlap, in units of log probability
_OVERLAP_WEIGHT = 5.0
# Side chains with an atom closer than this to another residue are cut back to CB
_SEVERE = 2.0
_MIN_PROBABILITY = 1e-6
_ROUNDS = 10


@functools.cache
def _library():
    "{type: (atom names, positions (rotamers, atoms, 3), frequencies)}"
    with open(_LIBRARY) as stream:
        data = json.load(stream)
    library = {}
    for name, rotamers in data.items():
        atoms = list(rotamers[0]["atoms"])
        positions = np.array([[r["atoms"][a] for a in atoms] for r in rotamers])
        frequencies = np.array([r["probability"] for r in rotamers])
        library[name] = atoms, positions, frequencies
    return library


def conformations(name, n, ca, c):
    """Side-chain atom names, conformations (rotamers, atoms, 3) and library
    frequencies (rotamers,) of a residue type on N, CA and C"""
    if name not in _library():
        return [], np.zeros((1, 0, 3)), np.ones(1)
    atoms, positions, frequencies = _library()[name]
    rotation, origin = frame(n, ca, c)
    return atoms, positions @ rotation.T + origin, frequencies


def _polar(names):
    return np.array([name[:1] in ("N", "O") for name in names])


class _Residue:
    "A residue's main-chain atoms, rotamers and their costs, and the chosen rotamer"

    def __init__(self, index, residue, build_as, probabilities):
        self.index = index
        self.residue = residue
        self.main = {
            a: np.array(residue.find_atom(a, "*").pos.tolist())
            for a in ("N", "CA", "C")
        }
        n, ca, c = self.main["N"], self.main["CA"], self.main["C"]
        self.atoms, self.options, frequencies = conformations(build_as, n, ca, c)
        if probabilities is None:
            probabilities = frequencies
        self.costs = -np.log(np.maximum(probabilities, _MIN_PROBABILITY))
        self.previous = None  # the residue before, if PRO's ring bonds to it

    @property
    def beyond(self):
        "Indices of the side-chain atoms beyond CB"
        return [i for i, a in enumerate(self.atoms) if a != "CB"]

    def choose(self, option):
        "Use a rotamer, or cut the side chain back to CB (None)"
        self.option = option
        if option is None:
            keep = 1 if self.atoms[:1] == ["CB"] else 0
            self.names, self.side = self.atoms[:keep], self.side[:keep]
        else:
            self.names, self.side = self.atoms, self.options[option]

    def fixed(self):
        "Names and positions (atoms, 3) of the main-chain atoms and CB"
        names = ["N", "CA", "C", "O"]
        positions = [self.main[a] for a in names]
        if "CB" in self.atoms:
            names.append("CB")
            positions.append(self.options[0][self.atoms.index("CB")])
        return names, np.array(positions)

    def excluded(self):
        "Indices of the residues whose fixed atoms the side chain can bond to"
        return {self.index} | ({self.previous} if self.previous is not None else set())


def _bonded(residue, following):
    "Whether the following residue (or None) is bonded to the residue"
    if following is None:
        return False
    return np.sum((following.main["N"] - residue.main["C"]) ** 2) < 2.2


def _oxygen(residue, following):
    "Carbonyl O (as in Clipper), placed from the next residue's N if bonded"
    n, ca, c = residue.main["N"], residue.main["CA"], residue.main["C"]
    if _bonded(residue, following):
        n2 = following.main["N"]
        v0 = (ca - c) / np.linalg.norm(ca - c)
        v1 = np.cross(v0, np.cross(v0, n2 - c))
        return c - 0.637 * v0 + 1.064 * v1 / np.linalg.norm(v1)
    return torsion_position((n, ca, c), 1.24, 2.09, -0.58)


def _overlaps(distances, polar_a, polar_b):
    "Overlaps (A) of atoms at distances, broadcasting with whether each is N or O"
    contact = np.where(polar_a & polar_b, _POLAR_CONTACT, _CONTACT)
    return np.maximum(contact - distances, 0.0)


def _fixed_penalties(residues, sym):
    "Penalties (rotamers,) of each residue's rotamers for overlapping fixed atoms"
    names, positions, owners = [], [], []
    for r in residues:
        n, p = r.fixed()
        names += n
        positions.append(p)
        owners += [r.index] * len(n)
    polar, owners = _polar(names), np.array(owners)
    sym.index(np.concatenate(positions), _CONTACT)
    penalties = []
    for r in residues:
        beyond = r.beyond
        penalty = np.zeros(len(r.costs))
        if beyond:
            side = r.options[:, beyond]
            queries, found, d = sym.within(side.reshape(-1, 3), _CONTACT)
            keep = ~np.isin(owners[found], list(r.excluded()))
            queries, found, d = queries[keep], found[keep], d[keep]
            side_polar = np.tile(_polar([r.atoms[i] for i in beyond]), len(side))
            overlap = _overlaps(d, side_polar[queries], polar[found])
            np.add.at(penalty, queries // len(beyond), _OVERLAP_WEIGHT * overlap)
        penalties.append(penalty)
    return penalties


def _pair_penalties(residues, sym):
    """Penalties (rotamers of i, rotamers of j) for each pair of residues whose side
    chains overlap, as pairs[i][j]"""
    centres, reach, sides = [], [], []
    for r in residues:
        side = r.options[:, r.beyond].reshape(-1, 3)
        # The centre of a sphere about all the rotamers' atoms beyond CB
        centre = side.mean(axis=0) if len(side) else r.main["CA"]
        centres.append(centre)
        sides.append(side)
        reach.append(np.linalg.norm(side - centre, axis=-1).max(initial=0.0))
    centres, reach = np.array(centres), np.array(reach)
    pairs = [{} for _ in residues]
    if reach.max() == 0:
        return pairs
    sym.index(centres, 2 * reach.max() + _CONTACT)
    for i, r in enumerate(residues):
        if reach[i] == 0:
            continue
        polar_i = _polar([r.atoms[k] for k in r.beyond])
        radius = reach[i] + reach.max() + _CONTACT
        for j, transform in sym.near(centres[i], radius):
            if j <= i or reach[j] == 0:
                continue
            gap = np.linalg.norm(sym.apply(transform, centres[j]) - centres[i])
            if gap > reach[i] + reach[j] + _CONTACT:
                continue
            d = cdist(sides[i], sym.apply(transform, sides[j]))
            if d.min() >= _CONTACT:
                continue
            s = residues[j]
            a, b = len(r.beyond), len(s.beyond)
            d = d.reshape(len(r.costs), a, len(s.costs), b)
            polar_j = _polar([s.atoms[k] for k in s.beyond])
            overlap = _overlaps(d, polar_i[:, None, None], polar_j).sum((1, 3))
            if overlap.any():
                matrix = _OVERLAP_WEIGHT * overlap
                pairs[i][j] = pairs[i].get(j, 0) + matrix
                pairs[j][i] = pairs[j].get(i, 0) + matrix.T
    return pairs


def _choose(residues, unary, pairs):
    "Rotamer indices approximately minimising the unary costs plus pair penalties"
    choice = [int(np.argmin(cost)) for cost in unary]

    def given_others(i, exclude=None):
        cost = unary[i].copy()
        for j, matrix in pairs[i].items():
            if j != exclude:
                cost += matrix[:, choice[j]]
        return cost

    for _ in range(_ROUNDS):
        changed = False
        for i in range(len(residues)):
            if pairs[i]:
                best = int(np.argmin(given_others(i)))
                changed |= best != choice[i]
                choice[i] = best
        for i in range(len(residues)):
            for j, matrix in pairs[i].items():
                if j < i or matrix[choice[i], choice[j]] == 0:
                    continue
                cost = given_others(i, j)[:, None] + given_others(j, i)[None] + matrix
                a, b = np.unravel_index(np.argmin(cost), cost.shape)
                if cost[a, b] < cost[choice[i], choice[j]]:
                    choice[i], choice[j] = int(a), int(b)
                    changed = True
        if not changed:
            break
    return choice


def _cut_severe(residues, sym):
    """Cut back to CB side chains with an atom within _SEVERE of another residue,
    keeping the more probable of two clashing side chains"""
    names, positions, owners = [], [], []
    for r in residues:
        n, p = r.fixed()
        beyond = r.beyond if r.option is not None else []
        names += n + [r.names[k] for k in beyond]
        positions += [p, r.side[beyond].reshape(-1, 3)]
        owners += [(r.index, False)] * len(n) + [(r.index, True)] * len(beyond)
    sym.index(np.concatenate(positions), _SEVERE)
    clashes = []
    for r in residues:
        if r.option is None or not r.beyond:
            continue
        _, found, _ = sym.within(r.side[r.beyond], _SEVERE)
        for k in found:
            j, is_side = owners[k]
            if j in r.excluded():
                continue
            clashes.append((r.index, j, is_side))
    for i, j, is_side in clashes:
        first, second = residues[i], residues[j]
        if first.option is None or (is_side and second.option is None):
            continue
        if is_side and second.costs[second.option] > first.costs[first.option]:
            second.choose(None)
        else:
            first.choose(None)


def _write(r):
    "Replace the atoms of the gemmi residue other than N, CA and C"
    residue = r.residue
    # Copied before the residue changes, which can move its atoms in memory
    ca = residue.find_atom("CA", "*")
    b_iso, occ = ca.b_iso, ca.occ
    for k in reversed(range(len(residue))):
        if residue[k].name not in ("N", "CA", "C"):
            del residue[k]
    for name, position in [("O", r.main["O"]), *zip(r.names, r.side)]:
        atom = gemmi.Atom()
        atom.name = name
        atom.element = gemmi.Element("Se" if name.startswith("SE") else name[0])
        atom.pos = gemmi.Position(*position)
        atom.b_iso = b_iso
        atom.occ = occ
        residue.add_atom(atom)


def _residues(chains, probabilities, unknown):
    "_Residues for chains of gemmi residues, with one row of probabilities each"
    built, row = [], 0
    for residues in chains:
        chain = []
        for residue in residues:
            build_as = residue.name
            if build_as != "GLY" and build_as not in _library():
                build_as = unknown
            given = side.rotamer_probabilities(probabilities[row][None], build_as)
            given = None if given is None else given[0]
            chain.append(_Residue(row, residue, build_as, given))
            row += 1
        for k, r in enumerate(chain):
            following = chain[k + 1] if k + 1 < len(chain) else None
            r.main["O"] = _oxygen(r, following)
            if _bonded(r, following) and following.residue.name == "PRO":
                following.previous = r.index
        built += chain
    return built


def build(
    structure: gemmi.Structure,
    grid: gemmi.FloatGrid = None,
    *,
    probabilities=None,
    unknown="ALA",
):
    """Add side chains and carbonyl oxygens to a model of N, CA and C atoms, in place.

    The side network's probabilities for the residues are predicted from the grid if
    not given. Residues of types without rotamers (e.g. UNK) are built as unknown.
    """
    chains, rows = [], []
    for chain in structure[0]:
        residues = []
        for residue in chain:
            atoms = [residue.find_atom(a, "*") for a in ("N", "CA", "C")]
            if all(atoms):
                residues.append(residue)
                rows.append([a.pos.tolist() for a in atoms])
        chains.append(residues)
    if probabilities is None:
        if grid is None:
            raise ValueError("Either a grid or probabilities are needed")
        probabilities = side.probabilities(grid, np.array(rows).reshape(-1, 3, 3))
    if len(probabilities) != len(rows):
        raise ValueError("The probabilities don't match the residues")
    built = _residues(chains, probabilities, unknown)
    if not built:
        return
    sym = Symmetry(structure.cell, structure.find_spacegroup())
    unary = [r.costs + p for r, p in zip(built, _fixed_penalties(built, sym))]
    pairs = _pair_penalties(built, sym)
    for r, option in zip(built, _choose(built, unary, pairs)):
        r.choose(option)
    _cut_severe(built, sym)
    for r in built:
        _write(r)
