import collections
import dataclasses
import functools
import math
import re
import gemmi
from .contents import AsuContents, Polymer, PolymerType
from .monlib import chemcomp


def solvent_fraction(contents: AsuContents, mtz: gemmi.Mtz) -> float:
    volume = _contents_volume(contents)
    asu_volume = mtz.cell.volume / len(mtz.spacegroup.operations())
    copies = contents.copies or _guess_copies(contents, mtz)
    return 1 - copies * volume / asu_volume


@functools.lru_cache(maxsize=None)
def _library_weight(code: str) -> float:
    return sum(atom.el.weight for atom in chemcomp(code).atoms)


@functools.lru_cache(maxsize=None)
def _library_volume(code: str) -> float:
    return sum(18 for atom in chemcomp(code).atoms if not atom.is_hydrogen())


def _polymer_weight(polymer: Polymer) -> float:
    codes = polymer.residue_codes(modified=False)
    total = sum(_library_weight(code) for code in codes)
    total -= _library_weight("HOH") * (len(codes) - 1)
    return total


def _polymer_volume(polymer: Polymer) -> float:
    density = 1.35 if polymer.type == PolymerType.PROTEIN else 2.0
    return _polymer_weight(polymer) / (density * 0.602214)


def _smiles_volume(smiles: str) -> float:
    atoms = re.findall(pattern="[A-Z][a-z]?", string=smiles)
    return 18 * len(atoms)


def _contents_volume(contents: AsuContents) -> float:
    return sum(
        item.volume * item.stoichiometry for item in _volume_components(contents)
    )


@dataclasses.dataclass
class _VolumeComponent:
    description: str
    stoichiometry: int
    stoichiometry_assumed: bool
    volume: float


def _volume_components(contents: AsuContents):
    for kind, polymers in (
        ("Protein", contents.proteins),
        ("RNA", contents.rnas),
        ("DNA", contents.dnas),
    ):
        for polymer in polymers:
            sequence = polymer.sequence
            description = f"{kind} with {len(sequence)} residues: "
            if len(sequence) > 9:
                description += f"{sequence[:3]}...{sequence[-3:]}"
            else:
                description += f"{sequence:9}"
            stoichiometry = polymer.stoichiometry or 1
            stoichiometry_assumed = polymer.stoichiometry is None
            volume = _polymer_volume(polymer)
            yield _VolumeComponent(
                description, stoichiometry, stoichiometry_assumed, volume
            )
    for carb in contents.carbs:
        description = "Carb:"
        stoichiometry = carb.stoichiometry or 1
        stoichiometry_assumed = carb.stoichiometry is None
        volume = 0
        length = 0
        for code, count in carb.codes.items():
            description += f" {count}x{code}"
            length += count
            if code in contents.smiles:
                volume += _smiles_volume(contents.smiles[code]) * count
            else:
                volume += _library_volume(code) * count
        volume -= _library_volume("HOH") * length
        yield _VolumeComponent(
            description, stoichiometry, stoichiometry_assumed, volume
        )
    for ligand in contents.ligands:
        description = "Ligand: " + ligand.code
        stoichiometry = ligand.stoichiometry or 1
        stoichiometry_assumed = ligand.stoichiometry is None
        if ligand.code in contents.smiles:
            volume = _smiles_volume(contents.smiles[ligand.code])
        else:
            volume = _library_volume(ligand.code)
        yield _VolumeComponent(
            description, stoichiometry, stoichiometry_assumed, volume
        )


@dataclasses.dataclass
class _CopiesOption:
    copies: int
    solvent: float
    probability: float


def _copies_options(contents: AsuContents, mtz: gemmi.Mtz) -> list:
    options = []
    nucleic_acids = contents.rnas + contents.dnas
    mwp = sum(_polymer_weight(p) * (p.stoichiometry or 1) for p in contents.proteins)
    mwn = sum(_polymer_weight(n) * (n.stoichiometry or 1) for n in nucleic_acids)
    asu_volume = mtz.cell.volume / len(mtz.spacegroup.operations())
    contents_volume = _contents_volume(contents)
    resolution = mtz.resolution_high()
    total_probability = 0
    for copies in range(1, 60):
        solvent = 1 - copies * contents_volume / asu_volume
        probability = _probability(mwp, mwn, copies, asu_volume, resolution)
        if solvent < 0:
            break
        options.append(_CopiesOption(copies, solvent, probability))
        total_probability += probability
    for option in options:
        option.probability /= total_probability
    return options


def _guess_copies(contents: AsuContents, mtz: gemmi.Mtz) -> int:
    options = _copies_options(contents, mtz)
    if len(options) == 0:
        raise ValueError("Contents are too big to fit into the asymmetric unit")
    chosen = max(options, key=lambda option: option.probability)
    return chosen.copies


_MatthewsProbabilitySetting = collections.namedtuple(
    "_MatthewsProbabilitySetting", ["rbin", "p0", "vmbar", "w", "a", "s"]
)


_MATTHEWS_PROBABILITY_SETTINGS = [
    _MatthewsProbabilitySetting(1.199, 0.085, 2.052, 0.213, 28.38, 0.953),
    _MatthewsProbabilitySetting(1.501, 0.312, 2.102, 0.214, 102.7, 0.807),
    _MatthewsProbabilitySetting(1.650, 0.400, 2.122, 0.220, 187.5, 0.775),
    _MatthewsProbabilitySetting(1.801, 0.503, 2.132, 0.225, 339.3, 0.702),
    _MatthewsProbabilitySetting(1.901, 0.597, 2.140, 0.226, 434.1, 0.648),
    _MatthewsProbabilitySetting(2.001, 0.729, 2.155, 0.231, 540.5, 0.640),
    _MatthewsProbabilitySetting(2.201, 1.052, 2.171, 0.236, 686.2, 0.635),
    _MatthewsProbabilitySetting(2.401, 1.781, 2.182, 0.241, 767.9, 0.589),
    _MatthewsProbabilitySetting(2.601, 2.852, 2.191, 0.242, 835.9, 0.584),
    _MatthewsProbabilitySetting(2.801, 3.386, 2.192, 0.244, 856.9, 0.542),
    _MatthewsProbabilitySetting(3.101, 3.841, 2.205, 0.244, 854.0, 0.500),
    _MatthewsProbabilitySetting(3.501, 4.281, 2.211, 0.244, 849.6, 0.485),
    _MatthewsProbabilitySetting(5.001, 4.592, 2.210, 0.245, 846.7, 0.480),
    _MatthewsProbabilitySetting(5.001, 1.503, 2.256, 0.446, 136.6, 1.180),
    _MatthewsProbabilitySetting(5.001, 0.257, 2.324, 0.327, 47.10, 0.466),
]


def _probability(
    protein_mw: float,
    nucleic_mw: float,
    copies: int,
    asu_volume: float,
    resolution: float,
) -> float:
    total_mw = protein_mw + nucleic_mw
    matt = asu_volume / (total_mw * copies)
    if protein_mw > 0.9 * total_mw:
        for index in range(12):
            if resolution < _MATTHEWS_PROBABILITY_SETTINGS[index].rbin:
                break
    elif nucleic_mw > 0.8 * total_mw:
        index = 13
    else:
        index = 14
    _, p0, vmbar, w, a, s = _MATTHEWS_PROBABILITY_SETTINGS[index]
    z = (matt - vmbar) / w
    return p0 + a * (math.exp(-math.exp(-z) - z * s + 1))
