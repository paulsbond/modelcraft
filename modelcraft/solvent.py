import collections
import dataclasses
import math
import gemmi
from .contents import AsuContents
from .monlib import MonLib


def solvent_fraction(
    contents: AsuContents,
    cell: gemmi.UnitCell,
    spacegroup: gemmi.SpaceGroup,
    resolution: float,
) -> float:
    monlib = MonLib(contents.monomer_codes(), include_standard=True)
    asu_volume = cell.volume / len(spacegroup.operations())
    copies = contents.copies or _guess_copies(contents, cell, spacegroup, resolution)
    return 1 - copies * contents.volume(monlib) / asu_volume


@dataclasses.dataclass
class CopiesOption:
    copies: int
    solvent: float
    probability: float


def copies_options(
    contents: AsuContents,
    cell: gemmi.UnitCell,
    spacegroup: gemmi.SpaceGroup,
    resolution: float,
    monlib: MonLib,
) -> list:
    options = []
    nucleic_acids = contents.rnas + contents.dnas
    mwp = sum(p.weight(monlib) * (p.stoichiometry or 1) for p in contents.proteins)
    mwn = sum(n.weight(monlib) * (n.stoichiometry or 1) for n in nucleic_acids)
    asu_volume = cell.volume / len(spacegroup.operations())
    contents_volume = contents.volume(monlib)
    total_probability = 0
    for copies in range(1, 60):
        solvent = 1 - copies * contents_volume / asu_volume
        probability = _matthews_probability(mwp, mwn, copies, asu_volume, resolution)
        if solvent < 0:
            break
        options.append(CopiesOption(copies, solvent, probability))
        total_probability += probability
    for option in options:
        option.probability /= total_probability
    return options


def _guess_copies(
    contents: AsuContents,
    cell: gemmi.UnitCell,
    spacegroup: gemmi.SpaceGroup,
    resolution: float,
) -> int:
    options = copies_options(contents, cell, spacegroup, resolution)
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


def _matthews_probability(
    protein_mw: float,
    nucleic_mw: float,
    copies: int,
    asu_volume: float,
    resolution: float,
) -> float:
    total_mw = (protein_mw + nucleic_mw) * copies
    matthews = asu_volume / total_mw
    if protein_mw > 0.9 * total_mw:
        for index in range(12):
            if resolution < _MATTHEWS_PROBABILITY_SETTINGS[index].rbin:
                break
    elif nucleic_mw > 0.8 * total_mw:
        index = 13
    else:
        index = 14
    _, p0, vmbar, w, a, s = _MATTHEWS_PROBABILITY_SETTINGS[index]
    z = (matthews - vmbar) / w
    return p0 + a * (math.exp(-math.exp(-z) - z * s + 1))
