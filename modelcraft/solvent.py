from collections import namedtuple
from math import exp
import gemmi
from .contents import AsuContents


def solvent_fraction(contents: AsuContents, mtz: gemmi.Mtz) -> float:
    asu_volume = mtz.cell.volume / len(mtz.spacegroup.operations())
    copies = contents.copies or _guess_copies(contents, mtz)
    return 1 - copies * contents.volume() / asu_volume


class CopiesOption:
    def __init__(self, copies: int, solvent: float, probability: float):
        self.copies = copies
        self.solvent = solvent
        self.probability = probability


def copies_options(contents: AsuContents, mtz: gemmi.Mtz) -> list:
    options = []
    mwp = sum(p.weight() * (p.copies or 1) for p in contents.proteins)
    mwn = sum(n.weight() * (n.copies or 1) for n in contents.rnas + contents.dnas)
    asu_volume = mtz.cell.volume / len(mtz.spacegroup.operations())
    resolution = mtz.resolution_high()
    total_probability = 0
    for copies in range(1, 60):
        solvent = 1 - copies * contents.volume() / asu_volume
        probability = _probability(mwp, mwn, copies, asu_volume, resolution)
        if solvent < 0:
            break
        options.append(CopiesOption(copies, solvent, probability))
        total_probability += probability
    for option in options:
        option.probability /= total_probability
    return options


def _guess_copies(contents: AsuContents, mtz: gemmi.Mtz) -> int:
    options = copies_options(contents, mtz)
    if len(options) == 0:
        raise ValueError("Contents are too big to fit into the asymmetric unit")
    chosen = max(options, key=lambda option: option.probability)
    return chosen.copies


def _probability(
    protein_mw: float,
    nucleic_mw: float,
    copies: int,
    asu_volume: float,
    resolution: float,
) -> float:
    mw = protein_mw + nucleic_mw
    matt = asu_volume / (mw * copies)
    Setting = namedtuple("Setting", ["rbin", "p0", "vmbar", "w", "a", "s"])
    settings = [
        Setting(1.199, 0.085, 2.052, 0.213, 28.38, 0.953),
        Setting(1.501, 0.312, 2.102, 0.214, 102.7, 0.807),
        Setting(1.650, 0.400, 2.122, 0.220, 187.5, 0.775),
        Setting(1.801, 0.503, 2.132, 0.225, 339.3, 0.702),
        Setting(1.901, 0.597, 2.140, 0.226, 434.1, 0.648),
        Setting(2.001, 0.729, 2.155, 0.231, 540.5, 0.640),
        Setting(2.201, 1.052, 2.171, 0.236, 686.2, 0.635),
        Setting(2.401, 1.781, 2.182, 0.241, 767.9, 0.589),
        Setting(2.601, 2.852, 2.191, 0.242, 835.9, 0.584),
        Setting(2.801, 3.386, 2.192, 0.244, 856.9, 0.542),
        Setting(3.101, 3.841, 2.205, 0.244, 854.0, 0.500),
        Setting(3.501, 4.281, 2.211, 0.244, 849.6, 0.485),
        Setting(5.001, 4.592, 2.210, 0.245, 846.7, 0.480),
        Setting(5.001, 1.503, 2.256, 0.446, 136.6, 1.180),
        Setting(5.001, 0.257, 2.324, 0.327, 47.10, 0.466),
    ]
    if protein_mw > 0.9 * mw:
        for n in range(12):
            if resolution < settings[n].rbin:
                break
    elif nucleic_mw > 0.8 * mw:
        n = 13
    else:
        n = 14
    _, p0, vmbar, w, a, s = settings[n]
    z = (matt - vmbar) / w
    return p0 + a * (exp(-exp(-z) - z * s + 1))
