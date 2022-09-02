import argparse
import collections
import dataclasses
import math
import sys
import gemmi
from ..contents import AsuContents
from ..environ import setup_environ


def main(argument_list=None):
    setup_environ()
    if argument_list is None:
        argument_list = sys.argv[1:]
    parser = argparse.ArgumentParser()
    parser.add_argument("contents", help="Path to contents file")
    parser.add_argument("mtz", help="Path to MTZ file")
    args = parser.parse_args(argument_list)
    contents = AsuContents.from_file(args.contents)
    mtz = gemmi.read_mtz_file(args.mtz)

    cell = mtz.cell
    asu_volume = cell.volume / len(mtz.spacegroup.operations())
    print("## MTZ\n")
    print(
        "Cell        %.3f  %.3f  %.3f  %.2f  %.2f  %.2f"
        % (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    )
    print("Spacegroup ", mtz.spacegroup.hm)
    print("ASU Volume  %.0f" % asu_volume)
    print("Resolution  %.2f - %.2f" % (mtz.resolution_low(), mtz.resolution_high()))
    print("")

    print("## Components\n")
    print("| Description                                         | Volume   | Copies |")
    print("|-----------------------------------------------------|----------|--------|")
    for component in _components(contents):
        description = component.description
        copies = component.copies
        volume = component.volume
        print("| %51s | %8.0f | %6d |" % (description[:51], volume, copies))
    print("|-----------------------------------------------------|----------|--------|")
    print("| %51s | %8.0f |        |" % ("Total", contents.volume()))
    print("")

    options = _options(contents, mtz)
    print("## Multiples\n")
    if len(options) == 0:
        print("Contents are too big to fit into the asymmetric unit")
    else:
        print("| Multiple | Solvent Fraction | Probability |")
        print("|----------|------------------|-------------|")
        for option in options:
            print(
                "| %6d | %16.3f | %11.3f |"
                % (option.multiple, option.solvent, option.probability)
            )
    print("")


@dataclasses.dataclass
class _Component:
    description: str
    copies: int
    volume: float


def _components(contents: AsuContents):
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
            yield _Component(description, polymer.copies, polymer.volume())
    for carb in contents.carbs:
        description = "Carb:"
        for code, count in carb.codes.items():
            description += f" {count}x{code}"
        yield _Component(description, carb.copies, carb.volume())
    for ligand in contents.ligands:
        description = "Ligand: " + ligand.code
        yield _Component(description, ligand.copies, ligand.volume())


@dataclasses.dataclass
class _Option:
    multiple: int
    solvent: float
    probability: float


def _options(contents: AsuContents, mtz: gemmi.Mtz) -> list:
    options = []
    mwp = sum(p.weight() * p.copies for p in contents.proteins)
    mwn = sum(n.weight() * n.copies for n in contents.rnas + contents.dnas)
    asu_volume = mtz.cell.volume / len(mtz.spacegroup.operations())
    resolution = mtz.resolution_high()
    total_probability = 0
    for multiple in range(1, 60):
        solvent = 1 - multiple * contents.volume() / asu_volume
        probability = _probability(mwp, mwn, multiple, asu_volume, resolution)
        if solvent < 0:
            break
        options.append(_Option(multiple, solvent, probability))
        total_probability += probability
    for option in options:
        option.probability /= total_probability
    return options


def _probability(
    protein_mw: float,
    nucleic_mw: float,
    multiple: int,
    asu_volume: float,
    resolution: float,
) -> float:
    total_mw = protein_mw + nucleic_mw
    matt = asu_volume / (total_mw * multiple)
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


if __name__ == "__main__":
    main()
