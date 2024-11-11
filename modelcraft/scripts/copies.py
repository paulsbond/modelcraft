import argparse
import sys
import gemmi
from ..contents import AsuContents
from ..environ import setup_environ
from ..solvent import copies_options


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
    print("| Description                                  | Stoichiometry | Volume   |")
    print("|----------------------------------------------|---------------|----------|")
    for component in contents.components():
        description = component.description()
        stoichiometry = component.stoichiometry or 1
        assumed = "(assumed)" if component.stoichiometry is None else ""
        volume = component.volume
        print(
            "| %44s | %9s %3d | %8.0f |"
            % (description[:44], assumed, stoichiometry, volume)
        )
    print("|----------------------------------------------|---------------|----------|")
    print("| %44s |               | %8.0f |" % ("Total", contents.volume()))
    print("")

    options = copies_options(contents, mtz)
    print("## Copies\n")
    if len(options) == 0:
        print("Contents are too big to fit into the asymmetric unit")
    else:
        print("| Copies | Solvent Fraction | Probability |")
        print("|--------|------------------|-------------|")
        for option in options:
            print(
                "| %6d | %16.3f | %11.3f |"
                % (option.copies, option.solvent, option.probability)
            )
    print("")


if __name__ == "__main__":
    main()
