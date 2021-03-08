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
    parser.add_argument("contents")
    parser.add_argument("mtz")
    args = parser.parse_args(argument_list)
    contents = AsuContents(args.contents)
    mtz = gemmi.read_mtz_file(args.mtz)

    cell = mtz.cell
    asu_volume = cell.volume / len(mtz.spacegroup.operations())
    print("## Cell\n")
    print(
        "Cell        %.3f  %.3f  %.3f  %.2f  %.2f  %.2f"
        % (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    )
    print("Spacegroup ", mtz.spacegroup.hm)
    print("ASU Volume  %.0f" % asu_volume)
    print("")

    print("## Components\n")
    print("| Type    | Stoichiometry | Volume   |")
    print("|---------|---------------|----------|")
    for kind, items in (
        ("Protein", contents.proteins),
        ("RNA", contents.rnas),
        ("DNA", contents.dnas),
        ("Carb", contents.carbs),
        ("Ligand", contents.ligands),
    ):
        for item in items:
            copies = item.copies or 1
            assumed = "" if item.copies else "(assumed)"
            volume = item.volume() * copies
            print("| %7s | %9s %3d | %8.0f |" % (kind, assumed, copies, volume))
    print("|---------|---------------|----------|")
    print("|         |         Total | %8.0f |" % contents.volume())
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
