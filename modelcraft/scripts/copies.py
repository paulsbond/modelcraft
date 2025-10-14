import argparse
import sys

import gemmi

from ..contents import AsuContents
from ..environ import setup_environ
from ..monlib import MonLib
from ..solvent import copies_options


def main(argument_list=None):
    setup_environ()
    if argument_list is None:
        argument_list = sys.argv[1:]
    parser = argparse.ArgumentParser()
    parser.add_argument("contents", help="Path to contents file")
    parser.add_argument("mtz", help="Path to MTZ file")
    parser.add_argument("libin", help="Path to custom restraint dictionary")
    args = parser.parse_args(argument_list)

    contents = AsuContents.from_file(args.contents)
    mtz = gemmi.read_mtz_file(args.mtz)
    monlib = MonLib(contents.monomer_codes(), args.libin, include_standard=True)

    cell = mtz.cell
    asu_volume = cell.volume / len(mtz.spacegroup.operations())
    print("## MTZ\n")
    print(
        f"Cell        {cell.a:.3f}  {cell.b:.3f}  {cell.c:.3f}"
        f"  {cell.alpha:.2f}  {cell.beta:.2f}  {cell.gamma:.2f}"
    )
    print(f"Spacegroup  {mtz.spacegroup.hm}")
    print(f"ASU Volume  {asu_volume:.0f}")
    print(f"Resolution  {mtz.resolution_low():.2f} - {mtz.resolution_high():.2f}")
    print("")

    print("## Components\n")
    print("| Description                                  | Stoichiometry | Volume   |")
    print("|----------------------------------------------|---------------|----------|")
    for component in contents.components():
        stoichiometry = component.stoichiometry or 1
        assumed = "(assumed)" if component.stoichiometry is None else ""
        volume = component.volume(monlib)
        print(
            f"| {str(component)[:44]:44s} "
            f"| {assumed:9s} {stoichiometry:3d} "
            f"| {volume:8.0f} |"
        )
    print("|----------------------------------------------|---------------|----------|")
    print(f"| {'Total':44s} |               | {contents.volume(monlib):8.0f} |")
    print("")

    options = copies_options(
        contents, cell, mtz.spacegroup, mtz.resolution_high(), monlib
    )
    print("## Copies\n")
    if len(options) == 0:
        print("Contents are too big to fit into the asymmetric unit")
    else:
        print("| Copies | Solvent Fraction | Probability |")
        print("|--------|------------------|-------------|")
        for option in options:
            print(
                f"| {option.copies:6d} "
                f"| {option.solvent:16.3f} "
                f"| {option.probability:11.3f} |"
            )
    print("")


if __name__ == "__main__":
    main()
