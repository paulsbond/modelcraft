import argparse
import sys
import gemmi
from tabulate import tabulate
from ..environ import setup_environ
from ..reflections import DataItem
from ..validation import validate


def _parse_args(argument_list):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "structure",
        help="Input structure in PDB, mmCIF, mmJSON format",
    )
    parser.add_argument(
        "mtz",
        help="MTZ file from Refmac (with standard output column labels)",
    )
    parser.add_argument(
        "--format",
        default="table",
        choices=["table", "csv"],
        help="Print the results as a human-readable table or a CSV file",
    )
    parser.add_argument(
        "--libin",
        metavar="PATH",
        help="Path to a custom restraint dictionary in CIF format",
    )
    parser.add_argument(
        "--sort",
        action="store_true",
        help="Order the output with the worse scoring residues first",
    )
    parser.add_argument(
        "--model-index",
        type=int,
        default=0,
        metavar="N",
        help="Index of the model to analyse (with 0 being the first model)",
    )
    return parser.parse_args(argument_list or sys.argv[1:])


def main(argument_list=None):
    setup_environ()
    args = _parse_args(argument_list)

    structure = gemmi.read_structure(args.structure, format=gemmi.CoorFormat.Detect)
    mtz = gemmi.read_mtz_file(args.mtz)
    fphi_best = DataItem(mtz, "FWT,PHWT")
    fphi_diff = DataItem(mtz, "DELFWT,PHDELWT")
    fphi_calc = DataItem(mtz, "FC_ALL_LS,PHIC_ALL_LS")

    metrics = validate(
        structure, fphi_best, fphi_diff, fphi_calc, args.model_index, args.libin
    )
    if args.sort:
        metrics.sort_values("Score", ascending=True, inplace=True)

    if args.format == "table":
        metrics["Sig"] = metrics["Score"].apply(lambda x: "+" * min(5, -int(x)))
        print(tabulate(metrics, headers="keys", showindex=False, floatfmt=".1f"))
    else:
        print(metrics.to_csv(index=False))


if __name__ == "__main__":
    main()
