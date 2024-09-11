import argparse
import sys
import gemmi
from ..environ import setup_environ
from ..rszd import __doc__, per_residue_rszd
from ..reflections import DataItem


def main(argument_list=None):
    setup_environ()
    if argument_list is None:
        argument_list = sys.argv[1:]
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
        help="MTZ file with coefficients for the mFo-DFc difference map",
    )
    parser.add_argument(
        "--fphi-diff",
        default="DELFWT,PHDELWT",
        help="Column labels for F,PHI of the mFo-DFc difference map",
    )
    parser.add_argument(
        "--format",
        default="table",
        choices=["table", "csv"],
        help="Print the results as a human-readable table or a CSV file",
    )
    parser.add_argument(
        "--model-index",
        type=int,
        default=0,
        help="Index of the model to analyse (with 0 being the first model)",
    )
    args = parser.parse_args(argument_list)

    structure = gemmi.read_structure(args.structure, format=gemmi.CoorFormat.Detect)
    mtz = gemmi.read_mtz_file(args.mtz)
    fphi_diff = DataItem(mtz, args.fphi_diff)

    rszds = per_residue_rszd(structure, fphi_diff, args.model_index)

    if args.format == "table":
        print("| Chain | Number | ICode |  RSZD | Sig |")
        print("|-------|--------|-------|-------|-----|")
    else:
        print("chain,number,icode,rszd")
    for chain in structure[args.model_index]:
        for residue in chain:
            num = residue.seqid.num
            icode = residue.seqid.icode
            rszd = rszds[(chain.name, str(residue.seqid))]
            sig = "+" * min(int(rszd), 3)
            if args.format == "table":
                print(
                    f"| {chain.name:5s} | {num:6d} | {icode:5s} | {rszd:5.3f} | {sig:3s} |"
                )
            else:
                print(f"{chain.name},{num},{icode},{rszd:.3f}")


if __name__ == "__main__":
    main()
