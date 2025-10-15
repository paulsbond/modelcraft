"Add missing side chains to a protein model"

import argparse
import sys

import chapi
import gemmi

from ..environ import setup_environ


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
        help=(
            "MTZ file from Refmac with standard output column labels "
            "(the output MTZ from ModelCraft meets this requirement)."
        ),
    )
    parser.add_argument(
        "output",
        help="Path to write the output structure with added side chains",
    )
    parser.add_argument(
        "--model-index",
        type=int,
        default=0,
        metavar="N",
        help="Index of the model to analyse (with 0 being the first model)",
    )
    return parser.parse_args(argument_list or sys.argv[1:])


SIDE_CHAIN_ATOMS = {
    "ARG": {"CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"CG", "OD1", "ND2"},
    "ASP": {"CG", "OD1", "OD2"},
    "CYS": {"SG"},
    "GLN": {"CG", "CD", "OE1", "NE2"},
    "GLU": {"CG", "CD", "OE1", "OE2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"CG1", "CG2", "CD1"},
    "LEU": {"CG", "CD1", "CD2"},
    "LYS": {"CG", "CD", "CE", "NZ"},
    "MET": {"CG", "SD", "CE"},
    "MSE": {"CG", "SE", "CE"},
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"CG", "CD"},
    "SER": {"OG"},
    "THR": {"OG1", "CG2"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": {"CG1", "CG2"},
}


def has_full_side_chain(residue: gemmi.Residue) -> bool:
    "Check if a residue has all side chain atoms from gamma onwards."
    expected = SIDE_CHAIN_ATOMS.get(residue.name, set())
    built = {atom.name for atom in residue}
    return built > expected


def any_missing_side_chains(structure: gemmi.Structure) -> bool:
    "Check if any residue in a structure has missing side chain atoms."
    for chain in structure[0]:
        for residue in chain:
            if not has_full_side_chain(residue):
                return True
    return False


def main(argument_list=None):
    setup_environ()
    args = _parse_args(argument_list)
    structure = gemmi.read_structure(args.structure)
    if not any_missing_side_chains(structure):
        print("No missing side chains detected, no action taken")
        return
    mc = chapi.molecules_container_t(True)
    mc.set_use_gemmi(False)
    imol = mc.read_coordinates(args.structure)
    imap = mc.read_mtz(args.mtz, "FWT", "PHWT", "", False, False)
    mc.set_imol_refinement_map(imap)
    mc.set_use_torsion_restraints(True)
    mc.set_use_rama_plot_restraints(True)
    for chain in structure[0]:
        for residue in chain:
            if not has_full_side_chain(residue):
                num = residue.seqid.num
                icode = residue.seqid.icode
                icode = "" if icode == " " else icode
                mc.refine_residues(imol, chain.name, num, icode, "", "TRIPLE", 1000)
                mc.fill_partial_residue(imol, chain.name, num, icode)
                mc.auto_fit_rotamer(imol, chain.name, num, icode, "", imap)
                mc.refine_residues(imol, chain.name, num, icode, "", "TRIPLE", 1000)
    mc.write_coordinates(imol, args.output)


if __name__ == "__main__":
    main()
