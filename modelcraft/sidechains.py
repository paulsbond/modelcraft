from tempfile import NamedTemporaryFile
import chapi
import gemmi
from .reflections import DataItem, write_mtz
from .structure import write_mmcif


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
    return expected.issubset(built)


def build_missing_side_chains(
    structure: gemmi.Structure, fphi_best: DataItem
) -> gemmi.Structure:
    "Build missing side chains in a structure using Coot's CHAPI interface."
    mc = chapi.molecules_container_t(True)
    with NamedTemporaryFile(suffix=".cif") as temp_cif:
        write_mmcif(temp_cif.name, structure)
        imol = mc.read_coordinates(temp_cif.name)
    with NamedTemporaryFile(suffix=".mtz") as temp_mtz:
        write_mtz(temp_mtz.name, [fphi_best], ["FWT,PHWT"])
        imap = mc.read_mtz(temp_mtz.name, "FWT", "PHWT", "", False, False)
    mc.set_imol_refinement_map(imap)
    for chain in structure[0]:
        for residue in chain:
            if not has_full_side_chain(residue):
                num = residue.seqid.num
                icode = residue.seqid.icode
                mc.fill_partial_residue(imol, chain.name, num, icode)
                mc.auto_fit_rotamer(imol, chain.name, num, icode, "", imap)
    with NamedTemporaryFile(suffix=".cif") as temp_cif:
        mc.write_coordinates(imol, temp_cif.name)
        return gemmi.read_structure(temp_cif.name)
