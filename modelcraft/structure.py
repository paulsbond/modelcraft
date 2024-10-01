from typing import Iterator
import gemmi
from .monlib import atom_ids, in_library, is_protein, is_nucleic


def read_structure(path: str) -> gemmi.Structure:
    structure = gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
    structure.remove_empty_chains()
    structure.remove_hydrogens()
    # TODO: Currently altconfs appear in CIF auth_atom_id after sheetbend
    # TODO: Keep alternative conformations after problem is fixed
    structure.remove_alternative_conformations()
    _patch_names(structure)
    return structure


def consecutive_residues(chain: gemmi.Chain):
    "Iterate through lists of residues with consecutive seqnums (first conformer only)"
    consecutive = []
    last_seqnum = None
    for residue in chain.first_conformer():
        if last_seqnum is None or residue.seqid.num == last_seqnum + 1:
            consecutive.append(residue)
        else:
            yield consecutive
            consecutive = [residue]
        last_seqnum = residue.seqid.num
    if len(consecutive) > 0:
        yield consecutive


def contains_residue(structure: gemmi.Structure, name: str) -> bool:
    return any(residue.name == name for residue in _residues(structure))


def remove_residues(structure: gemmi.Structure, names) -> None:
    for model in structure:
        for chain in model:
            for i, residue in reversed(list(enumerate(chain))):
                if residue.name in names:
                    del chain[i]
    structure.remove_empty_chains()


def remove_non_library_atoms(structure: gemmi.Structure) -> None:
    for residue in _residues(structure):
        if in_library(residue.name):
            for i, atom in reversed(list(enumerate(residue))):
                if atom.name not in atom_ids(residue.name):
                    del residue[i]
    structure.remove_empty_chains()


def remove_non_protein(structure: gemmi.Structure) -> None:
    for model in structure:
        for chain in model:
            for subchain in chain.subchains():
                if subchain.check_polymer_type() != gemmi.PolymerType.PeptideL:
                    for i in reversed(range(len(subchain))):
                        del subchain[i]


def write_mmcif(path: str, structure: gemmi.Structure) -> None:
    groups = gemmi.MmcifOutputGroups(True)
    groups.title_keywords = False
    doc = structure.make_mmcif_document(groups)
    doc.write_file(path, gemmi.cif.Style.Aligned)


class ModelStats:
    def __init__(self, structure: gemmi.Structure):
        self.residues: int = 0
        self.protein: int = 0
        self.nucleic: int = 0
        self.waters: int = 0
        self.dummy_atoms: int = 0

        for residue in _residues(structure):
            if residue.name == "HOH":
                self.waters += 1
            elif residue.name == "DUM":
                self.dummy_atoms += 1
            else:
                self.residues += 1
                if is_protein(residue.name):
                    self.protein += 1
                if is_nucleic(residue.name):
                    self.nucleic += 1

    def __eq__(self, other):
        if isinstance(other, ModelStats):
            return (
                self.residues == other.residues
                and self.waters == other.waters
                and self.dummy_atoms == other.dummy_atoms
            )
        return NotImplemented

    def __ne__(self, other):
        equal = self.__eq__(other)
        return NotImplemented if equal is not NotImplemented else not equal


def _residues(structure: gemmi.Structure) -> Iterator[gemmi.Residue]:
    for model in structure:
        for chain in model:
            for residue in chain:
                yield residue


def _patch_names(structure: gemmi.Structure) -> None:
    residue_patches = {"SUL": "SO4"}
    atom_patches = {("HOH", "O1"): "O"}
    for residue in _residues(structure):
        residue.name = residue.name.strip()
        residue.name = residue_patches.get(residue.name, residue.name)
        for atom in residue:
            atom.name = atom.name.strip()
            atom.name = atom_patches.get((residue.name, atom.name), atom.name)


def are_connected(residue1: gemmi.Residue, residue2: gemmi.Residue) -> bool:
    if (
        is_protein(residue1.name)
        and is_protein(residue2.name)
        and "C" in residue1
        and "N" in residue2
    ):
        for atom1 in residue1["C"]:
            for atom2 in residue2["N"]:
                if atom1.pos.dist(atom2.pos) < 2.5:
                    return True
    if (
        is_nucleic(residue1.name)
        and is_nucleic(residue2.name)
        and "O3'" in residue1
        and "P" in residue2
    ):
        for atom1 in residue1["O3'"]:
            for atom2 in residue2["P"]:
                if atom1.pos.dist(atom2.pos) < 2.5:
                    return True
    return False


def remove_isolated_fragments(chain: gemmi.Chain, max_length: int):
    to_remove = []
    fragment = []
    for i, residue in enumerate(chain):
        if i > 0 and are_connected(chain[i - 1], residue):
            fragment.append(i)
        else:
            if len(fragment) <= max_length:
                to_remove.extend(fragment)
            fragment = [i]
    if len(fragment) <= max_length:
        to_remove.extend(fragment)
    for i in reversed(to_remove):
        del chain[i]
