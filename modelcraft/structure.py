import gemmi
from .monlib import atom_ids, in_library


def read_structure(path: str) -> gemmi.Structure:
    structure = gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
    structure.remove_empty_chains()
    structure.remove_hydrogens()
    # TODO: Currently altconfs appear in CIF auth_atom_id after sheetbend
    # TODO: Keep alternative conformations after problem is fixed
    structure.remove_alternative_conformations()
    return structure


def contains_residue(structure: gemmi.Structure, name: str) -> bool:
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.name == name:
                    return True
    return False


def remove_residues(structure: gemmi.Structure, names) -> None:
    for model in structure:
        for chain in model:
            for i, residue in reversed(list(enumerate(chain))):
                if residue.name in names:
                    del chain[i]
    structure.remove_empty_chains()


def remove_non_library_atoms(structure: gemmi.Structure) -> None:
    for model in structure:
        for chain in model:
            for residue in chain:
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
    doc.write_file(path)


class ModelStats:
    def __init__(self, structure: gemmi.Structure):
        self.residues: int = 0
        self.waters: int = 0
        self.dummy_atoms: int = 0

        model = structure[0]
        for chain in model:
            for residue in chain:
                if residue.name == "HOH":
                    self.waters += 1
                elif residue.name == "DUM":
                    self.dummy_atoms += 1
                else:
                    self.residues += 1

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
