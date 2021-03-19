import gemmi


def read_structure(path: str) -> gemmi.Structure:
    if (
        path[-4:] == ".cif"
        or path[-7:] == ".cif.gz"
        or path[-6:] == ".mmcif"
        or path[-9:] == ".mmcif.gz"
    ):
        document = gemmi.cif.read(path)
        block = document[0]  # Assume the first block is the structure
        structure = gemmi.make_structure_from_block(block)
    else:
        structure = gemmi.read_structure(path)
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


def write_mmcif(path: str, structure: gemmi.Structure) -> None:
    structure.make_mmcif_document().write_file(path)


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
