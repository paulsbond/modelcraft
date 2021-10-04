import gemmi


def read_structure(path: str) -> gemmi.Structure:
    structure = None
    errors = {}
    for format_ in (
        gemmi.CoorFormat.Pdb,
        gemmi.CoorFormat.Mmcif,
        gemmi.CoorFormat.Mmjson,
    ):
        try:
            structure = gemmi.read_structure(path, format=format_)
            break
        except (RuntimeError, ValueError) as error:
            errors[format_] = error
    if structure is None:
        message = "Unable to read structure"
        for format_, error in errors.items():
            message += f"\nError for {format_}:\n{error}"
        raise RuntimeError(message)
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
