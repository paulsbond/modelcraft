import gemmi

_KNOWN_PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "MSE",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "UNK",
    "VAL",
}


def _distance(atom1: gemmi.Atom, atom2: gemmi.Atom) -> float:
    x_2 = (atom1.pos[0] - atom2.pos[0]) ** 2
    y_2 = (atom1.pos[1] - atom2.pos[1]) ** 2
    z_2 = (atom1.pos[2] - atom2.pos[2]) ** 2
    return (x_2 + y_2 + z_2) ** 0.5


def _min_distance(atom_group1: gemmi.AtomGroup, atom_group2: gemmi.AtomGroup) -> float:
    distances = []
    for atom1 in atom_group1:
        for atom2 in atom_group2:
            distances.append(_distance(atom1, atom2))
    return min(distances)


def _is_protein(residue: gemmi.Residue) -> bool:
    return residue.name in _KNOWN_PROTEIN_RESIDUES and all(
        atom in residue for atom in ("N", "CA", "C")
    )


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
        self.sequenced_residues: int = 0
        self.fragments: int = 0
        self.longest_fragment: int = 0
        self.waters: int = 0
        self.dummy_atoms: int = 0

        current_fragment_length = 0
        model = structure[0]
        for chain in model:
            for i, residue in enumerate(chain):
                if residue.name == "HOH":
                    self.waters += 1
                elif residue.name == "DUM":
                    self.dummy_atoms += 1
                elif _is_protein(residue):
                    self.residues += 1
                    current_fragment_length += 1
                    if residue.name != "UNK":
                        self.sequenced_residues += 1
                    if i + 1 < len(chain):
                        next_residue = chain[i + 1]
                        if (
                            _is_protein(next_residue)
                            and _min_distance(residue["C"], next_residue["N"]) < 1.7
                        ):
                            continue
                    self.fragments += 1
                    self.longest_fragment = max(
                        self.longest_fragment, current_fragment_length
                    )
                    current_fragment_length = 0

    def __eq__(self, other):
        if isinstance(other, ModelStats):
            return (
                self.residues == other.residues
                and self.sequenced_residues == other.sequenced_residues
                and self.fragments == other.fragments
                and self.longest_fragment == other.longest_fragment
                and self.waters == other.waters
                and self.dummy_atoms == other.dummy_atoms
            )
        return NotImplemented

    def __ne__(self, other):
        equal = self.__eq__(other)
        return NotImplemented if equal is not NotImplemented else not equal
