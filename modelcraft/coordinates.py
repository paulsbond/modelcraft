import gemmi
import os

_known_protein_residues = {
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


def _distance(atom1, atom2):
    x2 = (atom1.pos[0] - atom2.pos[0]) ** 2
    y2 = (atom1.pos[1] - atom2.pos[1]) ** 2
    z2 = (atom1.pos[2] - atom2.pos[2]) ** 2
    return (x2 + y2 + z2) ** 0.5


def _min_distance(atom_group1, atom_group2):
    distances = []
    for atom1 in atom_group1:
        for atom2 in atom_group2:
            distances.append(_distance(atom1, atom2))
    return min(distances)


def _is_protein(residue):
    if residue.name not in _known_protein_residues:
        return False
    return all(atom in residue for atom in ("N", "CA", "C"))


class CoordinateFile:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.residues = 0
        self.sequenced_residues = 0
        self.fragments = 0
        self.longest_fragment = 0
        self.waters = 0
        self.dummys = 0
        self.exists = os.path.exists(path)
        if self.exists:
            self.get_stats()

    def get_stats(self):
        current_fragment_length = 0
        structure = gemmi.read_structure(self.path)
        for chain in structure[0]:
            for i, residue in enumerate(chain):
                if residue.name == "HOH":
                    self.waters += 1
                elif residue.name == "DUM":
                    self.dummys += 1
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
