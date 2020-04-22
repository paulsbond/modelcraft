from dataclasses import dataclass
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


@dataclass
class ModelStats:
    residues: int = 0
    sequenced_residues: int = 0
    fragments: int = 0
    longest_fragment: int = 0
    waters: int = 0
    dummy_atoms: int = 0


def model_stats(model: gemmi.Model) -> ModelStats:
    stats = ModelStats
    current_fragment_length = 0
    for chain in model:
        for i, residue in enumerate(chain):
            if residue.name == "HOH":
                stats.waters += 1
            elif residue.name == "DUM":
                stats.dummy_atoms += 1
            elif _is_protein(residue):
                stats.residues += 1
                current_fragment_length += 1
                if residue.name != "UNK":
                    stats.sequenced_residues += 1
                if i + 1 < len(chain):
                    next_residue = chain[i + 1]
                    if (
                        _is_protein(next_residue)
                        and _min_distance(residue["C"], next_residue["N"]) < 1.7
                    ):
                        continue
                stats.fragments += 1
                stats.longest_fragment = max(
                    stats.longest_fragment, current_fragment_length
                )
                current_fragment_length = 0
    return stats
