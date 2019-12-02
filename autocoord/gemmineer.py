import gemmi


_known_protein_residues = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "MSE", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "UNK", "VAL",
}


def _distance(atom1, atom2):
    x2 = (atom1.pos[0] - atom2.pos[0]) ** 2
    y2 = (atom1.pos[1] - atom2.pos[1]) ** 2
    z2 = (atom1.pos[2] - atom2.pos[2]) ** 2
    return (x2 + y2 + z2) ** 0.5


def _is_protein(residue):
    if residue.name not in _known_protein_residues:
        return False
    return all(atom in residue for atom in ("N", "CA", "C"))


def model_stats(xyzin):
    stats = {
        "residues_built": 0,
        "residues_sequenced": 0,
        "fragments_built": 0,
        "longest_fragment": 0,
    }
    current_fragment = 0
    structure = gemmi.read_structure(xyzin)
    for chain in structure[0]:
        for i, residue in enumerate(chain):
            if _is_protein(residue):
                stats["residues_built"] += 1
                current_fragment += 1
                if residue.name != "UNK":
                    stats["residues_sequenced"] += 1
                if i + 1 < len(chain):
                    next_residue = chain[i + 1]
                    if (_is_protein(next_residue) and
                            _distance(residue["C"], next_residue["N"]) < 1.7):
                        continue
                stats["fragments_built"] += 1
                if current_fragment > stats["longest_fragment"]:
                    stats["longest_fragment"] = current_fragment
                current_fragment = 0
    return stats
