#!/usr/bin/python3

import gemmi

_known_protein_residues = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "UNK", "VAL",
}


def distance(atom1, atom2):
    x2 = (atom1.pos[0] - atom2.pos[0]) ** 2
    y2 = (atom1.pos[1] - atom2.pos[1]) ** 2
    z2 = (atom1.pos[2] - atom2.pos[2]) ** 2
    return (x2 + y2 + z2) ** 0.5


def is_protein(residue):
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
            if is_protein(residue):
                stats["residues_built"] += 1
                current_fragment += 1
                if residue.name != "UNK":
                    stats["residues_sequenced"] += 1
                if i + 1 < len(chain):
                    next_residue = chain[i + 1]
                    if is_protein(next_residue) and distance(residue["C"], next_residue["N"]) < 1.7:
                        continue
                stats["fragments_built"] += 1
                if current_fragment > stats["longest_fragment"]:
                    stats["longest_fragment"] = current_fragment
                current_fragment = 0
    return stats


def fo_columns(mtz):
    def differ_by_sig(label1, label2):
        label1 = label1.lower().replace("sig", "")
        label2 = label2.lower().replace("sig", "")
        return label1 == label2
    for f in mtz.columns_with_type("F"):
        for sig in mtz.columns_with_type("Q"):
            if differ_by_sig(f.label, sig.label):
                yield "%s,%s" % (f.label, sig.label)


def free_columns(mtz):
    for column in mtz.columns_with_type("I"):
        if "free" in column.label.lower():
            yield column.label


def hl_columns(mtz):
    def belong_to_the_same_group(label1, label2):
        if len(label1) != len(label2):
            return False
        differences = 0
        for i in range(len(label1)):
            if label1[i] != label2[i]:
                differences += 1
        return differences == 1
    groups = []
    for column in mtz.columns_with_type("A"):
        grouped = False
        for group in groups:
            if belong_to_the_same_group(column.label, group[0]):
                grouped = True
                group.append(column.label)
                if len(group) == 4:
                    yield ",".join(sorted(group))
        if not grouped:
            groups.append([column.label])


def phifom_columns(mtz):

    def is_phi(column):
        if column.type != "P":
            return False
        label = column.label.lower()
        return "ph" in label and "del" not in label

    def is_fom(column):
        if column.type != "W":
            return False
        label = column.label.lower()
        return "fom" in label

    phis = [column for column in mtz.columns if is_phi(column)]
    foms = [column for column in mtz.columns if is_fom(column)]
    for phi in phis:
        for fom in foms:
            yield "%s,%s" % (phi.label, fom.label)
