from collections import namedtuple


Residue = namedtuple("Residue", ["code1", "code3", "name"])

PROTEIN = [
    Residue("A", "ALA", "Alanine"),
    Residue("C", "CYS", "Cysteine"),
    Residue("D", "ASP", "Aspartic Acid"),
    Residue("E", "GLU", "Glutamic Acid"),
    Residue("F", "PHE", "Phenylalanine"),
    Residue("G", "GLY", "Glycine"),
    Residue("H", "HIS", "Histidine"),
    Residue("I", "ILE", "Isoleucine"),
    Residue("K", "LYS", "Lysine"),
    Residue("L", "LEU", "Leucine"),
    Residue("M", "MET", "Methionine"),
    Residue("N", "ASN", "Asparagine"),
    Residue("P", "PRO", "Proline"),
    Residue("Q", "GLN", "Glutamine"),
    Residue("R", "ARG", "Arginine"),
    Residue("S", "SER", "Serine"),
    Residue("T", "THR", "Threonine"),
    Residue("V", "VAL", "Valine"),
    Residue("W", "TRP", "Tryptophan"),
    Residue("Y", "TYR", "Tyrosine"),
]

RNA = [
    Residue("A", "A", "Adenosine monophosphate"),
    Residue("C", "C", "Cytidine monophosphate"),
    Residue("G", "G", "Guanosine monophosphate"),
    Residue("U", "U", "Uridine monophosphate"),
]

DNA = [
    Residue("A", "DA", "Deoxyadenosine monophosphate"),
    Residue("C", "DC", "Deoxycytidine monophosphate"),
    Residue("G", "DG", "Deoxyguanosine monophosphate"),
    Residue("T", "DT", "Thymidine monophosphate"),
]
