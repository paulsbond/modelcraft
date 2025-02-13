PROTEIN_CODES = {
    "A": "ALA",
    "B": "ASX",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "O": "PYL",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "U": "SEC",
    "V": "VAL",
    "W": "TRP",
    "X": "UNK",
    "Y": "TYR",
    "Z": "GLX",
}


DNA_CODES = {
    "A": "DA",
    "C": "DC",
    "G": "DG",
    "I": "DI",
    "T": "DT",
    "U": "DU",
    "X": "DN",
}

RNA_CODES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "I": "I",
    "U": "U",
    "X": "N",
}


class PolymerType:
    PROTEIN = "PolymerType"
    DNA = "PolymerType"
    RNA = "PolymerType"

    def __init__(self, name: str, codes: dict[str, str]):
        self.name = name
        self.codes = codes

    def parse(self, sequence: str) -> list[str]:
        names = []
        for code in sequence:
            if code not in self.codes:
                raise ValueError(f"'{code}' is not a valid {self.name} code")
            names.append(self.codes[code])
        return names

    @classmethod
    def guess(cls, sequence: str):
        codes = set(sequence)
        if codes & set("DEFHKLMNPQRSVWY") or codes in ({"A"}, {"G"}):
            return cls.PROTEIN
        if "U" in codes:
            return cls.RNA
        if "T" in codes:
            return cls.DNA
        return cls.RNA


PolymerType.PROTEIN = PolymerType("Protein", PROTEIN_CODES)
PolymerType.DNA = PolymerType("DNA", DNA_CODES)
PolymerType.RNA = PolymerType("RNA", RNA_CODES)


PIR_CODES = {"D1", "DC", "DL", "F1", "N1", "N3", "P1", "RC", "RL", "XX"}


def sequences_in_file(contents: str) -> list:
    sequence = ""
    sequences = []
    skip_line = False
    skip_lines = False
    lines = contents.splitlines(keepends=False)
    for line in lines:
        if skip_line:
            skip_line = False
            continue
        if line[:1] == ">":
            if len(sequence) > 0:
                sequences.append(sequence)
            sequence = ""
            if line[1:3] in PIR_CODES and line[3] == ";":
                skip_line = True
            skip_lines = False
        elif line[:1] != ";" and not skip_lines:
            sequence += "".join(c for c in line if c.isalpha())
            if line[-1:] == "*":
                skip_lines = True
    if len(sequence) > 0:
        sequences.append(sequence)
    return sequences
