import enum
from .monomers import Monomers


class PolymerType(enum.Enum):
    PROTEIN = "PROTEIN"
    RNA = "RNA"
    DNA = "DNA"

    @classmethod
    def parse(cls, string: str) -> "PolymerType":
        if string.lower() in ("protein", "polypeptide(l)"):
            return cls.PROTEIN
        if string.lower() in ("rna", "polyribonucleotide"):
            return cls.RNA
        if string.lower() in ("dna", "polydeoxyribonucleotide"):
            return cls.DNA
        raise ValueError(f"Unknown polymer type: '{string}'")

    @classmethod
    def from_sequence(cls, sequence: str) -> "PolymerType":
        codes = set(sequence)
        if "U" in codes:
            return cls.RNA
        if codes & set("DEFHIKLMNPQRSVWY"):
            return cls.PROTEIN
        if codes == {"A"}:
            return cls.PROTEIN
        if codes == {"G"}:
            return cls.PROTEIN
        if "T" in codes:
            return cls.DNA
        return cls.RNA


class Polymer:
    def __init__(
        self,
        sequence: str,
        start: int = None,
        copies: int = None,
        polymer_type: PolymerType = None,
        modifications: list = None,
        monomers: Monomers = None,
    ):
        self.sequence = sequence.upper()
        self.start = start or 1
        self.copies = copies
        if polymer_type is None:
            self.type = PolymerType.from_sequence(self.sequence)
        else:
            self.type = polymer_type
        self.modifications = modifications or []
        self.monomers = monomers or Monomers()

    def __eq__(self, other) -> bool:
        if isinstance(other, Polymer):
            return (
                self.sequence == other.sequence
                and self.type == other.type
                and self.modifications == other.modifications
            )
        return NotImplemented

    @classmethod
    def from_json(cls, component: dict) -> "Polymer":
        return cls(
            sequence=component["sequence"],
            start=component.get("start"),
            copies=component.get("copies"),
            modifications=component.get("modifications"),
        )

    @classmethod
    def from_pdbe_molecule_dict(cls, mol: dict) -> "Polymer":
        return cls(
            sequence=mol["sequence"],
            copies=mol["number_of_copies"],
            polymer_type=PolymerType.parse(mol["molecule_type"]),
            modifications=_modifications_in_pdbe_molecule_dict(mol),
        )

    @classmethod
    def from_sequence_file(cls, path: str, polymer_type: PolymerType = None):
        sequence = ""
        with open(path) as stream:
            for line in stream:
                if line[0] == ">":
                    if len(sequence) > 0:
                        yield cls(sequence=sequence, polymer_type=polymer_type)
                    sequence = ""
                elif line[0] != ";":
                    sequence += "".join(c for c in line if c.isalpha())
        if len(sequence) > 0:
            yield cls(sequence=sequence, polymer_type=polymer_type)

    def to_json(self) -> dict:
        return {
            "sequence": self.sequence,
            "start": self.start,
            "copies": self.copies,
            "modifications": self.modifications,
        }

    def residue_codes(self, modified: bool = True) -> list:
        codes = [code1_to_code3(code1, self.type) for code1 in self.sequence]
        if modified:
            for mod in self.modifications:
                source, code = mod.split("->")
                if source[0].isdigit():
                    index = int(source) - self.start
                    codes[index] = code
                else:
                    for index, code1 in enumerate(self.sequence):
                        if code1 == source:
                            codes[index] = code
        return codes

    def weight(self, modified=False) -> float:
        codes = self.residue_codes(modified=modified)
        total = sum(self.monomers.weight(code) for code in codes)
        total -= self.monomers.weight("HOH") * (len(codes) - 1)
        return total

    def volume(self) -> float:
        density = 1.35 if self.type == PolymerType.PROTEIN else 2.0
        return self.weight() / (density * 0.602214)

    def is_selenomet(self) -> bool:
        return "M->MSE" in self.modifications


def _modifications_in_pdbe_molecule_dict(mol: dict) -> list:
    indices = {}
    for index, mod in mol["pdb_sequence_indices_with_multiple_residues"].items():
        code1 = mod["one_letter_code"]
        code3 = mod["three_letter_code"]
        if code3 not in ("DA", "DC", "DG", "DT"):
            key = code1, code3
            indices.setdefault(key, []).append(index)
    modifications = []
    for key in indices:
        code1, code3 = key
        total = mol["sequence"].count(code1)
        if code1 == "M" and mol["sequence"][0] == "M":
            total -= 1
        if len(indices[key]) >= total:
            modifications.append(f"{code1}->{code3}")
        else:
            modifications.extend(f"{index}->{code3}" for index in indices[key])
    return modifications


def code1_to_code3(code1: str, polymer_type: PolymerType) -> str:
    if polymer_type == PolymerType.PROTEIN:
        return {
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
            "Y": "TYR",
            "Z": "GLX",
        }.get(code1) or "UNK"
    if polymer_type == PolymerType.RNA:
        return {
            "A": "A",
            "C": "C",
            "G": "G",
            "I": "I",
            "U": "U",
        }.get(code1) or "N"
    return {
        "A": "DA",
        "C": "DC",
        "G": "DG",
        "I": "DI",
        "T": "DT",
        "U": "DU",
    }.get(code1) or "DN"
