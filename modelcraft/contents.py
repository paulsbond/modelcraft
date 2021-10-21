import enum
import json

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

RNA_CODES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "I": "I",
    "U": "U",
    "X": "N",
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

PIR_CODES = {"D1", "DC", "DL", "F1", "N1", "N3", "P1", "RC", "RL", "XX"}


class PolymerType(enum.Enum):
    PROTEIN = "PROTEIN"
    RNA = "RNA"
    DNA = "DNA"


class Polymer:
    def __init__(
        self,
        sequence: str,
        stoichiometry: int = None,
        polymer_type: PolymerType = None,
        modifications: list = None,
    ):
        self.sequence = sequence.upper()
        self.stoichiometry = stoichiometry
        self.type = polymer_type or guess_sequence_type(self.sequence)
        self.modifications = modifications or []

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
            stoichiometry=component.get("stoichiometry"),
            modifications=component.get("modifications"),
        )

    @classmethod
    def from_sequence_file(cls, path: str, polymer_type: PolymerType = None):
        with open(path) as stream:
            contents = stream.read()
            for sequence in sequences_in_file(contents=contents):
                yield cls(sequence=sequence, polymer_type=polymer_type)

    def to_json(self) -> dict:
        return {
            "sequence": self.sequence,
            "stoichiometry": self.stoichiometry,
            "modifications": self.modifications,
        }

    def residue_codes(self, modified: bool = True) -> list:
        codes = [code1_to_code3(code1, self.type) for code1 in self.sequence]
        if modified:
            for mod in self.modifications:
                source, code = mod.split("->")
                if source[0].isdigit():
                    index = int(source) - 1
                    codes[index] = code
                else:
                    for index, code1 in enumerate(self.sequence):
                        if code1 == source:
                            codes[index] = code
        return codes

    def is_selenomet(self) -> bool:
        return "M->MSE" in self.modifications


class Carb:
    def __init__(self, codes: dict, stoichiometry: int = None):
        self.codes = codes
        self.stoichiometry = stoichiometry

    def __eq__(self, other) -> bool:
        if isinstance(other, Carb):
            return self.codes == other.codes
        return NotImplemented

    @classmethod
    def from_json(cls, component: dict) -> "Carb":
        return cls(
            codes=component["codes"],
            stoichiometry=component.get("stoichiometry"),
        )

    def to_json(self) -> dict:
        return {"codes": self.codes, "stoichiometry": self.stoichiometry}


class Ligand:
    def __init__(self, code: str, stoichiometry: int = None):
        self.code = code
        self.stoichiometry = stoichiometry

    def __eq__(self, other) -> bool:
        if isinstance(other, Ligand):
            return self.code == other.code
        return NotImplemented

    @classmethod
    def from_json(cls, component: dict) -> "Ligand":
        return cls(
            code=component["code"],
            stoichiometry=component.get("stoichiometry"),
        )

    def to_json(self) -> dict:
        return {"code": self.code, "stoichiometry": self.stoichiometry}


class AsuContents:
    def __init__(
        self,
        copies: int = None,
        proteins: list = None,
        rnas: list = None,
        dnas: list = None,
        carbs: list = None,
        ligands: list = None,
        buffers: list = None,
        smiles: dict = None,
    ):
        self.copies = copies
        self.proteins = proteins or []
        self.rnas = rnas or []
        self.dnas = dnas or []
        self.carbs = carbs or []
        self.ligands = ligands or []
        self.buffers = buffers or []
        self.smiles = smiles or {}

    @classmethod
    def from_file(cls, path: str) -> "AsuContents":
        if path[-5:] == ".json":
            return cls.from_json_file(path)
        return cls.from_sequence_file(path)

    @classmethod
    def from_json_file(cls, path: str) -> "AsuContents":
        contents = cls()
        with open(path) as stream:
            contents_json = json.load(stream)
        contents.copies = contents_json.get("copies")
        for obj in contents_json.get("proteins") or []:
            polymer = Polymer.from_json(obj)
            polymer.type = PolymerType.PROTEIN
            contents.proteins.append(polymer)
        for obj in contents_json.get("rnas") or []:
            polymer = Polymer.from_json(obj)
            polymer.type = PolymerType.RNA
            contents.rnas.append(polymer)
        for obj in contents_json.get("dnas") or []:
            polymer = Polymer.from_json(obj)
            polymer.type = PolymerType.DNA
            contents.dnas.append(polymer)
        for obj in contents_json.get("carbs") or []:
            carb = Carb.from_json(obj)
            contents.carbs.append(carb)
        for obj in contents_json.get("ligands") or []:
            ligand = Ligand.from_json(obj)
            contents.ligands.append(ligand)
        contents.buffers = contents_json.get("buffers") or []
        contents.smiles = contents_json.get("smiles") or []
        return contents

    @classmethod
    def from_sequence_file(
        cls, path: str, polymer_type: PolymerType = None
    ) -> "AsuContents":
        contents = cls()
        for polymer in Polymer.from_sequence_file(path, polymer_type):
            contents.add_polymer(polymer)
        return contents

    def add_polymer(self, polymer: Polymer) -> None:
        if polymer.type == PolymerType.PROTEIN:
            self.proteins.append(polymer)
        if polymer.type == PolymerType.RNA:
            self.rnas.append(polymer)
        if polymer.type == PolymerType.DNA:
            self.dnas.append(polymer)

    def monomer_codes(self) -> set:
        codes = set()
        for polymer in self.proteins + self.rnas + self.dnas:
            codes.update(set(polymer.residue_codes(modified=True)))
        for carb in self.carbs:
            codes.update(set(carb.codes.keys()))
        for ligand in self.ligands:
            codes.add(ligand.code)
        codes.update(set(self.buffers))
        return codes

    def is_selenomet(self) -> bool:
        return len(self.proteins) > 0 and all(p.is_selenomet() for p in self.proteins)

    def to_json(self) -> list:
        return {
            "copies": self.copies,
            "proteins": [protein.to_json() for protein in self.proteins],
            "rnas": [rna.to_json() for rna in self.rnas],
            "dnas": [dna.to_json() for dna in self.dnas],
            "carbs": [carb.to_json() for carb in self.carbs],
            "ligands": [ligand.to_json() for ligand in self.ligands],
            "buffers": self.buffers,
            "smiles": self.smiles,
        }

    def write_json_file(self, path: str) -> None:
        with open(path, "w") as stream:
            json.dump(self.to_json(), stream, indent=2)

    def write_sequence_file(
        self, path: str, types: list = None, line_length: int = 60
    ) -> None:
        with open(path, "w") as stream:
            for polymer in self.proteins + self.rnas + self.dnas:
                if types is None or polymer.type in types:
                    stream.write(f">{polymer.type.value}\n")
                    for i in range(0, len(polymer.sequence), line_length):
                        stream.write(polymer.sequence[i : i + line_length] + "\n")


def code1_to_code3(code1: str, polymer_type: PolymerType) -> str:
    return {
        PolymerType.PROTEIN: PROTEIN_CODES.get(code1) or PROTEIN_CODES["X"],
        PolymerType.RNA: RNA_CODES.get(code1) or RNA_CODES["X"],
        PolymerType.DNA: DNA_CODES.get(code1) or DNA_CODES["X"],
    }[polymer_type]


def guess_sequence_type(sequence: str) -> PolymerType:
    codes = set(sequence)
    if "U" in codes:
        return PolymerType.RNA
    if codes & set("DEFHIKLMNPQRSVWY"):
        return PolymerType.PROTEIN
    if codes == {"A"}:
        return PolymerType.PROTEIN
    if codes == {"G"}:
        return PolymerType.PROTEIN
    if "T" in codes:
        return PolymerType.DNA
    return PolymerType.RNA


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
            if line[1:3] in PIR_CODES and line[3:4] == ";":
                skip_line = True
            skip_lines = False
        elif line[:1] != ";" and not skip_lines:
            sequence += "".join(c for c in line if c.isalpha())
            if line[-1:] == "*":
                skip_lines = True
    if len(sequence) > 0:
        sequences.append(sequence)
    return sequences
