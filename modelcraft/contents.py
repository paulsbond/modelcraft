import enum
import functools
import json
import math
from . import monlib
from . import pdbe

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

    def __str__(self) -> str:
        kind = {
            PolymerType.PROTEIN: "Protein",
            PolymerType.RNA: "RNA",
            PolymerType.DNA: "DNA",
        }[self.type]
        s = f"{kind} with {len(self.sequence)} residues: "
        if len(self.sequence) > 9:
            s += f"{self.sequence[:3]}...{self.sequence[-3:]}"
        else:
            s += f"{self.sequence:9}"
        return s

    @classmethod
    def from_json(cls, component: dict, polymer_type: PolymerType) -> "Polymer":
        return cls(
            sequence=component["sequence"],
            stoichiometry=component.get("stoichiometry"),
            polymer_type=polymer_type,
            modifications=component.get("modifications"),
        )

    @classmethod
    def from_pdbe(cls, mol: dict, polymer_type: PolymerType) -> "Polymer":
        mod_indices = {}
        for index, mod in mol["pdb_sequence_indices_with_multiple_residues"].items():
            code1 = mod["one_letter_code"]
            code3 = mod["three_letter_code"]
            if code3 not in ("DA", "DC", "DG", "DT"):
                key = code1, code3
                mod_indices.setdefault(key, []).append(index)
        modifications = []
        for (code1, code3), indices in mod_indices.items():
            total = mol["sequence"].count(code1)
            if code1 == "M" and mol["sequence"][0] == "M":
                total -= 1
            if len(indices) >= total:
                modifications.append(f"{code1}->{code3}")
            else:
                modifications.extend(f"{index}->{code3}" for index in indices)
        return cls(
            sequence=mol["sequence"],
            stoichiometry=mol["number_of_copies"],
            polymer_type=polymer_type,
            modifications=modifications,
        )

    @classmethod
    def from_sequence_file(cls, path: str, polymer_type: PolymerType = None):
        with open(path, encoding="utf-8") as stream:
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

    def weight(self) -> float:
        codes = self.residue_codes(modified=False)
        weight = sum(monlib.weight(code) for code in codes)
        weight -= monlib.weight("HOH") * (len(codes) - 1)
        return weight

    def volume(self) -> float:
        density = 1.35 if self.type == PolymerType.PROTEIN else 2.0
        return self.weight() / (density * 0.602214)


class Carb:
    def __init__(self, codes: dict, stoichiometry: int = None):
        self.codes = codes
        self.stoichiometry = stoichiometry

    def __eq__(self, other) -> bool:
        if isinstance(other, Carb):
            return self.codes == other.codes
        return NotImplemented

    def __str__(self) -> str:
        s = "Carb:"
        for code, count in self.codes.items():
            s += f" {count}x{code}"
        return s

    @classmethod
    def from_json(cls, component: dict) -> "Carb":
        return cls(
            codes=component["codes"],
            stoichiometry=component.get("stoichiometry"),
        )

    @classmethod
    def from_pdbe(cls, mol: dict) -> "Carb":
        codes = mol["carb_codes"]
        length = sum(codes.values())
        stoichiometry = mol["number_of_copies"] // length
        return cls(codes=codes, stoichiometry=stoichiometry)

    def to_json(self) -> dict:
        return {"codes": self.codes, "stoichiometry": self.stoichiometry}

    def volume(self) -> float:
        monomers = sum(self.codes.values())
        volume = sum(monlib.volume(code) for code in self.codes)
        volume -= monomers * monlib.volume("HOH")
        return volume


class Ligand:
    def __init__(self, code: str, stoichiometry: int = None):
        self.code = code
        self.stoichiometry = stoichiometry

    def __eq__(self, other) -> bool:
        if isinstance(other, Ligand):
            return self.code == other.code
        return NotImplemented

    def __str__(self) -> str:
        return f"Ligand: {self.code}"

    @classmethod
    def from_json(cls, component: dict) -> "Ligand":
        return cls(
            code=component["code"],
            stoichiometry=component.get("stoichiometry"),
        )

    @classmethod
    def from_pdbe(cls, mol: dict) -> "Ligand":
        return cls(code=mol["chem_comp_ids"][0], stoichiometry=mol["number_of_copies"])

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
    ):
        self.copies = copies
        self.proteins = proteins or []
        self.rnas = rnas or []
        self.dnas = dnas or []
        self.carbs = carbs or []
        self.ligands = ligands or []
        self.buffers = buffers or []

    @classmethod
    def from_file(cls, path: str) -> "AsuContents":
        if path[-5:] == ".json":
            return cls.from_json_file(path)
        return cls.from_sequence_file(path)

    @classmethod
    def from_json_file(cls, path: str) -> "AsuContents":
        with open(path, encoding="utf-8") as stream:
            contents = json.load(stream)
        return cls(
            copies=contents.get("copies"),
            proteins=[
                Polymer.from_json(obj, PolymerType.PROTEIN)
                for obj in contents.get("proteins", [])
            ],
            rnas=[
                Polymer.from_json(obj, PolymerType.RNA)
                for obj in contents.get("rnas", [])
            ],
            dnas=[
                Polymer.from_json(obj, PolymerType.DNA)
                for obj in contents.get("dnas", [])
            ],
            carbs=[Carb.from_json(obj) for obj in contents.get("carbs", [])],
            ligands=[Ligand.from_json(obj) for obj in contents.get("ligands", [])],
            buffers=contents.get("buffers", []),
        )

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

    @classmethod
    def from_pdbe(cls, entry_id: str) -> "AsuContents":
        contents = cls(copies=1)
        for mol in pdbe.molecule_dicts(entry_id):
            molecule_type = mol["molecule_type"].lower()
            if molecule_type.startswith("polypeptide"):
                protein = Polymer.from_pdbe(mol, PolymerType.PROTEIN)
                contents.proteins.append(protein)
            elif molecule_type.startswith("polyribonucleotide"):
                rna = Polymer.from_pdbe(mol, PolymerType.RNA)
                contents.rnas.append(rna)
            elif molecule_type.startswith("polydeoxyribonucleotide"):
                dna = Polymer.from_pdbe(mol, PolymerType.DNA)
                contents.dnas.append(dna)
            elif molecule_type.startswith("carbohydrate"):
                carb = Carb.from_pdbe(mol)
                contents.carbs.append(carb)
            elif molecule_type.startswith("bound"):
                ligand = _ligand_from_pdbe_molecule_dict(mol)
                if monlib.is_buffer(ligand.code) or ligand.code == "UNX":
                    contents.buffers.append(ligand.code)
                else:
                    contents.ligands.append(ligand)
        contents.divide_stoichiometry()
        return contents

    def components(self) -> list:
        return self.proteins + self.rnas + self.dnas + self.carbs + self.ligands

    def divide_stoichiometry(self):
        counts = []
        for component in self.components:
            if component.stoichiometry is not None:
                counts.append(component.stoichiometry)
        if len(counts) > 0:
            if len(counts) > 1:
                divisor = functools.reduce(math.gcd, counts)
            else:
                divisor = counts[0]
            if divisor > 1:
                self.copies *= divisor
                for component in self.components():
                    component.stoichiometry //= divisor

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

    def volume(self) -> float:
        return sum(c.volume() * c.stoichiometry for c in self.components())

    def to_json(self) -> list:
        return {
            "copies": self.copies,
            "proteins": [protein.to_json() for protein in self.proteins],
            "rnas": [rna.to_json() for rna in self.rnas],
            "dnas": [dna.to_json() for dna in self.dnas],
            "carbs": [carb.to_json() for carb in self.carbs],
            "ligands": [ligand.to_json() for ligand in self.ligands],
            "buffers": self.buffers,
        }

    def write_json_file(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as stream:
            json.dump(self.to_json(), stream, indent=2)

    def write_sequence_file(
        self, path: str, types: list = None, line_length: int = 60
    ) -> None:
        with open(path, "w", encoding="utf-8") as stream:
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
