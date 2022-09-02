import enum
import json
import gemmi
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


class PolymerType(enum.Enum):
    PROTEIN = "PROTEIN"
    RNA = "RNA"
    DNA = "DNA"


class Polymer:
    def __init__(
        self,
        sequence: str,
        copies: int,
        polymer_type: PolymerType,
        modifications: list = None,
    ):
        self.sequence = sequence.upper()
        self.copies = copies
        self.type = polymer_type
        self.modifications = modifications or []

    @classmethod
    def from_json(cls, component: dict, polymer_type: PolymerType) -> "Polymer":
        return cls(
            sequence=component["sequence"],
            copies=component["copies"],
            polymer_type=polymer_type,
            modifications=component.get("modifications"),
        )

    def to_json(self) -> dict:
        return {
            "sequence": self.sequence,
            "copies": self.copies,
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
    def __init__(self, codes: dict, copies: int):
        self.codes = codes
        self.copies = copies

    def __eq__(self, other) -> bool:
        if isinstance(other, Carb):
            return self.codes == other.codes
        return NotImplemented

    @classmethod
    def from_json(cls, component: dict) -> "Carb":
        return cls(codes=component["codes"], copies=component["copies"])

    def to_json(self) -> dict:
        return {"codes": self.codes, "copies": self.copies}

    def volume(self) -> float:
        monomers = sum(self.codes.values())
        volume = sum(monlib.volume(code) for code in self.codes)
        volume -= monomers * monlib.volume("HOH")
        return volume


class Ligand:
    def __init__(self, code: str, copies: int = None):
        self.code = code
        self.copies = copies

    def __eq__(self, other) -> bool:
        if isinstance(other, Ligand):
            return self.code == other.code
        return NotImplemented

    @classmethod
    def from_json(cls, component: dict) -> "Ligand":
        return cls(code=component["code"], copies=component["copies"])

    def to_json(self) -> dict:
        return {"code": self.code, "copies": self.copies}


class AsuContents:
    def __init__(self):
        self.proteins = []
        self.rnas = []
        self.dnas = []
        self.carbs = []
        self.ligands = []
        self.buffers = []

    @classmethod
    def from_file(cls, path: str) -> "AsuContents":
        contents = cls()
        with open(path) as stream:
            contents_json = json.load(stream)
        for obj in contents_json.get("proteins") or []:
            polymer = Polymer.from_json(obj, PolymerType.PROTEIN)
            polymer.type = PolymerType.PROTEIN
            contents.proteins.append(polymer)
        for obj in contents_json.get("rnas") or []:
            polymer = Polymer.from_json(obj, PolymerType.RNA)
            polymer.type = PolymerType.RNA
            contents.rnas.append(polymer)
        for obj in contents_json.get("dnas") or []:
            polymer = Polymer.from_json(obj, PolymerType.DNA)
            polymer.type = PolymerType.DNA
            contents.dnas.append(polymer)
        for obj in contents_json.get("carbs") or []:
            carb = Carb.from_json(obj)
            contents.carbs.append(carb)
        for obj in contents_json.get("ligands") or []:
            ligand = Ligand.from_json(obj)
            contents.ligands.append(ligand)
        contents.buffers = contents_json.get("buffers") or []
        return contents

    @classmethod
    def from_pdbe(cls, entry_id: str) -> "AsuContents":
        contents = cls()
        for mol in pdbe.molecule_dicts(entry_id):
            molecule_type = mol["molecule_type"].lower()
            if molecule_type.startswith("polypeptide"):
                protein = _polymer_from_pdbe_molecule_dict(mol, PolymerType.PROTEIN)
                contents.proteins.append(protein)
            elif molecule_type.startswith("polyribonucleotide"):
                rna = _polymer_from_pdbe_molecule_dict(mol, PolymerType.RNA)
                contents.rnas.append(rna)
            elif molecule_type.startswith("polydeoxyribonucleotide"):
                dna = _polymer_from_pdbe_molecule_dict(mol, PolymerType.DNA)
                contents.dnas.append(dna)
            elif molecule_type.startswith("carbohydrate"):
                carb = _carb_from_pdbe_molecule_dict(mol)
                contents.carbs.append(carb)
            elif molecule_type.startswith("bound"):
                ligand = _ligand_from_pdbe_molecule_dict(mol)
                if monlib.is_buffer(ligand.code):
                    contents.buffers.append(ligand.code)
                elif ligand.code not in ("UNL", "UNX"):
                    contents.ligands.append(ligand)
        return contents

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
        components = self.proteins + self.rnas + self.dnas + self.carbs + self.ligands
        return sum(component.volume() * component.copies for component in components)

    def solvent_fraction(
        self, cell: gemmi.UnitCell, spacegroup: gemmi.SpaceGroup
    ) -> float:
        asu_volume = cell.volume / len(spacegroup.operations())
        return 1 - self.volume() / asu_volume

    def to_json(self) -> list:
        return {
            "proteins": [protein.to_json() for protein in self.proteins],
            "rnas": [rna.to_json() for rna in self.rnas],
            "dnas": [dna.to_json() for dna in self.dnas],
            "carbs": [carb.to_json() for carb in self.carbs],
            "ligands": [ligand.to_json() for ligand in self.ligands],
            "buffers": self.buffers,
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


def _polymer_from_pdbe_molecule_dict(mol: dict, polymer_type: PolymerType) -> Polymer:
    mod_indices = {}
    for index, mod in mol["pdb_sequence_indices_with_multiple_residues"].items():
        code1 = mod["one_letter_code"]
        code3 = mod["three_letter_code"]
        if code3 not in ("DA", "DC", "DG", "DT"):
            key = code1, code3
            mod_indices.setdefault(key, []).append(index)
    modifications = []
    for key in mod_indices:
        code1, code3 = key
        total = mol["sequence"].count(code1)
        if code1 == "M" and mol["sequence"][0] == "M":
            total -= 1
        if len(mod_indices[key]) >= total:
            modifications.append(f"{code1}->{code3}")
        else:
            modifications.extend(f"{index}->{code3}" for index in mod_indices[key])
    return Polymer(
        sequence=mol["sequence"],
        copies=mol["number_of_copies"],
        polymer_type=polymer_type,
        modifications=modifications,
    )


def _carb_from_pdbe_molecule_dict(mol: dict) -> Carb:
    codes = mol["carb_codes"]
    length = sum(codes.values())
    copies = mol["number_of_copies"] // length
    return Carb(codes=codes, copies=copies)


def _ligand_from_pdbe_molecule_dict(mol: dict) -> Ligand:
    return Ligand(code=mol["chem_comp_ids"][0], copies=mol["number_of_copies"])
