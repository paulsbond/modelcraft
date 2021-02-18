from enum import Enum
from typing import Dict, Iterator, List, Optional
import json
import os
import modelcraft.residues as residues
import modelcraft.pdbe as pdbe


class PolymerType(Enum):
    PROTEIN = "PROTEIN"
    RNA = "RNA"
    DNA = "DNA"

    @classmethod
    def parse(cls, s: str) -> "PolymerType":
        if s.lower() in ("protein", "polypeptide(l)"):
            return cls.PROTEIN
        if s.lower() in ("rna", "polyribonucleotide"):
            return cls.RNA
        if s.lower() in ("dna", "polydeoxyribonucleotide"):
            return cls.DNA
        raise ValueError(f"Unknown polymer type: '{s}'")

    @classmethod
    def from_sequence(cls, sequence: str) -> "PolymerType":
        codes = set(sequence)
        if "U" in codes:
            return cls.RNA
        protein_codes = {residue.code1 for residue in residues.PROTEIN}
        unique_protein_codes = protein_codes - {"A", "C", "G", "T"}
        if codes & unique_protein_codes:
            return cls.PROTEIN
        if codes == {"A"}:
            return cls.PROTEIN
        if codes == {"G"}:
            return cls.PROTEIN
        if "T" in codes:
            return cls.DNA
        return cls.RNA


def modifications_in_pdbe_molecule_dict(mol: dict) -> List[str]:
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
        if len(indices[key]) == total:
            modifications.append(code3)
        else:
            modifications.extend(code3 + index for index in indices[key])
    return modifications


class Polymer:
    def __init__(
        self,
        sequence: str,
        start: Optional[int] = None,
        copies: Optional[int] = None,
        polymer_type: Optional[PolymerType] = None,
        modifications: Optional[List[str]] = None,
    ):
        self.sequence = sequence.upper()
        self.start = start or 1
        self.copies = copies
        if polymer_type is None:
            self.type = PolymerType.from_sequence(self.sequence)
        else:
            self.type = polymer_type
        self.modifications = modifications

    def __eq__(self, other) -> bool:
        if isinstance(other, Polymer):
            return (
                self.sequence == other.sequence
                and self.type == other.type
                and self.modifications == other.modifications
            )
        return NotImplemented

    @classmethod
    def from_component_json(cls, component: dict) -> "Polymer":
        return cls(
            sequence=component["sequence"],
            start=component.get("start"),
            copies=component.get("copies"),
            polymer_type=PolymerType.parse(component.get("type")),
            modifications=component.get("modifications"),
        )

    @classmethod
    def from_pdbe_molecule_dict(cls, mol: dict) -> "Polymer":
        return cls(
            sequence=mol["sequence"],
            start=mol["source"][0]["mappings"][0]["start"]["residue_number"],
            copies=mol["number_of_copies"],
            polymer_type=PolymerType.parse(mol["molecule_type"]),
            modifications=modifications_in_pdbe_molecule_dict(mol),
        )

    @classmethod
    def from_sequence_file(
        cls, path: str, polymer_type: Optional[PolymerType] = None
    ) -> Iterator["Polymer"]:
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

    def to_component_json(self) -> dict:
        return {
            "sequence": self.sequence,
            "type": self.type.value,
            "start": self.start,
            "copies": self.copies,
            "modifications": self.modifications,
        }


class Carb:
    def __init__(self, codes: Dict[str, int], copies: Optional[int] = None):
        self.codes = codes
        self.copies = copies

    def __eq__(self, other) -> bool:
        if isinstance(other, Ligand):
            return self.codes == other.codes and self.length == other.length
        return NotImplemented

    @classmethod
    def from_component_json(cls, component: dict) -> "Carb":
        return cls(codes=component["codes"], copies=component.get("copies"))

    @classmethod
    def from_pdbe_molecule_dict(cls, mol: dict) -> "Carb":
        codes = mol["carb_codes"]
        length = sum(codes.values())
        copies = mol["number_of_copies"] // length
        return cls(codes=codes, copies=copies)

    def to_component_json(self) -> dict:
        return {"codes": self.codes, "copies": self.copies}


class Ligand:
    def __init__(self, code: str, copies: Optional[int] = None):
        self.code = code
        self.copies = copies

    def __eq__(self, other) -> bool:
        if isinstance(other, Ligand):
            return self.code == other.code
        return NotImplemented

    @classmethod
    def from_component_json(cls, component: dict) -> "Ligand":
        return cls(code=component["code"], copies=component.get("copies"))

    @classmethod
    def from_pdbe_molecule_dict(cls, mol: dict) -> "Ligand":
        return cls(code=mol["chem_comp_ids"][0], copies=mol["number_of_copies"])

    def to_component_json(self) -> dict:
        return {"code": self.code, "copies": self.copies}


class AsuContents:
    def __init__(self, path_or_pdbid: Optional[str] = None):
        self.polymers: List[Polymer] = []
        self.carbs: List[Polymer] = []
        self.ligands: List[Ligand] = []
        if path_or_pdbid is not None:
            if len(path_or_pdbid) == 4:
                self.add_from_pdbid(path_or_pdbid)
            else:
                _, ext = os.path.splitext(path_or_pdbid)
                if ext == ".json":
                    self.add_from_json_file(path_or_pdbid)
                else:
                    self.add_from_sequence_file(path_or_pdbid)

    def add_from_sequence_file(
        self, path: str, polymer_type: Optional[PolymerType] = None
    ) -> None:
        for polymer in Polymer.from_sequence_file(path, polymer_type):
            self.polymers.append(polymer)

    def add_from_json_file(self, path: str) -> None:
        with open(path) as stream:
            components = json.load(stream)
        for component in components:
            if "sequence" in component:
                polymer = Polymer.from_component_json(component)
                self.polymers.append(polymer)
            elif "codes" in component:
                carb = Carb.from_component_json(component)
                self.carbs.append(carb)
            elif "code" in component:
                ligand = Ligand.from_component_json(component)
                self.ligands.append(ligand)

    def add_from_pdbid(self, pdbid: str) -> None:
        for mol in pdbe.molecules(pdbid):
            if "sequence" in mol:
                polymer = Polymer.from_pdbe_molecule_dict(mol)
                self.polymers.append(polymer)
            if mol["molecule_type"] == "carbohydrate polymer":
                carb = Carb.from_pdbe_molecule_dict(mol)
                self.carbs.append(carb)
            if mol["molecule_type"] == "bound":
                ligand = Ligand.from_pdbe_molecule_dict(mol)
                if ligand.code not in ("UNL", "UNX"):
                    self.ligands.append(ligand)

    def components_json(self) -> list:
        polymers = [polymer.to_component_json() for polymer in self.polymers]
        carbs = [carb.to_component_json() for carb in self.carbs]
        ligands = [ligand.to_component_json() for ligand in self.ligands]
        return polymers + carbs + ligands

    def write_sequence_file(
        self,
        path: str,
        polymer_type: Optional[PolymerType] = None,
        line_length: int = 60,
    ) -> None:
        with open(path, "w") as stream:
            for polymer in self.polymers:
                if polymer_type is None or polymer.type == polymer_type:
                    stream.write(f">{polymer.type.value}\n")
                    for i in range(0, len(polymer.sequence), line_length):
                        stream.write(polymer.sequence[i : i + line_length] + "\n")

    def write_json_file(self, path: str) -> None:
        with open(path, "w") as stream:
            json.dump(self.components_json(), stream, indent=2)
