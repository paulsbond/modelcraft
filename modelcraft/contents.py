from functools import reduce
from math import gcd
from typing import List, Optional
import json
import modelcraft.pdbe as pdbe
from .carb import Carb
from .ligand import Ligand
from .polymer import Polymer, PolymerType
from .residues import is_buffer


class AsuContents:
    def __init__(self, path_or_pdbid: Optional[str] = None):
        self.copies = None
        self.proteins: List[Polymer] = []
        self.rnas: List[Polymer] = []
        self.dnas: List[Polymer] = []
        self.carbs: List[Carb] = []
        self.ligands: List[Ligand] = []
        self.buffers: List[str] = []

        if path_or_pdbid is not None:
            if len(path_or_pdbid) == 4:
                self._from_pdbid(path_or_pdbid)
            else:
                if path_or_pdbid[-5:] == ".json":
                    self._from_json_file(path_or_pdbid)
                else:
                    self._from_sequence_file(path_or_pdbid)

    def add_polymer(self, polymer: Polymer) -> None:
        if polymer.type == PolymerType.PROTEIN:
            self.proteins.append(polymer)
        if polymer.type == PolymerType.RNA:
            self.rnas.append(polymer)
        if polymer.type == PolymerType.DNA:
            self.dnas.append(polymer)

    def _from_sequence_file(
        self, path: str, polymer_type: Optional[PolymerType] = None
    ) -> None:
        for polymer in Polymer.from_sequence_file(path, polymer_type):
            self.add_polymer(polymer)

    def _from_json_file(self, path: str) -> None:
        with open(path) as stream:
            contents = json.load(stream)
        self.copies = contents.get("copies")
        for obj in contents.get("proteins") or []:
            polymer = Polymer.from_json(obj)
            polymer.type = PolymerType.PROTEIN
            self.proteins.append(polymer)
        for obj in contents.get("rnas") or []:
            polymer = Polymer.from_json(obj)
            polymer.type = PolymerType.RNA
            self.rnas.append(polymer)
        for obj in contents.get("dnas") or []:
            polymer = Polymer.from_json(obj)
            polymer.type = PolymerType.DNA
            self.dnas.append(polymer)
        for obj in contents.get("carbs") or []:
            carb = Carb.from_json(obj)
            self.carbs.append(carb)
        for obj in contents.get("ligands") or []:
            ligand = Ligand.from_json(obj)
            self.ligands.append(ligand)
        self.buffers = contents.get("buffers") or []

    def _from_pdbid(self, pdbid: str) -> None:
        self.copies = 1
        for mol in pdbe.molecules(pdbid):
            if "sequence" in mol:
                polymer = Polymer.from_pdbe_molecule_dict(mol)
                self.add_polymer(polymer)
            if mol["molecule_type"] == "carbohydrate polymer":
                carb = Carb.from_pdbe_molecule_dict(mol)
                self.carbs.append(carb)
            if mol["molecule_type"] == "bound":
                ligand = Ligand.from_pdbe_molecule_dict(mol)
                if is_buffer(ligand.code):
                    self.buffers.append(ligand.code)
                elif ligand.code not in ("UNL", "UNX"):
                    self.ligands.append(ligand)
        self._divide_copies()

    def _divide_copies(self):
        copies = []
        for item in self.proteins + self.rnas + self.dnas + self.carbs + self.ligands:
            if item.copies is not None:
                copies.append(item.copies)
        divisor = copies[0] if len(copies) == 1 else reduce(gcd, copies)
        self.copies *= divisor
        for item in self.proteins + self.rnas + self.dnas + self.carbs + self.ligands:
            item.copies //= divisor

    def is_selenomet(self) -> bool:
        return len(self.proteins) > 0 and all(p.is_selenomet() for p in self.proteins)

    def volume(self) -> float:
        total = 0
        for item in self.proteins + self.rnas + self.dnas + self.carbs + self.ligands:
            total += item.volume() * (item.copies or 1)
        return total

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

    def write_sequence_file(
        self,
        path: str,
        polymer_type: Optional[PolymerType] = None,
        line_length: int = 60,
    ) -> None:
        with open(path, "w") as stream:
            for polymer in self.proteins + self.rnas + self.dnas:
                if polymer_type is None or polymer.type == polymer_type:
                    stream.write(f">{polymer.type.value}\n")
                    for i in range(0, len(polymer.sequence), line_length):
                        stream.write(polymer.sequence[i : i + line_length] + "\n")

    def write_json_file(self, path: str) -> None:
        with open(path, "w") as stream:
            json.dump(self.to_json(), stream, indent=2)
