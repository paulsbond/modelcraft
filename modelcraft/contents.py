from enum import Enum
from typing import Iterator, List, Optional
import json
import modelcraft.residues as residues
import modelcraft.pdbe as pdbe


class PolymerType(Enum):
    PROTEIN = "PROTEIN"
    RNA = "RNA"
    DNA = "DNA"

    @classmethod
    def parse(cls, s: str) -> "PolymerType":
        try:
            return cls[s.upper()] if s is not None else None
        except KeyError:
            raise KeyError(f"Unknown polymer type: '{s}'")

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


class Polymer:
    def __init__(
        self,
        sequence: str,
        start: int = 1,
        label: Optional[str] = None,
        copies: Optional[int] = None,
        polymer_type: Optional[PolymerType] = None,
        modifications: Optional[List[str]] = None,
    ):
        self.sequence = sequence.upper()
        self.start = start
        self.label = label
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
        return Polymer(
            sequence=component["sequence"],
            start=component.get("start"),
            label=component.get("label"),
            copies=component.get("copies"),
            polymer_type=PolymerType.parse(component.get("type")),
            modifications=component.get("modifications"),
        )

    @classmethod
    def from_sequence_file(
        cls, path: str, polymer_type: Optional[PolymerType] = None
    ) -> Iterator["Polymer"]:
        label = None
        sequence = ""
        with open(path) as stream:
            for line in stream:
                if line[0] == ">":
                    if len(sequence) > 0:
                        yield Polymer(
                            sequence=sequence, label=label, polymer_type=polymer_type
                        )
                    label = line[1:].strip() if line[1:].strip() else None
                    sequence = ""
                elif line[0] != ";":
                    sequence += "".join(c for c in line if c.isalpha())
        if len(sequence) > 0:
            yield Polymer(sequence=sequence, label=label, polymer_type=polymer_type)

    def to_component_json(self) -> dict:
        return {
            "sequence": self.sequence,
            "label": self.label,
            "copies": self.copies,
            "type": str(self.type),
            "modifications": self.modifications,
        }


class Ligand:
    def __init__(
        self, code: str, label: Optional[str] = None, copies: Optional[int] = None
    ):
        self.code = code
        self.label = label
        self.copies = copies

    def __eq__(self, other) -> bool:
        if isinstance(other, Ligand):
            return self.code == other.code
        return NotImplemented

    @classmethod
    def from_component_json(cls, component: dict) -> "Ligand":
        return Ligand(
            code=component["code"],
            label=component.get("label"),
            copies=component.get("copies"),
        )

    def to_component_json(self) -> dict:
        return {"code": self.code, "label": self.label, "copies": self.copies}


class AsuContents:
    def __init__(self):
        self.polymers: List[Polymer] = []
        self.ligands: List[Ligand] = []

    def add_from_sequence_file(
        self, path: str, polymer_type: Optional[PolymerType] = None
    ) -> None:
        for polymer in Polymer.from_sequence_file(path, polymer_type):
            if polymer not in self.polymers:
                self.polymers.append(polymer)

    def add_from_json_file(self, path: str) -> None:
        with open(path) as stream:
            components = json.load(stream)
        for component in components:
            if "sequence" in component:
                polymer = Polymer.from_component_json(component)
                if polymer not in self.polymers:
                    self.polymers.append(polymer)
            elif "code" in component:
                ligand = Ligand.from_component_json(component)
                if ligand not in self.ligands:
                    self.ligands.append(ligand)

    def add_from_pdbid(self, pdbid: str) -> None:
        molecules = pdbe.molecules(pdbid)
        pass

    def write_sequence_file(
        self,
        path: str,
        polymer_type: Optional[PolymerType] = None,
        line_length: int = 60,
    ) -> None:
        with open(path, "w") as stream:
            for polymer in self.polymers:
                if polymer_type is None or polymer.type == polymer_type:
                    stream.write(f"> {polymer.label or polymer.type}\n")
                    for i in range(0, len(polymer.sequence), line_length):
                        stream.write(polymer.sequence[i : i + line_length] + "\n")

    def write_json_file(self, path: str) -> None:
        polymers = [polymer.to_component_json() for polymer in self.polymers]
        ligands = [ligand.to_component_json() for ligand in self.ligands]
        components = polymers + ligands
        with open(path, "w") as stream:
            json.dump(components, stream, indent=2)
