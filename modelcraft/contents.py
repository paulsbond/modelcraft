from enum import Enum
from typing import Iterator, List, Optional
import modelcraft.residues as residues


class PolymerType(Enum):
    PROTEIN = "protein"
    RNA = "rna"
    DNA = "dna"

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
        label: str = "",
        copies: int = 1,
        polymer_type: Optional[PolymerType] = None,
    ):
        self.sequence = sequence.upper()
        self.label = label
        self.copies = copies
        if polymer_type is None:
            self.polymer_type = PolymerType.from_sequence(sequence)
        else:
            self.polymer_type = polymer_type

    def __eq__(self, other) -> bool:
        if not isinstance(other, Polymer):
            return NotImplemented
        return self.sequence == other.sequence

    @classmethod
    def from_sequence_file(
        cls, path: str, polymer_type: Optional[PolymerType] = None
    ) -> Iterator["Polymer"]:
        label = ""
        sequence = ""
        with open(path) as sequence_file:
            for line in sequence_file:
                if line[0] == ">":
                    if len(sequence) > 0:
                        yield Polymer(
                            sequence=sequence, label=label, polymer_type=polymer_type,
                        )
                    label = line[1:].split()[0]
                    sequence = ""
                else:
                    sequence += line.strip()
        if len(sequence) > 0:
            yield Polymer(sequence=sequence, label=label, polymer_type=polymer_type)


class Ligand:
    def __init__(self, code: str, copies: Optional[int] = None):
        self.code = code
        self.copies = copies


class AsuContents:
    def __init__(self, path: Optional[str] = None):
        self.polymers: List[Polymer] = []
        self.ligands: List[Ligand] = []
        if path is not None:
            self.polymers.extend(Polymer.from_sequence_file(path))

    def sequence_file_lines(
        self, polymer_type: Optional[PolymerType] = None, line_length: int = 60
    ) -> Iterator[str]:
        for polymer in self.polymers:
            if polymer_type is None or polymer.polymer_type == polymer_type:
                yield f"> {polymer.label}\n"
                for i in range(0, len(polymer.sequence), line_length):
                    yield polymer.sequence[i : i + line_length] + "\n"

    def write_sequence_file(
        self,
        path: str,
        polymer_type: Optional[PolymerType] = None,
        line_length: int = 60,
    ):
        with open(path, "w") as sequence_file:
            for line in self.sequence_file_lines(
                polymer_type=polymer_type, line_length=line_length
            ):
                sequence_file.write(line)
