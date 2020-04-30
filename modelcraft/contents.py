from dataclasses import dataclass
from enum import Enum
from typing import Iterator, Optional
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


@dataclass
class Ligand:
    code: str
    copies: Optional[int] = None


def read_sequence_file(
    path: str, polymer_type: Optional[PolymerType] = None
) -> Iterator[Polymer]:
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


# class AsuContents:
#     def __init__(self, relative=False):
#         self.polymers = []
#         self.ligands = []
#         self.heavy_atoms = []
#         self.relative = relative

#     def add_from_sequence_file(self, path, polymer_type="auto"):
#         for record in Bio.SeqIO.parse(path, "fasta"):
#             sequence = str(record.seq).upper()
#             if polymer_type == "auto":
#                 polymer_type = determine_polymer_type_from_sequence(sequence)
#             if polymer_type == "protein":
#                 self.polymers.append(Protein(sequence))
#             elif polymer_type == "rna":
#                 self.polymers.append(Rna(sequence))
#             elif polymer_type == "dna":
#                 self.polymers.append(Dna(sequence))
#             else:
#                 raise ValueError("Unknown polymer type: %s" % polymer_type)

#     def add_from_coordinate_file(self, path):
#         records = Bio.SeqIO.parse(path, "pdb-seqres")  # pdb-seqres cif-seqres pdb-atom cif-atom

#     def write_protein(self, path: str):
#         pass


# def print_sequence(sequence: str, line_length: int = 60) -> None:
#     for i in range(0, len(sequence), line_length):
#         print(sequence[i : i + line_length])
