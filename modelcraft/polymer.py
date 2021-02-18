from typing import Iterator, List, Optional
import enum
import modelcraft.residues as residues


class PolymerType(enum.Enum):
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
            copies=mol["number_of_copies"],
            polymer_type=PolymerType.parse(mol["molecule_type"]),
            modifications=_modifications_in_pdbe_molecule_dict(mol),
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


def _modifications_in_pdbe_molecule_dict(mol: dict) -> List[str]:
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
