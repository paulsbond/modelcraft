import enum
import dataclasses
from typing import List, Tuple
import gemmi


@dataclasses.dataclass(unsafe_hash=True)
class Clash:
    """Clash stores data relating to a clash between a protein chain and nucleic acid
    chain
    """
    pro_chain_len: int
    na_chain_len: int
    pro_key: Tuple[str, str]
    na_key: Tuple[str, str]


@dataclasses.dataclass(unsafe_hash=True)
class ClashZone:
    """
    A class representing the ClashZone.

    Attributes:
    - pro_keys (List[Tuple[str, str]]): A list of tuples representing the protein Clash keys.
    - na_keys (List[Tuple[str, str]]): A list of tuples representing the nucleic acid Clash keys.
    """
    pro_keys: List[Tuple[str, str]]
    na_keys: List[Tuple[str, str]]


class StructureType(enum.Enum):
    """Enumeration class representing different types of structures."""
    protein = 1
    nucleic_acid = 2

    def is_same(self, kind: gemmi.ResidueKind):
        if self.value == self.protein.value:
            if kind == gemmi.ResidueKind.AA:
                return True
            return False
        else:
            if kind in [gemmi.ResidueKind.DNA, gemmi.ResidueKind.RNA]:
                return True
            return False
