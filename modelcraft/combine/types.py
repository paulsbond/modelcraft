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
    pro_keys: List[Tuple[str, str]]
    na_keys: List[Tuple[str, str]]


class StructureType(enum.Enum):
    """Enumeration class representing different types of structures."""
    protein = 1
    nucleic_acid = 2

    def is_same(self, kind: gemmi.ResidueKind):
        print(self.value, kind)
        if self.value == self.protein:
            if kind == gemmi.ResidueKind.AA:
                return True
            else:
                return False
        else:
            if kind in [gemmi.ResidueKind.DNA, gemmi.ResidueKind.RNA]:
                return True
            else:
                return False
