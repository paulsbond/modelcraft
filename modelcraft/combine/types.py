import dataclasses
from typing import List, Tuple


@dataclasses.dataclass(unsafe_hash=True)
class Clash:
    """Clash stores data relating to a clash between a protein chain and nucleic acid
    chain
    """

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
