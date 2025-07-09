import random
from typing import Dict, Tuple
import gemmi


def random_id(length: int = 10) -> str:
    chars = "23456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz"
    return "".join(random.choice(chars) for _ in range(length))

def check_nucleic(name: str) -> bool:
    """Check if a residue is a nucleic acid using gemmi tables"""
    kind = gemmi.find_tabulated_residue(name)
    return kind.is_nucleic_acid()


def find_nucleic_acid_only_chains(structure) -> Dict[str, Tuple[str, str]]:
    """Find chains which contain only nucleic acids"""
    nucleic_acid_only_chains = {}
    for chain in structure[0]:
        nucleic_acid_chain = all(map(lambda residue: check_nucleic(residue.name), chain))
        if not nucleic_acid_chain:
            continue

        start_residue = chain[0].seqid.num
        end_residue = chain[-1].seqid.num
        nucleic_acid_only_chains[chain.name] = (str(start_residue), str(end_residue))
    return nucleic_acid_only_chains