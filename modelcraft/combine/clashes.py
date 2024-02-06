from typing import Tuple, List, Set
import gemmi
from .types import Clash, ClashZone


def identify_clashes(
    protein_structure: gemmi.Structure,
    nucleic_acid_structure: gemmi.Structure,
    search: gemmi.NeighborSearch,
) -> Set[Clash]:
    clashes = set()
    for chain in protein_structure[0]:
        for residue in chain:
            for atom in residue:
                near_atoms = search.find_atoms(atom.pos, alt="\0", radius=1)
                for near_atom in near_atoms:
                    near_chain: gemmi.Chain = nucleic_acid_structure[0][
                        near_atom.chain_idx
                    ]
                    near_residue: gemmi.Residue = near_chain[near_atom.residue_idx]
                    residue_info: gemmi.ResidueInfo = gemmi.find_tabulated_residue(
                        near_residue.name
                    )
                    if not residue_info.is_nucleic_acid():
                        continue

                    detected_clash = Clash(
                        pro_key=(chain.name, str(residue.seqid)),
                        na_key=(near_chain.name, str(near_residue.seqid)),
                    )
                    clashes.add(detected_clash)
    return clashes


def identify_clash_zones(clashes: Set[Clash]):
    def is_sequential(key1: Tuple[str, str], key2: Tuple[str, str]) -> bool:
        if abs(int(key1[1]) - int(key2[1])) == 1:
            return True
        return False

    clash_map = {}
    for clash in clashes:
        clash_map.setdefault(clash.na_key, []).append(clash.pro_key)
    sorted_keys = sorted(list(clash_map))
    sequential_na_clashes = []
    current_na_key = ()
    current_na_clash_zone = []
    for clash_key in sorted_keys:
        if not current_na_key:
            current_na_key = clash_key
            current_na_clash_zone.append(clash_key)
            continue
        if is_sequential(current_na_key, clash_key):
            current_na_clash_zone.append(clash_key)
            current_na_key = clash_key
        else:
            sequential_na_clashes.append(current_na_clash_zone)
            current_na_key = clash_key
            current_na_clash_zone = [clash_key]
    sequential_na_clashes.append(current_na_clash_zone)
    clashing_zones: List[ClashZone] = []
    for clash_zone in sequential_na_clashes:
        clashing_pro_keys = set()
        for key in clash_zone:
            clash_pro_values = clash_map.get(key)
            for value in clash_pro_values:
                clashing_pro_keys.add(value)
        zone = ClashZone(pro_keys=clashing_pro_keys, na_keys=clash_zone)
        clashing_zones.append(zone)
    return clashing_zones
