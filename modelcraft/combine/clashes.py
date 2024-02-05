import dataclasses
from typing import Tuple, List, Set
from modelcraft.combine.types import Clash, ClashZone
import gemmi


def identify_clashes(protein_structure: gemmi.Structure, nucleic_acid_structure: gemmi.Structure,
                     search: gemmi.NeighborSearch) -> Set[Clash]:
    """
    Identify Clashes

    Finds clashes between protein and nucleic acid structures using a neighbor search.

    Parameters:
    - protein_structure: gemmi.Structure - The protein structure to search for clashes.
    - nucleic_acid_structure: gemmi.Structure - The nucleic acid structure to search for clashes.
    - search: gemmi.NeighborSearch - The neighbor search used to find nearby atoms.

    Returns:
    - Set[Clash] - A set of Clash objects representing the detected clashes.

    Example usage:
    protein_structure = gemmi.read_structure('protein.pdb')
    nucleic_acid_structure = gemmi.read_structure('nucleic_acid.pdb')
    search = gemmi.NeighborSearch(protein_structure)
    clashes = identify_clashes(protein_structure, nucleic_acid_structure, search)
    for clash in clashes:
        print(clash)
    """
    clashes = set()

    for chain in protein_structure[0]:
        for residue in chain:
            for atom in residue:
                near_atoms = search.find_atoms(atom.pos, alt='\0', radius=1)
                for near_atom in near_atoms:
                    near_chain: gemmi.Chain = nucleic_acid_structure[0][near_atom.chain_idx]
                    near_residue: gemmi.Residue = near_chain[near_atom.residue_idx]
                    residue_info: gemmi.ResidueInfo = gemmi.find_tabulated_residue(near_residue.name)
                    if not residue_info.is_nucleic_acid():
                        continue

                    detected_clash = Clash(pro_chain_len=len(chain),
                                           na_chain_len=len(near_chain),
                                           pro_key=(chain.name, str(residue.seqid)),
                                           na_key=(near_chain.name, str(near_residue.seqid)))
                    clashes.add(detected_clash)

    return clashes


def identify_clash_zones(clashes: Set[Clash]):
    """
    Identify Clash Zones

    Identifies and groups clashes based on sequential keys.

    Parameters:
    - clashes: A set of Clash objects representing the clashes to be analyzed.

    Returns:
    - List[ClashZone]: A list of ClashZone objects representing the identified clash zones.

    """

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


def inflate_bfactors(residue: gemmi.Residue, value: float) -> gemmi.Residue:
    """
    Inflates the B-factors of all atoms in the given residue to the specified value.

    :param residue: The residue to inflate the B-factors for.
    :type residue: gemmi.Residue

    :param value: The value to set the B-factors to.
    :type value: float

    :return: The modified residue with inflated B-factors.
    :rtype: gemmi.Residue
    """
    for atom in residue:
        atom.b_iso = value


def extract_residue(clash_key: Tuple[str, str], structure: gemmi.Structure) -> gemmi.Residue:
    """
    Extracts a specific residue from a given gemmi.Structure based on the provided clash_key.

    Parameters:
    - clash_key (Tuple[str, str]): A tuple containing the chain name and residue seqid for the desired residue.
    - structure (gemmi.Structure): The gemmi.Structure object from which the residue should be extracted.

    Returns:
    - gemmi.Residue: The extracted residue.

    Example usage:
    clash_key = ("A", "123")
    structure = gemmi.Structure()
    residue = extract_residue(clash_key, structure)
    """
    chain_name = clash_key[0]
    residue_seqid = clash_key[1]
    chain = structure[0].find_chain(chain_name)
    residue_group = chain[residue_seqid]
    residue = residue_group[0]
    return residue


def inflate_bfactors_on_clash(protein_structure: gemmi.Structure,
                              nucleic_acid_structure: gemmi.Structure) -> gemmi.Structure:
    clashes: Set[Clash] = identify_clashes(protein_structure, nucleic_acid_structure)

    for chain in protein_structure[0]:
        for idx, residue in enumerate(chain):
            if (chain.name, str(residue.seqid)) in clashes:
                pro_residue = extract_residue((chain.name, residue.seqid), protein_structure)
                inflate_bfactors(pro_residue, 100.0)
                chain[idx] = pro_residue

    for chain in nucleic_acid_structure[0]:
        for idx, residue in enumerate(chain):
            if (chain.name, str(residue.seqid)) in clashes:
                na_residue = extract_residue((chain.name, residue.seqid), nucleic_acid_structure)
                inflate_bfactors(na_residue, 100.0)
                chain[idx] = na_residue

    structure = gemmi.Structure()
    structure.add_model(protein_structure[0])
    structure.add_model(nucleic_acid_structure[0])
    return structure
