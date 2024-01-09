import dataclasses
from typing import Tuple, List, Set

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


def identify_clashes(protein_structure: gemmi.Structure, nucleic_acid_structure: gemmi.Structure, search: gemmi.NeighborSearch) -> Set[Clash]:
    # pro_neighbour_search = gemmi.NeighborSearch(
    #     protein_structure[0], protein_structure.cell, 1.5).populate()

    clashes = set()

    for chain in protein_structure[0]:
        for residue in chain:
            for atom in residue:
                near_atoms = search.find_atoms(atom.pos, alt='\0', radius=1)
                for near_atom in near_atoms:
                    near_atom_cra = near_atom.to_cra(nucleic_acid_structure[0])
                    detected_clash = Clash(pro_chain_len=len(chain),
                                           na_chain_len=len(near_atom_cra.chain),
                                           pro_key=(chain.name, str(residue.seqid)),
                                           na_key=(near_atom_cra.chain.name,
                                                   str(near_atom_cra.residue.seqid)),
                                           )
                    clashes.add(detected_clash)

    return clashes


def inflate_bfactors(residue: gemmi.Residue, value: float) -> gemmi.Residue:
    for atom in residue:
        atom.b_iso = value


def extract_residue(clash_key: Tuple[str, str], structure: gemmi.Structure) -> gemmi.Residue:
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
