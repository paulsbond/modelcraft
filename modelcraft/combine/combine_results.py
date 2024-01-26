import enum
from typing import Set, List
import gemmi
from modelcraft.jobs.refmac import RefmacResult
from modelcraft.combine.clashes import identify_clashes, identify_clash_zones
from modelcraft.combine.statistics import calculate_stats_per_residue, score_from_zone
from modelcraft.combine.types import Clash, ClashZone, StructureType


def combine(pipeline, buccaneer_results: RefmacResult, nautilus_results: RefmacResult):
    """
    Combines the results of the from Buccaneer and Nautilus.

    :param pipeline: The current pipeline object.
    :param buccaneer_results: The result object from Buccaneer.
    :param nautilus_results: The result object from Nautilus.
    """
    pro_neighbour_search = gemmi.NeighborSearch(
        buccaneer_results.structure[0], buccaneer_results.structure.cell, 1.5).populate()

    na_neighbour_search = gemmi.NeighborSearch(
        nautilus_results.structure[0], nautilus_results.structure.cell, 1.5).populate()

    na_stats = calculate_stats_per_residue(fphi_calc=nautilus_results.fphi_calc, fphi_best=pipeline.current_fphi_best,
                                           fphi_diff=nautilus_results.fphi_diff, search=na_neighbour_search,
                                           structure=nautilus_results.structure)

    pro_stats = calculate_stats_per_residue(fphi_calc=buccaneer_results.fphi_calc, fphi_best=pipeline.current_fphi_best,
                                            fphi_diff=buccaneer_results.fphi_diff, search=pro_neighbour_search,
                                            structure=buccaneer_results.structure)

    clashes: Set[Clash] = identify_clashes(buccaneer_results.structure,
                                           nautilus_results.structure,
                                           search=na_neighbour_search)

    clash_zones: List[ClashZone] = identify_clash_zones(clashes)

    to_remove = set()

    for clash_zone in clash_zones:
        pro_total_score = score_from_zone(zone=clash_zone.pro_keys, stats=pro_stats, structure=buccaneer_results.structure)
        na_total_score = score_from_zone(zone=clash_zone.na_keys, stats=na_stats, structure=nautilus_results.structure)

        if na_total_score > pro_total_score:
            for pro_key in clash_zone.pro_keys:
                to_remove.add((StructureType.protein, *pro_key))
        if pro_total_score > na_total_score:
            for na_key in clash_zone.na_keys:
                to_remove.add((StructureType.nucleic_acid, *na_key))

    combined_structure = rebuild_model(to_remove, buccaneer_structure=buccaneer_results.structure,
                                       nautilus_structure=nautilus_results.structure)
    return combined_structure


def rebuild_model(to_remove: Set, buccaneer_structure: gemmi.Structure, nautilus_structure: gemmi.Structure):
    """
    Rebuilds a model by combining two structures while removing specified residues.

    Parameters: to_remove (Set): A set containing tuples of residues to be removed. Each tuple should have the format
    (structure_type, chain_name, residue_id).
    buccaneer_structure (gemmi.Structure): The structure containing
    residues from the Buccaneer model.
    nautilus_structure (gemmi.Structure): The structure containing residues from
    the Nautilus model.

    Returns:
        combined_structure (gemmi.Structure): The combined structure with the removed residues.

    """
    combined_structure = gemmi.Structure()
    combined_structure.cell = buccaneer_structure.cell
    combined_structure.spacegroup_hm = buccaneer_structure.spacegroup_hm
    combined_model = gemmi.Model(buccaneer_structure[0].name)
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    for n_ch, chain in enumerate(buccaneer_structure[0]):
        to_add_chain = gemmi.Chain(alphabet[n_ch])
        for residue in chain:
            if (StructureType.protein, chain.name, str(residue.seqid)) in to_remove:
                continue

            to_add_chain.add_residue(residue)

        combined_model.add_chain(to_add_chain)

    for n_ch, chain in enumerate(nautilus_structure[0]):
        to_add_chain = gemmi.Chain(alphabet[len(combined_model) + n_ch])
        for residue in chain:
            if (StructureType.nucleic_acid, chain.name, str(residue.seqid)) in to_remove:
                continue

            # Only add nucleic acid, otherwise, protein chains will be duplicated
            residue_kind: gemmi.ResidueInfo = gemmi.find_tabulated_residue(residue.name)
            if residue_kind.is_nucleic_acid():
                to_add_chain.add_residue(residue)

        combined_model.add_chain(to_add_chain)

    combined_structure.add_model(combined_model)
    return combined_structure
