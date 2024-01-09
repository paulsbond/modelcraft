from typing import Set, Tuple
import gemmi
from modelcraft.jobs.refmac import RefmacResult
from modelcraft.combine.clashes import identify_clashes, Clash
from modelcraft.combine.statistics import calculate_stats_per_residue, score_from_key


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

    na_stats = calculate_stats_per_residue(fphi_calc=nautilus_results.fphi_calc,
                                           fphi_best=pipeline.current_fphi_best,
                                           search=na_neighbour_search,
                                           structure=nautilus_results.structure)

    pro_stats = calculate_stats_per_residue(fphi_calc=buccaneer_results.fphi_calc,
                                            fphi_best=pipeline.current_fphi_best,
                                            search=pro_neighbour_search,
                                            structure=buccaneer_results.structure)

    clashes: Set[Clash] = identify_clashes(buccaneer_results.structure,
                                           nautilus_results.structure,
                                           search=na_neighbour_search)

    to_remove = set()

    for clash in clashes:
        pro_total_score = score_from_key(key=clash.pro_key, stats=pro_stats, structure=buccaneer_results.structure)
        na_total_score = score_from_key(key=clash.na_key, stats=na_stats, structure=nautilus_results.structure)

        if na_total_score > pro_total_score:
            to_remove.add((0, *clash.pro_key))
        if pro_total_score > na_total_score:
            to_remove.add((1, *clash.na_key))

    combined_structure = rebuild_model(to_remove, buccaneer_structure=buccaneer_results.structure, nautilus_structure=nautilus_results.structure)
    return combined_structure


def rebuild_model(to_remove: Set, buccaneer_structure: gemmi.Structure, nautilus_structure: gemmi.Structure):
    combined_structure = gemmi.Structure()
    combined_structure.cell = buccaneer_structure.cell
    combined_structure.spacegroup_hm = buccaneer_structure.spacegroup_hm
    combined_model = gemmi.Model(buccaneer_structure[0].name)

    for n_ch, chain in enumerate(buccaneer_structure[0]):
        to_add_chain = gemmi.Chain(str(n_ch))
        for residue in chain:
            if (0, chain.name, str(residue.seqid)) in to_remove:
                continue

            to_add_chain.add_residue(residue)

        combined_model.add_chain(to_add_chain)

    for n_ch, chain in enumerate(nautilus_structure[0]):
        to_add_chain = gemmi.Chain(str(len(combined_model) + (n_ch)))
        for residue in chain:
            if (1, chain.name, str(residue.seqid)) in to_remove:
                continue

            to_add_chain.add_residue(residue)

        combined_model.add_chain(to_add_chain)

    combined_structure.add_model(combined_model)
    return combined_structure
