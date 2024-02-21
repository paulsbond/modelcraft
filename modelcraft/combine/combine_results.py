from typing import Set, List
import gemmi
from ..jobs.refmac import RefmacResult
from .clashes import identify_clashes, identify_clash_zones
from .statistics import calculate_stats_per_residue, score_from_zone
from .types import Clash, ClashZone


def combine(buccaneer: RefmacResult, nautilus: RefmacResult):
    pro_neighbour_search = gemmi.NeighborSearch(
        buccaneer.structure[0], buccaneer.structure.cell, 5
    ).populate()

    na_neighbour_search = gemmi.NeighborSearch(
        nautilus.structure[0], nautilus.structure.cell, 5
    ).populate()

    na_stats = calculate_stats_per_residue(
        fphi_diff=nautilus.fphi_diff,
        search=na_neighbour_search,
        structure=nautilus.structure,
    )

    pro_stats = calculate_stats_per_residue(
        fphi_diff=buccaneer.fphi_diff,
        search=pro_neighbour_search,
        structure=buccaneer.structure,
    )

    clashes: Set[Clash] = identify_clashes(
        buccaneer.structure,
        nautilus.structure,
        search=na_neighbour_search,
    )

    clash_zones: List[ClashZone] = identify_clash_zones(clashes, nautilus.structure)

    protein_to_remove = set()
    nucleic_acid_to_remove = set()

    for clash_zone in clash_zones:
        pro_total_score = score_from_zone(
            zone=clash_zone.pro_keys,
            stats=pro_stats,
            structure=buccaneer.structure,
        )
        na_total_score = score_from_zone(
            zone=clash_zone.na_keys,
            stats=na_stats,
            structure=nautilus.structure,
        )

        if na_total_score > pro_total_score:
            for pro_key in clash_zone.pro_keys:
                protein_to_remove.add(pro_key)

        if pro_total_score > na_total_score:
            for na_key in clash_zone.na_keys:
                nucleic_acid_to_remove.add(na_key)

    combined_structure = rebuild_model(
        protein_to_remove=protein_to_remove,
        nucleic_acid_to_remove=nucleic_acid_to_remove,
        buccaneer_structure=buccaneer.structure,
        nautilus_structure=nautilus.structure,
    )
    return combined_structure


def rebuild_model(
    protein_to_remove: Set,
    nucleic_acid_to_remove: Set,
    buccaneer_structure: gemmi.Structure,
    nautilus_structure: gemmi.Structure,
):
    combined_structure = gemmi.Structure()
    combined_structure.cell = buccaneer_structure.cell
    combined_structure.spacegroup_hm = buccaneer_structure.spacegroup_hm
    combined_model = gemmi.Model(buccaneer_structure[0].name)

    for chain in buccaneer_structure[0]:
        to_add_chain = gemmi.Chain(chain.name)
        for residue in chain:
            if (chain.name, str(residue.seqid)) in protein_to_remove:
                continue

            to_add_chain.add_residue(residue)

        combined_model.add_chain(to_add_chain, unique_name=True)

    for n_ch, chain in enumerate(nautilus_structure[0]):
        to_add_chain = gemmi.Chain(chain.name)
        for residue in chain:
            if (chain.name, str(residue.seqid)) in nucleic_acid_to_remove:
                continue

            # Only add nucleic acid, otherwise, protein chains will be duplicated
            residue_kind: gemmi.ResidueInfo = gemmi.find_tabulated_residue(residue.name)
            if residue_kind.is_nucleic_acid():
                to_add_chain.add_residue(residue)

        combined_model.add_chain(to_add_chain, unique_name=True)

    combined_structure.add_model(combined_model)
    return combined_structure
