import dataclasses
from typing import Dict, List, Set, Tuple
import gemmi
from .jobs.refmac import RefmacResult
from .reflections import DataItem


@dataclasses.dataclass(unsafe_hash=True)
class Clash:
    pro_key: Tuple[str, str]
    na_key: Tuple[str, str]


@dataclasses.dataclass(unsafe_hash=True)
class ClashZone:
    pro_keys: List[Tuple[str, str]]
    na_keys: List[Tuple[str, str]]


def combine_results(buccaneer: RefmacResult, nautilus: RefmacResult) -> gemmi.Structure:
    protein_search = gemmi.NeighborSearch(
        buccaneer.structure[0], buccaneer.structure.cell, 5
    ).populate()
    nucleic_search = gemmi.NeighborSearch(
        nautilus.structure[0], nautilus.structure.cell, 5
    ).populate()

    na_stats = calculate_stats_per_residue(
        fphi_diff=nautilus.fphi_diff,
        search=nucleic_search,
        structure=nautilus.structure,
    )
    pro_stats = calculate_stats_per_residue(
        fphi_diff=buccaneer.fphi_diff,
        search=protein_search,
        structure=buccaneer.structure,
    )

    clashes = identify_clashes(
        buccaneer.structure, nautilus.structure, search=nucleic_search
    )
    clash_zones = identify_clash_zones(clashes, nautilus.structure)

    protein_to_remove = set()
    nucleic_acid_to_remove = set()

    for clash_zone in clash_zones:
        pro_total_score = score_from_zone(
            zone=clash_zone.pro_keys, stats=pro_stats, structure=buccaneer.structure
        )
        na_total_score = score_from_zone(
            zone=clash_zone.na_keys, stats=na_stats, structure=nautilus.structure
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
    protein_to_remove: set,
    nucleic_acid_to_remove: set,
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
            if (chain.name, str(residue.seqid)) not in protein_to_remove:
                to_add_chain.add_residue(residue)
        if len(to_add_chain) > 0:
            combined_model.add_chain(to_add_chain, unique_name=True)

    for chain in nautilus_structure[0]:
        to_add_chain = gemmi.Chain(chain.name)
        for residue in chain:
            if (chain.name, str(residue.seqid)) in nucleic_acid_to_remove:
                continue
            # Only add nucleic acid, otherwise, protein chains will be duplicated
            residue_kind = gemmi.find_tabulated_residue(residue.name)
            if residue_kind.is_nucleic_acid():
                to_add_chain.add_residue(residue)
        if len(to_add_chain) > 0:
            combined_model.add_chain(to_add_chain, unique_name=True)

    combined_structure.add_model(combined_model)
    return combined_structure


def extract_residue(
    clash_key: Tuple[str, str], structure: gemmi.Structure
) -> gemmi.Residue:
    chain_name = clash_key[0]
    residue_seqid = clash_key[1]
    chain = structure[0].find_chain(chain_name)
    residue_group = chain[residue_seqid]
    return residue_group[0]


def identify_clashes(
    protein_structure: gemmi.Structure,
    nucleic_structure: gemmi.Structure,
    search: gemmi.NeighborSearch,
) -> Set[Clash]:
    clashes = set()
    for chain in protein_structure[0]:
        for residue in chain:
            for atom in residue:
                near_atoms = search.find_atoms(atom.pos, alt="\0", radius=1)
                for near_atom in near_atoms:
                    near_chain = nucleic_structure[0][near_atom.chain_idx]
                    near_residue = near_chain[near_atom.residue_idx]
                    residue_info = gemmi.find_tabulated_residue(near_residue.name)
                    if residue_info.is_nucleic_acid():
                        clash = Clash(
                            pro_key=(chain.name, str(residue.seqid)),
                            na_key=(near_chain.name, str(near_residue.seqid)),
                        )
                        clashes.add(clash)
    return clashes


def identify_clash_zones(clashes: Set[Clash], nautilus_structure: gemmi.Structure):
    def is_sequential(key1: Tuple[str, str], key2: Tuple[str, str]) -> bool:
        key1_residue = extract_residue(key1, nautilus_structure)
        key2_residue = extract_residue(key2, nautilus_structure)

        key1_o3 = key1_residue.sole_atom("O3'")
        key1_p = key1_residue.sole_atom("P")

        key2_o3 = key2_residue.sole_atom("O3'")
        key2_p = key2_residue.sole_atom("P")

        delta_1 = key1_o3.pos - key2_p.pos
        delta_2 = key2_o3.pos - key1_p.pos

        return delta_1.length() < 2 or delta_2.length() < 2

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
        elif is_sequential(current_na_key, clash_key):
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
        clashing_pro_keys = sorted(clashing_pro_keys)
        zone = ClashZone(pro_keys=clashing_pro_keys, na_keys=clash_zone)
        clashing_zones.append(zone)

    return clashing_zones


def calculate_stats_per_residue(
    fphi_diff: DataItem,
    search: gemmi.NeighborSearch,
    structure: gemmi.Structure,
) -> Dict[Tuple[str, str], float]:
    density = fphi_diff.transform_f_phi_to_map(fphi_diff.label(0), fphi_diff.label(1))
    stats = {}
    for point in density.masked_asu():
        position = density.point_to_position(point)
        mark = search.find_nearest_atom(position)
        if mark is not None:
            cra = mark.to_cra(structure[0])
            key = (cra.chain.name, str(cra.residue.seqid))
            stats.setdefault(key, [0, 0])
            if point.value > 0:
                stats[key][0] += point.value
            stats[key][1] += 1
    return {k: v[0] / v[1] if v[1] > 0 else 0 for k, v in stats.items()}


def score_from_zone(
    zone: List[Tuple[str, str]],
    stats: Dict[Tuple[str, str], float],
    structure: gemmi.Structure,
) -> float:

    score_sum = 0

    for index, key in enumerate(zone):
        chain_name, residue_seq_id = key
        chain = structure[0].find_chain(chain_name)
        residue = chain[residue_seq_id][0]

        count = 1
        score = -stats.get(key, 0)

        if index == 0:
            previous_residue = chain.previous_residue(residue)
            if not previous_residue:
                continue

            previous_residue_key = (chain_name, str(previous_residue.seqid))
            if previous_residue_key in stats:
                score -= stats.get(previous_residue_key, 0)
                count += 1

        if index == len(zone) - 1:
            next_residue = chain.next_residue(residue)
            if not next_residue:
                continue

            next_residue_key = (chain_name, str(next_residue.seqid))
            if next_residue_key in stats:
                score -= stats.get(next_residue_key, 0)
                count += 1

        score_sum += score / count

    return score_sum / len(zone) if zone else 0
