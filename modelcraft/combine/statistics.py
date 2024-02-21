from typing import Dict, Tuple, List
import gemmi
from ..reflections import DataItem


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

        score = 0.0
        count = 1
        score -= stats.get(key, 0)

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

    if zone:
        return score_sum / len(zone)
    return 0
