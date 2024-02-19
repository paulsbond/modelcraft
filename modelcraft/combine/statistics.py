from typing import Dict, Tuple, List
import gemmi
import numpy as np
import scipy
from ..reflections import DataItem


def calculate_stats_per_residue(
    fphi_diff: DataItem,
    search: gemmi.NeighborSearch,
    structure: gemmi.Structure,
) -> Dict[Tuple[str, str], Dict[str, float]]:
    difference_map: gemmi.FloatGrid = fphi_diff.transform_f_phi_to_map(
        fphi_diff.label(0), fphi_diff.label(1)
    )
    difference_scores = {}
    for point in difference_map.masked_asu():
        position = difference_map.point_to_position(point)
        mark = search.find_nearest_atom(position)
        if mark is not None:
            cra = mark.to_cra(structure[0])
            key = (cra.chain.name, str(cra.residue.seqid))

            difference_value = difference_map.get_value(point.u, point.v, point.w)
            difference_scores.setdefault(key, [0, 0])
            if difference_value > 0:
                difference_scores[key][0] += difference_value
            difference_scores[key][1] += 1

    stats = {}
    for key, value in difference_scores.items():
        stats.setdefault(key, 0)
        stats[key] = value[0] / value[1] if value[1] > 0 else 0
        
    return stats


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

        if index == len(zone)-1:
            next_residue = chain.next_residue(residue)
            if not next_residue:
                continue

            next_residue_key = (chain_name, str(next_residue.seqid))
            if next_residue_key in stats:
                score -= stats.get(next_residue_key, 0)
                count += 1

        score_sum += score/count

    if zone:
        return score_sum / len(zone)
    return 0

