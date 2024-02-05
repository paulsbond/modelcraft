from typing import Dict, Tuple, List
import gemmi
import numpy as np
import scipy
import modelcraft


def calculate_stats_per_residue(
    fphi_diff: modelcraft.DataItem,
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
    correlations = {}
    for key, value in difference_scores.items():
        correlations[key]["difference"] = value[0] / value[1] if value[1] > 0 else 0
    # Convert Difference Score to Z Score
    differences = [v["difference"] for v in correlations.values()]
    z_differences = scipy.stats.zscore(differences)
    for index, (k, v) in enumerate(correlations.items()):
        z_differences[index] = 0.5 * scipy.special.erfc(
            (1 / np.sqrt(2)) * z_differences[index]
        )
        z_differences[index] = scipy.stats.norm.cdf(-abs(z_differences[index])) * 2
        correlations[k]["z_difference"] = z_differences[index]
    return correlations


def score_from_zone(
    zone: List[Tuple[str, str]],
    stats: Dict[Tuple[str, str], Dict[str, float]],
    structure: gemmi.Structure,
) -> float:
    score_sum = 0
    for key in zone:
        score_sum += score_from_key(key, stats, structure)
    if zone:
        return score_sum / len(zone)
    return 0


def score_from_key(
    key: Tuple[str, str],
    stats: Dict[Tuple[str, str], Dict[str, float]],
    structure: gemmi.Structure,
) -> float:
    chain_name, residue_seq_id = key
    chain = structure[0].find_chain(chain_name)
    residue = chain[residue_seq_id][0]

    previous_residue = chain.previous_residue(residue)
    next_residue = chain.next_residue(residue)

    score = 0.0
    count = 1
    score -= stats[key].get("z_difference")

    if previous_residue:
        previous_residue_key = (chain_name, str(previous_residue.seqid))
        if previous_residue_key in stats:
            score -= stats[previous_residue_key].get("z_difference")
            count += 1

    if next_residue:
        next_residue_key = (chain_name, str(next_residue.seqid))
        if next_residue_key in stats:
            score -= stats[next_residue_key].get("z_difference")
            count += 1

    return score / count
