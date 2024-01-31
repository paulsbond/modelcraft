import dataclasses
import random
from typing import Dict, Tuple, List

import gemmi
import numpy as np
import scipy

import modelcraft

from modelcraft.combine.clashes import extract_residue
from modelcraft.combine.types import StructureType


def calculate_average_b_factor(residue: gemmi.Residue) -> float:
    """
    Calculate the average B-factor for a given residue.

    :param residue: gemmi.Residue object representing a residue
    :type residue: gemmi.Residue

    :return: the average B-factor for the residue
    :rtype: float
    """
    return np.average([atom.b_iso for atom in residue])


def calculate_stats_per_residue(fphi_calc: modelcraft.DataItem, fphi_best: modelcraft.DataItem,
                                fphi_diff: modelcraft.DataItem, search: gemmi.NeighborSearch,
                                structure: gemmi.Structure) -> Dict[
    Tuple[str, str], Dict[str, float]]:
    """Calculate RSCC for each residue in a structure

    RSCC Calculation - Paul Bond - https://github.com/paulsbond/modelcraft/blob/dev/modelcraft/rscc.py

    Args:
        fphi_calc (modelcraft.DataItem): calculated phases, from DensityCalculator() or refinement
        fphi_best (modelcraft.DataItem): best phases from refinement (before model building)
        fphi_diff (modelcraft.DataItem): difference phases from refinement (before model building)

        search (gemmi.NeighborSearch): constructed neighbor search object
        structure (gemmi.Structure): structure to calculate

    Returns:
        Dict[Tuple[str,str], Dict[str, float]: Dictionary of RSCC for each clashing residue in each chain
        In format {
            ('A', '1') : {'rscc': 0.5, 'z_rscc': 1.0, 'avg_b_factor': 0.5, 'z_b_factor': 1.0, 'z_difference': 0.4}
        }
    """
    calculated_map: gemmi.FloatGrid = fphi_calc.transform_f_phi_to_map(fphi_calc.label(0), fphi_calc.label(1))
    best_map: gemmi.FloatGrid = fphi_best.transform_f_phi_to_map(fphi_best.label(0), fphi_best.label(1))
    difference_map: gemmi.FloatGrid = fphi_diff.transform_f_phi_to_map(fphi_diff.label(0), fphi_diff.label(1))

    # Calculate RSCC
    residue_pairs = {}
    difference_scores = {}

    for point in calculated_map.masked_asu():
        position = calculated_map.point_to_position(point)
        mark = search.find_nearest_atom(position)
        if mark is not None:
            cra = mark.to_cra(structure[0])
            key = (cra.chain.name, str(cra.residue.seqid))

            value1 = point.value
            value2 = best_map.get_value(point.u, point.v, point.w)
            residue_pairs.setdefault(key, []).append((value1, value2))

            difference_value = difference_map.get_value(point.u, point.v, point.w)
            difference_scores.setdefault(key, [0, 0])
            if difference_value > 0:
                difference_scores[key][0] += difference_value
            difference_scores[key][1] += 1

    correlations = {}
    for key, pairs in residue_pairs.items():
        if len(pairs) > 1:
            values1, values2 = zip(*pairs)
            correlations[key] = {"rscc": np.corrcoef(values1, values2)[0, 1]}

        correlations[key]["difference"] = difference_scores[key][0] / difference_scores[key][1]

    # Convert RSCC to Z Score
    rsccs = [v["rscc"] for v in correlations.values()]
    z_rsccs = scipy.stats.zscore(rsccs)
    for index, (k, v) in enumerate(correlations.items()):
        z_rsccs[index] = 0.5 * scipy.special.erfc((1 / np.sqrt(2)) * z_rsccs[index])
        z_rsccs[index] = scipy.stats.norm.cdf(-abs(z_rsccs[index])) * 2
        correlations[k]["z_rscc"] = z_rsccs[index]

    # Calculate average B Factor
    for index, (k, v) in enumerate(correlations.items()):
        residue: gemmi.Residue = extract_residue(k, structure)

        residue_avg_b_factor = calculate_average_b_factor(residue)
        correlations[k]["b_factor"] = residue_avg_b_factor

    # Convert B Factor to Z Score
    b_factors = [v["b_factor"] for v in correlations.values()]
    z_bfactors = scipy.stats.zscore(b_factors)
    for index, (k, v) in enumerate(correlations.items()):
        z_bfactors[index] = 0.5 * scipy.special.erfc((1 / np.sqrt(2)) * z_bfactors[index])
        z_bfactors[index] = scipy.stats.norm.cdf(-abs(z_bfactors[index])) * 2
        correlations[k]["z_b_factor"] = z_bfactors[index]

    # Convert Difference Score to Z Score
    differences = [v["difference"] for v in correlations.values()]
    z_differences = scipy.stats.zscore(differences)
    for index, (k, v) in enumerate(correlations.items()):
        z_differences[index] = 0.5 * scipy.special.erfc((1 / np.sqrt(2)) * z_differences[index])
        z_differences[index] = scipy.stats.norm.cdf(-abs(z_differences[index])) * 2
        correlations[k]["z_difference"] = z_differences[index]

    return correlations


def score_from_zone(zone: List[Tuple[str, str]], stats: Dict[Tuple[str, str], Dict[str, float]],
                    structure: gemmi.Structure) -> float:
    """

    Calculate the score from a given zone.

    Parameters:
    - zone (List[Tuple[str, str]]): A list of tuples representing the zone.
    - stats (Dict[Tuple[str, str], Dict[str, float]]): A dictionary mapping keys to statistics.
    - structure (gemmi.Structure): The structure object.

    Note: If there are no residues in the zone the score is 0

    Returns:
    - float: The calculated score from the zone.

    """
    score_sum = 0
    for key in zone:
        score_sum += score_from_key(key, stats, structure)

    if zone:
        return score_sum / len(zone)
    else:
        return 0


def score_from_key(key: Tuple[str, str], stats: Dict[Tuple[str, str], Dict[str, float]],
                   structure: gemmi.Structure) -> float:
    """
    Calculates the score for a given key based on the statistics and structure.

    :param key: A tuple of chain name and residue sequence ID.
    :type key: str
    :param stats: A dictionary containing the statistical data for residues.
    :type stats: Dict[Tuple[str, str], Dict[str, float]]
    :param structure: The structure object containing the chain and residue information.
    :type structure: gemmi.Structure
    :return: The calculated score for the given key.
    :rtype: float
    """
    chain_name, residue_seq_id = key
    chain = structure[0].find_chain(chain_name)
    residue = chain[residue_seq_id][0]

    previous_residue = chain.previous_residue(residue)
    next_residue = chain.next_residue(residue)

    score = 0.0
    count = 1
    score += score_from_statistics(stats[key])

    if previous_residue:
        previous_residue_key = (chain_name, str(previous_residue.seqid))
        if previous_residue_key in stats:
            score += score_from_statistics(stats[previous_residue_key])
            count += 1

    if next_residue:
        next_residue_key = (chain_name, str(next_residue.seqid))
        if next_residue_key in stats:
            score += score_from_statistics(stats[next_residue_key])
            count += 1

    return score / count


def score_from_statistics(statistics: Dict[str, float]) -> float:
    """
    Calculate the score from the given statistics.

    Parameters:
        statistics (Dict[str, float]): A dictionary containing the statistics.
            - The key "z_rscc" represents the z-score for RSCC.
            - The key "z_bfactor" represents the z-score for B-factor.

    Returns:
        float: The calculated score obtained by multiplying the z-score for RSCC
            with the z-score for B-factor. If either of the z-scores is missing,
            the score will be 0.0.

    Example:
        statistics = {"z_rscc": 2.1, "z_bfactor": 2.0}
        score = score_from_statistics(statistics)
        print(score)
        # Output: 1.05
    """
    z_rscc = statistics.get("z_rscc")
    z_bfactor = statistics.get("z_b_factor")
    z_difference = statistics.get("z_difference")

    if z_rscc is None or z_bfactor is None:
        return 0.0
    return z_rscc - z_bfactor - z_difference
