import dataclasses
import gemmi
from .jobs.refmac import RefmacResult


@dataclasses.dataclass(unsafe_hash=True)
class _Clash:
    key1: tuple
    key2: tuple


@dataclasses.dataclass(unsafe_hash=True)
class _ClashZone:
    keys1: list
    keys2: list


def combine_results(result1: RefmacResult, result2: RefmacResult) -> gemmi.Structure:
    scores1 = _residue_scores(result1)
    scores2 = _residue_scores(result2)

    structure1 = result1.structure
    structure2 = result2.structure
    clashes = identify_clashes(structure1, structure2)
    clash_zones = identify_clash_zones(clashes, structure2)

    to_remove1 = set()
    to_remove2 = set()
    for clash_zone in clash_zones:
        total_score1 = score_from_zone(clash_zone.keys1, scores1, structure1)
        total_score2 = score_from_zone(clash_zone.keys2, scores2, structure2)
        if total_score2 > total_score1:
            for key1 in clash_zone.keys1:
                to_remove1.add(key1)
        else:
            for key2 in clash_zone.keys2:
                to_remove2.add(key2)

    combined_structure = rebuild_model(structure1, structure2, to_remove1, to_remove2)
    return combined_structure


def _residue_scores(result: RefmacResult) -> dict:
    search = gemmi.NeighborSearch(result.structure, 5).populate()
    fphi = result.fphi_diff
    density = fphi.transform_f_phi_to_map(fphi.label(0), fphi.label(1))
    stats = {}
    for point in density.masked_asu():
        position = density.point_to_position(point)
        mark = search.find_nearest_atom(position)
        if mark is not None:
            cra = mark.to_cra(result.structure[0])
            key = (cra.chain.name, str(cra.residue.seqid))
            stats.setdefault(key, [0, 0])
            if point.value > 0:
                stats[key][0] += point.value
            stats[key][1] += 1
    return {k: -v[0] / v[1] if v[1] > 0 else 0 for k, v in stats.items()}


def identify_clashes(structure1: gemmi.Structure, structure2: gemmi.Structure) -> set:
    search = gemmi.NeighborSearch(structure2[0], structure2.cell, 1).populate()
    clashes = set()
    for chain1 in structure1[0]:
        for residue1 in chain1:
            for atom1 in residue1:
                for atom2 in search.find_atoms(atom1.pos, alt="\0", radius=1):
                    chain2 = structure2[0][atom2.chain_idx]
                    residue2 = chain2[atom2.residue_idx]
                    if gemmi.find_tabulated_residue(residue2.name).is_nucleic_acid():
                        clash = _Clash(
                            key1=(chain1.name, str(residue1.seqid)),
                            key2=(chain2.name, str(residue2.seqid)),
                        )
                        clashes.add(clash)
    return clashes


def identify_clash_zones(clashes: set, structure2: gemmi.Structure):
    clash_map = {}
    for clash in clashes:
        clash_map.setdefault(clash.key2, []).append(clash.key1)

    zones2 = []
    zone2 = []
    previous_key2 = None
    for key2 in sorted(list(clash_map)):
        if previous_key2 is None or is_sequential(structure2, previous_key2, key2):
            zone2.append(key2)
        else:
            zones2.append(zone2)
            zone2 = [key2]
        previous_key2 = key2
    zones2.append(zone2)

    zones = []
    for zone2 in zones2:
        zone1 = set()
        for key2 in zone2:
            for key1 in clash_map.get(key2):
                zone1.add(key1)
        zone1 = sorted(zone1)
        zone = _ClashZone(keys1=zone1, keys2=zone2)
        zones.append(zone)
    return zones


def is_sequential(structure: gemmi.Structure, key1: tuple, key2: tuple) -> bool:
    key1_residue = extract_residue(key1, structure)
    key2_residue = extract_residue(key2, structure)
    key1_o3 = key1_residue.sole_atom("O3'").pos
    key1_p = key1_residue.sole_atom("P").pos
    key2_o3 = key2_residue.sole_atom("O3'").pos
    key2_p = key2_residue.sole_atom("P").pos
    return key1_o3.dist(key2_p) < 2 or key2_o3.dist(key1_p) < 2


def extract_residue(clash_key: tuple, structure: gemmi.Structure) -> gemmi.Residue:
    chain_name = clash_key[0]
    residue_seqid = clash_key[1]
    chain = structure[0].find_chain(chain_name)
    residue_group = chain[residue_seqid]
    return residue_group[0]


def score_from_zone(zone: list, scores: dict, structure: gemmi.Structure) -> float:
    score_sum = 0

    for index, key in enumerate(zone):
        chain_name, residue_seq_id = key
        chain = structure[0].find_chain(chain_name)
        residue = chain[residue_seq_id][0]

        score = scores.get(key, 0)
        count = 1

        if index == 0:
            previous_residue = chain.previous_residue(residue)
            if not previous_residue:
                continue

            previous_residue_key = (chain_name, str(previous_residue.seqid))
            if previous_residue_key in scores:
                score += scores.get(previous_residue_key, 0)
                count += 1

        if index == len(zone) - 1:
            next_residue = chain.next_residue(residue)
            if not next_residue:
                continue

            next_residue_key = (chain_name, str(next_residue.seqid))
            if next_residue_key in scores:
                score += scores.get(next_residue_key, 0)
                count += 1

        score_sum += score / count

    return score_sum / len(zone) if zone else 0


def rebuild_model(
    structure1: gemmi.Structrue,
    structure2: gemmi.Structrue,
    to_remove1: set,
    to_remove2: set,
):
    combined_structure = gemmi.Structure()
    combined_structure.cell = structure1.cell
    combined_structure.spacegroup_hm = structure1.spacegroup_hm
    combined_model = gemmi.Model(structure1[0].name)

    for chain in structure1[0]:
        to_add_chain = gemmi.Chain(chain.name)
        for residue in chain:
            if (chain.name, str(residue.seqid)) not in to_remove1:
                to_add_chain.add_residue(residue)
        if len(to_add_chain) > 0:
            combined_model.add_chain(to_add_chain, unique_name=True)

    for chain in structure2[0]:
        to_add_chain = gemmi.Chain(chain.name)
        for residue in chain:
            if (chain.name, str(residue.seqid)) in to_remove2:
                continue
            # Only add nucleic acid, otherwise, protein chains will be duplicated
            residue_kind = gemmi.find_tabulated_residue(residue.name)
            if residue_kind.is_nucleic_acid():
                to_add_chain.add_residue(residue)
        if len(to_add_chain) > 0:
            combined_model.add_chain(to_add_chain, unique_name=True)

    combined_structure.add_model(combined_model)
    return combined_structure
