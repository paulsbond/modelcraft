import dataclasses
import gemmi
from .jobs.refmac import RefmacResult


def combine_results(result1: RefmacResult, result2: RefmacResult) -> gemmi.Structure:
    residues1 = _my_residues(result1.structure)
    residues2 = _my_residues(result2.structure)
    _assign_scores(residues1, _scores(result1))
    _assign_scores(residues2, _scores(result2))
    clashes = _clashes(result1.structure, result2.structure)
    _assign_clashes(residues1, residues2, clashes)

    # to_remove1 = set()
    # to_remove2 = set()
    # for clash_zone in clash_zones:
    #     total_score1 = score_from_zone(clash_zone.keys1, scores1, structure1)
    #     total_score2 = score_from_zone(clash_zone.keys2, scores2, structure2)
    #     if total_score2 > total_score1:
    #         for key1 in clash_zone.keys1:
    #             to_remove1.add(key1)
    #     else:
    #         for key2 in clash_zone.keys2:
    #             to_remove2.add(key2)

    # combined_structure = rebuild_model(structure1, structure2, to_remove1, to_remove2)
    # return combined_structure


def _key(chain: gemmi.Chain, residue: gemmi.Residue) -> tuple:
    return (chain.name, residue.seqid.num, residue.seqid.icode)


@dataclasses.dataclass
class _MyResidue:
    def __init__(self, chain: gemmi.Chain, residue: gemmi.Residue):
        self.chain = chain
        self.residue = residue
        self.next = None
        self.prev = None
        self.score = None
        self.clashing = set()
        self.zone = None

    def __hash__(self) -> int:
        return hash(_key(self.chain, self.residue))


def _my_residues(structure: gemmi.Structure) -> dict:
    my_residues = {}
    for chain in structure[0]:
        previous = None
        for residue in chain:
            my_residue = _MyResidue(chain, residue)
            my_residues[_key(chain, residue)] = my_residue
            if previous is not None and _are_neighbours(previous, my_residue):
                my_residue.prev = previous
                previous.next = my_residue
            previous = my_residue
    return my_residues


def _are_neighbours(res1: _MyResidue, res2: _MyResidue) -> bool:
    res1 = res1.residue
    res2 = res2.residue
    info1 = gemmi.find_tabulated_residue(res1.name)
    info2 = gemmi.find_tabulated_residue(res2.name)
    return (
        info1.is_amino_acid()
        and info2.is_amino_acid()
        and "C" in res1
        and "N" in res2
        and res1["C"][0].pos.dist(res2["N"][0].pos) < 2
    ) or (
        info1.is_nucleic_acid()
        and info2.is_nucleic_acid()
        and "O3'" in res1
        and "P" in res2
        and res1["O3'"][0].pos.dist(res2["P"][0].pos) < 2
    )


def _scores(result: RefmacResult) -> dict:
    search = gemmi.NeighborSearch(result.structure, max_radius=5).populate()
    diff_density = result.fphi_diff.transform_f_phi_to_map(
        result.fphi_diff.label(0),
        result.fphi_diff.label(1),
        sample_rate=(result.fphi_diff.resolution_high() / 0.7),  # 0.7A grid spacing
    )
    sum_counts = {}
    for point in diff_density.masked_asu():
        position = diff_density.point_to_position(point)
        mark = search.find_nearest_atom(position)
        if mark is not None:
            cra = mark.to_cra(result.structure[0])
            key = _key(cra.chain, cra.residue)
            sum_counts.setdefault(key, [0, 0])
            if point.value > 0:
                sum_counts[key][0] += point.value
            sum_counts[key][1] += 1
    return {key: -sum / count for key, (sum, count) in sum_counts.items()}


def _assign_scores(residues: dict, scores: dict) -> None:
    for key, score in scores.items():
        residues[key].score = score


def _clashes(structure1: gemmi.Structure, structure2: gemmi.Structure) -> set:
    search = gemmi.NeighborSearch(structure2, max_radius=1).populate()
    clashes = set()
    for chain1 in structure1[0]:
        for residue1 in chain1:
            key1 = _key(chain1, residue1)
            for atom1 in residue1:
                for mark in search.find_atoms(atom1.pos, radius=1):
                    cra = mark.to_cra(structure2[0])
                    key2 = _key(cra.chain, cra.residue)
                    clashes.add((key1, key2))
    return clashes


def _assign_clashes(residues1: dict, residues2: dict, clashes: set) -> None:
    for key1, key2 in clashes:
        residues1[key1].clashing.add(residues2[key2])
        residues2[key2].clashing.add(residues1[key1])


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
