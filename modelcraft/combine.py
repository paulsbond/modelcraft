import gemmi
from .jobs.refmac import RefmacResult
from .monlib import is_protein, is_nucleic


def combine_results(buccaneer: RefmacResult, nautilus: RefmacResult) -> gemmi.Structure:
    structure = buccaneer.structure.clone()
    for i, chain in reversed(list(enumerate(structure[0]))):
        if _is_nucleic_chain(chain):
            del structure[0][i]
    chains_to_add, clashing_to_remove = _resolve_clashes(structure, buccaneer, nautilus)
    for chain in structure[0]:
        protein = _is_protein_chain(chain)
        any_removed = False
        for i, residue in reversed(list(enumerate(chain))):
            if _key(chain, residue) in clashing_to_remove:
                del chain[i]
                any_removed = True
        if protein and any_removed:
            _remove_isolated_fragments(chain, _are_joined_protein, min_length=6)
    structure.remove_empty_chains()
    for chain in chains_to_add:
        structure[0].add_chain(chain, unique_name=True)
    return structure


def _resolve_clashes(
    structure: gemmi.Structure, buccaneer: RefmacResult, nautilus: RefmacResult
) -> tuple:
    chains_to_add = []
    clashing_to_remove = set()
    bucc_scores = _scores(buccaneer)
    naut_scores = _scores(nautilus)
    search = gemmi.NeighborSearch(structure, max_radius=1).populate()
    for chain in nautilus.structure[0]:
        if _is_nucleic_chain(chain):
            residues_to_remove = set()
            for keys, clashing in _clashing_zones(chain, search, structure):
                if _mean_score(keys, naut_scores) < _mean_score(clashing, bucc_scores):
                    clashing_to_remove |= clashing
                else:
                    residues_to_remove |= keys
            if residues_to_remove:
                for i, residue in reversed(list(enumerate(chain))):
                    if _key(chain, residue) in residues_to_remove:
                        del chain[i]
                _remove_isolated_fragments(chain, _are_joined_nucleic, min_length=2)
            if len(chain) > 0:
                chains_to_add.append(chain)
    return chains_to_add, clashing_to_remove


def _scores(result: RefmacResult) -> dict:
    search = gemmi.NeighborSearch(result.structure, max_radius=5).populate()
    diff_density = result.fphi_diff.transform_f_phi_to_map(
        result.fphi_diff.label(0),
        result.fphi_diff.label(1),
        sample_rate=(result.fphi_diff.resolution_high() / 1.0),  # 1.0A grid spacing
    )
    sum_counts = {}
    for point in diff_density.masked_asu():
        position = diff_density.point_to_position(point)
        mark = search.find_nearest_atom(position)
        if mark is not None:
            cra = mark.to_cra(result.structure[0])
            key = _key(cra.chain, cra.residue)
            sum_counts.setdefault(key, [0, 0])
            sum_counts[key][0] += abs(point.value)
            sum_counts[key][1] += 1
    return {key: sum / count for key, (sum, count) in sum_counts.items()}


def _key(chain: gemmi.Chain, residue: gemmi.Residue) -> tuple:
    return (chain.name, residue.seqid.num, residue.seqid.icode)


def _is_nucleic_chain(chain: gemmi.Chain) -> bool:
    return len(chain) > 1 and all(is_nucleic(res.name) for res in chain)


def _is_protein_chain(chain: gemmi.Chain) -> bool:
    return len(chain) > 1 and all(is_protein(res.name) for res in chain)


def _clashing_zones(
    chain: gemmi.Chain, search: gemmi.NeighborSearch, structure: gemmi.Structure
):
    keys = set()
    clashing_keys = set()
    for residue in chain:
        residue_clashes = set()
        for atom in residue:
            for mark in search.find_atoms(atom.pos, radius=1):
                cra = mark.to_cra(structure[0])
                residue_clashes.add(_key(cra.chain, cra.residue))
        if residue_clashes:
            keys.add(_key(chain, residue))
            clashing_keys |= residue_clashes
        elif keys:
            yield keys, clashing_keys
            keys = set()
            clashing_keys = set()
    if keys:
        yield keys, clashing_keys


def _mean_score(keys: set, scores: dict) -> float:
    return sum(scores[key] for key in keys) / len(keys)


def _remove_isolated_fragments(chain: gemmi.Chain, are_joined, min_length: int):
    to_remove = []
    fragment = []
    for i, residue in enumerate(chain):
        if i > 0 and are_joined(chain[i - 1], residue):
            fragment.append(i)
        else:
            if len(fragment) < min_length:
                to_remove.extend(fragment)
            fragment = [i]
    if len(fragment) < min_length:
        to_remove.extend(fragment)
    for i in reversed(to_remove):
        del chain[i]


def _are_joined_protein(res1: gemmi.Residue, res2: gemmi.Residue) -> bool:
    return "C" in res1 and "N" in res2 and res1["C"][0].pos.dist(res2["N"][0].pos) < 2.5


def _are_joined_nucleic(res1: gemmi.Residue, res2: gemmi.Residue) -> bool:
    return (
        "O3'" in res1
        and "P" in res2
        and res1["O3'"][0].pos.dist(res2["P"][0].pos) < 2.5
    )
