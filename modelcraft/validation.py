from os import environ
import gemmi
import numpy as np
from .jobs.refmac import RefmacResult
from .utils import modified_zscore


def validate(result: RefmacResult) -> dict:
    den_score = density_score(result)
    geom_score = geometry_score(result)
    return {key: (den_score[key] + geom_score[key]) / 2 for key in den_score.keys()}


def density_score(result: RefmacResult) -> dict:
    """
    Score residues based on their fit to the density.
    A higher score is better.
    """
    best_map = result.fphi_best.map(spacing=1.0)
    diff_map = result.fphi_diff.map(size=best_map.shape)
    calc_map = result.fphi_calc.map(size=best_map.shape)

    best_score = {}
    for chain in result.structure[0]:
        for residue in chain:
            key = residue_key(chain, residue)
            best_score[key] = density_fit(residue, best_map)

    search = gemmi.NeighborSearch(result.structure, max_radius=3)
    search.populate(include_h=False)

    best_values = {}
    diff_values = {}
    calc_values = {}
    for point in best_map.masked_asu():
        position = best_map.point_to_position(point)
        mark = search.find_nearest_atom(position, radius=3)
        if mark is not None:
            cra = mark.to_cra(result.structure[0])
            key = residue_key(cra.chain, cra.residue)
            best_value = point.value
            diff_value = diff_map.get_value(point.u, point.v, point.w)
            calc_value = calc_map.get_value(point.u, point.v, point.w)
            best_values.setdefault(key, []).append(best_value)
            diff_values.setdefault(key, []).append(diff_value)
            calc_values.setdefault(key, []).append(calc_value)

    rscc_score = {}
    diff_score = {}
    for key in best_values.keys():
        rscc_score[key] = np.corrcoef(best_values[key], calc_values[key])[0, 1]
        diff_score[key] = -np.mean(np.abs(diff_values[key]))

    rscc_score = _dictmod(rscc_score, modified_zscore)
    diff_score = _dictmod(diff_score, modified_zscore)
    best_score = _dictmod(best_score, modified_zscore)
    best_score = _dictmod(best_score, lambda values: -np.abs(values))

    return {
        key: (rscc_score[key] + diff_score[key] + best_score[key]) / 3
        for key in rscc_score.keys()
    }


def geometry_score(structure: gemmi.Structure) -> dict:
    """
    Calculates Z-scores for bond lengths, angles, torsions and planes.
    Returns the absolue negative of the worse Z-score for each residue
    (i.e. 0 is the best score and more negative scores are worse).
    """
    structure.assign_serial_numbers()
    residue_names = structure[0].get_all_residue_names()
    monlib = gemmi.read_monomer_lib(environ["CLIBD_MON"], residue_names)
    topo = gemmi.prepare_topology(structure, monlib)
    atom_scores = {}
    for bond in topo.bonds:
        z = abs(bond.calculate_z())
        for atom in bond.atoms:
            atom_scores[atom.serial] = max(atom_scores.get(atom.serial, 0), z)
    for angle in topo.angles:
        z = abs(angle.calculate_z())
        for atom in angle.atoms:
            atom_scores[atom.serial] = max(atom_scores.get(atom.serial, 0), z)
    for torsion in topo.torsions:
        z = abs(torsion.calculate_z())
        for atom in torsion.atoms:
            atom_scores[atom.serial] = max(atom_scores.get(atom.serial, 0), z)
    for plane in topo.planes:
        best_plane = gemmi.find_best_plane(plane.atoms)
        for atom in plane.atoms:
            z = gemmi.get_distance_from_plane(atom.pos, best_plane) / plane.restr.esd
            atom_scores[atom.serial] = max(atom_scores.get(atom.serial, 0), z)
    residue_scores = {}
    for chain in structure[0]:
        for residue in chain:
            worst = max(atom_scores.get(atom.serial, 0) for atom in residue)
            residue_scores[residue_key(chain, residue)] = -worst
    return residue_scores


def residue_key(chain: gemmi.Chain, residue: gemmi.Residue) -> tuple:
    return (chain.name, residue.seqid.num, residue.seqid.icode)


def density_fit(residue: gemmi.Residue, density: gemmi.FloatGrid) -> float:
    """
    Quick density fit score using the mean value at atomic centres
    (weighted by the atomic number of the atom).
    """
    sum = 0
    count = 0
    for atom in residue:
        value = density.interpolate_value(atom.pos)
        sum += value / atom.element.atomic_number
        count += 1
    return sum / count


def _dictmod(d: dict, f) -> dict:
    return {k: v for k, v in zip(d.keys(), f(d.values()))}
