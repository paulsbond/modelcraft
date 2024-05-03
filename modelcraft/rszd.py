"""
Real-space difference density Z-score (RSZD)
Implementation adapted from:
Tickle, I. J. (2012). Acta Cryst. D68, 454-467.
https://doi.org/10.1107/S0907444911035918
"""

import gemmi
from .reflections import DataItem


def rszd(structure: gemmi.Structure, fphi_diff: DataItem, model_index: int = 0) -> dict:
    dmin = fphi_diff.resolution_high()
    radii = {}
    for ci, chain in enumerate(structure[model_index]):
        for ri, residue in enumerate(chain):
            for ai, atom in enumerate(residue):
                radii[(ci, ri, ai)] = max_radius(dmin, atom)

    search_radius = max(radii.values())
    search = gemmi.NeighborSearch(structure, search_radius, model_index)
    search.populate(include_h=False)

    diff_den = fphi_diff.transform_f_phi_to_map(fphi_diff.label(0), fphi_diff.label(1))
    diff_den.normalize()  # TODO: Estimation of the standard uncertainty in Δρ

    # TODO: Calculate RSZD for main chain and side chain separately
    residue_values = {}
    for point in diff_den.masked_asu():
        position = diff_den.point_to_position(point)
        mark = search.find_nearest_atom(position)
        if mark is not None:
            distance = mark.pos.dist(position)  # TODO: Handle image_idx
            if distance < radii[(mark.chain_idx, mark.residue_idx, mark.atom_idx)]:
                cra = mark.to_cra(structure[model_index])
                key = (cra.chain.name, str(cra.residue.seqid))
                residue_values.setdefault(key, []).append(point.value)


def max_radius(dmin: float, atom: gemmi.Atom):
    return 1.5  # TODO: Calculate from dmin, element and B-factor
