"""
Real-space difference density Z-score (RSZD).
Implementation adapted from:
Tickle, I. J. (2012). Acta Cryst. D68, 454-467.
https://doi.org/10.1107/S0907444911035918
"""

import gemmi
import scipy.stats
from .reflections import DataItem


def per_residue_rszd(
    structure: gemmi.Structure, fphi_diff: DataItem, model_index: int = 0
) -> dict:
    dmin = fphi_diff.resolution_high()
    radii = {}
    for ci, chain in enumerate(structure[model_index]):
        for ri, residue in enumerate(chain):
            for ai, atom in enumerate(residue):
                radii[(ci, ri, ai, atom.altloc)] = max_radius(dmin, atom)

    search_radius = max(radii.values())
    search = gemmi.NeighborSearch(structure, search_radius, model_index)
    search.populate(include_h=False)

    diff_den = fphi_diff.transform_f_phi_to_map(fphi_diff.label(0), fphi_diff.label(1))
    diff_den.normalize()  # TODO: Estimation of the standard uncertainty in Δρ

    # TODO: Calculate RSZD for main chain and side chain separately
    values_dict = {}
    for point in diff_den.masked_asu():
        position = diff_den.point_to_position(point)
        mark = search.find_nearest_atom(position, search_radius)
        if mark is not None:
            nearest = structure.cell.find_nearest_pbc_image(position, mark.pos, 0)
            key = (mark.chain_idx, mark.residue_idx, mark.atom_idx, mark.altloc)
            if nearest.dist() < radii[key]:
                cra = mark.to_cra(structure[model_index])
                key = (cra.chain.name, str(cra.residue.seqid))
                values_dict.setdefault(key, []).append(point.value)

    rszd_dict = {}
    for key, values in values_dict.items():
        # TODO: Statistically independent difference density values from resampling
        s = sum(x * x for x in values)
        n = len(values)
        p = 1 - scipy.stats.chi2.cdf(s, n)
        z = abs(scipy.stats.norm.ppf(p / 2))
        rszd_dict[key] = z
    return rszd_dict


def max_radius(dmin: float, atom: gemmi.Atom):
    return 1.5  # TODO: Calculate from dmin, element and B-factor
