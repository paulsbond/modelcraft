import gemmi
import numpy as np
from modelcraft.reflections import DataItem


def per_residue_rscc(
    structure: gemmi.Structure, fphi: DataItem, radius: float = 1.5
) -> dict:
    calculator = gemmi.DensityCalculatorX()
    calculator.d_min = fphi.resolution_high()
    calculator.set_grid_cell_and_spacegroup(structure)
    calculator.put_model_density_on_grid(structure[0])
    density = fphi.transform_f_phi_to_map(
        fphi.label(0), fphi.label(1), exact_size=calculator.grid.shape
    )
    search = gemmi.NeighborSearch(structure, max_radius=radius)
    search.populate(include_h=False)
    residue_pairs = {}
    for point in density.masked_asu():
        position = density.point_to_position(point)
        mark = search.find_nearest_atom(position)
        if mark is not None:
            cra = mark.to_cra(structure[0])
            key = (cra.chain.name, str(cra.residue.seqid))
            value1 = point.value
            value2 = calculator.grid.get_value(point.u, point.v, point.w)
            residue_pairs.setdefault(key, []).append((value1, value2))
    correlations = {}
    for key, pairs in residue_pairs.items():
        if len(pairs) > 1:
            values1, values2 = zip(*pairs)
            correlations[key] = np.corrcoef(values1, values2)[0, 1]
    return correlations
