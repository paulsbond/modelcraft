from collections import defaultdict

import gemmi
import numpy as np
import pandas as pd

from .geometry import per_residue_geometry_rmsz
from .jobs.refmac import RefmacResult
from .monlib import MonLib
from .reflections import DataItem
from .utils import modified_zscore


def validate(
    structure: gemmi.Structure,
    fphi_best: DataItem,
    fphi_diff: DataItem,
    fphi_calc: DataItem,
    monlib: MonLib,
    model_index: int = 0,
) -> pd.DataFrame:
    best_map = fphi_best.map(spacing=1.0)
    diff_map = fphi_diff.map(size=best_map.shape)
    calc_map = fphi_calc.map(size=best_map.shape)

    bfac = _bfac(structure[model_index])
    rscc, diff = _rscc_diff(structure, best_map, diff_map, calc_map, model_index)
    geom = per_residue_geometry_rmsz(structure, monlib, model_index)

    data = {
        "Chain": [],
        "SeqId": [],
        "Name": [],
        "BFac": [],
        "RSCC": [],
        "Diff": [],
        "Geom": [],
    }
    for chain in structure[model_index]:
        for residue in chain:
            if monlib.is_protein(residue.name):
                key = (chain.name, str(residue.seqid))
                data["Chain"].append(chain.name)
                data["SeqId"].append(str(residue.seqid))
                data["Name"].append(residue.name)
                data["BFac"].append(bfac[key])
                data["RSCC"].append(rscc[key])
                data["Diff"].append(diff[key])
                data["Geom"].append(geom[key])

    df = pd.DataFrame(data)
    df["BFac"] = -modified_zscore(df["BFac"])
    df["RSCC"] = modified_zscore(df["RSCC"])
    df["Diff"] = -modified_zscore(df["Diff"])
    df["Geom"] = -modified_zscore(df["Geom"])
    density_score = (df["BFac"] + df["RSCC"] + df["Diff"]) / 3
    df["Score"] = modified_zscore((density_score + df["Geom"]) / 2)
    return df


def validate_refmac(result: RefmacResult, monlib: MonLib) -> pd.DataFrame:
    return validate(
        result.structure, result.fphi_best, result.fphi_diff, result.fphi_calc, monlib
    )


def _bfac(model: gemmi.Model) -> dict:
    return {
        (chain.name, str(residue.seqid)): np.mean([a.b_iso for a in residue])
        for chain in model
        for residue in chain
    }


def _rscc_diff(
    structure: gemmi.Structure,
    best_map: gemmi.FloatGrid,
    diff_map: gemmi.FloatGrid,
    calc_map: gemmi.FloatGrid,
    model_index: int,
) -> dict:
    search = gemmi.NeighborSearch(structure, max_radius=3, model_index=model_index)
    search.populate(include_h=False)
    best_values = defaultdict(list)
    diff_values = defaultdict(list)
    calc_values = defaultdict(list)
    for point in best_map.masked_asu():
        position = best_map.point_to_position(point)
        mark = search.find_nearest_atom(position, radius=3)
        if mark is not None:
            cra = mark.to_cra(structure[model_index])
            key = (cra.chain.name, str(cra.residue.seqid))
            best_values[key].append(point.value)
            diff_values[key].append(diff_map.get_value(point.u, point.v, point.w))
            calc_values[key].append(calc_map.get_value(point.u, point.v, point.w))
    rscc = {}
    diff = {}
    for key, best_value in best_values.items():
        rscc[key] = np.corrcoef(best_value, calc_values[key])[0, 1]
        diff[key] = np.sqrt(np.mean(np.square(diff_values[key])))
    return rscc, diff
