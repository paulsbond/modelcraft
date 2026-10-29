import collections

import gemmi
import numpy as np

from .monlib import MonLib


def per_residue_geometry_rmsz(
    structure: gemmi.Structure, monlib: MonLib, model_index: int = 0
) -> dict:
    atom_zs = _atom_zs(structure, monlib, model_index)
    rv = {}
    for chain in structure[model_index]:
        for residue in chain:
            zs = np.concatenate([atom_zs.get(atom.serial, []) for atom in residue])
            rmsz = np.sqrt(np.mean(np.square(zs))) if len(zs) > 0 else np.nan
            rv[(chain.name, str(residue.seqid))] = rmsz
    return rv


def _atom_zs(structure: gemmi.Structure, monlib: MonLib, model_index: int) -> dict:
    structure.assign_serial_numbers()
    topo = gemmi.prepare_topology(structure, monlib, model_index)
    atom_zs = collections.defaultdict(list)
    for bond in topo.bonds:
        z = bond.calculate_z()
        for atom in bond.atoms:
            atom_zs[atom.serial].append(z)
    for angle in topo.angles:
        z = angle.calculate_z()
        for atom in angle.atoms:
            atom_zs[atom.serial].append(z)
    for torsion in topo.torsions:
        if torsion.restr.esd > 0:  # Some torsions are only restrained by planes
            z = torsion.calculate_z()
            for atom in torsion.atoms:
                atom_zs[atom.serial].append(z)
    for plane in topo.planes:
        best_plane = gemmi.find_best_plane(plane.atoms)
        for atom in plane.atoms:
            z = gemmi.get_distance_from_plane(atom.pos, best_plane) / plane.restr.esd
            atom_zs[atom.serial].append(z)
    return atom_zs
