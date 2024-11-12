import collections
import os
import gemmi
import numpy as np


def topology(
    structure: gemmi.Structure, model_index: int = 0, libin: str = ""
) -> gemmi.Topo:
    monlib = gemmi.read_monomer_lib(
        monomer_dir=os.environ["CLIBD_MON"],
        resnames=structure[model_index].get_all_residue_names(),
        libin=libin or "",
    )
    return gemmi.prepare_topology(
        st=structure,
        monlib=monlib,
        model_index=model_index,
    )


def per_residue_geometry_rmsz(
    structure: gemmi.Structure, model_index: int = 0, libin: str = ""
) -> dict:
    atom_zs = _atom_zs(structure, model_index, libin)
    rv = {}
    for chain in structure[model_index]:
        for residue in chain:
            zs = np.concatenate([atom_zs.get(atom.serial, []) for atom in residue])
            rmsz = np.sqrt(np.mean(np.square(zs))) if len(zs) > 0 else np.nan
            rv[(chain.name, str(residue.seqid))] = rmsz
    return rv


def _atom_zs(structure: gemmi.Structure, model_index, libin) -> dict:
    structure.assign_serial_numbers()
    topo = topology(structure, model_index, libin)
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
        if torsion.restr.esd > 0:
            z = torsion.calculate_z()
            for atom in torsion.atoms:
                atom_zs[atom.serial].append(z)
    for plane in topo.planes:
        best_plane = gemmi.find_best_plane(plane.atoms)
        for atom in plane.atoms:
            z = gemmi.get_distance_from_plane(atom.pos, best_plane) / plane.restr.esd
            atom_zs[atom.serial].append(z)
    return atom_zs
