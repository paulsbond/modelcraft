import math
import os
import gemmi
from modelcraft.monlib import in_library


def rmsz(structure: gemmi.Structure) -> float:
    codes = [code for code in structure[0].get_all_residue_names() if in_library(code)]
    monlib = gemmi.read_monomer_lib(os.environ["CLIBD_MON"], codes)
    devnull = open(os.devnull, "w")
    topology = gemmi.prepare_topology(structure, monlib, warnings=devnull)
    num_of_squares = 0
    sum_of_squares = 0.0
    for bond in topology.bonds:
        num_of_squares += 1
        sum_of_squares += bond.calculate_z() ** 2
    for angle in topology.angles:
        num_of_squares += 1
        sum_of_squares += angle.calculate_z() ** 2
    for torsion in topology.torsions:
        if torsion.restr.esd > 0:
            num_of_squares += 1
            sum_of_squares += torsion.calculate_z() ** 2
    for plane in topology.planes:
        best_plane = gemmi.find_best_plane(plane.atoms)
        max_z = 0
        for atom in plane.atoms:
            distance = gemmi.get_distance_from_plane(atom.pos, best_plane)
            max_z = max(distance / plane.restr.esd, max_z)
        num_of_squares += 1
        sum_of_squares += max_z ** 2
    return math.sqrt(sum_of_squares / num_of_squares)
