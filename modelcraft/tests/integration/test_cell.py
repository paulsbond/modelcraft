import urllib.request
import gemmi
from modelcraft.cell import max_distortion, remove_scale, update_cell
from modelcraft.structure import read_structure
from . import in_temp_directory


@in_temp_directory
def test_1ana():
    url = "https://ftp.ebi.ac.uk/pub/databases/pdb_versioned/data/entries/"
    url += "an/pdb_00001ana/pdb_00001ana_xyz_v1-2.cif.gz"
    urllib.request.urlretrieve(url, "1ana.cif.gz")
    structure = read_structure("1ana.cif.gz")
    assert structure.cell.parameters == (41.1, 41.1, 26.7, 90, 90, 90)
    assert structure.cell.fractionalization_matrix.tolist()[0][0] == 0
    remove_scale(structure)
    assert structure.cell.fractionalization_matrix.tolist()[0][0] != 0
    new_parameters = (41, 41, 27, 90, 90, 90)
    new_cell = gemmi.UnitCell(*new_parameters)
    assert max_distortion(structure.cell, new_cell) < 0.05
    old_pos = list(structure[0][0][0][0].pos)
    update_cell(structure, new_cell)
    assert structure.cell.parameters == new_parameters
    new_pos = list(structure[0][0][0][0].pos)
    assert old_pos != new_pos
