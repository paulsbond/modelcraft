import os
import shutil
import uuid
import gemmi
from modelcraft.model import write_mmcif
from modelcraft.tests import data_path


def test_read_write_structure():
    directory = "test_%s" % uuid.uuid4()
    os.mkdir(directory)
    os.chdir(directory)
    structure = gemmi.read_structure(data_path("1kv9_model.pdb"))
    write_mmcif("structure.cif", structure)
    structure = gemmi.read_structure("structure.cif")
    os.chdir("..")
    shutil.rmtree(directory)
