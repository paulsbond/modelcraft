from ...jobs.libg import Libg
from ...structure import read_structure
from . import ccp4_path


def test_libg():
    pdb_path = ccp4_path("lib", "data", "nautilus_lib.pdb")
    structure = read_structure(pdb_path)
    Libg(structure).run()
