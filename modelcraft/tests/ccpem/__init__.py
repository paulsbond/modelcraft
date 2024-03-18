from os import environ
from pathlib import Path


def _test_data(path: str) -> str:
    return str(Path(environ["CCPEM"], "lib/py2/ccpem/src/ccpem_core/test_data", path))


def halfmap1_path() -> str:
    return _test_data("map/mrc/3488_run_half1_class001_unfil.mrc")


def halfmap2_path() -> str:
    return _test_data("map/mrc/3488_run_half2_class001_unfil.mrc")


def density_path() -> str:
    return _test_data("map/mrc/emd_3488.map")


def mask_path() -> str:
    return _test_data("map/mrc/emd_3488_mask.mrc")


def structure_path() -> str:
    return _test_data("pdb/5ni1.cif")


def sequence_path() -> str:
    return _test_data("seq/5ni1.fasta")
