import functools
import os
import gemmi


def _path(code: str) -> str:
    return os.path.join(os.environ["CLIBD_MON"], code[0].lower(), f"{code.upper()}.cif")


@functools.lru_cache(maxsize=None)
def atom_ids(code: str) -> set:
    return {atom.id for atom in chemcomp(code).atoms}


@functools.lru_cache(maxsize=None)
def chemcomp(code: str) -> gemmi.ChemComp:
    doc = gemmi.cif.read(_path(code))
    return gemmi.make_chemcomp_from_block(doc[-1])


@functools.lru_cache(maxsize=None)
def in_library(code: str) -> bool:
    return os.path.exists(_path(code))
