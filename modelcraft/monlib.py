import functools
import os
import pathlib
import gemmi


@functools.lru_cache(maxsize=None)
def _buffers() -> set:
    path = pathlib.Path(os.environ["CCP4"], "share", "pisa", "agents.dat")
    agents = set()
    with path.open() as stream:
        for line in stream:
            if line[0] != "#" and "," in line:
                code = line.split(",")[0]
                agents.add(code)
    return agents


@functools.lru_cache(maxsize=None)
def _path(code: str) -> pathlib.Path:
    directory = pathlib.Path(os.environ["CLIBD_MON"], code[0].lower())
    single = directory / f"{code.upper()}.cif"
    double = directory / f"{code.upper()}_{code.upper()}.cif"
    return double if double.exists() else single


@functools.lru_cache(maxsize=None)
def atom_ids(code: str) -> set:
    return {atom.id for atom in chemcomp(code).atoms}


@functools.lru_cache(maxsize=None)
def chemcomp(code: str) -> gemmi.ChemComp:
    doc = gemmi.cif.read(str(_path(code)))
    return gemmi.make_chemcomp_from_block(doc[-1])


@functools.lru_cache(maxsize=None)
def in_library(code: str) -> bool:
    return _path(code).exists()


@functools.lru_cache(maxsize=None)
def is_buffer(code: str) -> float:
    return code.upper() in _buffers()


@functools.lru_cache(maxsize=None)
def volume(code: str) -> float:
    return sum(18 for atom in chemcomp(code).atoms if not atom.is_hydrogen())


@functools.lru_cache(maxsize=None)
def weight(code: str) -> float:
    return sum(atom.el.weight for atom in chemcomp(code).atoms)
