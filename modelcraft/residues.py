from functools import lru_cache
from typing import Set
import os
import gemmi


@lru_cache(maxsize=None)
def is_buffer(name: str) -> float:
    return name in _buffer_codes()


@lru_cache(maxsize=None)
def _buffer_codes() -> Set[str]:
    path = os.path.join(os.environ["CCP4"], "share", "pisa", "agents.dat")
    agents = set()
    with open(path) as stream:
        for line in stream:
            if line[0] != "#" and "," in line:
                code = line.split(",")[0]
                agents.add(code)
    return agents


@lru_cache(maxsize=None)
def weight(name: str) -> float:
    return sum(atom.el.weight for atom in _chemcomp(name).atoms)


@lru_cache(maxsize=None)
def volume(name: str) -> float:
    "18 cubic Angstroms per non-hydrogen"
    return sum(18 for atom in _chemcomp(name).atoms if not atom.is_hydrogen())


@lru_cache(maxsize=None)
def _chemcomp(name: str) -> gemmi.ChemComp:
    filename = name.upper() + ".cif"
    path = os.path.join(os.environ["CLIBD"], "monomers", name[0].lower(), filename)
    block = gemmi.cif.read(path)[-1]
    return gemmi.make_chemcomp_from_block(block)
