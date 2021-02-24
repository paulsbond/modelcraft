from functools import lru_cache
from math import pi
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
def weight(name: str, include_hydrogens: bool = True) -> float:
    weight = 0
    for atom in _chemcomp(name).atoms:
        if include_hydrogens or not atom.is_hydrogen():
            weight += atom.el.weight
    return weight


@lru_cache(maxsize=None)
def vdw_volume(name: str, include_hydrogrens: bool = False) -> float:
    volume = 0
    for atom in _chemcomp(name).atoms:
        if include_hydrogrens or not atom.is_hydrogen():
            volume += 4 / 3 * pi * atom.el.vdw_r ** 3
    return volume


@lru_cache(maxsize=None)
def _chemcomp(name: str) -> gemmi.ChemComp:
    filename = name.upper() + ".cif"
    path = os.path.join(os.environ["CLIBD"], "monomers", name[0].lower(), filename)
    block = gemmi.cif.read(path)[-1]
    return gemmi.make_chemcomp_from_block(block)
