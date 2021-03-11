import functools
import os
import gemmi


class Monomers:
    def __init__(self, library: str = None):
        self._chemcomps = {}
        self._library = library or os.path.join(os.environ("CLIBD"), "monomers")

    def add(self, code: str, chemcomp: gemmi.ChemComp):
        self._chemcomps[code.upper()] = chemcomp

    def weight(self, code: str) -> float:
        return _weight(self.chemcomp(code))

    def volume(self, code: str) -> float:
        return _volume(self.chemcomp(code))

    def chemcomp(self, code: str) -> gemmi.ChemComp:
        code = code.upper()
        if code in self._chemcomps:
            return self._chemcomps[code]
        path = os.path.join(self._library, code[0].lower(), code + ".cif")
        if not os.path.exists(path):
            raise RuntimeError(f"Monomer {code} not found at {path}")
        block = gemmi.cif.read(path)[-1]
        chemcomp = gemmi.make_chemcomp_from_block(block)
        self._chemcomps[code] = chemcomp
        return chemcomp


@functools.lru_cache(maxsize=None)
def is_buffer(code: str) -> float:
    return code.upper() in _buffers()


@functools.lru_cache(maxsize=None)
def _weight(chemcomp: gemmi.ChemComp) -> float:
    return sum(atom.el.weight for atom in chemcomp.atoms)


@functools.lru_cache(maxsize=None)
def _volume(chemcomp: gemmi.ChemComp) -> float:
    return sum(18 for atom in chemcomp if not atom.is_hydrogen())


@functools.lru_cache(maxsize=None)
def _buffers() -> set:
    path = os.path.join(os.environ["CCP4"], "share", "pisa", "agents.dat")
    agents = set()
    with open(path) as stream:
        for line in stream:
            if line[0] != "#" and "," in line:
                code = line.split(",")[0]
                agents.add(code)
    return agents
