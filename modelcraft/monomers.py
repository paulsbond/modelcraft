import functools
import os
import gemmi


class MonomerLibrary:
    def __init__(self, path: str = None):
        self._path = path or os.path.join(os.environ["CLIBD"], "monomers")
        self._chemcomps = {}

    def add(self, code: str, chemcomp: gemmi.ChemComp):
        self._chemcomps[code.upper()] = chemcomp

    def chemcomp(self, code: str) -> gemmi.ChemComp:
        code = code.upper()
        if code in self._chemcomps:
            return self._chemcomps[code]
        path = os.path.join(self._path, code[0].lower(), code + ".cif")
        if not os.path.exists(path):
            raise RuntimeError(f"Monomer {code} not found at {path}")
        block = gemmi.cif.read(path)[-1]
        chemcomp = gemmi.make_chemcomp_from_block(block)
        self._chemcomps[code] = chemcomp
        return chemcomp

    def weight(self, code: str) -> float:
        return _weight(self.chemcomp(code))

    def volume(self, code: str) -> float:
        return _volume(self.chemcomp(code))


# def weight(self, modified=False) -> float:
#     codes = self.residue_codes(modified=modified)
#     total = sum(self.monomers.weight(code) for code in codes)
#     total -= self.monomers.weight("HOH") * (len(codes) - 1)
#     return total

# def volume(self) -> float:
#     density = 1.35 if self.type == PolymerType.PROTEIN else 2.0
#     return self.weight(modified=False) / (density * 0.602214)

# def volume(self) -> float:
#     length = 0
#     total = 0
#     for code, copies in self.codes.items():
#         length += copies
#         total += monomers.volume(code) * copies
#     total -= monomers.volume("HOH") * length
#     return total

# def volume(self) -> float:
#     total = 0
#     for item in self.proteins + self.rnas + self.dnas + self.carbs + self.ligands:
#         total += item.volume() * (item.copies or 1)
#     return total


@functools.lru_cache(maxsize=None)
def _weight(chemcomp: gemmi.ChemComp) -> float:
    return sum(atom.el.weight for atom in chemcomp.atoms)


@functools.lru_cache(maxsize=None)
def _volume(chemcomp: gemmi.ChemComp) -> float:
    return sum(18 for atom in chemcomp if not atom.is_hydrogen())
