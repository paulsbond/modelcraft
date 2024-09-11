import functools
import os
import gemmi


@functools.lru_cache(maxsize=None)
def _path(code: str) -> str:
    directory = os.path.join(os.environ["CLIBD_MON"], code[0].lower())
    single = os.path.join(directory, f"{code.upper()}.cif")
    double = os.path.join(directory, f"{code.upper()}_{code.upper()}.cif")
    return double if os.path.exists(double) else single


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


@functools.lru_cache(maxsize=None)
def group(code: str) -> gemmi.ChemComp.Group:
    if in_library(code):
        doc = gemmi.cif.read(_path(code))
        monlib = gemmi.MonLib()
        monlib.read_monomer_doc(doc)
        return monlib.monomers[code].group
    return None


@functools.lru_cache(maxsize=None)
def is_protein(code: str) -> bool:
    return group(code) in {
        gemmi.ChemComp.Group.Peptide,
        gemmi.ChemComp.Group.PPeptide,
        gemmi.ChemComp.Group.MPeptide,
    }


@functools.lru_cache(maxsize=None)
def is_nucleic(code: str) -> bool:
    return group(code) in {
        gemmi.ChemComp.Group.Dna,
        gemmi.ChemComp.Group.Rna,
        gemmi.ChemComp.Group.DnaRna,
    }
