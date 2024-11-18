import functools
import os
import pathlib
import gemmi


@functools.cache
def _buffers() -> set:
    path = pathlib.Path(os.environ["CCP4"], "share", "pisa", "agents.dat")
    agents = {"UNX"}
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line[0] != "#" and "," in line:
                code = line.split(",")[0]
                agents.add(code)
    return agents


@functools.cache
def _path(code: str) -> pathlib.Path:
    directory = pathlib.Path(os.environ["CLIBD_MON"], code[0].lower())
    single = directory / f"{code.upper()}.cif"
    double = directory / f"{code.upper()}_{code.upper()}.cif"
    return double if double.exists() else single


@functools.cache
def atom_ids(code: str) -> set:
    return {atom.id for atom in chemcomp(code).atoms}


@functools.cache
def chemcomp(code: str) -> gemmi.ChemComp:
    doc = gemmi.cif.read(str(_path(code)))
    return gemmi.make_chemcomp_from_block(doc[-1])


@functools.cache
def in_library(code: str) -> bool:
    return _path(code).exists()


@functools.cache
def is_buffer(code: str) -> float:
    return code.upper() in _buffers()


@functools.cache
def volume(code: str) -> float:
    if in_library(code):
        return sum(18 for atom in chemcomp(code).atoms if not atom.is_hydrogen())
    return 0


@functools.cache
def weight(code: str) -> float:
    if in_library(code):
        return sum(atom.el.weight for atom in chemcomp(code).atoms)
    return 0


@functools.cache
def group(code: str) -> gemmi.ChemComp.Group:
    if in_library(code):
        doc = gemmi.cif.read(str(_path(code)))
        monlib = gemmi.MonLib()
        monlib.read_monomer_doc(doc)
        return monlib.monomers[code].group
    return None


@functools.cache
def is_protein(code: str) -> bool:
    return group(code) in {
        gemmi.ChemComp.Group.Peptide,
        gemmi.ChemComp.Group.PPeptide,
        gemmi.ChemComp.Group.MPeptide,
    }


@functools.cache
def is_nucleic(code: str) -> bool:
    return group(code) in {
        gemmi.ChemComp.Group.Dna,
        gemmi.ChemComp.Group.Rna,
        gemmi.ChemComp.Group.DnaRna,
    }
