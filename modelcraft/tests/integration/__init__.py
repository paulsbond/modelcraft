import functools
import os
import gemmi
from modelcraft.contents import AsuContents, Ligand, Polymer, PolymerType
from modelcraft.jobs.refmac import Refmac
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure


def ccp4_path(*paths: str) -> str:
    return os.path.join(os.environ["CCP4"], *paths)


@functools.lru_cache(maxsize=None)
def insulin_fsigf():
    mtz_path = ccp4_path("examples", "data", "insulin.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    return DataItem(mtz, "F,SIGF")


@functools.lru_cache(maxsize=None)
def insulin_freer():
    mtz_path = ccp4_path("examples", "data", "insulin.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    return DataItem(mtz, "FreeR_flag")


@functools.lru_cache(maxsize=None)
def insulin_structure():
    pdb_path = ccp4_path("examples", "data", "insulin.pdb")
    return read_structure(pdb_path)


@functools.lru_cache(maxsize=None)
def insulin_refmac():
    fsigf = insulin_fsigf()
    freer = insulin_freer()
    structure = insulin_structure()
    refmac = Refmac(structure=structure, fsigf=fsigf, freer=freer, cycles=0)
    return refmac.run()


@functools.lru_cache(maxsize=None)
def insulin_contents():
    contents = AsuContents()
    chain_a = Polymer(
        sequence="GIVEQCCASVCSLYQLENYCN",
        polymer_type=PolymerType.PROTEIN,
    )
    chain_b = Polymer(
        sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKA",
        polymer_type=PolymerType.PROTEIN,
    )
    contents.add_polymer(chain_a)
    contents.add_polymer(chain_b)
    return contents


@functools.lru_cache(maxsize=None)
def pdb1rxf_contents():
    sequence = (
        "MDTTVPTFSLAELQQGLHQDEFRRCLRDKGLFYLTDCGLTDTELKSAKDLVIDFFEHGSE"
        "AEKRAVTSPVPTMRRGFTGLESESTAQITNTGSYSDYSMCYSMGTADNLFPSGDFERIWT"
        "QYFDRQYTASRAVAREVLRATGTEPDGGVEAFLDCEPLLRFRYFPQVPEHRSAEEQPLRM"
        "APHYDLSMVTLIQQTPCANGFVSLQAEVGGAFTDLPYRPDAVLVFCGAIATLVTGGQVKA"
        "PRHHVAAPRRDQIAGSSRTSSVFFLRPNADFTFSVPLARECGFDVSLDGETATFQDWIGG"
        "NYVNIRRTSKA"
    )
    protein = Polymer(sequence=sequence, polymer_type=PolymerType.PROTEIN)
    ligand = Ligand(code="FE")
    contents = AsuContents()
    contents.proteins.append(protein)
    contents.ligands.append(ligand)
    contents.copies = 1
    return contents
