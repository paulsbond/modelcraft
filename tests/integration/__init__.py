import functools
import os
import shutil
import gemmi
from modelcraft.contents import AsuContents, Polymer, PolymerType
from modelcraft.jobs import Refmac
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure


def ccp4_path(*paths: str) -> str:
    if "CCP4" not in os.environ:
        raise EnvironmentError("CCP4 environment not set")
    return os.path.join(os.environ["CCP4"], *paths)


def remove_logs():
    shutil.rmtree("modelcraft-logs", ignore_errors=True)


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
    return Refmac(structure, fsigf, freer, cycles=0)


@functools.lru_cache(maxsize=None)
def insulin_contents():
    contents = AsuContents()
    a = Polymer(
        sequence="GIVEQCCASVCSLYQLENYCN",
        label="A",
        copies=1,
        polymer_type=PolymerType.PROTEIN,
    )
    b = Polymer(
        sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKA",
        label="B",
        copies=1,
        polymer_type=PolymerType.PROTEIN,
    )
    contents.polymers = [a, b]
    return contents


@functools.lru_cache(maxsize=None)
def pdb1rxf_contents():
    contents = AsuContents()
    a = Polymer(
        sequence=(
            "MDTTVPTFSLAELQQGLHQDEFRRCLRDKGLFYLTDCGLTDTELKSAKDLVIDFFEHGSE"
            "AEKRAVTSPVPTMRRGFTGLESESTAQITNTGSYSDYSMCYSMGTADNLFPSGDFERIWT"
            "QYFDRQYTASRAVAREVLRATGTEPDGGVEAFLDCEPLLRFRYFPQVPEHRSAEEQPLRM"
            "APHYDLSMVTLIQQTPCANGFVSLQAEVGGAFTDLPYRPDAVLVFCGAIATLVTGGQVKA"
            "PRHHVAAPRRDQIAGSSRTSSVFFLRPNADFTFSVPLARECGFDVSLDGETATFQDWIGG"
            "NYVNIRRTSKA"
        ),
        label="A",
        copies=1,
        polymer_type=PolymerType.PROTEIN,
    )
    contents.polymers = [a]
    return contents
