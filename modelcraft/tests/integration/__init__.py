import functools
import os
import shutil
import uuid
import urllib.request
import gemmi
from modelcraft.contents import AsuContents, Ligand, Polymer, PolymerType
from modelcraft.jobs.refmac import RefmacXray
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure


def ccp4_path(*paths: str) -> str:
    return os.path.join(os.environ["CCP4"], *paths)


def i2_demo_path(*paths: str) -> str:
    ccp4_7_dir = ccp4_path("share", "ccp4i2", "demo_data")
    ccp4_8_dir = ccp4_path("lib", "python3.7", "site-packages", "ccp4i2", "demo_data")
    if os.path.exists(ccp4_7_dir):
        return os.path.join(ccp4_7_dir, *paths)
    return os.path.join(ccp4_8_dir, *paths)


def in_temp_directory(func):
    def wrapper():
        tmp_dir = "tmp%s" % uuid.uuid4()
        os.mkdir(tmp_dir)
        os.chdir(tmp_dir)
        try:
            func()
        except Exception as exception:
            os.chdir("..")
            raise exception
        os.chdir("..")
        shutil.rmtree(tmp_dir)

    return wrapper


def pdbe_download(filename: str) -> None:
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    urllib.request.urlretrieve(url, filename)


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
    refmac = RefmacXray(structure=structure, fsigf=fsigf, freer=freer, cycles=0)
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
