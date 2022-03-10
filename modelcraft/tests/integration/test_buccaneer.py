from modelcraft.jobs.buccaneer import Buccaneer, _known_structure_ids
from modelcraft.structure import ModelStats, read_structure
from . import (
    in_temp_directory,
    insulin_fsigf,
    insulin_freer,
    insulin_refmac,
    insulin_contents,
    pdbe_download,
)


def test_insulin():
    fsigf = insulin_fsigf()
    freer = insulin_freer()
    refmac = insulin_refmac()
    contents = insulin_contents()
    buccaneer = Buccaneer(
        contents=contents,
        fsigf=fsigf,
        phases=refmac.abcd,
        freer=freer,
        cycles=1,
    ).run()
    stats = ModelStats(buccaneer.structure)
    assert stats.residues > 0


def _known_structure_ids_test(pdb_id: str, expected: list):
    # Testing with PDB format so the function only has access to author IDs
    filename = f"pdb{pdb_id}.ent"
    pdbe_download(filename)
    structure = read_structure(filename)
    assert _known_structure_ids(structure) == expected


@in_temp_directory
def test_1o6a_known_structure_ids():
    "Modified MSE residues but no ligands"
    _known_structure_ids_test("1o6a", [])


@in_temp_directory
def test_1bkx_known_structure_ids():
    "Modified SEP and TPO residues but no ligands"
    _known_structure_ids_test("1bkx", [])


@in_temp_directory
def test_5oah_known_structure_ids():
    "95B ligands (with CA atoms) in each chain"
    _known_structure_ids_test("5oah", ["/A/402/*/:1.0", "/B/402/*/:1.0"])
