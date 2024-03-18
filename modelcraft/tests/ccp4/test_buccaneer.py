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
    assert buccaneer.completeness_res > 0
    assert buccaneer.completeness_chn > 0
    assert buccaneer.chains_built > 0
    assert buccaneer.fragments_built > 0
    assert buccaneer.residues_built > 0
    assert buccaneer.residues_sequenced > 0
    assert buccaneer.residues_unique > 0
    assert buccaneer.longest_fragment > 0
    assert buccaneer.seconds > 0
    stats = ModelStats(buccaneer.structure)
    assert stats.residues > 0


def _known_structure_ids_test(pdb_id: str, expected: list):
    # Testing with PDB format so the function only has access to author IDs
    filename = f"pdb{pdb_id}.ent"
    pdbe_download(filename)
    structure = read_structure(filename)
    assert list(_known_structure_ids(structure)) == expected


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


@in_temp_directory
def test_4kui_known_structure_ids():
    "Modified residue at the start of the chain"
    _known_structure_ids_test("4kui", [])


@in_temp_directory
def test_5e1k_known_structure_ids():
    "Two METs next to each other and 5 calcium ions"
    _known_structure_ids_test(
        "5e1k",
        [
            "/A/201/*/:1.0",
            "/A/202/*/:1.0",
            "/A/203/*/:1.0",
            "/A/204/*/:1.0",
            "/A/205/*/:1.0",
        ],
    )
