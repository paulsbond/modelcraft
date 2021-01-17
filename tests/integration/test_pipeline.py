import json
import os
import pytest
import shutil
import uuid
from modelcraft.__main__ import main
from modelcraft.contents import PolymerType
from modelcraft.reflections import write_mtz
from modelcraft.structure import contains_residue, read_structure
from tests.integration import (
    ccp4_path,
    insulin_contents,
    insulin_fsigf,
    insulin_freer,
    insulin_refmac,
    pdb1rxf_contents,
)


def test_insulin_from_phases():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    fsigf = insulin_fsigf()
    freer = insulin_freer()
    refmac = insulin_refmac()
    write_mtz("data.mtz", [fsigf, freer, refmac.abcd])
    contents = insulin_contents()
    contents.write_sequence_file("sequence.fasta")
    args = []
    args += ["--hklin", "data.mtz"]
    args += ["--seqin", "sequence.fasta"]
    args += ["--cycles", "1"]
    with pytest.raises(SystemExit):
        main(args)
    with open("modelcraft.json") as report_file:
        report = json.load(report_file)
    assert report["real_time"]["total"] > 0
    assert report["termination_reason"] == "Normal"
    os.chdir("..")
    shutil.rmtree(tmp_dir)


def test_1rxf_from_model():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    hklin = ccp4_path("examples", "data", "1rxf.mtz")
    xyzin = ccp4_path("examples", "data", "1rxf_randomise.pdb")
    structure = read_structure(xyzin)
    structure.remove_alternative_conformations()
    structure.write_minimal_pdb("model.pdb")
    contents = pdb1rxf_contents()
    contents.write_sequence_file("1rxf.fasta", polymer_type=PolymerType.PROTEIN)
    args = []
    args += ["--hklin", hklin]
    args += ["--amplitudes", "F,SIGF"]
    args += ["--freerflag", "FreeR_flag"]
    args += ["--seqin", "1rxf.fasta"]
    args += ["--model", "model.pdb"]
    args += ["--cycles", "2"]
    with pytest.raises(SystemExit):
        main(args)
    with open("modelcraft.json") as report_file:
        report = json.load(report_file)
    assert report["real_time"]["total"] > 0
    assert report["termination_reason"] == "Normal"
    structure = read_structure("modelcraft.cif")
    assert contains_residue(structure, "FE")
    os.chdir("..")
    shutil.rmtree(tmp_dir)
