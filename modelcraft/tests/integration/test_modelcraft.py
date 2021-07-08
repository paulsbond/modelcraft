import json
import os
import shutil
import uuid
import pytest
from modelcraft.scripts.modelcraft import main
from modelcraft.reflections import write_mtz
from modelcraft.structure import contains_residue, read_structure
from . import (
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
    contents.write_json_file("contents.json")
    args = []
    args += ["--data", "data.mtz"]
    args += ["--contents", "contents.json"]
    args += ["--cycles", "1"]
    with pytest.raises(SystemExit):
        main(args)
    with open("modelcraft.json") as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["termination_reason"] == "Normal"
    os.chdir("..")
    shutil.rmtree(tmp_dir)


def test_1rxf_from_model():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    hklin = ccp4_path("examples", "data", "1rxf.mtz")
    xyzin = ccp4_path("examples", "data", "1rxf_randomise.pdb")
    contents = pdb1rxf_contents()
    contents.write_sequence_file("sequence.fasta")
    args = []
    args += ["--data", hklin]
    args += ["--observations", "I,SIGI"]
    args += ["--freerflag", "FreeR_flag"]
    args += ["--contents", "sequence.fasta"]
    args += ["--model", xyzin]
    args += ["--cycles", "2"]
    with pytest.raises(SystemExit):
        main(args)
    with open("modelcraft.json") as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["termination_reason"] == "Normal"
    structure = read_structure("modelcraft.cif")
    assert contains_residue(structure, "FE")
    os.chdir("..")
    shutil.rmtree(tmp_dir)
