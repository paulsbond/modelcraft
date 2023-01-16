import json
import os
import pytest
from modelcraft.scripts.modelcraft import main
from modelcraft.reflections import write_mtz
from modelcraft.structure import contains_residue, read_structure, write_mmcif
from . import (
    ccp4_path,
    in_temp_directory,
    insulin_contents,
    insulin_fsigf,
    insulin_freer,
    insulin_refmac,
    pdb1rxf_contents,
)


@in_temp_directory
def test_insulin_from_phases():
    fsigf = insulin_fsigf()
    freer = insulin_freer()
    refmac = insulin_refmac()
    write_mtz("data.mtz", [fsigf, freer, refmac.abcd])
    contents = insulin_contents()
    contents.write_json_file("contents.json")
    args = ["xray"]
    args += ["--data", "data.mtz"]
    args += ["--contents", "contents.json"]
    args += ["--cycles", "1"]
    with pytest.raises(SystemExit):
        main(args)
    with open(os.path.join("modelcraft", "modelcraft.json")) as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["termination_reason"] == "Normal"


@in_temp_directory
def test_1rxf_from_model():
    hklin = ccp4_path("examples", "data", "1rxf.mtz")
    xyzin = ccp4_path("examples", "data", "1rxf_randomise.pdb")
    contents = pdb1rxf_contents()
    contents.write_sequence_file("sequence.fasta")
    args = ["xray"]
    args += ["--data", hklin]
    args += ["--observations", "I,SIGI"]
    args += ["--freerflag", "FreeR_flag"]
    args += ["--contents", "sequence.fasta"]
    args += ["--model", xyzin]
    args += ["--cycles", "2"]
    with pytest.raises(SystemExit):
        main(args)
    with open(os.path.join("modelcraft", "modelcraft.json")) as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["termination_reason"] == "Normal"
    structure = read_structure(os.path.join("modelcraft", "modelcraft.cif"))
    assert contains_residue(structure, "FE")
