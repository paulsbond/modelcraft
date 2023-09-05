import json
import os
import pytest
from modelcraft.scripts.modelcraft import main
from ..ccp4 import in_temp_directory
from . import (
    density_path,
    halfmap1_path,
    halfmap2_path,
    mask_path,
    sequence_path,
    structure_path,
)


@in_temp_directory
def test_3488_model():
    args = ["em"]
    args += ["--contents", sequence_path()]
    args += ["--map", density_path()]
    args += ["--resolution", "3.2"]
    args += ["--blur", "10"]
    args += ["--cycles", "1"]
    args += ["--model", structure_path()]
    with pytest.raises(SystemExit):
        main(args)
    with open(os.path.join("modelcraft", "modelcraft.json")) as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["final"]["fsc"] > 0.6
    assert report["termination_reason"] == "Normal"


@in_temp_directory
def test_3488_halfmaps():
    args = ["em"]
    args += ["--contents", sequence_path()]
    args += ["--map", halfmap1_path(), halfmap2_path()]
    args += ["--mask", mask_path()]
    args += ["--resolution", "3.2"]
    args += ["--cycles", "1"]
    with pytest.raises(SystemExit):
        main(args)
    with open(os.path.join("modelcraft", "modelcraft.json")) as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["termination_reason"] == "Normal"
