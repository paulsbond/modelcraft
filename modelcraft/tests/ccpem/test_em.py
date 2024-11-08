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
def test_3488_single_map():
    args = ["em"]
    args += ["--contents", sequence_path()]
    args += ["--single-map", density_path()]
    args += ["--resolution", "3.2"]
    args += ["--cycles", "1"]
    args += ["--model", structure_path()]
    args += ["--mask", "auto"]
    with pytest.raises(SystemExit):
        main(args)
    report_path = os.path.join("modelcraft", "modelcraft.json")
    with open(report_path, encoding="utf-8") as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["final"]["fsc"] > 0.6
    assert report["termination_reason"] == "Normal"


@in_temp_directory
def test_3488_half_maps():
    args = ["em"]
    args += ["--contents", sequence_path()]
    args += ["--half-maps", halfmap1_path(), halfmap2_path()]
    args += ["--build-map", density_path()]
    args += ["--mask", mask_path()]
    args += ["--resolution", "3.2"]
    args += ["--cycles", "1"]
    with pytest.raises(SystemExit):
        main(args)
    report_path = os.path.join("modelcraft", "modelcraft.json")
    with open(report_path, encoding="utf-8") as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["termination_reason"] == "Normal"
