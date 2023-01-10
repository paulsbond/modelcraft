import json
import os
import pytest
from modelcraft.scripts.modelcraft import main
from ..integration import in_temp_directory
from . import density_path, sequence_path


@in_temp_directory
def test_3488():
    args = ["em"]
    args += ["--contents", sequence_path()]
    args += ["--map", density_path()]
    args += ["--resolution", "3.2"]
    args += ["--cycles", "1"]
    with pytest.raises(SystemExit):
        main(args)
    with open(os.path.join("modelcraft", "modelcraft.json")) as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["termination_reason"] == "Normal"
