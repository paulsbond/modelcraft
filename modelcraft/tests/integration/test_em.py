import json
import os
import urllib.request
import pytest
from modelcraft.scripts.modelcraft import main
from modelcraft.scripts.contents import _entry_contents
from . import in_temp_directory


@in_temp_directory
def test_7dy0():
    url_base = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures"
    map_name = "emd_30913.map.gz"
    url = f"{url_base}/EMD-30913/map/{map_name}"
    urllib.request.urlretrieve(url, map_name)
    contents = _entry_contents("7dy0")
    contents.write_json_file("contents.json")
    args = ["em"]
    args += ["--contents", "contents.json"]
    args += ["--map", map_name]
    args += ["--resolution", "1.93"]
    args += ["--cycles", "1"]
    with pytest.raises(SystemExit):
        main(args)
    with open(os.path.join("modelcraft", "modelcraft.json")) as report_file:
        report = json.load(report_file)
    assert report["seconds"]["total"] > 0
    assert report["termination_reason"] == "Normal"
