from typing import List
import subprocess
from datetime import datetime
import os
from dataclasses import dataclass


def get_version(modelcraft_exec: str) -> str:
    p = subprocess.run([modelcraft_exec, "--version"], stdout=subprocess.PIPE, encoding="UTF-8")
    return p.stdout.rstrip("\n")


def create_test_directory(modelcraft_exec: str, base_dir: str, pdb: str):
    version = get_version(modelcraft_exec)
    pdb = pdb.split("/")[-1]
    now = datetime.now()
    formatted_date = now.strftime("%Y-%m-%d-%H:%M")
    test_dir = os.path.join(base_dir, version, formatted_date, pdb)
    os.makedirs(test_dir, exist_ok=True)
    return test_dir


def get_test_data(data_dir: str):
    return [os.path.abspath(x.path) for x in os.scandir(data_dir)]


def run_test_worker(modelcraft_exec: str, test_folder: str, work_dir: str):
    contents_path = os.path.join(test_folder, "contents.json")
    mtz_path = os.path.join(test_folder, "cpm_hklout.mtz")
    log_path = os.path.join(work_dir, "main.log")

    with open(log_path, "w") as log_file:
        subprocess.run([
            modelcraft_exec,
            "xray",
            "--contents", contents_path,
            "--data", mtz_path,
            "--phases", "phasematch.F_phi.phi,FOM",
            "--cycles", "10",
            "--overwrite-directory"],
            cwd=work_dir,
            stdout=log_file)


def main():
    modelcraft_sources: List[str] = [
        "/Users/dialpuri/Development/modelcraft/pyenv/bin/modelcraft"
    ]

    for modelcraft_exec in modelcraft_sources:
        for pdb in get_test_data("data"):
            work_dir = create_test_directory(modelcraft_exec=modelcraft_exec, base_dir="runs", pdb=pdb)
            run_test_worker(modelcraft_exec, pdb, work_dir)
            exit()


if __name__ == "__main__":
    main()
    # generate_contents("data")


def generate_contents(test_dir: str):
    for test_data in os.scandir(test_dir):
        subprocess.run(["/Users/dialpuri/Development/modelcraft/pyenv/bin/modelcraft-contents", test_data.name,
                        os.path.join(test_data.path, "contents.json")])
