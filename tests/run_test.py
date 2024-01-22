from typing import List, Tuple
import subprocess
from datetime import datetime
import os
from dataclasses import dataclass
from tqdm import tqdm
from multiprocessing import Pool 


def get_version(modelcraft_exec: str) -> str:
    p = subprocess.run([modelcraft_exec, "--version"], stdout=subprocess.PIPE, encoding="UTF-8")
    return p.stdout.rstrip("\n")


def create_test_directory(modelcraft_exec: str, base_dir: str, pdb: str, date: str):
    version = get_version(modelcraft_exec)
    pdb = pdb.split("/")[-1]
    
    test_dir = os.path.join(base_dir, version, date, pdb)
    os.makedirs(test_dir, exist_ok=True)
    return test_dir


def get_test_data(data_dir: str):
    return [os.path.abspath(x.path) for x in os.scandir(data_dir)]


#def run_test_worker(modelcraft_exec: str, test_folder: str, work_dir: str):
def run_test_worker(data: Tuple[str, str, str]):
    modelcraft_exec, test_folder, work_dir = data
    contents_path = os.path.join(test_folder, "contents.json")
    mtz_path = os.path.join(test_folder, "cpm_hklout.mtz")
    model_path = os.path.join(test_folder, "symmatched_xyzout.pdb")

    log_path = os.path.join(work_dir, "main.log")

    with open(log_path, "w") as log_file:
        subprocess.run([
            modelcraft_exec,
            "xray",
            "--contents", contents_path,
            "--model", model_path,
            "--data", mtz_path,
            "--cycles", "10",
            "--overwrite-directory"],
            cwd=work_dir,
            stdout=log_file)


def main():
    modelcraft_sources: List[str] = [
    # "/home/jordan/dev/modelcraft_combine/pyenv3.2/bin/modelcraft",
	"/home/jordan/dev/modelcraft_combine/pyenv/bin/modelcraft"
    ]

    now = datetime.now()
    formatted_date = now.strftime("%Y-%m-%d-%H:%M")

    for modelcraft_exec in modelcraft_sources:

        test_data = []
        for pdb in get_test_data("data"):
            work_dir = create_test_directory(modelcraft_exec, "runs", pdb, formatted_date)
            test_data.append((modelcraft_exec, pdb, work_dir))

        with Pool(12) as pool_:
            x = list(tqdm(pool_.imap_unordered(run_test_worker, test_data),total=len(test_data)))

#        for pdb in tqdm(get_test_data("data")):
 #           work_dir = create_test_directory(modelcraft_exec=modelcraft_exec, base_dir="runs", pdb=pdb)
  #          run_test_worker(modelcraft_exec, pdb, work_dir)
            


if __name__ == "__main__":
    main()
    # generate_contents("data")


def generate_contents(test_dir: str):
    for test_data in os.scandir(test_dir):
        subprocess.run(["/Users/dialpuri/Development/modelcraft/pyenv/bin/modelcraft-contents", test_data.name,
                        os.path.join(test_data.path, "contents.json")])
