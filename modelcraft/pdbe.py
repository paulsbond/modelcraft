import functools
import math
import os
import re
import time
import requests
from .contents import AsuContents, Carb, Ligand, Polymer, PolymerType


def pdbe_molecules(pdbid: str) -> list:
    pdbid = pdbid.lower()
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/" + pdbid
    time.sleep(2)
    print("Requesting PDBe molecule data for", pdbid)
    response = requests.get(url)
    if response.status_code != 200:
        raise ConnectionError(response.text)
    mols = response.json()[pdbid]
    if any(mol["molecule_type"] == "carbohydrate polymer" for mol in mols):
        codes = _carb_codes(pdbid)
        for mol in mols:
            mol["carb_codes"] = codes.get(mol["entity_id"])
    return mols


def _carb_codes(pdbid: str) -> dict:
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    query = "pdb_id:" + pdbid
    filter_list = "carb_compound_id_entity"
    request_data = {"q": query, "fl": filter_list, "wt": "json"}
    time.sleep(2)
    response = requests.post(url, data=request_data)
    if response.status_code != 200:
        raise ConnectionError(response.text)
    docs = response.json()["response"]["docs"]
    codes = {}
    for doc in docs:
        for line in doc["carb_compound_id_entity"]:
            match = re.match(r"(.+)\((\d+)\)_(\d+)", line)
            code, copies, entity = match.groups()
            codes.setdefault(int(entity), {})[code] = int(copies)
    return codes


def carb_from_pdbe_molecule_dict(mol: dict) -> Carb:
    codes = mol["carb_codes"]
    length = sum(codes.values())
    copies = mol["number_of_copies"] // length
    return Carb(codes=codes, copies=copies)


def contents_from_pdbid(self, pdbid: str) -> AsuContents:
    contents = AsuContents()
    contents.copies = 1
    for mol in pdbe_molecules(pdbid):
        if "sequence" in mol:
            polymer = polymer_from_pdbe_molecule_dict(mol)
            contents.add_polymer(polymer)
        if mol["molecule_type"] == "carbohydrate polymer":
            carb = carb_from_pdbe_molecule_dict(mol)
            contents.carbs.append(carb)
        if mol["molecule_type"] == "bound":
            ligand = ligand_from_pdbe_molecule_dict(mol)
            if is_buffer(ligand.code):
                self.buffers.append(ligand.code)
            elif ligand.code not in ("UNL", "UNX"):
                self.ligands.append(ligand)
    _divide_copies(contents)


def polymer_from_pdbe_molecule_dict(mol: dict) -> Polymer:
    polymer_type = {
        "polypeptide(l)": PolymerType.PROTEIN,
        "polyribonucleotide": PolymerType.RNA,
        "polydeoxyribonucleotide": PolymerType.DNA,
    }[mol["molecule_type"]]
    return Polymer(
        sequence=mol["sequence"],
        copies=mol["number_of_copies"],
        polymer_type=polymer_type,
        modifications=_modifications_in_pdbe_molecule_dict(mol),
    )


def ligand_from_pdbe_molecule_dict(mol: dict) -> Ligand:
    return Ligand(code=mol["chem_comp_ids"][0], copies=mol["number_of_copies"])


def _modifications_in_pdbe_molecule_dict(mol: dict) -> list:
    indices = {}
    for index, mod in mol["pdb_sequence_indices_with_multiple_residues"].items():
        code1 = mod["one_letter_code"]
        code3 = mod["three_letter_code"]
        if code3 not in ("DA", "DC", "DG", "DT"):
            key = code1, code3
            indices.setdefault(key, []).append(index)
    modifications = []
    for key in indices:
        code1, code3 = key
        total = mol["sequence"].count(code1)
        if code1 == "M" and mol["sequence"][0] == "M":
            total -= 1
        if len(indices[key]) >= total:
            modifications.append(f"{code1}->{code3}")
        else:
            modifications.extend(f"{index}->{code3}" for index in indices[key])
    return modifications


@functools.lru_cache(maxsize=None)
def is_buffer(code: str) -> float:
    return code.upper() in _buffers()


@functools.lru_cache(maxsize=None)
def _buffers() -> set:
    path = os.path.join(os.environ["CCP4"], "share", "pisa", "agents.dat")
    agents = set()
    with open(path) as stream:
        for line in stream:
            if line[0] != "#" and "," in line:
                code = line.split(",")[0]
                agents.add(code)
    return agents


def _divide_copies(contents: AsuContents):
    copies = []
    for item in (
        contents.proteins
        + contents.rnas
        + contents.dnas
        + contents.carbs
        + contents.ligands
    ):
        if item.copies is not None:
            copies.append(item.copies)
    divisor = copies[0] if len(copies) == 1 else functools.reduce(math.gcd, copies)
    contents.copies *= divisor
    for item in (
        contents.proteins
        + contents.rnas
        + contents.dnas
        + contents.carbs
        + contents.ligands
    ):
        item.copies //= divisor
