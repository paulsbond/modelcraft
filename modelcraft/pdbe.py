import requests
import re
import time


def molecules(pdbid: str) -> list:
    pdbid = pdbid.lower()
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/" + pdbid
    time.sleep(1)
    print("Requesting PDBe molecule data for", pdbid)
    response = requests.get(url)
    if response.status_code != 200:
        raise ConnectionError(response.text)
    mols = response.json()[pdbid]
    if any(mol["molecule_type"] == "carbohydrate polymer" for mol in mols):
        codes = carb_codes(pdbid)
        for mol in mols:
            mol["carb_codes"] = codes.get(mol["entity_id"])
    return mols


def carb_codes(pdbid: str) -> dict:
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    query = "pdb_id:" + pdbid
    filter_list = "carb_compound_id_entity"
    request_data = {"q": query, "fl": filter_list, "wt": "json"}
    time.sleep(1)
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
