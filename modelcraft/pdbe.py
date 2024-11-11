import multiprocessing
import re
import requests


_MULTIPROCESSING_LOCK = multiprocessing.Lock()
_SERVER = "https://www.ebi.ac.uk/pdbe"


def _response_json(url, data=None):
    with _MULTIPROCESSING_LOCK:
        print("Requesting:", url)
        if data is None:
            response = requests.get(url, timeout=30)
        else:
            response = requests.post(url, data=data, timeout=30)
        response.raise_for_status()
    return response.json()


def molecule_dicts(entry_id: str) -> list:
    entry_id = _superceeding_entry(entry_id)
    url = _SERVER + "/api/pdb/entry/molecules/" + entry_id
    response = _response_json(url)
    mols = response[entry_id]
    if any(mol["molecule_type"] == "carbohydrate polymer" for mol in mols):
        codes = _carb_codes(entry_id)
        for mol in mols:
            mol["carb_codes"] = codes.get(mol["entity_id"])
    return mols


def _superceeding_entry(entry_id: str) -> str:
    entry_id = entry_id.lower()
    url = _SERVER + "/api/pdb/entry/status/" + entry_id
    response = _response_json(url)
    superceded_by = response[entry_id][0].get("superceded_by", [])
    return entry_id if len(superceded_by) == 0 else superceded_by[-1]


def _carb_codes(entry: str) -> dict:
    url = _SERVER + "/search/pdb/select?"
    query = "pdb_id:" + entry
    filter_list = "carb_compound_id_entity"
    request_data = {"q": query, "fl": filter_list, "wt": "json"}
    response = _response_json(url, data=request_data)
    docs = response["response"]["docs"]
    codes = {}
    for doc in docs:
        for line in doc["carb_compound_id_entity"]:
            match = re.match(r"(.+)\((\d+)\)_(\d+)", line)
            code, copies, entity = match.groups()
            codes.setdefault(int(entity), {})[code] = int(copies)
    return codes
