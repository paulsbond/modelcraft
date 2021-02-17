import requests
import time


def molecules(pdbid: str) -> list:
    pdbid = pdbid.lower()
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/" + pdbid
    time.sleep(1)
    print("Requesting PDBe molecule data for", pdbid)
    response = requests.get(url)
    if response.status_code != 200:
        raise ConnectionError("Could not retrieve PDBe molecule data for " + pdbid)
    return response.json()[pdbid]
