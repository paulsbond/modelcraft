import requests
import time


def molecules(pdbid: str) -> list:
    pdbid = pdbid.lower()
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/" + pdbid
    time.sleep(1)
    response = requests.get(url)
    if response.status_code != 200:
        raise ConnectionError("Could not reqreive entry data for " + pdbid)
    return response.json()[pdbid]
