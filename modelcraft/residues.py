from functools import lru_cache
from typing import Set
import os


@lru_cache(maxsize=None)
def is_buffer(name: str) -> float:
    return name in _buffer_codes()


@lru_cache(maxsize=None)
def _buffer_codes() -> Set[str]:
    path = os.path.join(os.environ["CCP4"], "share", "pisa", "agents.dat")
    agents = set()
    with open(path) as stream:
        for line in stream:
            if line[0] != "#" and "," in line:
                code = line.split(",")[0]
                agents.add(code)
    return agents
