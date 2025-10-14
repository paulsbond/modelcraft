from random import choice
from string import ascii_letters, digits

import numpy as np


def modified_zscore(a) -> np.ndarray:
    "Z-score using median and median absolute deviation"
    median = np.median(a)
    deviations = a - median
    mad = np.median(np.abs(deviations))
    return np.zeros_like(a) if mad == 0 else 0.67449 * deviations / mad


def puid(length: int = 10) -> str:
    "Probably unique identifier"
    chars = ascii_letters + digits
    return "".join(choice(chars) for _ in range(length))
