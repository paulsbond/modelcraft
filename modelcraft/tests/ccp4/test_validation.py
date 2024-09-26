import math
import numpy as np
from ...validation import validate_refmac
from . import insulin_refmac


def test_insulin():
    refmac = insulin_refmac()
    metrics = validate_refmac(refmac)
    assert len(metrics) > 0
    assert math.isclose(np.median(metrics["BFac"]), 0)
    assert math.isclose(np.median(metrics["RSCC"]), 0)
    assert math.isclose(np.median(metrics["Diff"]), 0)
    assert math.isclose(np.median(metrics["Geom"]), 0)
    assert math.isclose(np.median(metrics["Score"]), 0)
