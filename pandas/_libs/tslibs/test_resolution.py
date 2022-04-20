from pandas._libs.tslibs import Resolution, get_resolution
import numpy as np


def test_get_resolution_nano():
    # don't return the fallback RESO_DAY
    arr = np.array([1])
    res = get_resolution(arr)
    assert res == Resolution.RESO_NS
