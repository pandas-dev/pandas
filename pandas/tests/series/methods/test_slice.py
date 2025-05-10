import gc
import pandas as pd
from pandas import Series


class Dummy:
    def __init__(self, val):
        self.val = val


def count_dummy_instances():
    gc.collect()
    return sum(1 for obj in gc.get_objects() if isinstance(obj, Dummy))


def test_slicing_releases_dummy_instances():
    """Ensure that slicing a Series does not retain references to all original Dummy instances."""
    NDATA = 100_000
    baseline = count_dummy_instances()
    a = Series([Dummy(i) for i in range(NDATA)])
    a = a[-1:]
    gc.collect()
    after = count_dummy_instances()
    retained = after - baseline
    assert (
        retained <= 1
    ), f"{retained} Dummy instances were retained; expected at most 1"
