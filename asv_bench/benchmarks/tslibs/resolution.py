"""
ipython analogue:

tr = TimeResolution()
mi = pd.MultiIndex.from_product(tr.params[:-1] + ([str(x) for x in tr.params[-1]],))
df = pd.DataFrame(np.nan, index=mi, columns=["mean", "stdev"])

for unit in tr.params[0]:
    for size in tr.params[1]:
        for tz in tr.params[2]:
            tr.setup(unit, size, tz)
            key = (unit, size, str(tz))
            print(key)

            val = %timeit -o tr.time_get_resolution(unit, size, tz)

            df.loc[key] = (val.average, val.stdev)

"""
import numpy as np

try:
    from pandas._libs.tslibs import get_resolution
except ImportError:
    from pandas._libs.tslibs.resolution import get_resolution

from .tslib import (
    _sizes,
    _tzs,
    tzlocal_obj,
)


class TimeResolution:
    params = (
        ["D", "h", "m", "s", "us", "ns"],
        _sizes,
        _tzs,
    )
    param_names = ["unit", "size", "tz"]

    def setup(self, unit, size, tz):
        if size == 10 ** 6 and tz is tzlocal_obj:
            # tzlocal is cumbersomely slow, so skip to keep runtime in check
            raise NotImplementedError

        arr = np.random.randint(0, 10, size=size, dtype="i8")
        arr = arr.view(f"M8[{unit}]").astype("M8[ns]").view("i8")
        self.i8data = arr

    def time_get_resolution(self, unit, size, tz):
        get_resolution(self.i8data, tz)
