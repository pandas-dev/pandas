"""
ipython analogue:

tr = TimeIntsToPydatetime()
mi = pd.MultiIndex.from_product(
    tr.params[:-1] + ([str(x) for x in tr.params[-1]],)
)
df = pd.DataFrame(np.nan, index=mi, columns=["mean", "stdev"])
for box in tr.params[0]:
    for size in tr.params[1]:
        for tz in tr.params[2]:
            tr.setup(box, size, tz)
            key = (box, size, str(tz))
            print(key)
            val = %timeit -o tr.time_ints_to_pydatetime(box, size, tz)
            df.loc[key] = (val.average, val.stdev)
"""
from datetime import (
    timedelta,
    timezone,
)

from dateutil.tz import (
    gettz,
    tzlocal,
)
import numpy as np
import pytz

try:
    from pandas._libs.tslibs import ints_to_pydatetime
except ImportError:
    from pandas._libs.tslib import ints_to_pydatetime

tzlocal_obj = tzlocal()
_tzs = [
    None,
    timezone.utc,
    timezone(timedelta(minutes=60)),
    pytz.timezone("US/Pacific"),
    gettz("Asia/Tokyo"),
    tzlocal_obj,
]
_sizes = [0, 1, 100, 10 ** 4, 10 ** 6]


class TimeIntsToPydatetime:
    params = (
        ["time", "date", "datetime", "timestamp"],
        _sizes,
        _tzs,
    )
    param_names = ["box", "size", "tz"]
    # TODO: fold? freq?

    def setup(self, box, size, tz):
        if box == "date" and tz is not None:
            # tz is ignored, so avoid running redundant benchmarks
            raise NotImplementedError  # skip benchmark
        if size == 10 ** 6 and tz is _tzs[-1]:
            # This is cumbersomely-slow, so skip to trim runtime
            raise NotImplementedError  # skip benchmark

        arr = np.random.randint(0, 10, size=size, dtype="i8")
        self.i8data = arr

    def time_ints_to_pydatetime(self, box, size, tz):
        ints_to_pydatetime(self.i8data, tz, box=box)
