try:
    from pandas._libs.tslibs import (
        is_date_array_normalized,
        normalize_i8_timestamps,
    )
except ImportError:
    from pandas._libs.tslibs.conversion import (
        normalize_i8_timestamps,
        is_date_array_normalized,
    )

import pandas as pd

from .tslib import (
    _sizes,
    _tzs,
    tzlocal_obj,
)


class Normalize:
    params = [
        _sizes,
        _tzs,
    ]
    param_names = ["size", "tz"]

    def setup(self, size, tz):
        # use an array that will have is_date_array_normalized give True,
        #  so we do not short-circuit early.
        dti = pd.date_range("2016-01-01", periods=10, tz=tz).repeat(size // 10)
        self.i8data = dti.asi8

        if size == 10 ** 6 and tz is tzlocal_obj:
            # tzlocal is cumbersomely slow, so skip to keep runtime in check
            raise NotImplementedError

    def time_normalize_i8_timestamps(self, size, tz):
        normalize_i8_timestamps(self.i8data, tz)

    def time_is_date_array_normalized(self, size, tz):
        # TODO: cases with different levels of short-circuiting
        is_date_array_normalized(self.i8data, tz)
