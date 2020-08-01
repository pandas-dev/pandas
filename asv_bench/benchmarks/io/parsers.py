import numpy as np

try:
    from pandas._libs.tslibs.parsing import (
        _does_string_look_like_datetime,
        concat_date_cols,
    )
except ImportError:
    # Avoid whole benchmark suite import failure on asv (currently 0.4)
    pass


class DoesStringLookLikeDatetime:

    params = (["2Q2005", "0.0", "10000"],)
    param_names = ["value"]

    def setup(self, value):
        self.objects = [value] * 1000000

    def time_check_datetimes(self, value):
        for obj in self.objects:
            _does_string_look_like_datetime(obj)


class ConcatDateCols:

    params = ([1234567890, "AAAA"], [1, 2])
    param_names = ["value", "dim"]

    def setup(self, value, dim):
        count_elem = 10000
        if dim == 1:
            self.object = (np.array([value] * count_elem),)
        if dim == 2:
            self.object = (
                np.array([value] * count_elem),
                np.array([value] * count_elem),
            )

    def time_check_concat(self, value, dim):
        concat_date_cols(self.object)
