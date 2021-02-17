"""
Benchmarks for code in pandas/_libs, excluding pandas/_libs/tslibs,
which has its own directory
"""
import numpy as np

from pandas._libs.lib import (
    is_list_like,
    is_scalar,
)

from pandas import (
    NA,
    NaT,
)

# TODO: share with something in pd._testing?
scalars = [
    0,
    1.0,
    1 + 2j,
    True,
    "foo",
    b"bar",
    None,
    np.datetime64(123, "ns"),
    np.timedelta64(123, "ns"),
    NaT,
    NA,
]
zero_dims = [np.array("123")]
listlikes = [np.array([1, 2, 3]), {0: 1}, {1, 2, 3}, [1, 2, 3], (1, 2, 3)]


class ScalarListLike:
    params = scalars + zero_dims + listlikes

    def time_is_list_like(self, param):
        is_list_like(param)

    def time_is_scalar(self, param):
        is_scalar(param)
