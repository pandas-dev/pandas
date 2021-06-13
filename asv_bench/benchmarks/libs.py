"""
Benchmarks for code in pandas/_libs, excluding pandas/_libs/tslibs,
which has its own directory.

If a PR does not edit anything in _libs/, then it is unlikely that thes
benchmarks will be affected.
"""
import numpy as np

from pandas._libs.lib import (
    infer_dtype,
    is_list_like,
    is_scalar,
)

from pandas import (
    NA,
    NaT,
)

from .pandas_vb_common import (
    lib,
    tm,
)

try:
    from pandas.util import cache_readonly
except ImportError:
    from pandas.util.decorators import cache_readonly


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


class FastZip:
    def setup(self):
        N = 10000
        K = 10
        key1 = tm.makeStringIndex(N).values.repeat(K)
        key2 = tm.makeStringIndex(N).values.repeat(K)
        col_array = np.vstack([key1, key2, np.random.randn(N * K)])
        col_array2 = col_array.copy()
        col_array2[:, :10000] = np.nan
        self.col_array_list = list(col_array)

    def time_lib_fast_zip(self):
        lib.fast_zip(self.col_array_list)


class InferDtype:
    param_names = ["dtype"]
    data_dict = {
        "np-object": np.array([1] * 100000, dtype="O"),
        "py-object": [1] * 100000,
        "np-null": np.array([1] * 50000 + [np.nan] * 50000),
        "py-null": [1] * 50000 + [None] * 50000,
        "np-int": np.array([1] * 100000, dtype=int),
        "np-floating": np.array([1.0] * 100000, dtype=float),
        "empty": [],
        "bytes": [b"a"] * 100000,
    }
    params = list(data_dict.keys())

    def time_infer_dtype_skipna(self, dtype):
        infer_dtype(self.data_dict[dtype], skipna=True)

    def time_infer_dtype(self, dtype):
        infer_dtype(self.data_dict[dtype], skipna=False)


class CacheReadonly:
    def setup(self):
        class Foo:
            @cache_readonly
            def prop(self):
                return 5

        self.obj = Foo()

    def time_cache_readonly(self):
        self.obj.prop
