import numpy as np

from pandas import DatetimeIndex, Index, MultiIndex, Series, Timestamp
import pandas._testing as tm


def no_change(arr):
    return arr


def list_of_str(arr):
    return list(arr.astype(str))


def gen_of_str(arr):
    return (x for x in arr.astype(str))


def arr_dict(arr):
    return dict(zip(range(len(arr)), arr))


def list_of_tuples(arr):
    return [(i, -i) for i in arr]


def gen_of_tuples(arr):
    return ((i, -i) for i in arr)


def list_of_lists(arr):
    return [[i, -i] for i in arr]


def list_of_tuples_with_none(arr):
    return [(i, -i) for i in arr][:-1] + [None]


def list_of_lists_with_none(arr):
    return [[i, -i] for i in arr][:-1] + [None]


class SeriesConstructors:

    param_names = ["data_fmt", "with_index", "dtype"]
    params = [
        [
            no_change,
            list,
            list_of_str,
            gen_of_str,
            arr_dict,
            list_of_tuples,
            gen_of_tuples,
            list_of_lists,
            list_of_tuples_with_none,
            list_of_lists_with_none,
        ],
        [False, True],
        ["float", "int"],
    ]

    # Generators get exhausted on use, so run setup before every call
    number = 1
    repeat = (3, 250, 10)

    def setup(self, data_fmt, with_index, dtype):
        if data_fmt in (gen_of_str, gen_of_tuples) and with_index:
            raise NotImplementedError(
                "Series constructors do not support using generators with indexes"
            )
        N = 10 ** 4
        if dtype == "float":
            arr = np.random.randn(N)
        else:
            arr = np.arange(N)
        self.data = data_fmt(arr)
        self.index = np.arange(N) if with_index else None

    def time_series_constructor(self, data_fmt, with_index, dtype):
        Series(self.data, index=self.index)


class SeriesDtypesConstructors:
    def setup(self):
        N = 10 ** 4
        self.arr = np.random.randn(N)
        self.arr_str = np.array(["foo", "bar", "baz"], dtype=object)
        self.s = Series(
            [Timestamp("20110101"), Timestamp("20120101"), Timestamp("20130101")]
            * N
            * 10
        )

    def time_index_from_array_string(self):
        Index(self.arr_str)

    def time_index_from_array_floats(self):
        Index(self.arr)

    def time_dtindex_from_series(self):
        DatetimeIndex(self.s)

    def time_dtindex_from_index_with_series(self):
        Index(self.s)


class MultiIndexConstructor:
    def setup(self):
        N = 10 ** 4
        self.iterables = [tm.makeStringIndex(N), range(20)]

    def time_multiindex_from_iterables(self):
        MultiIndex.from_product(self.iterables)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
