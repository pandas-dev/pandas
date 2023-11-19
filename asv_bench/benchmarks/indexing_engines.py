"""
Benchmarks in this file depend mostly on code in _libs/

We have to created masked arrays to test the masked engine though. The
array is unpacked on the Cython level.

If a PR does not edit anything in _libs, it is very unlikely that benchmarks
in this file will be affected.
"""

import numpy as np

from pandas._libs import index as libindex

from pandas.core.arrays import BaseMaskedArray


def _get_numeric_engines():
    engine_names = [
        ("Int64Engine", np.int64),
        ("Int32Engine", np.int32),
        ("Int16Engine", np.int16),
        ("Int8Engine", np.int8),
        ("UInt64Engine", np.uint64),
        ("UInt32Engine", np.uint32),
        ("UInt16engine", np.uint16),
        ("UInt8Engine", np.uint8),
        ("Float64Engine", np.float64),
        ("Float32Engine", np.float32),
    ]
    return [
        (getattr(libindex, engine_name), dtype)
        for engine_name, dtype in engine_names
        if hasattr(libindex, engine_name)
    ]


def _get_masked_engines():
    engine_names = [
        ("MaskedInt64Engine", "Int64"),
        ("MaskedInt32Engine", "Int32"),
        ("MaskedInt16Engine", "Int16"),
        ("MaskedInt8Engine", "Int8"),
        ("MaskedUInt64Engine", "UInt64"),
        ("MaskedUInt32Engine", "UInt32"),
        ("MaskedUInt16engine", "UInt16"),
        ("MaskedUInt8Engine", "UInt8"),
        ("MaskedFloat64Engine", "Float64"),
        ("MaskedFloat32Engine", "Float32"),
    ]
    return [
        (getattr(libindex, engine_name), dtype)
        for engine_name, dtype in engine_names
        if hasattr(libindex, engine_name)
    ]


class NumericEngineIndexing:
    params = [
        _get_numeric_engines(),
        ["monotonic_incr", "monotonic_decr", "non_monotonic"],
        [True, False],
        [10**5, 2 * 10**6],  # 2e6 is above SIZE_CUTOFF
    ]
    param_names = ["engine_and_dtype", "index_type", "unique", "N"]

    def setup(self, engine_and_dtype, index_type, unique, N):
        engine, dtype = engine_and_dtype

        if index_type == "monotonic_incr":
            if unique:
                arr = np.arange(N * 3, dtype=dtype)
            else:
                arr = np.array([1, 2, 3], dtype=dtype).repeat(N)
        elif index_type == "monotonic_decr":
            if unique:
                arr = np.arange(N * 3, dtype=dtype)[::-1]
            else:
                arr = np.array([3, 2, 1], dtype=dtype).repeat(N)
        else:
            assert index_type == "non_monotonic"
            if unique:
                arr = np.empty(N * 3, dtype=dtype)
                arr[:N] = np.arange(N * 2, N * 3, dtype=dtype)
                arr[N:] = np.arange(N * 2, dtype=dtype)
            else:
                arr = np.array([1, 2, 3], dtype=dtype).repeat(N)

        self.data = engine(arr)
        # code belows avoids populating the mapping etc. while timing.
        self.data.get_loc(2)

        self.key_middle = arr[len(arr) // 2]
        self.key_early = arr[2]

    def time_get_loc(self, engine_and_dtype, index_type, unique, N):
        self.data.get_loc(self.key_early)

    def time_get_loc_near_middle(self, engine_and_dtype, index_type, unique, N):
        # searchsorted performance may be different near the middle of a range
        #  vs near an endpoint
        self.data.get_loc(self.key_middle)


class MaskedNumericEngineIndexing:
    params = [
        _get_masked_engines(),
        ["monotonic_incr", "monotonic_decr", "non_monotonic"],
        [True, False],
        [10**5, 2 * 10**6],  # 2e6 is above SIZE_CUTOFF
    ]
    param_names = ["engine_and_dtype", "index_type", "unique", "N"]

    def setup(self, engine_and_dtype, index_type, unique, N):
        engine, dtype = engine_and_dtype
        dtype = dtype.lower()

        if index_type == "monotonic_incr":
            if unique:
                arr = np.arange(N * 3, dtype=dtype)
            else:
                arr = np.array([1, 2, 3], dtype=dtype).repeat(N)
            mask = np.zeros(N * 3, dtype=np.bool_)
        elif index_type == "monotonic_decr":
            if unique:
                arr = np.arange(N * 3, dtype=dtype)[::-1]
            else:
                arr = np.array([3, 2, 1], dtype=dtype).repeat(N)
            mask = np.zeros(N * 3, dtype=np.bool_)
        else:
            assert index_type == "non_monotonic"
            if unique:
                arr = np.zeros(N * 3, dtype=dtype)
                arr[:N] = np.arange(N * 2, N * 3, dtype=dtype)
                arr[N:] = np.arange(N * 2, dtype=dtype)

            else:
                arr = np.array([1, 2, 3], dtype=dtype).repeat(N)
            mask = np.zeros(N * 3, dtype=np.bool_)
            mask[-1] = True

        self.data = engine(BaseMaskedArray(arr, mask))
        # code belows avoids populating the mapping etc. while timing.
        self.data.get_loc(2)

        self.key_middle = arr[len(arr) // 2]
        self.key_early = arr[2]

    def time_get_loc(self, engine_and_dtype, index_type, unique, N):
        self.data.get_loc(self.key_early)

    def time_get_loc_near_middle(self, engine_and_dtype, index_type, unique, N):
        # searchsorted performance may be different near the middle of a range
        #  vs near an endpoint
        self.data.get_loc(self.key_middle)


class ObjectEngineIndexing:
    params = [("monotonic_incr", "monotonic_decr", "non_monotonic")]
    param_names = ["index_type"]

    def setup(self, index_type):
        N = 10**5
        values = list("a" * N + "b" * N + "c" * N)
        arr = {
            "monotonic_incr": np.array(values, dtype=object),
            "monotonic_decr": np.array(list(reversed(values)), dtype=object),
            "non_monotonic": np.array(list("abc") * N, dtype=object),
        }[index_type]

        self.data = libindex.ObjectEngine(arr)
        # code belows avoids populating the mapping etc. while timing.
        self.data.get_loc("b")

    def time_get_loc(self, index_type):
        self.data.get_loc("b")
