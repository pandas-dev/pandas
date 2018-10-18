import numpy as np

from pandas._libs import index as li


def _get_numeric_engines():
    engine_names = [
        ('Int64Engine', np.int64), ('Int32Engine', np.int32),
        ('Int16Engine', np.int16), ('Int8Engine', np.int8),
        ('UInt64Engine', np.uint64), ('UInt32Engine', np.uint32),
        ('UInt16engine', np.uint16), ('UInt8Engine', np.uint8),
        ('Float64Engine', np.float64), ('Float32Engine', np.float32),
    ]
    return [(getattr(li, engine_name), dtype)
            for engine_name, dtype in engine_names if hasattr(li, engine_name)]


class NumericEngineIndexing(object):

    params = [_get_numeric_engines(),
              ['monotonic_incr', 'monotonic_decr', 'non_monotonic'],
              ]
    param_names = ['engine_and_dtype', 'index_type']

    def setup(self, engine_and_dtype, index_type):
        engine, dtype = engine_and_dtype
        N = 10**5
        values = list([1] * N + [2] * N + [3] * N)
        arr = {
            'monotonic_incr': np.array(values, dtype=dtype),
            'monotonic_decr': np.array(list(reversed(values)),
                                       dtype=dtype),
            'non_monotonic': np.array([1, 2, 3] * N, dtype=dtype),
        }[index_type]

        self.data = engine(lambda: arr, len(arr))
        # code belows avoids populating the mapping etc. while timing.
        self.data.get_loc(2)

    def time_get_loc(self, engine_and_dtype, index_type):
        self.data.get_loc(2)


class ObjectEngineIndexing(object):

    params = [('monotonic_incr', 'monotonic_decr', 'non_monotonic')]
    param_names = ['index_type']

    def setup(self, index_type):
        N = 10**5
        values = list('a' * N + 'b' * N + 'c' * N)
        arr = {
            'monotonic_incr': np.array(values, dtype=object),
            'monotonic_decr': np.array(list(reversed(values)), dtype=object),
            'non_monotonic': np.array(list('abc') * N, dtype=object),
        }[index_type]

        self.data = li.ObjectEngine(lambda: arr, len(arr))
        # code belows avoids populating the mapping etc. while timing.
        self.data.get_loc('b')

    def time_get_loc(self, index_type):
        self.data.get_loc('b')
