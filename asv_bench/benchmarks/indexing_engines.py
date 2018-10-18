import numpy as np

from pandas._libs.index import (Int64Engine, UInt64Engine, Float64Engine,
                                ObjectEngine)


class NumericEngineIndexing(object):

    goal_time = 0.2
    params = [[Int64Engine, UInt64Engine, Float64Engine],
              [np.int64, np.uint64, np.float64],
              ['monotonic_incr', 'monotonic_decr', 'non_monotonic'],
              ]
    param_names = ['engine', 'dtype', 'index_type']

    def setup(self, engine, dtype, index_type):
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

    def time_get_loc(self, engine, dtype, index_type):
        self.data.get_loc(2)


class ObjectEngineIndexing(object):

    goal_time = 0.2
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

        self.data = ObjectEngine(lambda: arr, len(arr))
        # code belows avoids populating the mapping etc. while timing.
        self.data.get_loc('b')

    def time_get_loc(self, index_type):
        self.data.get_loc('b')
