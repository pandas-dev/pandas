import numpy as np

from pandas._libs.index import (Int64Engine, Int32Engine,
                                Int16Engine, Int8Engine,
                                UInt64Engine, UInt32Engine,
                                UInt16Engine, UInt8Engine,
                                Float64Engine, Float32Engine,
                                ObjectEngine,
                                )


class NumericEngineIndexing(object):

    goal_time = 0.2
    params = [[Int64Engine, Int32Engine, Int16Engine, Int8Engine,
               UInt64Engine, UInt32Engine, UInt16Engine, UInt8Engine,
               Float64Engine, Float32Engine,
               ],
              ['monotonic_incr', 'monotonic_decr', 'non_monotonic']
              ]
    param_names = ['engine', 'index_type']

    def setup(self, engine, index_type):
        N = 10**5
        values = list([1] * N + [2] * N + [3] * N)
        array_ = {
            'monotonic_incr': np.array(values, dtype=engine._dtype),
            'monotonic_decr': np.array(list(reversed(values)),
                                       dtype=engine._dtype),
            'non_monotonic': np.array([1, 2, 3] * N, dtype=engine._dtype),
        }[index_type]

        self.data = engine(lambda: array_, len(array_))

    def time_get_loc(self, engine, index_type):
        self.data.get_loc(2)


class ObjectEngineIndexing(object):

    goal_time = 0.2
    params = [[ObjectEngine],
              ['monotonic_incr', 'monotonic_decr', 'non_monotonic']
              ]
    param_names = ['engine', 'index_type']

    def setup(self, engine, index_type):
        N = 10**5
        values = list('a' * N + 'b' * N + 'c' * N)
        array_ = {
            'monotonic_incr': np.array(values, dtype=engine._dtype),
            'monotonic_decr': np.array(list(reversed(values)),
                                       dtype=engine._dtype),
            'non_monotonic': np.array(list('abc') * N, dtype=engine._dtype),
        }[index_type]

        self.data = engine(lambda: array_, len(array_))

    def time_get_loc(self, engine, index_type):
        self.data.get_loc(2)
