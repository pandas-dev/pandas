import numpy as np

import pandas as pd


class BooleanArray:
    def setup(self):
        self.values_bool = np.array([True, False, True, False])
        self.values_float = np.array([1.0, 0.0, 1.0, 0.0])
        self.values_integer = np.array([1, 0, 1, 0])
        self.values_integer_like = [1, 0, 1, 0]
        self.data = np.array([True, False, True, False])
        self.mask = np.array([False, False, True, False])

    def time_constructor(self):
        pd.arrays.BooleanArray(self.data, self.mask)

    def time_from_bool_array(self):
        pd.array(self.values_bool, dtype="boolean")

    def time_from_integer_array(self):
        pd.array(self.values_integer, dtype="boolean")

    def time_from_integer_like(self):
        pd.array(self.values_integer_like, dtype="boolean")

    def time_from_float_array(self):
        pd.array(self.values_float, dtype="boolean")


class IntegerArray:
    def setup(self):
        self.values_integer = np.array([1, 0, 1, 0])
        self.data = np.array([1, 2, 3, 4], dtype="int64")
        self.mask = np.array([False, False, True, False])

    def time_constructor(self):
        pd.arrays.IntegerArray(self.data, self.mask)

    def time_from_integer_array(self):
        pd.array(self.values_integer, dtype="Int64")


class Isna:
    def setup(self):
        N = 10 ** 8
        data = np.random.default_rng().choice(np.ones(2 ** 20), size=N)
        self.pandas_array = pd.array(data, dtype=int)
        self.integer_array = pd.array(data, dtype="Int64")

    def time_isna_integer_array(self):
        pd.isna(self.integer_array)

    def time_isna_pandas_array(self):
        pd.isna(self.pandas_array)
