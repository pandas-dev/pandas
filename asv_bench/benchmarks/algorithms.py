import numpy as np
import pandas as pd


class algorithm(object):
    goal_time = 0.2

    def setup(self):
        N = 100000

        self.int_unique = pd.Int64Index(np.arange(N * 5))
        # cache is_unique
        self.int_unique.is_unique

        self.int = pd.Int64Index(np.arange(N).repeat(5))
        self.float = pd.Float64Index(np.random.randn(N).repeat(5))

        # Convenience naming.
        self.checked_add = pd.core.nanops._checked_add_with_arr

        self.arr = np.arange(1000000)
        self.arrpos = np.arange(1000000)
        self.arrneg = np.arange(-1000000, 0)
        self.arrmixed = np.array([1, -1]).repeat(500000)

    def time_int_factorize(self):
        self.int.factorize()

    def time_float_factorize(self):
        self.int.factorize()

    def time_int_unique_duplicated(self):
        self.int_unique.duplicated()

    def time_int_duplicated(self):
        self.int.duplicated()

    def time_float_duplicated(self):
        self.float.duplicated()

    def time_add_overflow_pos_scalar(self):
        self.checked_add(self.arr, 1)

    def time_add_overflow_neg_scalar(self):
        self.checked_add(self.arr, -1)

    def time_add_overflow_zero_scalar(self):
        self.checked_add(self.arr, 0)

    def time_add_overflow_pos_arr(self):
        self.checked_add(self.arr, self.arrpos)

    def time_add_overflow_neg_arr(self):
        self.checked_add(self.arr, self.arrneg)

    def time_add_overflow_mixed_arr(self):
        self.checked_add(self.arr, self.arrmixed)
