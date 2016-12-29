import numpy as np
import pandas as pd
from pandas.util import testing as tm


class Algorithms(object):
    goal_time = 0.2

    def setup(self):
        N = 100000
        np.random.seed(1234)

        self.int_unique = pd.Int64Index(np.arange(N * 5))
        # cache is_unique
        self.int_unique.is_unique

        self.int = pd.Int64Index(np.arange(N).repeat(5))
        self.float = pd.Float64Index(np.random.randn(N).repeat(5))

        # Convenience naming.
        self.checked_add = pd.core.algorithms.checked_add_with_arr

        self.arr = np.arange(1000000)
        self.arrpos = np.arange(1000000)
        self.arrneg = np.arange(-1000000, 0)
        self.arrmixed = np.array([1, -1]).repeat(500000)
        self.strings = tm.makeStringIndex(100000)

        self.arr_nan = np.random.choice([True, False], size=1000000)
        self.arrmixed_nan = np.random.choice([True, False], size=1000000)

        # match
        self.uniques = tm.makeStringIndex(1000).values
        self.all = self.uniques.repeat(10)

    def time_factorize_string(self):
        self.strings.factorize()

    def time_factorize_int(self):
        self.int.factorize()

    def time_factorize_float(self):
        self.int.factorize()

    def time_duplicated_int_unique(self):
        self.int_unique.duplicated()

    def time_duplicated_int(self):
        self.int.duplicated()

    def time_duplicated_float(self):
        self.float.duplicated()

    def time_match_strings(self):
        pd.match(self.all, self.uniques)

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

    def time_add_overflow_first_arg_nan(self):
        self.checked_add(self.arr, self.arrmixed, arr_mask=self.arr_nan)

    def time_add_overflow_second_arg_nan(self):
        self.checked_add(self.arr, self.arrmixed, b_mask=self.arrmixed_nan)

    def time_add_overflow_both_arg_nan(self):
        self.checked_add(self.arr, self.arrmixed, arr_mask=self.arr_nan,
                         b_mask=self.arrmixed_nan)


class Hashing(object):
    goal_time = 0.2

    def setup(self):
        N = 100000

        self.df = pd.DataFrame(
            {'A': pd.Series(tm.makeStringIndex(100).take(
                np.random.randint(0, 100, size=N))),
             'B': pd.Series(tm.makeStringIndex(10000).take(
                 np.random.randint(0, 10000, size=N))),
             'D': np.random.randn(N),
             'E': np.arange(N),
             'F': pd.date_range('20110101', freq='s', periods=N),
             'G': pd.timedelta_range('1 day', freq='s', periods=N),
             })
        self.df['C'] = self.df['B'].astype('category')
        self.df.iloc[10:20] = np.nan

    def time_frame(self):
        self.df.hash()

    def time_series_int(self):
        self.df.E.hash()

    def time_series_string(self):
        self.df.B.hash()

    def time_series_categorical(self):
        self.df.C.hash()
