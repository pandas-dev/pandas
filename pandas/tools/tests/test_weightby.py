import numpy as np
import pandas as pd

from pandas import DataFrame, Series
from pandas.util import testing as tm
from pandas.core import common as com


class TestWeightsby(tm.TestCase):

    def setUp(self):
        self.df = DataFrame({'A': [0.25, 0.25, 0.25, 0.25],
                             'B': [1, 2, 3, 4]})
        self.df2 = DataFrame({'A': [1, 2, 3, 4],
                              'B': [1, 2, 3, 4]})
        self.df3 = DataFrame({'A': [1, 2, 3, 4],
                              'B': [1, 2, 3, 4],
                              'C': [1, 1, 2, 2]})

    @property
    def rs(self):
        # always return the same starting random state object
        return com._random_state(1234)

    def test_basic(self):

        for f in ['sum', 'mean']:
            weights = (self.df[['A']] / self.df.A.sum()).values
            result = getattr(self.df, f)(weights='A')
            expected = getattr(self.df[['B']] * weights, f)()
            tm.assert_series_equal(result, expected)

            weights2 = (self.df2[['A']] / self.df2.A.sum()).values
            result = getattr(self.df2, f)(weights='A')
            expected = getattr(self.df2[['B']] * weights2, f)()
            tm.assert_series_equal(result, expected)

        for f in ['kurt', 'skew', 'sem']:
            weights = (self.df[['A']] / self.df.A.sum()).values
            result = getattr(self.df, f)(weights='A')
            expected = getattr(self.df[['B']] * weights, f)()
            # tm.assert_series_equal(result, expected)

            weights2 = (self.df2[['A']] / self.df2.A.sum()).values
            result = getattr(self.df2, f)(weights='A')
            expected = getattr(self.df2[['B']] * weights2, f)()
            # tm.assert_series_equal(result, expected)

        for f in ['std', 'var']:

            weights = (self.df[['A']] / self.df.A.sum()).values
            result = getattr(self.df, f)(weights='A', ddof=2)
            expected = getattr(self.df[['B']] * weights, f)(ddof=2)
            # tm.assert_series_equal(result, expected)

            weights2 = (self.df2[['A']] / self.df2.A.sum()).values
            result = getattr(self.df2, f)(weights='A', ddof=2)
            expected = getattr(self.df2[['B']] * weights2, f)(ddof=2)
            # tm.assert_series_equal(result, expected)

    def test_groupby(self):

        for f in ['mean', 'sum']:

            weights = (self.df3['A'] / self.df3.A.sum()).values
            result = getattr(self.df3.groupby('C'), f)(weights='A')
            adj = self.df3.assign(A=self.df3.A * weights,
                                  B=self.df3.B * weights)
            expected = getattr(adj.groupby('C'), f)()
            tm.assert_frame_equal(result, expected)

            weights = (self.df3['A'] / self.df3.A.sum()).values
            result = getattr(self.df3.groupby('C').B, f)(weights='A')
            adj = self.df3.assign(B=self.df3.B * weights)
            expected = getattr(adj.groupby('C').B, f)()
            tm.assert_series_equal(result, expected)

    def test_unsupported(self):
        for f in ['first', 'median', 'min', 'max', 'prod']:

            def func():
                getattr(self.df, f)(weights='A')
            self.assertRaises(TypeError, func)

    def test_panel_unsupported(self):
        panel = pd.Panel(items=[0, 1, 2], major_axis=[2, 3, 4],
                         minor_axis=[3, 4, 5])
        with tm.assertRaises(NotImplementedError):
            panel.sum(weights='weight_column')

    def test_weights_validation(self):
        o = DataFrame(np.random.randn(10, 10))

        # Weight length must be right
        with tm.assertRaises(ValueError):
            o.sample(n=3, random_state=self.rs, weights=[0, 1])

        with tm.assertRaises(ValueError):
            bad_weights = [0.5] * 11
            o.sample(n=3, random_state=self.rs, weights=bad_weights)

        # Check won't accept negative weights
        with tm.assertRaises(ValueError):
            bad_weights = [-0.1] * 10
            o.sample(n=3, random_state=self.rs, weights=bad_weights)

        # Check inf and -inf throw errors:
        with tm.assertRaises(ValueError):
            weights_with_inf = [0.1] * 10
            weights_with_inf[0] = np.inf
            o.sample(n=3, random_state=self.rs, weights=weights_with_inf)

        with tm.assertRaises(ValueError):
            weights_with_ninf = [0.1] * 10
            weights_with_ninf[0] = -np.inf
            o.sample(n=3, random_state=self.rs, weights=weights_with_ninf)

        # All zeros raises errors
        zero_weights = [0] * 10
        with tm.assertRaises(ValueError):
            o.sample(n=3, random_state=self.rs, weights=zero_weights)

        # All missing weights
        nan_weights = [np.nan] * 10
        with tm.assertRaises(ValueError):
            o.sample(n=3, random_state=self.rs, weights=nan_weights)

        # Check np.nan are replaced by zeros.
        weights_with_nan = [np.nan] * 10
        weights_with_nan[5] = 0.5
        result = o.sample(n=1, random_state=self.rs, weights=weights_with_nan)
        expected = o.iloc[5:6]
        tm.assert_frame_equal(result, expected)

        # Check None are also replaced by zeros.
        weights_with_None = [None] * 10
        weights_with_None[5] = 0.5
        result = o.sample(n=1, random_state=self.rs, weights=weights_with_None)
        expected = o.iloc[5:6]
        tm.assert_frame_equal(result, expected)

    def test_weights_strings(self):
        # Fixes issue: 2419
        # additional specific object based tests

        # A few dataframe test with degenerate weights.
        easy_weight_list = [0] * 10
        easy_weight_list[5] = 1

        df = pd.DataFrame({'col1': range(10, 20),
                           'col2': range(20, 30),
                           'colString': ['a'] * 10,
                           'easyweights': easy_weight_list})
        result = df.sample(n=1, random_state=self.rs, weights='easyweights')
        expected = df[['col1', 'col2', 'colString']].iloc[5:6]
        tm.assert_frame_equal(result, expected)

        # Ensure proper error if string given as weight for Series, panel, or
        # DataFrame with axis = 1.
        s = Series(range(10))
        with tm.assertRaises(ValueError):
            s.sample(n=3, random_state=self.rs, weights='weight_column')

        with tm.assertRaises(ValueError):
            df.sample(n=1, random_state=self.rs,
                      weights='weight_column', axis=1)

        # Check weighting key error
        with tm.assertRaises(KeyError):
            df.sample(n=3, random_state=self.rs,
                      weights='not_a_real_column_name')

        # Check that re-normalizes weights that don't sum to one.
        weights_less_than_1 = [0] * 10
        weights_less_than_1[0] = 0.5
        result = df.sample(n=1, random_state=self.rs,
                           weights=weights_less_than_1)
        expected = df.iloc[[0]]
        tm.assert_frame_equal(result, expected)

    def test_weights_axis(self):

        # Test axis argument
        df = pd.DataFrame({'col1': range(10), 'col2': ['a'] * 10})
        second_column_weight = [0, 1]
        result = df.sample(n=1, random_state=self.rs,
                           weights=second_column_weight, axis=1)
        tm.assert_frame_equal(result, df[['col2']])

        # Different axis arg types
        result = df.sample(n=1, random_state=self.rs,
                           weights=second_column_weight, axis='columns')
        tm.assert_frame_equal(result, df[['col2']])

        weight = [0] * 10
        weight[5] = 0.5
        result = df.sample(n=1, random_state=self.rs,
                           weights=weight, axis='index')
        expected = df.iloc[5:6]
        tm.assert_frame_equal(result, expected)

        # Test weight length compared to correct axis
        with tm.assertRaises(ValueError):
            df.sample(n=1, random_state=self.rs, weights=[0.5] * 10, axis=1)

        # Check weights with axis = 1
        easy_weight_list = [0] * 3
        easy_weight_list[2] = 1

        df = pd.DataFrame({'col1': range(10, 20),
                           'col2': range(20, 30),
                           'colString': ['a'] * 10})
        result = df.sample(n=1, random_state=self.rs,
                           weights=easy_weight_list, axis=1)
        expected = df[['colString']]
        tm.assert_frame_equal(result, expected)

        # Test that function aligns weights with frame
        df = DataFrame(
            {'col1': [5, 6, 7],
             'col2': ['a', 'b', 'c'], }, index=[9, 5, 3])
        s = Series([1, 0, 0], index=[3, 5, 9])
        result = df.sample(1, random_state=self.rs, weights=s)
        tm.assert_frame_equal(result, df.loc[[3]])

        # Weights have index values to be dropped because not in
        # sampled DataFrame
        s2 = Series([0.001, 0, 10000], index=[3, 5, 10])
        result = df.sample(1, random_state=self.rs, weights=s2)
        tm.assert_frame_equal(result, df.loc[[3]])

        # Weights have empty values to be filed with zeros
        s3 = Series([0.01, 0], index=[3, 5])
        result = df.sample(1, random_state=self.rs, weights=s3)
        tm.assert_frame_equal(result, df.loc[[3]])

        # No overlap in weight and sampled DataFrame indices
        s4 = Series([1, 0], index=[1, 2])
        with tm.assertRaises(ValueError):
            df.sample(1, random_state=self.rs, weights=s4)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
