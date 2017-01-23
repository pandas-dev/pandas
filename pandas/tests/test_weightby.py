import numpy as np
import pandas as pd

from pandas import DataFrame, Series
from pandas.util import testing as tm
from pandas.core import common as com


class TestWeights(tm.TestCase):

    def setUp(self):
        self.df = DataFrame({'A': [0.25, 0.25, 0.25, 0.25],
                             'B': [1, 2, 3, 4]})
        self.df2 = DataFrame({'A': [1, 2, 3, 4],
                              'B': [1, 2, 3, 4]})

    def test_basic(self):

        for f in ['sum', 'mean']:
            weights = self.df[['A']] / self.df.A.sum()
            result = getattr(self.df.weightby('A'), f)()
            expected = getattr(self.df[['B']] * weights.values, f)()
            tm.assert_series_equal(result, expected)

            weights2 = self.df2[['A']] / self.df2.A.sum()
            result = getattr(self.df2.weightby('A'), f)()
            expected = getattr(self.df2[['B']] * weights2.values, f)()
            tm.assert_series_equal(result, expected)

        for f in ['kurt', 'skew', 'sem']:
            weights = self.df[['A']] / self.df.A.sum()
            result = getattr(self.df.weightby('A'), f)()
            expected = getattr(self.df[['B']] * weights.values, f)()
            # tm.assert_series_equal(result, expected)

            weights2 = self.df2[['A']] / self.df2.A.sum()
            result = getattr(self.df2.weightby('A'), f)()
            expected = getattr(self.df2[['B']] * weights2.values, f)()
            # tm.assert_series_equal(result, expected)

        for f in ['std', 'var']:

            weights = self.df[['A']] / self.df.A.sum()
            result = getattr(self.df.weightby('A'), f)(ddof=2)
            expected = getattr(self.df[['B']] * weights.values, f)(ddof=2)
            # tm.assert_series_equal(result, expected)

            weights2 = self.df2[['A']] / self.df2.A.sum()
            result = getattr(self.df2.weightby('A'), f)(ddof=2)
            expected = getattr(self.df2[['B']] * weights2.values, f)(ddof=2)
            # tm.assert_series_equal(result, expected)

    def test_gotitem(self):

        result = self.df.weightby('A')['B'].sum()
        expected = self.df.weightby('A').sum()['B']
        self.assertEqual(result, expected)

        result = self.df.weightby('A').B.sum()
        self.assertEqual(result, expected)

        result = self.df['B'].weightby(self.df['A']).sum()
        self.assertEqual(result, expected)

    def test_sample_deprecation(self):
        rs = com._random_state(1234)
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = self.df.sample(2, random_state=rs, weights='A')

        expected = self.df.iloc[[0, 2]][['B']]
        tm.assert_frame_equal(result, expected)

    def test_unsupported(self):
        for f in ['first', 'median', 'min', 'max', 'prod']:

            def func():
                getattr(self.df.weightby('A'), f)()
            self.assertRaises(AttributeError, func)

    def test_panel_unsupported(self):
        panel = pd.Panel(items=[0, 1, 2], major_axis=[2, 3, 4],
                         minor_axis=[3, 4, 5])
        with tm.assertRaises(AttributeError):
            panel.weightby('weight_column')

    def test_weights_validation(self):
        o = DataFrame(np.random.randn(10, 10))

        # Weight length must be right
        with tm.assertRaises(ValueError):
            o.weightby([0, 1]).sample(n=3)

        with tm.assertRaises(ValueError):
            bad_weights = [0.5] * 11
            o.weightby(bad_weights).sample(n=3)

        with tm.assertRaises(ValueError):
            bad_weight_series = Series([0, 0, 0.2])
            o.weightby(bad_weight_series).sample(n=4)

        # Check won't accept negative weights
        with tm.assertRaises(ValueError):
            bad_weights = [-0.1] * 10
            o.weightby(bad_weights).sample(n=3)

        # Check inf and -inf throw errors:
        with tm.assertRaises(ValueError):
            weights_with_inf = [0.1] * 10
            weights_with_inf[0] = np.inf
            o.weightby(weights_with_inf).sample(n=3)

        with tm.assertRaises(ValueError):
            weights_with_ninf = [0.1] * 10
            weights_with_ninf[0] = -np.inf
            o.weightby(weights_with_ninf).sample(n=3)

        # All zeros raises errors
        zero_weights = [0] * 10
        with tm.assertRaises(ValueError):
            o.weightby(zero_weights).sample(n=3)

        # All missing weights
        nan_weights = [np.nan] * 10
        with tm.assertRaises(ValueError):
            o.weightby(nan_weights).sample(n=3)

        # Check np.nan are replaced by zeros.
        weights_with_nan = [np.nan] * 10
        weights_with_nan[5] = 0.5
        tm.assert_frame_equal(
            o.weightby(weights_with_nan, axis=0).sample(n=1), o.iloc[5:6])

        # Check None are also replaced by zeros.
        weights_with_None = [None] * 10
        weights_with_None[5] = 0.5
        tm.assert_frame_equal(
            o.weightby(weights_with_None, axis=0).sample(n=1), o.iloc[5:6])

    def test_weights_strings(sel):
        # Fixes issue: 2419
        # additional specific object based tests

        # A few dataframe test with degenerate weights.
        easy_weight_list = [0] * 10
        easy_weight_list[5] = 1

        df = pd.DataFrame({'col1': range(10, 20),
                           'col2': range(20, 30),
                           'colString': ['a'] * 10,
                           'easyweights': easy_weight_list})
        result = df.weightby('easyweights').sample(n=1)
        expected = df.iloc[5:6, 0:-1]
        tm.assert_frame_equal(result, expected)

        # Ensure proper error if string given as weight for Series, panel, or
        # DataFrame with axis = 1.
        s = Series(range(10))
        with tm.assertRaises(ValueError):
            s.weightby('weight_column').sample(n=3)

        with tm.assertRaises(ValueError):
            df.weightby('weight_column', axis=1).sample(n=1)

        # Check weighting key error
        with tm.assertRaises(KeyError):
            df.weightby('not_a_real_column_name').sample(n=3)

        # Check that re-normalizes weights that don't sum to one.
        weights_less_than_1 = [0] * 10
        weights_less_than_1[0] = 0.5
        tm.assert_frame_equal(
            df.weightby(weights_less_than_1).sample(n=1), df.iloc[:1])

    def test_weights_axis(sel):

        # Test axis argument
        df = pd.DataFrame({'col1': range(10), 'col2': ['a'] * 10})
        second_column_weight = [0, 1]
        result = df.weightby(second_column_weight, axis=1).sample(n=1)
        tm.assert_frame_equal(result, df[['col2']])

        # Different axis arg types
        result = df.weightby(second_column_weight, axis='columns').sample(n=1)
        tm.assert_frame_equal(result, df[['col2']])

        weight = [0] * 10
        weight[5] = 0.5
        tm.assert_frame_equal(df.weightby(weight, axis='index').sample(n=1),
                              df.iloc[5:6])

        # Test weight length compared to correct axis
        with tm.assertRaises(ValueError):
            df.weightby([0.5] * 10, axis=1).sample(n=1)

        # Check weights with axis = 1
        easy_weight_list = [0] * 3
        easy_weight_list[2] = 1

        df = pd.DataFrame({'col1': range(10, 20),
                           'col2': range(20, 30),
                           'colString': ['a'] * 10})
        result = df.weightby(easy_weight_list, axis=1).sample(n=1)
        tm.assert_frame_equal(result, df[['colString']])

        # Test that function aligns weights with frame
        df = DataFrame(
            {'col1': [5, 6, 7],
             'col2': ['a', 'b', 'c'], }, index=[9, 5, 3])
        s = Series([1, 0, 0], index=[3, 5, 9])
        result = df.weightby(s).sample(1)
        tm.assert_frame_equal(result, df.loc[[3]])

        # Weights have index values to be dropped because not in
        # sampled DataFrame
        s2 = Series([0.001, 0, 10000], index=[3, 5, 10])
        result = df.weightby(s2).sample(1)
        tm.assert_frame_equal(result, df.loc[[3]])

        # Weights have empty values to be filed with zeros
        s3 = Series([0.01, 0], index=[3, 5])
        result = df.weightby(s3).sample(1)
        tm.assert_frame_equal(result, df.loc[[3]])

        # No overlap in weight and sampled DataFrame indices
        s4 = Series([1, 0], index=[1, 2])
        with tm.assertRaises(ValueError):
            df.weightby(s4).sample(1)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
