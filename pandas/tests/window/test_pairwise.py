import warnings

import pytest

from pandas import DataFrame, Series
from pandas.core.sorting import safe_sort
import pandas.util.testing as tm


class TestPairwise:

    # GH 7738
    df1s = [
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[0, 1]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[1, 0]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[1, 1]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=["C", "C"]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[1.0, 0]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=[0.0, 1]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1]], columns=["C", 1]),
        DataFrame([[2.0, 4.0], [1.0, 2.0], [5.0, 2.0], [8.0, 1.0]], columns=[1, 0.0]),
        DataFrame([[2, 4.0], [1, 2.0], [5, 2.0], [8, 1.0]], columns=[0, 1.0]),
        DataFrame([[2, 4], [1, 2], [5, 2], [8, 1.0]], columns=[1.0, "X"]),
    ]
    df2 = DataFrame(
        [[None, 1, 1], [None, 1, 2], [None, 3, 2], [None, 8, 1]],
        columns=["Y", "Z", "X"],
    )
    s = Series([1, 1, 3, 8])

    def compare(self, result, expected):

        # since we have sorted the results
        # we can only compare non-nans
        result = result.dropna().values
        expected = expected.dropna().values

        tm.assert_numpy_array_equal(result, expected, check_dtype=False)

    @pytest.mark.parametrize("f", [lambda x: x.cov(), lambda x: x.corr()])
    def test_no_flex(self, f):

        # DataFrame methods (which do not call _flex_binary_moment())

        results = [f(df) for df in self.df1s]
        for (df, result) in zip(self.df1s, results):
            tm.assert_index_equal(result.index, df.columns)
            tm.assert_index_equal(result.columns, df.columns)
        for i, result in enumerate(results):
            if i > 0:
                self.compare(result, results[0])

    @pytest.mark.parametrize(
        "f",
        [
            lambda x: x.expanding().cov(pairwise=True),
            lambda x: x.expanding().corr(pairwise=True),
            lambda x: x.rolling(window=3).cov(pairwise=True),
            lambda x: x.rolling(window=3).corr(pairwise=True),
            lambda x: x.ewm(com=3).cov(pairwise=True),
            lambda x: x.ewm(com=3).corr(pairwise=True),
        ],
    )
    def test_pairwise_with_self(self, f):

        # DataFrame with itself, pairwise=True
        # note that we may construct the 1st level of the MI
        # in a non-monotonic way, so compare accordingly
        results = []
        for i, df in enumerate(self.df1s):
            result = f(df)
            tm.assert_index_equal(result.index.levels[0], df.index, check_names=False)
            tm.assert_numpy_array_equal(
                safe_sort(result.index.levels[1]), safe_sort(df.columns.unique())
            )
            tm.assert_index_equal(result.columns, df.columns)
            results.append(df)

        for i, result in enumerate(results):
            if i > 0:
                self.compare(result, results[0])

    @pytest.mark.parametrize(
        "f",
        [
            lambda x: x.expanding().cov(pairwise=False),
            lambda x: x.expanding().corr(pairwise=False),
            lambda x: x.rolling(window=3).cov(pairwise=False),
            lambda x: x.rolling(window=3).corr(pairwise=False),
            lambda x: x.ewm(com=3).cov(pairwise=False),
            lambda x: x.ewm(com=3).corr(pairwise=False),
        ],
    )
    def test_no_pairwise_with_self(self, f):

        # DataFrame with itself, pairwise=False
        results = [f(df) for df in self.df1s]
        for (df, result) in zip(self.df1s, results):
            tm.assert_index_equal(result.index, df.index)
            tm.assert_index_equal(result.columns, df.columns)
        for i, result in enumerate(results):
            if i > 0:
                self.compare(result, results[0])

    @pytest.mark.parametrize(
        "f",
        [
            lambda x, y: x.expanding().cov(y, pairwise=True),
            lambda x, y: x.expanding().corr(y, pairwise=True),
            lambda x, y: x.rolling(window=3).cov(y, pairwise=True),
            lambda x, y: x.rolling(window=3).corr(y, pairwise=True),
            lambda x, y: x.ewm(com=3).cov(y, pairwise=True),
            lambda x, y: x.ewm(com=3).corr(y, pairwise=True),
        ],
    )
    def test_pairwise_with_other(self, f):

        # DataFrame with another DataFrame, pairwise=True
        results = [f(df, self.df2) for df in self.df1s]
        for (df, result) in zip(self.df1s, results):
            tm.assert_index_equal(result.index.levels[0], df.index, check_names=False)
            tm.assert_numpy_array_equal(
                safe_sort(result.index.levels[1]), safe_sort(self.df2.columns.unique())
            )
        for i, result in enumerate(results):
            if i > 0:
                self.compare(result, results[0])

    @pytest.mark.parametrize(
        "f",
        [
            lambda x, y: x.expanding().cov(y, pairwise=False),
            lambda x, y: x.expanding().corr(y, pairwise=False),
            lambda x, y: x.rolling(window=3).cov(y, pairwise=False),
            lambda x, y: x.rolling(window=3).corr(y, pairwise=False),
            lambda x, y: x.ewm(com=3).cov(y, pairwise=False),
            lambda x, y: x.ewm(com=3).corr(y, pairwise=False),
        ],
    )
    def test_no_pairwise_with_other(self, f):

        # DataFrame with another DataFrame, pairwise=False
        results = [
            f(df, self.df2) if df.columns.is_unique else None for df in self.df1s
        ]
        for (df, result) in zip(self.df1s, results):
            if result is not None:
                with warnings.catch_warnings(record=True):
                    warnings.simplefilter("ignore", RuntimeWarning)
                    # we can have int and str columns
                    expected_index = df.index.union(self.df2.index)
                    expected_columns = df.columns.union(self.df2.columns)
                tm.assert_index_equal(result.index, expected_index)
                tm.assert_index_equal(result.columns, expected_columns)
            else:
                with pytest.raises(ValueError, match="'arg1' columns are not unique"):
                    f(df, self.df2)
                with pytest.raises(ValueError, match="'arg2' columns are not unique"):
                    f(self.df2, df)

    @pytest.mark.parametrize(
        "f",
        [
            lambda x, y: x.expanding().cov(y),
            lambda x, y: x.expanding().corr(y),
            lambda x, y: x.rolling(window=3).cov(y),
            lambda x, y: x.rolling(window=3).corr(y),
            lambda x, y: x.ewm(com=3).cov(y),
            lambda x, y: x.ewm(com=3).corr(y),
        ],
    )
    def test_pairwise_with_series(self, f):

        # DataFrame with a Series
        results = [f(df, self.s) for df in self.df1s] + [
            f(self.s, df) for df in self.df1s
        ]
        for (df, result) in zip(self.df1s, results):
            tm.assert_index_equal(result.index, df.index)
            tm.assert_index_equal(result.columns, df.columns)
        for i, result in enumerate(results):
            if i > 0:
                self.compare(result, results[0])
