"""
Tests for DataFrame.mask; tests DataFrame.where as a side-effect.
"""

import numpy as np

from pandas import (
    NA,
    DataFrame,
    Series,
    StringDtype,
    isna,
)
import pandas._testing as tm


class TestDataFrameMask:
    def test_mask(self):
        df = DataFrame(np.random.randn(5, 3))
        cond = df > 0

        rs = df.where(cond, np.nan)
        tm.assert_frame_equal(rs, df.mask(df <= 0))
        tm.assert_frame_equal(rs, df.mask(~cond))

        other = DataFrame(np.random.randn(5, 3))
        rs = df.where(cond, other)
        tm.assert_frame_equal(rs, df.mask(df <= 0, other=other))
        tm.assert_frame_equal(rs, df.mask(~cond, other=other))

        # see GH#21891
        df = DataFrame([1, 2])
        res = df.mask([[True], [False]])

        exp = DataFrame([np.nan, 2])
        tm.assert_frame_equal(res, exp)

    def test_mask_inplace(self):
        # GH#8801
        df = DataFrame(np.random.randn(5, 3))
        cond = df > 0

        rdf = df.copy()

        return_value = rdf.where(cond, inplace=True)
        assert return_value is None
        tm.assert_frame_equal(rdf, df.where(cond))
        tm.assert_frame_equal(rdf, df.mask(~cond))

        rdf = df.copy()
        return_value = rdf.where(cond, -df, inplace=True)
        assert return_value is None
        tm.assert_frame_equal(rdf, df.where(cond, -df))
        tm.assert_frame_equal(rdf, df.mask(~cond, other=-df))

    def test_mask_edge_case_1xN_frame(self):
        # GH#4071
        df = DataFrame([[1, 2]])
        res = df.mask(DataFrame([[True, False]]))
        expec = DataFrame([[np.nan, 2]])
        tm.assert_frame_equal(res, expec)

    def test_mask_callable(self):
        # GH#12533
        df = DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = df.mask(lambda x: x > 4, other=lambda x: x + 1)
        exp = DataFrame([[1, 2, 3], [4, 6, 7], [8, 9, 10]])
        tm.assert_frame_equal(result, exp)
        tm.assert_frame_equal(result, df.mask(df > 4, other=df + 1))

        # return ndarray and scalar
        result = df.mask(lambda x: (x % 2 == 0).values, other=lambda x: 99)
        exp = DataFrame([[1, 99, 3], [99, 5, 99], [7, 99, 9]])
        tm.assert_frame_equal(result, exp)
        tm.assert_frame_equal(result, df.mask(df % 2 == 0, other=99))

        # chain
        result = (df + 2).mask(lambda x: x > 8, other=lambda x: x + 10)
        exp = DataFrame([[3, 4, 5], [6, 7, 8], [19, 20, 21]])
        tm.assert_frame_equal(result, exp)
        tm.assert_frame_equal(result, (df + 2).mask((df + 2) > 8, other=(df + 2) + 10))

    def test_mask_dtype_bool_conversion(self):
        # GH#3733
        df = DataFrame(data=np.random.randn(100, 50))
        df = df.where(df > 0)  # create nans
        bools = df > 0
        mask = isna(df)
        expected = bools.astype(object).mask(mask)
        result = bools.mask(mask)
        tm.assert_frame_equal(result, expected)

    def test_mask_pos_args_deprecation(self):
        # https://github.com/pandas-dev/pandas/issues/41485
        df = DataFrame(np.random.randn(5, 5))
        cond = df > 0
        other = DataFrame(np.random.randn(5, 3))

        msg = (
            r"Starting with Pandas version 2\.0 all arguments of mask except for the "
            r"arguments 'self' and 'cond' will be keyword-only"
        )
        with tm.assert_produces_warning(FutureWarning, match=msg):
            df.mask(cond, other)

        # tm.assert_frame_equal(df.mask(cond, other), df.mask(cond, other=other))


def test_mask_try_cast_deprecated(frame_or_series):

    obj = DataFrame(np.random.randn(4, 3))
    if frame_or_series is not DataFrame:
        obj = obj[0]

    mask = obj > 0

    with tm.assert_produces_warning(FutureWarning):
        # try_cast keyword deprecated
        obj.mask(mask, other=-1, try_cast=True)


def test_mask_stringdtype():
    # GH 40824
    df = DataFrame(
        {"A": ["foo", "bar", "baz", NA]},
        index=["id1", "id2", "id3", "id4"],
        dtype=StringDtype(),
    )
    filtered_df = DataFrame(
        {"A": ["this", "that"]}, index=["id2", "id3"], dtype=StringDtype()
    )
    filter_ser = Series([False, True, True, False])
    result = df.mask(filter_ser, other=filtered_df)

    expected = DataFrame(
        {"A": [NA, "this", "that", NA]},
        index=["id1", "id2", "id3", "id4"],
        dtype=StringDtype(),
    )
    tm.assert_frame_equal(result, expected)