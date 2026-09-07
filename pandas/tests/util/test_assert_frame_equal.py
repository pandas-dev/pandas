import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


def _assert_frame_equal_both(a, b, **kwargs):
    """
    Check that two DataFrame equal.

    This check is performed commutatively.

    Parameters
    ----------
    a : DataFrame
        The first DataFrame to compare.
    b : DataFrame
        The second DataFrame to compare.
    kwargs : dict
        The arguments passed to `tm.assert_frame_equal`.
    """
    tm.assert_frame_equal(a, b, **kwargs)
    tm.assert_frame_equal(b, a, **kwargs)


@pytest.mark.parametrize(
    "levels,codes,arrays",
    [
        (
            [["USD"], np.array([np.nan], dtype=np.float64)],
            [[0], [0]],
            [["USD"], np.array([np.nan], dtype=np.float64)],
        ),
        (
            [["USD"], np.array([np.nan, 1.0], dtype=np.float64)],
            [[0, 0], [0, 1]],
            [["USD", "USD"], np.array([np.nan, 1.0], dtype=np.float64)],
        ),
        (
            [np.array([np.nan], dtype=np.float64), ["cash"]],
            [[0], [0]],
            [np.array([np.nan], dtype=np.float64), ["cash"]],
        ),
    ],
    ids=["all-na-level", "na-and-valid-level", "na-first-level"],
)
def test_assert_frame_equal_multiindex_columns_with_na(levels, codes, arrays):
    # GH#54521
    left_columns = pd.MultiIndex(
        levels=levels,
        codes=codes,
        names=["Currency", "Collateral"],
        verify_integrity=False,
    )
    right_columns = pd.MultiIndex.from_arrays(
        arrays,
        names=["Currency", "Collateral"],
    )
    values = np.full((2, len(left_columns)), 200.0)
    left = pd.DataFrame(values, columns=left_columns)
    right = pd.DataFrame(values, columns=right_columns)

    _assert_frame_equal_both(left, right)


def test_assert_frame_equal_groupby_multiindex_columns_with_na():
    # GH#54521
    date1 = pd.Timestamp("2023-01-01")
    date2 = pd.Timestamp("2023-02-01")
    df = pd.DataFrame(
        {
            "Currency": ["USD", "USD", "USD"],
            "Collateral": [None, None, None],
            "Payment": [date1, date1, date2],
            "Cashflow": [100.0, 100.0, 200.0],
        }
    )
    result = (
        df.groupby(["Currency", "Collateral", "Payment"], dropna=False)
        .sum()
        .unstack([0, 1])
        .droplevel(0, axis=1)
    )
    expected = pd.DataFrame(
        [200.0, 200.0],
        columns=pd.MultiIndex.from_arrays(
            [
                pd.Index(["USD"], name="Currency"),
                pd.Index([None], name="Collateral", dtype="float64"),
            ]
        ),
        index=pd.DatetimeIndex([date1, date2], name="Payment"),
    )

    _assert_frame_equal_both(result, expected)


@pytest.mark.parametrize("check_like", [True, False])
def test_frame_equal_row_order_mismatch(check_like, frame_or_series):
    df1 = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=["a", "b", "c"])
    df2 = pd.DataFrame({"A": [3, 2, 1], "B": [6, 5, 4]}, index=["c", "b", "a"])

    if not check_like:  # Do not ignore row-column orderings.
        msg = f"{frame_or_series.__name__}.index are different"
        with pytest.raises(AssertionError, match=msg):
            tm.assert_frame_equal(
                df1, df2, check_like=check_like, obj=frame_or_series.__name__
            )
    else:
        _assert_frame_equal_both(
            df1, df2, check_like=check_like, obj=frame_or_series.__name__
        )


@pytest.mark.parametrize(
    "df1,df2",
    [
        ({"A": [1, 2, 3]}, {"A": [1, 2, 3, 4]}),
        ({"A": [1, 2, 3], "B": [4, 5, 6]}, {"A": [1, 2, 3]}),
    ],
)
def test_frame_equal_shape_mismatch(df1, df2, frame_or_series):
    df1 = pd.DataFrame(df1)
    df2 = pd.DataFrame(df2)
    msg = f"{frame_or_series.__name__} are different"

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2, obj=frame_or_series.__name__)


@pytest.mark.parametrize(
    "df1,df2,msg",
    [
        # Index
        (
            pd.DataFrame({"a": [1, 2], "c": ["l1", "l2"]}).set_index("a"),
            pd.DataFrame({"a": [1.0, 2.0], "c": ["l1", "l2"]}).set_index("a"),
            "DataFrame\\.index are different",
        ),
        # MultiIndex
        (
            pd.DataFrame({"a": [1, 2], "b": [2.1, 1.5], "c": ["l1", "l2"]}).set_index(
                ["a", "b"]
            ),
            pd.DataFrame(
                {"a": [1.0, 2.0], "b": [2.1, 1.5], "c": ["l1", "l2"]}
            ).set_index(["a", "b"]),
            "DataFrame\\.index level \\[0\\] are different",
        ),
    ],
)
def test_frame_equal_index_dtype_mismatch(df1, df2, msg, check_index_type):
    kwargs = {"check_index_type": check_index_type}

    if check_index_type:
        with pytest.raises(AssertionError, match=msg):
            tm.assert_frame_equal(df1, df2, **kwargs)
    else:
        tm.assert_frame_equal(df1, df2, **kwargs)


def test_empty_dtypes(check_dtype):
    columns = ["col1", "col2"]
    df1 = pd.DataFrame(columns=columns)
    df2 = pd.DataFrame(columns=columns)

    kwargs = {"check_dtype": check_dtype}
    df1["col1"] = df1["col1"].astype("int64")

    if check_dtype:
        msg = r"Attributes of DataFrame\..* are different"
        with pytest.raises(AssertionError, match=msg):
            tm.assert_frame_equal(df1, df2, **kwargs)
    else:
        tm.assert_frame_equal(df1, df2, **kwargs)


@pytest.mark.parametrize("check_like", [True, False])
def test_frame_equal_index_mismatch(check_like, frame_or_series, using_infer_string):
    if using_infer_string:
        dtype = "str"
    else:
        dtype = "object"
    msg = f"""{frame_or_series.__name__}\\.index are different

{frame_or_series.__name__}\\.index values are different \\(33\\.33333 %\\)
\\[left\\]:  Index\\(\\['a', 'b', 'c'\\], dtype='{dtype}'\\)
\\[right\\]: Index\\(\\['a', 'b', 'd'\\], dtype='{dtype}'\\)
At positional index 2, first diff: 'c' != 'd'"""

    df1 = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=["a", "b", "c"])
    df2 = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=["a", "b", "d"])

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(
            df1, df2, check_like=check_like, obj=frame_or_series.__name__
        )


@pytest.mark.parametrize("check_like", [True, False])
def test_frame_equal_columns_mismatch(check_like, frame_or_series, using_infer_string):
    if using_infer_string:
        dtype = "str"
    else:
        dtype = "object"
    msg = f"""{frame_or_series.__name__}\\.columns are different

{frame_or_series.__name__}\\.columns values are different \\(50\\.0 %\\)
\\[left\\]:  Index\\(\\['A', 'B'\\], dtype='{dtype}'\\)
\\[right\\]: Index\\(\\['A', 'b'\\], dtype='{dtype}'\\)"""

    df1 = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=["a", "b", "c"])
    df2 = pd.DataFrame({"A": [1, 2, 3], "b": [4, 5, 6]}, index=["a", "b", "c"])

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(
            df1, df2, check_like=check_like, obj=frame_or_series.__name__
        )


def test_frame_equal_block_mismatch(frame_or_series):
    obj = frame_or_series.__name__
    msg = f"""{obj}\\.iloc\\[:, 1\\] \\(column name="B"\\) are different

{obj}\\.iloc\\[:, 1\\] \\(column name="B"\\) values are different \\(33\\.33333 %\\)
\\[index\\]: \\[0, 1, 2\\]
\\[left\\]:  \\[4, 5, 6\\]
\\[right\\]: \\[4, 5, 7\\]"""

    df1 = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
    df2 = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 7]})

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2, obj=obj)


@pytest.mark.parametrize(
    "df1,df2,msg",
    [
        (
            {"A": ["á", "à", "ä"], "E": ["é", "è", "ë"]},
            {"A": ["á", "à", "ä"], "E": ["é", "è", "e̊"]},
            """{obj}\\.iloc\\[:, 1\\] \\(column name="E"\\) are different

{obj}\\.iloc\\[:, 1\\] \\(column name="E"\\) values are different \\(33\\.33333 %\\)
\\[index\\]: \\[0, 1, 2\\]
\\[left\\]:  \\[é, è, ë\\]
\\[right\\]: \\[é, è, e̊\\]""",
        ),
        (
            {"A": ["á", "à", "ä"], "E": ["é", "è", "ë"]},
            {"A": ["a", "a", "a"], "E": ["e", "e", "e"]},
            """{obj}\\.iloc\\[:, 0\\] \\(column name="A"\\) are different

{obj}\\.iloc\\[:, 0\\] \\(column name="A"\\) values are different \\(100\\.0 %\\)
\\[index\\]: \\[0, 1, 2\\]
\\[left\\]:  \\[á, à, ä\\]
\\[right\\]: \\[a, a, a\\]""",
        ),
    ],
)
def test_frame_equal_unicode(df1, df2, msg, frame_or_series):
    # see gh-20503
    #
    # Test ensures that `tm.assert_frame_equals` raises the right exception
    # when comparing DataFrames containing differing unicode objects.
    df1 = pd.DataFrame(df1)
    df2 = pd.DataFrame(df2)
    msg = msg.format(obj=frame_or_series.__name__)
    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2, obj=frame_or_series.__name__)


def test_assert_frame_equal_extension_dtype_mismatch():
    # https://github.com/pandas-dev/pandas/issues/32747
    left = pd.DataFrame({"a": [1, 2, 3]}, dtype="Int64")
    right = left.astype(int)

    msg = (
        "Attributes of DataFrame\\.iloc\\[:, 0\\] "
        '\\(column name="a"\\) are different\n\n'
        'Attribute "dtype" are different\n'
        "\\[left\\]:  Int64\n"
        "\\[right\\]: int[32|64]"
    )

    tm.assert_frame_equal(left, right, check_dtype=False)

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(left, right, check_dtype=True)


def test_assert_frame_equal_interval_dtype_mismatch():
    # https://github.com/pandas-dev/pandas/issues/32747
    left = pd.DataFrame({"a": [pd.Interval(0, 1)]}, dtype="interval")
    right = left.astype(object)

    msg = (
        "Attributes of DataFrame\\.iloc\\[:, 0\\] "
        '\\(column name="a"\\) are different\n\n'
        'Attribute "dtype" are different\n'
        "\\[left\\]:  interval\\[int64, right\\]\n"
        "\\[right\\]: object"
    )

    tm.assert_frame_equal(left, right, check_dtype=False)

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(left, right, check_dtype=True)


def test_assert_frame_equal_ignore_extension_dtype_mismatch():
    # https://github.com/pandas-dev/pandas/issues/35715
    left = pd.DataFrame({"a": [1, 2, 3]}, dtype="Int64")
    right = pd.DataFrame({"a": [1, 2, 3]}, dtype="Int32")
    tm.assert_frame_equal(left, right, check_dtype=False)


def test_assert_frame_equal_ignore_extension_dtype_mismatch_cross_class():
    # https://github.com/pandas-dev/pandas/issues/35715
    left = pd.DataFrame({"a": [1, 2, 3]}, dtype="Int64")
    right = pd.DataFrame({"a": [1, 2, 3]}, dtype="int64")
    tm.assert_frame_equal(left, right, check_dtype=False)


@pytest.mark.parametrize(
    "dtype", ["timedelta64[ns]", "datetime64[ns, UTC]", "Period[D]"]
)
def test_assert_frame_equal_datetime_like_dtype_mismatch(dtype):
    df1 = pd.DataFrame({"a": []}, dtype=dtype)
    df2 = pd.DataFrame({"a": []})
    tm.assert_frame_equal(df1, df2, check_dtype=False)


def test_allows_duplicate_labels():
    left = pd.DataFrame()
    right = pd.DataFrame().set_flags(allows_duplicate_labels=False)
    tm.assert_frame_equal(left, left)
    tm.assert_frame_equal(right, right)
    tm.assert_frame_equal(left, right, check_flags=False)
    tm.assert_frame_equal(right, left, check_flags=False)

    with pytest.raises(AssertionError, match="<Flags"):
        tm.assert_frame_equal(left, right)

    with pytest.raises(AssertionError, match="<Flags"):
        tm.assert_frame_equal(left, right)


def test_assert_frame_equal_columns_mixed_dtype():
    # GH#39168
    df = pd.DataFrame([[0, 1, 2]], columns=["foo", "bar", 42], index=[1, "test", 2])
    tm.assert_frame_equal(df, df, check_like=True)


def test_frame_equal_extension_dtype(frame_or_series, any_numeric_ea_dtype):
    # GH#39410
    obj = frame_or_series([1, 2], dtype=any_numeric_ea_dtype)
    tm.assert_equal(obj, obj, check_exact=True)


@pytest.mark.parametrize("indexer", [(0, 1), (1, 0)])
def test_frame_equal_mixed_dtypes(frame_or_series, any_numeric_ea_dtype, indexer):
    dtypes = (any_numeric_ea_dtype, "int64")
    obj1 = frame_or_series([1, 2], dtype=dtypes[indexer[0]])
    obj2 = frame_or_series([1, 2], dtype=dtypes[indexer[1]])
    tm.assert_equal(obj1, obj2, check_exact=True, check_dtype=False)


def test_assert_frame_equal_check_like_different_indexes():
    # GH#39739
    df1 = pd.DataFrame(index=pd.Index([], dtype="object"))
    df2 = pd.DataFrame(index=pd.RangeIndex(start=0, stop=0, step=1))
    with pytest.raises(AssertionError, match="DataFrame.index are different"):
        tm.assert_frame_equal(df1, df2, check_like=True)


def test_assert_frame_equal_checking_allow_dups_flag():
    # GH#45554
    left = pd.DataFrame([[1, 2], [3, 4]])
    left.flags.allows_duplicate_labels = False

    right = pd.DataFrame([[1, 2], [3, 4]])
    right.flags.allows_duplicate_labels = True
    tm.assert_frame_equal(left, right, check_flags=False)

    with pytest.raises(AssertionError, match="allows_duplicate_labels"):
        tm.assert_frame_equal(left, right, check_flags=True)


def test_assert_frame_equal_check_like_categorical_midx():
    # GH#48975
    left = pd.DataFrame(
        [[1], [2], [3]],
        index=pd.MultiIndex.from_arrays(
            [
                pd.Categorical(["a", "b", "c"]),
                pd.Categorical(["a", "b", "c"]),
            ]
        ),
    )
    right = pd.DataFrame(
        [[3], [2], [1]],
        index=pd.MultiIndex.from_arrays(
            [
                pd.Categorical(["c", "b", "a"]),
                pd.Categorical(["c", "b", "a"]),
            ]
        ),
    )
    tm.assert_frame_equal(left, right, check_like=True)


def test_assert_frame_equal_ea_column_definition_in_exception_mask():
    # GH#50323
    df1 = pd.DataFrame({"a": pd.Series([pd.NA, 1], dtype="Int64")})
    df2 = pd.DataFrame({"a": pd.Series([1, 1], dtype="Int64")})

    msg = r'DataFrame.iloc\[:, 0\] \(column name="a"\) NA mask values are different'
    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2)


def test_assert_frame_equal_ea_column_definition_in_exception():
    # GH#50323
    df1 = pd.DataFrame({"a": pd.Series([pd.NA, 1], dtype="Int64")})
    df2 = pd.DataFrame({"a": pd.Series([pd.NA, 2], dtype="Int64")})

    msg = r'DataFrame.iloc\[:, 0\] \(column name="a"\) values are different'
    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2)

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2, check_exact=True)


def test_assert_frame_equal_ts_column():
    # GH#50323
    df1 = pd.DataFrame({"a": [pd.Timestamp("2019-12-31"), pd.Timestamp("2020-12-31")]})
    df2 = pd.DataFrame({"a": [pd.Timestamp("2020-12-31"), pd.Timestamp("2020-12-31")]})

    msg = r'DataFrame.iloc\[:, 0\] \(column name="a"\) values are different'
    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2)


def test_assert_frame_equal_set():
    # GH#51727
    df1 = pd.DataFrame({"set_column": [{1, 2, 3}, {4, 5, 6}]})
    df2 = pd.DataFrame({"set_column": [{1, 2, 3}, {4, 5, 6}]})
    tm.assert_frame_equal(df1, df2)


def test_assert_frame_equal_set_mismatch():
    # GH#51727
    df1 = pd.DataFrame({"set_column": [{1, 2, 3}, {4, 5, 6}]})
    df2 = pd.DataFrame({"set_column": [{1, 2, 3}, {4, 5, 7}]})

    msg = r'DataFrame.iloc\[:, 0\] \(column name="set_column"\) values are different'
    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2)


def test_datetimelike_compat_deprecated():
    # GH#55638
    df = pd.DataFrame({"a": [1]})

    msg = "the 'check_datetimelike_compat' keyword is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        tm.assert_frame_equal(df, df, check_datetimelike_compat=True)
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        tm.assert_frame_equal(df, df, check_datetimelike_compat=False)

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        tm.assert_series_equal(df["a"], df["a"], check_datetimelike_compat=True)
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        tm.assert_series_equal(df["a"], df["a"], check_datetimelike_compat=False)


def test_by_blocks_deprecated():
    # GH#65911
    df = pd.DataFrame({"a": [1]})

    msg = "the 'by_blocks' keyword is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        tm.assert_frame_equal(df, df, by_blocks=True)
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        tm.assert_frame_equal(df, df, by_blocks=False)


def test_assert_frame_equal_int_near_bounds():
    # GH#40719 - integer comparisons near int64 bounds should be exact by default
    min_val = np.iinfo(np.int64).min
    df1 = pd.DataFrame({"B": [min_val]}, dtype=np.int64)
    df2 = pd.DataFrame({"B": [min_val + 1]}, dtype=np.int64)

    msg = r'DataFrame.iloc\[:, 0\] \(column name="B"\) values are different'
    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2)


@pytest.mark.parametrize("na_value", [pd.NA, np.nan, None])
def test_assert_frame_equal_nested_df_na(na_value):
    # GH#43022
    inner = pd.DataFrame({"a": [1, na_value]})
    df1 = pd.DataFrame({"df": [inner]})
    df2 = pd.DataFrame({"df": [inner]})
    tm.assert_frame_equal(df1, df2)


def test_assert_frame_equal_check_freq_columns():
    # GH#51920 a freq mismatch on datetimelike columns is being introduced via
    #  a deprecation: it warns by default, raises only with check_freq=True
    cols = pd.date_range("2016-01-01", periods=2)
    df1 = pd.DataFrame([[1, 2]], columns=cols)
    df2 = pd.DataFrame([[1, 2]], columns=cols._with_freq(None))

    warn_msg = "will check the 'freq' attribute"
    with tm.assert_produces_warning(Pandas4Warning, match=warn_msg):
        tm.assert_frame_equal(df1, df2)

    raise_msg = 'Attribute "freq" are different'
    with pytest.raises(AssertionError, match=raise_msg):
        tm.assert_frame_equal(df1, df2, check_freq=True)

    with tm.assert_produces_warning(None):
        tm.assert_frame_equal(df1, df2, check_freq=False)


@pytest.mark.parametrize(
    "left,right",
    [
        (
            pd.DataFrame({"a": [pd.array([1, 2, 3])]}),
            pd.DataFrame({"a": [np.array([1, 2, 3])]}),
        ),
        (
            pd.DataFrame({"a": [[1, 2, 3]]}),
            pd.DataFrame({"a": [np.array([1, 2, 3])]}),
        ),
    ],
    ids=["extensionarray-vs-ndarray", "list-vs-ndarray"],
)
def test_assert_frame_equal_nested_arraylike_type_mismatch_check_exact(left, right):
    # GH#63904
    msg = r'DataFrame.iloc\[:, 0\] \(column name="a"\) are different'

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(left, right, check_exact=True)

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(right, left, check_exact=True)


def test_assert_frame_equal_check_like_check_freq():
    # GH#51920 sorting a shuffled DatetimeIndex does not restore its freq, so
    #  the freq check is skipped with check_like=True
    idx = pd.date_range("2020-01-01", periods=3, freq="D")
    df = pd.DataFrame({"a": [1, 2, 3]}, index=idx)
    with tm.assert_produces_warning(None):
        tm.assert_frame_equal(df.iloc[[2, 0, 1]], df, check_like=True)


def test_assert_frame_equal_first_diff_quotes_strings():
    # GH#37488 string and numeric labels have identical str(), so the first-diff
    #  line was unreadable without quoting
    df1 = pd.DataFrame({"1": [1, 2, 3], "2": [1, 5, 9]})
    df2 = pd.DataFrame({1: [1, 2, 3], 2: [1, 5, 9]})
    msg = "At positional index 0, first diff: '1' != 1"

    with pytest.raises(AssertionError, match=msg):
        tm.assert_frame_equal(df1, df2, check_column_type=False)


def test_assert_frame_equal_check_freq_multiindex_level():
    # GH#66761 a freq mismatch in a MultiIndex level was not checked before
    #  the check_freq deprecation, so it warns rather than raising
    dates = pd.date_range("2012-01-01", periods=3)
    left = pd.DataFrame(
        {"a": [1, 2, 3]}, index=pd.MultiIndex.from_arrays([dates, [1, 2, 3]])
    )
    right = pd.DataFrame(
        {"a": [1, 2, 3]},
        index=pd.MultiIndex.from_arrays([dates._with_freq(None), [1, 2, 3]]),
    )

    warn_msg = "will check the 'freq' attribute"
    with tm.assert_produces_warning(Pandas4Warning, match=warn_msg):
        tm.assert_frame_equal(left, right)
    with tm.assert_produces_warning(Pandas4Warning, match=warn_msg):
        tm.assert_equal(left, right)
    with tm.assert_produces_warning(Pandas4Warning, match=warn_msg):
        tm.assert_equal(left.index, right.index)

    raise_msg = 'Attribute "freq" are different'
    with pytest.raises(AssertionError, match=raise_msg):
        tm.assert_frame_equal(left, right, check_freq=True)

    with tm.assert_produces_warning(None):
        tm.assert_frame_equal(left, right, check_freq=False)


def test_assert_frame_equal_check_freq_multiindex_repeated_codes():
    # GH#66761 a level whose codes repeat loses its freq to get_level_values, so
    #  neither side has a freq to compare and check_freq=True passes; the
    #  deprecation must not warn where the future default accepts
    dates = pd.date_range("2012-01-01", periods=3)
    left = pd.DataFrame(
        {"a": range(6)}, index=pd.MultiIndex.from_product([dates, ["a", "b"]])
    )
    right = pd.DataFrame(
        {"a": range(6)},
        index=pd.MultiIndex.from_product([dates._with_freq(None), ["a", "b"]]),
    )

    with tm.assert_produces_warning(None):
        tm.assert_frame_equal(left, right)
    tm.assert_frame_equal(left, right, check_freq=True)


def test_assert_frame_equal_check_freq_categorical_column():
    # GH#66761 the categories of a Categorical column go through the
    #  deprecation even when the frame's own index resolves check_freq to the
    #  long-standing hard check
    dates = pd.date_range("2012-01-01", periods=3)
    nofreq = dates._with_freq(None)
    left = pd.DataFrame({"a": pd.Categorical(dates, categories=dates)}, index=dates)
    right = pd.DataFrame({"a": pd.Categorical(nofreq, categories=nofreq)}, index=dates)

    warn_msg = "will check the 'freq' attribute"
    with tm.assert_produces_warning(Pandas4Warning, match=warn_msg):
        tm.assert_frame_equal(left, right)

    raise_msg = 'Attribute "freq" are different'
    with pytest.raises(AssertionError, match=raise_msg):
        tm.assert_frame_equal(left, right, check_freq=True)

    with tm.assert_produces_warning(None):
        tm.assert_frame_equal(left, right, check_freq=False)


def test_assert_frame_equal_check_freq_categorical_index():
    # GH#66761 a freq mismatch in Categorical categories was not checked before
    #  the check_freq deprecation, so it warns rather than raising
    dates = pd.date_range("2012-01-01", periods=3)
    left = pd.DataFrame({"a": [1, 2, 3]}, index=pd.CategoricalIndex(dates))
    right = pd.DataFrame(
        {"a": [1, 2, 3]}, index=pd.CategoricalIndex(dates._with_freq(None))
    )

    warn_msg = "will check the 'freq' attribute"
    with tm.assert_produces_warning(Pandas4Warning, match=warn_msg):
        tm.assert_frame_equal(left, right)

    raise_msg = 'Attribute "freq" are different'
    with pytest.raises(AssertionError, match=raise_msg):
        tm.assert_frame_equal(left, right, check_freq=True)

    with tm.assert_produces_warning(None):
        tm.assert_frame_equal(left, right, check_freq=False)
