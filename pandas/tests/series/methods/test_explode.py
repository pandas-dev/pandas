import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_basic():
    s = pd.Series([[0, 1, 2], np.nan, [], (3, 4)], index=list("abcd"), name="foo")
    result = s.explode()
    expected = pd.Series(
        [0, 1, 2, np.nan, np.nan, 3, 4], index=list("aaabcdd"), dtype=object, name="foo"
    )
    tm.assert_series_equal(result, expected)


def test_mixed_type():
    s = pd.Series(
        [[0, 1, 2], np.nan, None, np.array([]), pd.Series(["a", "b"])], name="foo"
    )
    result = s.explode()
    expected = pd.Series(
        [0, 1, 2, np.nan, None, np.nan, "a", "b"],
        index=[0, 0, 0, 1, 2, 3, 4, 4],
        dtype=object,
        name="foo",
    )
    tm.assert_series_equal(result, expected)


def test_empty():
    s = pd.Series(dtype=object)
    result = s.explode()
    expected = s.copy()
    tm.assert_series_equal(result, expected)


def test_nested_lists():
    s = pd.Series([[[1, 2, 3]], [1, 2], 1])
    result = s.explode()
    expected = pd.Series([[1, 2, 3], 1, 2, 1], index=[0, 1, 1, 2])
    tm.assert_series_equal(result, expected)


def test_multi_index():
    s = pd.Series(
        [[0, 1, 2], np.nan, [], (3, 4)],
        name="foo",
        index=pd.MultiIndex.from_product([list("ab"), range(2)], names=["foo", "bar"]),
    )
    result = s.explode()
    index = pd.MultiIndex.from_tuples(
        [("a", 0), ("a", 0), ("a", 0), ("a", 1), ("b", 0), ("b", 1), ("b", 1)],
        names=["foo", "bar"],
    )
    expected = pd.Series(
        [0, 1, 2, np.nan, np.nan, 3, 4], index=index, dtype=object, name="foo"
    )
    tm.assert_series_equal(result, expected)


def test_large():
    s = pd.Series([range(256)]).explode()
    result = s.explode()
    tm.assert_series_equal(result, s)


def test_invert_array():
    df = pd.DataFrame({"a": pd.date_range("20190101", periods=3, tz="UTC")})

    listify = df.apply(lambda x: x.array, axis=1)
    result = listify.explode()
    tm.assert_series_equal(result, df["a"].rename())


@pytest.mark.parametrize(
    "data", [[1, 2, 3], pd.date_range("2019", periods=3, tz="UTC")]
)
def test_non_object_dtype(data):
    ser = pd.Series(data)
    result = ser.explode()
    tm.assert_series_equal(result, ser)


def test_typical_usecase():
    df = pd.DataFrame(
        [{"var1": "a,b,c", "var2": 1}, {"var1": "d,e,f", "var2": 2}],
        columns=["var1", "var2"],
    )
    exploded = df.var1.str.split(",").explode()
    result = df[["var2"]].join(exploded)
    expected = pd.DataFrame(
        {"var2": [1, 1, 1, 2, 2, 2], "var1": list("abcdef")},
        columns=["var2", "var1"],
        index=[0, 0, 0, 1, 1, 1],
    )
    tm.assert_frame_equal(result, expected)


def test_nested_EA():
    # a nested EA array
    s = pd.Series(
        [
            pd.date_range("20170101", periods=3, tz="UTC"),
            pd.date_range("20170104", periods=3, tz="UTC"),
        ]
    )
    result = s.explode()
    expected = pd.Series(
        pd.date_range("20170101", periods=6, tz="UTC"), index=[0, 0, 0, 1, 1, 1]
    )
    tm.assert_series_equal(result, expected)


def test_duplicate_index():
    # GH 28005
    s = pd.Series([[1, 2], [3, 4]], index=[0, 0])
    result = s.explode()
    expected = pd.Series([1, 2, 3, 4], index=[0, 0, 0, 0], dtype=object)
    tm.assert_series_equal(result, expected)


def test_ignore_index():
    # GH 34932
    s = pd.Series([[1, 2], [3, 4]])
    result = s.explode(ignore_index=True)
    expected = pd.Series([1, 2, 3, 4], index=[0, 1, 2, 3], dtype=object)
    tm.assert_series_equal(result, expected)


def test_explode_sets():
    # https://github.com/pandas-dev/pandas/issues/35614
    s = pd.Series([{"a", "b", "c"}], index=[1])
    result = s.explode().sort_values()
    expected = pd.Series(["a", "b", "c"], index=[1, 1, 1])
    tm.assert_series_equal(result, expected)


def test_explode_scalars_can_ignore_index():
    # https://github.com/pandas-dev/pandas/issues/40487
    s = pd.Series([1, 2, 3], index=["a", "b", "c"])
    result = s.explode(ignore_index=True)
    expected = pd.Series([1, 2, 3])
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ignore_index", [True, False])
@pytest.mark.parametrize("list_type", ["list_", "large_list"])
def test_explode_pyarrow_list_type(ignore_index, list_type):
    # GH 53602, 61091
    pa = pytest.importorskip("pyarrow")

    data = [
        [None, None],
        [1],
        [],
        [2, 3],
        None,
    ]
    ser = pd.Series(data, dtype=pd.ArrowDtype(getattr(pa, list_type)(pa.int64())))
    result = ser.explode(ignore_index=ignore_index)
    expected = pd.Series(
        data=[None, None, 1, None, 2, 3, None],
        index=None if ignore_index else [0, 0, 1, 2, 3, 3, 4],
        dtype=pd.ArrowDtype(pa.int64()),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ignore_index", [True, False])
def test_explode_pyarrow_non_list_type(ignore_index):
    pa = pytest.importorskip("pyarrow")
    data = [1, 2, 3]
    ser = pd.Series(data, dtype=pd.ArrowDtype(pa.int64()))
    result = ser.explode(ignore_index=ignore_index)
    expected = pd.Series([1, 2, 3], dtype="int64[pyarrow]", index=[0, 1, 2])
    tm.assert_series_equal(result, expected)


def test_explode_preserves_datetime_unit():
    # Create datetime64[ms] array manually
    dt64_ms = np.array(
        [
            "2020-01-01T00:00:00.000",
            "2020-01-01T01:00:00.000",
            "2020-01-01T02:00:00.000",
        ],
        dtype="datetime64[ms]",
    )
    s = pd.Series([dt64_ms])

    # Explode the Series
    result = s.explode()

    # Ensure the dtype (including unit) is preserved
    assert result.dtype == dt64_ms.dtype, (
        f"Expected dtype {dt64_ms.dtype}, got {result.dtype}"
    )


def test_single_column_explode_preserves_datetime_unit():
    # Use freq in ms since unit='ms'
    rng = pd.date_range("2020-01-01T00:00:00Z", periods=3, freq="3600000ms", unit="ms")
    s = pd.Series([rng])
    result = s.explode()
    assert result.dtype == rng.dtype


def test_multi_column_explode_preserves_datetime_unit():
    rng1 = pd.date_range("2020-01-01", periods=2, freq="3600000ms", unit="ms")
    rng2 = pd.date_range("2020-01-01", periods=2, freq="3600000ms", unit="ms")
    df = pd.DataFrame({"A": [rng1], "B": [rng2]})
    result = df.explode(["A", "B"])
    assert result["A"].dtype == rng1.dtype
    assert result["B"].dtype == rng2.dtype
