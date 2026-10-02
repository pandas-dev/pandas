from __future__ import annotations

from datetime import (
    date,
    time,
)
from decimal import Decimal
import re

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")

from pandas.core.arrays.arrow.array import ArrowExtensionArray


@pytest.mark.parametrize("pa_type", [pa.int64(), pa.float64()], ids=str)
@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize("frame", [True, False])
@pytest.mark.parametrize("chunked", [True, False])
def test_groupby_cumsum_with_nulls(pa_type, skipna, frame, chunked):
    # GH#62477: filling nulls must not overwrite non-null values on Windows.
    values = [4, 2, 3, None, 1, 2, 3]
    if chunked:
        arr = pa.chunked_array([values[:2], values[2:5], values[5:]], type=pa_type)
    else:
        arr = pa.array(values, type=pa_type)
    ser = pd.Series(ArrowExtensionArray(arr), name="a")
    obj = pd.DataFrame({"a": ser, "i": 1}) if frame else ser

    result = obj.groupby([1] * len(ser)).transform("cumsum", skipna=skipna)

    expected_values = (
        [4, 6, 9, None, 10, 12, 15] if skipna else [4, 6, 9, None, None, None, None]
    )
    expected = pd.Series(expected_values, dtype=ArrowDtype(pa_type), name="a")
    if frame:
        expected = pd.DataFrame({"a": expected, "i": range(1, len(ser) + 1)})
        tm.assert_frame_equal(result, expected)
    else:
        tm.assert_series_equal(result, expected)


class TestGroupbyAggPyArrowNative:
    """Tests for PyArrow-native groupby aggregations on decimal and string types."""

    @pytest.mark.parametrize(
        "agg_func, expected",
        [
            ("sum", [Decimal("1"), Decimal("5"), Decimal("4")]),
            ("prod", [Decimal("0"), Decimal("6"), Decimal("4")]),
            ("min", [Decimal("0"), Decimal("2"), Decimal("4")]),
            ("max", [Decimal("1"), Decimal("3"), Decimal("4")]),
            ("mean", [Decimal("0.5"), Decimal("2.5"), Decimal("4")]),
            ("count", [2, 2, 1]),
        ],
    )
    def test_groupby_decimal_aggregations(self, agg_func, expected):
        # PyArrow-native decimal groupby returns the correct values.
        values = [Decimal(str(i)) for i in range(5)]
        ser = pd.Series(values, dtype=ArrowDtype(pa.decimal128(10, 2)))
        # groups: 1 -> [0, 1], 2 -> [2, 3], 3 -> [4]
        result = ser.groupby([1, 1, 2, 2, 3]).agg(agg_func)
        assert result.index.tolist() == [1, 2, 3]
        assert isinstance(result.dtype, ArrowDtype)
        # Decimal equality is scale-insensitive (Decimal("1") == Decimal("1.00"))
        assert result.tolist() == expected

    @pytest.mark.parametrize(
        "agg_func, expected",
        [
            ("var", 0.5),
            ("std", 0.5**0.5),
            ("sem", 0.5),
        ],
    )
    def test_groupby_decimal_variance_aggregations(self, agg_func, expected):
        # std/var/sem on decimal return float64; a single-element group is NA.
        values = [Decimal(str(i)) for i in range(5)]
        ser = pd.Series(values, dtype=ArrowDtype(pa.decimal128(10, 2)))
        # groups: 1 -> [0, 1], 2 -> [2, 3], 3 -> [4] (single element -> NA)
        result = ser.groupby([1, 1, 2, 2, 3]).agg(agg_func)
        assert result.dtype == ArrowDtype(pa.float64())
        assert result.iloc[0] == pytest.approx(expected)
        assert result.iloc[1] == pytest.approx(expected)
        assert pd.isna(result.iloc[2])

    @pytest.mark.parametrize(
        "agg_func, expected",
        [
            ("min", ["a", "c", "e"]),
            ("max", ["b", "d", "e"]),
            ("count", [2, 2, 1]),
        ],
    )
    @pytest.mark.parametrize("dtype", [pa.string(), pa.large_string()])
    def test_groupby_string_aggregations(self, dtype, agg_func, expected):
        # PyArrow-native string groupby returns the correct values.
        ser = pd.Series(list("abcde"), dtype=ArrowDtype(dtype))
        # groups: 1 -> [a, b], 2 -> [c, d], 3 -> [e]
        result = ser.groupby([1, 1, 2, 2, 3]).agg(agg_func)
        assert result.index.tolist() == [1, 2, 3]
        assert isinstance(result.dtype, ArrowDtype)
        assert result.tolist() == expected

    @pytest.mark.parametrize("how", ["min", "max"])
    @pytest.mark.parametrize(
        "dtype",
        [
            pd.StringDtype("pyarrow", na_value=np.nan),
            pd.StringDtype("pyarrow", na_value=pd.NA),
            ArrowDtype(pa.string()),
            ArrowDtype(pa.large_string()),
        ],
        ids=["str[pyarrow]", "string[pyarrow]", "ArrowDtype", "ArrowDtype_large"],
    )
    def test_groupby_string_dtypes_min_max(self, dtype, how):
        # GH#63416 every PyArrow-backed string dtype takes the same path
        ser = pd.Series(["b", "a", "d", "c"], dtype=dtype)
        result = getattr(ser.groupby([1, 1, 2, 2]), how)()
        expected = ["a", "c"] if how == "min" else ["b", "d"]
        assert result.dtype == ser.dtype
        assert result.tolist() == expected

    @pytest.mark.parametrize("how", ["min", "max"])
    @pytest.mark.parametrize(
        "dtype",
        [
            pd.StringDtype("pyarrow", na_value=np.nan),
            pd.StringDtype("pyarrow", na_value=pd.NA),
            ArrowDtype(pa.string()),
            ArrowDtype(pa.large_string()),
        ],
        ids=["str[pyarrow]", "string[pyarrow]", "ArrowDtype", "ArrowDtype_large"],
    )
    def test_groupby_string_dtypes_skipna_false(self, dtype, how):
        # GH#63416 a group containing NA used to aggregate to a value, so the
        # NA was silently ignored; masked dtypes already returned NA here
        ser = pd.Series(["b", None, "d", "c"], dtype=dtype)
        result = getattr(ser.groupby([1, 1, 2, 2]), how)(skipna=False)
        assert result.dtype == ser.dtype
        assert pd.isna(result.iloc[0])
        assert result.iloc[1] == ("c" if how == "min" else "d")

    @pytest.mark.parametrize("how", ["min", "max"])
    @pytest.mark.parametrize(
        "dtype",
        [
            pd.StringDtype("pyarrow", na_value=np.nan),
            pd.StringDtype("pyarrow", na_value=pd.NA),
            ArrowDtype(pa.string()),
            ArrowDtype(pa.large_string()),
        ],
        ids=["str[pyarrow]", "string[pyarrow]", "ArrowDtype", "ArrowDtype_large"],
    )
    def test_groupby_string_dtypes_min_count(self, dtype, how):
        # GH#63416 min_count used to be ignored, so a group with fewer non-NA
        # values than min_count still got a value
        ser = pd.Series(["b", None, "d", "c"], dtype=dtype)
        result = getattr(ser.groupby([1, 1, 2, 2]), how)(min_count=2)
        assert result.dtype == ser.dtype
        assert pd.isna(result.iloc[0])
        assert result.iloc[1] == ("c" if how == "min" else "d")

    @pytest.mark.parametrize(
        "dtype,values,expected,agg_func",
        [
            (
                pa.decimal128(10, 2),
                [Decimal("1.0"), None, Decimal("3.0"), None],
                [Decimal("1.0"), Decimal("3.0")],
                "min",
            ),
            (pa.string(), ["a", None, "c", None], ["a", "c"], "min"),
            (pa.string(), ["a", None, "c", None], ["a", "c"], "max"),
        ],
    )
    def test_groupby_with_nulls(self, dtype, values, expected, agg_func):
        # Test groupby with null values.
        ser = pd.Series(values, dtype=ArrowDtype(dtype))
        result = ser.groupby([1, 1, 2, 2]).agg(agg_func)
        assert len(result) == 2
        assert result.iloc[0] == expected[0]
        assert result.iloc[1] == expected[1]

    @pytest.mark.parametrize(
        "values,keys,expected_na",
        [
            # Multiple values per group - sem is computable
            ([0, 1, 2, 3], [1, 1, 2, 2], [False, False]),
            # Single value per group - sem is NA (stddev undefined)
            ([1, 2], [1, 2], [True, True]),
            # All nulls in group 2 - sem is NA for that group
            ([1, 2, None, None], [1, 1, 2, 2], [False, True]),
        ],
    )
    def test_groupby_sem(self, values, keys, expected_na):
        # Test that sem returns float64 and handles edge cases correctly.
        ser = pd.Series(
            [Decimal(str(v)) if v is not None else None for v in values],
            dtype=ArrowDtype(pa.decimal128(10, 2)),
        )
        result = ser.groupby(keys).sem()
        assert result.dtype == ArrowDtype(pa.float64())
        assert pd.isna(result).tolist() == expected_na

    @pytest.mark.parametrize(
        "values,keys,expected_na",
        [
            # Group 1 has 2 values >= min_count, Group 2 has 1 < min_count
            ([0, 1, 2], [1, 1, 2], [False, True]),
            # With nulls: min_count uses non-null count, not group size
            # Group 1: 1 non-null < min_count=2, Group 2: 2 non-null >= min_count
            ([1, None, 2, 3, None], [1, 1, 2, 2, 2], [True, False]),
        ],
    )
    @pytest.mark.parametrize("agg_func", ["sum", "prod", "min", "max"])
    def test_groupby_min_count(self, agg_func, values, keys, expected_na):
        # Test min_count parameter with and without nulls.
        ser = pd.Series(
            [Decimal(str(v)) if v is not None else None for v in values],
            dtype=ArrowDtype(pa.decimal128(10, 2)),
        )
        result = ser.groupby(keys).agg(agg_func, min_count=2)
        assert pd.isna(result).tolist() == expected_na

    @pytest.mark.parametrize(
        "agg_func,default_value",
        [
            ("sum", 0),
            ("prod", 1),
        ],
    )
    def test_groupby_missing_groups(self, agg_func, default_value):
        # Test that missing groups get identity values.
        values = [Decimal(str(i)) for i in range(4)]
        ser = pd.Series(values, dtype=ArrowDtype(pa.decimal128(10, 2)))
        keys = pd.Categorical([0, 0, 2, 2], categories=[0, 1, 2])
        result = ser.groupby(keys, observed=False).agg(agg_func)
        assert len(result) == 3
        assert result.iloc[1] == Decimal(str(default_value))

    @pytest.mark.parametrize("agg_func", ["min", "max"])
    @pytest.mark.parametrize(
        "values, dtype",
        [
            # ordered like "badc" so the assertions below hold for both types
            (
                [Decimal("1"), Decimal("0"), Decimal("3"), Decimal("2")],
                ArrowDtype(pa.decimal128(10, 2)),
            ),
            (list("badc"), ArrowDtype(pa.string())),
            (list("badc"), pd.StringDtype("pyarrow", na_value=np.nan)),
        ],
        ids=["decimal", "ArrowDtype", "str[pyarrow]"],
    )
    def test_groupby_missing_groups_min_max(self, values, dtype, agg_func):
        # GH#63416 min and max have no identity element, so an unobserved
        # group is NA rather than filled
        ser = pd.Series(values, dtype=dtype)
        keys = pd.Categorical([0, 0, 2, 2], categories=[0, 1, 2])
        result = getattr(ser.groupby(keys, observed=False), agg_func)()
        assert len(result) == 3
        assert pd.isna(result.iloc[1])
        assert result.iloc[0] == (values[0] if agg_func == "max" else values[1])
        assert result.iloc[2] == (values[2] if agg_func == "max" else values[3])

    @pytest.mark.parametrize(
        "dropna, expected_len",
        [
            (True, 2),
            (False, 3),
        ],
    )
    def test_groupby_dropna(self, dropna, expected_len):
        # Test that NA keys are excluded when dropna=True.
        values = [Decimal(str(i)) for i in range(6)]
        ser = pd.Series(values, dtype=ArrowDtype(pa.decimal128(10, 2)))
        result = ser.groupby([1, 1, None, 2, 2, None], dropna=dropna).sum()
        assert len(result) == expected_len
        assert result.iloc[0] == Decimal("1.0")  # 0 + 1
        assert result.iloc[1] == Decimal("7.0")  # 3 + 4
        if not dropna:
            assert result.iloc[2] == Decimal("7.0")  # 2 + 5 (NA group)

    # with dropna=False the NA key sorts between the 1 and 2 groups
    @pytest.mark.parametrize(
        "dropna, expected", [(True, ["a", "d"]), (False, ["a", "c", "d"])]
    )
    @pytest.mark.parametrize(
        "dtype",
        [
            pd.StringDtype("pyarrow", na_value=np.nan),
            pd.StringDtype("pyarrow", na_value=pd.NA),
            ArrowDtype(pa.string()),
        ],
        ids=["str[pyarrow]", "string[pyarrow]", "ArrowDtype"],
    )
    def test_groupby_string_dropna(self, dtype, dropna, expected):
        # GH#63416 rows with an NA key are dropped before aggregating unless
        # they form their own group
        ser = pd.Series(list("baedc"), dtype=dtype)
        result = ser.groupby([1, 1, None, 2, None], dropna=dropna).min()
        assert result.dtype == ser.dtype
        assert result.tolist() == expected

    @pytest.mark.parametrize(
        "how", ["sum", "prod", "min", "max", "mean", "std", "var", "sem"]
    )
    def test_groupby_skipna_false(self, how):
        # GH#63416 with skipna=False, a group containing a null aggregates to NA
        values = [Decimal("1"), None, Decimal("3"), Decimal("4")]
        ser = pd.Series(values, dtype=ArrowDtype(pa.decimal128(10, 2)))
        result = getattr(ser.groupby([1, 1, 2, 2]), how)(skipna=False)
        # group 1 contains a null -> NA; group 2 has no nulls -> a real value
        assert result.iloc[0] is pd.NA
        assert result.iloc[1] is not pd.NA

    @pytest.mark.parametrize("agg_func", ["sum", "prod"])
    @pytest.mark.parametrize(
        "pa_type", [pa.decimal128(10, 2), pa.decimal256(40, 2)], ids=str
    )
    def test_groupby_sum_prod_widen_decimal_precision(self, pa_type, agg_func):
        # GH#63416 sum and prod widen to the maximum precision, matching
        # Series.sum and Series.prod, so a group result that needs more digits
        # than the input precision does not overflow
        ser = pd.Series(
            [Decimal("1"), Decimal("2"), Decimal("3")], dtype=ArrowDtype(pa_type)
        )
        result = getattr(ser.groupby([1, 1, 2]), agg_func)()
        if pa.types.is_decimal128(pa_type):
            expected_type = pa.decimal128(38, pa_type.scale)
        else:
            expected_type = pa.decimal256(76, pa_type.scale)
        assert result.dtype.pyarrow_dtype == expected_type
        # other reductions keep the input type
        assert ser.groupby([1, 1, 2]).min().dtype == ser.dtype

    @pytest.mark.parametrize(
        "agg_func, expected", [("sum", Decimal("1998")), ("prod", Decimal("998001"))]
    )
    def test_groupby_sum_prod_no_precision_overflow(self, agg_func, expected):
        # GH#63416 the group result no longer has to fit the input precision
        ser = pd.Series(
            [Decimal("999"), Decimal("999")], dtype=ArrowDtype(pa.decimal128(3, 0))
        )
        result = getattr(ser.groupby([1, 1]), agg_func)()
        assert result.iloc[0] == expected
        assert result.iloc[0] == getattr(ser, agg_func)()

    @pytest.mark.parametrize(
        "how, pa_type, values, expected",
        [
            # product needs 39 digits, one more than decimal128 can hold
            (
                "prod",
                pa.decimal128(20, 0),
                [Decimal(10**19), Decimal(12 * 10**18)],
                Decimal(12 * 10**37),
            ),
            # sum needs 39 digits, one more than decimal128 can hold
            (
                "sum",
                pa.decimal128(38, 0),
                [Decimal(3 * 10**37)] * 5,
                Decimal(15 * 10**37),
            ),
            # same, but the input already has the maximum precision, so the
            # widening cast is a no-op and cannot report the overflow
            (
                "prod",
                pa.decimal128(38, 0),
                [Decimal(10**19), Decimal(12 * 10**18)],
                Decimal(12 * 10**37),
            ),
        ],
    )
    def test_groupby_sum_prod_exceeds_max_precision(
        self, how, pa_type, values, expected
    ):
        # GH#63416 a group result that does not fit the maximum decimal
        # precision falls back to a wider type inferred from the values
        # instead of raising
        ser = pd.Series(values, dtype=ArrowDtype(pa_type))
        result = getattr(ser.groupby([1] * len(values)), how)()
        assert result.dtype == ArrowDtype(pa.decimal256(39, 0))
        assert result.iloc[0] == expected

    @pytest.mark.xfail(
        reason="PyArrow's product wraps silently once the result exceeds int256, "
        "so the group result is a wrong (negative) value; this predates the "
        "PyArrow-native path and is unchanged by it"
    )
    def test_groupby_prod_exceeds_int256(self):
        # GH#63416 a product that does not fit int256 cannot be represented by
        # any decimal type, and PyArrow reports no error for the overflow
        ser = pd.Series(
            [Decimal(10**39), Decimal(10**39)], dtype=ArrowDtype(pa.decimal256(40, 0))
        )
        result = ser.groupby([1, 1]).prod()
        assert result.iloc[0] == Decimal(10**78)

    @pytest.mark.parametrize("how", ["std", "sem"])
    def test_groupby_std_sem_supported(self, how):
        # GH#63416 these used to raise NotImplementedError on decimal
        ser = pd.Series(
            [Decimal("1"), Decimal("2"), Decimal("3"), Decimal("5")],
            dtype=ArrowDtype(pa.decimal128(10, 2)),
        )
        result = getattr(ser.groupby([1, 1, 2, 2]), how)()
        expected = getattr(pd.Series([1.0, 2.0, 3.0, 5.0]).groupby([1, 1, 2, 2]), how)()
        tm.assert_series_equal(result.astype("float64"), expected)

    @pytest.mark.parametrize("how", ["var", "std", "sem"])
    @pytest.mark.parametrize("ddof", [0, 1, 2])
    def test_groupby_decimal_ddof(self, how, ddof):
        # GH#63416 ddof is forwarded to PyArrow; ddof >= the group size is NA
        values = [Decimal("1"), Decimal("2"), Decimal("3"), Decimal("5")]
        ser = pd.Series(values, dtype=ArrowDtype(pa.decimal128(10, 2)))
        result = getattr(ser.groupby([1, 1, 2, 2]), how)(ddof=ddof)
        expected = getattr(pd.Series([1.0, 2.0, 3.0, 5.0]).groupby([1, 1, 2, 2]), how)(
            ddof=ddof
        )
        tm.assert_series_equal(result.astype("float64"), expected)

    @pytest.mark.parametrize("how", ["sum", "min", "max"])
    @pytest.mark.parametrize(
        "dtype",
        [
            ArrowDtype(pa.decimal128(10, 2)),
            ArrowDtype(pa.string()),
            pd.StringDtype("pyarrow", na_value=np.nan),
        ],
        ids=["decimal", "ArrowDtype", "str[pyarrow]"],
    )
    def test_groupby_empty(self, dtype, how):
        # GH#63416 an empty input has no groups to place results into
        ser = pd.Series([], dtype=dtype)
        result = getattr(ser.groupby([]), how)()
        assert len(result) == 0

    def test_groupby_dataframe_decimal_and_string(self):
        # GH#63416 both column types take the native path in one aggregation
        df = pd.DataFrame(
            {
                "key": [1, 1, 2, 2],
                "dec": pd.array(
                    [Decimal(str(i)) for i in range(4)],
                    dtype=ArrowDtype(pa.decimal128(10, 2)),
                ),
                "string": pd.array(list("badc"), dtype=ArrowDtype(pa.string())),
            }
        )
        result = df.groupby("key").min()
        assert result["dec"].dtype == df["dec"].dtype
        assert result["string"].dtype == df["string"].dtype
        assert result["dec"].tolist() == [Decimal("0"), Decimal("2")]
        assert result["string"].tolist() == ["a", "c"]

    @pytest.mark.parametrize(
        "dtype",
        [
            pd.StringDtype("pyarrow", na_value=np.nan),
            pd.StringDtype("pyarrow", na_value=pd.NA),
            ArrowDtype(pa.string()),
        ],
        ids=["str[pyarrow]", "string[pyarrow]", "ArrowDtype"],
    )
    def test_groupby_string_sum_falls_back(self, dtype):
        # GH#63416 PyArrow has no string sum, so it goes to the fallback path
        # and concatenates
        ser = pd.Series(["b", "a", "d", "c"], dtype=dtype)
        result = ser.groupby([1, 1, 2, 2]).sum()
        assert result.dtype == ser.dtype
        assert result.tolist() == ["ba", "dc"]

    @pytest.mark.parametrize("how", ["mean", "std", "var", "sem", "prod"])
    def test_groupby_string_unsupported_ops_raise(self, how):
        # GH#63416 the native path must not make these ops start working
        ser = pd.Series(["b", "a", "d", "c"], dtype="string[pyarrow]")
        with pytest.raises(TypeError, match=f"does not support operation '{how}'"):
            getattr(ser.groupby([1, 1, 2, 2]), how)()

    # ordered like "badc" for every type
    _dates = [date(2020, 1, 2), date(2020, 1, 1), date(2021, 1, 4), date(2021, 1, 3)]
    _times = [time(1), time(0), time(3), time(2)]
    _bytes = [b"b", b"a", b"d", b"c"]
    _temporal_binary_cases = [
        pytest.param(pa_type, values, id=str(pa_type))
        for pa_type, values in [
            (pa.date32(), _dates),
            (pa.date64(), _dates),
            (pa.time32("s"), _times),
            (pa.time32("ms"), _times),
            (pa.time64("us"), _times),
            (pa.time64("ns"), _times),
            (pa.binary(), _bytes),
            (pa.large_binary(), _bytes),
            (pa.binary(1), _bytes),
        ]
    ]

    @pytest.mark.parametrize("how", ["min", "max"])
    @pytest.mark.parametrize("pa_type, values", _temporal_binary_cases)
    def test_groupby_temporal_binary_min_max(self, pa_type, values, how):
        # GH#66626
        ser = pd.Series([*values, None], dtype=ArrowDtype(pa_type))
        keys = pd.Categorical([0, 0, 2, 2, 2], categories=[0, 1, 2])
        result = getattr(ser.groupby(keys, observed=False), how)()
        expected = pd.Series(
            [values[0], None, values[2]]
            if how == "max"
            else [values[1], None, values[3]],
            index=pd.CategoricalIndex([0, 1, 2]),
            dtype=ArrowDtype(pa_type),
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("how", ["min", "max"])
    @pytest.mark.parametrize("pa_type, values", _temporal_binary_cases)
    def test_groupby_temporal_binary_all_na_keeps_dtype(self, pa_type, values, how):
        # GH#66626 an all-NA result keeps the column's dtype
        ser = pd.Series([values[0], None, values[1], None], dtype=ArrowDtype(pa_type))
        result = getattr(ser.groupby([0, 0, 1, 1]), how)(skipna=False)
        expected = pd.Series([None, None], index=[0, 1], dtype=ArrowDtype(pa_type))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("how", ["min", "max"])
    @pytest.mark.parametrize("pa_type, values", _temporal_binary_cases)
    def test_groupby_temporal_binary_min_max_min_count(self, pa_type, values, how):
        # GH#66626 min_count is respected
        ser = pd.Series(
            [values[0], None, values[2], values[3]], dtype=ArrowDtype(pa_type)
        )
        result = getattr(ser.groupby([0, 0, 1, 1]), how)(min_count=2)
        expected = pd.Series(
            [None, values[2] if how == "max" else values[3]],
            index=[0, 1],
            dtype=ArrowDtype(pa_type),
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("pa_type, values", _temporal_binary_cases)
    def test_groupby_temporal_binary_min_chunked_unsorted(self, pa_type, values):
        # GH#66626
        arr = pa.chunked_array(
            [pa.array(values[:2], pa_type), pa.array(values[2:], pa_type)]
        )
        ser = pd.Series(ArrowExtensionArray(arr))
        result = ser.groupby([2, None, 1, 2], sort=False).min()
        expected = pd.Series(
            [values[0], values[2]], index=[2, 1], dtype=ArrowDtype(pa_type)
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype",
        [
            ArrowDtype(pa.string()),
            pd.StringDtype("pyarrow"),
            ArrowDtype(pa.binary()),
            ArrowDtype(pa.large_binary()),
        ],
        ids=["string", "str[pyarrow]", "binary", "large_binary"],
    )
    def test_groupby_min_count_many_groups(self, dtype):
        # GH#66626 a chunked aggregate with PyArrow < 24 (GH#64320)
        ngroups = 40_000
        values = [f"v{i}" for i in range(2 * ngroups)]
        if dtype == ArrowDtype(pa.binary()) or dtype == ArrowDtype(pa.large_binary()):
            values = [value.encode() for value in values]
        ser = pd.Series(values, dtype=dtype)
        keys = np.arange(2 * ngroups) % ngroups
        result = ser.groupby(keys).min(min_count=1)
        result.array._pa_array.validate(full=True)
        tm.assert_series_equal(result, ser.groupby(keys).min())

    @pytest.mark.parametrize("min_count", [0, 1, 3])
    @pytest.mark.parametrize("skipna", [True, False])
    @pytest.mark.parametrize("how", ["min", "max"])
    def test_groupby_null_min_max(self, how, skipna, min_count):
        # GH#66626
        ser = pd.Series(ArrowExtensionArray(pa.nulls(4)))
        keys = pd.Categorical([0, 0, 2, 2], categories=[0, 1, 2])
        result = getattr(ser.groupby(keys, observed=False), how)(
            skipna=skipna, min_count=min_count
        )
        expected = pd.Series(
            ArrowExtensionArray(pa.nulls(3)), index=pd.CategoricalIndex([0, 1, 2])
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("how", ["min", "max"])
    def test_groupby_time64_ns_min_count_keeps_nanoseconds(self, how):
        # GH#66626 nanoseconds survive min_count
        arr = pa.array([1, 2, 3], pa.int64()).cast(pa.time64("ns"))
        ser = pd.Series(ArrowExtensionArray(arr))
        result = getattr(ser.groupby([0, 0, 1]), how)(min_count=2)
        expected = pd.Series(
            ArrowExtensionArray(
                pa.array([1 if how == "min" else 2, None]).cast(pa.time64("ns"))
            )
        )
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("frame", [True, False])
@pytest.mark.parametrize("how", ["any", "all", "std", "sem", "idxmin", "idxmax"])
@pytest.mark.parametrize(
    "arr",
    [
        pa.array([date(2020, 1, 2), date(2020, 1, 1), None, date(2020, 1, 3)]),
        pa.array([time(2), time(1), None, time(3)]),
        pa.array([b"b", b"a", None, b"c"]),
        pa.array(["b", "a", None, "c"]),
        pa.array(["b", "a", None, "c"]).dictionary_encode(),
        pa.array([[2], [1], None, [3]]),
    ],
    ids=lambda arr: str(arr.type),
)
def test_groupby_unsupported_op_raises_typeerror(arr, how, frame):
    # GH#69717 used to raise NotImplementedError, mostly with no message
    ser = pd.Series(ArrowExtensionArray(arr))
    obj = ser.to_frame() if frame else ser
    msg = f"{how} is not supported for {re.escape(str(ser.dtype))} dtype"
    with pytest.raises(TypeError, match=msg):
        getattr(obj.groupby([0, 0, 1, 1]), how)()


@pytest.mark.parametrize("op_name", ["var", "std", "sem", "mean"])
@pytest.mark.parametrize("dtype", ["int64[pyarrow]", "float64[pyarrow]"])
def test_groupby_cython_agg_pyarrow_dtype_retention(op_name, dtype):
    # GH#54627
    arr = pd.array([1, 2, 3, 4], dtype=dtype)
    df = pd.DataFrame({"key": ["a", "a", "b", "b"], "col": arr})
    grouped = df.groupby("key")
    expected_dtype = ArrowDtype(pa.float64())

    result = getattr(grouped, op_name)()
    assert result["col"].dtype == expected_dtype

    result = grouped.aggregate(op_name)
    assert result["col"].dtype == expected_dtype

    result = getattr(grouped["col"], op_name)()
    assert result.dtype == expected_dtype

    result = grouped["col"].aggregate(op_name)
    assert result.dtype == expected_dtype


def test_groupby_series_size_returns_pa_int(data):
    # GH 54132
    ser = pd.Series(data[:3], index=["a", "a", "b"])
    result = ser.groupby(level=0).size()
    expected = pd.Series([2, 1], dtype="int64[pyarrow]", index=["a", "b"])
    tm.assert_series_equal(result, expected)


def test_groupby_count_return_arrow_dtype(data_missing):
    df = pd.DataFrame({"A": [1, 1], "B": data_missing, "C": data_missing})
    result = df.groupby("A").count()
    expected = pd.DataFrame(
        [[1, 1]],
        index=pd.Index([1], name="A"),
        columns=["B", "C"],
        dtype="int64[pyarrow]",
    )
    tm.assert_frame_equal(result, expected)
