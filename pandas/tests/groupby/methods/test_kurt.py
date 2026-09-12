import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


def test_groupby_kurt_equivalence():
    # GH#40139
    # Test that the groupby kurt method (which uses libgroupby.group_kurt)
    #  matches the results of operating group-by-group (which uses nanops.nankurt)
    nrows = 1000
    ngroups = 3
    ncols = 2
    nan_frac = 0.05

    arr = np.random.default_rng(2).standard_normal((nrows, ncols))
    arr[np.random.default_rng(2).random(nrows) < nan_frac] = np.nan

    df = pd.DataFrame(arr)
    grps = np.random.default_rng(2).integers(0, ngroups, size=nrows)
    gb = df.groupby(grps)

    result = gb.kurt()

    grpwise = [grp.kurt().to_frame(i).T for i, grp in gb]
    expected = pd.concat(grpwise, axis=0)
    expected.index = expected.index.astype("int64")  # 32bit builds
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "dtype",
    [
        pytest.param("float64[pyarrow]", marks=td.skip_if_no("pyarrow")),
        "Float64",
    ],
)
def test_groupby_kurt_arrow_float64(dtype):
    # GH#40139
    # Test groupby.kurt() with float64[pyarrow] and Float64 dtypes
    df = pd.DataFrame(
        {
            "x": [1.0, pd.NA, 3.2, 4.8, 2.3, 1.9, 8.9],
            "y": [1.6, 3.3, 3.2, 6.8, 1.3, 2.9, 9.0],
        },
        dtype=dtype,
    )
    gb = df.groupby(by=lambda x: 0)

    result = gb.kurt()
    expected = pd.DataFrame({"x": [2.1644713], "y": [0.1513969]}, dtype=dtype)
    tm.assert_almost_equal(result, expected)


def test_groupby_kurt_noskipna():
    # GH#40139
    # Test groupby.kurt() with skipna = False
    df = pd.DataFrame(
        {
            "x": [1.0, np.nan, 3.2, 4.8, 2.3, 1.9, 8.9],
            "y": [1.6, 3.3, 3.2, 6.8, 1.3, 2.9, 9.0],
        }
    )
    gb = df.groupby(by=lambda x: 0)

    result = gb.kurt(skipna=False)
    expected = pd.DataFrame({"x": [np.nan], "y": [0.1513969]})
    tm.assert_almost_equal(result, expected)


def test_groupby_kurt_all_ones():
    # GH#40139
    # Test groupby.kurt() with constant values
    df = pd.DataFrame(
        {
            "x": [1.0] * 10,
        }
    )
    gb = df.groupby(by=lambda x: 0)

    result = gb.kurt(skipna=False)
    expected = pd.DataFrame(
        {
            "x": [np.nan],  # Same behavior as pd.DataFrame.kurt()
        }
    )
    tm.assert_almost_equal(result, expected)


@pytest.mark.parametrize("bias", [True, False])
def test_groupby_kurt_bias(bias):
    sp_stats = pytest.importorskip("scipy.stats")

    df = pd.DataFrame({"g": ["a", "a", "a", "a", "a"], "v": [1.0, 2.0, 2.0, 3.0, 10.0]})
    result = df.groupby("g")["v"].kurt(bias=bias)
    expected = sp_stats.kurtosis(df["v"], bias=bias)
    tm.assert_almost_equal(result.iloc[0], expected)


@pytest.mark.parametrize("bias", [True, False])
def test_groupby_kurt_bias_mixed_group_sizes(bias):
    # GH#54556, GH#66659: each group is reduced independently. With bias=False a
    # group with fewer than 4 observations is NaN; with bias=True the nobs-based
    # NaN gate in calc_kurt no longer applies, so a small non-degenerate group
    # computes a real value while a zero-variance group is still NaN.
    sp_stats = pytest.importorskip("scipy.stats")

    df = pd.DataFrame(
        {
            "g": ["small"] * 3 + ["const"] * 3 + ["big"] * 5,
            "v": [1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 1.0, 2.0, 2.0, 3.0, 10.0],
        }
    )
    result = df.groupby("g")["v"].kurt(bias=bias)

    small = sp_stats.kurtosis([1.0, 2.0, 3.0], bias=True) if bias else np.nan
    expected = pd.Series(
        [sp_stats.kurtosis([1.0, 2.0, 2.0, 3.0, 10.0], bias=bias), np.nan, small],
        index=pd.Index(["big", "const", "small"], name="g"),
        name="v",
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("dtype", ["Float64", "Int64"])
@pytest.mark.parametrize("group", ["allna", "one", "const"])
@pytest.mark.parametrize("bias", [True, False])
def test_groupby_kurt_bias_degenerate_group_is_na(bias, group, dtype, using_nan_is_na):
    # GH#54556, GH#66659: a group with zero second moment has undefined kurt
    #  for either bias, so a nullable result reports NA, not a bare NaN.
    values = {
        "allna": [pd.NA] * 4,
        "one": [5, pd.NA, pd.NA, pd.NA],
        "const": [4, 4, 4, 4],
    }[group]
    df = pd.DataFrame({"g": [group] * 4, "v": pd.array(values, dtype=dtype)})

    result = df.groupby("g")["v"].kurt(bias=bias)

    expected = pd.Series(
        pd.array([pd.NA], dtype="Float64"),
        index=pd.Index([group], name="g"),
        name="v",
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("bias", [True, False])
def test_groupby_kurt_bias_unmasked_nan_stays_nan(bias, using_nan_is_na):
    # GH#66659: a NaN the mask marks as not-NA is a value, not missingness,
    #  so it propagates unmasked instead of becoming NA.
    values = np.array([1.0, 2.0, 3.0, 4.0, np.nan], dtype="float64")
    mask = np.zeros(5, dtype="bool")
    ser = pd.Series(pd.arrays.FloatingArray(values, mask))
    df = pd.DataFrame({"g": [1] * 5, "v": ser})

    result = df.groupby("g")["v"].kurt(bias=bias)

    expected = pd.Series(
        pd.arrays.FloatingArray(
            np.array([np.nan]), np.array([using_nan_is_na], dtype="bool")
        ),
        index=pd.Index([1], name="g"),
        name="v",
    )
    tm.assert_series_equal(result, expected)
