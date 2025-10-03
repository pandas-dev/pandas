import numpy as np
import pandas as pd
import pytest


def test_groupby_sum_single_group_large_int64_matches_df_sum():
    df = pd.DataFrame({"gb": ["A", "A"], "val": pd.Series([14, 2**60], dtype="int64")})
    got = df.groupby("gb")["val"].sum().iloc[0]
    exp = df["val"].sum()
    assert got == exp
    assert df["val"].dtype == "int64"
    assert df.groupby("gb")["val"].sum().dtype == "int64"


def test_groupby_sum_multi_groups_matches_series_sum_int64():
    vals = pd.Series([2**60, 14, 2**60 + 3, 7], dtype="int64")
    gb = pd.Series(["A", "A", "B", "B"])
    df = pd.DataFrame({"gb": gb, "val": vals})

    got = df.groupby("gb")["val"].sum()
    exp = pd.Series(
        {"A": vals.iloc[:2].sum(), "B": vals.iloc[2:].sum()},
        dtype="int64",
    )
    exp.index.name = "gb"
    exp.name = "val"       # <- this aligns the Series name with got

    pd.testing.assert_series_equal(got, exp)


@pytest.mark.parametrize(
    "dtype, big, small",
    [
        ("int64", 2**60, 123),
        ("uint64", np.uint64(2**60), np.uint64(123)),
    ],
)
def test_groupby_sum_preserves_dtype_no_float_cast(dtype, big, small):
    df = pd.DataFrame(
        {"gb": ["A", "A", "B"], "val": pd.Series([big, small, big], dtype=dtype)}
    )
    out = df.groupby("gb")["val"].sum()
    assert out.dtype.name == dtype
    assert out.loc["A"] == pd.Series([big, small], dtype=dtype).sum()
    assert out.loc["B"] == big


def test_groupby_sum_nullable_uint64_min_count_behavior():
    s = pd.Series([pd.NA, np.uint64(2**60)], dtype="UInt64")
    df = pd.DataFrame({"gb": ["A", "A"], "val": s})

    out_na = df.groupby("gb")["val"].sum(min_count=2)
    assert out_na.dtype.name == "UInt64"
    assert out_na.iloc[0] is pd.NA

    out_ok = df.groupby("gb")["val"].sum(min_count=1)
    assert out_ok.dtype.name == "UInt64"
    assert out_ok.iloc[0] == np.uint64(2**60)


def test_groupby_sum_nullable_all_na_respects_min_count():
    s = pd.Series([pd.NA, pd.NA], dtype="Int64")
    df = pd.DataFrame({"gb": ["A", "A"], "val": s})
    out = df.groupby("gb")["val"].sum(min_count=1)
    assert out.dtype.name == "Int64"
    assert out.iloc[0] is pd.NA


def test_groupby_sum_dataframe_multiple_integer_columns_preserve_dtypes():
    # int64 + uint64 columns; ensure values and dtypes preserved
    df = pd.DataFrame(
        {
            "gb": ["A", "A", "B"],
            "i64": pd.Series([2**60, 5, 7], dtype="int64"),
            "u64": pd.Series(
                [np.uint64(10), np.uint64(2**54), np.uint64(3)],
                dtype="uint64",
            ),
        }
    )

    got = df.groupby("gb")[["i64", "u64"]].sum()

    exp = pd.DataFrame(
        {
            "i64": pd.Series(
                {
                    "A": pd.Series([2**60, 5], dtype="int64").sum(),
                    "B": pd.Series([7], dtype="int64").sum(),
                },
                dtype="int64",
            ),
            "u64": pd.Series(
                {
                    "A": pd.Series(
                        [np.uint64(10), np.uint64(2**54)], dtype="uint64"
                    ).sum(),
                    "B": pd.Series([np.uint64(3)], dtype="uint64").sum(),
                },
                dtype="uint64",
            ),
        }
    )
    exp.index.name = "gb"  # align index name with groupby result

    pd.testing.assert_frame_equal(got, exp)
    assert got["i64"].dtype == "int64"
    assert got["u64"].dtype == "uint64"


def test_groupby_sum_dataframe_nullable_integers_min_count_by_column():
    # Nullable Int64 / UInt64 with missing values; verify per-column min_count behavior
    df = pd.DataFrame(
        {
            "gb": ["A", "A", "A", "B"],
            "I": pd.Series([pd.NA, 2**60, pd.NA, 5], dtype="Int64"),
            "U": pd.Series([pd.NA, np.uint64(7), pd.NA, np.uint64(2)], dtype="UInt64"),
        }
    )

    out_na = df.groupby("gb")[["I", "U"]].sum(min_count=2)
    assert out_na.loc["A", "I"] is pd.NA
    assert out_na.loc["A", "U"] is pd.NA
    assert out_na.loc["B", "I"] is pd.NA
    assert out_na.loc["B", "U"] is pd.NA
    assert out_na["I"].dtype.name == "Int64"
    assert out_na["U"].dtype.name == "UInt64"

    out_ok = df.groupby("gb")[["I", "U"]].sum(min_count=1)
    assert out_ok["I"].dtype.name == "Int64"
    assert out_ok["U"].dtype.name == "UInt64"
    assert out_ok.loc["A", "I"] == 2**60
    assert out_ok.loc["A", "U"] == np.uint64(7)
    assert out_ok.loc["B", "I"] == 5
    assert out_ok.loc["B", "U"] == np.uint64(2)
