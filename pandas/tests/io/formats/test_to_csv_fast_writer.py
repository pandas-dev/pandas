"""
Differential tests for the fast to_csv chunk writer
(libwriters.write_csv_chunk): output must be byte-identical to the legacy
csv.writer path for every dialect the fast writer claims.
"""

import csv

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame

from pandas.io.formats.csvs import CSVFormatter


def _legacy_to_csv(monkeypatch, df, **kwargs):
    with monkeypatch.context() as ctx:
        ctx.setattr(CSVFormatter, "_use_fast_writer", property(lambda self: False))
        return df.to_csv(**kwargs)


def assert_fast_matches_legacy(monkeypatch, df, **kwargs):
    # get_values_for_csv casts some ExtensionArrays (e.g. IntervalArray) to
    # str, which can emit a spurious "invalid value encountered in cast"
    # RuntimeWarning on some platforms; it is shared by both writers and
    # irrelevant to the byte comparison here.
    with np.errstate(invalid="ignore"):
        result = df.to_csv(**kwargs)
        expected = _legacy_to_csv(monkeypatch, df, **kwargs)
    assert result == expected


@pytest.fixture
def mixed_frame():
    return DataFrame(
        {
            "f8": [1.5, np.nan, -np.inf, 1e-4, 9.999999999999999e-05, 1e16, -0.0],
            "f4": np.array([0.1, np.nan, 1e6, 999999.94, 1e-4, 5e-324, 2.5], "f4"),
            "i8": [0, -1, 2**63 - 1, -(2**63), 3, 4, 5],
            "u8": np.array([0, 1, 2**64 - 1, 2, 3, 4, 5], dtype=np.uint64),
            "i4": np.array([1, 2, 3, 4, 5, 6, 7], dtype=np.int32),
            "b": [True, False] * 3 + [True],
            "dt": pd.to_datetime(
                [
                    "2020-01-01",
                    "NaT",
                    "2020-01-02 03:04:05",
                    "2020-01-02 03:04:05.123",
                    "2020-01-02 03:04:05.123456",
                    "2020-01-02 03:04:05.123456789",
                    "1677-09-21 00:12:43.145224193",
                ],
                format="ISO8601",
            ),
            "s": ["plain", "with,comma", 'with"quote', "with\nnewline", "", "x", "y"],
            "o": np.array(
                ["a", None, np.nan, 3, 2.5, pd.Timestamp("2020-01-01"), ("t", "u")],
                dtype=object,
            ),
        }
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"index": False},
        {"header": False},
        {"sep": ";"},
        {"sep": "."},
        {"sep": " "},
        {"sep": "0"},
        {"quotechar": "'"},
        {"lineterminator": "\r\n"},
        {"lineterminator": "\r"},
        {"na_rep": "NA"},
        {"na_rep": "N,A"},
        {"na_rep": 'N"A'},
        {"chunksize": 2},
        {"chunksize": 3, "na_rep": "NULL", "sep": ";"},
        {"float_format": "%.3f"},
        {"decimal": ","},
        {"date_format": "%Y/%m/%d %H:%M"},
    ],
)
def test_fast_writer_matches_legacy_options(monkeypatch, mixed_frame, kwargs):
    assert_fast_matches_legacy(monkeypatch, mixed_frame, **kwargs)


@pytest.mark.parametrize(
    "data",
    [
        pd.array([1, None, -3], dtype="Int64"),
        pd.array([1.5, None, -3.25], dtype="Float64"),
        pd.array([True, None, False], dtype="boolean"),
        pd.array(["x", "y,z", None], dtype="str"),
        pd.Categorical(["x", "y,z", None]),
        pd.Categorical(pd.to_datetime(["2020-01-01", "NaT", "2021-02-03"])),
        pd.interval_range(0.0, 3.0, periods=3),
        pd.period_range("2020-01", periods=3, freq="M"),
        pd.to_timedelta(["1 days 02:03:04", "NaT", "-1 days"]),
        pd.to_datetime(
            ["2020-01-01 01:02:03", "NaT", "2021-01-01"], format="ISO8601"
        ).tz_localize("US/Pacific"),
        np.array([1 + 2j, np.nan + 0j, -3.5j]),
        np.array([1.5, np.nan, 2.5], dtype=np.float16),
        np.array(["0045-01-02", "0999-06-07"], dtype="M8[s]"),
        np.array(["0045-01-02T03:04:05", "0999-06-07T08:09:10"], dtype="M8[s]"),
        pd.DatetimeIndex([pd.NaT, pd.NaT]).as_unit("us"),
    ],
)
def test_fast_writer_matches_legacy_dtypes(monkeypatch, data):
    df = DataFrame({"a": data})
    assert_fast_matches_legacy(monkeypatch, df)
    assert_fast_matches_legacy(monkeypatch, df, na_rep="NAVAL")


@pytest.mark.parametrize(
    "index",
    [
        pd.Index([0.5, 1e16, np.nan], name="fi"),
        pd.Index(["a,b", 'c"d', "e"]),
        pd.date_range("2020-01-01", periods=3, freq="12h"),
        pd.RangeIndex(10, 40, 10),
        pd.MultiIndex.from_tuples([("a", 1), ("b", 2), ("c,x", 3)]),
        pd.MultiIndex.from_arrays([["p", "q,r", "s"]], names=["only"]),
        pd.CategoricalIndex(["x", "y", None]),
    ],
)
def test_fast_writer_matches_legacy_index(monkeypatch, index):
    df = DataFrame({"val": [1.5, 2.5, np.nan]}, index=index)
    assert_fast_matches_legacy(monkeypatch, df)
    assert_fast_matches_legacy(monkeypatch, df, chunksize=2)


def test_fast_writer_float_bit_patterns(monkeypatch):
    # shortest-round-trip float formatting must match numpy's str() exactly
    rng = np.random.default_rng(20260701)
    f8 = rng.integers(0, 2**64, 10_000, dtype=np.uint64).view(np.float64)
    assert_fast_matches_legacy(monkeypatch, DataFrame({"a": f8}), index=False)
    f4 = rng.integers(0, 2**32, 10_000, dtype=np.uint32).view(np.float32)
    assert_fast_matches_legacy(monkeypatch, DataFrame({"a": f4}), index=False)


def test_fast_writer_float_boundaries(monkeypatch):
    # numpy's fixed/scientific cutoff is value-based: [1e-4, 1e16) for
    # float64, [1e-4, 1e6) for float32
    specials = [
        0.0,
        -0.0,
        np.inf,
        -np.inf,
        np.nan,
        5e-324,
        2.2250738585072014e-308,
        1.7976931348623157e308,
        1e-4,
        np.nextafter(1e-4, 0),
        9.999999999999999e-05,
        1e16,
        np.nextafter(1e16, 0),
        1e15,
        0.1,
        1 / 3,
        1e100,
        1e-100,
    ]
    assert_fast_matches_legacy(monkeypatch, DataFrame({"a": specials}))
    f4_specials = np.array(
        [0.0, -0.0, 1e-4, 1e6, 999999.94, 16777217.0, 1e-45, 3.4028235e38],
        dtype=np.float32,
    )
    assert_fast_matches_legacy(monkeypatch, DataFrame({"a": f4_specials}))


def test_fast_writer_lone_empty_field(monkeypatch):
    # csv.writer quotes an empty field iff it is the row's only field
    # (lineterminator is pinned so the expected string is platform-independent)
    df = DataFrame({"a": ["", "x", ""]})
    assert_fast_matches_legacy(monkeypatch, df, index=False)
    assert df.to_csv(index=False, lineterminator="\n") == 'a\n""\nx\n""\n'
    # ... including for NaN rendered with the default na_rep=""
    df = DataFrame({"a": [np.nan, 1.0]})
    assert_fast_matches_legacy(monkeypatch, df, index=False)
    assert df.to_csv(index=False, lineterminator="\n") == 'a\n""\n1.0\n'


@pytest.mark.parametrize("lineterminator", ["\n", "\r\n", "\r"])
def test_fast_writer_embedded_newline_quoting(monkeypatch, lineterminator):
    # csv.writer quotes a field containing \r or \n; up to CPython 3.11.1 it
    # did so only when the char was part of the line terminator. The fast
    # writer must match the running module for each line terminator.
    df = DataFrame({"a": ["plain", "has\nlf", "has\rcr", "has\r\ncrlf"]})
    assert_fast_matches_legacy(monkeypatch, df, lineterminator=lineterminator)
    assert_fast_matches_legacy(
        monkeypatch, df, lineterminator=lineterminator, na_rep="a\nb"
    )


def test_fast_writer_long_na_rep_truncation(monkeypatch):
    # the legacy float path assigns na_rep into a '<U32' array, truncating
    # longer values; the fast writer diverts to that path to match
    df = DataFrame({"a": [1.5, np.nan]})
    for length in [32, 33, 40]:
        assert_fast_matches_legacy(monkeypatch, df, na_rep="M" * length)


def test_fast_writer_non_str_na_rep(monkeypatch):
    df = DataFrame({"a": [1.5, np.nan], "b": [pd.NaT, pd.Timestamp("2020-01-01")]})
    assert_fast_matches_legacy(monkeypatch, df, na_rep=999)


def test_fast_writer_surrogate_fallback(monkeypatch):
    # PyUnicode_AsUTF8AndSize rejects lone surrogates; the chunk falls back
    # to the legacy writer, which passes them through
    df = DataFrame({"a": ["ok"] * 3 + ["bad\ud800esc"] + ["ok"] * 3}, dtype=object)
    assert_fast_matches_legacy(monkeypatch, df)
    assert_fast_matches_legacy(monkeypatch, df, chunksize=2)
    assert "\ud800" in df.to_csv()


def test_fast_writer_dt64_per_chunk_resolution(monkeypatch):
    # the fractional-digits and dates-only decisions are made per chunk,
    # matching the legacy conversion
    df = DataFrame(
        {
            "a": pd.to_datetime(
                [
                    "2020-01-01",
                    "2020-01-02",
                    "2020-01-03 00:00:00.5",
                    "2020-01-04",
                    "2020-01-05 00:00:00.000001",
                    "2020-01-06",
                ],
                format="ISO8601",
            )
        }
    )
    for chunksize in [1, 2, 3, None]:
        assert_fast_matches_legacy(monkeypatch, df, chunksize=chunksize)
    result = df.to_csv(chunksize=2).splitlines()
    # rows 0-1 render dates-only, rows 2-3 with (ms-resolution) fractional
    # seconds, rows 4-5 with us-resolution fractional seconds
    assert result[1] == "0,2020-01-01"
    assert result[3] == "2,2020-01-03 00:00:00.500"
    assert result[4] == "3,2020-01-04 00:00:00.000"
    assert result[5] == "4,2020-01-05 00:00:00.000001"


def test_fast_writer_used_for_default_dialect(monkeypatch, mixed_frame):
    from pandas._libs import writers as libwriters

    calls = []
    real = libwriters.write_csv_chunk

    def spy(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(libwriters, "write_csv_chunk", spy)
    mixed_frame.to_csv()
    assert calls

    # non-default dialects divert to the legacy writer
    calls.clear()
    mixed_frame.to_csv(quoting=csv.QUOTE_ALL)
    mixed_frame.to_csv(quoting=csv.QUOTE_NONNUMERIC)
    mixed_frame.to_csv(escapechar="\\", doublequote=False)
    mixed_frame.to_csv(sep="\N{FULLWIDTH COMMA}")
    assert not calls


def test_fast_writer_quoting_dialects_match(monkeypatch, mixed_frame):
    # diverted dialects still round-trip through the same legacy code
    assert_fast_matches_legacy(monkeypatch, mixed_frame, quoting=csv.QUOTE_ALL)
    assert_fast_matches_legacy(
        monkeypatch, mixed_frame, escapechar="\\", doublequote=False
    )
