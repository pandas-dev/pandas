import zoneinfo

import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io._util import arrow_table_to_pandas

pa = pytest.importorskip("pyarrow")
pytz = pytest.importorskip("pytz")


def test_arrow_table_to_pandas_normalize_timezones():
    df = pd.DataFrame(
        {"ts": pd.date_range("2024-03-01", periods=2, tz="America/New_York")},
        index=pd.date_range("2024-01-01", periods=2, tz="US/Eastern"),
    )
    expected = df.copy()
    expected.index = expected.index._with_freq(None)

    table = pa.Table.from_pandas(df)
    result = arrow_table_to_pandas(table)

    tm.assert_frame_equal(result, expected)
    assert isinstance(result["ts"].dtype.tz, zoneinfo.ZoneInfo)
    assert isinstance(result.index.tz, zoneinfo.ZoneInfo)


def test_arrow_table_to_pandas_normalize_timezones_columns():
    df = pd.DataFrame(
        [[1, 2], [3, 4]],
        columns=pd.date_range("2024-02-01", periods=2, tz="Europe/Berlin"),
    )

    table = pa.Table.from_pandas(df)
    result = arrow_table_to_pandas(table)

    tm.assert_frame_equal(result, df)
    assert isinstance(result.columns.tz, zoneinfo.ZoneInfo)


def test_arrow_table_to_pandas_normalize_timezones_multiindex():
    df = pd.DataFrame(
        {"ts": pd.date_range("2024-03-01", periods=2, tz="America/New_York")},
    ).set_index("ts", append=True, drop=False)
    expected = df.copy()

    table = pa.Table.from_pandas(df)
    result = arrow_table_to_pandas(table)

    tm.assert_frame_equal(result, expected)
    assert isinstance(result["ts"].dtype.tz, zoneinfo.ZoneInfo)
    assert isinstance(result.index.get_level_values("ts").tz, zoneinfo.ZoneInfo)
