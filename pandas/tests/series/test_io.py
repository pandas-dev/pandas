import collections
from datetime import datetime
from io import StringIO

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series
import pandas.util.testing as tm
from pandas.util.testing import (
    assert_almost_equal,
    assert_frame_equal,
    assert_series_equal,
    ensure_clean,
)

from pandas.io.common import _get_handle


class TestSeriesToCSV:
    def read_csv(self, path, **kwargs):
        params = dict(squeeze=True, index_col=0, header=None, parse_dates=True)
        params.update(**kwargs)

        header = params.get("header")
        out = pd.read_csv(path, **params)

        if header is None:
            out.name = out.index.name = None

        return out

    @pytest.mark.parametrize("arg", ["path", "header", "both"])
    def test_to_csv_deprecation(self, arg, datetime_series):
        # see gh-19715
        with ensure_clean() as path:
            if arg == "path":
                kwargs = dict(path=path, header=False)
            elif arg == "header":
                kwargs = dict(path_or_buf=path)
            else:  # Both discrepancies match.
                kwargs = dict(path=path)

            with tm.assert_produces_warning(FutureWarning):
                datetime_series.to_csv(**kwargs)

                # Make sure roundtrip still works.
                ts = self.read_csv(path)
                assert_series_equal(datetime_series, ts, check_names=False)

    def test_from_csv(self, datetime_series, string_series):

        with ensure_clean() as path:
            datetime_series.to_csv(path, header=False)
            ts = self.read_csv(path)
            assert_series_equal(datetime_series, ts, check_names=False)

            assert ts.name is None
            assert ts.index.name is None

            # see gh-10483
            datetime_series.to_csv(path, header=True)
            ts_h = self.read_csv(path, header=0)
            assert ts_h.name == "ts"

            string_series.to_csv(path, header=False)
            series = self.read_csv(path)
            assert_series_equal(string_series, series, check_names=False)

            assert series.name is None
            assert series.index.name is None

            string_series.to_csv(path, header=True)
            series_h = self.read_csv(path, header=0)
            assert series_h.name == "series"

            with open(path, "w") as outfile:
                outfile.write("1998-01-01|1.0\n1999-01-01|2.0")

            series = self.read_csv(path, sep="|")
            check_series = Series(
                {datetime(1998, 1, 1): 1.0, datetime(1999, 1, 1): 2.0}
            )
            assert_series_equal(check_series, series)

            series = self.read_csv(path, sep="|", parse_dates=False)
            check_series = Series({"1998-01-01": 1.0, "1999-01-01": 2.0})
            assert_series_equal(check_series, series)

    def test_to_csv(self, datetime_series):
        import io

        with ensure_clean() as path:
            datetime_series.to_csv(path, header=False)

            with io.open(path, newline=None) as f:
                lines = f.readlines()
            assert lines[1] != "\n"

            datetime_series.to_csv(path, index=False, header=False)
            arr = np.loadtxt(path)
            assert_almost_equal(arr, datetime_series.values)

    def test_to_csv_unicode_index(self):
        buf = StringIO()
        s = Series(["\u05d0", "d2"], index=["\u05d0", "\u05d1"])

        s.to_csv(buf, encoding="UTF-8", header=False)
        buf.seek(0)

        s2 = self.read_csv(buf, index_col=0, encoding="UTF-8")
        assert_series_equal(s, s2)

    def test_to_csv_float_format(self):

        with ensure_clean() as filename:
            ser = Series([0.123456, 0.234567, 0.567567])
            ser.to_csv(filename, float_format="%.2f", header=False)

            rs = self.read_csv(filename)
            xp = Series([0.12, 0.23, 0.57])
            assert_series_equal(rs, xp)

    def test_to_csv_list_entries(self):
        s = Series(["jack and jill", "jesse and frank"])

        split = s.str.split(r"\s+and\s+")

        buf = StringIO()
        split.to_csv(buf, header=False)

    def test_to_csv_path_is_none(self):
        # GH 8215
        # Series.to_csv() was returning None, inconsistent with
        # DataFrame.to_csv() which returned string
        s = Series([1, 2, 3])
        csv_str = s.to_csv(path_or_buf=None, header=False)
        assert isinstance(csv_str, str)

    @pytest.mark.parametrize(
        "s,encoding",
        [
            (
                Series([0.123456, 0.234567, 0.567567], index=["A", "B", "C"], name="X"),
                None,
            ),
            # GH 21241, 21118
            (Series(["abc", "def", "ghi"], name="X"), "ascii"),
            (Series(["123", "你好", "世界"], name="中文"), "gb2312"),
            (Series(["123", "Γειά σου", "Κόσμε"], name="Ελληνικά"), "cp737"),
        ],
    )
    def test_to_csv_compression(self, s, encoding, compression):

        with ensure_clean() as filename:

            s.to_csv(filename, compression=compression, encoding=encoding, header=True)
            # test the round trip - to_csv -> read_csv
            result = pd.read_csv(
                filename,
                compression=compression,
                encoding=encoding,
                index_col=0,
                squeeze=True,
            )
            assert_series_equal(s, result)

            # test the round trip using file handle - to_csv -> read_csv
            f, _handles = _get_handle(
                filename, "w", compression=compression, encoding=encoding
            )
            with f:
                s.to_csv(f, encoding=encoding, header=True)
            result = pd.read_csv(
                filename,
                compression=compression,
                encoding=encoding,
                index_col=0,
                squeeze=True,
            )
            assert_series_equal(s, result)

            # explicitly ensure file was compressed
            with tm.decompress_file(filename, compression) as fh:
                text = fh.read().decode(encoding or "utf8")
                assert s.name in text

            with tm.decompress_file(filename, compression) as fh:
                assert_series_equal(
                    s, pd.read_csv(fh, index_col=0, squeeze=True, encoding=encoding)
                )

    def test_to_csv_interval_index(self):
        # GH 28210
        s = Series(["foo", "bar", "baz"], index=pd.interval_range(0, 3))

        with ensure_clean("__tmp_to_csv_interval_index__.csv") as path:
            s.to_csv(path, header=False)
            result = self.read_csv(path, index_col=0, squeeze=True)

            # can't roundtrip intervalindex via read_csv so check string repr (GH 23595)
            expected = s.copy()
            expected.index = expected.index.astype(str)

            assert_series_equal(result, expected)


class TestSeriesIO:
    def test_to_frame(self, datetime_series):
        datetime_series.name = None
        rs = datetime_series.to_frame()
        xp = pd.DataFrame(datetime_series.values, index=datetime_series.index)
        assert_frame_equal(rs, xp)

        datetime_series.name = "testname"
        rs = datetime_series.to_frame()
        xp = pd.DataFrame(
            dict(testname=datetime_series.values), index=datetime_series.index
        )
        assert_frame_equal(rs, xp)

        rs = datetime_series.to_frame(name="testdifferent")
        xp = pd.DataFrame(
            dict(testdifferent=datetime_series.values), index=datetime_series.index
        )
        assert_frame_equal(rs, xp)

    def test_timeseries_periodindex(self):
        # GH2891
        from pandas import period_range

        prng = period_range("1/1/2011", "1/1/2012", freq="M")
        ts = Series(np.random.randn(len(prng)), prng)
        new_ts = tm.round_trip_pickle(ts)
        assert new_ts.index.freq == "M"

    def test_pickle_preserve_name(self):
        for n in [777, 777.0, "name", datetime(2001, 11, 11), (1, 2)]:
            unpickled = self._pickle_roundtrip_name(tm.makeTimeSeries(name=n))
            assert unpickled.name == n

    def test_pickle_categorical_ordered_from_sentinel(self):
        # GH 27295: can remove test when _ordered_from_sentinel is removed (GH 26403)
        s = Series(["a", "b", "c", "a"], dtype="category")
        result = tm.round_trip_pickle(s)
        result = result.astype("category")

        tm.assert_series_equal(result, s)
        assert result.dtype._ordered_from_sentinel is False

    def _pickle_roundtrip_name(self, obj):

        with ensure_clean() as path:
            obj.to_pickle(path)
            unpickled = pd.read_pickle(path)
            return unpickled

    def test_to_frame_expanddim(self):
        # GH 9762

        class SubclassedSeries(Series):
            @property
            def _constructor_expanddim(self):
                return SubclassedFrame

        class SubclassedFrame(DataFrame):
            pass

        s = SubclassedSeries([1, 2, 3], name="X")
        result = s.to_frame()
        assert isinstance(result, SubclassedFrame)
        expected = SubclassedFrame({"X": [1, 2, 3]})
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "mapping", (dict, collections.defaultdict(list), collections.OrderedDict)
    )
    def test_to_dict(self, mapping, datetime_series):
        # GH16122
        tm.assert_series_equal(
            Series(datetime_series.to_dict(mapping), name="ts"), datetime_series
        )
        from_method = Series(datetime_series.to_dict(collections.Counter))
        from_constructor = Series(collections.Counter(datetime_series.items()))
        tm.assert_series_equal(from_method, from_constructor)
