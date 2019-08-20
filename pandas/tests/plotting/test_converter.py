from datetime import date, datetime
import subprocess
import sys

import numpy as np
import pytest

import pandas._config.config as cf

from pandas.compat.numpy import np_datetime64_compat

from pandas import Index, Period, Series, Timestamp, date_range
import pandas.util.testing as tm

from pandas.plotting import (
    deregister_matplotlib_converters,
    register_matplotlib_converters,
)
from pandas.tseries.offsets import Day, Micro, Milli, Second

try:
    from pandas.plotting._matplotlib import converter
except ImportError:
    # try / except, rather than skip, to avoid internal refactoring
    # causing an improprer skip
    pass

pytest.importorskip("matplotlib.pyplot")


def test_initial_warning():
    code = (
        "import pandas as pd; import matplotlib.pyplot as plt; "
        "s = pd.Series(1, pd.date_range('2000', periods=12)); "
        "fig, ax = plt.subplots(); "
        "ax.plot(s.index, s.values)"
    )
    call = [sys.executable, "-c", code]
    out = subprocess.check_output(call, stderr=subprocess.STDOUT).decode()
    assert "Using an implicitly" in out


def test_registry_mpl_resets():
    # Check that Matplotlib converters are properly reset (see issue #27481)
    code = (
        "import matplotlib.units as units; "
        "import matplotlib.dates as mdates; "
        "n_conv = len(units.registry); "
        "import pandas as pd; "
        "pd.plotting.register_matplotlib_converters(); "
        "pd.plotting.deregister_matplotlib_converters(); "
        "assert len(units.registry) == n_conv"
    )
    call = [sys.executable, "-c", code]
    subprocess.check_output(call)


def test_timtetonum_accepts_unicode():
    assert converter.time2num("00:01") == converter.time2num("00:01")


class TestRegistration:
    def test_register_by_default(self):
        # Run in subprocess to ensure a clean state
        code = (
            "'import matplotlib.units; "
            "import pandas as pd; "
            "units = dict(matplotlib.units.registry); "
            "assert pd.Timestamp in units)'"
        )
        call = [sys.executable, "-c", code]
        assert subprocess.check_call(call) == 0

    def test_warns(self):
        plt = pytest.importorskip("matplotlib.pyplot")
        s = Series(range(12), index=date_range("2017", periods=12))
        _, ax = plt.subplots()

        # Set to the "warning" state, in case this isn't the first test run
        converter._WARN = True
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False) as w:
            ax.plot(s.index, s.values)
            plt.close()

        assert len(w) == 1
        assert "Using an implicitly registered datetime converter" in str(w[0])

    def test_registering_no_warning(self):
        plt = pytest.importorskip("matplotlib.pyplot")
        s = Series(range(12), index=date_range("2017", periods=12))
        _, ax = plt.subplots()

        # Set to the "warn" state, in case this isn't the first test run
        converter._WARN = True
        register_matplotlib_converters()
        with tm.assert_produces_warning(None) as w:
            ax.plot(s.index, s.values)

        assert len(w) == 0

    def test_pandas_plots_register(self):
        pytest.importorskip("matplotlib.pyplot")
        s = Series(range(12), index=date_range("2017", periods=12))
        # Set to the "warn" state, in case this isn't the first test run
        converter._WARN = True
        with tm.assert_produces_warning(None) as w:
            s.plot()

        assert len(w) == 0

    def test_matplotlib_formatters(self):
        units = pytest.importorskip("matplotlib.units")
        assert Timestamp in units.registry

        ctx = cf.option_context("plotting.matplotlib.register_converters", False)
        with ctx:
            assert Timestamp not in units.registry

        assert Timestamp in units.registry

    def test_option_no_warning(self):
        pytest.importorskip("matplotlib.pyplot")
        ctx = cf.option_context("plotting.matplotlib.register_converters", False)
        plt = pytest.importorskip("matplotlib.pyplot")
        s = Series(range(12), index=date_range("2017", periods=12))
        _, ax = plt.subplots()

        converter._WARN = True
        # Test without registering first, no warning
        with ctx:
            with tm.assert_produces_warning(None) as w:
                ax.plot(s.index, s.values)

        assert len(w) == 0

        # Now test with registering
        converter._WARN = True
        register_matplotlib_converters()
        with ctx:
            with tm.assert_produces_warning(None) as w:
                ax.plot(s.index, s.values)

        assert len(w) == 0

    def test_registry_resets(self):
        units = pytest.importorskip("matplotlib.units")
        dates = pytest.importorskip("matplotlib.dates")

        # make a copy, to reset to
        original = dict(units.registry)

        try:
            # get to a known state
            units.registry.clear()
            date_converter = dates.DateConverter()
            units.registry[datetime] = date_converter
            units.registry[date] = date_converter

            register_matplotlib_converters()
            assert units.registry[date] is not date_converter
            deregister_matplotlib_converters()
            assert units.registry[date] is date_converter

        finally:
            # restore original stater
            units.registry.clear()
            for k, v in original.items():
                units.registry[k] = v

    def test_old_import_warns(self):
        with tm.assert_produces_warning(FutureWarning) as w:
            from pandas.tseries import converter

            converter.register()

        assert len(w)
        assert "pandas.plotting.register_matplotlib_converters" in str(w[0].message)


class TestDateTimeConverter:
    def setup_method(self, method):
        self.dtc = converter.DatetimeConverter()
        self.tc = converter.TimeFormatter(None)

    def test_convert_accepts_unicode(self):
        r1 = self.dtc.convert("12:22", None, None)
        r2 = self.dtc.convert("12:22", None, None)
        assert r1 == r2, "DatetimeConverter.convert should accept unicode"

    def test_conversion(self):
        rs = self.dtc.convert(["2012-1-1"], None, None)[0]
        xp = datetime(2012, 1, 1).toordinal()
        assert rs == xp

        rs = self.dtc.convert("2012-1-1", None, None)
        assert rs == xp

        rs = self.dtc.convert(date(2012, 1, 1), None, None)
        assert rs == xp

        rs = self.dtc.convert(datetime(2012, 1, 1).toordinal(), None, None)
        assert rs == xp

        rs = self.dtc.convert("2012-1-1", None, None)
        assert rs == xp

        rs = self.dtc.convert(Timestamp("2012-1-1"), None, None)
        assert rs == xp

        # also testing datetime64 dtype (GH8614)
        rs = self.dtc.convert(np_datetime64_compat("2012-01-01"), None, None)
        assert rs == xp

        rs = self.dtc.convert(
            np_datetime64_compat("2012-01-01 00:00:00+0000"), None, None
        )
        assert rs == xp

        rs = self.dtc.convert(
            np.array(
                [
                    np_datetime64_compat("2012-01-01 00:00:00+0000"),
                    np_datetime64_compat("2012-01-02 00:00:00+0000"),
                ]
            ),
            None,
            None,
        )
        assert rs[0] == xp

        # we have a tz-aware date (constructed to that when we turn to utc it
        # is the same as our sample)
        ts = Timestamp("2012-01-01").tz_localize("UTC").tz_convert("US/Eastern")
        rs = self.dtc.convert(ts, None, None)
        assert rs == xp

        rs = self.dtc.convert(ts.to_pydatetime(), None, None)
        assert rs == xp

        rs = self.dtc.convert(Index([ts - Day(1), ts]), None, None)
        assert rs[1] == xp

        rs = self.dtc.convert(Index([ts - Day(1), ts]).to_pydatetime(), None, None)
        assert rs[1] == xp

    def test_conversion_float(self):
        decimals = 9

        rs = self.dtc.convert(Timestamp("2012-1-1 01:02:03", tz="UTC"), None, None)
        xp = converter.dates.date2num(Timestamp("2012-1-1 01:02:03", tz="UTC"))
        tm.assert_almost_equal(rs, xp, decimals)

        rs = self.dtc.convert(
            Timestamp("2012-1-1 09:02:03", tz="Asia/Hong_Kong"), None, None
        )
        tm.assert_almost_equal(rs, xp, decimals)

        rs = self.dtc.convert(datetime(2012, 1, 1, 1, 2, 3), None, None)
        tm.assert_almost_equal(rs, xp, decimals)

    def test_conversion_outofbounds_datetime(self):
        # 2579
        values = [date(1677, 1, 1), date(1677, 1, 2)]
        rs = self.dtc.convert(values, None, None)
        xp = converter.dates.date2num(values)
        tm.assert_numpy_array_equal(rs, xp)
        rs = self.dtc.convert(values[0], None, None)
        xp = converter.dates.date2num(values[0])
        assert rs == xp

        values = [datetime(1677, 1, 1, 12), datetime(1677, 1, 2, 12)]
        rs = self.dtc.convert(values, None, None)
        xp = converter.dates.date2num(values)
        tm.assert_numpy_array_equal(rs, xp)
        rs = self.dtc.convert(values[0], None, None)
        xp = converter.dates.date2num(values[0])
        assert rs == xp

    @pytest.mark.parametrize(
        "time,format_expected",
        [
            (0, "00:00"),  # time2num(datetime.time.min)
            (86399.999999, "23:59:59.999999"),  # time2num(datetime.time.max)
            (90000, "01:00"),
            (3723, "01:02:03"),
            (39723.2, "11:02:03.200"),
        ],
    )
    def test_time_formatter(self, time, format_expected):
        # issue 18478
        result = self.tc(time)
        assert result == format_expected

    def test_dateindex_conversion(self):
        decimals = 9

        for freq in ("B", "L", "S"):
            dateindex = tm.makeDateIndex(k=10, freq=freq)
            rs = self.dtc.convert(dateindex, None, None)
            xp = converter.dates.date2num(dateindex._mpl_repr())
            tm.assert_almost_equal(rs, xp, decimals)

    def test_resolution(self):
        def _assert_less(ts1, ts2):
            val1 = self.dtc.convert(ts1, None, None)
            val2 = self.dtc.convert(ts2, None, None)
            if not val1 < val2:
                raise AssertionError("{0} is not less than {1}.".format(val1, val2))

        # Matplotlib's time representation using floats cannot distinguish
        # intervals smaller than ~10 microsecond in the common range of years.
        ts = Timestamp("2012-1-1")
        _assert_less(ts, ts + Second())
        _assert_less(ts, ts + Milli())
        _assert_less(ts, ts + Micro(50))

    def test_convert_nested(self):
        inner = [Timestamp("2017-01-01"), Timestamp("2017-01-02")]
        data = [inner, inner]
        result = self.dtc.convert(data, None, None)
        expected = [self.dtc.convert(x, None, None) for x in data]
        assert (np.array(result) == expected).all()


class TestPeriodConverter:
    def setup_method(self, method):
        self.pc = converter.PeriodConverter()

        class Axis:
            pass

        self.axis = Axis()
        self.axis.freq = "D"

    def test_convert_accepts_unicode(self):
        r1 = self.pc.convert("2012-1-1", None, self.axis)
        r2 = self.pc.convert("2012-1-1", None, self.axis)
        assert r1 == r2

    def test_conversion(self):
        rs = self.pc.convert(["2012-1-1"], None, self.axis)[0]
        xp = Period("2012-1-1").ordinal
        assert rs == xp

        rs = self.pc.convert("2012-1-1", None, self.axis)
        assert rs == xp

        rs = self.pc.convert([date(2012, 1, 1)], None, self.axis)[0]
        assert rs == xp

        rs = self.pc.convert(date(2012, 1, 1), None, self.axis)
        assert rs == xp

        rs = self.pc.convert([Timestamp("2012-1-1")], None, self.axis)[0]
        assert rs == xp

        rs = self.pc.convert(Timestamp("2012-1-1"), None, self.axis)
        assert rs == xp

        rs = self.pc.convert(np_datetime64_compat("2012-01-01"), None, self.axis)
        assert rs == xp

        rs = self.pc.convert(
            np_datetime64_compat("2012-01-01 00:00:00+0000"), None, self.axis
        )
        assert rs == xp

        rs = self.pc.convert(
            np.array(
                [
                    np_datetime64_compat("2012-01-01 00:00:00+0000"),
                    np_datetime64_compat("2012-01-02 00:00:00+0000"),
                ]
            ),
            None,
            self.axis,
        )
        assert rs[0] == xp

    def test_integer_passthrough(self):
        # GH9012
        rs = self.pc.convert([0, 1], None, self.axis)
        xp = [0, 1]
        assert rs == xp

    def test_convert_nested(self):
        data = ["2012-1-1", "2012-1-2"]
        r1 = self.pc.convert([data, data], None, self.axis)
        r2 = [self.pc.convert(data, None, self.axis) for _ in range(2)]
        assert r1 == r2
