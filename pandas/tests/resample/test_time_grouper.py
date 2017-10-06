# pylint: disable=E1101

from datetime import datetime
from warnings import catch_warnings

import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import DataFrame, Panel, Series
from pandas.compat import zip
from pandas.core.indexes.datetimes import date_range
from pandas.core.resample import TimeGrouper
from pandas.util.testing import assert_frame_equal, assert_series_equal


class TestTimeGrouper(object):

    def setup_method(self, method):
        self.ts = Series(np.random.randn(1000),
                         index=date_range('1/1/2000', periods=1000))

    def test_apply(self):
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            grouper = pd.TimeGrouper(freq='A', label='right', closed='right')

        grouped = self.ts.groupby(grouper)

        f = lambda x: x.sort_values()[-3:]

        applied = grouped.apply(f)
        expected = self.ts.groupby(lambda x: x.year).apply(f)

        applied.index = applied.index.droplevel(0)
        expected.index = expected.index.droplevel(0)
        assert_series_equal(applied, expected)

    def test_count(self):
        self.ts[::3] = np.nan

        expected = self.ts.groupby(lambda x: x.year).count()

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            grouper = pd.TimeGrouper(freq='A', label='right', closed='right')
        result = self.ts.groupby(grouper).count()
        expected.index = result.index
        assert_series_equal(result, expected)

        result = self.ts.resample('A').count()
        expected.index = result.index
        assert_series_equal(result, expected)

    def test_numpy_reduction(self):
        result = self.ts.resample('A', closed='right').prod()

        expected = self.ts.groupby(lambda x: x.year).agg(np.prod)
        expected.index = result.index

        assert_series_equal(result, expected)

    def test_apply_iteration(self):
        # #2300
        N = 1000
        ind = pd.date_range(start="2000-01-01", freq="D", periods=N)
        df = DataFrame({'open': 1, 'close': 2}, index=ind)
        tg = TimeGrouper('M')

        _, grouper, _ = tg._get_grouper(df)

        # Errors
        grouped = df.groupby(grouper, group_keys=False)
        f = lambda df: df['close'] / df['open']

        # it works!
        result = grouped.apply(f)
        tm.assert_index_equal(result.index, df.index)

    def test_panel_aggregation(self):
        ind = pd.date_range('1/1/2000', periods=100)
        data = np.random.randn(2, len(ind), 4)

        with catch_warnings(record=True):
            wp = Panel(data, items=['Item1', 'Item2'], major_axis=ind,
                       minor_axis=['A', 'B', 'C', 'D'])

            tg = TimeGrouper('M', axis=1)
            _, grouper, _ = tg._get_grouper(wp)
            bingrouped = wp.groupby(grouper)
            binagg = bingrouped.mean()

            def f(x):
                assert (isinstance(x, Panel))
                return x.mean(1)

            result = bingrouped.agg(f)
            tm.assert_panel_equal(result, binagg)

    def test_fails_on_no_datetime_index(self):
        index_names = ('Int64Index', 'Index', 'Float64Index', 'MultiIndex')
        index_funcs = (tm.makeIntIndex,
                       tm.makeUnicodeIndex, tm.makeFloatIndex,
                       lambda m: tm.makeCustomIndex(m, 2))
        n = 2
        for name, func in zip(index_names, index_funcs):
            index = func(n)
            df = DataFrame({'a': np.random.randn(n)}, index=index)
            with tm.assert_raises_regex(TypeError,
                                        "Only valid with "
                                        "DatetimeIndex, TimedeltaIndex "
                                        "or PeriodIndex, but got an "
                                        "instance of %r" % name):
                df.groupby(TimeGrouper('D'))

    def test_aaa_group_order(self):
        # GH 12840
        # check TimeGrouper perform stable sorts
        n = 20
        data = np.random.randn(n, 4)
        df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        df['key'] = [datetime(2013, 1, 1), datetime(2013, 1, 2),
                     datetime(2013, 1, 3), datetime(2013, 1, 4),
                     datetime(2013, 1, 5)] * 4
        grouped = df.groupby(TimeGrouper(key='key', freq='D'))

        tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 1)),
                              df[::5])
        tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 2)),
                              df[1::5])
        tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 3)),
                              df[2::5])
        tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 4)),
                              df[3::5])
        tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 5)),
                              df[4::5])

    def test_aggregate_normal(self):
        # check TimeGrouper's aggregation is identical as normal groupby

        n = 20
        data = np.random.randn(n, 4)
        normal_df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        normal_df['key'] = [1, 2, 3, 4, 5] * 4

        dt_df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        dt_df['key'] = [datetime(2013, 1, 1), datetime(2013, 1, 2),
                        datetime(2013, 1, 3), datetime(2013, 1, 4),
                        datetime(2013, 1, 5)] * 4

        normal_grouped = normal_df.groupby('key')
        dt_grouped = dt_df.groupby(TimeGrouper(key='key', freq='D'))

        for func in ['min', 'max', 'prod', 'var', 'std', 'mean']:
            expected = getattr(normal_grouped, func)()
            dt_result = getattr(dt_grouped, func)()
            expected.index = date_range(start='2013-01-01', freq='D',
                                        periods=5, name='key')
            assert_frame_equal(expected, dt_result)

        for func in ['count', 'sum']:
            expected = getattr(normal_grouped, func)()
            expected.index = date_range(start='2013-01-01', freq='D',
                                        periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_frame_equal(expected, dt_result)

        # GH 7453
        for func in ['size']:
            expected = getattr(normal_grouped, func)()
            expected.index = date_range(start='2013-01-01', freq='D',
                                        periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_series_equal(expected, dt_result)

        # GH 7453
        for func in ['first', 'last']:
            expected = getattr(normal_grouped, func)()
            expected.index = date_range(start='2013-01-01', freq='D',
                                        periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_frame_equal(expected, dt_result)

        # if TimeGrouper is used included, 'nth' doesn't work yet

        """
        for func in ['nth']:
            expected = getattr(normal_grouped, func)(3)
            expected.index = date_range(start='2013-01-01',
                                        freq='D', periods=5, name='key')
            dt_result = getattr(dt_grouped, func)(3)
            assert_frame_equal(expected, dt_result)
        """

    def test_aggregate_with_nat(self):
        # check TimeGrouper's aggregation is identical as normal groupby

        n = 20
        data = np.random.randn(n, 4).astype('int64')
        normal_df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        normal_df['key'] = [1, 2, np.nan, 4, 5] * 4

        dt_df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        dt_df['key'] = [datetime(2013, 1, 1), datetime(2013, 1, 2), pd.NaT,
                        datetime(2013, 1, 4), datetime(2013, 1, 5)] * 4

        normal_grouped = normal_df.groupby('key')
        dt_grouped = dt_df.groupby(TimeGrouper(key='key', freq='D'))

        for func in ['min', 'max', 'sum', 'prod']:
            normal_result = getattr(normal_grouped, func)()
            dt_result = getattr(dt_grouped, func)()
            pad = DataFrame([[np.nan, np.nan, np.nan, np.nan]], index=[3],
                            columns=['A', 'B', 'C', 'D'])
            expected = normal_result.append(pad)
            expected = expected.sort_index()
            expected.index = date_range(start='2013-01-01', freq='D',
                                        periods=5, name='key')
            assert_frame_equal(expected, dt_result)

        for func in ['count']:
            normal_result = getattr(normal_grouped, func)()
            pad = DataFrame([[0, 0, 0, 0]], index=[3],
                            columns=['A', 'B', 'C', 'D'])
            expected = normal_result.append(pad)
            expected = expected.sort_index()
            expected.index = date_range(start='2013-01-01', freq='D',
                                        periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_frame_equal(expected, dt_result)

        for func in ['size']:
            normal_result = getattr(normal_grouped, func)()
            pad = Series([0], index=[3])
            expected = normal_result.append(pad)
            expected = expected.sort_index()
            expected.index = date_range(start='2013-01-01', freq='D',
                                        periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_series_equal(expected, dt_result)
            # GH 9925
            assert dt_result.index.name == 'key'

            # if NaT is included, 'var', 'std', 'mean', 'first','last'
            # and 'nth' doesn't work yet
