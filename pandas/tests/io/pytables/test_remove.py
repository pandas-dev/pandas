import pytest
from warnings import catch_warnings

import numpy as np
from pandas import Index
from pandas.tests.io.pytables.common import ensure_clean_store, _maybe_remove
import pandas.util.testing as tm
from pandas.util.testing import assert_panel_equal


def test_remove():
    with ensure_clean_store() as store:
        ts = tm.makeTimeSeries()
        df = tm.makeDataFrame()
        store['a'] = ts
        store['b'] = df
        _maybe_remove(store, 'a')
        assert len(store) == 1
        tm.assert_frame_equal(df, store['b'])

        _maybe_remove(store, 'b')
        assert len(store) == 0

        # nonexistence
        pytest.raises(KeyError, store.remove, 'a_nonexistent_store')

        # pathing
        store['a'] = ts
        store['b/foo'] = df
        _maybe_remove(store, 'foo')
        _maybe_remove(store, 'b/foo')
        assert len(store) == 1

        store['a'] = ts
        store['b/foo'] = df
        _maybe_remove(store, 'b')
        assert len(store) == 1

        # __delitem__
        store['a'] = ts
        store['b'] = df
        del store['a']
        del store['b']
        assert len(store) == 0


def test_remove_where():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            # non-existance
            crit1 = 'index>foo'
            pytest.raises(KeyError, store.remove, 'a', [crit1])

            # try to remove non-table (with crit)
            # non-table ok (where = None)
            wp = tm.makePanel(30)
            store.put('wp', wp, format='table')
            store.remove('wp', ["minor_axis=['A', 'D']"])
            rs = store.select('wp')
            expected = wp.reindex(minor_axis=['B', 'C'])
            assert_panel_equal(rs, expected)

            # empty where
            _maybe_remove(store, 'wp')
            store.put('wp', wp, format='table')

            # deleted number (entire table)
            n = store.remove('wp', [])
            assert n == 120

            # non - empty where
            _maybe_remove(store, 'wp')
            store.put('wp', wp, format='table')
            pytest.raises(ValueError, store.remove, 'wp', ['foo'])


def test_remove_startstop():
    # GH #4835 and #6177
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = tm.makePanel(30)

            # start
            _maybe_remove(store, 'wp1')
            store.put('wp1', wp, format='t')
            n = store.remove('wp1', start=32)
            assert n == 120 - 32
            result = store.select('wp1')
            expected = wp.reindex(major_axis=wp.major_axis[:32 // 4])
            assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp2')
            store.put('wp2', wp, format='t')
            n = store.remove('wp2', start=-32)
            assert n == 32
            result = store.select('wp2')
            expected = wp.reindex(major_axis=wp.major_axis[:-32 // 4])
            assert_panel_equal(result, expected)

            # stop
            _maybe_remove(store, 'wp3')
            store.put('wp3', wp, format='t')
            n = store.remove('wp3', stop=32)
            assert n == 32
            result = store.select('wp3')
            expected = wp.reindex(major_axis=wp.major_axis[32 // 4:])
            assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp4')
            store.put('wp4', wp, format='t')
            n = store.remove('wp4', stop=-32)
            assert n == 120 - 32
            result = store.select('wp4')
            expected = wp.reindex(major_axis=wp.major_axis[-32 // 4:])
            assert_panel_equal(result, expected)

            # start n stop
            _maybe_remove(store, 'wp5')
            store.put('wp5', wp, format='t')
            n = store.remove('wp5', start=16, stop=-16)
            assert n == 120 - 32
            result = store.select('wp5')
            expected = wp.reindex(
                major_axis=(wp.major_axis[:16 // 4]
                            .union(wp.major_axis[-16 // 4:])))
            assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp6')
            store.put('wp6', wp, format='t')
            n = store.remove('wp6', start=16, stop=16)
            assert n == 0
            result = store.select('wp6')
            expected = wp.reindex(major_axis=wp.major_axis)
            assert_panel_equal(result, expected)

            # with where
            _maybe_remove(store, 'wp7')

            # TODO: unused?
            date = wp.major_axis.take(np.arange(0, 30, 3))  # noqa

            crit = 'major_axis=date'
            store.put('wp7', wp, format='t')
            n = store.remove('wp7', where=[crit], stop=80)
            assert n == 28
            result = store.select('wp7')
            expected = wp.reindex(major_axis=wp.major_axis.difference(
                wp.major_axis[np.arange(0, 20, 3)]))
            assert_panel_equal(result, expected)


def test_remove_crit():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = tm.makePanel(30)

            # group row removal
            _maybe_remove(store, 'wp3')
            date4 = wp.major_axis.take([0, 1, 2, 4, 5, 6, 8, 9, 10])
            crit4 = 'major_axis=date4'
            store.put('wp3', wp, format='t')
            n = store.remove('wp3', where=[crit4])
            assert n == 36

            result = store.select('wp3')
            expected = wp.reindex(
                major_axis=wp.major_axis.difference(date4))
            assert_panel_equal(result, expected)

            # upper half
            _maybe_remove(store, 'wp')
            store.put('wp', wp, format='table')
            date = wp.major_axis[len(wp.major_axis) // 2]

            crit1 = 'major_axis>date'
            crit2 = "minor_axis=['A', 'D']"
            n = store.remove('wp', where=[crit1])
            assert n == 56

            n = store.remove('wp', where=[crit2])
            assert n == 32

            result = store['wp']
            expected = wp.truncate(after=date).reindex(minor=['B', 'C'])
            assert_panel_equal(result, expected)

            # individual row elements
            _maybe_remove(store, 'wp2')
            store.put('wp2', wp, format='table')

            date1 = wp.major_axis[1:3]
            crit1 = 'major_axis=date1'
            store.remove('wp2', where=[crit1])
            result = store.select('wp2')
            expected = wp.reindex(
                major_axis=wp.major_axis.difference(date1))
            assert_panel_equal(result, expected)

            date2 = wp.major_axis[5]
            crit2 = 'major_axis=date2'
            store.remove('wp2', where=[crit2])
            result = store['wp2']
            expected = wp.reindex(
                major_axis=(wp.major_axis
                            .difference(date1)
                            .difference(Index([date2]))
                            ))
            assert_panel_equal(result, expected)

            date3 = [wp.major_axis[7], wp.major_axis[9]]
            crit3 = 'major_axis=date3'
            store.remove('wp2', where=[crit3])
            result = store['wp2']
            expected = wp.reindex(major_axis=wp.major_axis
                                  .difference(date1)
                                  .difference(Index([date2]))
                                  .difference(Index(date3)))
            assert_panel_equal(result, expected)

            # corners
            _maybe_remove(store, 'wp4')
            store.put('wp4', wp, format='table')
            n = store.remove(
                'wp4', where="major_axis>wp.major_axis[-1]")
            result = store.select('wp4')
            assert_panel_equal(result, wp)
