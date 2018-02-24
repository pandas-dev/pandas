import pytest
from warnings import catch_warnings

import numpy as np
from pandas import Index, Panel, date_range, Timestamp
from .common import _check_roundtrip, ensure_clean_store, _maybe_remove
import pandas.util.testing as tm


def test_wide():
    with catch_warnings(record=True):
        wp = tm.makePanel()
        _check_roundtrip(wp, tm.assert_panel_equal)


def test_wide_table_dups():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = tm.makePanel()
            store.put('panel', wp, format='table')
            store.put('panel', wp, format='table', append=True)

            recons = store['panel']
            tm.assert_panel_equal(recons, wp)


def test_long():
    def _check(left, right):
        tm.assert_panel_equal(left.to_panel(), right.to_panel())

    with catch_warnings(record=True):
        wp = tm.makePanel()
        _check_roundtrip(wp.to_frame(), _check)


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
            tm.assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp2')
            store.put('wp2', wp, format='t')
            n = store.remove('wp2', start=-32)
            assert n == 32
            result = store.select('wp2')
            expected = wp.reindex(major_axis=wp.major_axis[:-32 // 4])
            tm.assert_panel_equal(result, expected)

            # stop
            _maybe_remove(store, 'wp3')
            store.put('wp3', wp, format='t')
            n = store.remove('wp3', stop=32)
            assert n == 32
            result = store.select('wp3')
            expected = wp.reindex(major_axis=wp.major_axis[32 // 4:])
            tm.assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp4')
            store.put('wp4', wp, format='t')
            n = store.remove('wp4', stop=-32)
            assert n == 120 - 32
            result = store.select('wp4')
            expected = wp.reindex(major_axis=wp.major_axis[-32 // 4:])
            tm.assert_panel_equal(result, expected)

            # start n stop
            _maybe_remove(store, 'wp5')
            store.put('wp5', wp, format='t')
            n = store.remove('wp5', start=16, stop=-16)
            assert n == 120 - 32
            result = store.select('wp5')
            expected = wp.reindex(
                major_axis=(wp.major_axis[:16 // 4]
                            .union(wp.major_axis[-16 // 4:])))
            tm.assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp6')
            store.put('wp6', wp, format='t')
            n = store.remove('wp6', start=16, stop=16)
            assert n == 0
            result = store.select('wp6')
            expected = wp.reindex(major_axis=wp.major_axis)
            tm.assert_panel_equal(result, expected)

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
            tm.assert_panel_equal(result, expected)


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
            tm.assert_panel_equal(result, expected)

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
            tm.assert_panel_equal(result, expected)

            # individual row elements
            _maybe_remove(store, 'wp2')
            store.put('wp2', wp, format='table')

            date1 = wp.major_axis[1:3]
            crit1 = 'major_axis=date1'
            store.remove('wp2', where=[crit1])
            result = store.select('wp2')
            expected = wp.reindex(
                major_axis=wp.major_axis.difference(date1))
            tm.assert_panel_equal(result, expected)

            date2 = wp.major_axis[5]
            crit2 = 'major_axis=date2'
            store.remove('wp2', where=[crit2])
            result = store['wp2']
            expected = wp.reindex(
                major_axis=(wp.major_axis
                            .difference(date1)
                            .difference(Index([date2]))
                            ))
            tm.assert_panel_equal(result, expected)

            date3 = [wp.major_axis[7], wp.major_axis[9]]
            crit3 = 'major_axis=date3'
            store.remove('wp2', where=[crit3])
            result = store['wp2']
            expected = wp.reindex(major_axis=wp.major_axis
                                  .difference(date1)
                                  .difference(Index([date2]))
                                  .difference(Index(date3)))
            tm.assert_panel_equal(result, expected)

            # corners
            _maybe_remove(store, 'wp4')
            store.put('wp4', wp, format='table')
            n = store.remove(
                'wp4', where="major_axis>wp.major_axis[-1]")
            result = store.select('wp4')
            tm.assert_panel_equal(result, wp)


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
            tm.assert_panel_equal(rs, expected)

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


def test_terms():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = tm.makePanel()
            wpneg = Panel.fromDict({-1: tm.makeDataFrame(),
                                    0: tm.makeDataFrame(),
                                    1: tm.makeDataFrame()})

            store.put('wp', wp, format='table')
            store.put('wpneg', wpneg, format='table')

            # panel
            result = store.select(
                'wp',
                "major_axis<'20000108' and minor_axis=['A', 'B']")
            expected = wp.truncate(
                after='20000108').reindex(minor=['A', 'B'])
            tm.assert_panel_equal(result, expected)

            # with deprecation
            result = store.select(
                'wp', where=("major_axis<'20000108' "
                             "and minor_axis=['A', 'B']"))
            expected = wp.truncate(
                after='20000108').reindex(minor=['A', 'B'])
            tm.assert_panel_equal(result, expected)

        with catch_warnings(record=True):
            # valid terms
            terms = [('major_axis=20121114'),
                     ('major_axis>20121114'),
                     (("major_axis=['20121114', '20121114']"),),
                     ('major_axis=datetime.datetime(2012, 11, 14)'),
                     'major_axis> 20121114',
                     'major_axis >20121114',
                     'major_axis > 20121114',
                     (("minor_axis=['A', 'B']"),),
                     (("minor_axis=['A', 'B']"),),
                     ((("minor_axis==['A', 'B']"),),),
                     (("items=['ItemA', 'ItemB']"),),
                     ('items=ItemA'),
                     ]

            for t in terms:
                store.select('wp', t)

            with tm.assert_raises_regex(
                    TypeError, 'Only named functions are supported'):
                store.select(
                    'wp',
                    'major_axis == (lambda x: x)("20130101")')

        with catch_warnings(record=True):
            # check USub node parsing
            res = store.select('wpneg', 'items == -1')
            expected = Panel({-1: wpneg[-1]})
            tm.assert_panel_equal(res, expected)

            with tm.assert_raises_regex(NotImplementedError,
                                        'Unary addition '
                                        'not supported'):
                store.select('wpneg', 'items == +1')


def test_term_compat():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = Panel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                       major_axis=date_range('1/1/2000', periods=5),
                       minor_axis=['A', 'B', 'C', 'D'])
            store.append('wp', wp)

            result = store.select(
                'wp', where=("major_axis>20000102 "
                             "and minor_axis=['A', 'B']"))
            expected = wp.loc[:, wp.major_axis >
                              Timestamp('20000102'), ['A', 'B']]
            tm.assert_panel_equal(result, expected)

            store.remove('wp', 'major_axis>20000103')
            result = store.select('wp')
            expected = wp.loc[:, wp.major_axis <= Timestamp('20000103'), :]
            tm.assert_panel_equal(result, expected)

    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = Panel(np.random.randn(2, 5, 4),
                       items=['Item1', 'Item2'],
                       major_axis=date_range('1/1/2000', periods=5),
                       minor_axis=['A', 'B', 'C', 'D'])
            store.append('wp', wp)

            # stringified datetimes
            result = store.select(
                'wp', 'major_axis>datetime.datetime(2000, 1, 2)')
            expected = wp.loc[:, wp.major_axis > Timestamp('20000102')]
            tm.assert_panel_equal(result, expected)

            result = store.select(
                'wp', 'major_axis>datetime.datetime(2000, 1, 2)')
            expected = wp.loc[:, wp.major_axis > Timestamp('20000102')]
            tm.assert_panel_equal(result, expected)

            result = store.select(
                'wp',
                "major_axis=[datetime.datetime(2000, 1, 2, 0, 0), "
                "datetime.datetime(2000, 1, 3, 0, 0)]")
            expected = wp.loc[:, [Timestamp('20000102'),
                                  Timestamp('20000103')]]
            tm.assert_panel_equal(result, expected)

            result = store.select(
                'wp', "minor_axis=['A', 'B']")
            expected = wp.loc[:, :, ['A', 'B']]
            tm.assert_panel_equal(result, expected)
