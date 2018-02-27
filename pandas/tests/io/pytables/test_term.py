import pytest
from warnings import catch_warnings

import numpy as np
from pandas import DataFrame, Panel, date_range, Timestamp
from pandas.tests.io.pytables.common import (ensure_clean_store,
                                             ensure_clean_path)

import pandas.util.testing as tm
from pandas.util.testing import assert_panel_equal
from pandas.io.pytables import read_hdf, Term


def test_invalid_terms():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            df = tm.makeTimeDataFrame()
            df['string'] = 'foo'
            df.loc[0:4, 'string'] = 'bar'
            wp = tm.makePanel()

            store.put('df', df, format='table')
            store.put('wp', wp, format='table')

            # some invalid terms
            pytest.raises(ValueError, store.select,
                          'wp', "minor=['A', 'B']")
            pytest.raises(ValueError, store.select,
                          'wp', ["index=['20121114']"])
            pytest.raises(ValueError, store.select, 'wp', [
                "index=['20121114', '20121114']"])
            pytest.raises(TypeError, Term)

            # more invalid
            pytest.raises(
                ValueError, store.select, 'df', 'df.index[3]')
            pytest.raises(SyntaxError, store.select, 'df', 'index>')
            pytest.raises(
                ValueError, store.select, 'wp',
                "major_axis<'20000108' & minor_axis['A', 'B']")

    # from the docs
    with ensure_clean_path() as path:
        dfq = DataFrame(np.random.randn(10, 4), columns=list(
            'ABCD'), index=date_range('20130101', periods=10))
        dfq.to_hdf(path, 'dfq', format='table', data_columns=True)

        # check ok
        read_hdf(path, 'dfq',
                 where="index>Timestamp('20130104') & columns=['A', 'B']")
        read_hdf(path, 'dfq', where="A>0 or C>0")

    # catch the invalid reference
    with ensure_clean_path() as path:
        dfq = DataFrame(np.random.randn(10, 4), columns=list(
            'ABCD'), index=date_range('20130101', periods=10))
        dfq.to_hdf(path, 'dfq', format='table')

        pytest.raises(ValueError, read_hdf, path,
                      'dfq', where="A>0 or C>0")


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
            assert_panel_equal(result, expected)

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
            assert_panel_equal(result, expected)

            store.remove('wp', 'major_axis>20000103')
            result = store.select('wp')
            expected = wp.loc[:, wp.major_axis <= Timestamp('20000103'), :]
            assert_panel_equal(result, expected)

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
            assert_panel_equal(result, expected)

            result = store.select(
                'wp', 'major_axis>datetime.datetime(2000, 1, 2)')
            expected = wp.loc[:, wp.major_axis > Timestamp('20000102')]
            assert_panel_equal(result, expected)

            result = store.select(
                'wp',
                "major_axis=[datetime.datetime(2000, 1, 2, 0, 0), "
                "datetime.datetime(2000, 1, 3, 0, 0)]")
            expected = wp.loc[:, [Timestamp('20000102'),
                                  Timestamp('20000103')]]
            assert_panel_equal(result, expected)

            result = store.select(
                'wp', "minor_axis=['A', 'B']")
            expected = wp.loc[:, :, ['A', 'B']]
            assert_panel_equal(result, expected)
