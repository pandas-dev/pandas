# pylint: disable=W0612,E1101


from datetime import datetime
import os
import operator
import unittest
import nose

import numpy as np

from pandas import DataFrame, Index, isnull, notnull, pivot, MultiIndex
from pandas.core.datetools import bday
from pandas.core.frame import group_agg
from pandas.core.panel import Panel
from pandas.core.series import remove_na
import pandas.core.common as com
import pandas.core.panel as panelmod
from pandas.util import py3compat
from pandas.io.parsers import (ExcelFile, ExcelWriter)

from pandas.util.testing import (assert_panel_equal,
                                 assert_frame_equal,
                                 assert_series_equal,
                                 assert_almost_equal)
import pandas.core.panel as panelm
import pandas.util.testing as tm

class PanelTests(object):
    panel = None

    def test_pickle(self):
        import cPickle
        pickled = cPickle.dumps(self.panel)
        unpickled = cPickle.loads(pickled)
        assert_frame_equal(unpickled['ItemA'], self.panel['ItemA'])

    def test_cumsum(self):
        cumsum = self.panel.cumsum()
        assert_frame_equal(cumsum['ItemA'], self.panel['ItemA'].cumsum())

class SafeForLongAndSparse(object):

    def test_repr(self):
        foo = repr(self.panel)

    def test_iter(self):
        tm.equalContents(list(self.panel), self.panel.items)

    def test_count(self):
        f = lambda s: notnull(s).sum()
        self._check_stat_op('count', f, obj=self.panel, has_skipna=False)

    def test_sum(self):
        self._check_stat_op('sum', np.sum)

    def test_mean(self):
        self._check_stat_op('mean', np.mean)

    def test_prod(self):
        self._check_stat_op('prod', np.prod)

    def test_median(self):
        def wrapper(x):
            if isnull(x).any():
                return np.nan
            return np.median(x)

        self._check_stat_op('median', wrapper)

    def test_min(self):
        self._check_stat_op('min', np.min)

    def test_max(self):
        self._check_stat_op('max', np.max)

    def test_skew(self):
        from scipy.stats import skew
        def this_skew(x):
            if len(x) < 3:
                return np.nan
            return skew(x, bias=False)
        self._check_stat_op('skew', this_skew)

    # def test_mad(self):
    #     f = lambda x: np.abs(x - x.mean()).mean()
    #     self._check_stat_op('mad', f)

    def test_var(self):
        def alt(x):
            if len(x) < 2:
                return np.nan
            return np.var(x, ddof=1)
        self._check_stat_op('var', alt)

    def test_std(self):
        def alt(x):
            if len(x) < 2:
                return np.nan
            return np.std(x, ddof=1)
        self._check_stat_op('std', alt)

    # def test_skew(self):
    #     from scipy.stats import skew

    #     def alt(x):
    #         if len(x) < 3:
    #             return np.nan
    #         return skew(x, bias=False)

    #     self._check_stat_op('skew', alt)

    def _check_stat_op(self, name, alternative, obj=None, has_skipna=True):
        if obj is None:
            obj = self.panel

            # # set some NAs
            # obj.ix[5:10] = np.nan
            # obj.ix[15:20, -2:] = np.nan

        f = getattr(obj, name)

        if has_skipna:
            def skipna_wrapper(x):
                nona = remove_na(x)
                if len(nona) == 0:
                    return np.nan
                return alternative(nona)

            def wrapper(x):
                return alternative(np.asarray(x))

            for i in range(obj.ndim):
                result = f(axis=i, skipna=False)
                assert_frame_equal(result, obj.apply(wrapper, axis=i))
        else:
            skipna_wrapper = alternative
            wrapper = alternative

        for i in range(obj.ndim):
            result = f(axis=i)
            assert_frame_equal(result, obj.apply(skipna_wrapper, axis=i))

        self.assertRaises(Exception, f, axis=obj.ndim)

class SafeForSparse(object):

    @classmethod
    def assert_panel_equal(cls, x, y):
        assert_panel_equal(x, y)

    def test_get_axis(self):
        assert(self.panel._get_axis(0) is self.panel.items)
        assert(self.panel._get_axis(1) is self.panel.major_axis)
        assert(self.panel._get_axis(2) is self.panel.minor_axis)

    def test_set_axis(self):
        new_items = Index(np.arange(len(self.panel.items)))
        new_major = Index(np.arange(len(self.panel.major_axis)))
        new_minor = Index(np.arange(len(self.panel.minor_axis)))

        # ensure propagate to potentially prior-cached items too
        item = self.panel['ItemA']
        self.panel.items = new_items

        if hasattr(self.panel, '_item_cache'):
            self.assert_('ItemA' not in self.panel._item_cache)
        self.assert_(self.panel.items is new_items)

        item = self.panel[0]
        self.panel.major_axis = new_major
        self.assert_(self.panel[0].index is new_major)
        self.assert_(self.panel.major_axis is new_major)

        item = self.panel[0]
        self.panel.minor_axis = new_minor
        self.assert_(self.panel[0].columns is new_minor)
        self.assert_(self.panel.minor_axis is new_minor)

    def test_get_axis_number(self):
        self.assertEqual(self.panel._get_axis_number('items'), 0)
        self.assertEqual(self.panel._get_axis_number('major'), 1)
        self.assertEqual(self.panel._get_axis_number('minor'), 2)

    def test_get_axis_name(self):
        self.assertEqual(self.panel._get_axis_name(0), 'items')
        self.assertEqual(self.panel._get_axis_name(1), 'major_axis')
        self.assertEqual(self.panel._get_axis_name(2), 'minor_axis')

    def test_get_plane_axes(self):
        # what to do here?

        index, columns = self.panel._get_plane_axes('items')
        index, columns = self.panel._get_plane_axes('major_axis')
        index, columns = self.panel._get_plane_axes('minor_axis')
        index, columns = self.panel._get_plane_axes(0)

    def test_truncate(self):
        dates = self.panel.major_axis
        start, end = dates[1], dates[5]

        trunced = self.panel.truncate(start, end, axis='major')
        expected = self.panel['ItemA'].truncate(start, end)

        assert_frame_equal(trunced['ItemA'], expected)

        trunced = self.panel.truncate(before=start, axis='major')
        expected = self.panel['ItemA'].truncate(before=start)

        assert_frame_equal(trunced['ItemA'], expected)

        trunced = self.panel.truncate(after=end, axis='major')
        expected = self.panel['ItemA'].truncate(after=end)

        assert_frame_equal(trunced['ItemA'], expected)

        # XXX test other axes

    def test_arith(self):
        self._test_op(self.panel, operator.add)
        self._test_op(self.panel, operator.sub)
        self._test_op(self.panel, operator.mul)
        self._test_op(self.panel, operator.truediv)
        self._test_op(self.panel, operator.floordiv)
        self._test_op(self.panel, operator.pow)

        self._test_op(self.panel, lambda x, y: y + x)
        self._test_op(self.panel, lambda x, y: y - x)
        self._test_op(self.panel, lambda x, y: y * x)
        self._test_op(self.panel, lambda x, y: y / x)
        self._test_op(self.panel, lambda x, y: y ** x)

        self.assertRaises(Exception, self.panel.__add__, self.panel['ItemA'])

    @staticmethod
    def _test_op(panel, op):
        result = op(panel, 1)
        assert_frame_equal(result['ItemA'], op(panel['ItemA'], 1))

    def test_keys(self):
        tm.equalContents(self.panel.keys(), self.panel.items)

    def test_iteritems(self):
        """Test panel.iteritems(), aka panel.iterkv()"""
        # just test that it works
        for k, v in self.panel.iterkv():
            pass

        self.assertEqual(len(list(self.panel.iterkv())),
                         len(self.panel.items))

    def test_combineFrame(self):
        def check_op(op, name):
            # items
            df = self.panel['ItemA']

            func = getattr(self.panel, name)

            result = func(df, axis='items')

            assert_frame_equal(result['ItemB'], op(self.panel['ItemB'], df))

            # major
            xs = self.panel.major_xs(self.panel.major_axis[0])
            result = func(xs, axis='major')

            idx = self.panel.major_axis[1]

            assert_frame_equal(result.major_xs(idx),
                               op(self.panel.major_xs(idx), xs))

            # minor
            xs = self.panel.minor_xs(self.panel.minor_axis[0])
            result = func(xs, axis='minor')

            idx = self.panel.minor_axis[1]

            assert_frame_equal(result.minor_xs(idx),
                               op(self.panel.minor_xs(idx), xs))

        check_op(operator.add, 'add')
        check_op(operator.sub, 'subtract')
        check_op(operator.mul, 'multiply')
        if py3compat.PY3:
            check_op(operator.truediv, 'divide')
        else:
            check_op(operator.div, 'divide')

    def test_combinePanel(self):
        result = self.panel.add(self.panel)
        self.assert_panel_equal(result, self.panel * 2)

    def test_neg(self):
        self.assert_panel_equal(-self.panel, self.panel * -1)

    def test_select(self):
        p = self.panel

        # select items
        result = p.select(lambda x: x in ('ItemA', 'ItemC'), axis='items')
        expected = p.reindex(items=['ItemA', 'ItemC'])
        self.assert_panel_equal(result, expected)

        # select major_axis
        result = p.select(lambda x: x >= datetime(2000, 1, 15), axis='major')
        new_major = p.major_axis[p.major_axis >= datetime(2000, 1, 15)]
        expected = p.reindex(major=new_major)
        self.assert_panel_equal(result, expected)

        # select minor_axis
        result = p.select(lambda x: x in ('D', 'A'), axis=2)
        expected = p.reindex(minor=['A', 'D'])
        self.assert_panel_equal(result, expected)

        # corner case, empty thing
        result = p.select(lambda x: x in ('foo',), axis='items')
        self.assert_panel_equal(result, p.reindex(items=[]))

    def test_get_value(self):
        for item in self.panel.items:
            for mjr in self.panel.major_axis[::2]:
                for mnr in self.panel.minor_axis:
                    result = self.panel.get_value(item, mjr, mnr)
                    expected = self.panel[item][mnr][mjr]
                    assert_almost_equal(result, expected)

    def test_abs(self):
        result = self.panel.abs()
        expected = np.abs(self.panel)
        self.assert_panel_equal(result, expected)

        df = self.panel['ItemA']
        result = df.abs()
        expected = np.abs(df)
        assert_frame_equal(result, expected)

        s = df['A']
        result = s.abs()
        expected = np.abs(s)
        assert_series_equal(result, expected)

class CheckIndexing(object):


    def test_getitem(self):
        self.assertRaises(Exception, self.panel.__getitem__, 'ItemQ')

    def test_delitem_and_pop(self):
        expected = self.panel['ItemA']
        result = self.panel.pop('ItemA')
        assert_frame_equal(expected, result)
        self.assert_('ItemA' not in self.panel.items)

        del self.panel['ItemB']
        self.assert_('ItemB' not in self.panel.items)
        self.assertRaises(Exception, self.panel.__delitem__, 'ItemB')

        values = np.empty((3, 3, 3))
        values[0] = 0
        values[1] = 1
        values[2] = 2

        panel = Panel(values, range(3), range(3), range(3))

        # did we delete the right row?

        panelc = panel.copy()
        del panelc[0]
        assert_frame_equal(panelc[1], panel[1])
        assert_frame_equal(panelc[2], panel[2])

        panelc = panel.copy()
        del panelc[1]
        assert_frame_equal(panelc[0], panel[0])
        assert_frame_equal(panelc[2], panel[2])

        panelc = panel.copy()
        del panelc[2]
        assert_frame_equal(panelc[1], panel[1])
        assert_frame_equal(panelc[0], panel[0])

    def test_setitem(self):
        # LongPanel with one item
        lp = self.panel.filter(['ItemA', 'ItemB']).to_frame()
        self.assertRaises(Exception, self.panel.__setitem__,
                          'ItemE', lp)

        # DataFrame
        df = self.panel['ItemA'][2:].filter(items=['A', 'B'])
        self.panel['ItemF'] = df
        self.panel['ItemE'] = df

        df2 = self.panel['ItemF']

        assert_frame_equal(df, df2.reindex(index=df.index,
                                           columns=df.columns))

        # scalar
        self.panel['ItemG'] = 1
        self.panel['ItemE'] = True
        self.assert_(self.panel['ItemG'].values.dtype == np.int64)
        self.assert_(self.panel['ItemE'].values.dtype == np.bool_)

        # object dtype
        self.panel['ItemQ'] = 'foo'
        self.assert_(self.panel['ItemQ'].values.dtype == np.object_)

        # boolean dtype
        self.panel['ItemP'] = self.panel['ItemA'] > 0
        self.assert_(self.panel['ItemP'].values.dtype == np.bool_)

    def test_setitem_ndarray(self):
        from pandas import DateRange, datetools

        timeidx = DateRange(start=datetime(2009,1,1),
                            end=datetime(2009,12,31),
                            offset=datetools.MonthEnd())
        lons_coarse = np.linspace(-177.5, 177.5, 72)
        lats_coarse = np.linspace(-87.5, 87.5, 36)
        P = Panel(items=timeidx, major_axis=lons_coarse, minor_axis=lats_coarse)
        data = np.random.randn(72*36).reshape((72,36))
        key = datetime(2009,2,28)
        P[key] = data

        assert_almost_equal(P[key].values, data)

    def test_major_xs(self):
        ref = self.panel['ItemA']

        idx = self.panel.major_axis[5]
        xs = self.panel.major_xs(idx)

        assert_series_equal(xs['ItemA'], ref.xs(idx))

        # not contained
        idx = self.panel.major_axis[0] - bday
        self.assertRaises(Exception, self.panel.major_xs, idx)

    def test_major_xs_mixed(self):
        self.panel['ItemD'] = 'foo'
        xs = self.panel.major_xs(self.panel.major_axis[0])
        self.assert_(xs['ItemA'].dtype == np.float64)
        self.assert_(xs['ItemD'].dtype == np.object_)

    def test_minor_xs(self):
        ref = self.panel['ItemA']

        idx = self.panel.minor_axis[1]
        xs = self.panel.minor_xs(idx)

        assert_series_equal(xs['ItemA'], ref[idx])

        # not contained
        self.assertRaises(Exception, self.panel.minor_xs, 'E')

    def test_minor_xs_mixed(self):
        self.panel['ItemD'] = 'foo'

        xs = self.panel.minor_xs('D')
        self.assert_(xs['ItemA'].dtype == np.float64)
        self.assert_(xs['ItemD'].dtype == np.object_)

    def test_xs(self):
        itemA = self.panel.xs('ItemA', axis=0)
        expected = self.panel['ItemA']
        assert_frame_equal(itemA, expected)

        # not view by default
        itemA.values[:] = np.nan
        self.assert_(not np.isnan(self.panel['ItemA'].values).all())

        # but can get view
        itemA_view = self.panel.xs('ItemA', axis=0, copy=False)
        itemA_view.values[:] = np.nan
        self.assert_(np.isnan(self.panel['ItemA'].values).all())

        # mixed-type
        self.panel['strings'] = 'foo'
        self.assertRaises(Exception, self.panel.xs, 'D', axis=2,
                          copy=False)

    def test_getitem_fancy_labels(self):
        p = self.panel

        items = p.items[[1, 0]]
        dates = p.major_axis[::2]
        cols = ['D', 'C', 'F']

        # all 3 specified
        assert_panel_equal(p.ix[items, dates, cols],
                           p.reindex(items=items, major=dates, minor=cols))

        # 2 specified
        assert_panel_equal(p.ix[:, dates, cols],
                           p.reindex(major=dates, minor=cols))

        assert_panel_equal(p.ix[items, :, cols],
                           p.reindex(items=items, minor=cols))

        assert_panel_equal(p.ix[items, dates, :],
                           p.reindex(items=items, major=dates))

        # only 1
        assert_panel_equal(p.ix[items, :, :],
                           p.reindex(items=items))

        assert_panel_equal(p.ix[:, dates, :],
                           p.reindex(major=dates))

        assert_panel_equal(p.ix[:, :, cols],
                           p.reindex(minor=cols))

    def test_getitem_fancy_slice(self):
        pass

    def test_getitem_fancy_ints(self):
        pass

    def test_getitem_fancy_xs(self):
        p = self.panel
        item = 'ItemB'
        date = p.major_axis[5]
        col = 'C'

        # get DataFrame
        # item
        assert_frame_equal(p.ix[item], p[item])
        assert_frame_equal(p.ix[item, :], p[item])
        assert_frame_equal(p.ix[item, :, :], p[item])

        # major axis, axis=1
        assert_frame_equal(p.ix[:, date], p.major_xs(date))
        assert_frame_equal(p.ix[:, date, :], p.major_xs(date))

        # minor axis, axis=2
        assert_frame_equal(p.ix[:, :, 'C'], p.minor_xs('C'))

        # get Series
        assert_series_equal(p.ix[item, date], p[item].ix[date])
        assert_series_equal(p.ix[item, date, :], p[item].ix[date])
        assert_series_equal(p.ix[item, :, col], p[item][col])
        assert_series_equal(p.ix[:, date, col], p.major_xs(date).ix[col])

    def test_getitem_fancy_xs_check_view(self):
        item = 'ItemB'
        date = self.panel.major_axis[5]
        col = 'C'

        # make sure it's always a view
        NS = slice(None, None)

        # DataFrames
        comp = assert_frame_equal
        self._check_view(item, comp)
        self._check_view((item, NS), comp)
        self._check_view((item, NS, NS), comp)
        self._check_view((NS, date), comp)
        self._check_view((NS, date, NS), comp)
        self._check_view((NS, NS, 'C'), comp)

        # Series
        comp = assert_series_equal
        self._check_view((item, date), comp)
        self._check_view((item, date, NS), comp)
        self._check_view((item, NS, 'C'), comp)
        self._check_view((NS, date, 'C'), comp)

    def _check_view(self, indexer, comp):
        cp = self.panel.copy()
        obj = cp.ix[indexer]
        obj.values[:] = 0
        self.assert_((obj.values == 0).all())
        comp(cp.ix[indexer].reindex_like(obj), obj)

    def test_get_value(self):
        for item in self.panel.items:
            for mjr in self.panel.major_axis[::2]:
                for mnr in self.panel.minor_axis:
                    result = self.panel.get_value(item, mjr, mnr)
                    expected = self.panel[item][mnr][mjr]
                    assert_almost_equal(result, expected)

    def test_set_value(self):
        for item in self.panel.items:
            for mjr in self.panel.major_axis[::2]:
                for mnr in self.panel.minor_axis:
                    self.panel.set_value(item, mjr, mnr, 1.)
                    assert_almost_equal(self.panel[item][mnr][mjr], 1.)

        # resize
        res = self.panel.set_value('ItemE', 'foo', 'bar', 1.5)
        self.assert_(isinstance(res, Panel))
        self.assert_(res is not self.panel)
        self.assertEqual(res.get_value('ItemE', 'foo', 'bar'), 1.5)

        res3 = self.panel.set_value('ItemE', 'foobar', 'baz', 5)
        self.assert_(com.is_float_dtype(res3['ItemE'].values))

class TestPanel(unittest.TestCase, PanelTests, CheckIndexing,
                SafeForLongAndSparse,
                SafeForSparse):

    @classmethod
    def assert_panel_equal(cls,x, y):
        assert_panel_equal(x, y)

    def setUp(self):
        self.panel = tm.makePanel()
        tm.add_nans(self.panel)

    def test_constructor(self):
        # with BlockManager
        wp = Panel(self.panel._data)
        self.assert_(wp._data is self.panel._data)

        wp = Panel(self.panel._data, copy=True)
        self.assert_(wp._data is not self.panel._data)
        assert_panel_equal(wp, self.panel)

        # strings handled prop
        wp = Panel([[['foo', 'foo', 'foo',],
                         ['foo', 'foo', 'foo']]])
        self.assert_(wp.values.dtype == np.object_)

        vals = self.panel.values

        # no copy
        wp = Panel(vals)
        self.assert_(wp.values is vals)

        # copy
        wp = Panel(vals, copy=True)
        self.assert_(wp.values is not vals)

    def test_constructor_cast(self):
        casted = Panel(self.panel._data, dtype=int)
        casted2 = Panel(self.panel.values, dtype=int)

        exp_values = self.panel.values.astype(int)
        assert_almost_equal(casted.values, exp_values)
        assert_almost_equal(casted2.values, exp_values)

        # can't cast
        data = [[['foo', 'bar', 'baz']]]
        self.assertRaises(ValueError, Panel, data, dtype=float)

    def test_constructor_empty_panel(self):
        empty = Panel()
        self.assert_(len(empty.items) == 0)
        self.assert_(len(empty.major_axis) == 0)
        self.assert_(len(empty.minor_axis) == 0)

    def test_constructor_observe_dtype(self):
        # GH #411
        panel = Panel(items=range(3), major_axis=range(3),
                      minor_axis=range(3), dtype='O')
        self.assert_(panel.values.dtype == np.object_)

    def test_consolidate(self):
        self.assert_(self.panel._data.is_consolidated())

        self.panel['foo'] = 1.
        self.assert_(not self.panel._data.is_consolidated())

        panel = self.panel.consolidate()
        self.assert_(panel._data.is_consolidated())

    def test_ctor_dict(self):
        itema = self.panel['ItemA']
        itemb = self.panel['ItemB']

        d = {'A' : itema, 'B' : itemb[5:]}
        d2 = {'A' : itema._series, 'B' : itemb[5:]._series}
        d3 = {'A' : DataFrame(itema._series),
              'B' : DataFrame(itemb[5:]._series)}

        wp = Panel.from_dict(d)
        wp2 = Panel.from_dict(d2) # nested Dict
        wp3 = Panel.from_dict(d3)
        self.assert_(wp.major_axis.equals(self.panel.major_axis))
        assert_panel_equal(wp, wp2)

        # intersect
        wp = Panel.from_dict(d, intersect=True)
        self.assert_(wp.major_axis.equals(itemb.index[5:]))

        # use constructor
        assert_panel_equal(Panel(d), Panel.from_dict(d))
        assert_panel_equal(Panel(d2), Panel.from_dict(d2))
        assert_panel_equal(Panel(d3), Panel.from_dict(d3))

        # cast
        result = Panel(d, dtype=int)
        expected = Panel(dict((k, v.astype(int)) for k, v in d.iteritems()))

    def test_constructor_dict_mixed(self):
        data = dict((k, v.values) for k, v in self.panel.iterkv())
        result = Panel(data)
        exp_major = Index(np.arange(len(self.panel.major_axis)))
        self.assert_(result.major_axis.equals(exp_major))

        result = Panel(data, items=self.panel.items,
                       major_axis=self.panel.major_axis,
                       minor_axis=self.panel.minor_axis)
        assert_panel_equal(result, self.panel)

        data['ItemC'] = self.panel['ItemC']
        result = Panel(data)
        assert_panel_equal(result, self.panel)

        # corner, blow up
        data['ItemB'] = data['ItemB'][:-1]
        self.assertRaises(Exception, Panel, data)

        data['ItemB'] = self.panel['ItemB'].values[:, :-1]
        self.assertRaises(Exception, Panel, data)

    def test_constructor_resize(self):
        data = self.panel._data
        items = self.panel.items[:-1]
        major = self.panel.major_axis[:-1]
        minor = self.panel.minor_axis[:-1]

        result = Panel(data, items=items, major_axis=major,
                       minor_axis=minor)
        expected = self.panel.reindex(items=items, major=major, minor=minor)
        assert_panel_equal(result, expected)

        result = Panel(data, items=items, major_axis=major)
        expected = self.panel.reindex(items=items, major=major)
        assert_panel_equal(result, expected)

        result = Panel(data, items=items)
        expected = self.panel.reindex(items=items)
        assert_panel_equal(result, expected)

        result = Panel(data, minor_axis=minor)
        expected = self.panel.reindex(minor=minor)
        assert_panel_equal(result, expected)

    def test_from_dict_mixed_orient(self):
        df = tm.makeDataFrame()
        df['foo'] = 'bar'

        data = {'k1' : df,
                'k2' : df}

        panel = Panel.from_dict(data, orient='minor')

        self.assert_(panel['foo'].values.dtype == np.object_)
        self.assert_(panel['A'].values.dtype == np.float64)

    def test_values(self):
        self.assertRaises(Exception, Panel, np.random.randn(5, 5, 5),
                          range(5), range(5), range(4))

    def test_conform(self):
        df = self.panel['ItemA'][:-5].filter(items=['A', 'B'])
        conformed = self.panel.conform(df)

        assert(conformed.index.equals(self.panel.major_axis))
        assert(conformed.columns.equals(self.panel.minor_axis))

    def test_reindex(self):
        ref = self.panel['ItemB']

        # items
        result = self.panel.reindex(items=['ItemA', 'ItemB'])
        assert_frame_equal(result['ItemB'], ref)

        # major
        new_major = list(self.panel.major_axis[:10])
        result = self.panel.reindex(major=new_major)
        assert_frame_equal(result['ItemB'], ref.reindex(index=new_major))

        # raise exception put both major and major_axis
        self.assertRaises(Exception, self.panel.reindex,
                          major_axis=new_major, major=new_major)

        # minor
        new_minor = list(self.panel.minor_axis[:2])
        result = self.panel.reindex(minor=new_minor)
        assert_frame_equal(result['ItemB'], ref.reindex(columns=new_minor))

        result = self.panel.reindex(items=self.panel.items,
                                    major=self.panel.major_axis,
                                    minor=self.panel.minor_axis)

        assert(result.items is self.panel.items)
        assert(result.major_axis is self.panel.major_axis)
        assert(result.minor_axis is self.panel.minor_axis)

        self.assertRaises(Exception, self.panel.reindex)

        # with filling
        smaller_major = self.panel.major_axis[::5]
        smaller = self.panel.reindex(major=smaller_major)

        larger = smaller.reindex(major=self.panel.major_axis,
                                 method='pad')

        assert_frame_equal(larger.major_xs(self.panel.major_axis[1]),
                           smaller.major_xs(smaller_major[0]))

        # don't necessarily copy
        result = self.panel.reindex(major=self.panel.major_axis, copy=False)
        self.assert_(result is self.panel)

    def test_reindex_like(self):
        # reindex_like
        smaller = self.panel.reindex(items=self.panel.items[:-1],
                                     major=self.panel.major_axis[:-1],
                                     minor=self.panel.minor_axis[:-1])
        smaller_like = self.panel.reindex_like(smaller)
        assert_panel_equal(smaller, smaller_like)

    def test_take(self):
        # axis == 0
        result = self.panel.take([2, 0, 1], axis=0)
        expected = self.panel.reindex(items=['ItemC', 'ItemA', 'ItemB'])
        assert_panel_equal(result, expected)

        # axis >= 1
        result = self.panel.take([3, 0, 1, 2], axis=2)
        expected = self.panel.reindex(minor=['D', 'A', 'B', 'C'])
        assert_panel_equal(result, expected)

        self.assertRaises(Exception, self.panel.take, [3, -1, 1, 2], axis=2)
        self.assertRaises(Exception, self.panel.take, [4, 0, 1, 2], axis=2)

    def test_sort_index(self):
        import random

        ritems = list(self.panel.items)
        rmajor = list(self.panel.major_axis)
        rminor = list(self.panel.minor_axis)
        random.shuffle(ritems)
        random.shuffle(rmajor)
        random.shuffle(rminor)

        random_order = self.panel.reindex(items=ritems)
        sorted_panel = random_order.sort_index(axis=0)
        assert_panel_equal(sorted_panel, self.panel)

        # descending
        random_order = self.panel.reindex(items=ritems)
        sorted_panel = random_order.sort_index(axis=0, ascending=False)
        assert_panel_equal(sorted_panel,
                           self.panel.reindex(items=self.panel.items[::-1]))

        random_order = self.panel.reindex(major=rmajor)
        sorted_panel = random_order.sort_index(axis=1)
        assert_panel_equal(sorted_panel, self.panel)

        random_order = self.panel.reindex(minor=rminor)
        sorted_panel = random_order.sort_index(axis=2)
        assert_panel_equal(sorted_panel, self.panel)

    def test_fillna(self):
        filled = self.panel.fillna(0)
        self.assert_(np.isfinite(filled.values).all())

        filled = self.panel.fillna(method='backfill')
        assert_frame_equal(filled['ItemA'],
                           self.panel['ItemA'].fillna(method='backfill'))

        empty = self.panel.reindex(items=[])
        filled = empty.fillna(0)
        assert_panel_equal(filled, empty)

    def test_swapaxes(self):
        result = self.panel.swapaxes('items', 'minor')
        self.assert_(result.items is self.panel.minor_axis)

        result = self.panel.swapaxes('items', 'major')
        self.assert_(result.items is self.panel.major_axis)

        result = self.panel.swapaxes('major', 'minor')
        self.assert_(result.major_axis is self.panel.minor_axis)

        # this should also work
        result = self.panel.swapaxes(0, 1)
        self.assert_(result.items is self.panel.major_axis)

        # this should also work
        self.assertRaises(Exception, self.panel.swapaxes, 'items', 'items')

    def test_to_frame(self):
        # filtered
        filtered = self.panel.to_frame()
        expected = self.panel.to_frame().dropna(how='any')
        assert_frame_equal(filtered, expected)

        # unfiltered
        unfiltered = self.panel.to_frame(filter_observations=False)
        assert_panel_equal(unfiltered.to_panel(), self.panel)

        # names
        self.assertEqual(unfiltered.index.names, ['major', 'minor'])

    def test_to_frame_mixed(self):
        panel = self.panel.fillna(0)
        panel['str'] = 'foo'
        panel['bool'] = panel['ItemA'] > 0

        lp = panel.to_frame()
        wp = lp.to_panel()
        self.assertEqual(wp['bool'].values.dtype, np.bool_)
        assert_frame_equal(wp['bool'], panel['bool'])

    def test_filter(self):
        pass

    def test_apply(self):
        pass

    def test_compound(self):
        compounded = self.panel.compound()

        assert_series_equal(compounded['ItemA'],
                            (1 + self.panel['ItemA']).product(0) - 1)

    def test_shift(self):
        # major
        idx = self.panel.major_axis[0]
        idx_lag = self.panel.major_axis[1]

        shifted = self.panel.shift(1)

        assert_frame_equal(self.panel.major_xs(idx),
                           shifted.major_xs(idx_lag))

        # minor
        idx = self.panel.minor_axis[0]
        idx_lag = self.panel.minor_axis[1]

        shifted = self.panel.shift(1, axis='minor')

        assert_frame_equal(self.panel.minor_xs(idx),
                           shifted.minor_xs(idx_lag))

        self.assertRaises(Exception, self.panel.shift, 1, axis='items')

    def test_repr_empty(self):
        empty = Panel()
        repr(empty)

    def test_rename(self):
        mapper = {
            'ItemA' : 'foo',
            'ItemB' : 'bar',
            'ItemC' : 'baz'
        }

        renamed = self.panel.rename_axis(mapper, axis=0)
        exp = Index(['foo', 'bar', 'baz'])
        self.assert_(renamed.items.equals(exp))

        renamed = self.panel.rename_axis(str.lower, axis=2)
        exp = Index(['a', 'b', 'c', 'd'])
        self.assert_(renamed.minor_axis.equals(exp))

        # don't copy
        renamed_nocopy = self.panel.rename_axis(mapper, axis=0, copy=False)
        renamed_nocopy['foo'] = 3.
        self.assert_((self.panel['ItemA'].values == 3).all())

    def test_get_attr(self):
        assert_frame_equal(self.panel['ItemA'], self.panel.ItemA)

    def test_group_agg(self):
        values = np.ones((10, 2)) * np.arange(10).reshape((10, 1))
        bounds = np.arange(5) * 2
        f = lambda x: x.mean(axis=0)

        agged = group_agg(values, bounds, f)

        assert(agged[1][0] == 2.5)
        assert(agged[2][0] == 4.5)

        # test a function that doesn't aggregate
        f2 = lambda x: np.zeros((2,2))
        self.assertRaises(Exception, group_agg, values, bounds, f2)

    def test_from_frame_level1_unsorted(self):
        tuples = [('MSFT', 3), ('MSFT', 2), ('AAPL', 2),
                  ('AAPL', 1), ('MSFT', 1)]
        midx = MultiIndex.from_tuples(tuples)
        df = DataFrame(np.random.rand(5,4), index=midx)
        p = df.to_panel()
        assert_frame_equal(p.minor_xs(2), df.ix[:,2].sort_index())

    def test_to_excel(self):
        try:
            import xlwt
            import xlrd
            import openpyxl
        except ImportError:
            raise nose.SkipTest

        path = '__tmp__.xlsx'
        self.panel.to_excel(path)
        reader = ExcelFile(path)
        for item, df in self.panel.iteritems():
            recdf = reader.parse(str(item),index_col=0)
            assert_frame_equal(df, recdf)

class TestLongPanel(unittest.TestCase):
    """
    LongPanel no longer exists, but...
    """

    def setUp(self):
        panel = tm.makePanel()
        tm.add_nans(panel)

        self.panel = panel.to_frame()
        self.unfiltered_panel = panel.to_frame(filter_observations=False)

    def test_ops_differently_indexed(self):
        # trying to set non-identically indexed panel
        wp = self.panel.to_panel()
        wp2 = wp.reindex(major=wp.major_axis[:-1])
        lp2 = wp2.to_frame()

        result = self.panel + lp2
        assert_frame_equal(result.reindex(lp2.index), lp2 * 2)

        # careful, mutation
        self.panel['foo'] = lp2['ItemA']
        assert_series_equal(self.panel['foo'].reindex(lp2.index),
                            lp2['ItemA'])

    def test_ops_scalar(self):
        result = self.panel.mul(2)
        expected = DataFrame.__mul__(self.panel, 2)
        assert_frame_equal(result, expected)

    def test_combineFrame(self):
        wp = self.panel.to_panel()
        result = self.panel.add(wp['ItemA'].stack(), axis=0)
        assert_frame_equal(result.to_panel()['ItemA'], wp['ItemA'] * 2)

    def test_combinePanel(self):
        wp = self.panel.to_panel()
        result = self.panel.add(self.panel)
        wide_result = result.to_panel()
        assert_frame_equal(wp['ItemA'] * 2, wide_result['ItemA'])

        # one item
        result = self.panel.add(self.panel.filter(['ItemA']))

    def test_combine_scalar(self):
        result = self.panel.mul(2)
        expected = DataFrame(self.panel._data) * 2
        assert_frame_equal(result, expected)

    def test_combine_series(self):
        s = self.panel['ItemA'][:10]
        result = self.panel.add(s, axis=0)
        expected = DataFrame.add(self.panel, s, axis=0)
        assert_frame_equal(result, expected)

        s = self.panel.ix[5]
        result = self.panel + s
        expected = DataFrame.add(self.panel, s, axis=1)
        assert_frame_equal(result, expected)

    def test_operators(self):
        wp = self.panel.to_panel()
        result = (self.panel + 1).to_panel()
        assert_frame_equal(wp['ItemA'] + 1, result['ItemA'])

    def test_sort(self):
        def is_sorted(arr):
            return (arr[1:] > arr[:-1]).any()

        sorted_minor = self.panel.sortlevel(level=1)
        self.assert_(is_sorted(sorted_minor.index.labels[1]))

        sorted_major = sorted_minor.sortlevel(level=0)
        self.assert_(is_sorted(sorted_major.index.labels[0]))

    def test_to_string(self):
        from cStringIO import StringIO

        buf = StringIO()
        self.panel.to_string(buf)

    def test_truncate(self):
        dates = self.panel.index.levels[0]
        start, end = dates[1], dates[5]

        trunced = self.panel.truncate(start, end).to_panel()
        expected = self.panel.to_panel()['ItemA'].truncate(start, end)

        assert_frame_equal(trunced['ItemA'], expected)

        trunced = self.panel.truncate(before=start).to_panel()
        expected = self.panel.to_panel()['ItemA'].truncate(before=start)

        assert_frame_equal(trunced['ItemA'], expected)

        trunced = self.panel.truncate(after=end).to_panel()
        expected = self.panel.to_panel()['ItemA'].truncate(after=end)

        assert_frame_equal(trunced['ItemA'], expected)

        # truncate on dates that aren't in there
        wp = self.panel.to_panel()
        new_index = wp.major_axis[::5]

        wp2 = wp.reindex(major=new_index)

        lp2 = wp2.to_frame()
        lp_trunc = lp2.truncate(wp.major_axis[2], wp.major_axis[-2])

        wp_trunc = wp2.truncate(wp.major_axis[2], wp.major_axis[-2])

        assert_panel_equal(wp_trunc, lp_trunc.to_panel())

        # throw proper exception
        self.assertRaises(Exception, lp2.truncate, wp.major_axis[-2],
                          wp.major_axis[2])

    def test_axis_dummies(self):
        from pandas.core.reshape import make_axis_dummies

        minor_dummies = make_axis_dummies(self.panel, 'minor')
        self.assertEqual(len(minor_dummies.columns),
                         len(self.panel.index.levels[1]))

        major_dummies = make_axis_dummies(self.panel, 'major')
        self.assertEqual(len(major_dummies.columns),
                         len(self.panel.index.levels[0]))

        mapping = {'A' : 'one',
                   'B' : 'one',
                   'C' : 'two',
                   'D' : 'two'}

        transformed = make_axis_dummies(self.panel, 'minor',
                                                  transform=mapping.get)
        self.assertEqual(len(transformed.columns), 2)
        self.assert_(np.array_equal(transformed.columns, ['one', 'two']))

        # TODO: test correctness

    def test_get_dummies(self):
        from pandas.core.reshape import make_column_dummies, make_axis_dummies

        self.panel['Label'] = self.panel.index.labels[1]
        minor_dummies = make_axis_dummies(self.panel, 'minor')
        dummies = make_column_dummies(self.panel, 'Label')
        self.assert_(np.array_equal(dummies.values, minor_dummies.values))

    def test_apply(self):
        # ufunc
        applied = self.panel.apply(np.sqrt)
        self.assert_(assert_almost_equal(applied.values,
                                         np.sqrt(self.panel.values)))

    def test_mean(self):
        means = self.panel.mean(level='minor')

        # test versus Panel version
        wide_means = self.panel.to_panel().mean('major')
        assert_frame_equal(means, wide_means)

    def test_sum(self):
        sums = self.panel.sum(level='minor')

        # test versus Panel version
        wide_sums = self.panel.to_panel().sum('major')
        assert_frame_equal(sums, wide_sums)

    def test_count(self):
        index = self.panel.index

        major_count = self.panel.count(level=0)['ItemA']
        labels = index.labels[0]
        for i, idx in enumerate(index.levels[0]):
            self.assertEqual(major_count[i], (labels == i).sum())

        minor_count = self.panel.count(level=1)['ItemA']
        labels = index.labels[1]
        for i, idx in enumerate(index.levels[1]):
            self.assertEqual(minor_count[i], (labels == i).sum())

    def test_join(self):
        lp1 = self.panel.filter(['ItemA', 'ItemB'])
        lp2 = self.panel.filter(['ItemC'])

        joined = lp1.join(lp2)

        self.assertEqual(len(joined.columns), 3)

        self.assertRaises(Exception, lp1.join,
                          self.panel.filter(['ItemB', 'ItemC']))

    def test_pivot(self):
        from pandas.core.reshape import _slow_pivot

        one, two, three = (np.array([1, 2, 3, 4, 5]),
                           np.array(['a', 'b', 'c', 'd', 'e']),
                           np.array([1, 2, 3, 5, 4.]))
        df = pivot(one, two, three)
        self.assertEqual(df['a'][1], 1)
        self.assertEqual(df['b'][2], 2)
        self.assertEqual(df['c'][3], 3)
        self.assertEqual(df['d'][4], 5)
        self.assertEqual(df['e'][5], 4)
        assert_frame_equal(df, _slow_pivot(one, two, three))

        # weird overlap, TODO: test?
        a, b, c = (np.array([1, 2, 3, 4, 4]),
                   np.array(['a', 'a', 'a', 'a', 'a']),
                   np.array([1., 2., 3., 4., 5.]))
        self.assertRaises(Exception, pivot, a, b, c)

        # corner case, empty
        df = pivot(np.array([]), np.array([]), np.array([]))

def test_monotonic():
    pos = np.array([1, 2, 3, 5])

    assert panelm._monotonic(pos)

    neg = np.array([1, 2, 3, 4, 3])

    assert not panelm._monotonic(neg)

    neg2 = np.array([5, 1, 2, 3, 4, 5])

    assert not panelm._monotonic(neg2)

def test_panel_index():
    index = panelm.panel_index([1,2,3,4], [1,2,3])
    expected = MultiIndex.from_arrays([np.tile([1,2,3,4], 3),
                                       np.repeat([1,2,3], 4)])
    assert(index.equals(expected))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
