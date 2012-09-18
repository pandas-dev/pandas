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
from pandas.core.fdpanel import FDPanel
from pandas.core.series import remove_na
import pandas.core.common as com
import pandas.core.panel as panelmod
from pandas.util import py3compat
from pandas.io.parsers import (ExcelFile, ExcelWriter)

from pandas.util.testing import (assert_panel_equal,
                                 assert_fdpanel_equal,
                                 assert_frame_equal,
                                 assert_series_equal,
                                 assert_almost_equal)
import pandas.util.testing as tm

def add_nans(fdp):
    for l, label in enumerate(fdp.labels):
        panel = fdp[label]
        tm.add_nans(panel)

class SafeForLongAndSparse(object):

    def test_repr(self):
        foo = repr(self.fdpanel)

    def test_iter(self):
        tm.equalContents(list(self.fdpanel), self.fdpanel.labels)

    def test_count(self):
        f = lambda s: notnull(s).sum()
        self._check_stat_op('count', f, obj=self.fdpanel, has_skipna=False)

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
            obj = self.fdpanel

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
                assert_panel_equal(result, obj.apply(wrapper, axis=i))
        else:
            skipna_wrapper = alternative
            wrapper = alternative

        for i in range(obj.ndim):
            result = f(axis=i)
            assert_panel_equal(result, obj.apply(skipna_wrapper, axis=i))

        self.assertRaises(Exception, f, axis=obj.ndim)

class SafeForSparse(object):

    @classmethod
    def assert_panel_equal(cls, x, y):
        assert_panel_equal(x, y)

    @classmethod
    def assert_fdpanel_equal(cls, x, y):
        assert_fdpanel_equal(x, y)

    def test_get_axis(self):
        assert(self.fdpanel._get_axis(0) is self.fdpanel.labels)
        assert(self.fdpanel._get_axis(1) is self.fdpanel.items)
        assert(self.fdpanel._get_axis(2) is self.fdpanel.major_axis)
        assert(self.fdpanel._get_axis(3) is self.fdpanel.minor_axis)

    def test_set_axis(self):
        new_labels = Index(np.arange(len(self.fdpanel.labels)))
        new_items  = Index(np.arange(len(self.fdpanel.items)))
        new_major  = Index(np.arange(len(self.fdpanel.major_axis)))
        new_minor  = Index(np.arange(len(self.fdpanel.minor_axis)))

        # ensure propagate to potentially prior-cached items too
        label = self.fdpanel['l1']
        self.fdpanel.labels = new_labels

        if hasattr(self.fdpanel, '_item_cache'):
            self.assert_('l1' not in self.fdpanel._item_cache)
        self.assert_(self.fdpanel.labels is new_labels)

        self.fdpanel.major_axis = new_major
        self.assert_(self.fdpanel[0].major_axis is new_major)
        self.assert_(self.fdpanel.major_axis is new_major)

        self.fdpanel.minor_axis = new_minor
        self.assert_(self.fdpanel[0].minor_axis is new_minor)
        self.assert_(self.fdpanel.minor_axis is new_minor)

    def test_get_axis_number(self):
        self.assertEqual(self.fdpanel._get_axis_number('labels'), 0)
        self.assertEqual(self.fdpanel._get_axis_number('items'), 1)
        self.assertEqual(self.fdpanel._get_axis_number('major'), 2)
        self.assertEqual(self.fdpanel._get_axis_number('minor'), 3)

    def test_get_axis_name(self):
        self.assertEqual(self.fdpanel._get_axis_name(0), 'labels')
        self.assertEqual(self.fdpanel._get_axis_name(1), 'items')
        self.assertEqual(self.fdpanel._get_axis_name(2), 'major_axis')
        self.assertEqual(self.fdpanel._get_axis_name(3), 'minor_axis')

    #def test_get_plane_axes(self):
    #    # what to do here?

    #    index, columns = self.panel._get_plane_axes('items')
    #    index, columns = self.panel._get_plane_axes('major_axis')
    #    index, columns = self.panel._get_plane_axes('minor_axis')
    #    index, columns = self.panel._get_plane_axes(0)

    def test_truncate(self):
        raise nose.SkipTest

        #dates = self.panel.major_axis
        #start, end = dates[1], dates[5]

        #trunced = self.panel.truncate(start, end, axis='major')
        #expected = self.panel['ItemA'].truncate(start, end)

        #assert_frame_equal(trunced['ItemA'], expected)

        #trunced = self.panel.truncate(before=start, axis='major')
        #expected = self.panel['ItemA'].truncate(before=start)

        #assert_frame_equal(trunced['ItemA'], expected)

        #trunced = self.panel.truncate(after=end, axis='major')
        #expected = self.panel['ItemA'].truncate(after=end)

        #assert_frame_equal(trunced['ItemA'], expected)

        # XXX test other axes

    def test_arith(self):
        self._test_op(self.fdpanel, operator.add)
        self._test_op(self.fdpanel, operator.sub)
        self._test_op(self.fdpanel, operator.mul)
        self._test_op(self.fdpanel, operator.truediv)
        self._test_op(self.fdpanel, operator.floordiv)
        self._test_op(self.fdpanel, operator.pow)

        self._test_op(self.fdpanel, lambda x, y: y + x)
        self._test_op(self.fdpanel, lambda x, y: y - x)
        self._test_op(self.fdpanel, lambda x, y: y * x)
        self._test_op(self.fdpanel, lambda x, y: y / x)
        self._test_op(self.fdpanel, lambda x, y: y ** x)

        self.assertRaises(Exception, self.fdpanel.__add__, self.fdpanel['l1'])

    @staticmethod
    def _test_op(fdpanel, op):
        result = op(fdpanel, 1)
        assert_panel_equal(result['l1'], op(fdpanel['l1'], 1))

    def test_keys(self):
        tm.equalContents(self.fdpanel.keys(), self.fdpanel.labels)

    def test_iteritems(self):
        """Test fdpanel.iteritems(), aka fdpanel.iterkv()"""
        # just test that it works
        for k, v in self.fdpanel.iterkv():
            pass

        self.assertEqual(len(list(self.fdpanel.iterkv())),
                         len(self.fdpanel.labels))

    def test_combineFDPanel(self):
        result = self.fdpanel.add(self.fdpanel)
        self.assert_fdpanel_equal(result, self.fdpanel * 2)

    def test_neg(self):
        self.assert_fdpanel_equal(-self.fdpanel, self.fdpanel * -1)

    def test_select(self):
        p = self.fdpanel

        # select labels
        result = p.select(lambda x: x in ('l1', 'l3'), axis='labels')
        expected = p.reindex(labels=['l1','l3'])
        self.assert_fdpanel_equal(result, expected)

        # select items
        result = p.select(lambda x: x in ('ItemA', 'ItemC'), axis='items')
        expected = p.reindex(items=['ItemA', 'ItemC'])
        self.assert_fdpanel_equal(result, expected)

        # select major_axis
        result = p.select(lambda x: x >= datetime(2000, 1, 15), axis='major')
        new_major = p.major_axis[p.major_axis >= datetime(2000, 1, 15)]
        expected = p.reindex(major=new_major)
        self.assert_fdpanel_equal(result, expected)

        # select minor_axis
        result = p.select(lambda x: x in ('D', 'A'), axis=3)
        expected = p.reindex(minor=['A', 'D'])
        self.assert_fdpanel_equal(result, expected)

        # corner case, empty thing
        result = p.select(lambda x: x in ('foo',), axis='items')
        self.assert_fdpanel_equal(result, p.reindex(items=[]))

    def test_get_value(self):
        for item in self.panel.items:
            for mjr in self.panel.major_axis[::2]:
                for mnr in self.panel.minor_axis:
                    result = self.panel.get_value(item, mjr, mnr)
                    expected = self.panel[item][mnr][mjr]
                    assert_almost_equal(result, expected)

    def test_abs(self):
        result = self.fdpanel.abs()
        expected = np.abs(self.fdpanel)
        self.assert_fdpanel_equal(result, expected)

        p = self.fdpanel['l1']
        result = p.abs()
        expected = np.abs(p)
        assert_panel_equal(result, expected)

        df = p['ItemA']
        result = df.abs()
        expected = np.abs(df)
        assert_frame_equal(result, expected)

class CheckIndexing(object):


    def test_getitem(self):
        self.assertRaises(Exception, self.fdpanel.__getitem__, 'ItemQ')

    def test_delitem_and_pop(self):
        expected = self.fdpanel['l2']
        result   = self.fdpanel.pop('l2')
        assert_panel_equal(expected, result)
        self.assert_('l2' not in self.fdpanel.labels)

        del self.fdpanel['l3']
        self.assert_('l3' not in self.fdpanel.labels)
        self.assertRaises(Exception, self.fdpanel.__delitem__, 'l3')

        values = np.empty((4, 4, 4, 4))
        values[0] = 0
        values[1] = 1
        values[2] = 2
        values[3] = 3

        fdpanel = FDPanel(values, range(4), range(4), range(4), range(4))

        # did we delete the right row?

        fdpanelc = fdpanel.copy()
        del fdpanelc[0]
        assert_panel_equal(fdpanelc[1], fdpanel[1])
        assert_panel_equal(fdpanelc[2], fdpanel[2])
        assert_panel_equal(fdpanelc[3], fdpanel[3])

        fdpanelc = fdpanel.copy()
        del fdpanelc[1]
        assert_panel_equal(fdpanelc[0], fdpanel[0])
        assert_panel_equal(fdpanelc[2], fdpanel[2])
        assert_panel_equal(fdpanelc[3], fdpanel[3])

        fdpanelc = fdpanel.copy()
        del fdpanelc[2]
        assert_panel_equal(fdpanelc[1], fdpanel[1])
        assert_panel_equal(fdpanelc[0], fdpanel[0])
        assert_panel_equal(fdpanelc[3], fdpanel[3])

        fdpanelc = fdpanel.copy()
        del fdpanelc[3]
        assert_panel_equal(fdpanelc[1], fdpanel[1])
        assert_panel_equal(fdpanelc[2], fdpanel[2])
        assert_panel_equal(fdpanelc[0], fdpanel[0])

    def test_setitem(self):
        ## LongPanel with one item
        #lp = self.panel.filter(['ItemA', 'ItemB']).to_frame()
        #self.assertRaises(Exception, self.panel.__setitem__,
        #                  'ItemE', lp)

        # Panel
        p = Panel(dict(ItemA = self.fdpanel['l1']['ItemA'][2:].filter(items=['A', 'B'])))
        self.fdpanel['l4'] = p
        self.fdpanel['l5'] = p

        p2 = self.fdpanel['l4']

        assert_panel_equal(p, p2.reindex(items      = p.items,
                                         major_axis = p.major_axis,
                                         minor_axis = p.minor_axis))

        # scalar
        self.fdpanel['lG'] = 1
        self.fdpanel['lE'] = True
        self.assert_(self.fdpanel['lG'].values.dtype == np.int64)
        self.assert_(self.fdpanel['lE'].values.dtype == np.bool_)

        # object dtype
        self.fdpanel['lQ'] = 'foo'
        self.assert_(self.fdpanel['lQ'].values.dtype == np.object_)

        # boolean dtype
        self.fdpanel['lP'] = self.fdpanel['l1'] > 0
        self.assert_(self.fdpanel['lP'].values.dtype == np.bool_)

    def test_setitem_ndarray(self):
        raise nose.SkipTest
    #    from pandas import DateRange, datetools

    #    timeidx = DateRange(start=datetime(2009,1,1),
    #                        end=datetime(2009,12,31),
    #                        offset=datetools.MonthEnd())
    #    lons_coarse = np.linspace(-177.5, 177.5, 72)
    #    lats_coarse = np.linspace(-87.5, 87.5, 36)
    #    P = Panel(items=timeidx, major_axis=lons_coarse, minor_axis=lats_coarse)
    #    data = np.random.randn(72*36).reshape((72,36))
    #    key = datetime(2009,2,28)
    #    P[key] = data#

    #    assert_almost_equal(P[key].values, data)

    def test_major_xs(self):
        ref = self.fdpanel['l1']['ItemA']

        idx = self.fdpanel.major_axis[5]
        xs  = self.fdpanel.major_xs(idx)
        
        assert_series_equal(xs['l1'].T['ItemA'], ref.xs(idx))

        # not contained
        idx = self.fdpanel.major_axis[0] - bday
        self.assertRaises(Exception, self.fdpanel.major_xs, idx)

    def test_major_xs_mixed(self):
        self.fdpanel['l4'] = 'foo'
        xs = self.fdpanel.major_xs(self.fdpanel.major_axis[0])
        self.assert_(xs['l1']['A'].dtype == np.float64)
        self.assert_(xs['l4']['A'].dtype == np.object_)

    def test_minor_xs(self):
        ref = self.fdpanel['l1']['ItemA']

        idx = self.fdpanel.minor_axis[1]
        xs = self.fdpanel.minor_xs(idx)

        assert_series_equal(xs['l1'].T['ItemA'], ref[idx])
            
        # not contained
        self.assertRaises(Exception, self.fdpanel.minor_xs, 'E')

    def test_minor_xs_mixed(self):
        self.fdpanel['l4'] = 'foo'

        xs = self.fdpanel.minor_xs('D')
        self.assert_(xs['l1'].T['ItemA'].dtype == np.float64)
        self.assert_(xs['l4'].T['ItemA'].dtype == np.object_)

    def test_xs(self):
        l1 = self.fdpanel.xs('l1', axis=0)
        expected = self.fdpanel['l1']
        assert_panel_equal(l1, expected)

        # not view by default
        l1.values[:] = np.nan
        self.assert_(not np.isnan(self.fdpanel['l1'].values).all())

        # but can get view
        l1_view = self.fdpanel.xs('l1', axis=0, copy=False)
        l1_view.values[:] = np.nan
        self.assert_(np.isnan(self.fdpanel['l1'].values).all())

        # mixed-type
        self.fdpanel['strings'] = 'foo'
        self.assertRaises(Exception, self.fdpanel.xs, 'D', axis=2,
                          copy=False)

    def test_getitem_fancy_labels(self):
        fdp = self.fdpanel

        labels = fdp.labels[[1, 0]]
        items  = fdp.items[[1, 0]]
        dates  = fdp.major_axis[::2]
        cols   = ['D', 'C', 'F']

        # all 4 specified
        assert_fdpanel_equal(fdp.ix[labels, items, dates, cols],
                             fdp.reindex(labels=labels, items=items, major=dates, minor=cols))

        # 3 specified
        assert_fdpanel_equal(fdp.ix[:, items, dates, cols],
                             fdp.reindex(items=items, major=dates, minor=cols))

        # 2 specified
        assert_fdpanel_equal(fdp.ix[:, :, dates, cols],
                             fdp.reindex(major=dates, minor=cols))

        assert_fdpanel_equal(fdp.ix[:, items, :, cols],
                             fdp.reindex(items=items, minor=cols))

        assert_fdpanel_equal(fdp.ix[:, items, dates, :],
                             fdp.reindex(items=items, major=dates))

        # only 1
        assert_fdpanel_equal(fdp.ix[:, items, :, :],
                             fdp.reindex(items=items))

        assert_fdpanel_equal(fdp.ix[:, :, dates, :],
                             fdp.reindex(major=dates))
        
        assert_fdpanel_equal(fdp.ix[:, :, :, cols],
                             fdp.reindex(minor=cols))

    def test_getitem_fancy_slice(self):
        pass

    def test_getitem_fancy_ints(self):
        pass

    def test_getitem_fancy_xs(self):
        raise nose.SkipTest
        #self.assertRaises(NotImplementedError, self.fdpanel.major_xs)
        #self.assertRaises(NotImplementedError, self.fdpanel.minor_xs)

    def test_getitem_fancy_xs_check_view(self):
        raise nose.SkipTest
    #    item = 'ItemB'
    #    date = self.panel.major_axis[5]
    #    col = 'C'

    #    # make sure it's always a view
    #    NS = slice(None, None)

    #    # DataFrames
    #    comp = assert_frame_equal
    #    self._check_view(item, comp)
    #    self._check_view((item, NS), comp)
    #    self._check_view((item, NS, NS), comp)
    #    self._check_view((NS, date), comp)
    #    self._check_view((NS, date, NS), comp)
    #    self._check_view((NS, NS, 'C'), comp)

    #    # Series
    #    comp = assert_series_equal
    #    self._check_view((item, date), comp)
    #    self._check_view((item, date, NS), comp)
    #    self._check_view((item, NS, 'C'), comp)
    #    self._check_view((NS, date, 'C'), comp)#

    #def _check_view(self, indexer, comp):
    #    cp = self.panel.copy()
    #    obj = cp.ix[indexer]
    #    obj.values[:] = 0
    #    self.assert_((obj.values == 0).all())
    #    comp(cp.ix[indexer].reindex_like(obj), obj)

    def test_get_value(self):
        for label in self.fdpanel.labels:
            for item in self.fdpanel.items:
                for mjr in self.fdpanel.major_axis[::2]:
                    for mnr in self.fdpanel.minor_axis:
                        result   = self.fdpanel.get_value(label, item, mjr, mnr)
                        expected = self.fdpanel[label][item][mnr][mjr]
                        assert_almost_equal(result, expected)
                        
    def test_set_value(self):
        for label in self.fdpanel.labels:
            for item in self.fdpanel.items:
                for mjr in self.fdpanel.major_axis[::2]:
                    for mnr in self.fdpanel.minor_axis:
                        self.fdpanel.set_value(label, item, mjr, mnr, 1.)
                        assert_almost_equal(self.fdpanel[label][item][mnr][mjr], 1.)

        # resize
        res = self.fdpanel.set_value('l4', 'ItemE', 'foo', 'bar', 1.5)
        self.assert_(isinstance(res, FDPanel))
        self.assert_(res is not self.fdpanel)
        self.assertEqual(res.get_value('l4', 'ItemE', 'foo', 'bar'), 1.5)

        res3 = self.fdpanel.set_value('l4', 'ItemE', 'foobar', 'baz', 5)
        self.assert_(com.is_float_dtype(res3['l4'].values))

class TestFDPanel(unittest.TestCase, CheckIndexing, SafeForSparse, SafeForLongAndSparse):

    @classmethod
    def assert_fdpanel_equal(cls,x, y):
        assert_fdpanel_equal(x, y)

    def setUp(self):
        self.fdpanel = tm.makeFDPanel()
        add_nans(self.fdpanel)

    def test_constructor(self):
        # with BlockManager
        fdp = FDPanel(self.fdpanel._data)
        self.assert_(fdp._data is self.fdpanel._data)

        fdp = FDPanel(self.fdpanel._data, copy=True)
        self.assert_(fdp._data is not self.fdpanel._data)
        assert_fdpanel_equal(fdp, self.fdpanel)

        # strings handled prop
        #fdp = FDPanel([[['foo', 'foo', 'foo',],
        #                 ['foo', 'foo', 'foo']]])
        #self.assert_(wp.values.dtype == np.object_)

        vals = self.fdpanel.values

        # no copy
        fdp = FDPanel(vals)
        self.assert_(fdp.values is vals)

        # copy
        fdp = FDPanel(vals, copy=True)
        self.assert_(fdp.values is not vals)

    def test_constructor_cast(self):
        zero_filled = self.fdpanel.fillna(0)

        casted = FDPanel(zero_filled._data, dtype=int)
        casted2 = FDPanel(zero_filled.values, dtype=int)

        exp_values = zero_filled.values.astype(int)
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
        self.assert_(self.fdpanel._data.is_consolidated())

        self.fdpanel['foo'] = 1.
        self.assert_(not self.fdpanel._data.is_consolidated())

        fdpanel = self.fdpanel.consolidate()
        self.assert_(fdpanel._data.is_consolidated())

    def test_ctor_dict(self):
        l1 = self.fdpanel['l1']
        l2 = self.fdpanel['l2']

        d  = {'A' : l1, 'B' : l2.ix[['ItemB'],:,:] }
        #d2 = {'A' : itema._series, 'B' : itemb[5:]._series}
        #d3 = {'A' : DataFrame(itema._series),
        #      'B' : DataFrame(itemb[5:]._series)}

        fdp = FDPanel(d)
        #wp2 = Panel.from_dict(d2) # nested Dict
        #wp3 = Panel.from_dict(d3)
        #self.assert_(wp.major_axis.equals(self.panel.major_axis))
        assert_panel_equal(fdp['A'], self.fdpanel['l1'])
        assert_frame_equal(fdp.ix['B','ItemB',:,:], self.fdpanel.ix['l2',['ItemB'],:,:]['ItemB'])

        # intersect
        #wp = Panel.from_dict(d, intersect=True)
        #self.assert_(wp.major_axis.equals(itemb.index[5:]))

        # use constructor
        #assert_panel_equal(Panel(d), Panel.from_dict(d))
        #assert_panel_equal(Panel(d2), Panel.from_dict(d2))
        #assert_panel_equal(Panel(d3), Panel.from_dict(d3))

        # cast
        #dcasted = dict((k, v.reindex(wp.major_axis).fillna(0))
        #               for k, v in d.iteritems())
        #result = Panel(dcasted, dtype=int)
        #expected = Panel(dict((k, v.astype(int))
        #                      for k, v in dcasted.iteritems()))
        #assert_panel_equal(result, expected)

    def test_constructor_dict_mixed(self):
        data = dict((k, v.values) for k, v in self.fdpanel.iterkv())
        result = FDPanel(data)
        exp_major = Index(np.arange(len(self.fdpanel.major_axis)))
        self.assert_(result.major_axis.equals(exp_major))

        result = FDPanel(data, 
                         labels     = self.fdpanel.labels,
                         items      = self.fdpanel.items,
                         major_axis = self.fdpanel.major_axis,
                         minor_axis = self.fdpanel.minor_axis)
        assert_fdpanel_equal(result, self.fdpanel)

        data['l2'] = self.fdpanel['l2']
        result = FDPanel(data)
        assert_fdpanel_equal(result, self.fdpanel)

        # corner, blow up
        data['l2'] = data['l2']['ItemB']
        self.assertRaises(Exception, FDPanel, data)

        data['l2'] = self.fdpanel['l2'].values[:, :, :-1]
        self.assertRaises(Exception, FDPanel, data)

    def test_constructor_resize(self):
        data  = self.fdpanel._data
        labels= self.fdpanel.labels[:-1]
        items = self.fdpanel.items[:-1]
        major = self.fdpanel.major_axis[:-1]
        minor = self.fdpanel.minor_axis[:-1]

        result = FDPanel(data, labels=labels, items=items, major_axis=major, minor_axis=minor)
        expected = self.fdpanel.reindex(labels=labels, items=items, major=major, minor=minor)
        assert_fdpanel_equal(result, expected)

        result = FDPanel(data, items=items, major_axis=major)
        expected = self.fdpanel.reindex(items=items, major=major)
        assert_fdpanel_equal(result, expected)

        result = FDPanel(data, items=items)
        expected = self.fdpanel.reindex(items=items)
        assert_fdpanel_equal(result, expected)

        result = FDPanel(data, minor_axis=minor)
        expected = self.fdpanel.reindex(minor=minor)
        assert_fdpanel_equal(result, expected)

    def test_from_dict_mixed_orient(self):
        raise nose.SkipTest
    #    df = tm.makeDataFrame()
    #    df['foo'] = 'bar'

    #    data = {'k1' : df,
    #            'k2' : df}

    #    panel = Panel.from_dict(data, orient='minor')

    #    self.assert_(panel['foo'].values.dtype == np.object_)
    #    self.assert_(panel['A'].values.dtype == np.float64)

    def test_values(self):
        self.assertRaises(Exception, Panel, np.random.randn(5, 5, 5),
                          range(5), range(5), range(4))

    def test_conform(self):
        p = self.fdpanel['l1'].filter(items=['ItemA', 'ItemB'])
        conformed = self.fdpanel.conform(p)

        assert(conformed.items.equals(self.fdpanel.labels))
        assert(conformed.major_axis.equals(self.fdpanel.major_axis))
        assert(conformed.minor_axis.equals(self.fdpanel.minor_axis))

    def test_reindex(self):
        ref = self.fdpanel['l2']

        # labels
        result = self.fdpanel.reindex(labels=['l1','l2'])
        assert_panel_equal(result['l2'], ref)

        # items
        result = self.fdpanel.reindex(items=['ItemA', 'ItemB'])
        assert_frame_equal(result['l2']['ItemB'], ref['ItemB'])

        # major
        new_major = list(self.fdpanel.major_axis[:10])
        result = self.fdpanel.reindex(major=new_major)
        assert_frame_equal(result['l2']['ItemB'], ref['ItemB'].reindex(index=new_major))

        # raise exception put both major and major_axis
        self.assertRaises(Exception, self.fdpanel.reindex,
                          major_axis=new_major, major=new_major)

        # minor
        new_minor = list(self.fdpanel.minor_axis[:2])
        result = self.fdpanel.reindex(minor=new_minor)
        assert_frame_equal(result['l2']['ItemB'], ref['ItemB'].reindex(columns=new_minor))
        
        result = self.fdpanel.reindex(labels=self.fdpanel.labels,
                                      items =self.fdpanel.items,
                                      major =self.fdpanel.major_axis,
                                      minor =self.fdpanel.minor_axis)
        
        assert(result.labels is self.fdpanel.labels)
        assert(result.items is self.fdpanel.items)
        assert(result.major_axis is self.fdpanel.major_axis)
        assert(result.minor_axis is self.fdpanel.minor_axis)

        self.assertRaises(Exception, self.fdpanel.reindex)

        # with filling
        smaller_major = self.fdpanel.major_axis[::5]
        smaller = self.fdpanel.reindex(major=smaller_major)

        larger = smaller.reindex(major=self.fdpanel.major_axis,
                                 method='pad')

        assert_panel_equal(larger.ix[:,:,self.fdpanel.major_axis[1],:],
                           smaller.ix[:,:,smaller_major[0],:])
        
        # don't necessarily copy
        result = self.fdpanel.reindex(major=self.fdpanel.major_axis, copy=False)
        self.assert_(result is self.fdpanel)

    def test_reindex_like(self):
        # reindex_like
        smaller = self.fdpanel.reindex(labels=self.fdpanel.labels[:-1],
                                       items =self.fdpanel.items[:-1],
                                       major =self.fdpanel.major_axis[:-1],
                                       minor =self.fdpanel.minor_axis[:-1])
        smaller_like = self.fdpanel.reindex_like(smaller)
        assert_fdpanel_equal(smaller, smaller_like)

    def test_take(self):
        raise nose.SkipTest
    #    # axis == 0
    #    result = self.panel.take([2, 0, 1], axis=0)
    #    expected = self.panel.reindex(items=['ItemC', 'ItemA', 'ItemB'])
    #    assert_panel_equal(result, expected)#

    #    # axis >= 1
    #    result = self.panel.take([3, 0, 1, 2], axis=2)
    #    expected = self.panel.reindex(minor=['D', 'A', 'B', 'C'])
    #    assert_panel_equal(result, expected)

    #    self.assertRaises(Exception, self.panel.take, [3, -1, 1, 2], axis=2)
    #    self.assertRaises(Exception, self.panel.take, [4, 0, 1, 2], axis=2)

    def test_sort_index(self):
        import random

        rlabels= list(self.fdpanel.labels)
        ritems = list(self.fdpanel.items)
        rmajor = list(self.fdpanel.major_axis)
        rminor = list(self.fdpanel.minor_axis)
        random.shuffle(rlabels)
        random.shuffle(ritems)
        random.shuffle(rmajor)
        random.shuffle(rminor)

        random_order = self.fdpanel.reindex(labels=rlabels)
        sorted_fdpanel = random_order.sort_index(axis=0)
        assert_fdpanel_equal(sorted_fdpanel, self.fdpanel)

        # descending
        #random_order = self.panel.reindex(items=ritems)
        #sorted_panel = random_order.sort_index(axis=0, ascending=False)
        #assert_panel_equal(sorted_panel,
        #                   self.panel.reindex(items=self.panel.items[::-1]))

        #random_order = self.panel.reindex(major=rmajor)
        #sorted_panel = random_order.sort_index(axis=1)
        #assert_panel_equal(sorted_panel, self.panel)

        #random_order = self.panel.reindex(minor=rminor)
        #sorted_panel = random_order.sort_index(axis=2)
        #assert_panel_equal(sorted_panel, self.panel)

    def test_fillna(self):
        filled = self.fdpanel.fillna(0)
        self.assert_(np.isfinite(filled.values).all())

        filled = self.fdpanel.fillna(method='backfill')
        assert_panel_equal(filled['l1'],
                           self.fdpanel['l1'].fillna(method='backfill'))

        fdpanel = self.fdpanel.copy()
        fdpanel['str'] = 'foo'

        filled = fdpanel.fillna(method='backfill')
        assert_panel_equal(filled['l1'],
                           fdpanel['l1'].fillna(method='backfill'))

        empty = self.fdpanel.reindex(labels=[])
        filled = empty.fillna(0)
        assert_fdpanel_equal(filled, empty)

    def test_swapaxes(self):
        result = self.fdpanel.swapaxes('labels','items')
        self.assert_(result.items is self.fdpanel.labels)

        result = self.fdpanel.swapaxes('labels','minor')
        self.assert_(result.labels is self.fdpanel.minor_axis)

        result = self.fdpanel.swapaxes('items', 'minor')
        self.assert_(result.items is self.fdpanel.minor_axis)

        result = self.fdpanel.swapaxes('items', 'major')
        self.assert_(result.items is self.fdpanel.major_axis)

        result = self.fdpanel.swapaxes('major', 'minor')
        self.assert_(result.major_axis is self.fdpanel.minor_axis)

        # this should also work
        result = self.fdpanel.swapaxes(0, 1)
        self.assert_(result.labels is self.fdpanel.items)

        # this should also work
        self.assertRaises(Exception, self.fdpanel.swapaxes, 'items', 'items')

    def test_to_frame(self):
        raise nose.SkipTest
    #    # filtered
    #    filtered = self.panel.to_frame()
    #    expected = self.panel.to_frame().dropna(how='any')
    #    assert_frame_equal(filtered, expected)

    #    # unfiltered
    #    unfiltered = self.panel.to_frame(filter_observations=False)
    #    assert_panel_equal(unfiltered.to_panel(), self.panel)

    #    # names
    #    self.assertEqual(unfiltered.index.names, ['major', 'minor'])

    def test_to_frame_mixed(self):
        raise nose.SkipTest
    #    panel = self.panel.fillna(0)
    #    panel['str'] = 'foo'
    #    panel['bool'] = panel['ItemA'] > 0

    #    lp = panel.to_frame()
    #    wp = lp.to_panel()
    #    self.assertEqual(wp['bool'].values.dtype, np.bool_)
    #    assert_frame_equal(wp['bool'], panel['bool'])

    def test_filter(self):
        pass

    def test_apply(self):
        pass

    def test_compound(self):
        raise nose.SkipTest
    #    compounded = self.panel.compound()

    #    assert_series_equal(compounded['ItemA'],
    #                        (1 + self.panel['ItemA']).product(0) - 1)

    def test_shift(self):
        raise nose.SkipTest
    #    # major
    #    idx = self.panel.major_axis[0]
    #    idx_lag = self.panel.major_axis[1]

    #    shifted = self.panel.shift(1)

    #    assert_frame_equal(self.panel.major_xs(idx),
    #                       shifted.major_xs(idx_lag))

    #    # minor
    #    idx = self.panel.minor_axis[0]
    #    idx_lag = self.panel.minor_axis[1]

    #    shifted = self.panel.shift(1, axis='minor')

    #    assert_frame_equal(self.panel.minor_xs(idx),
    #                       shifted.minor_xs(idx_lag))

    #    self.assertRaises(Exception, self.panel.shift, 1, axis='items')

    def test_multiindex_get(self):
        raise nose.SkipTest
    #    ind = MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1), ('b',2)],
    #                                 names=['first', 'second'])
    #    wp = Panel(np.random.random((4,5,5)),
    #                                items=ind,
    #                                major_axis=np.arange(5),
    #                                minor_axis=np.arange(5))
    #    f1 = wp['a']
    #    f2 = wp.ix['a']
    #    assert_panel_equal(f1, f2)

    #    self.assert_((f1.items == [1, 2]).all())
    #    self.assert_((f2.items == [1, 2]).all())

    #    ind = MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)],
    #                                 names=['first', 'second'])

    def test_multiindex_blocks(self):
        raise nose.SkipTest
    #    ind = MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)],
    #                                 names=['first', 'second'])
    #    wp = Panel(self.panel._data)
    #    wp.items = ind
    #    f1 = wp['a']
    #    self.assert_((f1.items == [1, 2]).all())

    #    f1 = wp[('b',1)]
    #    self.assert_((f1.columns == ['A', 'B', 'C', 'D']).all())

    def test_repr_empty(self):
        empty = FDPanel()
        repr(empty)

    def test_rename(self):
        mapper = {
            'l1' : 'foo',
            'l2' : 'bar',
            'l3' : 'baz'
        }

        renamed = self.fdpanel.rename_axis(mapper, axis=0)
        exp = Index(['foo', 'bar', 'baz'])
        self.assert_(renamed.labels.equals(exp))

        renamed = self.fdpanel.rename_axis(str.lower, axis=3)
        exp = Index(['a', 'b', 'c', 'd'])
        self.assert_(renamed.minor_axis.equals(exp))

        # don't copy
        renamed_nocopy = self.fdpanel.rename_axis(mapper, axis=0, copy=False)
        renamed_nocopy['foo'] = 3.
        self.assert_((self.fdpanel['l1'].values == 3).all())

    def test_get_attr(self):
        assert_panel_equal(self.fdpanel['l1'], self.fdpanel.l1)

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
        raise nose.SkipTest
    #    tuples = [('MSFT', 3), ('MSFT', 2), ('AAPL', 2),
    #              ('AAPL', 1), ('MSFT', 1)]
    #    midx = MultiIndex.from_tuples(tuples)
    #    df = DataFrame(np.random.rand(5,4), index=midx)
    #    p = df.to_panel()
    #    assert_frame_equal(p.minor_xs(2), df.ix[:,2].sort_index())

    def test_to_excel(self):
        raise nose.SkipTest
    #    try:
    #        import xlwt
    #        import xlrd
    #        import openpyxl
    #    except ImportError:
    #        raise nose.SkipTest

    #    for ext in ['xls', 'xlsx']:
    #        path = '__tmp__.' + ext
    #        self.panel.to_excel(path)
    #        reader = ExcelFile(path)
    #        for item, df in self.panel.iteritems():
    #            recdf = reader.parse(str(item),index_col=0)
    #            assert_frame_equal(df, recdf)
    #        os.remove(path)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
