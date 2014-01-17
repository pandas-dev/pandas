from datetime import datetime
from pandas.compat import range, lrange
import os
import operator
import nose

import numpy as np

from pandas import Series, DataFrame, Index, isnull, notnull, pivot, MultiIndex
from pandas.core.datetools import bday
from pandas.core.frame import group_agg
from pandas.core.panel import Panel
from pandas.core.panel4d import Panel4D
from pandas.core.series import remove_na
import pandas.core.common as com
import pandas.core.panel as panelmod
from pandas import compat

from pandas.util.testing import (assert_panel_equal,
                                 assert_panel4d_equal,
                                 assert_frame_equal,
                                 assert_series_equal,
                                 assert_almost_equal)
import pandas.util.testing as tm
import pandas.compat as compat


def add_nans(panel4d):
    for l, label in enumerate(panel4d.labels):
        panel = panel4d[label]
        tm.add_nans(panel)


class SafeForLongAndSparse(object):

    _multiprocess_can_split_ = True

    def test_repr(self):
        foo = repr(self.panel4d)

    def test_iter(self):
        tm.equalContents(list(self.panel4d), self.panel4d.labels)

    def test_count(self):
        f = lambda s: notnull(s).sum()
        self._check_stat_op('count', f, obj=self.panel4d, has_skipna=False)

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
        try:
            from scipy.stats import skew
        except ImportError:
            raise nose.SkipTest("no scipy.stats.skew")

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
            obj = self.panel4d

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

    _multiprocess_can_split_ = True

    @classmethod
    def assert_panel_equal(cls, x, y):
        assert_panel_equal(x, y)

    @classmethod
    def assert_panel4d_equal(cls, x, y):
        assert_panel4d_equal(x, y)

    def test_get_axis(self):
        assert(self.panel4d._get_axis(0) is self.panel4d.labels)
        assert(self.panel4d._get_axis(1) is self.panel4d.items)
        assert(self.panel4d._get_axis(2) is self.panel4d.major_axis)
        assert(self.panel4d._get_axis(3) is self.panel4d.minor_axis)

    def test_set_axis(self):
        new_labels = Index(np.arange(len(self.panel4d.labels)))
        new_items = Index(np.arange(len(self.panel4d.items)))
        new_major = Index(np.arange(len(self.panel4d.major_axis)))
        new_minor = Index(np.arange(len(self.panel4d.minor_axis)))

        # ensure propagate to potentially prior-cached items too
        label = self.panel4d['l1']
        self.panel4d.labels = new_labels

        if hasattr(self.panel4d, '_item_cache'):
            self.assert_('l1' not in self.panel4d._item_cache)
        self.assert_(self.panel4d.labels is new_labels)

        self.panel4d.major_axis = new_major
        self.assert_(self.panel4d[0].major_axis is new_major)
        self.assert_(self.panel4d.major_axis is new_major)

        self.panel4d.minor_axis = new_minor
        self.assert_(self.panel4d[0].minor_axis is new_minor)
        self.assert_(self.panel4d.minor_axis is new_minor)

    def test_get_axis_number(self):
        self.assertEqual(self.panel4d._get_axis_number('labels'), 0)
        self.assertEqual(self.panel4d._get_axis_number('items'), 1)
        self.assertEqual(self.panel4d._get_axis_number('major'), 2)
        self.assertEqual(self.panel4d._get_axis_number('minor'), 3)

    def test_get_axis_name(self):
        self.assertEqual(self.panel4d._get_axis_name(0), 'labels')
        self.assertEqual(self.panel4d._get_axis_name(1), 'items')
        self.assertEqual(self.panel4d._get_axis_name(2), 'major_axis')
        self.assertEqual(self.panel4d._get_axis_name(3), 'minor_axis')

    def test_arith(self):
        self._test_op(self.panel4d, operator.add)
        self._test_op(self.panel4d, operator.sub)
        self._test_op(self.panel4d, operator.mul)
        self._test_op(self.panel4d, operator.truediv)
        self._test_op(self.panel4d, operator.floordiv)
        self._test_op(self.panel4d, operator.pow)

        self._test_op(self.panel4d, lambda x, y: y + x)
        self._test_op(self.panel4d, lambda x, y: y - x)
        self._test_op(self.panel4d, lambda x, y: y * x)
        self._test_op(self.panel4d, lambda x, y: y / x)
        self._test_op(self.panel4d, lambda x, y: y ** x)

        self.assertRaises(Exception, self.panel4d.__add__, self.panel4d['l1'])

    @staticmethod
    def _test_op(panel4d, op):
        result = op(panel4d, 1)
        assert_panel_equal(result['l1'], op(panel4d['l1'], 1))

    def test_keys(self):
        tm.equalContents(list(self.panel4d.keys()), self.panel4d.labels)

    def test_iteritems(self):
        """Test panel4d.iteritems()"""

        self.assertEqual(len(list(compat.iteritems(self.panel4d))),
                         len(self.panel4d.labels))

    def test_combinePanel4d(self):
        result = self.panel4d.add(self.panel4d)
        self.assert_panel4d_equal(result, self.panel4d * 2)

    def test_neg(self):
        self.assert_panel4d_equal(-self.panel4d, self.panel4d * -1)

    def test_select(self):
        p = self.panel4d

        # select labels
        result = p.select(lambda x: x in ('l1', 'l3'), axis='labels')
        expected = p.reindex(labels=['l1', 'l3'])
        self.assert_panel4d_equal(result, expected)

        # select items
        result = p.select(lambda x: x in ('ItemA', 'ItemC'), axis='items')
        expected = p.reindex(items=['ItemA', 'ItemC'])
        self.assert_panel4d_equal(result, expected)

        # select major_axis
        result = p.select(lambda x: x >= datetime(2000, 1, 15), axis='major')
        new_major = p.major_axis[p.major_axis >= datetime(2000, 1, 15)]
        expected = p.reindex(major=new_major)
        self.assert_panel4d_equal(result, expected)

        # select minor_axis
        result = p.select(lambda x: x in ('D', 'A'), axis=3)
        expected = p.reindex(minor=['A', 'D'])
        self.assert_panel4d_equal(result, expected)

        # corner case, empty thing
        result = p.select(lambda x: x in ('foo',), axis='items')
        self.assert_panel4d_equal(result, p.reindex(items=[]))

    def test_get_value(self):
        for item in self.panel.items:
            for mjr in self.panel.major_axis[::2]:
                for mnr in self.panel.minor_axis:
                    result = self.panel.get_value(item, mjr, mnr)
                    expected = self.panel[item][mnr][mjr]
                    assert_almost_equal(result, expected)

    def test_abs(self):
        result = self.panel4d.abs()
        expected = np.abs(self.panel4d)
        self.assert_panel4d_equal(result, expected)

        p = self.panel4d['l1']
        result = p.abs()
        expected = np.abs(p)
        assert_panel_equal(result, expected)

        df = p['ItemA']
        result = df.abs()
        expected = np.abs(df)
        assert_frame_equal(result, expected)


class CheckIndexing(object):

    _multiprocess_can_split_ = True

    def test_getitem(self):
        self.assertRaises(Exception, self.panel4d.__getitem__, 'ItemQ')

    def test_delitem_and_pop(self):
        expected = self.panel4d['l2']
        result = self.panel4d.pop('l2')
        assert_panel_equal(expected, result)
        self.assert_('l2' not in self.panel4d.labels)

        del self.panel4d['l3']
        self.assert_('l3' not in self.panel4d.labels)
        self.assertRaises(Exception, self.panel4d.__delitem__, 'l3')

        values = np.empty((4, 4, 4, 4))
        values[0] = 0
        values[1] = 1
        values[2] = 2
        values[3] = 3

        panel4d = Panel4D(values, lrange(4), lrange(4), lrange(4), lrange(4))

        # did we delete the right row?

        panel4dc = panel4d.copy()
        del panel4dc[0]
        assert_panel_equal(panel4dc[1], panel4d[1])
        assert_panel_equal(panel4dc[2], panel4d[2])
        assert_panel_equal(panel4dc[3], panel4d[3])

        panel4dc = panel4d.copy()
        del panel4dc[1]
        assert_panel_equal(panel4dc[0], panel4d[0])
        assert_panel_equal(panel4dc[2], panel4d[2])
        assert_panel_equal(panel4dc[3], panel4d[3])

        panel4dc = panel4d.copy()
        del panel4dc[2]
        assert_panel_equal(panel4dc[1], panel4d[1])
        assert_panel_equal(panel4dc[0], panel4d[0])
        assert_panel_equal(panel4dc[3], panel4d[3])

        panel4dc = panel4d.copy()
        del panel4dc[3]
        assert_panel_equal(panel4dc[1], panel4d[1])
        assert_panel_equal(panel4dc[2], panel4d[2])
        assert_panel_equal(panel4dc[0], panel4d[0])

    def test_setitem(self):
        ## LongPanel with one item
        # lp = self.panel.filter(['ItemA', 'ItemB']).to_frame()
        # self.assertRaises(Exception, self.panel.__setitem__,
        #                  'ItemE', lp)

        # Panel
        p = Panel(dict(
            ItemA=self.panel4d['l1']['ItemA'][2:].filter(items=['A', 'B'])))
        self.panel4d['l4'] = p
        self.panel4d['l5'] = p

        p2 = self.panel4d['l4']

        assert_panel_equal(p, p2.reindex(items=p.items,
                                         major_axis=p.major_axis,
                                         minor_axis=p.minor_axis))

        # scalar
        self.panel4d['lG'] = 1
        self.panel4d['lE'] = True
        self.assert_(self.panel4d['lG'].values.dtype == np.int64)
        self.assert_(self.panel4d['lE'].values.dtype == np.bool_)

        # object dtype
        self.panel4d['lQ'] = 'foo'
        self.assert_(self.panel4d['lQ'].values.dtype == np.object_)

        # boolean dtype
        self.panel4d['lP'] = self.panel4d['l1'] > 0
        self.assert_(self.panel4d['lP'].values.dtype == np.bool_)

    def test_comparisons(self):
        p1 = tm.makePanel4D()
        p2 = tm.makePanel4D()

        tp = p1.reindex(labels=p1.labels + ['foo'])
        p = p1[p1.labels[0]]

        def test_comp(func):
            result = func(p1, p2)
            self.assert_(np.array_equal(result.values,
                                        func(p1.values, p2.values)))

            # versus non-indexed same objs
            self.assertRaises(Exception, func, p1, tp)

            # versus different objs
            self.assertRaises(Exception, func, p1, p)

            result3 = func(self.panel4d, 0)
            self.assert_(np.array_equal(result3.values,
                                        func(self.panel4d.values, 0)))

        test_comp(operator.eq)
        test_comp(operator.ne)
        test_comp(operator.lt)
        test_comp(operator.gt)
        test_comp(operator.ge)
        test_comp(operator.le)

    def test_setitem_ndarray(self):
        raise nose.SkipTest("skipping for now")
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
        ref = self.panel4d['l1']['ItemA']

        idx = self.panel4d.major_axis[5]
        xs = self.panel4d.major_xs(idx)

        assert_series_equal(xs['l1'].T['ItemA'], ref.xs(idx))

        # not contained
        idx = self.panel4d.major_axis[0] - bday
        self.assertRaises(Exception, self.panel4d.major_xs, idx)

    def test_major_xs_mixed(self):
        self.panel4d['l4'] = 'foo'
        xs = self.panel4d.major_xs(self.panel4d.major_axis[0])
        self.assert_(xs['l1']['A'].dtype == np.float64)
        self.assert_(xs['l4']['A'].dtype == np.object_)

    def test_minor_xs(self):
        ref = self.panel4d['l1']['ItemA']

        idx = self.panel4d.minor_axis[1]
        xs = self.panel4d.minor_xs(idx)

        assert_series_equal(xs['l1'].T['ItemA'], ref[idx])

        # not contained
        self.assertRaises(Exception, self.panel4d.minor_xs, 'E')

    def test_minor_xs_mixed(self):
        self.panel4d['l4'] = 'foo'

        xs = self.panel4d.minor_xs('D')
        self.assert_(xs['l1'].T['ItemA'].dtype == np.float64)
        self.assert_(xs['l4'].T['ItemA'].dtype == np.object_)

    def test_xs(self):
        l1 = self.panel4d.xs('l1', axis=0)
        expected = self.panel4d['l1']
        assert_panel_equal(l1, expected)

        # not view by default
        l1.values[:] = np.nan
        self.assert_(not np.isnan(self.panel4d['l1'].values).all())

        # but can get view
        l1_view = self.panel4d.xs('l1', axis=0, copy=False)
        l1_view.values[:] = np.nan
        self.assert_(np.isnan(self.panel4d['l1'].values).all())

        # mixed-type
        self.panel4d['strings'] = 'foo'
        self.assertRaises(Exception, self.panel4d.xs, 'D', axis=2,
                          copy=False)

    def test_getitem_fancy_labels(self):
        panel4d = self.panel4d

        labels = panel4d.labels[[1, 0]]
        items = panel4d.items[[1, 0]]
        dates = panel4d.major_axis[::2]
        cols = ['D', 'C', 'F']

        # all 4 specified
        assert_panel4d_equal(panel4d.ix[labels, items, dates, cols],
                             panel4d.reindex(labels=labels, items=items, major=dates, minor=cols))

        # 3 specified
        assert_panel4d_equal(panel4d.ix[:, items, dates, cols],
                             panel4d.reindex(items=items, major=dates, minor=cols))

        # 2 specified
        assert_panel4d_equal(panel4d.ix[:, :, dates, cols],
                             panel4d.reindex(major=dates, minor=cols))

        assert_panel4d_equal(panel4d.ix[:, items, :, cols],
                             panel4d.reindex(items=items, minor=cols))

        assert_panel4d_equal(panel4d.ix[:, items, dates, :],
                             panel4d.reindex(items=items, major=dates))

        # only 1
        assert_panel4d_equal(panel4d.ix[:, items, :, :],
                             panel4d.reindex(items=items))

        assert_panel4d_equal(panel4d.ix[:, :, dates, :],
                             panel4d.reindex(major=dates))

        assert_panel4d_equal(panel4d.ix[:, :, :, cols],
                             panel4d.reindex(minor=cols))

    def test_getitem_fancy_slice(self):
        pass

    def test_getitem_fancy_ints(self):
        pass

    def test_getitem_fancy_xs(self):
        raise nose.SkipTest("skipping for now")
        # self.assertRaises(NotImplementedError, self.panel4d.major_xs)
        # self.assertRaises(NotImplementedError, self.panel4d.minor_xs)

    def test_get_value(self):
        for label in self.panel4d.labels:
            for item in self.panel4d.items:
                for mjr in self.panel4d.major_axis[::2]:
                    for mnr in self.panel4d.minor_axis:
                        result = self.panel4d.get_value(
                            label, item, mjr, mnr)
                        expected = self.panel4d[label][item][mnr][mjr]
                        assert_almost_equal(result, expected)

    def test_set_value(self):
        for label in self.panel4d.labels:
            for item in self.panel4d.items:
                for mjr in self.panel4d.major_axis[::2]:
                    for mnr in self.panel4d.minor_axis:
                        self.panel4d.set_value(label, item, mjr, mnr, 1.)
                        assert_almost_equal(
                            self.panel4d[label][item][mnr][mjr], 1.)

        # resize
        res = self.panel4d.set_value('l4', 'ItemE', 'foo', 'bar', 1.5)
        tm.assert_isinstance(res, Panel4D)
        self.assert_(res is not self.panel4d)
        self.assertEqual(res.get_value('l4', 'ItemE', 'foo', 'bar'), 1.5)

        res3 = self.panel4d.set_value('l4', 'ItemE', 'foobar', 'baz', 5)
        self.assert_(com.is_float_dtype(res3['l4'].values))


class TestPanel4d(tm.TestCase, CheckIndexing, SafeForSparse,
                  SafeForLongAndSparse):

    _multiprocess_can_split_ = True

    @classmethod
    def assert_panel4d_equal(cls, x, y):
        assert_panel4d_equal(x, y)

    def setUp(self):
        self.panel4d = tm.makePanel4D(nper=8)
        add_nans(self.panel4d)

    def test_constructor(self):
        # with BlockManager
        panel4d = Panel4D(self.panel4d._data)
        self.assert_(panel4d._data is self.panel4d._data)

        panel4d = Panel4D(self.panel4d._data, copy=True)
        self.assert_(panel4d._data is not self.panel4d._data)
        assert_panel4d_equal(panel4d, self.panel4d)

        # strings handled prop
        # panel4d = Panel4D([[['foo', 'foo', 'foo',],
        #                 ['foo', 'foo', 'foo']]])
        # self.assert_(wp.values.dtype == np.object_)

        vals = self.panel4d.values

        # no copy
        panel4d = Panel4D(vals)
        self.assert_(panel4d.values is vals)

        # copy
        panel4d = Panel4D(vals, copy=True)
        self.assert_(panel4d.values is not vals)

    def test_constructor_cast(self):
        zero_filled = self.panel4d.fillna(0)

        casted = Panel4D(zero_filled._data, dtype=int)
        casted2 = Panel4D(zero_filled.values, dtype=int)

        exp_values = zero_filled.values.astype(int)
        assert_almost_equal(casted.values, exp_values)
        assert_almost_equal(casted2.values, exp_values)

        casted = Panel4D(zero_filled._data, dtype=np.int32)
        casted2 = Panel4D(zero_filled.values, dtype=np.int32)

        exp_values = zero_filled.values.astype(np.int32)
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
        panel = Panel(items=lrange(3), major_axis=lrange(3),
                      minor_axis=lrange(3), dtype='O')
        self.assert_(panel.values.dtype == np.object_)

    def test_consolidate(self):
        self.assert_(self.panel4d._data.is_consolidated())

        self.panel4d['foo'] = 1.
        self.assert_(not self.panel4d._data.is_consolidated())

        panel4d = self.panel4d.consolidate()
        self.assert_(panel4d._data.is_consolidated())

    def test_ctor_dict(self):
        l1 = self.panel4d['l1']
        l2 = self.panel4d['l2']

        d = {'A': l1, 'B': l2.ix[['ItemB'], :, :]}
        # d2 = {'A' : itema._series, 'B' : itemb[5:]._series}
        # d3 = {'A' : DataFrame(itema._series),
        #      'B' : DataFrame(itemb[5:]._series)}

        panel4d = Panel4D(d)
        # wp2 = Panel.from_dict(d2) # nested Dict
        # wp3 = Panel.from_dict(d3)
        # self.assert_(wp.major_axis.equals(self.panel.major_axis))
        assert_panel_equal(panel4d['A'], self.panel4d['l1'])
        assert_frame_equal(panel4d.ix['B', 'ItemB', :, :],
                           self.panel4d.ix['l2', ['ItemB'], :, :]['ItemB'])

        # intersect
        # wp = Panel.from_dict(d, intersect=True)
        # self.assert_(wp.major_axis.equals(itemb.index[5:]))

        # use constructor
        # assert_panel_equal(Panel(d), Panel.from_dict(d))
        # assert_panel_equal(Panel(d2), Panel.from_dict(d2))
        # assert_panel_equal(Panel(d3), Panel.from_dict(d3))

        # cast
        # dcasted = dict((k, v.reindex(wp.major_axis).fillna(0))
        #               for k, v in d.iteritems())
        # result = Panel(dcasted, dtype=int)
        # expected = Panel(dict((k, v.astype(int))
        #                      for k, v in dcasted.iteritems()))
        # assert_panel_equal(result, expected)

    def test_constructor_dict_mixed(self):
        data = dict((k, v.values) for k, v in compat.iteritems(self.panel4d))
        result = Panel4D(data)
        exp_major = Index(np.arange(len(self.panel4d.major_axis)))
        self.assert_(result.major_axis.equals(exp_major))

        result = Panel4D(data,
                         labels=self.panel4d.labels,
                         items=self.panel4d.items,
                         major_axis=self.panel4d.major_axis,
                         minor_axis=self.panel4d.minor_axis)
        assert_panel4d_equal(result, self.panel4d)

        data['l2'] = self.panel4d['l2']
        result = Panel4D(data)
        assert_panel4d_equal(result, self.panel4d)

        # corner, blow up
        data['l2'] = data['l2']['ItemB']
        self.assertRaises(Exception, Panel4D, data)

        data['l2'] = self.panel4d['l2'].values[:, :, :-1]
        self.assertRaises(Exception, Panel4D, data)

    def test_constructor_resize(self):
        data = self.panel4d._data
        labels = self.panel4d.labels[:-1]
        items = self.panel4d.items[:-1]
        major = self.panel4d.major_axis[:-1]
        minor = self.panel4d.minor_axis[:-1]

        result = Panel4D(data, labels=labels, items=items,
                         major_axis=major, minor_axis=minor)
        expected = self.panel4d.reindex(
            labels=labels, items=items, major=major, minor=minor)
        assert_panel4d_equal(result, expected)

        result = Panel4D(data, items=items, major_axis=major)
        expected = self.panel4d.reindex(items=items, major=major)
        assert_panel4d_equal(result, expected)

        result = Panel4D(data, items=items)
        expected = self.panel4d.reindex(items=items)
        assert_panel4d_equal(result, expected)

        result = Panel4D(data, minor_axis=minor)
        expected = self.panel4d.reindex(minor=minor)
        assert_panel4d_equal(result, expected)

    def test_from_dict_mixed_orient(self):
        raise nose.SkipTest("skipping for now")
    #    df = tm.makeDataFrame()
    #    df['foo'] = 'bar'

    #    data = {'k1' : df,
    #            'k2' : df}

    #    panel = Panel.from_dict(data, orient='minor')

    #    self.assert_(panel['foo'].values.dtype == np.object_)
    #    self.assert_(panel['A'].values.dtype == np.float64)

    def test_values(self):
        self.assertRaises(Exception, Panel, np.random.randn(5, 5, 5),
                          lrange(5), lrange(5), lrange(4))

    def test_conform(self):
        p = self.panel4d['l1'].filter(items=['ItemA', 'ItemB'])
        conformed = self.panel4d.conform(p)

        assert(conformed.items.equals(self.panel4d.labels))
        assert(conformed.major_axis.equals(self.panel4d.major_axis))
        assert(conformed.minor_axis.equals(self.panel4d.minor_axis))

    def test_reindex(self):
        ref = self.panel4d['l2']

        # labels
        result = self.panel4d.reindex(labels=['l1', 'l2'])
        assert_panel_equal(result['l2'], ref)

        # items
        result = self.panel4d.reindex(items=['ItemA', 'ItemB'])
        assert_frame_equal(result['l2']['ItemB'], ref['ItemB'])

        # major
        new_major = list(self.panel4d.major_axis[:10])
        result = self.panel4d.reindex(major=new_major)
        assert_frame_equal(
            result['l2']['ItemB'], ref['ItemB'].reindex(index=new_major))

        # raise exception put both major and major_axis
        self.assertRaises(Exception, self.panel4d.reindex,
                          major_axis=new_major, major=new_major)

        # minor
        new_minor = list(self.panel4d.minor_axis[:2])
        result = self.panel4d.reindex(minor=new_minor)
        assert_frame_equal(
            result['l2']['ItemB'], ref['ItemB'].reindex(columns=new_minor))

        result = self.panel4d.reindex(labels=self.panel4d.labels,
                                      items=self.panel4d.items,
                                      major=self.panel4d.major_axis,
                                      minor=self.panel4d.minor_axis)

        # don't necessarily copy
        result = self.panel4d.reindex()
        assert_panel4d_equal(result,self.panel4d)
        self.assert_((result is self.panel4d) == False)

        # with filling
        smaller_major = self.panel4d.major_axis[::5]
        smaller = self.panel4d.reindex(major=smaller_major)

        larger = smaller.reindex(major=self.panel4d.major_axis,
                                 method='pad')

        assert_panel_equal(larger.ix[:, :, self.panel4d.major_axis[1], :],
                           smaller.ix[:, :, smaller_major[0], :])

        # don't necessarily copy
        result = self.panel4d.reindex(
            major=self.panel4d.major_axis, copy=False)
        assert_panel4d_equal(result,self.panel4d)
        self.assert_((result is self.panel4d) == True)

    def test_not_hashable(self):
        p4D_empty = Panel4D()
        self.assertRaises(TypeError, hash, p4D_empty)
        self.assertRaises(TypeError, hash, self.panel4d)

    def test_reindex_like(self):
        # reindex_like
        smaller = self.panel4d.reindex(labels=self.panel4d.labels[:-1],
                                       items=self.panel4d.items[:-1],
                                       major=self.panel4d.major_axis[:-1],
                                       minor=self.panel4d.minor_axis[:-1])
        smaller_like = self.panel4d.reindex_like(smaller)
        assert_panel4d_equal(smaller, smaller_like)

    def test_take(self):
        raise nose.SkipTest("skipping for now")

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

        rlabels = list(self.panel4d.labels)
        ritems = list(self.panel4d.items)
        rmajor = list(self.panel4d.major_axis)
        rminor = list(self.panel4d.minor_axis)
        random.shuffle(rlabels)
        random.shuffle(ritems)
        random.shuffle(rmajor)
        random.shuffle(rminor)

        random_order = self.panel4d.reindex(labels=rlabels)
        sorted_panel4d = random_order.sort_index(axis=0)
        assert_panel4d_equal(sorted_panel4d, self.panel4d)

        # descending
        # random_order = self.panel.reindex(items=ritems)
        # sorted_panel = random_order.sort_index(axis=0, ascending=False)
        # assert_panel_equal(sorted_panel,
        #                   self.panel.reindex(items=self.panel.items[::-1]))

        # random_order = self.panel.reindex(major=rmajor)
        # sorted_panel = random_order.sort_index(axis=1)
        # assert_panel_equal(sorted_panel, self.panel)

        # random_order = self.panel.reindex(minor=rminor)
        # sorted_panel = random_order.sort_index(axis=2)
        # assert_panel_equal(sorted_panel, self.panel)

    def test_fillna(self):
        self.assert_(not np.isfinite(self.panel4d.values).all())
        filled = self.panel4d.fillna(0)
        self.assert_(np.isfinite(filled.values).all())

        self.assertRaises(NotImplementedError, self.panel4d.fillna, method='pad')

    def test_swapaxes(self):
        result = self.panel4d.swapaxes('labels', 'items')
        self.assert_(result.items is self.panel4d.labels)

        result = self.panel4d.swapaxes('labels', 'minor')
        self.assert_(result.labels is self.panel4d.minor_axis)

        result = self.panel4d.swapaxes('items', 'minor')
        self.assert_(result.items is self.panel4d.minor_axis)

        result = self.panel4d.swapaxes('items', 'major')
        self.assert_(result.items is self.panel4d.major_axis)

        result = self.panel4d.swapaxes('major', 'minor')
        self.assert_(result.major_axis is self.panel4d.minor_axis)

        # this should also work
        result = self.panel4d.swapaxes(0, 1)
        self.assert_(result.labels is self.panel4d.items)

        # this works, but return a copy
        result = self.panel4d.swapaxes('items', 'items')
        assert_panel4d_equal(self.panel4d,result)
        self.assert_(id(self.panel4d) != id(result))

    def test_to_frame(self):
        raise nose.SkipTest("skipping for now")
    #    # filtered
    #    filtered = self.panel.to_frame()
    #    expected = self.panel.to_frame().dropna(how='any')
    #    assert_frame_equal(filtered, expected)

    #    # unfiltered
    #    unfiltered = self.panel.to_frame(filter_observations=False)
    #    assert_panel_equal(unfiltered.to_panel(), self.panel)

    #    # names
    #    self.assertEqual(unfiltered.index.names, ('major', 'minor'))

    def test_to_frame_mixed(self):
        raise nose.SkipTest("skipping for now")
    #    panel = self.panel.fillna(0)
    #    panel['str'] = 'foo'
    #    panel['bool'] = panel['ItemA'] > 0

    #    lp = panel.to_frame()
    #    wp = lp.to_panel()
    #    self.assertEqual(wp['bool'].values.dtype, np.bool_)
    #    assert_frame_equal(wp['bool'], panel['bool'])

    def test_update(self):

        p4d = Panel4D([[[[1.5, np.nan, 3.],
                         [1.5, np.nan, 3.],
                         [1.5, np.nan, 3.],
                         [1.5, np.nan, 3.]],
                        [[1.5, np.nan, 3.],
                         [1.5, np.nan, 3.],
                         [1.5, np.nan, 3.],
                         [1.5, np.nan, 3.]]]])

        other = Panel4D([[[[3.6, 2., np.nan]],
                          [[np.nan, np.nan, 7]]]])

        p4d.update(other)

        expected = Panel4D([[[[3.6, 2, 3.],
                              [1.5, np.nan, 3.],
                              [1.5, np.nan, 3.],
                              [1.5, np.nan, 3.]],
                             [[1.5, np.nan, 7],
                              [1.5, np.nan, 3.],
                              [1.5, np.nan, 3.],
                              [1.5, np.nan, 3.]]]])

        assert_panel4d_equal(p4d, expected)

    def test_filter(self):
        raise nose.SkipTest("skipping for now")

    def test_apply(self):
        raise nose.SkipTest("skipping for now")

    def test_dtypes(self):

        result = self.panel4d.dtypes
        expected = Series(np.dtype('float64'),index=self.panel4d.labels)
        assert_series_equal(result, expected)

    def test_compound(self):
        raise nose.SkipTest("skipping for now")
    #    compounded = self.panel.compound()

    #    assert_series_equal(compounded['ItemA'],
    #                        (1 + self.panel['ItemA']).product(0) - 1)

    def test_shift(self):
        raise nose.SkipTest("skipping for now")
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
        raise nose.SkipTest("skipping for now")
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
        raise nose.SkipTest("skipping for now")
    #    ind = MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)],
    #                                 names=['first', 'second'])
    #    wp = Panel(self.panel._data)
    #    wp.items = ind
    #    f1 = wp['a']
    #    self.assert_((f1.items == [1, 2]).all())

    #    f1 = wp[('b',1)]
    #    self.assert_((f1.columns == ['A', 'B', 'C', 'D']).all())

    def test_repr_empty(self):
        empty = Panel4D()
        repr(empty)

    def test_rename(self):
        mapper = {
            'l1': 'foo',
            'l2': 'bar',
            'l3': 'baz'
        }

        renamed = self.panel4d.rename_axis(mapper, axis=0)
        exp = Index(['foo', 'bar', 'baz'])
        self.assert_(renamed.labels.equals(exp))

        renamed = self.panel4d.rename_axis(str.lower, axis=3)
        exp = Index(['a', 'b', 'c', 'd'])
        self.assert_(renamed.minor_axis.equals(exp))

        # don't copy
        renamed_nocopy = self.panel4d.rename_axis(mapper, axis=0, copy=False)
        renamed_nocopy['foo'] = 3.
        self.assert_((self.panel4d['l1'].values == 3).all())

    def test_get_attr(self):
        assert_panel_equal(self.panel4d['l1'], self.panel4d.l1)

    def test_group_agg(self):
        values = np.ones((10, 2)) * np.arange(10).reshape((10, 1))
        bounds = np.arange(5) * 2
        f = lambda x: x.mean(axis=0)

        agged = group_agg(values, bounds, f)

        assert(agged[1][0] == 2.5)
        assert(agged[2][0] == 4.5)

        # test a function that doesn't aggregate
        f2 = lambda x: np.zeros((2, 2))
        self.assertRaises(Exception, group_agg, values, bounds, f2)

    def test_from_frame_level1_unsorted(self):
        raise nose.SkipTest("skipping for now")

    def test_to_excel(self):
        raise nose.SkipTest("skipping for now")


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure',
                         '--with-timer'],
                   exit=False)
