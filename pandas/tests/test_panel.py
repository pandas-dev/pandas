# -*- coding: utf-8 -*-
# pylint: disable=W0612,E1101
from collections import OrderedDict
from datetime import datetime

import numpy as np
import pytest

from pandas.compat import lrange

from pandas import DataFrame, MultiIndex, Series, date_range, notna
import pandas.core.panel as panelm
from pandas.core.panel import Panel
import pandas.util.testing as tm
from pandas.util.testing import (
    assert_almost_equal, assert_frame_equal, assert_series_equal,
    makeCustomDataframe as mkdf, makeMixedDataFrame)

from pandas.tseries.offsets import MonthEnd


@pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
class PanelTests(object):
    panel = None

    def not_hashable(self):
        c_empty = Panel()
        c = Panel(Panel([[[1]]]))
        pytest.raises(TypeError, hash, c_empty)
        pytest.raises(TypeError, hash, c)


@pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
class SafeForSparse(object):

    # issue 7692
    def test_raise_when_not_implemented(self):
        p = Panel(np.arange(3 * 4 * 5).reshape(3, 4, 5),
                  items=['ItemA', 'ItemB', 'ItemC'],
                  major_axis=date_range('20130101', periods=4),
                  minor_axis=list('ABCDE'))
        d = p.sum(axis=1).iloc[0]
        ops = ['add', 'sub', 'mul', 'truediv',
               'floordiv', 'div', 'mod', 'pow']
        for op in ops:
            with pytest.raises(NotImplementedError):
                getattr(p, op)(d, axis=0)


@pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
class CheckIndexing(object):

    def test_delitem_and_pop(self):

        values = np.empty((3, 3, 3))
        values[0] = 0
        values[1] = 1
        values[2] = 2

        panel = Panel(values, lrange(3), lrange(3), lrange(3))

        # did we delete the right row?

        panelc = panel.copy()
        del panelc[0]
        tm.assert_frame_equal(panelc[1], panel[1])
        tm.assert_frame_equal(panelc[2], panel[2])

        panelc = panel.copy()
        del panelc[1]
        tm.assert_frame_equal(panelc[0], panel[0])
        tm.assert_frame_equal(panelc[2], panel[2])

        panelc = panel.copy()
        del panelc[2]
        tm.assert_frame_equal(panelc[1], panel[1])
        tm.assert_frame_equal(panelc[0], panel[0])

    def test_setitem(self):
        # bad shape
        p = Panel(np.random.randn(4, 3, 2))
        msg = (r"shape of value must be \(3, 2\), "
               r"shape of given object was \(4, 2\)")
        with pytest.raises(ValueError, match=msg):
            p[0] = np.random.randn(4, 2)

    def test_setitem_ndarray(self):
        timeidx = date_range(start=datetime(2009, 1, 1),
                             end=datetime(2009, 12, 31),
                             freq=MonthEnd())
        lons_coarse = np.linspace(-177.5, 177.5, 72)
        lats_coarse = np.linspace(-87.5, 87.5, 36)
        P = Panel(items=timeidx, major_axis=lons_coarse,
                  minor_axis=lats_coarse)
        data = np.random.randn(72 * 36).reshape((72, 36))
        key = datetime(2009, 2, 28)
        P[key] = data

        assert_almost_equal(P[key].values, data)

    def test_set_minor_major(self):
        # GH 11014
        df1 = DataFrame(['a', 'a', 'a', np.nan, 'a', np.nan])
        df2 = DataFrame([1.0, np.nan, 1.0, np.nan, 1.0, 1.0])
        panel = Panel({'Item1': df1, 'Item2': df2})

        newminor = notna(panel.iloc[:, :, 0])
        panel.loc[:, :, 'NewMinor'] = newminor
        assert_frame_equal(panel.loc[:, :, 'NewMinor'],
                           newminor.astype(object))

        newmajor = notna(panel.iloc[:, 0, :])
        panel.loc[:, 'NewMajor', :] = newmajor
        assert_frame_equal(panel.loc[:, 'NewMajor', :],
                           newmajor.astype(object))

    def test_getitem_fancy_slice(self):
        pass

    def test_ix_setitem_slice_dataframe(self):
        a = Panel(items=[1, 2, 3], major_axis=[11, 22, 33],
                  minor_axis=[111, 222, 333])
        b = DataFrame(np.random.randn(2, 3), index=[111, 333],
                      columns=[1, 2, 3])

        a.loc[:, 22, [111, 333]] = b

        assert_frame_equal(a.loc[:, 22, [111, 333]], b)

    def test_ix_align(self):
        from pandas import Series
        b = Series(np.random.randn(10), name=0)
        b.sort_values()
        df_orig = Panel(np.random.randn(3, 10, 2))
        df = df_orig.copy()

        df.loc[0, :, 0] = b
        assert_series_equal(df.loc[0, :, 0].reindex(b.index), b)

        df = df_orig.swapaxes(0, 1)
        df.loc[:, 0, 0] = b
        assert_series_equal(df.loc[:, 0, 0].reindex(b.index), b)

        df = df_orig.swapaxes(1, 2)
        df.loc[0, 0, :] = b
        assert_series_equal(df.loc[0, 0, :].reindex(b.index), b)

    def test_ix_frame_align(self):
        # GH3830, panel assignent by values/frame
        for dtype in ['float64', 'int64']:

            panel = Panel(np.arange(40).reshape((2, 4, 5)),
                          items=['a1', 'a2'], dtype=dtype)
            df1 = panel.iloc[0]
            df2 = panel.iloc[1]

            tm.assert_frame_equal(panel.loc['a1'], df1)
            tm.assert_frame_equal(panel.loc['a2'], df2)

            # Assignment by Value Passes for 'a2'
            panel.loc['a2'] = df1.values
            tm.assert_frame_equal(panel.loc['a1'], df1)
            tm.assert_frame_equal(panel.loc['a2'], df1)

            # Assignment by DataFrame Ok w/o loc 'a2'
            panel['a2'] = df2
            tm.assert_frame_equal(panel.loc['a1'], df1)
            tm.assert_frame_equal(panel.loc['a2'], df2)

            # Assignment by DataFrame Fails for 'a2'
            panel.loc['a2'] = df2
            tm.assert_frame_equal(panel.loc['a1'], df1)
            tm.assert_frame_equal(panel.loc['a2'], df2)

    def test_logical_with_nas(self):
        d = Panel({'ItemA': {'a': [np.nan, False]},
                   'ItemB': {'a': [True, True]}})

        result = d['ItemA'] | d['ItemB']
        expected = DataFrame({'a': [np.nan, True]})
        assert_frame_equal(result, expected)

        # this is autodowncasted here
        result = d['ItemA'].fillna(False) | d['ItemB']
        expected = DataFrame({'a': [True, True]})
        assert_frame_equal(result, expected)


@pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
class TestPanel(PanelTests, CheckIndexing, SafeForSparse):

    def test_constructor_cast(self):
        # can't cast
        data = [[['foo', 'bar', 'baz']]]
        pytest.raises(ValueError, Panel, data, dtype=float)

    def test_constructor_empty_panel(self):
        empty = Panel()
        assert len(empty.items) == 0
        assert len(empty.major_axis) == 0
        assert len(empty.minor_axis) == 0

    def test_constructor_observe_dtype(self):
        # GH #411
        panel = Panel(items=lrange(3), major_axis=lrange(3),
                      minor_axis=lrange(3), dtype='O')
        assert panel.values.dtype == np.object_

    def test_constructor_dtypes(self):
        # GH #797

        def _check_dtype(panel, dtype):
            for i in panel.items:
                assert panel[i].values.dtype.name == dtype

        # only nan holding types allowed here
        for dtype in ['float64', 'float32', 'object']:
            panel = Panel(items=lrange(2), major_axis=lrange(10),
                          minor_axis=lrange(5), dtype=dtype)
            _check_dtype(panel, dtype)

        for dtype in ['float64', 'float32', 'int64', 'int32', 'object']:
            panel = Panel(np.array(np.random.randn(2, 10, 5), dtype=dtype),
                          items=lrange(2),
                          major_axis=lrange(10),
                          minor_axis=lrange(5), dtype=dtype)
            _check_dtype(panel, dtype)

        for dtype in ['float64', 'float32', 'int64', 'int32', 'object']:
            panel = Panel(np.array(np.random.randn(2, 10, 5), dtype='O'),
                          items=lrange(2),
                          major_axis=lrange(10),
                          minor_axis=lrange(5), dtype=dtype)
            _check_dtype(panel, dtype)

        for dtype in ['float64', 'float32', 'int64', 'int32', 'object']:
            panel = Panel(
                np.random.randn(2, 10, 5),
                items=lrange(2), major_axis=lrange(10),
                minor_axis=lrange(5),
                dtype=dtype)
            _check_dtype(panel, dtype)

        for dtype in ['float64', 'float32', 'int64', 'int32', 'object']:
            df1 = DataFrame(np.random.randn(2, 5),
                            index=lrange(2), columns=lrange(5))
            df2 = DataFrame(np.random.randn(2, 5),
                            index=lrange(2), columns=lrange(5))
            panel = Panel.from_dict({'a': df1, 'b': df2}, dtype=dtype)
            _check_dtype(panel, dtype)

    def test_constructor_fails_with_not_3d_input(self):
        msg = "The number of dimensions required is 3"
        with pytest.raises(ValueError, match=msg):
                Panel(np.random.randn(10, 2))

    def test_ctor_orderedDict(self):
        keys = list(set(np.random.randint(0, 5000, 100)))[
            :50]  # unique random int  keys
        d = OrderedDict([(k, mkdf(10, 5)) for k in keys])
        p = Panel(d)
        assert list(p.items) == keys

        p = Panel.from_dict(d)
        assert list(p.items) == keys

    def test_from_dict_mixed_orient(self):
        df = tm.makeDataFrame()
        df['foo'] = 'bar'

        data = {'k1': df, 'k2': df}

        panel = Panel.from_dict(data, orient='minor')

        assert panel['foo'].values.dtype == np.object_
        assert panel['A'].values.dtype == np.float64

    def test_constructor_error_msgs(self):
        msg = (r"Shape of passed values is \(3, 4, 5\), "
               r"indices imply \(4, 5, 5\)")
        with pytest.raises(ValueError, match=msg):
            Panel(np.random.randn(3, 4, 5),
                  lrange(4), lrange(5), lrange(5))

        msg = (r"Shape of passed values is \(3, 4, 5\), "
               r"indices imply \(5, 4, 5\)")
        with pytest.raises(ValueError, match=msg):
            Panel(np.random.randn(3, 4, 5),
                  lrange(5), lrange(4), lrange(5))

        msg = (r"Shape of passed values is \(3, 4, 5\), "
               r"indices imply \(5, 5, 4\)")
        with pytest.raises(ValueError, match=msg):
            Panel(np.random.randn(3, 4, 5),
                  lrange(5), lrange(5), lrange(4))

    def test_apply_slabs(self):
        # with multi-indexes
        # GH7469
        index = MultiIndex.from_tuples([('one', 'a'), ('one', 'b'), (
            'two', 'a'), ('two', 'b')])
        dfa = DataFrame(np.array(np.arange(12, dtype='int64')).reshape(
            4, 3), columns=list("ABC"), index=index)
        dfb = DataFrame(np.array(np.arange(10, 22, dtype='int64')).reshape(
            4, 3), columns=list("ABC"), index=index)
        p = Panel({'f': dfa, 'g': dfb})
        result = p.apply(lambda x: x.sum(), axis=0)

        # on windows this will be in32
        result = result.astype('int64')
        expected = p.sum(0)
        assert_frame_equal(result, expected)

    def test_apply_no_or_zero_ndim(self):
        # GH10332
        self.panel = Panel(np.random.rand(5, 5, 5))

        result_int = self.panel.apply(lambda df: 0, axis=[1, 2])
        result_float = self.panel.apply(lambda df: 0.0, axis=[1, 2])
        result_int64 = self.panel.apply(
            lambda df: np.int64(0), axis=[1, 2])
        result_float64 = self.panel.apply(lambda df: np.float64(0.0),
                                          axis=[1, 2])

        expected_int = expected_int64 = Series([0] * 5)
        expected_float = expected_float64 = Series([0.0] * 5)

        assert_series_equal(result_int, expected_int)
        assert_series_equal(result_int64, expected_int64)
        assert_series_equal(result_float, expected_float)
        assert_series_equal(result_float64, expected_float64)

    def test_fillna(self):
        # limit not implemented when only value is specified
        p = Panel(np.random.randn(3, 4, 5))
        p.iloc[0:2, 0:2, 0:2] = np.nan
        pytest.raises(NotImplementedError,
                      lambda: p.fillna(999, limit=1))

    def test_to_frame_multi_major(self):
        idx = MultiIndex.from_tuples(
            [(1, 'one'), (1, 'two'), (2, 'one'), (2, 'two')])
        df = DataFrame([[1, 'a', 1], [2, 'b', 1],
                        [3, 'c', 1], [4, 'd', 1]],
                       columns=['A', 'B', 'C'], index=idx)
        wp = Panel({'i1': df, 'i2': df})
        expected_idx = MultiIndex.from_tuples(
            [
                (1, 'one', 'A'), (1, 'one', 'B'),
                (1, 'one', 'C'), (1, 'two', 'A'),
                (1, 'two', 'B'), (1, 'two', 'C'),
                (2, 'one', 'A'), (2, 'one', 'B'),
                (2, 'one', 'C'), (2, 'two', 'A'),
                (2, 'two', 'B'), (2, 'two', 'C')
            ],
            names=[None, None, 'minor'])
        expected = DataFrame({'i1': [1, 'a', 1, 2, 'b', 1, 3,
                                     'c', 1, 4, 'd', 1],
                              'i2': [1, 'a', 1, 2, 'b',
                                     1, 3, 'c', 1, 4, 'd', 1]},
                             index=expected_idx)
        result = wp.to_frame()
        assert_frame_equal(result, expected)

        wp.iloc[0, 0].iloc[0] = np.nan  # BUG on setting. GH #5773
        result = wp.to_frame()
        assert_frame_equal(result, expected[1:])

        idx = MultiIndex.from_tuples(
            [(1, 'two'), (1, 'one'), (2, 'one'), (np.nan, 'two')])
        df = DataFrame([[1, 'a', 1], [2, 'b', 1],
                        [3, 'c', 1], [4, 'd', 1]],
                       columns=['A', 'B', 'C'], index=idx)
        wp = Panel({'i1': df, 'i2': df})
        ex_idx = MultiIndex.from_tuples([(1, 'two', 'A'), (1, 'two', 'B'),
                                         (1, 'two', 'C'),
                                         (1, 'one', 'A'),
                                         (1, 'one', 'B'),
                                         (1, 'one', 'C'),
                                         (2, 'one', 'A'),
                                         (2, 'one', 'B'),
                                         (2, 'one', 'C'),
                                         (np.nan, 'two', 'A'),
                                         (np.nan, 'two', 'B'),
                                         (np.nan, 'two', 'C')],
                                        names=[None, None, 'minor'])
        expected.index = ex_idx
        result = wp.to_frame()
        assert_frame_equal(result, expected)

    def test_to_frame_multi_major_minor(self):
        cols = MultiIndex(levels=[['C_A', 'C_B'], ['C_1', 'C_2']],
                          codes=[[0, 0, 1, 1], [0, 1, 0, 1]])
        idx = MultiIndex.from_tuples([(1, 'one'), (1, 'two'), (2, 'one'), (
            2, 'two'), (3, 'three'), (4, 'four')])
        df = DataFrame([[1, 2, 11, 12], [3, 4, 13, 14],
                        ['a', 'b', 'w', 'x'],
                        ['c', 'd', 'y', 'z'], [-1, -2, -3, -4],
                        [-5, -6, -7, -8]], columns=cols, index=idx)
        wp = Panel({'i1': df, 'i2': df})

        exp_idx = MultiIndex.from_tuples(
            [(1, 'one', 'C_A', 'C_1'), (1, 'one', 'C_A', 'C_2'),
             (1, 'one', 'C_B', 'C_1'), (1, 'one', 'C_B', 'C_2'),
             (1, 'two', 'C_A', 'C_1'), (1, 'two', 'C_A', 'C_2'),
             (1, 'two', 'C_B', 'C_1'), (1, 'two', 'C_B', 'C_2'),
             (2, 'one', 'C_A', 'C_1'), (2, 'one', 'C_A', 'C_2'),
             (2, 'one', 'C_B', 'C_1'), (2, 'one', 'C_B', 'C_2'),
             (2, 'two', 'C_A', 'C_1'), (2, 'two', 'C_A', 'C_2'),
             (2, 'two', 'C_B', 'C_1'), (2, 'two', 'C_B', 'C_2'),
             (3, 'three', 'C_A', 'C_1'), (3, 'three', 'C_A', 'C_2'),
             (3, 'three', 'C_B', 'C_1'), (3, 'three', 'C_B', 'C_2'),
             (4, 'four', 'C_A', 'C_1'), (4, 'four', 'C_A', 'C_2'),
             (4, 'four', 'C_B', 'C_1'), (4, 'four', 'C_B', 'C_2')],
            names=[None, None, None, None])
        exp_val = [[1, 1], [2, 2], [11, 11], [12, 12],
                   [3, 3], [4, 4],
                   [13, 13], [14, 14], ['a', 'a'],
                   ['b', 'b'], ['w', 'w'],
                   ['x', 'x'], ['c', 'c'], ['d', 'd'], [
                       'y', 'y'], ['z', 'z'],
                   [-1, -1], [-2, -2], [-3, -3], [-4, -4],
                   [-5, -5], [-6, -6],
                   [-7, -7], [-8, -8]]
        result = wp.to_frame()
        expected = DataFrame(exp_val, columns=['i1', 'i2'], index=exp_idx)
        assert_frame_equal(result, expected)

    def test_to_frame_multi_drop_level(self):
        idx = MultiIndex.from_tuples([(1, 'one'), (2, 'one'), (2, 'two')])
        df = DataFrame({'A': [np.nan, 1, 2]}, index=idx)
        wp = Panel({'i1': df, 'i2': df})
        result = wp.to_frame()
        exp_idx = MultiIndex.from_tuples(
            [(2, 'one', 'A'), (2, 'two', 'A')],
            names=[None, None, 'minor'])
        expected = DataFrame({'i1': [1., 2], 'i2': [1., 2]}, index=exp_idx)
        assert_frame_equal(result, expected)

    def test_panel_dups(self):

        # GH 4960
        # duplicates in an index

        # items
        data = np.random.randn(5, 100, 5)
        no_dup_panel = Panel(data, items=list("ABCDE"))
        panel = Panel(data, items=list("AACDE"))

        expected = no_dup_panel['A']
        result = panel.iloc[0]
        assert_frame_equal(result, expected)

        expected = no_dup_panel['E']
        result = panel.loc['E']
        assert_frame_equal(result, expected)

        # major
        data = np.random.randn(5, 5, 5)
        no_dup_panel = Panel(data, major_axis=list("ABCDE"))
        panel = Panel(data, major_axis=list("AACDE"))

        expected = no_dup_panel.loc[:, 'A']
        result = panel.iloc[:, 0]
        assert_frame_equal(result, expected)

        expected = no_dup_panel.loc[:, 'E']
        result = panel.loc[:, 'E']
        assert_frame_equal(result, expected)

        # minor
        data = np.random.randn(5, 100, 5)
        no_dup_panel = Panel(data, minor_axis=list("ABCDE"))
        panel = Panel(data, minor_axis=list("AACDE"))

        expected = no_dup_panel.loc[:, :, 'A']
        result = panel.iloc[:, :, 0]
        assert_frame_equal(result, expected)

        expected = no_dup_panel.loc[:, :, 'E']
        result = panel.loc[:, :, 'E']
        assert_frame_equal(result, expected)

    def test_filter(self):
        pass

    def test_shift(self):
        # mixed dtypes #6959
        data = [('item ' + ch, makeMixedDataFrame())
                for ch in list('abcde')]
        data = dict(data)
        mixed_panel = Panel.from_dict(data, orient='minor')
        shifted = mixed_panel.shift(1)
        assert_series_equal(mixed_panel.dtypes, shifted.dtypes)

    def test_numpy_round(self):
        values = [[[-3.2, 2.2], [0, -4.8213], [3.123, 123.12],
                   [-1566.213, 88.88], [-12, 94.5]],
                  [[-5.82, 3.5], [6.21, -73.272], [-9.087, 23.12],
                   [272.212, -99.99], [23, -76.5]]]
        p = Panel(values, items=['Item1', 'Item2'],
                  major_axis=date_range('1/1/2000', periods=5),
                  minor_axis=['A', 'B'])

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.round(p, out=p)

    # removing Panel before NumPy enforces, so just ignore
    @pytest.mark.filterwarnings("ignore:Using a non-tuple:FutureWarning")
    def test_multiindex_get(self):
        ind = MultiIndex.from_tuples(
            [('a', 1), ('a', 2), ('b', 1), ('b', 2)],
            names=['first', 'second'])
        wp = Panel(np.random.random((4, 5, 5)),
                   items=ind,
                   major_axis=np.arange(5),
                   minor_axis=np.arange(5))
        f1 = wp['a']
        f2 = wp.loc['a']

        assert (f1.items == [1, 2]).all()
        assert (f2.items == [1, 2]).all()

        MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)],
                               names=['first', 'second'])

    def test_repr_empty(self):
        empty = Panel()
        repr(empty)

    @pytest.mark.parametrize('bad_kwarg, exception, msg', [
        # errors must be 'ignore' or 'raise'
        ({'errors': 'something'}, ValueError, 'The parameter errors must.*'),
        ({'join': 'inner'}, NotImplementedError, 'Only left join is supported')
    ])
    def test_update_raise_bad_parameter(self, bad_kwarg, exception, msg):
        pan = Panel([[[1.5, np.nan, 3.]]])
        with pytest.raises(exception, match=msg):
            pan.update(pan, **bad_kwarg)

    def test_update_raise_on_overlap(self):
        pan = Panel([[[1.5, np.nan, 3.], [1.5, np.nan, 3.],
                      [1.5, np.nan, 3.],
                      [1.5, np.nan, 3.]],
                     [[1.5, np.nan, 3.], [1.5, np.nan, 3.],
                      [1.5, np.nan, 3.],
                      [1.5, np.nan, 3.]]])

        with pytest.raises(ValueError, match='Data overlaps'):
            pan.update(pan, errors='raise')

    @pytest.mark.parametrize('raise_conflict', [True, False])
    def test_update_deprecation(self, raise_conflict):
        pan = Panel([[[1.5, np.nan, 3.]]])
        other = Panel([[[]]])
        with tm.assert_produces_warning(FutureWarning):
            pan.update(other, raise_conflict=raise_conflict)


def test_panel_index():
    index = panelm.panel_index([1, 2, 3, 4], [1, 2, 3])
    expected = MultiIndex.from_arrays([np.tile([1, 2, 3, 4], 3),
                                       np.repeat([1, 2, 3], 4)],
                                      names=['time', 'panel'])
    tm.assert_index_equal(index, expected)


@pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
def test_panel_np_all():
    wp = Panel({"A": DataFrame({'b': [1, 2]})})
    result = np.all(wp)
    assert result == np.bool_(True)
