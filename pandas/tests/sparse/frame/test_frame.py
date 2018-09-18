# pylint: disable-msg=E1101,W0612

import operator

import pytest
from numpy import nan
import numpy as np
import pandas as pd

from pandas import Series, DataFrame, bdate_range, Panel
from pandas.core.indexes.datetimes import DatetimeIndex
from pandas.tseries.offsets import BDay
from pandas.util import testing as tm
from pandas.compat import lrange
from pandas import compat
from pandas.core.sparse import frame as spf

from pandas._libs.sparse import BlockIndex, IntIndex
from pandas.core.sparse.api import SparseSeries, SparseDataFrame, SparseArray
from pandas.tests.frame.test_api import SharedWithSparse


class TestSparseDataFrame(SharedWithSparse):
    klass = SparseDataFrame

    # SharedWithSparse tests use generic, klass-agnostic assertion
    _assert_frame_equal = staticmethod(tm.assert_sp_frame_equal)
    _assert_series_equal = staticmethod(tm.assert_sp_series_equal)

    def test_fill_value_when_combine_const(self):
        # GH12723
        dat = np.array([0, 1, np.nan, 3, 4, 5], dtype='float')
        df = SparseDataFrame({'foo': dat}, index=range(6))

        exp = df.fillna(0).add(2)
        res = df.add(2, fill_value=0)
        tm.assert_sp_frame_equal(res, exp)

    def test_values(self, empty_frame, float_frame):
        empty = empty_frame.values
        assert empty.shape == (0, 0)

        no_cols = SparseDataFrame(index=np.arange(10))
        mat = no_cols.values
        assert mat.shape == (10, 0)

        no_index = SparseDataFrame(columns=np.arange(10))
        mat = no_index.values
        assert mat.shape == (0, 10)

    def test_copy(self, float_frame):
        cp = float_frame.copy()
        assert isinstance(cp, SparseDataFrame)
        tm.assert_sp_frame_equal(cp, float_frame)

        # as of v0.15.0
        # this is now identical (but not is_a )
        assert cp.index.identical(float_frame.index)

    def test_constructor(self, float_frame, float_frame_int_kind,
                         float_frame_fill0):
        for col, series in compat.iteritems(float_frame):
            assert isinstance(series, SparseSeries)

        assert isinstance(float_frame_int_kind['A'].sp_index, IntIndex)

        # constructed zframe from matrix above
        assert float_frame_fill0['A'].fill_value == 0
        tm.assert_numpy_array_equal(pd.SparseArray([1., 2., 3., 4., 5., 6.]),
                                    float_frame_fill0['A'].values)
        tm.assert_numpy_array_equal(np.array([0., 0., 0., 0., 1., 2.,
                                              3., 4., 5., 6.]),
                                    float_frame_fill0['A'].to_dense().values)

        # construct no data
        sdf = SparseDataFrame(columns=np.arange(10), index=np.arange(10))
        for col, series in compat.iteritems(sdf):
            assert isinstance(series, SparseSeries)

        # construct from nested dict
        data = {}
        for c, s in compat.iteritems(float_frame):
            data[c] = s.to_dict()

        sdf = SparseDataFrame(data)
        tm.assert_sp_frame_equal(sdf, float_frame)

        # TODO: test data is copied from inputs

        # init dict with different index
        idx = float_frame.index[:5]
        cons = SparseDataFrame(
            float_frame, index=idx, columns=float_frame.columns,
            default_fill_value=float_frame.default_fill_value,
            default_kind=float_frame.default_kind, copy=True)
        reindexed = float_frame.reindex(idx)

        tm.assert_sp_frame_equal(cons, reindexed, exact_indices=False)

        # assert level parameter breaks reindex
        with pytest.raises(TypeError):
            float_frame.reindex(idx, level=0)

        repr(float_frame)

    def test_constructor_dict_order(self):
        # GH19018
        # initialization ordering: by insertion order if python>= 3.6, else
        # order by value
        d = {'b': [2, 3], 'a': [0, 1]}
        frame = SparseDataFrame(data=d)
        if compat.PY36:
            expected = SparseDataFrame(data=d, columns=list('ba'))
        else:
            expected = SparseDataFrame(data=d, columns=list('ab'))
        tm.assert_sp_frame_equal(frame, expected)

    def test_constructor_ndarray(self, float_frame):
        # no index or columns
        sp = SparseDataFrame(float_frame.values)

        # 1d
        sp = SparseDataFrame(float_frame['A'].values, index=float_frame.index,
                             columns=['A'])
        tm.assert_sp_frame_equal(sp, float_frame.reindex(columns=['A']))

        # raise on level argument
        pytest.raises(TypeError, float_frame.reindex, columns=['A'],
                      level=1)

        # wrong length index / columns
        with tm.assert_raises_regex(ValueError, "^Index length"):
            SparseDataFrame(float_frame.values, index=float_frame.index[:-1])

        with tm.assert_raises_regex(ValueError, "^Column length"):
            SparseDataFrame(float_frame.values,
                            columns=float_frame.columns[:-1])

    # GH 9272
    def test_constructor_empty(self):
        sp = SparseDataFrame()
        assert len(sp.index) == 0
        assert len(sp.columns) == 0

    def test_constructor_dataframe(self, float_frame):
        dense = float_frame.to_dense()
        sp = SparseDataFrame(dense)
        tm.assert_sp_frame_equal(sp, float_frame)

    def test_constructor_convert_index_once(self):
        arr = np.array([1.5, 2.5, 3.5])
        sdf = SparseDataFrame(columns=lrange(4), index=arr)
        assert sdf[0].index is sdf[1].index

    def test_constructor_from_series(self):

        # GH 2873
        x = Series(np.random.randn(10000), name='a')
        x = x.to_sparse(fill_value=0)
        assert isinstance(x, SparseSeries)
        df = SparseDataFrame(x)
        assert isinstance(df, SparseDataFrame)

        x = Series(np.random.randn(10000), name='a')
        y = Series(np.random.randn(10000), name='b')
        x2 = x.astype(float)
        x2.loc[:9998] = np.NaN
        # TODO: x_sparse is unused...fix
        x_sparse = x2.to_sparse(fill_value=np.NaN)  # noqa

        # Currently fails too with weird ufunc error
        # df1 = SparseDataFrame([x_sparse, y])

        y.loc[:9998] = 0
        # TODO: y_sparse is unsused...fix
        y_sparse = y.to_sparse(fill_value=0)  # noqa
        # without sparse value raises error
        # df2 = SparseDataFrame([x2_sparse, y])

    def test_constructor_from_dense_series(self):
        # GH 19393
        # series with name
        x = Series(np.random.randn(10000), name='a')
        result = SparseDataFrame(x)
        expected = x.to_frame().to_sparse()
        tm.assert_sp_frame_equal(result, expected)

        # series with no name
        x = Series(np.random.randn(10000))
        result = SparseDataFrame(x)
        expected = x.to_frame().to_sparse()
        tm.assert_sp_frame_equal(result, expected)

    def test_constructor_from_unknown_type(self):
        # GH 19393
        class Unknown(object):
            pass
        with pytest.raises(TypeError,
                           message='SparseDataFrame called with unknown type '
                                   '"Unknown" for data argument'):
            SparseDataFrame(Unknown())

    def test_constructor_preserve_attr(self):
        # GH 13866
        arr = pd.SparseArray([1, 0, 3, 0], dtype=np.int64, fill_value=0)
        assert arr.dtype == np.int64
        assert arr.fill_value == 0

        df = pd.SparseDataFrame({'x': arr})
        assert df['x'].dtype == np.int64
        assert df['x'].fill_value == 0

        s = pd.SparseSeries(arr, name='x')
        assert s.dtype == np.int64
        assert s.fill_value == 0

        df = pd.SparseDataFrame(s)
        assert df['x'].dtype == np.int64
        assert df['x'].fill_value == 0

        df = pd.SparseDataFrame({'x': s})
        assert df['x'].dtype == np.int64
        assert df['x'].fill_value == 0

    def test_constructor_nan_dataframe(self):
        # GH 10079
        trains = np.arange(100)
        thresholds = [10, 20, 30, 40, 50, 60]
        tuples = [(i, j) for i in trains for j in thresholds]
        index = pd.MultiIndex.from_tuples(tuples,
                                          names=['trains', 'thresholds'])
        matrix = np.empty((len(index), len(trains)))
        matrix.fill(np.nan)
        df = pd.DataFrame(matrix, index=index, columns=trains, dtype=float)
        result = df.to_sparse()
        expected = pd.SparseDataFrame(matrix, index=index, columns=trains,
                                      dtype=float)
        tm.assert_sp_frame_equal(result, expected)

    def test_type_coercion_at_construction(self):
        # GH 15682
        result = pd.SparseDataFrame(
            {'a': [1, 0, 0], 'b': [0, 1, 0], 'c': [0, 0, 1]}, dtype='uint8',
            default_fill_value=0)
        expected = pd.SparseDataFrame(
            {'a': pd.SparseSeries([1, 0, 0], dtype='uint8'),
             'b': pd.SparseSeries([0, 1, 0], dtype='uint8'),
             'c': pd.SparseSeries([0, 0, 1], dtype='uint8')},
            default_fill_value=0)
        tm.assert_sp_frame_equal(result, expected)

    def test_dtypes(self):
        df = DataFrame(np.random.randn(10000, 4))
        df.loc[:9998] = np.nan
        sdf = df.to_sparse()

        result = sdf.get_dtype_counts()
        expected = Series({'float64': 4})
        tm.assert_series_equal(result, expected)

    def test_shape(self, float_frame, float_frame_int_kind,
                   float_frame_fill0, float_frame_fill2):
        # see gh-10452
        assert float_frame.shape == (10, 4)
        assert float_frame_int_kind.shape == (10, 4)
        assert float_frame_fill0.shape == (10, 4)
        assert float_frame_fill2.shape == (10, 4)

    def test_str(self):
        df = DataFrame(np.random.randn(10000, 4))
        df.loc[:9998] = np.nan

        sdf = df.to_sparse()
        str(sdf)

    def test_array_interface(self, float_frame):
        res = np.sqrt(float_frame)
        dres = np.sqrt(float_frame.to_dense())
        tm.assert_frame_equal(res.to_dense(), dres)

    def test_pickle(self, float_frame, float_frame_int_kind, float_frame_dense,
                    float_frame_fill0, float_frame_fill0_dense,
                    float_frame_fill2, float_frame_fill2_dense):

        def _test_roundtrip(frame, orig):
            result = tm.round_trip_pickle(frame)
            tm.assert_sp_frame_equal(frame, result)
            tm.assert_frame_equal(result.to_dense(), orig, check_dtype=False)

        _test_roundtrip(SparseDataFrame(), DataFrame())
        _test_roundtrip(float_frame, float_frame_dense)
        _test_roundtrip(float_frame_int_kind, float_frame_dense)
        _test_roundtrip(float_frame_fill0, float_frame_fill0_dense)
        _test_roundtrip(float_frame_fill2, float_frame_fill2_dense)

    def test_dense_to_sparse(self):
        df = DataFrame({'A': [nan, nan, nan, 1, 2],
                        'B': [1, 2, nan, nan, nan]})
        sdf = df.to_sparse()
        assert isinstance(sdf, SparseDataFrame)
        assert np.isnan(sdf.default_fill_value)
        assert isinstance(sdf['A'].sp_index, BlockIndex)
        tm.assert_frame_equal(sdf.to_dense(), df)

        sdf = df.to_sparse(kind='integer')
        assert isinstance(sdf['A'].sp_index, IntIndex)

        df = DataFrame({'A': [0, 0, 0, 1, 2],
                        'B': [1, 2, 0, 0, 0]}, dtype=float)
        sdf = df.to_sparse(fill_value=0)
        assert sdf.default_fill_value == 0
        tm.assert_frame_equal(sdf.to_dense(), df)

    def test_density(self):
        df = SparseSeries([nan, nan, nan, 0, 1, 2, 3, 4, 5, 6])
        assert df.density == 0.7

        df = SparseDataFrame({'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
                              'B': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
                              'C': np.arange(10),
                              'D': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]})

        assert df.density == 0.75

    def test_sparse_to_dense(self):
        pass

    def test_sparse_series_ops(self, float_frame):
        self._check_frame_ops(float_frame)

    def test_sparse_series_ops_i(self, float_frame_int_kind):
        self._check_frame_ops(float_frame_int_kind)

    def test_sparse_series_ops_z(self, float_frame_fill0):
        self._check_frame_ops(float_frame_fill0)

    def test_sparse_series_ops_fill(self, float_frame_fill2):
        self._check_frame_ops(float_frame_fill2)

    def _check_frame_ops(self, frame):

        def _compare_to_dense(a, b, da, db, op):
            sparse_result = op(a, b)
            dense_result = op(da, db)

            fill = sparse_result.default_fill_value
            dense_result = dense_result.to_sparse(fill_value=fill)
            tm.assert_sp_frame_equal(sparse_result, dense_result,
                                     exact_indices=False)

            if isinstance(a, DataFrame) and isinstance(db, DataFrame):
                mixed_result = op(a, db)
                assert isinstance(mixed_result, SparseDataFrame)
                tm.assert_sp_frame_equal(mixed_result, sparse_result,
                                         exact_indices=False)

        opnames = ['add', 'sub', 'mul', 'truediv', 'floordiv']
        ops = [getattr(operator, name) for name in opnames]

        fidx = frame.index

        # time series operations

        series = [frame['A'], frame['B'], frame['C'], frame['D'],
                  frame['A'].reindex(fidx[:7]), frame['A'].reindex(fidx[::2]),
                  SparseSeries(
                      [], index=[])]

        for op in opnames:
            _compare_to_dense(frame, frame[::2], frame.to_dense(),
                              frame[::2].to_dense(), getattr(operator, op))

            # 2304, no auto-broadcasting
            for i, s in enumerate(series):
                f = lambda a, b: getattr(a, op)(b, axis='index')
                _compare_to_dense(frame, s, frame.to_dense(), s.to_dense(), f)

                # rops are not implemented
                # _compare_to_dense(s, frame, s.to_dense(),
                #                   frame.to_dense(), f)

                # cross-sectional operations
        series = [frame.xs(fidx[0]), frame.xs(fidx[3]), frame.xs(fidx[5]),
                  frame.xs(fidx[7]), frame.xs(fidx[5])[:2]]

        for op in ops:
            for s in series:
                _compare_to_dense(frame, s, frame.to_dense(), s, op)
                _compare_to_dense(s, frame, s, frame.to_dense(), op)

        # it works!
        result = frame + frame.loc[:, ['A', 'B']]  # noqa

    def test_op_corners(self, float_frame, empty_frame):
        empty = empty_frame + empty_frame
        assert empty.empty

        foo = float_frame + empty_frame
        assert isinstance(foo.index, DatetimeIndex)
        tm.assert_frame_equal(foo, float_frame * np.nan)

        foo = empty_frame + float_frame
        tm.assert_frame_equal(foo, float_frame * np.nan)

    def test_scalar_ops(self):
        pass

    def test_getitem(self):
        # 1585 select multiple columns
        sdf = SparseDataFrame(index=[0, 1, 2], columns=['a', 'b', 'c'])

        result = sdf[['a', 'b']]
        exp = sdf.reindex(columns=['a', 'b'])
        tm.assert_sp_frame_equal(result, exp)

        pytest.raises(Exception, sdf.__getitem__, ['a', 'd'])

    def test_iloc(self, float_frame):

        # GH 2227
        result = float_frame.iloc[:, 0]
        assert isinstance(result, SparseSeries)
        tm.assert_sp_series_equal(result, float_frame['A'])

        # preserve sparse index type. #2251
        data = {'A': [0, 1]}
        iframe = SparseDataFrame(data, default_kind='integer')
        tm.assert_class_equal(iframe['A'].sp_index,
                              iframe.iloc[:, 0].sp_index)

    def test_set_value(self, float_frame):

        # ok, as the index gets converted to object
        frame = float_frame.copy()
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            res = frame.set_value('foobar', 'B', 1.5)
        assert res.index.dtype == 'object'

        res = float_frame
        res.index = res.index.astype(object)

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            res = float_frame.set_value('foobar', 'B', 1.5)
        assert res is not float_frame
        assert res.index[-1] == 'foobar'
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            assert res.get_value('foobar', 'B') == 1.5

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            res2 = res.set_value('foobar', 'qux', 1.5)
        assert res2 is not res
        tm.assert_index_equal(res2.columns,
                              pd.Index(list(float_frame.columns) + ['qux']))
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            assert res2.get_value('foobar', 'qux') == 1.5

    def test_fancy_index_misc(self, float_frame):
        # axis = 0
        sliced = float_frame.iloc[-2:, :]
        expected = float_frame.reindex(index=float_frame.index[-2:])
        tm.assert_sp_frame_equal(sliced, expected)

        # axis = 1
        sliced = float_frame.iloc[:, -2:]
        expected = float_frame.reindex(columns=float_frame.columns[-2:])
        tm.assert_sp_frame_equal(sliced, expected)

    def test_getitem_overload(self, float_frame):
        # slicing
        sl = float_frame[:20]
        tm.assert_sp_frame_equal(sl,
                                 float_frame.reindex(float_frame.index[:20]))

        # boolean indexing
        d = float_frame.index[5]
        indexer = float_frame.index > d

        subindex = float_frame.index[indexer]
        subframe = float_frame[indexer]

        tm.assert_index_equal(subindex, subframe.index)
        pytest.raises(Exception, float_frame.__getitem__, indexer[:-1])

    def test_setitem(self, float_frame, float_frame_int_kind,
                     float_frame_dense,
                     float_frame_fill0, float_frame_fill0_dense,
                     float_frame_fill2, float_frame_fill2_dense):

        def _check_frame(frame, orig):
            N = len(frame)

            # insert SparseSeries
            frame['E'] = frame['A']
            assert isinstance(frame['E'], SparseSeries)
            tm.assert_sp_series_equal(frame['E'], frame['A'],
                                      check_names=False)

            # insert SparseSeries differently-indexed
            to_insert = frame['A'][::2]
            frame['E'] = to_insert
            expected = to_insert.to_dense().reindex(frame.index)
            result = frame['E'].to_dense()
            tm.assert_series_equal(result, expected, check_names=False)
            assert result.name == 'E'

            # insert Series
            frame['F'] = frame['A'].to_dense()
            assert isinstance(frame['F'], SparseSeries)
            tm.assert_sp_series_equal(frame['F'], frame['A'],
                                      check_names=False)

            # insert Series differently-indexed
            to_insert = frame['A'].to_dense()[::2]
            frame['G'] = to_insert
            expected = to_insert.reindex(frame.index)
            expected.name = 'G'
            tm.assert_series_equal(frame['G'].to_dense(), expected)

            # insert ndarray
            frame['H'] = np.random.randn(N)
            assert isinstance(frame['H'], SparseSeries)

            to_sparsify = np.random.randn(N)
            to_sparsify[N // 2:] = frame.default_fill_value
            frame['I'] = to_sparsify
            assert len(frame['I'].sp_values) == N // 2

            # insert ndarray wrong size
            pytest.raises(Exception, frame.__setitem__, 'foo',
                          np.random.randn(N - 1))

            # scalar value
            frame['J'] = 5
            assert len(frame['J'].sp_values) == N
            assert (frame['J'].sp_values == 5).all()

            frame['K'] = frame.default_fill_value
            assert len(frame['K'].sp_values) == 0

        _check_frame(float_frame, float_frame_dense)
        _check_frame(float_frame_int_kind, float_frame_dense)
        _check_frame(float_frame_fill0, float_frame_fill0_dense)
        _check_frame(float_frame_fill2, float_frame_fill2_dense)

    def test_setitem_corner(self, float_frame):
        float_frame['a'] = float_frame['B']
        tm.assert_sp_series_equal(float_frame['a'], float_frame['B'],
                                  check_names=False)

    def test_setitem_array(self, float_frame):
        arr = float_frame['B']

        float_frame['E'] = arr
        tm.assert_sp_series_equal(float_frame['E'], float_frame['B'],
                                  check_names=False)

        float_frame['F'] = arr[:-1]
        index = float_frame.index[:-1]
        tm.assert_sp_series_equal(float_frame['E'].reindex(index),
                                  float_frame['F'].reindex(index),
                                  check_names=False)

    def test_setitem_chained_no_consolidate(self):
        # https://github.com/pandas-dev/pandas/pull/19268
        # issuecomment-361696418
        # chained setitem used to cause consolidation
        sdf = pd.SparseDataFrame([[np.nan, 1], [2, np.nan]])
        with pd.option_context('mode.chained_assignment', None):
            sdf[0][1] = 2
        assert len(sdf._data.blocks) == 2

    def test_delitem(self, float_frame):
        A = float_frame['A']
        C = float_frame['C']

        del float_frame['B']
        assert 'B' not in float_frame
        tm.assert_sp_series_equal(float_frame['A'], A)
        tm.assert_sp_series_equal(float_frame['C'], C)

        del float_frame['D']
        assert 'D' not in float_frame

        del float_frame['A']
        assert 'A' not in float_frame

    def test_set_columns(self, float_frame):
        float_frame.columns = float_frame.columns
        pytest.raises(Exception, setattr, float_frame, 'columns',
                      float_frame.columns[:-1])

    def test_set_index(self, float_frame):
        float_frame.index = float_frame.index
        pytest.raises(Exception, setattr, float_frame, 'index',
                      float_frame.index[:-1])

    def test_append(self, float_frame):
        a = float_frame[:5]
        b = float_frame[5:]

        appended = a.append(b)
        tm.assert_sp_frame_equal(appended, float_frame, exact_indices=False)

        a = float_frame.iloc[:5, :3]
        b = float_frame.iloc[5:]
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            # Stacklevel is set for pd.concat, not append
            appended = a.append(b)
        tm.assert_sp_frame_equal(appended.iloc[:, :3], float_frame.iloc[:, :3],
                                 exact_indices=False)

        a = a[['B', 'C', 'A']].head(2)
        b = b.head(2)

        expected = pd.SparseDataFrame({
            "B": [0., 1, None, 3],
            "C": [0., 1, 5, 6],
            "A": [None, None, 2, 3],
            "D": [None, None, 5, None],
        }, index=a.index | b.index, columns=['B', 'C', 'A', 'D'])
        with tm.assert_produces_warning(None):
            appended = a.append(b, sort=False)

        tm.assert_frame_equal(appended, expected)

        with tm.assert_produces_warning(None):
            appended = a.append(b, sort=True)

        tm.assert_sp_frame_equal(appended, expected[['A', 'B', 'C', 'D']])

    def test_astype(self):
        sparse = pd.SparseDataFrame({'A': SparseArray([1, 2, 3, 4],
                                                      dtype=np.int64),
                                     'B': SparseArray([4, 5, 6, 7],
                                                      dtype=np.int64)})
        assert sparse['A'].dtype == np.int64
        assert sparse['B'].dtype == np.int64

        res = sparse.astype(np.float64)
        exp = pd.SparseDataFrame({'A': SparseArray([1., 2., 3., 4.],
                                                   fill_value=0.),
                                  'B': SparseArray([4., 5., 6., 7.],
                                                   fill_value=0.)},
                                 default_fill_value=np.nan)
        tm.assert_sp_frame_equal(res, exp)
        assert res['A'].dtype == np.float64
        assert res['B'].dtype == np.float64

        sparse = pd.SparseDataFrame({'A': SparseArray([0, 2, 0, 4],
                                                      dtype=np.int64),
                                     'B': SparseArray([0, 5, 0, 7],
                                                      dtype=np.int64)},
                                    default_fill_value=0)
        assert sparse['A'].dtype == np.int64
        assert sparse['B'].dtype == np.int64

        res = sparse.astype(np.float64)
        exp = pd.SparseDataFrame({'A': SparseArray([0., 2., 0., 4.],
                                                   fill_value=0.),
                                  'B': SparseArray([0., 5., 0., 7.],
                                                   fill_value=0.)},
                                 default_fill_value=0.)
        tm.assert_sp_frame_equal(res, exp)
        assert res['A'].dtype == np.float64
        assert res['B'].dtype == np.float64

    def test_astype_bool(self):
        sparse = pd.SparseDataFrame({'A': SparseArray([0, 2, 0, 4],
                                                      fill_value=0,
                                                      dtype=np.int64),
                                     'B': SparseArray([0, 5, 0, 7],
                                                      fill_value=0,
                                                      dtype=np.int64)},
                                    default_fill_value=0)
        assert sparse['A'].dtype == np.int64
        assert sparse['B'].dtype == np.int64

        res = sparse.astype(bool)
        exp = pd.SparseDataFrame({'A': SparseArray([False, True, False, True],
                                                   dtype=np.bool,
                                                   fill_value=False),
                                  'B': SparseArray([False, True, False, True],
                                                   dtype=np.bool,
                                                   fill_value=False)},
                                 default_fill_value=False)
        tm.assert_sp_frame_equal(res, exp)
        assert res['A'].dtype == np.bool
        assert res['B'].dtype == np.bool

    def test_fillna(self, float_frame_fill0, float_frame_fill0_dense):
        df = float_frame_fill0.reindex(lrange(5))
        dense = float_frame_fill0_dense.reindex(lrange(5))

        result = df.fillna(0)
        expected = dense.fillna(0)
        tm.assert_sp_frame_equal(result, expected.to_sparse(fill_value=0),
                                 exact_indices=False)
        tm.assert_frame_equal(result.to_dense(), expected)

        result = df.copy()
        result.fillna(0, inplace=True)
        expected = dense.fillna(0)

        tm.assert_sp_frame_equal(result, expected.to_sparse(fill_value=0),
                                 exact_indices=False)
        tm.assert_frame_equal(result.to_dense(), expected)

        result = df.copy()
        result = df['A']
        result.fillna(0, inplace=True)

        expected = dense['A'].fillna(0)
        # this changes internal SparseArray repr
        # tm.assert_sp_series_equal(result, expected.to_sparse(fill_value=0))
        tm.assert_series_equal(result.to_dense(), expected)

    def test_fillna_fill_value(self):
        df = pd.DataFrame({'A': [1, 0, 0], 'B': [np.nan, np.nan, 4]})

        sparse = pd.SparseDataFrame(df)
        tm.assert_frame_equal(sparse.fillna(-1).to_dense(),
                              df.fillna(-1), check_dtype=False)

        sparse = pd.SparseDataFrame(df, default_fill_value=0)
        tm.assert_frame_equal(sparse.fillna(-1).to_dense(),
                              df.fillna(-1), check_dtype=False)

    def test_sparse_frame_pad_backfill_limit(self):
        index = np.arange(10)
        df = DataFrame(np.random.randn(10, 4), index=index)
        sdf = df.to_sparse()

        result = sdf[:2].reindex(index, method='pad', limit=5)

        expected = sdf[:2].reindex(index).fillna(method='pad')
        expected = expected.to_dense()
        expected.values[-3:] = np.nan
        expected = expected.to_sparse()
        tm.assert_frame_equal(result, expected)

        result = sdf[-2:].reindex(index, method='backfill', limit=5)

        expected = sdf[-2:].reindex(index).fillna(method='backfill')
        expected = expected.to_dense()
        expected.values[:3] = np.nan
        expected = expected.to_sparse()
        tm.assert_frame_equal(result, expected)

    def test_sparse_frame_fillna_limit(self):
        index = np.arange(10)
        df = DataFrame(np.random.randn(10, 4), index=index)
        sdf = df.to_sparse()

        result = sdf[:2].reindex(index)
        result = result.fillna(method='pad', limit=5)

        expected = sdf[:2].reindex(index).fillna(method='pad')
        expected = expected.to_dense()
        expected.values[-3:] = np.nan
        expected = expected.to_sparse()
        tm.assert_frame_equal(result, expected)

        result = sdf[-2:].reindex(index)
        result = result.fillna(method='backfill', limit=5)

        expected = sdf[-2:].reindex(index).fillna(method='backfill')
        expected = expected.to_dense()
        expected.values[:3] = np.nan
        expected = expected.to_sparse()
        tm.assert_frame_equal(result, expected)

    def test_rename(self, float_frame):
        result = float_frame.rename(index=str)
        expected = SparseDataFrame(float_frame.values,
                                   index=float_frame.index.strftime(
                                       "%Y-%m-%d %H:%M:%S"),
                                   columns=list('ABCD'))
        tm.assert_sp_frame_equal(result, expected)

        result = float_frame.rename(columns=lambda x: '%s%d' % (x, 1))
        data = {'A1': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
                'B1': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
                'C1': np.arange(10, dtype=np.float64),
                'D1': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]}
        expected = SparseDataFrame(data, index=float_frame.index)
        tm.assert_sp_frame_equal(result, expected)

    def test_corr(self, float_frame):
        res = float_frame.corr()
        tm.assert_frame_equal(res, float_frame.to_dense().corr())

    def test_describe(self, float_frame):
        float_frame['foo'] = np.nan
        float_frame.get_dtype_counts()
        str(float_frame)
        desc = float_frame.describe()  # noqa

    def test_join(self, float_frame):
        left = float_frame.loc[:, ['A', 'B']]
        right = float_frame.loc[:, ['C', 'D']]
        joined = left.join(right)
        tm.assert_sp_frame_equal(joined, float_frame, exact_indices=False)

        right = float_frame.loc[:, ['B', 'D']]
        pytest.raises(Exception, left.join, right)

        with tm.assert_raises_regex(ValueError,
                                    'Other Series must have a name'):
            float_frame.join(Series(
                np.random.randn(len(float_frame)), index=float_frame.index))

    def test_reindex(self, float_frame, float_frame_int_kind,
                     float_frame_fill0, float_frame_fill2):

        def _check_frame(frame):
            index = frame.index
            sidx = index[::2]
            sidx2 = index[:5]  # noqa

            sparse_result = frame.reindex(sidx)
            dense_result = frame.to_dense().reindex(sidx)
            tm.assert_frame_equal(sparse_result.to_dense(), dense_result)

            tm.assert_frame_equal(frame.reindex(list(sidx)).to_dense(),
                                  dense_result)

            sparse_result2 = sparse_result.reindex(index)
            dense_result2 = dense_result.reindex(index)
            tm.assert_frame_equal(sparse_result2.to_dense(), dense_result2)

            # propagate CORRECT fill value
            tm.assert_almost_equal(sparse_result.default_fill_value,
                                   frame.default_fill_value)
            tm.assert_almost_equal(sparse_result['A'].fill_value,
                                   frame['A'].fill_value)

            # length zero
            length_zero = frame.reindex([])
            assert len(length_zero) == 0
            assert len(length_zero.columns) == len(frame.columns)
            assert len(length_zero['A']) == 0

            # frame being reindexed has length zero
            length_n = length_zero.reindex(index)
            assert len(length_n) == len(frame)
            assert len(length_n.columns) == len(frame.columns)
            assert len(length_n['A']) == len(frame)

            # reindex columns
            reindexed = frame.reindex(columns=['A', 'B', 'Z'])
            assert len(reindexed.columns) == 3
            tm.assert_almost_equal(reindexed['Z'].fill_value,
                                   frame.default_fill_value)
            assert np.isnan(reindexed['Z'].sp_values).all()

        _check_frame(float_frame)
        _check_frame(float_frame_int_kind)
        _check_frame(float_frame_fill0)
        _check_frame(float_frame_fill2)

        # with copy=False
        reindexed = float_frame.reindex(float_frame.index, copy=False)
        reindexed['F'] = reindexed['A']
        assert 'F' in float_frame

        reindexed = float_frame.reindex(float_frame.index)
        reindexed['G'] = reindexed['A']
        assert 'G' not in float_frame

    def test_reindex_fill_value(self, float_frame_fill0,
                                float_frame_fill0_dense):
        rng = bdate_range('20110110', periods=20)

        result = float_frame_fill0.reindex(rng, fill_value=0)
        exp = float_frame_fill0_dense.reindex(rng, fill_value=0)
        exp = exp.to_sparse(float_frame_fill0.default_fill_value)
        tm.assert_sp_frame_equal(result, exp)

    def test_reindex_method(self):

        sparse = SparseDataFrame(data=[[11., 12., 14.],
                                       [21., 22., 24.],
                                       [41., 42., 44.]],
                                 index=[1, 2, 4],
                                 columns=[1, 2, 4],
                                 dtype=float)

        # Over indices

        # default method
        result = sparse.reindex(index=range(6))
        expected = SparseDataFrame(data=[[nan, nan, nan],
                                         [11., 12., 14.],
                                         [21., 22., 24.],
                                         [nan, nan, nan],
                                         [41., 42., 44.],
                                         [nan, nan, nan]],
                                   index=range(6),
                                   columns=[1, 2, 4],
                                   dtype=float)
        tm.assert_sp_frame_equal(result, expected)

        # method='bfill'
        result = sparse.reindex(index=range(6), method='bfill')
        expected = SparseDataFrame(data=[[11., 12., 14.],
                                         [11., 12., 14.],
                                         [21., 22., 24.],
                                         [41., 42., 44.],
                                         [41., 42., 44.],
                                         [nan, nan, nan]],
                                   index=range(6),
                                   columns=[1, 2, 4],
                                   dtype=float)
        tm.assert_sp_frame_equal(result, expected)

        # method='ffill'
        result = sparse.reindex(index=range(6), method='ffill')
        expected = SparseDataFrame(data=[[nan, nan, nan],
                                         [11., 12., 14.],
                                         [21., 22., 24.],
                                         [21., 22., 24.],
                                         [41., 42., 44.],
                                         [41., 42., 44.]],
                                   index=range(6),
                                   columns=[1, 2, 4],
                                   dtype=float)
        tm.assert_sp_frame_equal(result, expected)

        # Over columns

        # default method
        result = sparse.reindex(columns=range(6))
        expected = SparseDataFrame(data=[[nan, 11., 12., nan, 14., nan],
                                         [nan, 21., 22., nan, 24., nan],
                                         [nan, 41., 42., nan, 44., nan]],
                                   index=[1, 2, 4],
                                   columns=range(6),
                                   dtype=float)
        tm.assert_sp_frame_equal(result, expected)

        # method='bfill'
        with pytest.raises(NotImplementedError):
            sparse.reindex(columns=range(6), method='bfill')

        # method='ffill'
        with pytest.raises(NotImplementedError):
            sparse.reindex(columns=range(6), method='ffill')

    def test_take(self, float_frame):
        result = float_frame.take([1, 0, 2], axis=1)
        expected = float_frame.reindex(columns=['B', 'A', 'C'])
        tm.assert_sp_frame_equal(result, expected)

    def test_to_dense(self, float_frame, float_frame_int_kind,
                      float_frame_dense,
                      float_frame_fill0, float_frame_fill0_dense,
                      float_frame_fill2, float_frame_fill2_dense):
        def _check(frame, orig):
            dense_dm = frame.to_dense()
            tm.assert_frame_equal(frame, dense_dm)
            tm.assert_frame_equal(dense_dm, orig, check_dtype=False)

        _check(float_frame, float_frame_dense)
        _check(float_frame_int_kind, float_frame_dense)
        _check(float_frame_fill0, float_frame_fill0_dense)
        _check(float_frame_fill2, float_frame_fill2_dense)

    @pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
    def test_stack_sparse_frame(self, float_frame, float_frame_int_kind,
                                float_frame_fill0, float_frame_fill2):
        with catch_warnings(record=True):

            def _check(frame):
                dense_frame = frame.to_dense()  # noqa

                wp = Panel.from_dict({'foo': frame})
                from_dense_lp = wp.to_frame()

                from_sparse_lp = spf.stack_sparse_frame(frame)

                tm.assert_numpy_array_equal(from_dense_lp.values,
                                            from_sparse_lp.values)

            _check(float_frame)
            _check(float_frame_int_kind)

            # for now
            pytest.raises(Exception, _check, float_frame_fill0)
            pytest.raises(Exception, _check, float_frame_fill2)

    def test_transpose(self, float_frame, float_frame_int_kind,
                       float_frame_dense,
                       float_frame_fill0, float_frame_fill0_dense,
                       float_frame_fill2, float_frame_fill2_dense):

        def _check(frame, orig):
            transposed = frame.T
            untransposed = transposed.T
            tm.assert_sp_frame_equal(frame, untransposed)

            tm.assert_frame_equal(frame.T.to_dense(), orig.T)
            tm.assert_frame_equal(frame.T.T.to_dense(), orig.T.T)
            tm.assert_sp_frame_equal(frame, frame.T.T, exact_indices=False)

        _check(float_frame, float_frame_dense)
        _check(float_frame_int_kind, float_frame_dense)
        _check(float_frame_fill0, float_frame_fill0_dense)
        _check(float_frame_fill2, float_frame_fill2_dense)

    def test_shift(self, float_frame, float_frame_int_kind, float_frame_dense,
                   float_frame_fill0, float_frame_fill0_dense,
                   float_frame_fill2, float_frame_fill2_dense):

        def _check(frame, orig):
            shifted = frame.shift(0)
            exp = orig.shift(0)
            tm.assert_frame_equal(shifted.to_dense(), exp)

            shifted = frame.shift(1)
            exp = orig.shift(1)
            tm.assert_frame_equal(shifted, exp)

            shifted = frame.shift(-2)
            exp = orig.shift(-2)
            tm.assert_frame_equal(shifted, exp)

            shifted = frame.shift(2, freq='B')
            exp = orig.shift(2, freq='B')
            exp = exp.to_sparse(frame.default_fill_value,
                                kind=frame.default_kind)
            tm.assert_frame_equal(shifted, exp)

            shifted = frame.shift(2, freq=BDay())
            exp = orig.shift(2, freq=BDay())
            exp = exp.to_sparse(frame.default_fill_value,
                                kind=frame.default_kind)
            tm.assert_frame_equal(shifted, exp)

        _check(float_frame, float_frame_dense)
        _check(float_frame_int_kind, float_frame_dense)
        _check(float_frame_fill0, float_frame_fill0_dense)
        _check(float_frame_fill2, float_frame_fill2_dense)

    def test_count(self, float_frame):
        dense_result = float_frame.to_dense().count()

        result = float_frame.count()
        tm.assert_series_equal(result, dense_result)

        result = float_frame.count(axis=None)
        tm.assert_series_equal(result, dense_result)

        result = float_frame.count(axis=0)
        tm.assert_series_equal(result, dense_result)

        result = float_frame.count(axis=1)
        dense_result = float_frame.to_dense().count(axis=1)

        # win32 don't check dtype
        tm.assert_series_equal(result, dense_result, check_dtype=False)

    def test_numpy_transpose(self):
        sdf = SparseDataFrame([1, 2, 3], index=[1, 2, 3], columns=['a'])
        result = np.transpose(np.transpose(sdf))
        tm.assert_sp_frame_equal(result, sdf)

        msg = "the 'axes' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.transpose, sdf, axes=1)

    def test_combine_first(self, float_frame):
        df = float_frame

        result = df[::2].combine_first(df)
        result2 = df[::2].combine_first(df.to_dense())

        expected = df[::2].to_dense().combine_first(df.to_dense())
        expected = expected.to_sparse(fill_value=df.default_fill_value)

        tm.assert_sp_frame_equal(result, result2)
        tm.assert_sp_frame_equal(result, expected)

    def test_combine_add(self, float_frame):
        df = float_frame.to_dense()
        df2 = df.copy()
        df2['C'][:3] = np.nan
        df['A'][:3] = 5.7

        result = df.to_sparse().add(df2.to_sparse(), fill_value=0)
        expected = df.add(df2, fill_value=0).to_sparse()
        tm.assert_sp_frame_equal(result, expected)

    def test_isin(self):
        sparse_df = DataFrame({'flag': [1., 0., 1.]}).to_sparse(fill_value=0.)
        xp = sparse_df[sparse_df.flag == 1.]
        rs = sparse_df[sparse_df.flag.isin([1.])]
        tm.assert_frame_equal(xp, rs)

    def test_sparse_pow_issue(self):
        # 2220
        df = SparseDataFrame({'A': [1.1, 3.3], 'B': [2.5, -3.9]})

        # note : no error without nan
        df = SparseDataFrame({'A': [nan, 0, 1]})

        # note that 2 ** df works fine, also df ** 1
        result = 1 ** df

        r1 = result.take([0], 1)['A']
        r2 = result['A']

        assert len(r2.sp_values) == len(r1.sp_values)

    def test_as_blocks(self):
        df = SparseDataFrame({'A': [1.1, 3.3], 'B': [nan, -3.9]},
                             dtype='float64')

        # deprecated 0.21.0
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            df_blocks = df.blocks
        assert list(df_blocks.keys()) == ['float64']
        tm.assert_frame_equal(df_blocks['float64'], df)

    @pytest.mark.xfail(reason='nan column names in _init_dict problematic '
                              '(GH#16894)',
                       strict=True)
    def test_nan_columnname(self):
        # GH 8822
        nan_colname = DataFrame(Series(1.0, index=[0]), columns=[nan])
        nan_colname_sparse = nan_colname.to_sparse()
        assert np.isnan(nan_colname_sparse.columns[0])

    def test_isna(self):
        # GH 8276
        df = pd.SparseDataFrame({'A': [np.nan, np.nan, 1, 2, np.nan],
                                 'B': [0, np.nan, np.nan, 2, np.nan]})

        res = df.isna()
        exp = pd.SparseDataFrame({'A': [True, True, False, False, True],
                                  'B': [False, True, True, False, True]},
                                 default_fill_value=True)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp)

        # if fill_value is not nan, True can be included in sp_values
        df = pd.SparseDataFrame({'A': [0, 0, 1, 2, np.nan],
                                 'B': [0, np.nan, 0, 2, np.nan]},
                                default_fill_value=0.)
        res = df.isna()
        assert isinstance(res, pd.SparseDataFrame)
        exp = pd.DataFrame({'A': [False, False, False, False, True],
                            'B': [False, True, False, False, True]})
        tm.assert_frame_equal(res.to_dense(), exp)

    def test_notna(self):
        # GH 8276
        df = pd.SparseDataFrame({'A': [np.nan, np.nan, 1, 2, np.nan],
                                 'B': [0, np.nan, np.nan, 2, np.nan]})

        res = df.notna()
        exp = pd.SparseDataFrame({'A': [False, False, True, True, False],
                                  'B': [True, False, False, True, False]},
                                 default_fill_value=False)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp)

        # if fill_value is not nan, True can be included in sp_values
        df = pd.SparseDataFrame({'A': [0, 0, 1, 2, np.nan],
                                 'B': [0, np.nan, 0, 2, np.nan]},
                                default_fill_value=0.)
        res = df.notna()
        assert isinstance(res, pd.SparseDataFrame)
        exp = pd.DataFrame({'A': [True, True, True, True, False],
                            'B': [True, False, True, True, False]})
        tm.assert_frame_equal(res.to_dense(), exp)


class TestSparseDataFrameArithmetic(object):

    def test_numeric_op_scalar(self):
        df = pd.DataFrame({'A': [nan, nan, 0, 1, ],
                           'B': [0, 1, 2, nan],
                           'C': [1., 2., 3., 4.],
                           'D': [nan, nan, nan, nan]})
        sparse = df.to_sparse()

        tm.assert_sp_frame_equal(sparse + 1, (df + 1).to_sparse())

    def test_comparison_op_scalar(self):
        # GH 13001
        df = pd.DataFrame({'A': [nan, nan, 0, 1, ],
                           'B': [0, 1, 2, nan],
                           'C': [1., 2., 3., 4.],
                           'D': [nan, nan, nan, nan]})
        sparse = df.to_sparse()

        # comparison changes internal repr, compare with dense
        res = sparse > 1
        assert isinstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), df > 1)

        res = sparse != 0
        assert isinstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), df != 0)


class TestSparseDataFrameAnalytics(object):

    def test_cumsum(self, float_frame):
        expected = SparseDataFrame(float_frame.to_dense().cumsum())

        result = float_frame.cumsum()
        tm.assert_sp_frame_equal(result, expected)

        result = float_frame.cumsum(axis=None)
        tm.assert_sp_frame_equal(result, expected)

        result = float_frame.cumsum(axis=0)
        tm.assert_sp_frame_equal(result, expected)

    def test_numpy_cumsum(self, float_frame):
        result = np.cumsum(float_frame)
        expected = SparseDataFrame(float_frame.to_dense().cumsum())
        tm.assert_sp_frame_equal(result, expected)

        msg = "the 'dtype' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.cumsum,
                               float_frame, dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.cumsum,
                               float_frame, out=result)

    def test_numpy_func_call(self, float_frame):
        # no exception should be raised even though
        # numpy passes in 'axis=None' or `axis=-1'
        funcs = ['sum', 'cumsum', 'var',
                 'mean', 'prod', 'cumprod',
                 'std', 'min', 'max']
        for func in funcs:
            getattr(np, func)(float_frame)

    @pytest.mark.xfail(reason='Wrong SparseBlock initialization (GH 17386)',
                       strict=True)
    def test_quantile(self):
        # GH 17386
        data = [[1, 1], [2, 10], [3, 100], [nan, nan]]
        q = 0.1

        sparse_df = SparseDataFrame(data)
        result = sparse_df.quantile(q)

        dense_df = DataFrame(data)
        dense_expected = dense_df.quantile(q)
        sparse_expected = SparseSeries(dense_expected)

        tm.assert_series_equal(result, dense_expected)
        tm.assert_sp_series_equal(result, sparse_expected)

    @pytest.mark.xfail(reason='Wrong SparseBlock initialization (GH 17386)',
                       strict=True)
    def test_quantile_multi(self):
        # GH 17386
        data = [[1, 1], [2, 10], [3, 100], [nan, nan]]
        q = [0.1, 0.5]

        sparse_df = SparseDataFrame(data)
        result = sparse_df.quantile(q)

        dense_df = DataFrame(data)
        dense_expected = dense_df.quantile(q)
        sparse_expected = SparseDataFrame(dense_expected)

        tm.assert_frame_equal(result, dense_expected)
        tm.assert_sp_frame_equal(result, sparse_expected)

    def test_assign_with_sparse_frame(self):
        # GH 19163
        df = pd.DataFrame({"a": [1, 2, 3]})
        res = df.to_sparse(fill_value=False).assign(newcol=False)
        exp = df.assign(newcol=False).to_sparse(fill_value=False)

        tm.assert_sp_frame_equal(res, exp)

        for column in res.columns:
            assert type(res[column]) is SparseSeries
