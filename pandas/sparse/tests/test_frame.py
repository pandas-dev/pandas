# pylint: disable-msg=E1101,W0612

import operator

from numpy import nan
import numpy as np
import pandas as pd

from pandas import Series, DataFrame, bdate_range, Panel
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.offsets import BDay
import pandas.util.testing as tm
from pandas.compat import lrange
from pandas import compat
import pandas.sparse.frame as spf

from pandas._sparse import BlockIndex, IntIndex
from pandas.sparse.api import SparseSeries, SparseDataFrame, SparseArray
from pandas.tests.frame.test_misc_api import SharedWithSparse


class TestSparseDataFrame(tm.TestCase, SharedWithSparse):

    klass = SparseDataFrame
    _multiprocess_can_split_ = True

    def setUp(self):
        self.data = {'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
                     'B': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
                     'C': np.arange(10, dtype=np.float64),
                     'D': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]}

        self.dates = bdate_range('1/1/2011', periods=10)

        self.orig = pd.DataFrame(self.data, index=self.dates)
        self.iorig = pd.DataFrame(self.data, index=self.dates)

        self.frame = SparseDataFrame(self.data, index=self.dates)
        self.iframe = SparseDataFrame(self.data, index=self.dates,
                                      default_kind='integer')

        values = self.frame.values.copy()
        values[np.isnan(values)] = 0

        self.zorig = pd.DataFrame(values, columns=['A', 'B', 'C', 'D'],
                                  index=self.dates)
        self.zframe = SparseDataFrame(values, columns=['A', 'B', 'C', 'D'],
                                      default_fill_value=0, index=self.dates)

        values = self.frame.values.copy()
        values[np.isnan(values)] = 2

        self.fill_orig = pd.DataFrame(values, columns=['A', 'B', 'C', 'D'],
                                      index=self.dates)
        self.fill_frame = SparseDataFrame(values, columns=['A', 'B', 'C', 'D'],
                                          default_fill_value=2,
                                          index=self.dates)

        self.empty = SparseDataFrame()

    def test_fill_value_when_combine_const(self):
        # GH12723
        dat = np.array([0, 1, np.nan, 3, 4, 5], dtype='float')
        df = SparseDataFrame({'foo': dat}, index=range(6))

        exp = df.fillna(0).add(2)
        res = df.add(2, fill_value=0)
        tm.assert_sp_frame_equal(res, exp)

    def test_as_matrix(self):
        empty = self.empty.as_matrix()
        self.assertEqual(empty.shape, (0, 0))

        no_cols = SparseDataFrame(index=np.arange(10))
        mat = no_cols.as_matrix()
        self.assertEqual(mat.shape, (10, 0))

        no_index = SparseDataFrame(columns=np.arange(10))
        mat = no_index.as_matrix()
        self.assertEqual(mat.shape, (0, 10))

    def test_copy(self):
        cp = self.frame.copy()
        tm.assertIsInstance(cp, SparseDataFrame)
        tm.assert_sp_frame_equal(cp, self.frame)

        # as of v0.15.0
        # this is now identical (but not is_a )
        self.assertTrue(cp.index.identical(self.frame.index))

    def test_constructor(self):
        for col, series in compat.iteritems(self.frame):
            tm.assertIsInstance(series, SparseSeries)

        tm.assertIsInstance(self.iframe['A'].sp_index, IntIndex)

        # constructed zframe from matrix above
        self.assertEqual(self.zframe['A'].fill_value, 0)
        tm.assert_numpy_array_equal(pd.SparseArray([1., 2., 3., 4., 5., 6.]),
                                    self.zframe['A'].values)
        tm.assert_numpy_array_equal(np.array([0., 0., 0., 0., 1., 2.,
                                              3., 4., 5., 6.]),
                                    self.zframe['A'].to_dense().values)

        # construct no data
        sdf = SparseDataFrame(columns=np.arange(10), index=np.arange(10))
        for col, series in compat.iteritems(sdf):
            tm.assertIsInstance(series, SparseSeries)

        # construct from nested dict
        data = {}
        for c, s in compat.iteritems(self.frame):
            data[c] = s.to_dict()

        sdf = SparseDataFrame(data)
        tm.assert_sp_frame_equal(sdf, self.frame)

        # TODO: test data is copied from inputs

        # init dict with different index
        idx = self.frame.index[:5]
        cons = SparseDataFrame(
            self.frame, index=idx, columns=self.frame.columns,
            default_fill_value=self.frame.default_fill_value,
            default_kind=self.frame.default_kind, copy=True)
        reindexed = self.frame.reindex(idx)

        tm.assert_sp_frame_equal(cons, reindexed, exact_indices=False)

        # assert level parameter breaks reindex
        with tm.assertRaises(TypeError):
            self.frame.reindex(idx, level=0)

        repr(self.frame)

    def test_constructor_ndarray(self):
        # no index or columns
        sp = SparseDataFrame(self.frame.values)

        # 1d
        sp = SparseDataFrame(self.data['A'], index=self.dates, columns=['A'])
        tm.assert_sp_frame_equal(sp, self.frame.reindex(columns=['A']))

        # raise on level argument
        self.assertRaises(TypeError, self.frame.reindex, columns=['A'],
                          level=1)

        # wrong length index / columns
        with tm.assertRaisesRegexp(ValueError, "^Index length"):
            SparseDataFrame(self.frame.values, index=self.frame.index[:-1])

        with tm.assertRaisesRegexp(ValueError, "^Column length"):
            SparseDataFrame(self.frame.values, columns=self.frame.columns[:-1])

    # GH 9272
    def test_constructor_empty(self):
        sp = SparseDataFrame()
        self.assertEqual(len(sp.index), 0)
        self.assertEqual(len(sp.columns), 0)

    def test_constructor_dataframe(self):
        dense = self.frame.to_dense()
        sp = SparseDataFrame(dense)
        tm.assert_sp_frame_equal(sp, self.frame)

    def test_constructor_convert_index_once(self):
        arr = np.array([1.5, 2.5, 3.5])
        sdf = SparseDataFrame(columns=lrange(4), index=arr)
        self.assertTrue(sdf[0].index is sdf[1].index)

    def test_constructor_from_series(self):

        # GH 2873
        x = Series(np.random.randn(10000), name='a')
        x = x.to_sparse(fill_value=0)
        tm.assertIsInstance(x, SparseSeries)
        df = SparseDataFrame(x)
        tm.assertIsInstance(df, SparseDataFrame)

        x = Series(np.random.randn(10000), name='a')
        y = Series(np.random.randn(10000), name='b')
        x2 = x.astype(float)
        x2.ix[:9998] = np.NaN
        # TODO: x_sparse is unused...fix
        x_sparse = x2.to_sparse(fill_value=np.NaN)  # noqa

        # Currently fails too with weird ufunc error
        # df1 = SparseDataFrame([x_sparse, y])

        y.ix[:9998] = 0
        # TODO: y_sparse is unsused...fix
        y_sparse = y.to_sparse(fill_value=0)  # noqa
        # without sparse value raises error
        # df2 = SparseDataFrame([x2_sparse, y])

    def test_constructor_preserve_attr(self):
        # GH 13866
        arr = pd.SparseArray([1, 0, 3, 0], dtype=np.int64, fill_value=0)
        self.assertEqual(arr.dtype, np.int64)
        self.assertEqual(arr.fill_value, 0)

        df = pd.SparseDataFrame({'x': arr})
        self.assertEqual(df['x'].dtype, np.int64)
        self.assertEqual(df['x'].fill_value, 0)

        s = pd.SparseSeries(arr, name='x')
        self.assertEqual(s.dtype, np.int64)
        self.assertEqual(s.fill_value, 0)

        df = pd.SparseDataFrame(s)
        self.assertEqual(df['x'].dtype, np.int64)
        self.assertEqual(df['x'].fill_value, 0)

        df = pd.SparseDataFrame({'x': s})
        self.assertEqual(df['x'].dtype, np.int64)
        self.assertEqual(df['x'].fill_value, 0)

    def test_dtypes(self):
        df = DataFrame(np.random.randn(10000, 4))
        df.ix[:9998] = np.nan
        sdf = df.to_sparse()

        result = sdf.get_dtype_counts()
        expected = Series({'float64': 4})
        tm.assert_series_equal(result, expected)

    def test_shape(self):
        # GH 10452
        self.assertEqual(self.frame.shape, (10, 4))
        self.assertEqual(self.iframe.shape, (10, 4))
        self.assertEqual(self.zframe.shape, (10, 4))
        self.assertEqual(self.fill_frame.shape, (10, 4))

    def test_str(self):
        df = DataFrame(np.random.randn(10000, 4))
        df.ix[:9998] = np.nan

        sdf = df.to_sparse()
        str(sdf)

    def test_array_interface(self):
        res = np.sqrt(self.frame)
        dres = np.sqrt(self.frame.to_dense())
        tm.assert_frame_equal(res.to_dense(), dres)

    def test_pickle(self):

        def _test_roundtrip(frame, orig):
            result = self.round_trip_pickle(frame)
            tm.assert_sp_frame_equal(frame, result)
            tm.assert_frame_equal(result.to_dense(), orig, check_dtype=False)

        _test_roundtrip(SparseDataFrame(), DataFrame())
        self._check_all(_test_roundtrip)

    def test_dense_to_sparse(self):
        df = DataFrame({'A': [nan, nan, nan, 1, 2],
                        'B': [1, 2, nan, nan, nan]})
        sdf = df.to_sparse()
        tm.assertIsInstance(sdf, SparseDataFrame)
        self.assertTrue(np.isnan(sdf.default_fill_value))
        tm.assertIsInstance(sdf['A'].sp_index, BlockIndex)
        tm.assert_frame_equal(sdf.to_dense(), df)

        sdf = df.to_sparse(kind='integer')
        tm.assertIsInstance(sdf['A'].sp_index, IntIndex)

        df = DataFrame({'A': [0, 0, 0, 1, 2],
                        'B': [1, 2, 0, 0, 0]}, dtype=float)
        sdf = df.to_sparse(fill_value=0)
        self.assertEqual(sdf.default_fill_value, 0)
        tm.assert_frame_equal(sdf.to_dense(), df)

    def test_density(self):
        df = SparseSeries([nan, nan, nan, 0, 1, 2, 3, 4, 5, 6])
        self.assertEqual(df.density, 0.7)

        df = SparseDataFrame({'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
                              'B': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
                              'C': np.arange(10),
                              'D': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]})

        self.assertEqual(df.density, 0.75)

    def test_sparse_to_dense(self):
        pass

    def test_sparse_series_ops(self):
        self._check_frame_ops(self.frame)

    def test_sparse_series_ops_i(self):
        self._check_frame_ops(self.iframe)

    def test_sparse_series_ops_z(self):
        self._check_frame_ops(self.zframe)

    def test_sparse_series_ops_fill(self):
        self._check_frame_ops(self.fill_frame)

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
                tm.assertIsInstance(mixed_result, SparseDataFrame)
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
        result = self.frame + self.frame.ix[:, ['A', 'B']]  # noqa

    def test_op_corners(self):
        empty = self.empty + self.empty
        self.assertTrue(empty.empty)

        foo = self.frame + self.empty
        tm.assertIsInstance(foo.index, DatetimeIndex)
        tm.assert_frame_equal(foo, self.frame * np.nan)

        foo = self.empty + self.frame
        tm.assert_frame_equal(foo, self.frame * np.nan)

    def test_scalar_ops(self):
        pass

    def test_getitem(self):
        # 1585 select multiple columns
        sdf = SparseDataFrame(index=[0, 1, 2], columns=['a', 'b', 'c'])

        result = sdf[['a', 'b']]
        exp = sdf.reindex(columns=['a', 'b'])
        tm.assert_sp_frame_equal(result, exp)

        self.assertRaises(Exception, sdf.__getitem__, ['a', 'd'])

    def test_icol(self):
        # 10711 deprecated

        # 2227
        result = self.frame.iloc[:, 0]
        self.assertTrue(isinstance(result, SparseSeries))
        tm.assert_sp_series_equal(result, self.frame['A'])

        # preserve sparse index type. #2251
        data = {'A': [0, 1]}
        iframe = SparseDataFrame(data, default_kind='integer')
        self.assertEqual(type(iframe['A'].sp_index),
                         type(iframe.iloc[:, 0].sp_index))

    def test_set_value(self):

        # ok as the index gets conver to object
        frame = self.frame.copy()
        res = frame.set_value('foobar', 'B', 1.5)
        self.assertEqual(res.index.dtype, 'object')

        res = self.frame
        res.index = res.index.astype(object)

        res = self.frame.set_value('foobar', 'B', 1.5)
        self.assertIsNot(res, self.frame)
        self.assertEqual(res.index[-1], 'foobar')
        self.assertEqual(res.get_value('foobar', 'B'), 1.5)

        res2 = res.set_value('foobar', 'qux', 1.5)
        self.assertIsNot(res2, res)
        self.assert_index_equal(res2.columns,
                                pd.Index(list(self.frame.columns) + ['qux']))
        self.assertEqual(res2.get_value('foobar', 'qux'), 1.5)

    def test_fancy_index_misc(self):
        # axis = 0
        sliced = self.frame.ix[-2:, :]
        expected = self.frame.reindex(index=self.frame.index[-2:])
        tm.assert_sp_frame_equal(sliced, expected)

        # axis = 1
        sliced = self.frame.ix[:, -2:]
        expected = self.frame.reindex(columns=self.frame.columns[-2:])
        tm.assert_sp_frame_equal(sliced, expected)

    def test_getitem_overload(self):
        # slicing
        sl = self.frame[:20]
        tm.assert_sp_frame_equal(sl, self.frame.reindex(self.frame.index[:20]))

        # boolean indexing
        d = self.frame.index[5]
        indexer = self.frame.index > d

        subindex = self.frame.index[indexer]
        subframe = self.frame[indexer]

        self.assert_index_equal(subindex, subframe.index)
        self.assertRaises(Exception, self.frame.__getitem__, indexer[:-1])

    def test_setitem(self):

        def _check_frame(frame, orig):
            N = len(frame)

            # insert SparseSeries
            frame['E'] = frame['A']
            tm.assertIsInstance(frame['E'], SparseSeries)
            tm.assert_sp_series_equal(frame['E'], frame['A'],
                                      check_names=False)

            # insert SparseSeries differently-indexed
            to_insert = frame['A'][::2]
            frame['E'] = to_insert
            expected = to_insert.to_dense().reindex(frame.index)
            result = frame['E'].to_dense()
            tm.assert_series_equal(result, expected, check_names=False)
            self.assertEqual(result.name, 'E')

            # insert Series
            frame['F'] = frame['A'].to_dense()
            tm.assertIsInstance(frame['F'], SparseSeries)
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
            tm.assertIsInstance(frame['H'], SparseSeries)

            to_sparsify = np.random.randn(N)
            to_sparsify[N // 2:] = frame.default_fill_value
            frame['I'] = to_sparsify
            self.assertEqual(len(frame['I'].sp_values), N // 2)

            # insert ndarray wrong size
            self.assertRaises(Exception, frame.__setitem__, 'foo',
                              np.random.randn(N - 1))

            # scalar value
            frame['J'] = 5
            self.assertEqual(len(frame['J'].sp_values), N)
            self.assertTrue((frame['J'].sp_values == 5).all())

            frame['K'] = frame.default_fill_value
            self.assertEqual(len(frame['K'].sp_values), 0)

        self._check_all(_check_frame)

    def test_setitem_corner(self):
        self.frame['a'] = self.frame['B']
        tm.assert_sp_series_equal(self.frame['a'], self.frame['B'],
                                  check_names=False)

    def test_setitem_array(self):
        arr = self.frame['B']

        self.frame['E'] = arr
        tm.assert_sp_series_equal(self.frame['E'], self.frame['B'],
                                  check_names=False)

        self.frame['F'] = arr[:-1]
        index = self.frame.index[:-1]
        tm.assert_sp_series_equal(self.frame['E'].reindex(index),
                                  self.frame['F'].reindex(index),
                                  check_names=False)

    def test_delitem(self):
        A = self.frame['A']
        C = self.frame['C']

        del self.frame['B']
        self.assertNotIn('B', self.frame)
        tm.assert_sp_series_equal(self.frame['A'], A)
        tm.assert_sp_series_equal(self.frame['C'], C)

        del self.frame['D']
        self.assertNotIn('D', self.frame)

        del self.frame['A']
        self.assertNotIn('A', self.frame)

    def test_set_columns(self):
        self.frame.columns = self.frame.columns
        self.assertRaises(Exception, setattr, self.frame, 'columns',
                          self.frame.columns[:-1])

    def test_set_index(self):
        self.frame.index = self.frame.index
        self.assertRaises(Exception, setattr, self.frame, 'index',
                          self.frame.index[:-1])

    def test_append(self):
        a = self.frame[:5]
        b = self.frame[5:]

        appended = a.append(b)
        tm.assert_sp_frame_equal(appended, self.frame, exact_indices=False)

        a = self.frame.ix[:5, :3]
        b = self.frame.ix[5:]
        appended = a.append(b)
        tm.assert_sp_frame_equal(appended.ix[:, :3], self.frame.ix[:, :3],
                                 exact_indices=False)

    def test_apply(self):
        applied = self.frame.apply(np.sqrt)
        tm.assertIsInstance(applied, SparseDataFrame)
        tm.assert_almost_equal(applied.values, np.sqrt(self.frame.values))

        applied = self.fill_frame.apply(np.sqrt)
        self.assertEqual(applied['A'].fill_value, np.sqrt(2))

        # agg / broadcast
        broadcasted = self.frame.apply(np.sum, broadcast=True)
        tm.assertIsInstance(broadcasted, SparseDataFrame)

        exp = self.frame.to_dense().apply(np.sum, broadcast=True)
        tm.assert_frame_equal(broadcasted.to_dense(), exp)

        self.assertIs(self.empty.apply(np.sqrt), self.empty)

        from pandas.core import nanops
        applied = self.frame.apply(np.sum)
        tm.assert_series_equal(applied,
                               self.frame.to_dense().apply(nanops.nansum))

    def test_apply_nonuq(self):
        orig = DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                         index=['a', 'a', 'c'])
        sparse = orig.to_sparse()
        res = sparse.apply(lambda s: s[0], axis=1)
        exp = orig.apply(lambda s: s[0], axis=1)
        # dtype must be kept
        self.assertEqual(res.dtype, np.int64)
        # ToDo: apply must return subclassed dtype
        self.assertIsInstance(res, pd.Series)
        tm.assert_series_equal(res.to_dense(), exp)

        # df.T breaks
        sparse = orig.T.to_sparse()
        res = sparse.apply(lambda s: s[0], axis=0)  # noqa
        exp = orig.T.apply(lambda s: s[0], axis=0)
        # TODO: no non-unique columns supported in sparse yet
        # tm.assert_series_equal(res.to_dense(), exp)

    def test_applymap(self):
        # just test that it works
        result = self.frame.applymap(lambda x: x * 2)
        tm.assertIsInstance(result, SparseDataFrame)

    def test_astype(self):
        sparse = pd.SparseDataFrame({'A': SparseArray([1, 2, 3, 4],
                                                      dtype=np.int64),
                                     'B': SparseArray([4, 5, 6, 7],
                                                      dtype=np.int64)})
        self.assertEqual(sparse['A'].dtype, np.int64)
        self.assertEqual(sparse['B'].dtype, np.int64)

        res = sparse.astype(np.float64)
        exp = pd.SparseDataFrame({'A': SparseArray([1., 2., 3., 4.],
                                                   fill_value=0.),
                                  'B': SparseArray([4., 5., 6., 7.],
                                                   fill_value=0.)},
                                 default_fill_value=np.nan)
        tm.assert_sp_frame_equal(res, exp)
        self.assertEqual(res['A'].dtype, np.float64)
        self.assertEqual(res['B'].dtype, np.float64)

        sparse = pd.SparseDataFrame({'A': SparseArray([0, 2, 0, 4],
                                                      dtype=np.int64),
                                     'B': SparseArray([0, 5, 0, 7],
                                                      dtype=np.int64)},
                                    default_fill_value=0)
        self.assertEqual(sparse['A'].dtype, np.int64)
        self.assertEqual(sparse['B'].dtype, np.int64)

        res = sparse.astype(np.float64)
        exp = pd.SparseDataFrame({'A': SparseArray([0., 2., 0., 4.],
                                                   fill_value=0.),
                                  'B': SparseArray([0., 5., 0., 7.],
                                                   fill_value=0.)},
                                 default_fill_value=0.)
        tm.assert_sp_frame_equal(res, exp)
        self.assertEqual(res['A'].dtype, np.float64)
        self.assertEqual(res['B'].dtype, np.float64)

    def test_astype_bool(self):
        sparse = pd.SparseDataFrame({'A': SparseArray([0, 2, 0, 4],
                                                      fill_value=0,
                                                      dtype=np.int64),
                                     'B': SparseArray([0, 5, 0, 7],
                                                      fill_value=0,
                                                      dtype=np.int64)},
                                    default_fill_value=0)
        self.assertEqual(sparse['A'].dtype, np.int64)
        self.assertEqual(sparse['B'].dtype, np.int64)

        res = sparse.astype(bool)
        exp = pd.SparseDataFrame({'A': SparseArray([False, True, False, True],
                                                   dtype=np.bool,
                                                   fill_value=False),
                                  'B': SparseArray([False, True, False, True],
                                                   dtype=np.bool,
                                                   fill_value=False)},
                                 default_fill_value=False)
        tm.assert_sp_frame_equal(res, exp)
        self.assertEqual(res['A'].dtype, np.bool)
        self.assertEqual(res['B'].dtype, np.bool)

    def test_fillna(self):
        df = self.zframe.reindex(lrange(5))
        dense = self.zorig.reindex(lrange(5))

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

    def test_rename(self):
        # just check this works
        renamed = self.frame.rename(index=str)  # noqa
        renamed = self.frame.rename(columns=lambda x: '%s%d' % (x, len(x)))  # noqa

    def test_corr(self):
        res = self.frame.corr()
        tm.assert_frame_equal(res, self.frame.to_dense().corr())

    def test_describe(self):
        self.frame['foo'] = np.nan
        self.frame.get_dtype_counts()
        str(self.frame)
        desc = self.frame.describe()  # noqa

    def test_join(self):
        left = self.frame.ix[:, ['A', 'B']]
        right = self.frame.ix[:, ['C', 'D']]
        joined = left.join(right)
        tm.assert_sp_frame_equal(joined, self.frame, exact_indices=False)

        right = self.frame.ix[:, ['B', 'D']]
        self.assertRaises(Exception, left.join, right)

        with tm.assertRaisesRegexp(ValueError,
                                   'Other Series must have a name'):
            self.frame.join(Series(
                np.random.randn(len(self.frame)), index=self.frame.index))

    def test_reindex(self):

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
            self.assertEqual(len(length_zero), 0)
            self.assertEqual(len(length_zero.columns), len(frame.columns))
            self.assertEqual(len(length_zero['A']), 0)

            # frame being reindexed has length zero
            length_n = length_zero.reindex(index)
            self.assertEqual(len(length_n), len(frame))
            self.assertEqual(len(length_n.columns), len(frame.columns))
            self.assertEqual(len(length_n['A']), len(frame))

            # reindex columns
            reindexed = frame.reindex(columns=['A', 'B', 'Z'])
            self.assertEqual(len(reindexed.columns), 3)
            tm.assert_almost_equal(reindexed['Z'].fill_value,
                                   frame.default_fill_value)
            self.assertTrue(np.isnan(reindexed['Z'].sp_values).all())

        _check_frame(self.frame)
        _check_frame(self.iframe)
        _check_frame(self.zframe)
        _check_frame(self.fill_frame)

        # with copy=False
        reindexed = self.frame.reindex(self.frame.index, copy=False)
        reindexed['F'] = reindexed['A']
        self.assertIn('F', self.frame)

        reindexed = self.frame.reindex(self.frame.index)
        reindexed['G'] = reindexed['A']
        self.assertNotIn('G', self.frame)

    def test_reindex_fill_value(self):
        rng = bdate_range('20110110', periods=20)

        result = self.zframe.reindex(rng, fill_value=0)
        exp = self.zorig.reindex(rng, fill_value=0)
        exp = exp.to_sparse(self.zframe.default_fill_value)
        tm.assert_sp_frame_equal(result, exp)

    def test_take(self):
        result = self.frame.take([1, 0, 2], axis=1)
        expected = self.frame.reindex(columns=['B', 'A', 'C'])
        tm.assert_sp_frame_equal(result, expected)

    def test_to_dense(self):
        def _check(frame, orig):
            dense_dm = frame.to_dense()
            tm.assert_frame_equal(frame, dense_dm)
            tm.assert_frame_equal(dense_dm, orig, check_dtype=False)

        self._check_all(_check)

    def test_stack_sparse_frame(self):
        def _check(frame):
            dense_frame = frame.to_dense()  # noqa

            wp = Panel.from_dict({'foo': frame})
            from_dense_lp = wp.to_frame()

            from_sparse_lp = spf.stack_sparse_frame(frame)

            self.assert_numpy_array_equal(from_dense_lp.values,
                                          from_sparse_lp.values)

        _check(self.frame)
        _check(self.iframe)

        # for now
        self.assertRaises(Exception, _check, self.zframe)
        self.assertRaises(Exception, _check, self.fill_frame)

    def test_transpose(self):

        def _check(frame, orig):
            transposed = frame.T
            untransposed = transposed.T
            tm.assert_sp_frame_equal(frame, untransposed)

            tm.assert_frame_equal(frame.T.to_dense(), orig.T)
            tm.assert_frame_equal(frame.T.T.to_dense(), orig.T.T)
            tm.assert_sp_frame_equal(frame, frame.T.T, exact_indices=False)

        self._check_all(_check)

    def test_shift(self):

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
            exp = exp.to_sparse(frame.default_fill_value)
            tm.assert_frame_equal(shifted, exp)

            shifted = frame.shift(2, freq=BDay())
            exp = orig.shift(2, freq=BDay())
            exp = exp.to_sparse(frame.default_fill_value)
            tm.assert_frame_equal(shifted, exp)

        self._check_all(_check)

    def test_count(self):
        dense_result = self.frame.to_dense().count()

        result = self.frame.count()
        tm.assert_series_equal(result, dense_result)

        result = self.frame.count(axis=None)
        tm.assert_series_equal(result, dense_result)

        result = self.frame.count(axis=0)
        tm.assert_series_equal(result, dense_result)

        result = self.frame.count(axis=1)
        dense_result = self.frame.to_dense().count(axis=1)

        # win32 don't check dtype
        tm.assert_series_equal(result, dense_result, check_dtype=False)

    def _check_all(self, check_func):
        check_func(self.frame, self.orig)
        check_func(self.iframe, self.iorig)
        check_func(self.zframe, self.zorig)
        check_func(self.fill_frame, self.fill_orig)

    def test_numpy_transpose(self):
        sdf = SparseDataFrame([1, 2, 3], index=[1, 2, 3], columns=['a'])
        result = np.transpose(np.transpose(sdf))
        tm.assert_sp_frame_equal(result, sdf)

        msg = "the 'axes' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.transpose, sdf, axes=1)

    def test_combine_first(self):
        df = self.frame

        result = df[::2].combine_first(df)
        result2 = df[::2].combine_first(df.to_dense())

        expected = df[::2].to_dense().combine_first(df.to_dense())
        expected = expected.to_sparse(fill_value=df.default_fill_value)

        tm.assert_sp_frame_equal(result, result2)
        tm.assert_sp_frame_equal(result, expected)

    def test_combine_add(self):
        df = self.frame.to_dense()
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
        result = 1**df

        r1 = result.take([0], 1)['A']
        r2 = result['A']

        self.assertEqual(len(r2.sp_values), len(r1.sp_values))

    def test_as_blocks(self):
        df = SparseDataFrame({'A': [1.1, 3.3], 'B': [nan, -3.9]},
                             dtype='float64')

        df_blocks = df.blocks
        self.assertEqual(list(df_blocks.keys()), ['float64'])
        tm.assert_frame_equal(df_blocks['float64'], df)

    def test_nan_columnname(self):
        # GH 8822
        nan_colname = DataFrame(Series(1.0, index=[0]), columns=[nan])
        nan_colname_sparse = nan_colname.to_sparse()
        self.assertTrue(np.isnan(nan_colname_sparse.columns[0]))

    def test_isnull(self):
        # GH 8276
        df = pd.SparseDataFrame({'A': [np.nan, np.nan, 1, 2, np.nan],
                                 'B': [0, np.nan, np.nan, 2, np.nan]})

        res = df.isnull()
        exp = pd.SparseDataFrame({'A': [True, True, False, False, True],
                                  'B': [False, True, True, False, True]},
                                 default_fill_value=True)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp)

        # if fill_value is not nan, True can be included in sp_values
        df = pd.SparseDataFrame({'A': [0, 0, 1, 2, np.nan],
                                 'B': [0, np.nan, 0, 2, np.nan]},
                                default_fill_value=0.)
        res = df.isnull()
        tm.assertIsInstance(res, pd.SparseDataFrame)
        exp = pd.DataFrame({'A': [False, False, False, False, True],
                            'B': [False, True, False, False, True]})
        tm.assert_frame_equal(res.to_dense(), exp)

    def test_isnotnull(self):
        # GH 8276
        df = pd.SparseDataFrame({'A': [np.nan, np.nan, 1, 2, np.nan],
                                 'B': [0, np.nan, np.nan, 2, np.nan]})

        res = df.isnotnull()
        exp = pd.SparseDataFrame({'A': [False, False, True, True, False],
                                  'B': [True, False, False, True, False]},
                                 default_fill_value=False)
        exp._default_fill_value = np.nan
        tm.assert_sp_frame_equal(res, exp)

        # if fill_value is not nan, True can be included in sp_values
        df = pd.SparseDataFrame({'A': [0, 0, 1, 2, np.nan],
                                 'B': [0, np.nan, 0, 2, np.nan]},
                                default_fill_value=0.)
        res = df.isnotnull()
        tm.assertIsInstance(res, pd.SparseDataFrame)
        exp = pd.DataFrame({'A': [True, True, True, True, False],
                            'B': [True, False, True, True, False]})
        tm.assert_frame_equal(res.to_dense(), exp)


class TestSparseDataFrameArithmetic(tm.TestCase):

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
        tm.assertIsInstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), df > 1)

        res = sparse != 0
        tm.assertIsInstance(res, pd.SparseDataFrame)
        tm.assert_frame_equal(res.to_dense(), df != 0)


class TestSparseDataFrameAnalytics(tm.TestCase):
    def setUp(self):
        self.data = {'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
                     'B': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
                     'C': np.arange(10, dtype=float),
                     'D': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]}

        self.dates = bdate_range('1/1/2011', periods=10)

        self.frame = SparseDataFrame(self.data, index=self.dates)

    def test_cumsum(self):
        expected = SparseDataFrame(self.frame.to_dense().cumsum())

        result = self.frame.cumsum()
        tm.assert_sp_frame_equal(result, expected)

        result = self.frame.cumsum(axis=None)
        tm.assert_sp_frame_equal(result, expected)

        result = self.frame.cumsum(axis=0)
        tm.assert_sp_frame_equal(result, expected)

    def test_numpy_cumsum(self):
        result = np.cumsum(self.frame)
        expected = SparseDataFrame(self.frame.to_dense().cumsum())
        tm.assert_sp_frame_equal(result, expected)

        msg = "the 'dtype' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.cumsum,
                              self.frame, dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.cumsum,
                              self.frame, out=result)

    def test_numpy_func_call(self):
        # no exception should be raised even though
        # numpy passes in 'axis=None' or `axis=-1'
        funcs = ['sum', 'cumsum', 'var',
                 'mean', 'prod', 'cumprod',
                 'std', 'min', 'max']
        for func in funcs:
            getattr(np, func)(self.frame)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
