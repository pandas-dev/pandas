# -*- coding: utf-8 -*-
# pylint: disable-msg=W0612,E1101
import sys
import nose
import itertools
import warnings
from datetime import datetime

from pandas.compat import range, lrange, lzip, StringIO, lmap, map
from pandas.tslib import NaT
from numpy import nan
from numpy.random import randn
import numpy as np

import pandas as pd
import pandas.core.common as com
from pandas import option_context
from pandas.core.indexing import _non_reducing_slice, _maybe_numeric_slice
from pandas.core.api import (DataFrame, Index, Series, Panel, isnull,
                             MultiIndex, Float64Index, Timestamp, Timedelta)
from pandas.util.testing import (assert_almost_equal, assert_series_equal,
                                 assert_frame_equal, assert_panel_equal,
                                 assert_attr_equal)
from pandas import concat, lib
from pandas.io.common import PerformanceWarning

import pandas.util.testing as tm
from pandas import date_range
from numpy.testing.decorators import slow

_verbose = False

#-------------------------------------------------------------------------------
# Indexing test cases


def _generate_indices(f, values=False):
    """ generate the indicies
          if values is True , use the axis values
                    is False, use the range
                    """

    axes = f.axes
    if values:
        axes = [ lrange(len(a)) for a in axes ]

    return itertools.product(*axes)

def _get_value(f, i, values=False):
    """ return the value for the location i """

    # check agains values
    if values:
        return f.values[i]

    # this is equiv of f[col][row].....
    #v = f
    #for a in reversed(i):
    #    v = v.__getitem__(a)
    #return v
    return f.ix[i]

def _get_result(obj, method, key, axis):
    """ return the result for this obj with this key and this axis """

    if isinstance(key, dict):
        key = key[axis]

    # use an artifical conversion to map the key as integers to the labels
    # so ix can work for comparisions
    if method == 'indexer':
        method = 'ix'
        key    = obj._get_axis(axis)[key]

    # in case we actually want 0 index slicing
    try:
        xp  = getattr(obj, method).__getitem__(_axify(obj,key,axis))
    except:
        xp  = getattr(obj, method).__getitem__(key)

    return xp

def _axify(obj, key, axis):
    # create a tuple accessor
    if axis is not None:
        axes = [ slice(None) ] * obj.ndim
        axes[axis] = key
        return tuple(axes)
    return k


def _mklbl(prefix,n):
    return ["%s%s" % (prefix,i)  for i in range(n)]

class TestIndexing(tm.TestCase):

    _multiprocess_can_split_ = True

    _objs = set(['series','frame','panel'])
    _typs = set(['ints','labels','mixed','ts','floats','empty'])

    def setUp(self):

        import warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)

        self.series_ints = Series(np.random.rand(4), index=lrange(0,8,2))
        self.frame_ints = DataFrame(np.random.randn(4, 4), index=lrange(0, 8, 2), columns=lrange(0,12,3))
        self.panel_ints = Panel(np.random.rand(4,4,4), items=lrange(0,8,2),major_axis=lrange(0,12,3),minor_axis=lrange(0,16,4))

        self.series_labels = Series(np.random.randn(4), index=list('abcd'))
        self.frame_labels  = DataFrame(np.random.randn(4, 4), index=list('abcd'), columns=list('ABCD'))
        self.panel_labels  = Panel(np.random.randn(4,4,4), items=list('abcd'), major_axis=list('ABCD'), minor_axis=list('ZYXW'))

        self.series_mixed  = Series(np.random.randn(4), index=[2, 4, 'null', 8])
        self.frame_mixed   = DataFrame(np.random.randn(4, 4), index=[2, 4, 'null', 8])
        self.panel_mixed   = Panel(np.random.randn(4,4,4), items=[2,4,'null',8])

        self.series_ts     = Series(np.random.randn(4), index=date_range('20130101', periods=4))
        self.frame_ts      = DataFrame(np.random.randn(4, 4), index=date_range('20130101', periods=4))
        self.panel_ts      = Panel(np.random.randn(4, 4, 4), items=date_range('20130101', periods=4))

        #self.series_floats = Series(np.random.randn(4), index=[1.00, 2.00, 3.00, 4.00])
        #self.frame_floats  = DataFrame(np.random.randn(4, 4), columns=[1.00, 2.00, 3.00, 4.00])
        #self.panel_floats  = Panel(np.random.rand(4,4,4), items = [1.00,2.00,3.00,4.00])

        self.frame_empty   = DataFrame({})
        self.series_empty  = Series({})
        self.panel_empty   = Panel({})

        # form agglomerates
        for o in self._objs:

            d = dict()
            for t in self._typs:
                d[t] = getattr(self,'%s_%s' % (o,t),None)

            setattr(self,o,d)

    def check_values(self, f, func, values = False):

        if f is None: return
        axes = f.axes
        indicies = itertools.product(*axes)

        for i in indicies:
            result = getattr(f,func)[i]

            # check agains values
            if values:
                expected = f.values[i]
            else:
                expected = f
                for a in reversed(i):
                    expected = expected.__getitem__(a)

            assert_almost_equal(result, expected)


    def check_result(self, name, method1, key1, method2, key2, typs = None, objs = None, axes = None, fails = None):


        def _eq(t, o, a, obj, k1, k2):
            """ compare equal for these 2 keys """

            if a is not None and a > obj.ndim-1:
                return

            def _print(result, error = None):
                if error is not None:
                    error = str(error)
                v = "%-16.16s [%-16.16s]: [typ->%-8.8s,obj->%-8.8s,key1->(%-4.4s),key2->(%-4.4s),axis->%s] %s" % (name,result,t,o,method1,method2,a,error or '')
                if _verbose:
                    com.pprint_thing(v)

            try:

                ### good debug location ###
                #if name == 'bool' and t == 'empty' and o == 'series' and method1 == 'loc':
                #    import pdb; pdb.set_trace()

                rs  = getattr(obj, method1).__getitem__(_axify(obj,k1,a))

                try:
                    xp = _get_result(obj,method2,k2,a)
                except:
                    result = 'no comp'
                    _print(result)
                    return

                try:
                    if np.isscalar(rs) and np.isscalar(xp):
                        self.assertEqual(rs, xp)
                    elif xp.ndim == 1:
                        assert_series_equal(rs,xp)
                    elif xp.ndim == 2:
                        assert_frame_equal(rs,xp)
                    elif xp.ndim == 3:
                        assert_panel_equal(rs,xp)
                    result = 'ok'
                except (AssertionError):
                    result = 'fail'

                # reverse the checks
                if fails is True:
                    if result == 'fail':
                        result = 'ok (fail)'

                if not result.startswith('ok'):
                    raise AssertionError(_print(result))

                _print(result)

            except AssertionError:
                raise
            except Exception as detail:

                # if we are in fails, the ok, otherwise raise it
                if fails is not None:
                    if isinstance(detail, fails):
                        result = 'ok (%s)' % type(detail).__name__
                        _print(result)
                        return

                result = type(detail).__name__
                raise AssertionError(_print(result, error = detail))

        if typs is None:
            typs = self._typs

        if objs is None:
            objs = self._objs

        if axes is not None:
            if not isinstance(axes,(tuple,list)):
                axes = [ axes ]
            else:
                axes = list(axes)
        else:
            axes = [ 0, 1, 2]

        # check
        for o in objs:
            if o not in self._objs:
                continue

            d = getattr(self,o)
            for a in axes:
                for t in typs:
                    if t not in self._typs:
                        continue

                    obj = d[t]
                    if obj is not None:
                        obj = obj.copy()

                        k2 = key2
                        _eq(t, o, a, obj, key1, k2)

    def test_indexer_caching(self):
        # GH5727
        # make sure that indexers are in the _internal_names_set
        n = 1000001
        arrays = [lrange(n), lrange(n)]
        index = MultiIndex.from_tuples(lzip(*arrays))
        s = Series(np.zeros(n), index=index)
        str(s)

        # setitem
        expected = Series(np.ones(n), index=index)
        s = Series(np.zeros(n), index=index)
        s[s==0] = 1
        assert_series_equal(s,expected)

    def test_at_and_iat_get(self):

        def _check(f, func, values = False):

            if f is not None:
                indicies = _generate_indices(f, values)
                for i in indicies:
                    result = getattr(f,func)[i]
                    expected = _get_value(f,i,values)
                    assert_almost_equal(result, expected)

        for o in self._objs:

            d = getattr(self,o)

            # iat
            _check(d['ints'],'iat', values=True)
            for f in [d['labels'],d['ts'],d['floats']]:
                if f is not None:
                    self.assertRaises(ValueError, self.check_values, f, 'iat')

            # at
            _check(d['ints'],  'at')
            _check(d['labels'],'at')
            _check(d['ts'],    'at')
            _check(d['floats'],'at')

    def test_at_and_iat_set(self):

        def _check(f, func, values = False):

            if f is not None:
                indicies = _generate_indices(f, values)
                for i in indicies:
                    getattr(f,func)[i] = 1
                    expected = _get_value(f,i,values)
                    assert_almost_equal(expected, 1)

        for t in self._objs:

            d = getattr(self,t)

            _check(d['ints'],'iat',values=True)
            for f in [d['labels'],d['ts'],d['floats']]:
                if f is not None:
                    self.assertRaises(ValueError, _check, f, 'iat')

            # at
            _check(d['ints'],  'at')
            _check(d['labels'],'at')
            _check(d['ts'],    'at')
            _check(d['floats'],'at')

    def test_at_iat_coercion(self):

        # as timestamp is not a tuple!
        dates = date_range('1/1/2000', periods=8)
        df = DataFrame(randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])
        s = df['A']

        result = s.at[dates[5]]
        xp     = s.values[5]
        self.assertEqual(result, xp)

        # GH 7729
        # make sure we are boxing the returns
        s = Series(['2014-01-01', '2014-02-02'], dtype='datetime64[ns]')
        expected = Timestamp('2014-02-02')

        for r in [ lambda : s.iat[1], lambda : s.iloc[1] ]:
            result = r()
            self.assertEqual(result, expected)

        s = Series(['1 days','2 days'], dtype='timedelta64[ns]')
        expected = Timedelta('2 days')

        for r in [ lambda : s.iat[1], lambda : s.iloc[1] ]:
            result = r()
            self.assertEqual(result, expected)

    def test_iat_invalid_args(self):
        pass

    def test_imethods_with_dups(self):

        # GH6493
        # iat/iloc with dups

        s = Series(range(5), index=[1,1,2,2,3], dtype='int64')
        result = s.iloc[2]
        self.assertEqual(result,2)
        result = s.iat[2]
        self.assertEqual(result,2)

        self.assertRaises(IndexError, lambda : s.iat[10])
        self.assertRaises(IndexError, lambda : s.iat[-10])

        result = s.iloc[[2,3]]
        expected = Series([2,3],[2,2],dtype='int64')
        assert_series_equal(result,expected)

        df = s.to_frame()
        result = df.iloc[2]
        expected = Series(2, index=[0], name=2)
        assert_series_equal(result, expected)

        result = df.iat[2,0]
        expected = 2
        self.assertEqual(result,2)

    def test_repeated_getitem_dups(self):
        # GH 5678
        # repeated gettitems on a dup index returing a ndarray
        df = DataFrame(np.random.random_sample((20,5)), index=['ABCDE'[x%5] for x in range(20)])
        expected = df.loc['A',0]
        result = df.loc[:,0].loc['A']
        assert_series_equal(result,expected)

    def test_iloc_exceeds_bounds(self):

        # GH6296
        # iloc should allow indexers that exceed the bounds
        df = DataFrame(np.random.random_sample((20,5)), columns=list('ABCDE'))
        expected = df

        # lists of positions should raise IndexErrror!
        with tm.assertRaisesRegexp(IndexError, 'positional indexers are out-of-bounds'):
            df.iloc[:,[0,1,2,3,4,5]]
        self.assertRaises(IndexError, lambda : df.iloc[[1,30]])
        self.assertRaises(IndexError, lambda : df.iloc[[1,-30]])
        self.assertRaises(IndexError, lambda : df.iloc[[100]])

        s = df['A']
        self.assertRaises(IndexError, lambda : s.iloc[[100]])
        self.assertRaises(IndexError, lambda : s.iloc[[-100]])

        # still raise on a single indexer
        with tm.assertRaisesRegexp(IndexError, 'single positional indexer is out-of-bounds'):
            df.iloc[30]
        self.assertRaises(IndexError, lambda : df.iloc[-30])

        # GH10779
        # single positive/negative indexer exceeding Series bounds should raise an IndexError
        with tm.assertRaisesRegexp(IndexError, 'single positional indexer is out-of-bounds'):
            s.iloc[30]
        self.assertRaises(IndexError, lambda : s.iloc[-30])

        # slices are ok
        result = df.iloc[:,4:10]  # 0 < start < len < stop
        expected = df.iloc[:,4:]
        assert_frame_equal(result,expected)

        result = df.iloc[:,-4:-10]  # stop < 0 < start < len
        expected = df.iloc[:,:0]
        assert_frame_equal(result,expected)

        result = df.iloc[:,10:4:-1]  # 0 < stop < len < start (down)
        expected = df.iloc[:,:4:-1]
        assert_frame_equal(result,expected)

        result = df.iloc[:,4:-10:-1]  # stop < 0 < start < len (down)
        expected = df.iloc[:,4::-1]
        assert_frame_equal(result,expected)

        result = df.iloc[:,-10:4]  # start < 0 < stop < len
        expected = df.iloc[:,:4]
        assert_frame_equal(result,expected)

        result = df.iloc[:,10:4]  # 0 < stop < len < start
        expected = df.iloc[:,:0]
        assert_frame_equal(result,expected)

        result = df.iloc[:,-10:-11:-1]  # stop < start < 0 < len (down)
        expected = df.iloc[:,:0]
        assert_frame_equal(result,expected)

        result = df.iloc[:,10:11]  # 0 < len < start < stop
        expected = df.iloc[:,:0]
        assert_frame_equal(result,expected)

        # slice bounds exceeding is ok
        result = s.iloc[18:30]
        expected = s.iloc[18:]
        assert_series_equal(result,expected)

        result = s.iloc[30:]
        expected = s.iloc[:0]
        assert_series_equal(result,expected)

        result = s.iloc[30::-1]
        expected = s.iloc[::-1]
        assert_series_equal(result,expected)

        # doc example
        def check(result,expected):
            str(result)
            result.dtypes
            assert_frame_equal(result,expected)

        dfl = DataFrame(np.random.randn(5,2),columns=list('AB'))
        check(dfl.iloc[:,2:3],DataFrame(index=dfl.index))
        check(dfl.iloc[:,1:3],dfl.iloc[:,[1]])
        check(dfl.iloc[4:6],dfl.iloc[[4]])

        self.assertRaises(IndexError, lambda : dfl.iloc[[4,5,6]])
        self.assertRaises(IndexError, lambda : dfl.iloc[:,4])

    def test_iloc_getitem_int(self):

        # integer
        self.check_result('integer', 'iloc', 2, 'ix', { 0 : 4, 1: 6, 2: 8 }, typs = ['ints'])
        self.check_result('integer', 'iloc', 2, 'indexer', 2, typs = ['labels','mixed','ts','floats','empty'], fails = IndexError)

    def test_iloc_getitem_neg_int(self):

        # neg integer
        self.check_result('neg int', 'iloc', -1, 'ix', { 0 : 6, 1: 9, 2: 12 }, typs = ['ints'])
        self.check_result('neg int', 'iloc', -1, 'indexer', -1, typs = ['labels','mixed','ts','floats','empty'], fails = IndexError)

    def test_iloc_getitem_list_int(self):

        # list of ints
        self.check_result('list int', 'iloc', [0,1,2], 'ix', { 0 : [0,2,4], 1 : [0,3,6], 2: [0,4,8] }, typs = ['ints'])
        self.check_result('list int', 'iloc', [2], 'ix', { 0 : [4], 1 : [6], 2: [8] }, typs = ['ints'])
        self.check_result('list int', 'iloc', [0,1,2], 'indexer', [0,1,2], typs = ['labels','mixed','ts','floats','empty'], fails = IndexError)

        # array of ints
        # (GH5006), make sure that a single indexer is returning the correct type
        self.check_result('array int', 'iloc', np.array([0,1,2]), 'ix', { 0 : [0,2,4], 1 : [0,3,6], 2: [0,4,8] }, typs = ['ints'])
        self.check_result('array int', 'iloc', np.array([2]), 'ix', { 0 : [4], 1 : [6], 2: [8] }, typs = ['ints'])
        self.check_result('array int', 'iloc', np.array([0,1,2]), 'indexer', [0,1,2], typs = ['labels','mixed','ts','floats','empty'], fails = IndexError)

    def test_iloc_getitem_neg_int_can_reach_first_index(self):
        # GH10547 and GH10779
        # negative integers should be able to reach index 0
        df = DataFrame({'A': [2, 3, 5], 'B': [7, 11, 13]})
        s = df['A']

        expected = df.iloc[0]
        result = df.iloc[-3]
        assert_series_equal(result, expected)

        expected = df.iloc[[0]]
        result = df.iloc[[-3]]
        assert_frame_equal(result, expected)

        expected = s.iloc[0]
        result = s.iloc[-3]
        self.assertEqual(result, expected)

        expected = s.iloc[[0]]
        result = s.iloc[[-3]]
        assert_series_equal(result, expected)

        # check the length 1 Series case highlighted in GH10547
        expected = pd.Series(['a'], index=['A'])
        result = expected.iloc[[-1]]
        assert_series_equal(result, expected)

    def test_iloc_getitem_dups(self):

        # no dups in panel (bug?)
        self.check_result('list int (dups)', 'iloc', [0,1,1,3], 'ix', { 0 : [0,2,2,6], 1 : [0,3,3,9] }, objs = ['series','frame'], typs = ['ints'])

        # GH 6766
        df1 = DataFrame([{'A':None, 'B':1},{'A':2, 'B':2}])
        df2 = DataFrame([{'A':3, 'B':3},{'A':4, 'B':4}])
        df = concat([df1, df2], axis=1)

        # cross-sectional indexing
        result = df.iloc[0,0]
        self.assertTrue(isnull(result))

        result = df.iloc[0,:]
        expected = Series([np.nan, 1, 3, 3], index=['A','B','A','B'], name=0)
        assert_series_equal(result,expected)

    def test_iloc_getitem_array(self):

        # array like
        s = Series(index=lrange(1,4))
        self.check_result('array like', 'iloc', s.index, 'ix', { 0 : [2,4,6], 1 : [3,6,9], 2: [4,8,12] }, typs = ['ints'])

    def test_iloc_getitem_bool(self):

        # boolean indexers
        b = [True,False,True,False,]
        self.check_result('bool', 'iloc', b, 'ix', b, typs = ['ints'])
        self.check_result('bool', 'iloc', b, 'ix', b, typs = ['labels','mixed','ts','floats','empty'], fails = IndexError)

    def test_iloc_getitem_slice(self):

        # slices
        self.check_result('slice', 'iloc', slice(1,3), 'ix', { 0 : [2,4], 1: [3,6], 2: [4,8] }, typs = ['ints'])
        self.check_result('slice', 'iloc', slice(1,3), 'indexer', slice(1,3), typs = ['labels','mixed','ts','floats','empty'], fails = IndexError)

    def test_iloc_getitem_slice_dups(self):

        df1 = DataFrame(np.random.randn(10,4),columns=['A','A','B','B'])
        df2 = DataFrame(np.random.randint(0,10,size=20).reshape(10,2),columns=['A','C'])

        # axis=1
        df = concat([df1,df2],axis=1)
        assert_frame_equal(df.iloc[:,:4],df1)
        assert_frame_equal(df.iloc[:,4:],df2)

        df = concat([df2,df1],axis=1)
        assert_frame_equal(df.iloc[:,:2],df2)
        assert_frame_equal(df.iloc[:,2:],df1)

        assert_frame_equal(df.iloc[:,0:3],concat([df2,df1.iloc[:,[0]]],axis=1))

        # axis=0
        df = concat([df,df],axis=0)
        assert_frame_equal(df.iloc[0:10,:2],df2)
        assert_frame_equal(df.iloc[0:10,2:],df1)
        assert_frame_equal(df.iloc[10:,:2],df2)
        assert_frame_equal(df.iloc[10:,2:],df1)

    def test_iloc_getitem_multiindex(self):

        arr = np.random.randn(3, 3)
        df = DataFrame(arr,
                       columns=[[2,2,4],[6,8,10]],
                       index=[[4,4,8],[8,10,12]])

        rs = df.iloc[2]
        xp = Series(arr[2],index=df.columns)
        assert_series_equal(rs, xp)

        rs = df.iloc[:,2]
        xp = Series(arr[:, 2],index=df.index)
        assert_series_equal(rs, xp)

        rs = df.iloc[2,2]
        xp = df.values[2,2]
        self.assertEqual(rs, xp)

        # for multiple items
        # GH 5528
        rs = df.iloc[[0,1]]
        xp = df.xs(4,drop_level=False)
        assert_frame_equal(rs,xp)

        tup = zip(*[['a','a','b','b'],['x','y','x','y']])
        index = MultiIndex.from_tuples(tup)
        df = DataFrame(np.random.randn(4, 4), index=index)
        rs = df.iloc[[2, 3]]
        xp = df.xs('b',drop_level=False)
        assert_frame_equal(rs,xp)

    def test_iloc_setitem(self):
        df = self.frame_ints

        df.iloc[1,1] = 1
        result = df.iloc[1,1]
        self.assertEqual(result, 1)

        df.iloc[:,2:3] = 0
        expected = df.iloc[:,2:3]
        result = df.iloc[:,2:3]
        assert_frame_equal(result, expected)

        # GH5771
        s = Series(0,index=[4,5,6])
        s.iloc[1:2] += 1
        expected = Series([0,1,0],index=[4,5,6])
        assert_series_equal(s, expected)

    def test_ix_loc_setitem_consistency(self):

        # GH 5771
        # loc with slice and series
        s = Series(0,index=[4,5,6])
        s.loc[4:5] += 1
        expected = Series([1,1,0],index=[4,5,6])
        assert_series_equal(s, expected)

        # GH 5928
        # chained indexing assignment
        df = DataFrame({'a' : [0,1,2] })
        expected = df.copy()
        expected.ix[[0,1,2],'a'] = -expected.ix[[0,1,2],'a']

        df['a'].ix[[0,1,2]] = -df['a'].ix[[0,1,2]]
        assert_frame_equal(df,expected)

        df = DataFrame({'a' : [0,1,2], 'b' :[0,1,2] })
        df['a'].ix[[0,1,2]] = -df['a'].ix[[0,1,2]].astype('float64') + 0.5
        expected = DataFrame({'a' : [0.5,-0.5,-1.5], 'b' : [0,1,2] })
        assert_frame_equal(df,expected)

        # GH 8607
        # ix setitem consistency
        df = DataFrame(
            {'timestamp':[1413840976, 1413842580, 1413760580],
             'delta':[1174, 904, 161],
             'elapsed':[7673, 9277, 1470]
             })
        expected = DataFrame(
            {'timestamp':pd.to_datetime([1413840976, 1413842580, 1413760580], unit='s'),
             'delta':[1174, 904, 161],
             'elapsed':[7673, 9277, 1470]
             })

        df2 = df.copy()
        df2['timestamp'] = pd.to_datetime(df['timestamp'], unit='s')
        assert_frame_equal(df2,expected)

        df2 = df.copy()
        df2.loc[:,'timestamp'] = pd.to_datetime(df['timestamp'], unit='s')
        assert_frame_equal(df2,expected)

        df2 = df.copy()
        df2.ix[:,2] = pd.to_datetime(df['timestamp'], unit='s')
        assert_frame_equal(df2,expected)

    def test_ix_loc_consistency(self):

        # GH 8613
        # some edge cases where ix/loc should return the same
        # this is not an exhaustive case

        def compare(result, expected):
            if lib.isscalar(expected):
                self.assertEqual(result, expected)
            else:
                self.assertTrue(expected.equals(result))

        # failure cases for .loc, but these work for .ix
        df = pd.DataFrame(np.random.randn(5,4), columns=list('ABCD'))
        for key in [ slice(1,3), tuple([slice(0,2),slice(0,2)]), tuple([slice(0,2),df.columns[0:2]]) ]:

            for index in [ tm.makeStringIndex, tm.makeUnicodeIndex,
                           tm.makeDateIndex, tm.makePeriodIndex, tm.makeTimedeltaIndex ]:
                df.index = index(len(df.index))
                df.ix[key]

                self.assertRaises(TypeError, lambda : df.loc[key])

        df = pd.DataFrame(np.random.randn(5,4), columns=list('ABCD'), index=pd.date_range('2012-01-01', periods=5))

        for key in [ '2012-01-03',
                     '2012-01-31',
                     slice('2012-01-03','2012-01-03'),
                     slice('2012-01-03','2012-01-04'),
                     slice('2012-01-03','2012-01-06',2),
                     slice('2012-01-03','2012-01-31'),
                     tuple([[True,True,True,False,True]]),
                     ]:

            # getitem

            # if the expected raises, then compare the exceptions
            try:
                expected = df.ix[key]
            except KeyError:
                self.assertRaises(KeyError, lambda : df.loc[key])
                continue

            result = df.loc[key]
            compare(result, expected)

            # setitem
            df1 = df.copy()
            df2 = df.copy()

            df1.ix[key] = 10
            df2.loc[key] = 10
            compare(df2, df1)

        # edge cases
        s = Series([1,2,3,4], index=list('abde'))

        result1 = s['a':'c']
        result2 = s.ix['a':'c']
        result3 = s.loc['a':'c']
        assert_series_equal(result1,result2)
        assert_series_equal(result1,result3)

        # now work rather than raising KeyError
        s = Series(range(5),[-2,-1,1,2,3])

        result1 = s.ix[-10:3]
        result2 = s.loc[-10:3]
        assert_series_equal(result1,result2)

        result1 = s.ix[0:3]
        result2 = s.loc[0:3]
        assert_series_equal(result1,result2)

    def test_setitem_multiindex(self):
        for index_fn in ('ix', 'loc'):
            def check(target, indexers, value, compare_fn, expected=None):
                fn = getattr(target, index_fn)
                fn.__setitem__(indexers, value)
                result = fn.__getitem__(indexers)
                if expected is None:
                    expected = value
                compare_fn(result, expected)
            # GH7190
            index = pd.MultiIndex.from_product([np.arange(0,100), np.arange(0, 80)], names=['time', 'firm'])
            t, n = 0, 2
            df = DataFrame(np.nan,columns=['A', 'w', 'l', 'a', 'x', 'X', 'd', 'profit'], index=index)
            check(
                target=df, indexers=((t,n), 'X'),
                value=0, compare_fn=self.assertEqual
            )

            df = DataFrame(-999,columns=['A', 'w', 'l', 'a', 'x', 'X', 'd', 'profit'], index=index)
            check(
                target=df, indexers=((t,n), 'X'),
                value=1, compare_fn=self.assertEqual
            )

            df = DataFrame(columns=['A', 'w', 'l', 'a', 'x', 'X', 'd', 'profit'], index=index)
            check(
                target=df, indexers=((t,n), 'X'),
                value=2, compare_fn=self.assertEqual
            )

            # GH 7218, assinging with 0-dim arrays
            df = DataFrame(-999,columns=['A', 'w', 'l', 'a', 'x', 'X', 'd', 'profit'], index=index)
            check(
                target=df, indexers=((t,n), 'X'),
                value=np.array(3), compare_fn=self.assertEqual,
                expected=3,
            )

            # GH5206
            df = pd.DataFrame(
                np.arange(25).reshape(5, 5), columns='A,B,C,D,E'.split(','),
                dtype=float
            )
            df['F'] = 99
            row_selection = df['A'] % 2 == 0
            col_selection = ['B', 'C']
            df.ix[row_selection, col_selection] = df['F']
            output = pd.DataFrame(99., index=[0, 2, 4], columns=['B', 'C'])
            assert_frame_equal(df.ix[row_selection, col_selection], output)
            check(
                target=df, indexers=(row_selection, col_selection),
                value=df['F'], compare_fn=assert_frame_equal,
                expected=output,
            )

            # GH11372
            idx = pd.MultiIndex.from_product([
                ['A', 'B', 'C'],
                pd.date_range('2015-01-01', '2015-04-01', freq='MS')
            ])
            cols = pd.MultiIndex.from_product([
                ['foo', 'bar'],
                pd.date_range('2016-01-01', '2016-02-01', freq='MS')
            ])
            df = pd.DataFrame(np.random.random((12, 4)), index=idx, columns=cols)
            subidx = pd.MultiIndex.from_tuples(
                [('A', pd.Timestamp('2015-01-01')), ('A', pd.Timestamp('2015-02-01'))]
            )
            subcols = pd.MultiIndex.from_tuples(
                [('foo', pd.Timestamp('2016-01-01')), ('foo', pd.Timestamp('2016-02-01'))]
            )
            vals = pd.DataFrame(np.random.random((2, 2)), index=subidx, columns=subcols)
            check(
                target=df, indexers=(subidx, subcols),
                value=vals, compare_fn=assert_frame_equal,
            )
            # set all columns
            vals = pd.DataFrame(np.random.random((2, 4)), index=subidx, columns=cols)
            check(
                target=df, indexers=(subidx, slice(None, None, None)),
                value=vals, compare_fn=assert_frame_equal,
            )
            # identity
            copy = df.copy()
            check(
                target=df, indexers=(df.index, df.columns),
                value=df, compare_fn=assert_frame_equal,
                expected=copy
            )

    def test_indexing_with_datetime_tz(self):

        # 8260
        # support datetime64 with tz

        idx = Index(date_range('20130101',periods=3,tz='US/Eastern'),
                    name='foo')
        dr = date_range('20130110',periods=3)
        df = DataFrame({'A' : idx, 'B' : dr})
        df['C'] = idx
        df.iloc[1,1] = pd.NaT
        df.iloc[1,2] = pd.NaT

        # indexing
        result = df.iloc[1]
        expected = Series([Timestamp('2013-01-02 00:00:00-0500', tz='US/Eastern'), np.nan, np.nan],
                          index=list('ABC'), dtype='object', name=1)
        assert_series_equal(result, expected)
        result = df.loc[1]
        expected = Series([Timestamp('2013-01-02 00:00:00-0500', tz='US/Eastern'), np.nan, np.nan],
                          index=list('ABC'), dtype='object', name=1)
        assert_series_equal(result, expected)

        # indexing - fast_xs
        df = DataFrame({'a': date_range('2014-01-01', periods=10, tz='UTC')})
        result = df.iloc[5]
        expected = Timestamp('2014-01-06 00:00:00+0000', tz='UTC', offset='D')
        self.assertEqual(result, expected)

        result = df.loc[5]
        self.assertEqual(result, expected)

        # indexing - boolean
        result = df[df.a > df.a[3]]
        expected = df.iloc[4:]
        assert_frame_equal(result, expected)

        # indexing - setting an element
        df = DataFrame( data = pd.to_datetime(['2015-03-30 20:12:32','2015-03-12 00:11:11']) ,columns=['time'] )
        df['new_col']=['new','old']
        df.time=df.set_index('time').index.tz_localize('UTC')
        v = df[df.new_col=='new'].set_index('time').index.tz_convert('US/Pacific')

        # trying to set a single element on a part of a different timezone
        def f():
            df.loc[df.new_col=='new','time'] = v
        self.assertRaises(ValueError, f)

        v = df.loc[df.new_col=='new','time'] + pd.Timedelta('1s')
        df.loc[df.new_col=='new','time'] = v
        assert_series_equal(df.loc[df.new_col=='new','time'],v)

    def test_loc_setitem_dups(self):

        # GH 6541
        df_orig = DataFrame({'me' : list('rttti'),
                             'foo': list('aaade'),
                             'bar': np.arange(5,dtype='float64')*1.34+2,
                             'bar2': np.arange(5,dtype='float64')*-.34+2}).set_index('me')

        indexer = tuple(['r',['bar','bar2']])
        df = df_orig.copy()
        df.loc[indexer]*=2.0
        assert_series_equal(df.loc[indexer],2.0*df_orig.loc[indexer])

        indexer = tuple(['r','bar'])
        df = df_orig.copy()
        df.loc[indexer]*=2.0
        self.assertEqual(df.loc[indexer],2.0*df_orig.loc[indexer])

        indexer = tuple(['t',['bar','bar2']])
        df = df_orig.copy()
        df.loc[indexer]*=2.0
        assert_frame_equal(df.loc[indexer],2.0*df_orig.loc[indexer])

    def test_iloc_setitem_dups(self):

        # GH 6766
        # iloc with a mask aligning from another iloc
        df1 = DataFrame([{'A':None, 'B':1},{'A':2, 'B':2}])
        df2 = DataFrame([{'A':3, 'B':3},{'A':4, 'B':4}])
        df = concat([df1, df2], axis=1)

        expected = df.fillna(3)
        expected['A'] = expected['A'].astype('float64')
        inds = np.isnan(df.iloc[:, 0])
        mask = inds[inds].index
        df.iloc[mask,0] = df.iloc[mask,2]
        assert_frame_equal(df, expected)

        # del a dup column across blocks
        expected = DataFrame({ 0 : [1,2], 1 : [3,4] })
        expected.columns=['B','B']
        del df['A']
        assert_frame_equal(df, expected)

        # assign back to self
        df.iloc[[0,1],[0,1]] = df.iloc[[0,1],[0,1]]
        assert_frame_equal(df, expected)

        # reversed x 2
        df.iloc[[1,0],[0,1]] = df.iloc[[1,0],[0,1]].reset_index(drop=True)
        df.iloc[[1,0],[0,1]] = df.iloc[[1,0],[0,1]].reset_index(drop=True)
        assert_frame_equal(df, expected)

    def test_chained_getitem_with_lists(self):

        # GH6394
        # Regression in chained getitem indexing with embedded list-like from 0.12
        def check(result, expected):
            tm.assert_numpy_array_equal(result,expected)
            tm.assertIsInstance(result, np.ndarray)


        df = DataFrame({'A': 5*[np.zeros(3)], 'B':5*[np.ones(3)]})
        expected = df['A'].iloc[2]
        result = df.loc[2,'A']
        check(result, expected)
        result2 = df.iloc[2]['A']
        check(result2, expected)
        result3 = df['A'].loc[2]
        check(result3, expected)
        result4 = df['A'].iloc[2]
        check(result4, expected)

    def test_loc_getitem_int(self):

        # int label
        self.check_result('int label', 'loc', 2, 'ix', 2, typs = ['ints'], axes = 0)
        self.check_result('int label', 'loc', 3, 'ix', 3, typs = ['ints'], axes = 1)
        self.check_result('int label', 'loc', 4, 'ix', 4, typs = ['ints'], axes = 2)
        self.check_result('int label', 'loc', 2, 'ix', 2, typs = ['label'], fails = KeyError)

    def test_loc_getitem_label(self):

        # label
        self.check_result('label', 'loc', 'c',    'ix', 'c',    typs = ['labels'], axes=0)
        self.check_result('label', 'loc', 'null', 'ix', 'null', typs = ['mixed'] , axes=0)
        self.check_result('label', 'loc', 8,      'ix', 8,      typs = ['mixed'] , axes=0)
        self.check_result('label', 'loc', Timestamp('20130102'), 'ix', 1, typs = ['ts'], axes=0)
        self.check_result('label', 'loc', 'c', 'ix', 'c', typs = ['empty'], fails = KeyError)

    def test_loc_getitem_label_out_of_range(self):

        # out of range label
        self.check_result('label range', 'loc', 'f', 'ix', 'f', typs = ['ints','labels','mixed','ts'], fails=KeyError)
        self.check_result('label range', 'loc', 'f', 'ix', 'f', typs = ['floats'], fails=TypeError)
        self.check_result('label range', 'loc', 20, 'ix', 20, typs = ['ints','labels','mixed'], fails=KeyError)
        self.check_result('label range', 'loc', 20, 'ix', 20, typs = ['ts'], axes=0, fails=TypeError)
        self.check_result('label range', 'loc', 20, 'ix', 20, typs = ['floats'], axes=0, fails=TypeError)

    def test_loc_getitem_label_list(self):

        # list of labels
        self.check_result('list lbl', 'loc', [0,2,4], 'ix', [0,2,4], typs = ['ints'], axes=0)
        self.check_result('list lbl', 'loc', [3,6,9], 'ix', [3,6,9], typs = ['ints'], axes=1)
        self.check_result('list lbl', 'loc', [4,8,12], 'ix', [4,8,12], typs = ['ints'], axes=2)
        self.check_result('list lbl', 'loc', ['a','b','d'], 'ix', ['a','b','d'], typs = ['labels'], axes=0)
        self.check_result('list lbl', 'loc', ['A','B','C'], 'ix', ['A','B','C'], typs = ['labels'], axes=1)
        self.check_result('list lbl', 'loc', ['Z','Y','W'], 'ix', ['Z','Y','W'], typs = ['labels'], axes=2)
        self.check_result('list lbl', 'loc', [2,8,'null'], 'ix', [2,8,'null'], typs = ['mixed'], axes=0)
        self.check_result('list lbl', 'loc', [Timestamp('20130102'),Timestamp('20130103')], 'ix',
                          [Timestamp('20130102'),Timestamp('20130103')], typs = ['ts'], axes=0)

        self.check_result('list lbl', 'loc', [0,1,2], 'indexer', [0,1,2], typs = ['empty'], fails = KeyError)
        self.check_result('list lbl', 'loc', [0,2,3], 'ix', [0,2,3], typs = ['ints'], axes=0, fails = KeyError)
        self.check_result('list lbl', 'loc', [3,6,7], 'ix', [3,6,7], typs = ['ints'], axes=1, fails = KeyError)
        self.check_result('list lbl', 'loc', [4,8,10], 'ix', [4,8,10], typs = ['ints'], axes=2, fails = KeyError)

        # fails
        self.check_result('list lbl', 'loc', [20,30,40], 'ix', [20,30,40], typs = ['ints'], axes=1, fails = KeyError)
        self.check_result('list lbl', 'loc', [20,30,40], 'ix', [20,30,40], typs = ['ints'], axes=2, fails = KeyError)

        # array like
        self.check_result('array like', 'loc', Series(index=[0,2,4]).index, 'ix', [0,2,4], typs = ['ints'], axes=0)
        self.check_result('array like', 'loc', Series(index=[3,6,9]).index, 'ix', [3,6,9], typs = ['ints'], axes=1)
        self.check_result('array like', 'loc', Series(index=[4,8,12]).index, 'ix', [4,8,12], typs = ['ints'], axes=2)

    def test_loc_getitem_bool(self):

        # boolean indexers
        b = [True,False,True,False]
        self.check_result('bool', 'loc', b, 'ix', b, typs = ['ints','labels','mixed','ts','floats'])
        self.check_result('bool', 'loc', b, 'ix', b, typs = ['empty'], fails = KeyError)

    def test_loc_getitem_int_slice(self):

        # ok
        self.check_result('int slice2', 'loc', slice(2,4), 'ix', [2,4], typs = ['ints'], axes = 0)
        self.check_result('int slice2', 'loc', slice(3,6), 'ix', [3,6], typs = ['ints'], axes = 1)
        self.check_result('int slice2', 'loc', slice(4,8), 'ix', [4,8], typs = ['ints'], axes = 2)

        # GH 3053
        # loc should treat integer slices like label slices
        from itertools import product

        index = MultiIndex.from_tuples([t for t in product([6,7,8], ['a', 'b'])])
        df = DataFrame(np.random.randn(6, 6), index, index)
        result = df.loc[6:8,:]
        expected = df.ix[6:8,:]
        assert_frame_equal(result,expected)

        index = MultiIndex.from_tuples([t for t in product([10, 20, 30], ['a', 'b'])])
        df = DataFrame(np.random.randn(6, 6), index, index)
        result = df.loc[20:30,:]
        expected = df.ix[20:30,:]
        assert_frame_equal(result,expected)

        # doc examples
        result = df.loc[10,:]
        expected = df.ix[10,:]
        assert_frame_equal(result,expected)

        result = df.loc[:,10]
        #expected = df.ix[:,10] (this fails)
        expected = df[10]
        assert_frame_equal(result,expected)

    def test_loc_to_fail(self):

        # GH3449
        df = DataFrame(np.random.random((3, 3)),
                       index=['a', 'b', 'c'],
                       columns=['e', 'f', 'g'])

        # raise a KeyError?
        self.assertRaises(KeyError, df.loc.__getitem__, tuple([[1, 2], [1, 2]]))

        # GH  7496
        # loc should not fallback

        s = Series()
        s.loc[1] = 1
        s.loc['a'] = 2

        self.assertRaises(KeyError, lambda : s.loc[-1])
        self.assertRaises(KeyError, lambda : s.loc[[-1, -2]])

        self.assertRaises(KeyError, lambda : s.loc[['4']])

        s.loc[-1] = 3
        result = s.loc[[-1,-2]]
        expected = Series([3,np.nan],index=[-1,-2])
        assert_series_equal(result, expected)

        s['a'] = 2
        self.assertRaises(KeyError, lambda : s.loc[[-2]])

        del s['a']
        def f():
            s.loc[[-2]] = 0
        self.assertRaises(KeyError, f)

        # inconsistency between .loc[values] and .loc[values,:]
        # GH 7999
        df = DataFrame([['a'],['b']],index=[1,2],columns=['value'])

        def f():
            df.loc[[3],:]
        self.assertRaises(KeyError, f)

        def f():
            df.loc[[3]]
        self.assertRaises(KeyError, f)

        # at should not fallback
        # GH 7814
        s = Series([1,2,3], index=list('abc'))
        result = s.at['a']
        self.assertEqual(result, 1)
        self.assertRaises(ValueError, lambda : s.at[0])

        df = DataFrame({'A' : [1,2,3]},index=list('abc'))
        result = df.at['a','A']
        self.assertEqual(result, 1)
        self.assertRaises(ValueError, lambda : df.at['a',0])

        s = Series([1,2,3], index=[3,2,1])
        result = s.at[1]
        self.assertEqual(result, 3)
        self.assertRaises(ValueError, lambda : s.at['a'])

        df = DataFrame({0 : [1,2,3]},index=[3,2,1])
        result = df.at[1,0]
        self.assertEqual(result, 3)
        self.assertRaises(ValueError, lambda : df.at['a',0])

    def test_loc_getitem_label_slice(self):

        # label slices (with ints)
        self.check_result('lab slice', 'loc', slice(1,3), 'ix', slice(1,3), typs = ['labels','mixed','empty','ts','floats'], fails=TypeError)

        # real label slices
        self.check_result('lab slice', 'loc', slice('a','c'), 'ix', slice('a','c'), typs = ['labels'], axes=0)
        self.check_result('lab slice', 'loc', slice('A','C'), 'ix', slice('A','C'), typs = ['labels'], axes=1)
        self.check_result('lab slice', 'loc', slice('W','Z'), 'ix', slice('W','Z'), typs = ['labels'], axes=2)

        self.check_result('ts  slice', 'loc', slice('20130102','20130104'), 'ix', slice('20130102','20130104'), typs = ['ts'], axes=0)
        self.check_result('ts  slice', 'loc', slice('20130102','20130104'), 'ix', slice('20130102','20130104'), typs = ['ts'], axes=1, fails=TypeError)
        self.check_result('ts  slice', 'loc', slice('20130102','20130104'), 'ix', slice('20130102','20130104'), typs = ['ts'], axes=2, fails=TypeError)

        self.check_result('mixed slice', 'loc', slice(2,8), 'ix', slice(2,8), typs = ['mixed'], axes=0, fails=TypeError)
        self.check_result('mixed slice', 'loc', slice(2,8), 'ix', slice(2,8), typs = ['mixed'], axes=1, fails=KeyError)
        self.check_result('mixed slice', 'loc', slice(2,8), 'ix', slice(2,8), typs = ['mixed'], axes=2, fails=KeyError)

        self.check_result('mixed slice', 'loc', slice(2,4,2), 'ix', slice(2,4,2), typs = ['mixed'], axes=0, fails=TypeError)

    def test_loc_general(self):

        df = DataFrame(np.random.rand(4,4),columns=['A','B','C','D'], index=['A','B','C','D'])

        # want this to work
        result = df.loc[:,"A":"B"].iloc[0:2,:]
        self.assertTrue((result.columns == ['A','B']).all() == True)
        self.assertTrue((result.index == ['A','B']).all() == True)

        # mixed type
        result = DataFrame({ 'a' : [Timestamp('20130101')], 'b' : [1] }).iloc[0]
        expected = Series([ Timestamp('20130101'), 1], index=['a','b'], name=0)
        assert_series_equal(result, expected)
        self.assertEqual(result.dtype, object)

    def test_loc_setitem_consistency(self):

        # GH 6149
        # coerce similary for setitem and loc when rows have a null-slice
        expected = DataFrame({ 'date': Series(0,index=range(5),dtype=np.int64),
                               'val' : Series(range(5),dtype=np.int64) })

        df = DataFrame({ 'date': date_range('2000-01-01','2000-01-5'),
                         'val' : Series(range(5),dtype=np.int64) })
        df.loc[:,'date'] = 0
        assert_frame_equal(df,expected)

        df = DataFrame({ 'date': date_range('2000-01-01','2000-01-5'),
                         'val' : Series(range(5),dtype=np.int64) })
        df.loc[:,'date'] = np.array(0,dtype=np.int64)
        assert_frame_equal(df,expected)

        df = DataFrame({ 'date': date_range('2000-01-01','2000-01-5'),
                         'val' : Series(range(5),dtype=np.int64) })
        df.loc[:,'date'] = np.array([0,0,0,0,0],dtype=np.int64)
        assert_frame_equal(df,expected)

        expected = DataFrame({ 'date': Series('foo',index=range(5)),
                               'val' : Series(range(5),dtype=np.int64) })
        df = DataFrame({ 'date': date_range('2000-01-01','2000-01-5'),
                         'val' : Series(range(5),dtype=np.int64) })
        df.loc[:,'date'] = 'foo'
        assert_frame_equal(df,expected)

        expected = DataFrame({ 'date': Series(1.0,index=range(5)),
                               'val' : Series(range(5),dtype=np.int64) })
        df = DataFrame({ 'date': date_range('2000-01-01','2000-01-5'),
                         'val' : Series(range(5),dtype=np.int64) })
        df.loc[:,'date'] = 1.0
        assert_frame_equal(df,expected)

        # empty (essentially noops)
        expected = DataFrame(columns=['x', 'y'])
        expected['x'] = expected['x'].astype(np.int64)
        df = DataFrame(columns=['x', 'y'])
        df.loc[:, 'x'] = 1
        assert_frame_equal(df,expected)

        df = DataFrame(columns=['x', 'y'])
        df['x'] = 1
        assert_frame_equal(df,expected)

        # .loc[:,column] setting with slice == len of the column
        # GH10408
        data = """Level_0,,,Respondent,Respondent,Respondent,OtherCat,OtherCat
Level_1,,,Something,StartDate,EndDate,Yes/No,SomethingElse
Region,Site,RespondentID,,,,,
Region_1,Site_1,3987227376,A,5/25/2015 10:59,5/25/2015 11:22,Yes,
Region_1,Site_1,3980680971,A,5/21/2015 9:40,5/21/2015 9:52,Yes,Yes
Region_1,Site_2,3977723249,A,5/20/2015 8:27,5/20/2015 8:41,Yes,
Region_1,Site_2,3977723089,A,5/20/2015 8:33,5/20/2015 9:09,Yes,No"""

        df = pd.read_csv(StringIO(data),header=[0,1], index_col=[0,1,2])
        df.loc[:,('Respondent','StartDate')] = pd.to_datetime(df.loc[:,('Respondent','StartDate')])
        df.loc[:,('Respondent','EndDate')] = pd.to_datetime(df.loc[:,('Respondent','EndDate')])
        df.loc[:,('Respondent','Duration')] = df.loc[:,('Respondent','EndDate')] - df.loc[:,('Respondent','StartDate')]

        df.loc[:,('Respondent','Duration')] = df.loc[:,('Respondent','Duration')].astype('timedelta64[s]')
        expected = Series([1380,720,840,2160.],index=df.index,name=('Respondent','Duration'))
        assert_series_equal(df[('Respondent','Duration')],expected)

    def test_loc_setitem_frame(self):
        df = self.frame_labels

        result = df.iloc[0,0]

        df.loc['a','A'] = 1
        result = df.loc['a','A']
        self.assertEqual(result, 1)

        result = df.iloc[0,0]
        self.assertEqual(result, 1)

        df.loc[:,'B':'D'] = 0
        expected = df.loc[:,'B':'D']
        result = df.ix[:,1:]
        assert_frame_equal(result, expected)

        # GH 6254
        # setting issue
        df = DataFrame(index=[3, 5, 4], columns=['A'])
        df.loc[[4, 3, 5], 'A'] = np.array([1, 2, 3],dtype='int64')
        expected = DataFrame(dict(A = Series([1,2,3],index=[4, 3, 5]))).reindex(index=[3,5,4])
        assert_frame_equal(df, expected)

        # GH 6252
        # setting with an empty frame
        keys1 = ['@' + str(i) for i in range(5)]
        val1 = np.arange(5,dtype='int64')

        keys2 = ['@' + str(i) for i in range(4)]
        val2 = np.arange(4,dtype='int64')

        index = list(set(keys1).union(keys2))
        df = DataFrame(index = index)
        df['A'] = nan
        df.loc[keys1, 'A'] = val1

        df['B'] = nan
        df.loc[keys2, 'B'] = val2

        expected = DataFrame(dict(A = Series(val1,index=keys1), B = Series(val2,index=keys2))).reindex(index=index)
        assert_frame_equal(df, expected)

        # GH 8669
        # invalid coercion of nan -> int
        df = DataFrame({'A' : [1,2,3], 'B' : np.nan })
        df.loc[df.B > df.A, 'B'] = df.A
        expected = DataFrame({'A' : [1,2,3], 'B' : np.nan})
        assert_frame_equal(df, expected)

        # GH 6546
        # setting with mixed labels
        df = DataFrame({1:[1,2],2:[3,4],'a':['a','b']})

        result = df.loc[0, [1,2]]
        expected = Series([1,3],index=[1,2],dtype=object, name=0)
        assert_series_equal(result, expected)

        expected = DataFrame({1:[5,2],2:[6,4],'a':['a','b']})
        df.loc[0, [1,2]] = [5,6]
        assert_frame_equal(df, expected)

    def test_loc_setitem_frame_multiples(self):
        # multiple setting
        df = DataFrame({ 'A' : ['foo','bar','baz'],
                         'B' : Series(range(3),dtype=np.int64) })
        rhs = df.loc[1:2]
        rhs.index = df.index[0:2]
        df.loc[0:1] = rhs
        expected = DataFrame({ 'A' : ['bar','baz','baz'],
                               'B' : Series([1,2,2],dtype=np.int64) })
        assert_frame_equal(df, expected)


        # multiple setting with frame on rhs (with M8)
        df = DataFrame({ 'date' : date_range('2000-01-01','2000-01-5'),
                         'val'  : Series(range(5),dtype=np.int64) })
        expected = DataFrame({ 'date' : [Timestamp('20000101'),Timestamp('20000102'),Timestamp('20000101'),
                                         Timestamp('20000102'),Timestamp('20000103')],
                               'val'  : Series([0,1,0,1,2],dtype=np.int64) })
        rhs = df.loc[0:2]
        rhs.index = df.index[2:5]
        df.loc[2:4] = rhs
        assert_frame_equal(df, expected)

    def test_iloc_getitem_frame(self):
        df = DataFrame(np.random.randn(10, 4), index=lrange(0, 20, 2), columns=lrange(0,8,2))

        result = df.iloc[2]
        exp = df.ix[4]
        assert_series_equal(result, exp)

        result = df.iloc[2,2]
        exp = df.ix[4,4]
        self.assertEqual(result, exp)

        # slice
        result = df.iloc[4:8]
        expected = df.ix[8:14]
        assert_frame_equal(result, expected)

        result = df.iloc[:,2:3]
        expected = df.ix[:,4:5]
        assert_frame_equal(result, expected)

        # list of integers
        result = df.iloc[[0,1,3]]
        expected = df.ix[[0,2,6]]
        assert_frame_equal(result, expected)

        result = df.iloc[[0,1,3],[0,1]]
        expected = df.ix[[0,2,6],[0,2]]
        assert_frame_equal(result, expected)

        # neg indicies
        result = df.iloc[[-1,1,3],[-1,1]]
        expected = df.ix[[18,2,6],[6,2]]
        assert_frame_equal(result, expected)

        # dups indicies
        result = df.iloc[[-1,-1,1,3],[-1,1]]
        expected = df.ix[[18,18,2,6],[6,2]]
        assert_frame_equal(result, expected)

        # with index-like
        s = Series(index=lrange(1,5))
        result = df.iloc[s.index]
        expected = df.ix[[2,4,6,8]]
        assert_frame_equal(result, expected)

        # try with labelled frame
        df = DataFrame(np.random.randn(10, 4), index=list('abcdefghij'), columns=list('ABCD'))

        result = df.iloc[1,1]
        exp = df.ix['b','B']
        self.assertEqual(result, exp)

        result = df.iloc[:,2:3]
        expected = df.ix[:,['C']]
        assert_frame_equal(result, expected)

        # negative indexing
        result = df.iloc[-1,-1]
        exp = df.ix['j','D']
        self.assertEqual(result, exp)

        # out-of-bounds exception
        self.assertRaises(IndexError, df.iloc.__getitem__, tuple([10,5]))

        # trying to use a label
        self.assertRaises(ValueError, df.iloc.__getitem__, tuple(['j','D']))

    def test_iloc_getitem_panel(self):

        # GH 7189
        p = Panel(np.arange(4*3*2).reshape(4,3,2),
                  items=['A','B','C','D'],
                  major_axis=['a','b','c'],
                  minor_axis=['one','two'])

        result = p.iloc[1]
        expected = p.loc['B']
        assert_frame_equal(result, expected)

        result = p.iloc[1,1]
        expected = p.loc['B','b']
        assert_series_equal(result, expected)

        result = p.iloc[1,1,1]
        expected = p.loc['B','b','two']
        self.assertEqual(result,expected)

        # slice
        result = p.iloc[1:3]
        expected = p.loc[['B','C']]
        assert_panel_equal(result, expected)

        result = p.iloc[:,0:2]
        expected = p.loc[:,['a','b']]
        assert_panel_equal(result, expected)

        # list of integers
        result = p.iloc[[0,2]]
        expected = p.loc[['A','C']]
        assert_panel_equal(result, expected)

        # neg indicies
        result = p.iloc[[-1,1],[-1,1]]
        expected = p.loc[['D','B'],['c','b']]
        assert_panel_equal(result, expected)

        # dups indicies
        result = p.iloc[[-1,-1,1],[-1,1]]
        expected = p.loc[['D','D','B'],['c','b']]
        assert_panel_equal(result, expected)

        # combined
        result = p.iloc[0,[True,True],[0,1]]
        expected = p.loc['A',['a','b'],['one','two']]
        assert_frame_equal(result, expected)

        # out-of-bounds exception
        self.assertRaises(IndexError, p.iloc.__getitem__, tuple([10,5]))
        def f():
            p.iloc[0,[True,True],[0,1,2]]
        self.assertRaises(IndexError, f)

        # trying to use a label
        self.assertRaises(ValueError, p.iloc.__getitem__, tuple(['j','D']))

        # GH
        p = Panel(np.random.rand(4,3,2), items=['A','B','C','D'], major_axis=['U','V','W'], minor_axis=['X','Y'])
        expected = p['A']

        result = p.iloc[0,:,:]
        assert_frame_equal(result, expected)

        result = p.iloc[0,[True,True,True],:]
        assert_frame_equal(result, expected)

        result = p.iloc[0,[True,True,True],[0,1]]
        assert_frame_equal(result, expected)

        def f():
            p.iloc[0,[True,True,True],[0,1,2]]
        self.assertRaises(IndexError, f)

        def f():
            p.iloc[0,[True,True,True],[2]]
        self.assertRaises(IndexError, f)

        # GH 7199
        # Panel with multi-index
        multi_index = pd.MultiIndex.from_tuples([('ONE', 'one'),
                                                 ('TWO', 'two'),
                                                 ('THREE', 'three')],
                                                names=['UPPER', 'lower'])

        simple_index = [x[0] for x in multi_index]
        wd1 = Panel(items=['First', 'Second'],
                    major_axis=['a', 'b', 'c', 'd'],
                    minor_axis=multi_index)

        wd2 = Panel(items=['First', 'Second'],
                    major_axis=['a', 'b', 'c', 'd'],
                    minor_axis=simple_index)

        expected1 = wd1['First'].iloc[[True, True, True, False], [0, 2]]
        result1 = wd1.iloc[0, [True, True, True, False], [0, 2]]  # WRONG
        assert_frame_equal(result1,expected1)

        expected2 = wd2['First'].iloc[[True, True, True, False], [0, 2]]
        result2 = wd2.iloc[0, [True, True, True, False], [0, 2]]
        assert_frame_equal(result2,expected2)

        expected1 = DataFrame(index=['a'],columns=multi_index,dtype='float64')
        result1 = wd1.iloc[0,[0],[0,1,2]]
        assert_frame_equal(result1,expected1)

        expected2 = DataFrame(index=['a'],columns=simple_index,dtype='float64')
        result2 = wd2.iloc[0,[0],[0,1,2]]
        assert_frame_equal(result2,expected2)

        # GH 7516
        mi = MultiIndex.from_tuples([(0,'x'), (1,'y'), (2,'z')])
        p = Panel(np.arange(3*3*3,dtype='int64').reshape(3,3,3), items=['a','b','c'], major_axis=mi, minor_axis=['u','v','w'])
        result = p.iloc[:, 1, 0]
        expected = Series([3,12,21],index=['a','b','c'], name='u')
        assert_series_equal(result,expected)

        result = p.loc[:, (1,'y'), 'u']
        assert_series_equal(result,expected)

    def test_iloc_getitem_doc_issue(self):

        # multi axis slicing issue with single block
        # surfaced in GH 6059

        arr = np.random.randn(6,4)
        index = date_range('20130101',periods=6)
        columns = list('ABCD')
        df = DataFrame(arr,index=index,columns=columns)

        # defines ref_locs
        df.describe()

        result = df.iloc[3:5,0:2]
        str(result)
        result.dtypes

        expected = DataFrame(arr[3:5,0:2],index=index[3:5],columns=columns[0:2])
        assert_frame_equal(result,expected)

        # for dups
        df.columns = list('aaaa')
        result = df.iloc[3:5,0:2]
        str(result)
        result.dtypes

        expected = DataFrame(arr[3:5,0:2],index=index[3:5],columns=list('aa'))
        assert_frame_equal(result,expected)

        # related
        arr = np.random.randn(6,4)
        index = list(range(0,12,2))
        columns = list(range(0,8,2))
        df = DataFrame(arr,index=index,columns=columns)

        df._data.blocks[0].mgr_locs
        result = df.iloc[1:5,2:4]
        str(result)
        result.dtypes
        expected = DataFrame(arr[1:5,2:4],index=index[1:5],columns=columns[2:4])
        assert_frame_equal(result,expected)

    def test_setitem_ndarray_1d(self):
        # GH5508

        # len of indexer vs length of the 1d ndarray
        df = DataFrame(index=Index(lrange(1,11)))
        df['foo'] = np.zeros(10, dtype=np.float64)
        df['bar'] = np.zeros(10, dtype=np.complex)

        # invalid
        def f():
            df.ix[2:5, 'bar'] = np.array([2.33j, 1.23+0.1j, 2.2])
        self.assertRaises(ValueError, f)

        # valid
        df.ix[2:5, 'bar'] = np.array([2.33j, 1.23+0.1j, 2.2, 1.0])

        result = df.ix[2:5, 'bar']
        expected = Series([2.33j, 1.23+0.1j, 2.2, 1.0], index=[2,3,4,5], name='bar')
        assert_series_equal(result, expected)

        # dtype getting changed?
        df = DataFrame(index=Index(lrange(1,11)))
        df['foo'] = np.zeros(10, dtype=np.float64)
        df['bar'] = np.zeros(10, dtype=np.complex)

        def f():
            df[2:5] = np.arange(1,4)*1j
        self.assertRaises(ValueError, f)

    def test_iloc_setitem_series(self):
        df = DataFrame(np.random.randn(10, 4), index=list('abcdefghij'), columns=list('ABCD'))

        df.iloc[1,1] = 1
        result = df.iloc[1,1]
        self.assertEqual(result, 1)

        df.iloc[:,2:3] = 0
        expected = df.iloc[:,2:3]
        result = df.iloc[:,2:3]
        assert_frame_equal(result, expected)

        s = Series(np.random.randn(10), index=lrange(0,20,2))

        s.iloc[1] = 1
        result = s.iloc[1]
        self.assertEqual(result, 1)

        s.iloc[:4] = 0
        expected = s.iloc[:4]
        result = s.iloc[:4]
        assert_series_equal(result, expected)

        s= Series([-1]*6)
        s.iloc[0::2]= [0,2,4]
        s.iloc[1::2]= [1,3,5]
        result  = s
        expected= Series([0,1,2,3,4,5])
        assert_series_equal(result, expected)

    def test_iloc_setitem_list_of_lists(self):

        # GH 7551
        # list-of-list is set incorrectly in mixed vs. single dtyped frames
        df = DataFrame(dict(A = np.arange(5,dtype='int64'), B = np.arange(5,10,dtype='int64')))
        df.iloc[2:4] = [[10,11],[12,13]]
        expected = DataFrame(dict(A = [0,1,10,12,4], B = [5,6,11,13,9]))
        assert_frame_equal(df, expected)

        df = DataFrame(dict(A = list('abcde'), B = np.arange(5,10,dtype='int64')))
        df.iloc[2:4] = [['x',11],['y',13]]
        expected = DataFrame(dict(A = ['a','b','x','y','e'], B = [5,6,11,13,9]))
        assert_frame_equal(df, expected)

    def test_iloc_getitem_multiindex(self):
        mi_labels = DataFrame(np.random.randn(4, 3), columns=[['i', 'i', 'j'],
                                                              ['A', 'A', 'B']],
                              index=[['i', 'i', 'j', 'k'], ['X', 'X', 'Y','Y']])

        mi_int    = DataFrame(np.random.randn(3, 3),
                              columns=[[2,2,4],[6,8,10]],
                              index=[[4,4,8],[8,10,12]])


        # the first row
        rs = mi_int.iloc[0]
        xp = mi_int.ix[4].ix[8]
        assert_series_equal(rs, xp, check_names=False)
        self.assertEqual(rs.name, (4, 8))
        self.assertEqual(xp.name, 8)

        # 2nd (last) columns
        rs = mi_int.iloc[:,2]
        xp = mi_int.ix[:,2]
        assert_series_equal(rs, xp)

        # corner column
        rs = mi_int.iloc[2,2]
        xp = mi_int.ix[:,2].ix[2]
        self.assertEqual(rs, xp)

        # this is basically regular indexing
        rs = mi_labels.iloc[2,2]
        xp = mi_labels.ix['j'].ix[:,'j'].ix[0,0]
        self.assertEqual(rs, xp)

    def test_loc_multiindex(self):

        mi_labels = DataFrame(np.random.randn(3, 3), columns=[['i', 'i', 'j'],
                                                              ['A', 'A', 'B']],
                              index=[['i', 'i', 'j'], ['X', 'X', 'Y']])

        mi_int    = DataFrame(np.random.randn(3, 3),
                              columns=[[2,2,4],[6,8,10]],
                              index=[[4,4,8],[8,10,12]])

        # the first row
        rs = mi_labels.loc['i']
        xp = mi_labels.ix['i']
        assert_frame_equal(rs, xp)

        # 2nd (last) columns
        rs = mi_labels.loc[:,'j']
        xp = mi_labels.ix[:,'j']
        assert_frame_equal(rs, xp)

        # corner column
        rs = mi_labels.loc['j'].loc[:,'j']
        xp = mi_labels.ix['j'].ix[:,'j']
        assert_frame_equal(rs,xp)

        # with a tuple
        rs = mi_labels.loc[('i','X')]
        xp = mi_labels.ix[('i','X')]
        assert_frame_equal(rs,xp)

        rs = mi_int.loc[4]
        xp = mi_int.ix[4]
        assert_frame_equal(rs,xp)

        # GH6788
        # multi-index indexer is None (meaning take all)
        attributes = ['Attribute' + str(i) for i in range(1)]
        attribute_values = ['Value' + str(i) for i in range(5)]

        index = MultiIndex.from_product([attributes,attribute_values])
        df = 0.1 * np.random.randn(10, 1 * 5) + 0.5
        df = DataFrame(df, columns=index)
        result = df[attributes]
        assert_frame_equal(result, df)

        # GH 7349
        # loc with a multi-index seems to be doing fallback
        df = DataFrame(np.arange(12).reshape(-1,1),index=pd.MultiIndex.from_product([[1,2,3,4],[1,2,3]]))

        expected = df.loc[([1,2],),:]
        result = df.loc[[1,2]]
        assert_frame_equal(result, expected)

        # GH 7399
        # incomplete indexers
        s = pd.Series(np.arange(15,dtype='int64'),MultiIndex.from_product([range(5), ['a', 'b', 'c']]))
        expected = s.loc[:, 'a':'c']

        result = s.loc[0:4, 'a':'c']
        assert_series_equal(result, expected)
        assert_series_equal(result, expected)

        result = s.loc[:4, 'a':'c']
        assert_series_equal(result, expected)
        assert_series_equal(result, expected)

        result = s.loc[0:, 'a':'c']
        assert_series_equal(result, expected)
        assert_series_equal(result, expected)

        # GH 7400
        # multiindexer gettitem with list of indexers skips wrong element
        s = pd.Series(np.arange(15,dtype='int64'),MultiIndex.from_product([range(5), ['a', 'b', 'c']]))
        expected = s.iloc[[6,7,8,12,13,14]]
        result = s.loc[2:4:2, 'a':'c']
        assert_series_equal(result, expected)

    def test_multiindex_perf_warn(self):

        if sys.version_info < (2, 7):
            raise nose.SkipTest('python version < 2.7')

        df = DataFrame({'jim':[0, 0, 1, 1],
                        'joe':['x', 'x', 'z', 'y'],
                        'jolie':np.random.rand(4)}).set_index(['jim', 'joe'])

        with tm.assert_produces_warning(PerformanceWarning, clear=[pd.core.index]):
            _ = df.loc[(1, 'z')]

        df = df.iloc[[2,1,3,0]]
        with tm.assert_produces_warning(PerformanceWarning):
            _ = df.loc[(0,)]

    @slow
    def test_multiindex_get_loc(self):  # GH7724, GH2646

        with warnings.catch_warnings(record=True):

            # test indexing into a multi-index before & past the lexsort depth
            from numpy.random import randint, choice, randn
            cols = ['jim', 'joe', 'jolie', 'joline', 'jolia']

            def validate(mi, df, key):
                mask = np.ones(len(df)).astype('bool')

                # test for all partials of this key
                for i, k in enumerate(key):
                    mask &= df.iloc[:, i] == k

                    if not mask.any():
                        self.assertNotIn(key[:i+1], mi.index)
                        continue

                    self.assertIn(key[:i+1], mi.index)
                    right = df[mask].copy()

                    if i + 1 != len(key):  # partial key
                        right.drop(cols[:i+1], axis=1, inplace=True)
                        right.set_index(cols[i+1:-1], inplace=True)
                        assert_frame_equal(mi.loc[key[:i+1]], right)

                    else:  # full key
                        right.set_index(cols[:-1], inplace=True)
                        if len(right) == 1:  # single hit
                            right = Series(right['jolia'].values,
                                           name=right.index[0], index=['jolia'])
                            assert_series_equal(mi.loc[key[:i+1]], right)
                        else:  # multi hit
                            assert_frame_equal(mi.loc[key[:i+1]], right)

            def loop(mi, df, keys):
                for key in keys:
                    validate(mi, df, key)

            n, m = 1000, 50

            vals = [randint(0, 10, n), choice(list('abcdefghij'), n),
                    choice(pd.date_range('20141009', periods=10).tolist(), n),
                    choice(list('ZYXWVUTSRQ'), n), randn(n)]
            vals = list(map(tuple, zip(*vals)))

            # bunch of keys for testing
            keys = [randint(0, 11, m), choice(list('abcdefghijk'), m),
                    choice(pd.date_range('20141009', periods=11).tolist(), m),
                    choice(list('ZYXWVUTSRQP'), m)]
            keys = list(map(tuple, zip(*keys)))
            keys += list(map(lambda t: t[:-1], vals[::n//m]))

            # covers both unique index and non-unique index
            df = pd.DataFrame(vals, columns=cols)
            a, b = pd.concat([df, df]), df.drop_duplicates(subset=cols[:-1])

            for frame in a, b:
                for i in range(5):  # lexsort depth
                    df = frame.copy() if i == 0 else frame.sort_values(by=cols[:i])
                    mi = df.set_index(cols[:-1])
                    assert not mi.index.lexsort_depth < i
                    loop(mi, df, keys)

    def test_series_getitem_multiindex(self):

        # GH 6018
        # series regression getitem with a multi-index

        s = Series([1,2,3])
        s.index = MultiIndex.from_tuples([(0,0),(1,1), (2,1)])

        result = s[:,0]
        expected = Series([1],index=[0])
        assert_series_equal(result,expected)

        result = s.ix[:,1]
        expected = Series([2,3],index=[1,2])
        assert_series_equal(result,expected)

        # xs
        result = s.xs(0,level=0)
        expected = Series([1],index=[0])
        assert_series_equal(result,expected)

        result = s.xs(1,level=1)
        expected = Series([2,3],index=[1,2])
        assert_series_equal(result,expected)

        # GH6258
        s = Series([1,3,4,1,3,4],
                   index=MultiIndex.from_product([list('AB'),
                                                  list(date_range('20130903',periods=3))]))
        result = s.xs('20130903',level=1)
        expected = Series([1,1],index=list('AB'))
        assert_series_equal(result,expected)

        # GH5684
        idx = MultiIndex.from_tuples([('a', 'one'), ('a', 'two'),
                                      ('b', 'one'), ('b', 'two')])
        s = Series([1, 2, 3, 4], index=idx)
        s.index.set_names(['L1', 'L2'], inplace=True)
        result = s.xs('one', level='L2')
        expected = Series([1, 3], index=['a', 'b'])
        expected.index.set_names(['L1'], inplace=True)
        assert_series_equal(result, expected)

    def test_ix_general(self):

        # ix general issues

        # GH 2817
        data = {'amount': {0: 700, 1: 600, 2: 222, 3: 333, 4: 444},
                'col': {0: 3.5, 1: 3.5, 2: 4.0, 3: 4.0, 4: 4.0},
                'year': {0: 2012, 1: 2011, 2: 2012, 3: 2012, 4: 2012}}
        df = DataFrame(data).set_index(keys=['col', 'year'])
        key = 4.0, 2012

        # emits a PerformanceWarning, ok
        with self.assert_produces_warning(PerformanceWarning):
            tm.assert_frame_equal(df.ix[key], df.iloc[2:])

        # this is ok
        df.sortlevel(inplace=True)
        res = df.ix[key]
        # col has float dtype, result should be Float64Index
        index = MultiIndex.from_arrays([[4.] * 3, [2012] * 3],
                                       names=['col', 'year'])
        expected = DataFrame({'amount': [222, 333, 444]}, index=index)
        tm.assert_frame_equal(res, expected)

    def test_ix_weird_slicing(self):
        ## http://stackoverflow.com/q/17056560/1240268
        df = DataFrame({'one' : [1, 2, 3, np.nan, np.nan], 'two' : [1, 2, 3, 4, 5]})
        df.ix[df['one']>1, 'two'] = -df['two']

        expected = DataFrame({'one': {0: 1.0, 1: 2.0, 2: 3.0, 3: nan, 4: nan},
                              'two': {0: 1, 1: -2, 2: -3, 3: 4, 4: 5}})
        assert_frame_equal(df, expected)

    def test_xs_multiindex(self):

        # GH2903
        columns = MultiIndex.from_tuples([('a', 'foo'), ('a', 'bar'), ('b', 'hello'), ('b', 'world')], names=['lvl0', 'lvl1'])
        df = DataFrame(np.random.randn(4, 4), columns=columns)
        df.sortlevel(axis=1,inplace=True)
        result = df.xs('a', level='lvl0', axis=1)
        expected = df.iloc[:,0:2].loc[:,'a']
        assert_frame_equal(result,expected)

        result = df.xs('foo', level='lvl1', axis=1)
        expected = df.iloc[:, 1:2].copy()
        expected.columns = expected.columns.droplevel('lvl1')
        assert_frame_equal(result, expected)

    def test_per_axis_per_level_getitem(self):

        # GH6134
        # example test case
        ix = MultiIndex.from_product([_mklbl('A',5),_mklbl('B',7),_mklbl('C',4),_mklbl('D',2)])
        df = DataFrame(np.arange(len(ix.get_values())),index=ix)

        result = df.loc[(slice('A1','A3'),slice(None), ['C1','C3']),:]
        expected = df.loc[[ tuple([a,b,c,d]) for a,b,c,d in df.index.values if (
            a == 'A1' or a == 'A2' or a == 'A3') and (c == 'C1' or c == 'C3')]]
        assert_frame_equal(result, expected)

        expected = df.loc[[ tuple([a,b,c,d]) for a,b,c,d in df.index.values if (
            a == 'A1' or a == 'A2' or a == 'A3') and (c == 'C1' or c == 'C2' or c == 'C3')]]
        result = df.loc[(slice('A1','A3'),slice(None), slice('C1','C3')),:]
        assert_frame_equal(result, expected)

        # test multi-index slicing with per axis and per index controls
        index = MultiIndex.from_tuples([('A',1),('A',2),('A',3),('B',1)],
                                       names=['one','two'])
        columns = MultiIndex.from_tuples([('a','foo'),('a','bar'),('b','foo'),('b','bah')],
                                         names=['lvl0', 'lvl1'])

        df = DataFrame(np.arange(16,dtype='int64').reshape(4, 4), index=index, columns=columns)
        df = df.sortlevel(axis=0).sortlevel(axis=1)

        # identity
        result = df.loc[(slice(None),slice(None)),:]
        assert_frame_equal(result, df)
        result = df.loc[(slice(None),slice(None)),(slice(None),slice(None))]
        assert_frame_equal(result, df)
        result = df.loc[:,(slice(None),slice(None))]
        assert_frame_equal(result, df)

        # index
        result = df.loc[(slice(None),[1]),:]
        expected = df.iloc[[0,3]]
        assert_frame_equal(result, expected)

        result = df.loc[(slice(None),1),:]
        expected = df.iloc[[0,3]]
        assert_frame_equal(result, expected)

        # columns
        result = df.loc[:,(slice(None),['foo'])]
        expected = df.iloc[:,[1,3]]
        assert_frame_equal(result, expected)

        # both
        result = df.loc[(slice(None),1),(slice(None),['foo'])]
        expected = df.iloc[[0,3],[1,3]]
        assert_frame_equal(result, expected)

        result = df.loc['A','a']
        expected = DataFrame(dict(bar = [1,5,9], foo = [0,4,8]),
                             index=Index([1,2,3],name='two'),
                             columns=Index(['bar','foo'],name='lvl1'))
        assert_frame_equal(result, expected)

        result = df.loc[(slice(None),[1,2]),:]
        expected = df.iloc[[0,1,3]]
        assert_frame_equal(result, expected)

        # multi-level series
        s = Series(np.arange(len(ix.get_values())),index=ix)
        result = s.loc['A1':'A3', :, ['C1','C3']]
        expected = s.loc[[ tuple([a,b,c,d]) for a,b,c,d in s.index.values if (
            a == 'A1' or a == 'A2' or a == 'A3') and (c == 'C1' or c == 'C3')]]
        assert_series_equal(result, expected)

        # boolean indexers
        result = df.loc[(slice(None),df.loc[:,('a','bar')]>5),:]
        expected = df.iloc[[2,3]]
        assert_frame_equal(result, expected)

        def f():
            df.loc[(slice(None),np.array([True,False])),:]
        self.assertRaises(ValueError, f)

        # ambiguous cases
        # these can be multiply interpreted (e.g. in this case
        # as df.loc[slice(None),[1]] as well
        self.assertRaises(KeyError, lambda : df.loc[slice(None),[1]])

        result = df.loc[(slice(None),[1]),:]
        expected = df.iloc[[0,3]]
        assert_frame_equal(result, expected)

        # not lexsorted
        self.assertEqual(df.index.lexsort_depth,2)
        df = df.sortlevel(level=1,axis=0)
        self.assertEqual(df.index.lexsort_depth,0)
        with tm.assertRaisesRegexp(KeyError, 'MultiIndex Slicing requires the index to be fully lexsorted tuple len \(2\), lexsort depth \(0\)'):
            df.loc[(slice(None),df.loc[:,('a','bar')]>5),:]

    def test_multiindex_slicers_non_unique(self):

        # GH 7106
        # non-unique mi index support
        df = DataFrame(dict(A = ['foo','foo','foo','foo'],
                            B = ['a','a','a','a'],
                            C = [1,2,1,3],
                            D = [1,2,3,4])).set_index(['A','B','C']).sortlevel()
        self.assertFalse(df.index.is_unique)
        expected = DataFrame(dict(A = ['foo','foo'],
                                  B = ['a','a'],
                                  C = [1,1],
                                  D = [1,3])).set_index(['A','B','C']).sortlevel()
        result = df.loc[(slice(None),slice(None),1),:]
        assert_frame_equal(result, expected)

        # this is equivalent of an xs expression
        result = df.xs(1,level=2,drop_level=False)
        assert_frame_equal(result, expected)

        df = DataFrame(dict(A = ['foo','foo','foo','foo'],
                            B = ['a','a','a','a'],
                            C = [1,2,1,2],
                            D = [1,2,3,4])).set_index(['A','B','C']).sortlevel()
        self.assertFalse(df.index.is_unique)
        expected = DataFrame(dict(A = ['foo','foo'],
                                  B = ['a','a'],
                                  C = [1,1],
                                  D = [1,3])).set_index(['A','B','C']).sortlevel()
        result = df.loc[(slice(None),slice(None),1),:]
        self.assertFalse(result.index.is_unique)
        assert_frame_equal(result, expected)

    def test_multiindex_slicers_datetimelike(self):

        # GH 7429
        # buggy/inconsistent behavior when slicing with datetime-like
        import datetime
        dates = [datetime.datetime(2012,1,1,12,12,12) + datetime.timedelta(days=i) for i in range(6)]
        freq = [1,2]
        index = MultiIndex.from_product([dates,freq], names=['date','frequency'])

        df = DataFrame(np.arange(6*2*4,dtype='int64').reshape(-1,4),index=index,columns=list('ABCD'))

        # multi-axis slicing
        idx = pd.IndexSlice
        expected = df.iloc[[0,2,4],[0,1]]
        result = df.loc[(slice(Timestamp('2012-01-01 12:12:12'),Timestamp('2012-01-03 12:12:12')),slice(1,1)), slice('A','B')]
        assert_frame_equal(result,expected)

        result = df.loc[(idx[Timestamp('2012-01-01 12:12:12'):Timestamp('2012-01-03 12:12:12')],idx[1:1]), slice('A','B')]
        assert_frame_equal(result,expected)

        result = df.loc[(slice(Timestamp('2012-01-01 12:12:12'),Timestamp('2012-01-03 12:12:12')),1), slice('A','B')]
        assert_frame_equal(result,expected)

        # with strings
        result = df.loc[(slice('2012-01-01 12:12:12','2012-01-03 12:12:12'),slice(1,1)), slice('A','B')]
        assert_frame_equal(result,expected)

        result = df.loc[(idx['2012-01-01 12:12:12':'2012-01-03 12:12:12'],1), idx['A','B']]
        assert_frame_equal(result,expected)


    def test_multiindex_slicers_edges(self):

        # GH 8132
        # various edge cases
        df = DataFrame({'A': ['A0'] * 5 + ['A1']*5 + ['A2']*5,
                        'B': ['B0','B0','B1','B1','B2'] * 3,
                        'DATE': ["2013-06-11",
                                 "2013-07-02",
                                 "2013-07-09",
                                 "2013-07-30",
                                 "2013-08-06",
                                 "2013-06-11",
                                 "2013-07-02",
                                 "2013-07-09",
                                 "2013-07-30",
                                 "2013-08-06",
                                 "2013-09-03",
                                 "2013-10-01",
                                 "2013-07-09",
                                 "2013-08-06",
                                 "2013-09-03"],
                        'VALUES': [22, 35, 14,  9,  4, 40, 18, 4, 2, 5, 1, 2, 3,4, 2]})

        df['DATE'] = pd.to_datetime(df['DATE'])
        df1 = df.set_index(['A', 'B', 'DATE'])
        df1 = df1.sortlevel()
        df2 = df.set_index('DATE')

        # A1 - Get all values under "A0" and "A1"
        result = df1.loc[(slice('A1')),:]
        expected = df1.iloc[0:10]
        assert_frame_equal(result, expected)

        # A2 - Get all values from the start to "A2"
        result = df1.loc[(slice('A2')),:]
        expected = df1
        assert_frame_equal(result, expected)

        # A3 - Get all values under "B1" or "B2"
        result = df1.loc[(slice(None),slice('B1','B2')),:]
        expected = df1.iloc[[2,3,4,7,8,9,12,13,14]]
        assert_frame_equal(result, expected)

        # A4 - Get all values between 2013-07-02 and 2013-07-09
        result = df1.loc[(slice(None),slice(None),slice('20130702','20130709')),:]
        expected = df1.iloc[[1,2,6,7,12]]
        assert_frame_equal(result, expected)

        # B1 - Get all values in B0 that are also under A0, A1 and A2
        result = df1.loc[(slice('A2'),slice('B0')),:]
        expected = df1.iloc[[0,1,5,6,10,11]]
        assert_frame_equal(result, expected)

        # B2 - Get all values in B0, B1 and B2 (similar to what #2 is doing for the As)
        result = df1.loc[(slice(None),slice('B2')),:]
        expected = df1
        assert_frame_equal(result, expected)

        # B3 - Get all values from B1 to B2 and up to 2013-08-06
        result = df1.loc[(slice(None),slice('B1','B2'),slice('2013-08-06')),:]
        expected = df1.iloc[[2,3,4,7,8,9,12,13]]
        assert_frame_equal(result, expected)

        # B4 - Same as A4 but the start of the date slice is not a key.
        #      shows indexing on a partial selection slice
        result = df1.loc[(slice(None),slice(None),slice('20130701','20130709')),:]
        expected = df1.iloc[[1,2,6,7,12]]
        assert_frame_equal(result, expected)

    def test_per_axis_per_level_doc_examples(self):

        # test index maker
        idx = pd.IndexSlice

        # from indexing.rst / advanced
        index = MultiIndex.from_product([_mklbl('A',4),
                                         _mklbl('B',2),
                                         _mklbl('C',4),
                                         _mklbl('D',2)])
        columns = MultiIndex.from_tuples([('a','foo'),('a','bar'),
                                          ('b','foo'),('b','bah')],
                                         names=['lvl0', 'lvl1'])
        df = DataFrame(np.arange(len(index)*len(columns),dtype='int64').reshape((len(index),len(columns))),
                       index=index,
                       columns=columns)
        result = df.loc[(slice('A1','A3'),slice(None), ['C1','C3']),:]
        expected = df.loc[[ tuple([a,b,c,d]) for a,b,c,d in df.index.values if (
            a == 'A1' or a == 'A2' or a == 'A3') and (c == 'C1' or c == 'C3')]]
        assert_frame_equal(result, expected)
        result = df.loc[idx['A1':'A3',:,['C1','C3']],:]
        assert_frame_equal(result, expected)

        result = df.loc[(slice(None),slice(None), ['C1','C3']),:]
        expected = df.loc[[ tuple([a,b,c,d]) for a,b,c,d in df.index.values if (
            c == 'C1' or c == 'C3')]]
        assert_frame_equal(result, expected)
        result = df.loc[idx[:,:,['C1','C3']],:]
        assert_frame_equal(result, expected)

        # not sorted
        def f():
            df.loc['A1',(slice(None),'foo')]
        self.assertRaises(KeyError, f)
        df = df.sortlevel(axis=1)

        # slicing
        df.loc['A1',(slice(None),'foo')]
        df.loc[(slice(None),slice(None), ['C1','C3']),(slice(None),'foo')]

        # setitem
        df.loc(axis=0)[:,:,['C1','C3']] = -10

    def test_loc_arguments(self):

        index = MultiIndex.from_product([_mklbl('A',4),
                                         _mklbl('B',2),
                                         _mklbl('C',4),
                                         _mklbl('D',2)])
        columns = MultiIndex.from_tuples([('a','foo'),('a','bar'),
                                          ('b','foo'),('b','bah')],
                                         names=['lvl0', 'lvl1'])
        df = DataFrame(np.arange(len(index)*len(columns),dtype='int64').reshape((len(index),len(columns))),
                       index=index,
                       columns=columns).sortlevel().sortlevel(axis=1)


        # axis 0
        result = df.loc(axis=0)['A1':'A3',:,['C1','C3']]
        expected = df.loc[[ tuple([a,b,c,d]) for a,b,c,d in df.index.values if (
            a == 'A1' or a == 'A2' or a == 'A3') and (c == 'C1' or c == 'C3')]]
        assert_frame_equal(result, expected)

        result = df.loc(axis='index')[:,:,['C1','C3']]
        expected = df.loc[[ tuple([a,b,c,d]) for a,b,c,d in df.index.values if (
            c == 'C1' or c == 'C3')]]
        assert_frame_equal(result, expected)

        # axis 1
        result = df.loc(axis=1)[:,'foo']
        expected = df.loc[:,(slice(None),'foo')]
        assert_frame_equal(result, expected)

        result = df.loc(axis='columns')[:,'foo']
        expected = df.loc[:,(slice(None),'foo')]
        assert_frame_equal(result, expected)

        # invalid axis
        def f():
            df.loc(axis=-1)[:,:,['C1','C3']]
        self.assertRaises(ValueError, f)

        def f():
            df.loc(axis=2)[:,:,['C1','C3']]
        self.assertRaises(ValueError, f)

        def f():
            df.loc(axis='foo')[:,:,['C1','C3']]
        self.assertRaises(ValueError, f)

    def test_per_axis_per_level_setitem(self):

        # test index maker
        idx = pd.IndexSlice

        # test multi-index slicing with per axis and per index controls
        index = MultiIndex.from_tuples([('A',1),('A',2),('A',3),('B',1)],
                                       names=['one','two'])
        columns = MultiIndex.from_tuples([('a','foo'),('a','bar'),('b','foo'),('b','bah')],
                                         names=['lvl0', 'lvl1'])

        df_orig = DataFrame(np.arange(16,dtype='int64').reshape(4, 4), index=index, columns=columns)
        df_orig = df_orig.sortlevel(axis=0).sortlevel(axis=1)

        # identity
        df = df_orig.copy()
        df.loc[(slice(None),slice(None)),:] = 100
        expected = df_orig.copy()
        expected.iloc[:,:] = 100
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc(axis=0)[:,:] = 100
        expected = df_orig.copy()
        expected.iloc[:,:] = 100
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[(slice(None),slice(None)),(slice(None),slice(None))] = 100
        expected = df_orig.copy()
        expected.iloc[:,:] = 100
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[:,(slice(None),slice(None))] = 100
        expected = df_orig.copy()
        expected.iloc[:,:] = 100
        assert_frame_equal(df, expected)

        # index
        df = df_orig.copy()
        df.loc[(slice(None),[1]),:] = 100
        expected = df_orig.copy()
        expected.iloc[[0,3]] = 100
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[(slice(None),1),:] = 100
        expected = df_orig.copy()
        expected.iloc[[0,3]] = 100
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc(axis=0)[:,1] = 100
        expected = df_orig.copy()
        expected.iloc[[0,3]] = 100
        assert_frame_equal(df, expected)

        # columns
        df = df_orig.copy()
        df.loc[:,(slice(None),['foo'])] = 100
        expected = df_orig.copy()
        expected.iloc[:,[1,3]] = 100
        assert_frame_equal(df, expected)

        # both
        df = df_orig.copy()
        df.loc[(slice(None),1),(slice(None),['foo'])] = 100
        expected = df_orig.copy()
        expected.iloc[[0,3],[1,3]] = 100
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[idx[:,1],idx[:,['foo']]] = 100
        expected = df_orig.copy()
        expected.iloc[[0,3],[1,3]] = 100
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc['A','a'] = 100
        expected = df_orig.copy()
        expected.iloc[0:3,0:2] = 100
        assert_frame_equal(df, expected)

        # setting with a list-like
        df = df_orig.copy()
        df.loc[(slice(None),1),(slice(None),['foo'])] = np.array([[100, 100], [100, 100]],dtype='int64')
        expected = df_orig.copy()
        expected.iloc[[0,3],[1,3]] = 100
        assert_frame_equal(df, expected)

        # not enough values
        df = df_orig.copy()
        def f():
            df.loc[(slice(None),1),(slice(None),['foo'])] = np.array([[100], [100, 100]],dtype='int64')
        self.assertRaises(ValueError, f)
        def f():
            df.loc[(slice(None),1),(slice(None),['foo'])] = np.array([100, 100, 100, 100],dtype='int64')
        self.assertRaises(ValueError, f)

        # with an alignable rhs
        df = df_orig.copy()
        df.loc[(slice(None),1),(slice(None),['foo'])] = df.loc[(slice(None),1),(slice(None),['foo'])] * 5
        expected = df_orig.copy()
        expected.iloc[[0,3],[1,3]] = expected.iloc[[0,3],[1,3]] * 5
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[(slice(None),1),(slice(None),['foo'])] *= df.loc[(slice(None),1),(slice(None),['foo'])]
        expected = df_orig.copy()
        expected.iloc[[0,3],[1,3]] *= expected.iloc[[0,3],[1,3]]
        assert_frame_equal(df, expected)

        rhs = df_orig.loc[(slice(None),1),(slice(None),['foo'])].copy()
        rhs.loc[:,('c','bah')] = 10
        df = df_orig.copy()
        df.loc[(slice(None),1),(slice(None),['foo'])] *= rhs
        expected = df_orig.copy()
        expected.iloc[[0,3],[1,3]] *= expected.iloc[[0,3],[1,3]]
        assert_frame_equal(df, expected)

    def test_multiindex_setitem(self):

        # GH 3738
        # setting with a multi-index right hand side
        arrays = [np.array(['bar', 'bar', 'baz', 'qux', 'qux', 'bar']),
                  np.array(['one', 'two', 'one', 'one', 'two', 'one']),
                  np.arange(0, 6, 1)]

        df_orig = pd.DataFrame(np.random.randn(6, 3),
                               index=arrays,
                               columns=['A', 'B', 'C']).sort_index()

        expected = df_orig.loc[['bar']]*2
        df = df_orig.copy()
        df.loc[['bar']] *= 2
        assert_frame_equal(df.loc[['bar']],expected)

        # raise because these have differing levels
        def f():
            df.loc['bar'] *= 2
        self.assertRaises(TypeError, f)

        # from SO
        #http://stackoverflow.com/questions/24572040/pandas-access-the-level-of-multiindex-for-inplace-operation
        df_orig = DataFrame.from_dict({'price': {
            ('DE', 'Coal', 'Stock'): 2,
            ('DE', 'Gas', 'Stock'): 4,
            ('DE', 'Elec', 'Demand'): 1,
            ('FR', 'Gas', 'Stock'): 5,
            ('FR', 'Solar', 'SupIm'): 0,
            ('FR', 'Wind', 'SupIm'): 0}})
        df_orig.index = MultiIndex.from_tuples(df_orig.index, names=['Sit', 'Com', 'Type'])

        expected = df_orig.copy()
        expected.iloc[[0,2,3]] *= 2

        idx = pd.IndexSlice
        df = df_orig.copy()
        df.loc[idx[:,:,'Stock'],:] *= 2
        assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[idx[:,:,'Stock'],'price'] *= 2
        assert_frame_equal(df, expected)

    def test_getitem_multiindex(self):

        # GH 5725
        # the 'A' happens to be a valid Timestamp so the doesn't raise the appropriate
        # error, only in PY3 of course!
        index = MultiIndex(levels=[['D', 'B', 'C'], [0, 26, 27, 37, 57, 67, 75, 82]],
                           labels=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2], [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                           names=['tag', 'day'])
        arr = np.random.randn(len(index),1)
        df = DataFrame(arr,index=index,columns=['val'])
        result = df.val['D']
        expected = Series(arr.ravel()[0:3],name='val',index=Index([26,37,57],name='day'))
        assert_series_equal(result,expected)

        def f():
            df.val['A']
        self.assertRaises(KeyError, f)

        def f():
            df.val['X']
        self.assertRaises(KeyError, f)

        # A is treated as a special Timestamp
        index = MultiIndex(levels=[['A', 'B', 'C'], [0, 26, 27, 37, 57, 67, 75, 82]],
                           labels=[[0, 0, 0, 1, 2, 2, 2, 2, 2, 2], [1, 3, 4, 6, 0, 2, 2, 3, 5, 7]],
                           names=['tag', 'day'])
        df = DataFrame(arr,index=index,columns=['val'])
        result = df.val['A']
        expected = Series(arr.ravel()[0:3],name='val',index=Index([26,37,57],name='day'))
        assert_series_equal(result,expected)

        def f():
            df.val['X']
        self.assertRaises(KeyError, f)


        # GH 7866
        # multi-index slicing with missing indexers
        s = pd.Series(np.arange(9,dtype='int64'),
                      index=pd.MultiIndex.from_product([['A','B','C'],['foo','bar','baz']],
                                                       names=['one','two'])
                      ).sortlevel()

        expected = pd.Series(np.arange(3,dtype='int64'),
                             index=pd.MultiIndex.from_product([['A'],['foo','bar','baz']],
                                                              names=['one','two'])
                             ).sortlevel()

        result = s.loc[['A']]
        assert_series_equal(result,expected)
        result = s.loc[['A','D']]
        assert_series_equal(result,expected)

        # not any values found
        self.assertRaises(KeyError, lambda : s.loc[['D']])

        # empty ok
        result = s.loc[[]]
        expected = s.iloc[[]]
        assert_series_equal(result,expected)

        idx = pd.IndexSlice
        expected = pd.Series([0,3,6],
                             index=pd.MultiIndex.from_product([['A','B','C'],['foo']],
                                                              names=['one','two'])
                             ).sortlevel()

        result = s.loc[idx[:,['foo']]]
        assert_series_equal(result,expected)
        result = s.loc[idx[:,['foo','bah']]]
        assert_series_equal(result,expected)

        # GH 8737
        # empty indexer
        multi_index = pd.MultiIndex.from_product((['foo', 'bar', 'baz'], ['alpha', 'beta']))
        df = DataFrame(np.random.randn(5, 6), index=range(5), columns=multi_index)
        df = df.sortlevel(0, axis=1)

        expected = DataFrame(index=range(5),columns=multi_index.reindex([])[0])
        result1 = df.loc[:, ([], slice(None))]
        result2 = df.loc[:, (['foo'], [])]
        assert_frame_equal(result1, expected)
        assert_frame_equal(result2, expected)

        # regression from < 0.14.0
        # GH 7914
        df = DataFrame([[np.mean, np.median],['mean','median']],
                       columns=MultiIndex.from_tuples([('functs','mean'),
                                                       ('functs','median')]),
                       index=['function', 'name'])
        result = df.loc['function',('functs','mean')]
        self.assertEqual(result,np.mean)

    def test_setitem_dtype_upcast(self):

        # GH3216
        df = DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df['c'] = np.nan
        self.assertEqual(df['c'].dtype, np.float64)

        df.ix[0,'c'] = 'foo'
        expected = DataFrame([{"a": 1, "c" : 'foo'}, {"a": 3, "b": 2, "c" : np.nan}])
        assert_frame_equal(df,expected)

        # GH10280
        df = DataFrame(np.arange(6,dtype='int64').reshape(2, 3),
                       index=list('ab'),
                       columns=['foo', 'bar', 'baz'])

        for val in [3.14, 'wxyz']:
            left = df.copy()
            left.loc['a', 'bar'] = val
            right = DataFrame([[0, val, 2], [3, 4, 5]],
                              index=list('ab'),
                              columns=['foo', 'bar', 'baz'])

            assert_frame_equal(left, right)
            self.assertTrue(com.is_integer_dtype(left['foo']))
            self.assertTrue(com.is_integer_dtype(left['baz']))

        left = DataFrame(np.arange(6,dtype='int64').reshape(2, 3) / 10.0,
                         index=list('ab'),
                         columns=['foo', 'bar', 'baz'])
        left.loc['a', 'bar'] = 'wxyz'

        right = DataFrame([[0, 'wxyz', .2], [.3, .4, .5]],
                          index=list('ab'),
                          columns=['foo', 'bar', 'baz'])

        assert_frame_equal(left, right)
        self.assertTrue(com.is_float_dtype(left['foo']))
        self.assertTrue(com.is_float_dtype(left['baz']))

    def test_setitem_iloc(self):


        # setitem with an iloc list
        df = DataFrame(np.arange(9).reshape((3, 3)), index=["A", "B", "C"], columns=["A", "B", "C"])
        df.iloc[[0,1],[1,2]]
        df.iloc[[0,1],[1,2]] += 100

        expected = DataFrame(np.array([0,101,102,3,104,105,6,7,8]).reshape((3, 3)), index=["A", "B", "C"], columns=["A", "B", "C"])
        assert_frame_equal(df,expected)

    def test_dups_fancy_indexing(self):

        # GH 3455
        from pandas.util.testing import makeCustomDataframe as mkdf
        df= mkdf(10, 3)
        df.columns = ['a','a','b']
        cols = ['b','a']
        result = df[['b','a']].columns
        expected = Index(['b','a','a'])
        self.assertTrue(result.equals(expected))

        # across dtypes
        df = DataFrame([[1,2,1.,2.,3.,'foo','bar']], columns=list('aaaaaaa'))
        df.head()
        str(df)
        result = DataFrame([[1,2,1.,2.,3.,'foo','bar']])
        result.columns = list('aaaaaaa')

        df_v  = df.iloc[:,4]
        res_v = result.iloc[:,4]

        assert_frame_equal(df,result)

        # GH 3561, dups not in selected order
        df = DataFrame({'test': [5,7,9,11], 'test1': [4.,5,6,7], 'other': list('abcd') }, index=['A', 'A', 'B', 'C'])
        rows = ['C', 'B']
        expected = DataFrame({'test' : [11,9], 'test1': [ 7., 6], 'other': ['d','c']},index=rows)
        result = df.ix[rows]
        assert_frame_equal(result, expected)

        result = df.ix[Index(rows)]
        assert_frame_equal(result, expected)

        rows = ['C','B','E']
        expected = DataFrame({'test' : [11,9,np.nan], 'test1': [7.,6,np.nan], 'other': ['d','c',np.nan]},index=rows)

        result = df.ix[rows]
        assert_frame_equal(result, expected)

        # see GH5553, make sure we use the right indexer
        rows = ['F','G','H','C','B','E']
        expected = DataFrame({'test' : [np.nan,np.nan,np.nan,11,9,np.nan],
                              'test1': [np.nan,np.nan,np.nan,7.,6,np.nan],
                              'other': [np.nan,np.nan,np.nan,'d','c',np.nan]},index=rows)
        result = df.ix[rows]
        assert_frame_equal(result, expected)

        # inconsistent returns for unique/duplicate indices when values are missing
        df = DataFrame(randn(4,3),index=list('ABCD'))
        expected = df.ix[['E']]

        dfnu = DataFrame(randn(5,3),index=list('AABCD'))
        result = dfnu.ix[['E']]
        assert_frame_equal(result, expected)

        # ToDo: check_index_type can be True after GH 11497

        # GH 4619; duplicate indexer with missing label
        df = DataFrame({"A": [0, 1, 2]})
        result = df.ix[[0, 8, 0]]
        expected = DataFrame({"A": [0, np.nan, 0]}, index=[0, 8, 0])
        assert_frame_equal(result, expected, check_index_type=False)

        df = DataFrame({"A": list('abc')})
        result = df.ix[[0,8,0]]
        expected = DataFrame({"A": ['a', np.nan, 'a']}, index=[0, 8, 0])
        assert_frame_equal(result, expected, check_index_type=False)

        # non unique with non unique selector
        df = DataFrame({'test': [5, 7, 9, 11]}, index=['A', 'A', 'B', 'C'])
        expected = DataFrame({'test' : [5, 7, 5, 7, np.nan]}, index=['A', 'A', 'A', 'A', 'E'])
        result = df.ix[['A', 'A', 'E']]
        assert_frame_equal(result, expected)

        # GH 5835
        # dups on index and missing values
        df = DataFrame(np.random.randn(5, 5), columns=['A', 'B', 'B', 'B', 'A'])

        expected = pd.concat([df.ix[:,['A','B']],DataFrame(np.nan,columns=['C'],index=df.index)],axis=1)
        result = df.ix[:,['A','B','C']]
        assert_frame_equal(result, expected)

        # GH 6504, multi-axis indexing
        df = DataFrame(np.random.randn(9,2), index=[1,1,1,2,2,2,3,3,3], columns=['a', 'b'])

        expected = df.iloc[0:6]
        result = df.loc[[1, 2]]
        assert_frame_equal(result, expected)

        expected = df
        result = df.loc[:,['a', 'b']]
        assert_frame_equal(result, expected)

        expected = df.iloc[0:6,:]
        result = df.loc[[1, 2], ['a', 'b']]
        assert_frame_equal(result, expected)

    def test_indexing_mixed_frame_bug(self):

        # GH3492
        df=DataFrame({'a':{1:'aaa',2:'bbb',3:'ccc'},'b':{1:111,2:222,3:333}})

        # this works, new column is created correctly
        df['test']=df['a'].apply(lambda x: '_' if x=='aaa' else x)

        # this does not work, ie column test is not changed
        idx=df['test']=='_'
        temp=df.ix[idx,'a'].apply(lambda x: '-----' if x=='aaa' else x)
        df.ix[idx,'test']=temp
        self.assertEqual(df.iloc[0,2], '-----')

        #if I look at df, then element [0,2] equals '_'. If instead I type df.ix[idx,'test'], I get '-----', finally by typing df.iloc[0,2] I get '_'.

    def test_multitype_list_index_access(self):
        #GH 10610
        df = pd.DataFrame(np.random.random((10, 5)), columns=["a"] + [20, 21, 22, 23])
        with self.assertRaises(IndexError):
            vals = df[[22, 26, -8]]
        self.assertEqual(df[21].shape[0], df.shape[0])

    def test_set_index_nan(self):

        # GH 3586
        df = DataFrame({'PRuid': {17: 'nonQC', 18: 'nonQC', 19: 'nonQC', 20: '10', 21: '11', 22: '12', 23: '13',
                                  24: '24', 25: '35', 26: '46', 27: '47', 28: '48', 29: '59', 30: '10'},
                        'QC': {17: 0.0, 18: 0.0, 19: 0.0, 20: nan, 21: nan, 22: nan, 23: nan, 24: 1.0, 25: nan,
                               26: nan, 27: nan, 28: nan, 29: nan, 30: nan},
                        'data': {17: 7.9544899999999998, 18: 8.0142609999999994, 19: 7.8591520000000008, 20: 0.86140349999999999,
                                 21: 0.87853110000000001, 22: 0.8427041999999999, 23: 0.78587700000000005, 24: 0.73062459999999996,
                                 25: 0.81668560000000001, 26: 0.81927080000000008, 27: 0.80705009999999999, 28: 0.81440240000000008,
                                 29: 0.80140849999999997, 30: 0.81307740000000006},
                        'year': {17: 2006, 18: 2007, 19: 2008, 20: 1985, 21: 1985, 22: 1985, 23: 1985,
                                 24: 1985, 25: 1985, 26: 1985, 27: 1985, 28: 1985, 29: 1985, 30: 1986}}).reset_index()

        result = df.set_index(['year','PRuid','QC']).reset_index().reindex(columns=df.columns)
        assert_frame_equal(result,df)

    def test_multi_nan_indexing(self):

        # GH 3588
        df = DataFrame({"a":['R1', 'R2', np.nan, 'R4'], 'b':["C1", "C2", "C3" , "C4"], "c":[10, 15, np.nan , 20]})
        result = df.set_index(['a','b'], drop=False)
        expected = DataFrame({"a":['R1', 'R2', np.nan, 'R4'], 'b':["C1", "C2", "C3" , "C4"], "c":[10, 15, np.nan , 20]},
                             index = [Index(['R1','R2',np.nan,'R4'],name='a'),Index(['C1','C2','C3','C4'],name='b')])
        assert_frame_equal(result,expected)


    def test_iloc_panel_issue(self):

        # GH 3617
        p = Panel(randn(4, 4, 4))

        self.assertEqual(p.iloc[:3, :3, :3].shape, (3,3,3))
        self.assertEqual(p.iloc[1, :3, :3].shape, (3,3))
        self.assertEqual(p.iloc[:3, 1, :3].shape, (3,3))
        self.assertEqual(p.iloc[:3, :3, 1].shape, (3,3))
        self.assertEqual(p.iloc[1, 1, :3].shape, (3,))
        self.assertEqual(p.iloc[1, :3, 1].shape, (3,))
        self.assertEqual(p.iloc[:3, 1, 1].shape, (3,))

    def test_panel_getitem(self):
        # GH4016, date selection returns a frame when a partial string selection
        ind = date_range(start="2000", freq="D", periods=1000)
        df = DataFrame(np.random.randn(len(ind), 5), index=ind, columns=list('ABCDE'))
        panel = Panel(dict([ ('frame_'+c,df) for c in list('ABC') ]))

        test2 = panel.ix[:, "2002":"2002-12-31"]
        test1 = panel.ix[:, "2002"]
        tm.assert_panel_equal(test1,test2)

        # GH8710
        # multi-element getting with a list
        panel = tm.makePanel()

        expected = panel.iloc[[0,1]]

        result = panel.loc[['ItemA','ItemB']]
        tm.assert_panel_equal(result,expected)

        result = panel.loc[['ItemA','ItemB'],:,:]
        tm.assert_panel_equal(result,expected)

        result = panel[['ItemA','ItemB']]
        tm.assert_panel_equal(result,expected)

        result = panel.loc['ItemA':'ItemB']
        tm.assert_panel_equal(result,expected)

        result = panel.ix['ItemA':'ItemB']
        tm.assert_panel_equal(result,expected)

        result = panel.ix[['ItemA','ItemB']]
        tm.assert_panel_equal(result,expected)

        # with an object-like
        # GH 9140
        class TestObject:
            def __str__(self):
                return "TestObject"

        obj = TestObject()

        p = Panel(np.random.randn(1,5,4), items=[obj],
                  major_axis = date_range('1/1/2000', periods=5),
                  minor_axis=['A', 'B', 'C', 'D'])

        expected = p.iloc[0]
        result = p[obj]
        tm.assert_frame_equal(result, expected)

    def test_panel_setitem(self):

        # GH 7763
        # loc and setitem have setting differences
        np.random.seed(0)
        index=range(3)
        columns = list('abc')

        panel = Panel({'A' : DataFrame(np.random.randn(3, 3), index=index, columns=columns),
                       'B' : DataFrame(np.random.randn(3, 3), index=index, columns=columns),
                       'C' : DataFrame(np.random.randn(3, 3), index=index, columns=columns)
                       })

        replace = DataFrame(np.eye(3,3), index=range(3), columns=columns)
        expected = Panel({ 'A' : replace, 'B' : replace, 'C' : replace })

        p = panel.copy()
        for idx in list('ABC'):
            p[idx] = replace
        tm.assert_panel_equal(p, expected)

        p = panel.copy()
        for idx in list('ABC'):
            p.loc[idx,:,:] = replace
        tm.assert_panel_equal(p, expected)


    def test_panel_setitem_with_multiindex(self):

        # 10360
        # failing with a multi-index
        arr = np.array([[[1,2,3],[0,0,0]],[[0,0,0],[0,0,0]]],dtype=np.float64)

        # reg index
        axes = dict(items=['A', 'B'], major_axis=[0, 1], minor_axis=['X', 'Y' ,'Z'])
        p1 = Panel(0., **axes)
        p1.iloc[0, 0, :] = [1, 2, 3]
        expected = Panel(arr, **axes)
        tm.assert_panel_equal(p1, expected)

        # multi-indexes
        axes['items'] = pd.MultiIndex.from_tuples([('A','a'), ('B','b')])
        p2 = Panel(0., **axes)
        p2.iloc[0, 0, :] = [1, 2, 3]
        expected = Panel(arr, **axes)
        tm.assert_panel_equal(p2, expected)

        axes['major_axis']=pd.MultiIndex.from_tuples([('A',1),('A',2)])
        p3 = Panel(0., **axes)
        p3.iloc[0, 0, :] = [1, 2, 3]
        expected = Panel(arr, **axes)
        tm.assert_panel_equal(p3, expected)

        axes['minor_axis']=pd.MultiIndex.from_product([['X'],range(3)])
        p4 = Panel(0., **axes)
        p4.iloc[0, 0, :] = [1, 2, 3]
        expected = Panel(arr, **axes)
        tm.assert_panel_equal(p4, expected)

        arr = np.array([[[1,0,0],[2,0,0]],[[0,0,0],[0,0,0]]],dtype=np.float64)
        p5 = Panel(0., **axes)
        p5.iloc[0, :, 0] = [1, 2]
        expected = Panel(arr, **axes)
        tm.assert_panel_equal(p5, expected)

    def test_panel_assignment(self):

        # GH3777
        wp = Panel(randn(2, 5, 4), items=['Item1', 'Item2'], major_axis=date_range('1/1/2000', periods=5), minor_axis=['A', 'B', 'C', 'D'])
        wp2 = Panel(randn(2, 5, 4), items=['Item1', 'Item2'], major_axis=date_range('1/1/2000', periods=5), minor_axis=['A', 'B', 'C', 'D'])
        expected = wp.loc[['Item1', 'Item2'], :, ['A', 'B']]

        def f():
            wp.loc[['Item1', 'Item2'], :, ['A', 'B']] = wp2.loc[['Item1', 'Item2'], :, ['A', 'B']]
        self.assertRaises(NotImplementedError, f)

        #wp.loc[['Item1', 'Item2'], :, ['A', 'B']] = wp2.loc[['Item1', 'Item2'], :, ['A', 'B']]
        #result = wp.loc[['Item1', 'Item2'], :, ['A', 'B']]
        #tm.assert_panel_equal(result,expected)

    def test_multiindex_assignment(self):

        # GH3777 part 2

        # mixed dtype
        df = DataFrame(np.random.randint(5,10,size=9).reshape(3, 3),
                       columns=list('abc'),
                       index=[[4,4,8],[8,10,12]])
        df['d'] = np.nan
        arr = np.array([0.,1.])

        df.ix[4,'d'] = arr
        assert_series_equal(df.ix[4,'d'],Series(arr,index=[8,10],name='d'))

        # single dtype
        df = DataFrame(np.random.randint(5,10,size=9).reshape(3, 3),
                       columns=list('abc'),
                       index=[[4,4,8],[8,10,12]])

        df.ix[4,'c'] = arr
        assert_series_equal(df.ix[4,'c'],Series(arr,index=[8,10],name='c',dtype='int64'))

        # scalar ok
        df.ix[4,'c'] = 10
        assert_series_equal(df.ix[4,'c'],Series(10,index=[8,10],name='c',dtype='int64'))

        # invalid assignments
        def f():
            df.ix[4,'c'] = [0,1,2,3]
        self.assertRaises(ValueError, f)

        def f():
            df.ix[4,'c'] = [0]
        self.assertRaises(ValueError, f)

        # groupby example
        NUM_ROWS = 100
        NUM_COLS = 10
        col_names = ['A'+num for num in map(str,np.arange(NUM_COLS).tolist())]
        index_cols = col_names[:5]

        df = DataFrame(np.random.randint(5, size=(NUM_ROWS,NUM_COLS)), dtype=np.int64, columns=col_names)
        df = df.set_index(index_cols).sort_index()
        grp = df.groupby(level=index_cols[:4])
        df['new_col'] = np.nan

        f_index = np.arange(5)
        def f(name,df2):
            return Series(np.arange(df2.shape[0]),name=df2.index.values[0]).reindex(f_index)
        new_df = pd.concat([ f(name,df2) for name, df2 in grp ],axis=1).T

        # we are actually operating on a copy here
        # but in this case, that's ok
        for name, df2 in grp:
            new_vals = np.arange(df2.shape[0])
            df.ix[name, 'new_col'] = new_vals

    def test_multi_assign(self):

        # GH 3626, an assignement of a sub-df to a df
        df = DataFrame({'FC':['a','b','a','b','a','b'],
                        'PF':[0,0,0,0,1,1],
                        'col1':lrange(6),
                        'col2':lrange(6,12)})
        df.ix[1,0]=np.nan
        df2 = df.copy()

        mask=~df2.FC.isnull()
        cols=['col1', 'col2']

        dft = df2 * 2
        dft.ix[3,3] = np.nan

        expected = DataFrame({'FC':['a',np.nan,'a','b','a','b'],
                              'PF':[0,0,0,0,1,1],
                              'col1':Series([0,1,4,6,8,10]),
                              'col2':[12,7,16,np.nan,20,22]})


        # frame on rhs
        df2.ix[mask, cols]= dft.ix[mask, cols]
        assert_frame_equal(df2,expected)

        df2.ix[mask, cols]= dft.ix[mask, cols]
        assert_frame_equal(df2,expected)

        # with an ndarray on rhs
        df2 = df.copy()
        df2.ix[mask, cols]= dft.ix[mask, cols].values
        assert_frame_equal(df2,expected)
        df2.ix[mask, cols]= dft.ix[mask, cols].values
        assert_frame_equal(df2,expected)

        # broadcasting on the rhs is required
        df = DataFrame(dict(A = [1,2,0,0,0],B=[0,0,0,10,11],C=[0,0,0,10,11],D=[3,4,5,6,7]))

        expected = df.copy()
        mask = expected['A'] == 0
        for col in ['A','B']:
            expected.loc[mask,col] = df['D']

        df.loc[df['A']==0,['A','B']] = df['D']
        assert_frame_equal(df,expected)

    def test_ix_assign_column_mixed(self):
        # GH #1142
        df = DataFrame(tm.getSeriesData())
        df['foo'] = 'bar'

        orig = df.ix[:, 'B'].copy()
        df.ix[:, 'B'] = df.ix[:, 'B'] + 1
        assert_series_equal(df.B, orig + 1)

        # GH 3668, mixed frame with series value
        df = DataFrame({'x':lrange(10), 'y':lrange(10,20),'z' : 'bar'})
        expected = df.copy()

        for i in range(5):
            indexer = i*2
            v = 1000 + i*200
            expected.ix[indexer, 'y'] = v
            self.assertEqual(expected.ix[indexer, 'y'], v)

        df.ix[df.x % 2 == 0, 'y'] = df.ix[df.x % 2 == 0, 'y'] * 100
        assert_frame_equal(df,expected)

        # GH 4508, making sure consistency of assignments
        df = DataFrame({'a':[1,2,3],'b':[0,1,2]})
        df.ix[[0,2,],'b'] = [100,-100]
        expected = DataFrame({'a' : [1,2,3], 'b' : [100,1,-100] })
        assert_frame_equal(df,expected)

        df = pd.DataFrame({'a': lrange(4) })
        df['b'] = np.nan
        df.ix[[1,3],'b'] = [100,-100]
        expected = DataFrame({'a' : [0,1,2,3], 'b' : [np.nan,100,np.nan,-100] })
        assert_frame_equal(df,expected)

        # ok, but chained assignments are dangerous
        # if we turn off chained assignement it will work
        with option_context('chained_assignment',None):
            df = pd.DataFrame({'a': lrange(4) })
            df['b'] = np.nan
            df['b'].ix[[1,3]] = [100,-100]
            assert_frame_equal(df,expected)

    def test_ix_get_set_consistency(self):

        # GH 4544
        # ix/loc get/set not consistent when
        # a mixed int/string index
        df = DataFrame(np.arange(16).reshape((4, 4)),
                       columns=['a', 'b', 8, 'c'],
                       index=['e', 7, 'f', 'g'])

        self.assertEqual(df.ix['e', 8], 2)
        self.assertEqual(df.loc['e', 8], 2)

        df.ix['e', 8] = 42
        self.assertEqual(df.ix['e', 8], 42)
        self.assertEqual(df.loc['e', 8], 42)

        df.loc['e', 8] = 45
        self.assertEqual(df.ix['e', 8], 45)
        self.assertEqual(df.loc['e', 8], 45)

    def test_setitem_list(self):

        # GH 6043
        # ix with a list
        df = DataFrame(index=[0,1], columns=[0])
        df.ix[1,0] = [1,2,3]
        df.ix[1,0] = [1,2]

        result = DataFrame(index=[0,1], columns=[0])
        result.ix[1,0] = [1,2]

        assert_frame_equal(result,df)

        # ix with an object
        class TO(object):
            def __init__(self, value):
                self.value = value
            def __str__(self):
                return "[{0}]".format(self.value)
            __repr__ = __str__
            def __eq__(self, other):
                return self.value == other.value
            def view(self):
                return self

        df = DataFrame(index=[0,1], columns=[0])
        df.ix[1,0] = TO(1)
        df.ix[1,0] = TO(2)

        result = DataFrame(index=[0,1], columns=[0])
        result.ix[1,0] = TO(2)

        assert_frame_equal(result,df)

        # remains object dtype even after setting it back
        df = DataFrame(index=[0,1], columns=[0])
        df.ix[1,0] = TO(1)
        df.ix[1,0] = np.nan
        result = DataFrame(index=[0,1], columns=[0])

        assert_frame_equal(result, df)

    def test_iloc_mask(self):

        # GH 3631, iloc with a mask (of a series) should raise
        df = DataFrame(lrange(5), list('ABCDE'), columns=['a'])
        mask = (df.a%2 == 0)
        self.assertRaises(ValueError, df.iloc.__getitem__, tuple([mask]))
        mask.index = lrange(len(mask))
        self.assertRaises(NotImplementedError, df.iloc.__getitem__, tuple([mask]))

        # ndarray ok
        result = df.iloc[np.array([True] * len(mask),dtype=bool)]
        assert_frame_equal(result,df)

        # the possibilities
        locs = np.arange(4)
        nums = 2**locs
        reps = lmap(bin, nums)
        df = DataFrame({'locs':locs, 'nums':nums}, reps)

        expected = {
            (None,'')     : '0b1100',
            (None,'.loc')  : '0b1100',
            (None,'.iloc') : '0b1100',
            ('index','')  : '0b11',
            ('index','.loc')  : '0b11',
            ('index','.iloc') : 'iLocation based boolean indexing cannot use an indexable as a mask',
            ('locs','')      : 'Unalignable boolean Series key provided',
            ('locs','.loc')   : 'Unalignable boolean Series key provided',
            ('locs','.iloc')  : 'iLocation based boolean indexing on an integer type is not available',
            }

        warnings.filterwarnings(action='ignore', category=UserWarning)
        result = dict()
        for idx in [None, 'index', 'locs']:
            mask = (df.nums>2).values
            if idx:
                mask = Series(mask, list(reversed(getattr(df, idx))))
            for method in ['', '.loc', '.iloc']:
                try:
                    if method:
                        accessor = getattr(df, method[1:])
                    else:
                        accessor = df
                    ans = str(bin(accessor[mask]['nums'].sum()))
                except Exception as e:
                    ans = str(e)

                key = tuple([idx,method])
                r = expected.get(key)
                if r != ans:
                    raise AssertionError("[%s] does not match [%s], received [%s]" %
                                         (key,ans,r))
        warnings.filterwarnings(action='always', category=UserWarning)

    def test_ix_slicing_strings(self):
        ##GH3836
        data = {'Classification': ['SA EQUITY CFD', 'bbb', 'SA EQUITY', 'SA SSF', 'aaa'],
                'Random': [1,2,3,4,5],
                'X': ['correct', 'wrong','correct', 'correct','wrong']}
        df = DataFrame(data)
        x = df[~df.Classification.isin(['SA EQUITY CFD', 'SA EQUITY', 'SA SSF'])]
        df.ix[x.index,'X'] = df['Classification']

        expected = DataFrame({'Classification': {0: 'SA EQUITY CFD', 1: 'bbb',
                                                2: 'SA EQUITY', 3: 'SA SSF', 4: 'aaa'},
                            'Random': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5},
                            'X': {0: 'correct', 1: 'bbb', 2: 'correct',
                            3: 'correct', 4: 'aaa'}})  # bug was 4: 'bbb'

        assert_frame_equal(df, expected)

    def test_non_unique_loc(self):
        ## GH3659
        ## non-unique indexer with loc slice
        ## https://groups.google.com/forum/?fromgroups#!topic/pydata/zTm2No0crYs

        # these are going to raise becuase the we are non monotonic
        df = DataFrame({'A' : [1,2,3,4,5,6], 'B' : [3,4,5,6,7,8]}, index = [0,1,0,1,2,3])
        self.assertRaises(KeyError, df.loc.__getitem__, tuple([slice(1,None)]))
        self.assertRaises(KeyError, df.loc.__getitem__, tuple([slice(0,None)]))
        self.assertRaises(KeyError, df.loc.__getitem__, tuple([slice(1,2)]))

        # monotonic are ok
        df = DataFrame({'A' : [1,2,3,4,5,6], 'B' : [3,4,5,6,7,8]}, index = [0,1,0,1,2,3]).sort_index(axis=0)
        result = df.loc[1:]
        expected = DataFrame({'A' : [2,4,5,6], 'B' : [4, 6,7,8]}, index = [1,1,2,3])
        assert_frame_equal(result,expected)

        result = df.loc[0:]
        assert_frame_equal(result,df)

        result = df.loc[1:2]
        expected = DataFrame({'A' : [2,4,5], 'B' : [4,6,7]}, index = [1,1,2])
        assert_frame_equal(result,expected)

    def test_loc_name(self):
        # GH 3880
        df = DataFrame([[1, 1], [1, 1]])
        df.index.name = 'index_name'
        result = df.iloc[[0, 1]].index.name
        self.assertEqual(result, 'index_name')

        result = df.ix[[0, 1]].index.name
        self.assertEqual(result, 'index_name')

        result = df.loc[[0, 1]].index.name
        self.assertEqual(result, 'index_name')

    def test_iloc_non_unique_indexing(self):

        #GH 4017, non-unique indexing (on the axis)
        df = DataFrame({'A' : [0.1] * 3000, 'B' : [1] * 3000})
        idx = np.array(lrange(30)) * 99
        expected = df.iloc[idx]

        df3 = pd.concat([df, 2*df, 3*df])
        result = df3.iloc[idx]

        assert_frame_equal(result, expected)

        df2 = DataFrame({'A' : [0.1] * 1000, 'B' : [1] * 1000})
        df2 = pd.concat([df2, 2*df2, 3*df2])

        sidx = df2.index.to_series()
        expected = df2.iloc[idx[idx<=sidx.max()]]

        new_list = []
        for r, s in expected.iterrows():
            new_list.append(s)
            new_list.append(s*2)
            new_list.append(s*3)

        expected = DataFrame(new_list)
        expected = pd.concat([ expected, DataFrame(index=idx[idx>sidx.max()]) ])
        result = df2.loc[idx]
        assert_frame_equal(result, expected, check_index_type=False)

    def test_mi_access(self):

        # GH 4145
        data = """h1 main  h3 sub  h5
0  a    A   1  A1   1
1  b    B   2  B1   2
2  c    B   3  A1   3
3  d    A   4  B2   4
4  e    A   5  B2   5
5  f    B   6  A2   6
"""

        df = pd.read_csv(StringIO(data),sep='\s+',index_col=0)
        df2 = df.set_index(['main', 'sub']).T.sort_index(1)
        index = Index(['h1','h3','h5'])
        columns = MultiIndex.from_tuples([('A','A1')],names=['main','sub'])
        expected = DataFrame([['a',1,1]],index=columns,columns=index).T

        result = df2.loc[:,('A','A1')]
        assert_frame_equal(result,expected)

        result = df2[('A','A1')]
        assert_frame_equal(result,expected)

        # GH 4146, not returning a block manager when selecting a unique index
        # from a duplicate index
        # as of 4879, this returns a Series (which is similar to what happens with a non-unique)
        expected = Series(['a',1,1], index=['h1','h3','h5'], name='A1')
        result = df2['A']['A1']
        assert_series_equal(result, expected)

        # selecting a non_unique from the 2nd level
        expected = DataFrame([['d',4,4],['e',5,5]],index=Index(['B2','B2'],name='sub'),columns=['h1','h3','h5'],).T
        result = df2['A']['B2']
        assert_frame_equal(result, expected)

    def test_non_unique_loc_memory_error(self):

        # GH 4280
        # non_unique index with a large selection triggers a memory error

        columns = list('ABCDEFG')
        def gen_test(l,l2):
            return pd.concat([ DataFrame(randn(l,len(columns)),index=lrange(l),columns=columns),
                               DataFrame(np.ones((l2,len(columns))),index=[0]*l2,columns=columns) ])


        def gen_expected(df,mask):
            l = len(mask)
            return pd.concat([
                df.take([0],convert=False),
                DataFrame(np.ones((l,len(columns))),index=[0]*l,columns=columns),
                df.take(mask[1:],convert=False) ])

        df = gen_test(900,100)
        self.assertFalse(df.index.is_unique)

        mask = np.arange(100)
        result = df.loc[mask]
        expected = gen_expected(df,mask)
        assert_frame_equal(result,expected)

        df = gen_test(900000,100000)
        self.assertFalse(df.index.is_unique)

        mask = np.arange(100000)
        result = df.loc[mask]
        expected = gen_expected(df,mask)
        assert_frame_equal(result,expected)

    def test_astype_assignment(self):

        # GH4312 (iloc)
        df_orig = DataFrame([['1','2','3','.4',5,6.,'foo']],columns=list('ABCDEFG'))

        df = df_orig.copy()
        df.iloc[:,0:2] = df.iloc[:,0:2].astype(np.int64)
        expected = DataFrame([[1,2,'3','.4',5,6.,'foo']],columns=list('ABCDEFG'))
        assert_frame_equal(df,expected)

        df = df_orig.copy()
        df.iloc[:,0:2] = df.iloc[:,0:2]._convert(datetime=True, numeric=True)
        expected =  DataFrame([[1,2,'3','.4',5,6.,'foo']],columns=list('ABCDEFG'))
        assert_frame_equal(df,expected)

        # GH5702 (loc)
        df = df_orig.copy()
        df.loc[:,'A'] = df.loc[:,'A'].astype(np.int64)
        expected = DataFrame([[1,'2','3','.4',5,6.,'foo']],columns=list('ABCDEFG'))
        assert_frame_equal(df,expected)

        df = df_orig.copy()
        df.loc[:,['B','C']] = df.loc[:,['B','C']].astype(np.int64)
        expected =  DataFrame([['1',2,3,'.4',5,6.,'foo']],columns=list('ABCDEFG'))
        assert_frame_equal(df,expected)

        # full replacements / no nans
        df = DataFrame({'A': [1., 2., 3., 4.]})
        df.iloc[:, 0] = df['A'].astype(np.int64)
        expected = DataFrame({'A': [1, 2, 3, 4]})
        assert_frame_equal(df,expected)

        df = DataFrame({'A': [1., 2., 3., 4.]})
        df.loc[:, 'A'] = df['A'].astype(np.int64)
        expected = DataFrame({'A': [1, 2, 3, 4]})
        assert_frame_equal(df,expected)

    def test_astype_assignment_with_dups(self):

        # GH 4686
        # assignment with dups that has a dtype change
        df = DataFrame(
            np.arange(3).reshape((1,3)),
            columns=pd.MultiIndex.from_tuples(
                [('A', '1'), ('B', '1'), ('A', '2')]
                ),
            dtype=object
            )
        index = df.index.copy()

        df['A'] = df['A'].astype(np.float64)
        result = df.get_dtype_counts().sort_index()
        expected = Series({ 'float64' : 2, 'object' : 1 }).sort_index()
        self.assertTrue(df.index.equals(index))

    def test_dups_loc(self):

        # GH4726
        # dup indexing with iloc/loc
        df = DataFrame([[1, 2, 'foo', 'bar', Timestamp('20130101')]],
                       columns=['a','a','a','a','a'], index=[1])
        expected = Series([1, 2, 'foo', 'bar', Timestamp('20130101')],
                          index=['a','a','a','a','a'], name=1)

        result = df.iloc[0]
        assert_series_equal(result, expected)

        result = df.loc[1]
        assert_series_equal(result, expected)

    def test_partial_setting(self):

        # GH2578, allow ix and friends to partially set

        ### series ###
        s_orig = Series([1,2,3])

        s = s_orig.copy()
        s[5] = 5
        expected = Series([1,2,3,5],index=[0,1,2,5])
        assert_series_equal(s,expected)

        s = s_orig.copy()
        s.loc[5] = 5
        expected = Series([1,2,3,5],index=[0,1,2,5])
        assert_series_equal(s,expected)

        s = s_orig.copy()
        s[5] = 5.
        expected = Series([1,2,3,5.],index=[0,1,2,5])
        assert_series_equal(s,expected)

        s = s_orig.copy()
        s.loc[5] = 5.
        expected = Series([1,2,3,5.],index=[0,1,2,5])
        assert_series_equal(s,expected)

        # iloc/iat raise
        s = s_orig.copy()
        def f():
            s.iloc[3] = 5.
        self.assertRaises(IndexError, f)
        def f():
            s.iat[3] = 5.
        self.assertRaises(IndexError, f)

        ### frame ###

        df_orig = DataFrame(np.arange(6).reshape(3,2),columns=['A','B'],dtype='int64')

        # iloc/iat raise
        df = df_orig.copy()
        def f():
            df.iloc[4,2] = 5.
        self.assertRaises(IndexError, f)
        def f():
            df.iat[4,2] = 5.
        self.assertRaises(IndexError, f)

        # row setting where it exists
        expected = DataFrame(dict({ 'A' : [0,4,4], 'B' : [1,5,5] }))
        df = df_orig.copy()
        df.iloc[1] = df.iloc[2]
        assert_frame_equal(df,expected)

        expected = DataFrame(dict({ 'A' : [0,4,4], 'B' : [1,5,5] }))
        df = df_orig.copy()
        df.loc[1] = df.loc[2]
        assert_frame_equal(df,expected)

        # like 2578, partial setting with dtype preservation
        expected = DataFrame(dict({ 'A' : [0,2,4,4], 'B' : [1,3,5,5] }))
        df = df_orig.copy()
        df.loc[3] = df.loc[2]
        assert_frame_equal(df,expected)

        # single dtype frame, overwrite
        expected = DataFrame(dict({ 'A' : [0,2,4], 'B' : [0,2,4] }))
        df = df_orig.copy()
        df.ix[:,'B'] = df.ix[:,'A']
        assert_frame_equal(df,expected)

        # mixed dtype frame, overwrite
        expected = DataFrame(dict({ 'A' : [0,2,4], 'B' : Series([0,2,4]) }))
        df = df_orig.copy()
        df['B'] = df['B'].astype(np.float64)
        df.ix[:,'B'] = df.ix[:,'A']
        assert_frame_equal(df,expected)

        # single dtype frame, partial setting
        expected = df_orig.copy()
        expected['C'] = df['A']
        df = df_orig.copy()
        df.ix[:,'C'] = df.ix[:,'A']
        assert_frame_equal(df,expected)

        # mixed frame, partial setting
        expected = df_orig.copy()
        expected['C'] = df['A']
        df = df_orig.copy()
        df.ix[:,'C'] = df.ix[:,'A']
        assert_frame_equal(df,expected)

        ### panel ###
        p_orig = Panel(np.arange(16).reshape(2,4,2),items=['Item1','Item2'],major_axis=pd.date_range('2001/1/12',periods=4),minor_axis=['A','B'],dtype='float64')

        # panel setting via item
        p_orig = Panel(np.arange(16).reshape(2,4,2),items=['Item1','Item2'],major_axis=pd.date_range('2001/1/12',periods=4),minor_axis=['A','B'],dtype='float64')
        expected = p_orig.copy()
        expected['Item3'] = expected['Item1']
        p = p_orig.copy()
        p.loc['Item3'] = p['Item1']
        assert_panel_equal(p,expected)

        # panel with aligned series
        expected = p_orig.copy()
        expected = expected.transpose(2,1,0)
        expected['C'] = DataFrame({ 'Item1' : [30,30,30,30], 'Item2' : [32,32,32,32] },index=p_orig.major_axis)
        expected = expected.transpose(2,1,0)
        p = p_orig.copy()
        p.loc[:,:,'C'] = Series([30,32],index=p_orig.items)
        assert_panel_equal(p,expected)

        # GH 8473
        dates = date_range('1/1/2000', periods=8)
        df_orig = DataFrame(np.random.randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])

        expected = pd.concat([df_orig,DataFrame({'A' : 7},index=[dates[-1]+1])])
        df = df_orig.copy()
        df.loc[dates[-1]+1, 'A'] = 7
        assert_frame_equal(df,expected)
        df = df_orig.copy()
        df.at[dates[-1]+1, 'A'] = 7
        assert_frame_equal(df,expected)

        expected = pd.concat([df_orig,DataFrame({0 : 7},index=[dates[-1]+1])],axis=1)
        df = df_orig.copy()
        df.loc[dates[-1]+1, 0] = 7
        assert_frame_equal(df,expected)
        df = df_orig.copy()
        df.at[dates[-1]+1, 0] = 7
        assert_frame_equal(df,expected)

    def test_partial_setting_mixed_dtype(self):

        # in a mixed dtype environment, try to preserve dtypes
        # by appending
        df = DataFrame([[True, 1],[False, 2]],
                       columns = ["female","fitness"])

        s = df.loc[1].copy()
        s.name = 2
        expected = df.append(s)

        df.loc[2] = df.loc[1]
        assert_frame_equal(df, expected)

        # columns will align
        df = DataFrame(columns=['A','B'])
        df.loc[0] = Series(1,index=range(4))
        assert_frame_equal(df,DataFrame(columns=['A','B'],index=[0]))

        # columns will align
        df = DataFrame(columns=['A','B'])
        df.loc[0] = Series(1,index=['B'])
        assert_frame_equal(df,DataFrame([[np.nan, 1]], columns=['A','B'],index=[0],dtype='float64'))

        # list-like must conform
        df = DataFrame(columns=['A','B'])
        def f():
            df.loc[0] = [1,2,3]
        self.assertRaises(ValueError, f)

        # these are coerced to float unavoidably (as its a list-like to begin)
        df = DataFrame(columns=['A','B'])
        df.loc[3] = [6,7]
        assert_frame_equal(df,DataFrame([[6,7]],index=[3],columns=['A','B'],dtype='float64'))

    def test_partial_setting_with_datetimelike_dtype(self):

        # GH9478
        # a datetimeindex alignment issue with partial setting
        df = pd.DataFrame(np.arange(6.).reshape(3,2), columns=list('AB'),
                          index=pd.date_range('1/1/2000', periods=3, freq='1H'))
        expected = df.copy()
        expected['C'] = [expected.index[0]] + [pd.NaT,pd.NaT]

        mask = df.A < 1
        df.loc[mask, 'C'] = df.loc[mask].index
        assert_frame_equal(df, expected)

    def test_loc_setitem_datetime(self):

        # GH 9516
        dt1 = Timestamp('20130101 09:00:00')
        dt2 = Timestamp('20130101 10:00:00')

        for conv in [lambda x: x, lambda x: x.to_datetime64(),
                     lambda x: x.to_pydatetime(), lambda x: np.datetime64(x)]:

            df = pd.DataFrame()
            df.loc[conv(dt1),'one'] = 100
            df.loc[conv(dt2),'one'] = 200

            expected = DataFrame({'one' : [100.0, 200.0]},index=[dt1, dt2])
            assert_frame_equal(df, expected)

    def test_series_partial_set(self):
        # partial set with new index
        # Regression from GH4825
        ser = Series([0.1, 0.2], index=[1, 2])

        # ToDo: check_index_type can be True after GH 11497

        # loc
        expected = Series([np.nan, 0.2, np.nan], index=[3, 2, 3])
        result = ser.loc[[3, 2, 3]]
        assert_series_equal(result, expected, check_index_type=False)

        # raises as nothing in in the index
        self.assertRaises(KeyError, lambda : ser.loc[[3, 3, 3]])

        expected = Series([0.2, 0.2, np.nan], index=[2, 2, 3])
        result = ser.loc[[2, 2, 3]]
        assert_series_equal(result, expected, check_index_type=False)

        expected = Series([0.3, np.nan, np.nan], index=[3, 4, 4])
        result = Series([0.1, 0.2, 0.3], index=[1, 2, 3]).loc[[3, 4, 4]]
        assert_series_equal(result, expected, check_index_type=False)

        expected = Series([np.nan, 0.3, 0.3], index=[5, 3, 3])
        result = Series([0.1, 0.2, 0.3, 0.4], index=[1, 2, 3, 4]).loc[[5, 3, 3]]
        assert_series_equal(result, expected, check_index_type=False)

        expected = Series([np.nan, 0.4, 0.4], index=[5, 4, 4])
        result = Series([0.1, 0.2, 0.3, 0.4], index=[1, 2, 3, 4]).loc[[5, 4, 4]]
        assert_series_equal(result, expected, check_index_type=False)

        expected = Series([0.4, np.nan, np.nan], index=[7, 2, 2])
        result = Series([0.1, 0.2, 0.3, 0.4], index=[4, 5, 6, 7]).loc[[7, 2, 2]]
        assert_series_equal(result, expected, check_index_type=False)

        expected = Series([0.4, np.nan, np.nan], index=[4, 5, 5])
        result = Series([0.1, 0.2, 0.3, 0.4], index=[1, 2, 3, 4]).loc[[4, 5, 5]]
        assert_series_equal(result, expected, check_index_type=False)

        # iloc
        expected = Series([0.2, 0.2, 0.1, 0.1], index=[2, 2, 1, 1])
        result = ser.iloc[[1, 1, 0, 0]]
        assert_series_equal(result, expected, check_index_type=False)

    def test_partial_set_invalid(self):

        # GH 4940
        # allow only setting of 'valid' values

        df = tm.makeTimeDataFrame()

        # don't allow not string inserts
        def f():
            df.loc[100.0, :] = df.ix[0]
        self.assertRaises(TypeError, f)
        def f():
            df.loc[100,:] = df.ix[0]
        self.assertRaises(TypeError, f)

        def f():
            df.ix[100.0, :] = df.ix[0]
        self.assertRaises(ValueError, f)
        def f():
            df.ix[100,:] = df.ix[0]
        self.assertRaises(ValueError, f)

        # allow object conversion here
        df.loc['a',:] = df.ix[0]

    def test_partial_set_empty(self):

        # GH5226

        # partially set with an empty object
        # series
        s = Series()
        s.loc[1] = 1
        assert_series_equal(s,Series([1],index=[1]))
        s.loc[3] = 3
        assert_series_equal(s,Series([1,3],index=[1,3]))

        s = Series()
        s.loc[1] = 1.
        assert_series_equal(s,Series([1.],index=[1]))
        s.loc[3] = 3.
        assert_series_equal(s,Series([1.,3.],index=[1,3]))

        s = Series()
        s.loc['foo'] = 1
        assert_series_equal(s,Series([1],index=['foo']))
        s.loc['bar'] = 3
        assert_series_equal(s,Series([1,3],index=['foo','bar']))
        s.loc[3] = 4
        assert_series_equal(s,Series([1,3,4],index=['foo','bar',3]))

        # partially set with an empty object
        # frame
        df = DataFrame()

        def f():
            df.loc[1] = 1
        self.assertRaises(ValueError, f)
        def f():
            df.loc[1] = Series([1],index=['foo'])
        self.assertRaises(ValueError, f)
        def f():
            df.loc[:,1] = 1
        self.assertRaises(ValueError, f)

        # these work as they don't really change
        # anything but the index
        # GH5632
        expected = DataFrame(columns=['foo'], index=pd.Index([], dtype='int64'))
        def f():
            df = DataFrame()
            df['foo'] = Series([], dtype='object')
            return df
        assert_frame_equal(f(), expected)
        def f():
            df = DataFrame()
            df['foo'] = Series(df.index)
            return df
        assert_frame_equal(f(), expected)
        def f():
            df = DataFrame()
            df['foo'] = df.index
            return df
        assert_frame_equal(f(), expected)

        expected = DataFrame(columns=['foo'], index=pd.Index([], dtype='int64'))
        expected['foo'] = expected['foo'].astype('float64')
        def f():
            df = DataFrame()
            df['foo'] = []
            return df
        assert_frame_equal(f(), expected)
        def f():
            df = DataFrame()
            df['foo'] = Series(range(len(df)))
            return df
        assert_frame_equal(f(), expected)
        def f():
            df = DataFrame()
            df['foo'] = range(len(df))
            return df
        assert_frame_equal(f(), expected)

        df = DataFrame()
        df2 = DataFrame()
        df2[1] = Series([1], index=['foo'])
        df.loc[:,1] = Series([1], index=['foo'])
        assert_frame_equal(df,DataFrame([[1]], index=['foo'], columns=[1]))
        assert_frame_equal(df,df2)

        # no index to start
        expected = DataFrame({ 0 : Series(1,index=range(4)) }, columns=['A','B',0])

        df = DataFrame(columns=['A','B'])
        df[0] = Series(1, index=range(4))
        df.dtypes
        str(df)
        assert_frame_equal(df,expected)

        df = DataFrame(columns=['A','B'])
        df.loc[:,0] = Series(1,index=range(4))
        df.dtypes
        str(df)
        assert_frame_equal(df,expected)

        # GH5720, GH5744
        # don't create rows when empty
        expected = DataFrame(columns=['A', 'B', 'New'], index=pd.Index([], dtype='int64'))
        expected['A'] = expected['A'].astype('int64')
        expected['B'] = expected['B'].astype('float64')
        expected['New'] = expected['New'].astype('float64')
        df = DataFrame({"A": [1, 2, 3], "B": [1.2, 4.2, 5.2]})
        y = df[df.A > 5]
        y['New'] = np.nan
        assert_frame_equal(y, expected)
        #assert_frame_equal(y,expected)

        expected = DataFrame(columns=['a', 'b', 'c c', 'd'])
        expected['d'] = expected['d'].astype('int64')
        df = DataFrame(columns=['a', 'b', 'c c'])
        df['d'] = 3
        assert_frame_equal(df, expected)
        assert_series_equal(df['c c'],Series(name='c c',dtype=object))

        # reindex columns is ok
        df = DataFrame({"A": [1, 2, 3], "B": [1.2, 4.2, 5.2]})
        y = df[df.A > 5]
        result = y.reindex(columns=['A','B','C'])
        expected = DataFrame(columns=['A','B','C'], index=pd.Index([], dtype='int64'))
        expected['A'] = expected['A'].astype('int64')
        expected['B'] = expected['B'].astype('float64')
        expected['C'] = expected['C'].astype('float64')
        assert_frame_equal(result,expected)

        # GH 5756
        # setting with empty Series
        df = DataFrame(Series())
        assert_frame_equal(df, DataFrame({ 0 : Series() }))

        df = DataFrame(Series(name='foo'))
        assert_frame_equal(df, DataFrame({ 'foo' : Series() }))

        # GH 5932
        # copy on empty with assignment fails
        df = DataFrame(index=[0])
        df = df.copy()
        df['a'] = 0
        expected = DataFrame(0,index=[0],columns=['a'])
        assert_frame_equal(df, expected)

        # GH 6171
        # consistency on empty frames
        df = DataFrame(columns=['x', 'y'])
        df['x'] = [1, 2]
        expected = DataFrame(dict(x = [1,2], y = [np.nan,np.nan]))
        assert_frame_equal(df, expected, check_dtype=False)

        df = DataFrame(columns=['x', 'y'])
        df['x'] = ['1', '2']
        expected = DataFrame(dict(x = ['1','2'], y = [np.nan,np.nan]),dtype=object)
        assert_frame_equal(df, expected)

        df = DataFrame(columns=['x', 'y'])
        df.loc[0, 'x'] = 1
        expected = DataFrame(dict(x = [1], y = [np.nan]))
        assert_frame_equal(df, expected, check_dtype=False)

    def test_cache_updating(self):
        # GH 4939, make sure to update the cache on setitem

        df = tm.makeDataFrame()
        df['A'] # cache series
        df.ix["Hello Friend"] = df.ix[0]
        self.assertIn("Hello Friend", df['A'].index)
        self.assertIn("Hello Friend", df['B'].index)

        panel = tm.makePanel()
        panel.ix[0] # get first item into cache
        panel.ix[:, :, 'A+1'] = panel.ix[:, :, 'A'] + 1
        self.assertIn("A+1", panel.ix[0].columns)
        self.assertIn("A+1", panel.ix[1].columns)

        # 5216
        # make sure that we don't try to set a dead cache
        a = np.random.rand(10, 3)
        df = DataFrame(a, columns=['x', 'y', 'z'])
        tuples = [(i, j) for i in range(5) for j in range(2)]
        index = MultiIndex.from_tuples(tuples)
        df.index = index

        # setting via chained assignment
        # but actually works, since everything is a view
        df.loc[0]['z'].iloc[0] = 1.
        result = df.loc[(0,0),'z']
        self.assertEqual(result, 1)

        # correct setting
        df.loc[(0,0),'z'] = 2
        result = df.loc[(0,0),'z']
        self.assertEqual(result, 2)

        # 10264
        df = DataFrame(np.zeros((5,5),dtype='int64'),columns=['a','b','c','d','e'],index=range(5))
        df['f'] = 0
        df.f.values[3] = 1
        y = df.iloc[np.arange(2,len(df))]
        df.f.values[3] = 2
        expected = DataFrame(np.zeros((5,6),dtype='int64'),columns=['a','b','c','d','e','f'],index=range(5))
        expected.at[3,'f'] = 2
        assert_frame_equal(df, expected)
        expected = Series([0,0,0,2,0],name='f')
        assert_series_equal(df.f, expected)

    def test_slice_consolidate_invalidate_item_cache(self):

        # this is chained assignment, but will 'work'
        with option_context('chained_assignment',None):

            # #3970
            df = DataFrame({ "aa":lrange(5), "bb":[2.2]*5})

            # Creates a second float block
            df["cc"] = 0.0

            # caches a reference to the 'bb' series
            df["bb"]

            # repr machinery triggers consolidation
            repr(df)

            # Assignment to wrong series
            df['bb'].iloc[0] = 0.17
            df._clear_item_cache()
            self.assertAlmostEqual(df['bb'][0], 0.17)

    def test_setitem_cache_updating(self):
        # GH 5424
        cont = ['one', 'two','three', 'four', 'five', 'six', 'seven']

        for do_ref in [False,False]:
            df = DataFrame({'a' : cont, "b":cont[3:]+cont[:3] ,'c' : np.arange(7)})

            # ref the cache
            if do_ref:
                df.ix[0,"c"]

            # set it
            df.ix[7,'c'] = 1

            self.assertEqual(df.ix[0,'c'], 0.0)
            self.assertEqual(df.ix[7,'c'], 1.0)

        # GH 7084
        # not updating cache on series setting with slices
        expected = DataFrame({'A': [600, 600, 600]}, index=date_range('5/7/2014', '5/9/2014'))
        out = DataFrame({'A': [0, 0, 0]}, index=date_range('5/7/2014', '5/9/2014'))
        df = DataFrame({'C': ['A', 'A', 'A'], 'D': [100, 200, 300]})

        #loop through df to update out
        six = Timestamp('5/7/2014')
        eix = Timestamp('5/9/2014')
        for ix, row in df.iterrows():
            out.loc[six:eix,row['C']] = out.loc[six:eix,row['C']] + row['D']

        assert_frame_equal(out, expected)
        assert_series_equal(out['A'], expected['A'])

        # try via a chain indexing
        # this actually works
        out = DataFrame({'A': [0, 0, 0]}, index=date_range('5/7/2014', '5/9/2014'))
        for ix, row in df.iterrows():
            v = out[row['C']][six:eix] + row['D']
            out[row['C']][six:eix] = v

        assert_frame_equal(out, expected)
        assert_series_equal(out['A'], expected['A'])

        out = DataFrame({'A': [0, 0, 0]}, index=date_range('5/7/2014', '5/9/2014'))
        for ix, row in df.iterrows():
            out.loc[six:eix,row['C']] += row['D']

        assert_frame_equal(out, expected)
        assert_series_equal(out['A'], expected['A'])

    def test_setitem_chained_setfault(self):

        # GH6026
        # setfaults under numpy 1.7.1 (ok on 1.8)
        data = ['right', 'left', 'left', 'left', 'right', 'left', 'timeout']
        mdata = ['right', 'left', 'left', 'left', 'right', 'left', 'none']

        df = DataFrame({'response': np.array(data)})
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        assert_frame_equal(df, DataFrame({'response': mdata }))

        recarray = np.rec.fromarrays([data], names=['response'])
        df = DataFrame(recarray)
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        assert_frame_equal(df, DataFrame({'response': mdata }))

        df = DataFrame({'response': data, 'response1' : data })
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        assert_frame_equal(df, DataFrame({'response': mdata, 'response1' : data }))

        # GH 6056
        expected = DataFrame(dict(A = [np.nan,'bar','bah','foo','bar']))
        df = DataFrame(dict(A = np.array(['foo','bar','bah','foo','bar'])))
        df['A'].iloc[0] = np.nan
        result = df.head()
        assert_frame_equal(result, expected)

        df = DataFrame(dict(A = np.array(['foo','bar','bah','foo','bar'])))
        df.A.iloc[0] = np.nan
        result = df.head()
        assert_frame_equal(result, expected)

    def test_detect_chained_assignment(self):

        pd.set_option('chained_assignment','raise')

        # work with the chain
        expected = DataFrame([[-5,1],[-6,3]],columns=list('AB'))
        df = DataFrame(np.arange(4).reshape(2,2),columns=list('AB'),dtype='int64')
        self.assertIsNone(df.is_copy)
        df['A'][0] = -5
        df['A'][1] = -6
        assert_frame_equal(df, expected)

        # test with the chaining
        df = DataFrame({ 'A' : Series(range(2),dtype='int64'), 'B' : np.array(np.arange(2,4),dtype=np.float64)})
        self.assertIsNone(df.is_copy)
        def f():
            df['A'][0] = -5
        self.assertRaises(com.SettingWithCopyError, f)
        def f():
            df['A'][1] = np.nan
        self.assertRaises(com.SettingWithCopyError, f)
        self.assertIsNone(df['A'].is_copy)

        # using a copy (the chain), fails
        df = DataFrame({ 'A' : Series(range(2),dtype='int64'), 'B' : np.array(np.arange(2,4),dtype=np.float64)})
        def f():
            df.loc[0]['A'] = -5
        self.assertRaises(com.SettingWithCopyError, f)

        # doc example
        df = DataFrame({'a' : ['one', 'one', 'two',
                               'three', 'two', 'one', 'six'],
                        'c' : Series(range(7),dtype='int64') })
        self.assertIsNone(df.is_copy)
        expected = DataFrame({'a' : ['one', 'one', 'two',
                                     'three', 'two', 'one', 'six'],
                              'c' : [42,42,2,3,4,42,6]})

        def f():
            indexer = df.a.str.startswith('o')
            df[indexer]['c'] = 42
        self.assertRaises(com.SettingWithCopyError, f)

        expected = DataFrame({'A':[111,'bbb','ccc'],'B':[1,2,3]})
        df = DataFrame({'A':['aaa','bbb','ccc'],'B':[1,2,3]})
        def f():
            df['A'][0] = 111
        self.assertRaises(com.SettingWithCopyError, f)
        def f():
            df.loc[0]['A'] = 111
        self.assertRaises(com.SettingWithCopyError, f)

        df.loc[0,'A'] = 111
        assert_frame_equal(df,expected)

        # make sure that is_copy is picked up reconstruction
        # GH5475
        df = DataFrame({"A": [1,2]})
        self.assertIsNone(df.is_copy)
        with tm.ensure_clean('__tmp__pickle') as path:
            df.to_pickle(path)
            df2 = pd.read_pickle(path)
            df2["B"] = df2["A"]
            df2["B"] = df2["A"]

        # a suprious raise as we are setting the entire column here
        # GH5597
        from string import ascii_letters as letters

        def random_text(nobs=100):
            df = []
            for i in range(nobs):
                idx= np.random.randint(len(letters), size=2)
                idx.sort()
                df.append([letters[idx[0]:idx[1]]])

            return DataFrame(df, columns=['letters'])

        df = random_text(100000)

        # always a copy
        x = df.iloc[[0,1,2]]
        self.assertIsNotNone(x.is_copy)
        x = df.iloc[[0,1,2,4]]
        self.assertIsNotNone(x.is_copy)

        # explicity copy
        indexer = df.letters.apply(lambda x : len(x) > 10)
        df = df.ix[indexer].copy()
        self.assertIsNone(df.is_copy)
        df['letters'] = df['letters'].apply(str.lower)

        # implicity take
        df = random_text(100000)
        indexer = df.letters.apply(lambda x : len(x) > 10)
        df = df.ix[indexer]
        self.assertIsNotNone(df.is_copy)
        df['letters'] = df['letters'].apply(str.lower)

        # implicity take 2
        df = random_text(100000)
        indexer = df.letters.apply(lambda x : len(x) > 10)
        df = df.ix[indexer]
        self.assertIsNotNone(df.is_copy)
        df.loc[:,'letters'] = df['letters'].apply(str.lower)

        # should be ok even though it's a copy!
        self.assertIsNone(df.is_copy)
        df['letters'] = df['letters'].apply(str.lower)
        self.assertIsNone(df.is_copy)

        df = random_text(100000)
        indexer = df.letters.apply(lambda x : len(x) > 10)
        df.ix[indexer,'letters'] = df.ix[indexer,'letters'].apply(str.lower)

        # an identical take, so no copy
        df = DataFrame({'a' : [1]}).dropna()
        self.assertIsNone(df.is_copy)
        df['a'] += 1

        # inplace ops
        # original from: http://stackoverflow.com/questions/20508968/series-fillna-in-a-multiindex-dataframe-does-not-fill-is-this-a-bug
        a = [12, 23]
        b = [123, None]
        c = [1234, 2345]
        d = [12345, 23456]
        tuples = [('eyes', 'left'), ('eyes', 'right'), ('ears', 'left'), ('ears', 'right')]
        events = {('eyes', 'left'): a, ('eyes', 'right'): b, ('ears', 'left'): c, ('ears', 'right'): d}
        multiind = MultiIndex.from_tuples(tuples, names=['part', 'side'])
        zed = DataFrame(events, index=['a', 'b'], columns=multiind)
        def f():
            zed['eyes']['right'].fillna(value=555, inplace=True)
        self.assertRaises(com.SettingWithCopyError, f)

        df = DataFrame(np.random.randn(10,4))
        s = df.iloc[:,0].sort_values()
        assert_series_equal(s,df.iloc[:,0].sort_values())
        assert_series_equal(s,df[0].sort_values())

        # false positives GH6025
        df = DataFrame ({'column1':['a', 'a', 'a'], 'column2': [4,8,9] })
        str(df)
        df['column1'] = df['column1'] + 'b'
        str(df)
        df = df [df['column2']!=8]
        str(df)
        df['column1'] = df['column1'] + 'c'
        str(df)

        # from SO: http://stackoverflow.com/questions/24054495/potential-bug-setting-value-for-undefined-column-using-iloc
        df = DataFrame(np.arange(0,9), columns=['count'])
        df['group'] = 'b'
        def f():
            df.iloc[0:5]['group'] = 'a'
        self.assertRaises(com.SettingWithCopyError, f)

        # mixed type setting
        # same dtype & changing dtype
        df = DataFrame(dict(A=date_range('20130101',periods=5),B=np.random.randn(5),C=np.arange(5,dtype='int64'),D=list('abcde')))

        def f():
            df.ix[2]['D'] = 'foo'
        self.assertRaises(com.SettingWithCopyError, f)
        def f():
            df.ix[2]['C'] = 'foo'
        self.assertRaises(com.SettingWithCopyError, f)
        def f():
            df['C'][2] = 'foo'
        self.assertRaises(com.SettingWithCopyError, f)

    def test_setting_with_copy_bug(self):

        # operating on a copy
        df = pd.DataFrame({'a': list(range(4)), 'b': list('ab..'), 'c': ['a', 'b', np.nan, 'd']})
        mask = pd.isnull(df.c)

        def f():
            df[['c']][mask] = df[['b']][mask]
        self.assertRaises(com.SettingWithCopyError, f)

        # invalid warning as we are returning a new object
        # GH 8730
        df1 = DataFrame({'x': Series(['a','b','c']), 'y': Series(['d','e','f'])})
        df2 = df1[['x']]

        # this should not raise
        df2['y'] = ['g', 'h', 'i']

    def test_detect_chained_assignment_warnings(self):

        # warnings
        with option_context('chained_assignment','warn'):
            df = DataFrame({'A':['aaa','bbb','ccc'],'B':[1,2,3]})
            with tm.assert_produces_warning(expected_warning=com.SettingWithCopyWarning):
                df.loc[0]['A'] = 111

    def test_float64index_slicing_bug(self):
        # GH 5557, related to slicing a float index
        ser = {256: 2321.0, 1: 78.0, 2: 2716.0, 3: 0.0, 4: 369.0, 5: 0.0, 6: 269.0, 7: 0.0, 8: 0.0, 9: 0.0, 10: 3536.0, 11: 0.0, 12: 24.0, 13: 0.0, 14: 931.0, 15: 0.0, 16: 101.0, 17: 78.0, 18: 9643.0, 19: 0.0, 20: 0.0, 21: 0.0, 22: 63761.0, 23: 0.0, 24: 446.0, 25: 0.0, 26: 34773.0, 27: 0.0, 28: 729.0, 29: 78.0, 30: 0.0, 31: 0.0, 32: 3374.0, 33: 0.0, 34: 1391.0, 35: 0.0, 36: 361.0, 37: 0.0, 38: 61808.0, 39: 0.0, 40: 0.0, 41: 0.0, 42: 6677.0, 43: 0.0, 44: 802.0, 45: 0.0, 46: 2691.0, 47: 0.0, 48: 3582.0, 49: 0.0, 50: 734.0, 51: 0.0, 52: 627.0, 53: 70.0, 54: 2584.0, 55: 0.0, 56: 324.0, 57: 0.0, 58: 605.0, 59: 0.0, 60: 0.0, 61: 0.0, 62: 3989.0, 63: 10.0, 64: 42.0, 65: 0.0, 66: 904.0, 67: 0.0, 68: 88.0, 69: 70.0, 70: 8172.0, 71: 0.0, 72: 0.0, 73: 0.0, 74: 64902.0, 75: 0.0, 76: 347.0, 77: 0.0, 78: 36605.0, 79: 0.0, 80: 379.0, 81: 70.0, 82: 0.0, 83: 0.0, 84: 3001.0, 85: 0.0, 86: 1630.0, 87: 7.0, 88: 364.0, 89: 0.0, 90: 67404.0, 91: 9.0, 92: 0.0, 93: 0.0, 94: 7685.0, 95: 0.0, 96: 1017.0, 97: 0.0, 98: 2831.0, 99: 0.0, 100: 2963.0, 101: 0.0, 102: 854.0, 103: 0.0, 104: 0.0, 105: 0.0, 106: 0.0, 107: 0.0, 108: 0.0, 109: 0.0, 110: 0.0, 111: 0.0, 112: 0.0, 113: 0.0, 114: 0.0, 115: 0.0, 116: 0.0, 117: 0.0, 118: 0.0, 119: 0.0, 120: 0.0, 121: 0.0, 122: 0.0, 123: 0.0, 124: 0.0, 125: 0.0, 126: 67744.0, 127: 22.0, 128: 264.0, 129: 0.0, 260: 197.0, 268: 0.0, 265: 0.0, 269: 0.0, 261: 0.0, 266: 1198.0, 267: 0.0, 262: 2629.0, 258: 775.0, 257: 0.0, 263: 0.0, 259: 0.0, 264: 163.0, 250: 10326.0, 251: 0.0, 252: 1228.0, 253: 0.0, 254: 2769.0, 255: 0.0}

        # smoke test for the repr
        s = Series(ser)
        result  = s.value_counts()
        str(result)

    def test_floating_index_doc_example(self):

        index = Index([1.5, 2, 3, 4.5, 5])
        s = Series(range(5),index=index)
        self.assertEqual(s[3], 2)
        self.assertEqual(s.ix[3], 2)
        self.assertEqual(s.loc[3], 2)
        self.assertEqual(s.iloc[3], 3)

    def test_floating_index(self):

        # related 236
        # scalar/slicing of a float index
        s = Series(np.arange(5), index=np.arange(5) * 2.5, dtype=np.int64)

        # label based slicing
        result1 = s[1.0:3.0]
        result2 = s.ix[1.0:3.0]
        result3 = s.loc[1.0:3.0]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)

        # exact indexing when found
        result1 = s[5.0]
        result2 = s.loc[5.0]
        result3 = s.ix[5.0]
        self.assertEqual(result1, result2)
        self.assertEqual(result1, result3)

        result1 = s[5]
        result2 = s.loc[5]
        result3 = s.ix[5]
        self.assertEqual(result1, result2)
        self.assertEqual(result1, result3)

        self.assertEqual(s[5.0], s[5])

        # value not found (and no fallbacking at all)

        # scalar integers
        self.assertRaises(KeyError, lambda : s.loc[4])
        self.assertRaises(KeyError, lambda : s.ix[4])
        self.assertRaises(KeyError, lambda : s[4])

        # fancy floats/integers create the correct entry (as nan)
        # fancy tests
        expected = Series([2, 0], index=Float64Index([5.0, 0.0]))
        for fancy_idx in [[5.0, 0.0], np.array([5.0, 0.0])]: # float
            assert_series_equal(s[fancy_idx], expected)
            assert_series_equal(s.loc[fancy_idx], expected)
            assert_series_equal(s.ix[fancy_idx], expected)

        expected = Series([2, 0], index=Index([5, 0], dtype='int64'))
        for fancy_idx in [[5, 0], np.array([5, 0])]: #int
            assert_series_equal(s[fancy_idx], expected)
            assert_series_equal(s.loc[fancy_idx], expected)
            assert_series_equal(s.ix[fancy_idx], expected)

        # all should return the same as we are slicing 'the same'
        result1 = s.loc[2:5]
        result2 = s.loc[2.0:5.0]
        result3 = s.loc[2.0:5]
        result4 = s.loc[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        # previously this did fallback indexing
        result1 = s[2:5]
        result2 = s[2.0:5.0]
        result3 = s[2.0:5]
        result4 = s[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        result1 = s.ix[2:5]
        result2 = s.ix[2.0:5.0]
        result3 = s.ix[2.0:5]
        result4 = s.ix[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        # combined test
        result1 = s.loc[2:5]
        result2 = s.ix[2:5]
        result3 = s[2:5]

        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)

        # list selection
        result1 = s[[0.0,5,10]]
        result2 = s.loc[[0.0,5,10]]
        result3 = s.ix[[0.0,5,10]]
        result4 = s.iloc[[0,2,4]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        result1 = s[[1.6,5,10]]
        result2 = s.loc[[1.6,5,10]]
        result3 = s.ix[[1.6,5,10]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series([np.nan,2,4],index=[1.6,5,10]))

        result1 = s[[0,1,2]]
        result2 = s.ix[[0,1,2]]
        result3 = s.loc[[0,1,2]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series([0.0,np.nan,np.nan],index=[0,1,2]))

        result1 = s.loc[[2.5, 5]]
        result2 = s.ix[[2.5, 5]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, Series([1,2],index=[2.5,5.0]))

        result1 = s[[2.5]]
        result2 = s.ix[[2.5]]
        result3 = s.loc[[2.5]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series([1],index=[2.5]))

    def test_scalar_indexer(self):
        # float indexing checked above

        def check_invalid(index, loc=None, iloc=None, ix=None, getitem=None):

            # related 236/4850
            # trying to access with a float index
            s = Series(np.arange(len(index)),index=index)

            if iloc is None:
                iloc = TypeError
            self.assertRaises(iloc, lambda : s.iloc[3.5])
            if loc is None:
                loc = TypeError
            self.assertRaises(loc, lambda : s.loc[3.5])
            if ix is None:
                ix = TypeError
            self.assertRaises(ix, lambda : s.ix[3.5])
            if getitem is None:
                getitem = TypeError
            self.assertRaises(getitem, lambda : s[3.5])

        for index in [ tm.makeStringIndex, tm.makeUnicodeIndex, tm.makeIntIndex,
                       tm.makeDateIndex, tm.makePeriodIndex ]:
            check_invalid(index())
        check_invalid(Index(np.arange(5) * 2.5),loc=KeyError, ix=KeyError, getitem=KeyError)

        def check_index(index, error):
            index = index()
            s = Series(np.arange(len(index)),index=index)

            # positional selection
            result1 = s[5]
            result2 = s[5.0]
            result3 = s.iloc[5]
            result4 = s.iloc[5.0]

            # by value
            self.assertRaises(error, lambda : s.loc[5])
            self.assertRaises(error, lambda : s.loc[5.0])

            # this is fallback, so it works
            result5 = s.ix[5]
            result6 = s.ix[5.0]

            self.assertEqual(result1, result2)
            self.assertEqual(result1, result3)
            self.assertEqual(result1, result4)
            self.assertEqual(result1, result5)
            self.assertEqual(result1, result6)

        # string-like
        for index in [ tm.makeStringIndex, tm.makeUnicodeIndex ]:
            check_index(index, KeyError)

        # datetimelike
        for index in [ tm.makeDateIndex, tm.makeTimedeltaIndex, tm.makePeriodIndex ]:
            check_index(index, TypeError)

        # exact indexing when found on IntIndex
        s = Series(np.arange(10),dtype='int64')

        result1 = s[5.0]
        result2 = s.loc[5.0]
        result3 = s.ix[5.0]
        result4 = s[5]
        result5 = s.loc[5]
        result6 = s.ix[5]
        self.assertEqual(result1, result2)
        self.assertEqual(result1, result3)
        self.assertEqual(result1, result4)
        self.assertEqual(result1, result5)
        self.assertEqual(result1, result6)

    def test_slice_indexer(self):

        def check_iloc_compat(s):
            # invalid type for iloc (but works with a warning)
            # check_stacklevel=False -> impossible to get it right for all
            # index types
            with self.assert_produces_warning(
                    FutureWarning, check_stacklevel=False):
                s.iloc[6.0:8]
            with self.assert_produces_warning(
                    FutureWarning, check_stacklevel=False):
                s.iloc[6.0:8.0]
            with self.assert_produces_warning(
                    FutureWarning, check_stacklevel=False):
                s.iloc[6:8.0]

        def check_slicing_positional(index):

            s = Series(np.arange(len(index))+10,index=index)

            # these are all positional
            result1 = s[2:5]
            result2 = s.ix[2:5]
            result3 = s.iloc[2:5]
            assert_series_equal(result1, result2)
            assert_series_equal(result1, result3)

            # loc will fail
            self.assertRaises(TypeError, lambda : s.loc[2:5])

            # make all float slicing fail
            self.assertRaises(TypeError, lambda : s[2.0:5])
            self.assertRaises(TypeError, lambda : s[2.0:5.0])
            self.assertRaises(TypeError, lambda : s[2:5.0])

            self.assertRaises(TypeError, lambda : s.ix[2.0:5])
            self.assertRaises(TypeError, lambda : s.ix[2.0:5.0])
            self.assertRaises(TypeError, lambda : s.ix[2:5.0])

            self.assertRaises(TypeError, lambda : s.loc[2.0:5])
            self.assertRaises(TypeError, lambda : s.loc[2.0:5.0])
            self.assertRaises(TypeError, lambda : s.loc[2:5.0])

            check_iloc_compat(s)

        # all index types except int, float
        for index in [ tm.makeStringIndex, tm.makeUnicodeIndex,
                       tm.makeDateIndex, tm.makeTimedeltaIndex, tm.makePeriodIndex ]:
            check_slicing_positional(index())

        ############
        # IntIndex #
        ############
        index = tm.makeIntIndex()
        s = Series(np.arange(len(index),dtype='int64')+10,index+5)

        # this is positional
        result1 = s[2:5]
        result4 = s.iloc[2:5]
        assert_series_equal(result1, result4)

        # these are all label based
        result2 = s.ix[2:5]
        result3 = s.loc[2:5]
        assert_series_equal(result2, result3)

        # float slicers on an int index
        expected = Series([11,12,13],index=[6,7,8])
        for method in [lambda x: x.loc, lambda x: x.ix]:
            result = method(s)[6.0:8.5]
            assert_series_equal(result, expected)

            result = method(s)[5.5:8.5]
            assert_series_equal(result, expected)

            result = method(s)[5.5:8.0]
            assert_series_equal(result, expected)

        # make all float slicing fail for [] with an int index
        self.assertRaises(TypeError, lambda : s[6.0:8])
        self.assertRaises(TypeError, lambda : s[6.0:8.0])
        self.assertRaises(TypeError, lambda : s[6:8.0])

        check_iloc_compat(s)

        ##############
        # FloatIndex #
        ##############
        s.index = s.index.astype('float64')

        # these are all value based
        result1 = s[6:8]
        result2 = s.ix[6:8]
        result3 = s.loc[6:8]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)

        # these are valid for all methods
        # these are treated like labels (e.g. the rhs IS included)
        def compare(slicers, expected):
            for method in [lambda x: x, lambda x: x.loc, lambda x: x.ix ]:
                for slices in slicers:

                    result = method(s)[slices]
                    assert_series_equal(result, expected)

        compare([slice(6.0,8),slice(6.0,8.0),slice(6,8.0)],
                s[(s.index>=6.0)&(s.index<=8)])
        compare([slice(6.5,8),slice(6.5,8.5)],
                s[(s.index>=6.5)&(s.index<=8.5)])
        compare([slice(6,8.5)],
                s[(s.index>=6.0)&(s.index<=8.5)])
        compare([slice(6.5,6.5)],
                s[(s.index>=6.5)&(s.index<=6.5)])

        check_iloc_compat(s)

    def test_set_ix_out_of_bounds_axis_0(self):
        df = pd.DataFrame(randn(2, 5), index=["row%s" % i for i in range(2)], columns=["col%s" % i for i in range(5)])
        self.assertRaises(ValueError, df.ix.__setitem__, (2, 0), 100)

    def test_set_ix_out_of_bounds_axis_1(self):
        df = pd.DataFrame(randn(5, 2), index=["row%s" % i for i in range(5)], columns=["col%s" % i for i in range(2)])
        self.assertRaises(ValueError, df.ix.__setitem__, (0 , 2), 100)

    def test_iloc_empty_list_indexer_is_ok(self):
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = mkdf(5, 2)
        # vertical empty
        assert_frame_equal(df.iloc[:, []], df.iloc[:, :0],
                           check_index_type=True, check_column_type=True)
        # horizontal empty
        assert_frame_equal(df.iloc[[], :], df.iloc[:0, :],
                           check_index_type=True, check_column_type=True)
        # horizontal empty
        assert_frame_equal(df.iloc[[]], df.iloc[:0, :],
                           check_index_type=True, check_column_type=True)

    def test_loc_empty_list_indexer_is_ok(self):
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = mkdf(5, 2)
        # vertical empty
        assert_frame_equal(df.loc[:, []], df.iloc[:, :0],
                           check_index_type=True, check_column_type=True)
        # horizontal empty
        assert_frame_equal(df.loc[[], :], df.iloc[:0, :],
                           check_index_type=True, check_column_type=True)
        # horizontal empty
        assert_frame_equal(df.loc[[]], df.iloc[:0, :],
                           check_index_type=True, check_column_type=True)

    def test_ix_empty_list_indexer_is_ok(self):
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = mkdf(5, 2)
        # vertical empty
        assert_frame_equal(df.ix[:, []], df.iloc[:, :0],
                           check_index_type=True, check_column_type=True)
        # horizontal empty
        assert_frame_equal(df.ix[[], :], df.iloc[:0, :],
                           check_index_type=True, check_column_type=True)
        # horizontal empty
        assert_frame_equal(df.ix[[]], df.iloc[:0, :],
                           check_index_type=True, check_column_type=True)

    def test_deprecate_float_indexers(self):

        # GH 4892
        # deprecate allowing float indexers that are equal to ints to be used
        # as indexers in non-float indices

        import warnings
        warnings.filterwarnings(action='error', category=FutureWarning)

        def check_index(index):
            i = index(5)

            for s in  [ Series(np.arange(len(i)),index=i), DataFrame(np.random.randn(len(i),len(i)),index=i,columns=i) ]:
                self.assertRaises(FutureWarning, lambda :
                                  s.iloc[3.0])

                # setting
                def f():
                    s.iloc[3.0] = 0
                self.assertRaises(FutureWarning, f)

            # fallsback to position selection ,series only
            s = Series(np.arange(len(i)),index=i)
            s[3]
            self.assertRaises(FutureWarning, lambda : s[3.0])

        for index in [ tm.makeStringIndex, tm.makeUnicodeIndex,
                       tm.makeDateIndex, tm.makeTimedeltaIndex, tm.makePeriodIndex ]:
            check_index(index)

        # ints
        i = index(5)
        for s in [ Series(np.arange(len(i))), DataFrame(np.random.randn(len(i),len(i)),index=i,columns=i) ]:
            self.assertRaises(FutureWarning, lambda :
                              s.iloc[3.0])

            # on some arch's this doesn't provide a warning (and thus raise)
            # and some it does
            try:
                s[3.0]
            except:
                pass

            # setting
            def f():
                s.iloc[3.0] = 0
            self.assertRaises(FutureWarning, f)

        # floats: these are all ok!
        i = np.arange(5.)

        for s in [ Series(np.arange(len(i)),index=i), DataFrame(np.random.randn(len(i),len(i)),index=i,columns=i) ]:
            with tm.assert_produces_warning(False):
                s[3.0]

            with tm.assert_produces_warning(False):
                s[3]

            self.assertRaises(FutureWarning, lambda :
                              s.iloc[3.0])

            with tm.assert_produces_warning(False):
                s.iloc[3]

            with tm.assert_produces_warning(False):
                s.loc[3.0]

            with tm.assert_produces_warning(False):
                s.loc[3]

            def f():
                s.iloc[3.0] = 0
            self.assertRaises(FutureWarning, f)

        # slices
        for index in [ tm.makeIntIndex, tm.makeFloatIndex,
                       tm.makeStringIndex, tm.makeUnicodeIndex,
                       tm.makeDateIndex, tm.makePeriodIndex ]:

            index = index(5)
            for s in [ Series(range(5),index=index), DataFrame(np.random.randn(5,2),index=index) ]:

                # getitem
                self.assertRaises(FutureWarning, lambda :
                                  s.iloc[3.0:4])
                self.assertRaises(FutureWarning, lambda :
                                  s.iloc[3.0:4.0])
                self.assertRaises(FutureWarning, lambda :
                                  s.iloc[3:4.0])

                # setitem
                def f():
                    s.iloc[3.0:4] = 0
                self.assertRaises(FutureWarning, f)
                def f():
                    s.iloc[3:4.0] = 0
                self.assertRaises(FutureWarning, f)
                def f():
                    s.iloc[3.0:4.0] = 0
                self.assertRaises(FutureWarning, f)

        warnings.filterwarnings(action='ignore', category=FutureWarning)

    def test_float_index_to_mixed(self):
        df = DataFrame({0.0: np.random.rand(10),
                        1.0: np.random.rand(10)})
        df['a'] = 10
        tm.assert_frame_equal(DataFrame({0.0: df[0.0],
                                         1.0: df[1.0],
                                         'a': [10] * 10}),
                              df)

    def test_duplicate_ix_returns_series(self):
        df = DataFrame(np.random.randn(3, 3), index=[0.1, 0.2, 0.2],
                       columns=list('abc'))
        r = df.ix[0.2, 'a']
        e = df.loc[0.2, 'a']
        tm.assert_series_equal(r, e)

    def test_float_index_non_scalar_assignment(self):
        df = DataFrame({'a': [1,2,3], 'b': [3,4,5]},index=[1.,2.,3.])
        df.loc[df.index[:2]] = 1
        expected = DataFrame({'a':[1,1,3],'b':[1,1,5]},index=df.index)
        tm.assert_frame_equal(expected, df)

        df = DataFrame({'a': [1,2,3], 'b': [3,4,5]},index=[1.,2.,3.])
        df2 = df.copy()
        df.loc[df.index] = df.loc[df.index]
        tm.assert_frame_equal(df,df2)

    def test_float_index_at_iat(self):
        s = pd.Series([1, 2, 3], index=[0.1, 0.2, 0.3])
        for el, item in s.iteritems():
            self.assertEqual(s.at[el], item)
        for i in range(len(s)):
            self.assertEqual(s.iat[i], i + 1)

    def test_rhs_alignment(self):
        # GH8258, tests that both rows & columns are aligned to what is
        # assigned to. covers both uniform data-type & multi-type cases
        def run_tests(df, rhs, right):
            # label, index, slice
            r, i, s = list('bcd'), [1, 2, 3], slice(1, 4)
            c, j, l = ['joe', 'jolie'], [1, 2], slice(1, 3)

            left = df.copy()
            left.loc[r, c] = rhs
            assert_frame_equal(left, right)

            left = df.copy()
            left.iloc[i, j] = rhs
            assert_frame_equal(left, right)

            left = df.copy()
            left.ix[s, l] = rhs
            assert_frame_equal(left, right)

            left = df.copy()
            left.ix[i, j] = rhs
            assert_frame_equal(left, right)

            left = df.copy()
            left.ix[r, c] = rhs
            assert_frame_equal(left, right)

        xs = np.arange(20).reshape(5, 4)
        cols = ['jim', 'joe', 'jolie', 'joline']
        df = pd.DataFrame(xs, columns=cols, index=list('abcde'))

        # right hand side; permute the indices and multiplpy by -2
        rhs = - 2 * df.iloc[3:0:-1, 2:0:-1]

        # expected `right` result; just multiply by -2
        right = df.copy()
        right.iloc[1:4, 1:3] *= -2

        # run tests with uniform dtypes
        run_tests(df, rhs, right)

        # make frames multi-type & re-run tests
        for frame in [df, rhs, right]:
            frame['joe'] = frame['joe'].astype('float64')
            frame['jolie'] = frame['jolie'].map('@{0}'.format)

        run_tests(df, rhs, right)

    def test_str_label_slicing_with_negative_step(self):
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            assert_series_equal(s.loc[l_slc], s.iloc[i_slc])

            if not idx.is_integer:
                # For integer indices, ix and plain getitem are position-based.
                assert_series_equal(s[l_slc], s.iloc[i_slc])
                assert_series_equal(s.ix[l_slc], s.iloc[i_slc])

        for idx in [_mklbl('A', 20), np.arange(20) + 100,
                    np.linspace(100, 150, 20)]:
            idx = Index(idx)
            s = Series(np.arange(20), index=idx)
            assert_slices_equivalent(SLC[idx[9]::-1], SLC[9::-1])
            assert_slices_equivalent(SLC[:idx[9]:-1], SLC[:8:-1])
            assert_slices_equivalent(SLC[idx[13]:idx[9]:-1], SLC[13:8:-1])
            assert_slices_equivalent(SLC[idx[9]:idx[13]:-1], SLC[:0])

    def test_multiindex_label_slicing_with_negative_step(self):
        s = Series(np.arange(20),
                   MultiIndex.from_product([list('abcde'), np.arange(4)]))
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            assert_series_equal(s.loc[l_slc], s.iloc[i_slc])
            assert_series_equal(s[l_slc], s.iloc[i_slc])
            assert_series_equal(s.ix[l_slc], s.iloc[i_slc])

        assert_slices_equivalent(SLC[::-1], SLC[::-1])

        assert_slices_equivalent(SLC['d'::-1], SLC[15::-1])
        assert_slices_equivalent(SLC[('d',)::-1], SLC[15::-1])

        assert_slices_equivalent(SLC[:'d':-1], SLC[:11:-1])
        assert_slices_equivalent(SLC[:('d',):-1], SLC[:11:-1])

        assert_slices_equivalent(SLC['d':'b':-1], SLC[15:3:-1])
        assert_slices_equivalent(SLC[('d',):'b':-1], SLC[15:3:-1])
        assert_slices_equivalent(SLC['d':('b',):-1], SLC[15:3:-1])
        assert_slices_equivalent(SLC[('d',):('b',):-1], SLC[15:3:-1])
        assert_slices_equivalent(SLC['b':'d':-1], SLC[:0])

        assert_slices_equivalent(SLC[('c', 2)::-1], SLC[10::-1])
        assert_slices_equivalent(SLC[:('c', 2):-1], SLC[:9:-1])
        assert_slices_equivalent(SLC[('e', 0):('c', 2):-1], SLC[16:9:-1])

    def test_slice_with_zero_step_raises(self):
        s = Series(np.arange(20), index=_mklbl('A', 20))
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: s[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: s.loc[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: s.ix[::0])

    def test_indexing_assignment_dict_already_exists(self):
        df = pd.DataFrame({'x': [1, 2, 6],
                           'y': [2, 2, 8],
                           'z': [-5, 0, 5]}).set_index('z')
        expected = df.copy()
        rhs = dict(x=9, y=99)
        df.loc[5] = rhs
        expected.loc[5] = [9, 99]
        tm.assert_frame_equal(df, expected)

    def test_indexing_dtypes_on_empty(self):
        # Check that .iloc and .ix return correct dtypes GH9983
        df = DataFrame({'a':[1,2,3],'b':['b','b2','b3']})
        df2 = df.ix[[],:]

        self.assertEqual(df2.loc[:,'a'].dtype, np.int64)
        assert_series_equal(df2.loc[:,'a'], df2.iloc[:,0])
        assert_series_equal(df2.loc[:,'a'], df2.ix[:,0])

    def test_range_in_series_indexing(self):
        # range can cause an indexing error
        # GH 11652
        for x in [5, 999999, 1000000]:
            s = pd.Series(index=range(x))
            s.loc[range(1)] = 42
            assert_series_equal(s.loc[range(1)],Series(42.0,index=[0]))

            s.loc[range(2)] = 43
            assert_series_equal(s.loc[range(2)],Series(43.0,index=[0,1]))

    @slow
    def test_large_dataframe_indexing(self):
        #GH10692
        result = DataFrame({'x': range(10**6)},dtype='int64')
        result.loc[len(result)] = len(result) + 1
        expected = DataFrame({'x': range(10**6 + 1)},dtype='int64')
        assert_frame_equal(result, expected)

    @slow
    def test_large_mi_dataframe_indexing(self):
        #GH10645
        result = MultiIndex.from_arrays([range(10**6), range(10**6)])
        assert(not (10**6, 0) in result)

    def test_non_reducing_slice(self):
        df = pd.DataFrame([[0, 1], [2, 3]])

        slices = [
            # pd.IndexSlice[:, :],
            pd.IndexSlice[:, 1],
            pd.IndexSlice[1, :],
            pd.IndexSlice[[1], [1]],
            pd.IndexSlice[1, [1]],
            pd.IndexSlice[[1], 1],
            pd.IndexSlice[1],
            pd.IndexSlice[1, 1],
            slice(None, None, None),
            [0, 1],
            np.array([0, 1]),
            pd.Series([0, 1])
        ]
        for slice_ in slices:
            tslice_ = _non_reducing_slice(slice_)
            self.assertTrue(isinstance(df.loc[tslice_], DataFrame))

    def test_list_slice(self):
        # like dataframe getitem
        slices = [['A'], pd.Series(['A']), np.array(['A'])]
        df = pd.DataFrame({'A': [1, 2], 'B': [3, 4]}, index=['A', 'B'])
        expected = pd.IndexSlice[:, ['A']]
        for subset in slices:
            result = _non_reducing_slice(subset)
            tm.assert_frame_equal(df.loc[result], df.loc[expected])

    def test_maybe_numeric_slice(self):
        df = pd.DataFrame({'A': [1, 2], 'B': ['c', 'd'], 'C': [True, False]})
        result = _maybe_numeric_slice(df, slice_=None)
        expected = pd.IndexSlice[:, ['A']]
        self.assertEqual(result, expected)

        result = _maybe_numeric_slice(df, None, include_bool=True)
        expected = pd.IndexSlice[:, ['A', 'C']]
        result = _maybe_numeric_slice(df, [1])
        expected = [1]
        self.assertEqual(result, expected)


class TestCategoricalIndex(tm.TestCase):

    def setUp(self):

        self.df = DataFrame({'A' : np.arange(6,dtype='int64'),
                             'B' : Series(list('aabbca')).astype('category',categories=list('cab')) }).set_index('B')
        self.df2 = DataFrame({'A' : np.arange(6,dtype='int64'),
                              'B' : Series(list('aabbca')).astype('category',categories=list('cabe')) }).set_index('B')
        self.df3 = DataFrame({'A' : np.arange(6,dtype='int64'),
                              'B' : Series([1,1,2,1,3,2]).astype('category',categories=[3,2,1],ordered=True) }).set_index('B')
        self.df4 = DataFrame({'A' : np.arange(6,dtype='int64'),
                              'B' : Series([1,1,2,1,3,2]).astype('category',categories=[3,2,1],ordered=False) }).set_index('B')


    def test_loc_scalar(self):

        result = self.df.loc['a']
        expected = DataFrame({'A' : [0,1,5],
                              'B' : Series(list('aaa')).astype('category',categories=list('cab')) }).set_index('B')
        assert_frame_equal(result, expected)


        df = self.df.copy()
        df.loc['a'] = 20
        expected = DataFrame({'A' : [20,20,2,3,4,20],
                              'B' : Series(list('aabbca')).astype('category',categories=list('cab')) }).set_index('B')
        assert_frame_equal(df, expected)

        # value not in the categories
        self.assertRaises(KeyError, lambda : df.loc['d'])

        def f():
            df.loc['d'] = 10
        self.assertRaises(TypeError, f)

        def f():
            df.loc['d','A'] = 10
        self.assertRaises(TypeError, f)

        def f():
            df.loc['d','C'] = 10
        self.assertRaises(TypeError, f)

    def test_loc_listlike(self):

        # list of labels
        result = self.df.loc[['c','a']]
        expected = self.df.iloc[[4,0,1,5]]
        assert_frame_equal(result, expected)

        # ToDo: check_index_type can be True after GH XXX

        result = self.df2.loc[['a','b','e']]
        exp_index = pd.CategoricalIndex(list('aaabbe'), categories=list('cabe'), name='B')
        expected = DataFrame({'A' : [0,1,5,2,3,np.nan]}, index=exp_index)
        assert_frame_equal(result, expected, check_index_type=False)

        # element in the categories but not in the values
        self.assertRaises(KeyError, lambda : self.df2.loc['e'])

        # assign is ok
        df = self.df2.copy()
        df.loc['e'] = 20
        result = df.loc[['a','b','e']]
        exp_index = pd.CategoricalIndex(list('aaabbe'), categories=list('cabe'), name='B')
        expected = DataFrame({'A' : [0,1,5,2,3,20]}, index=exp_index)
        assert_frame_equal(result, expected)

        df = self.df2.copy()
        result = df.loc[['a','b','e']]
        expected = DataFrame({'A' : [0,1,5,2,3,np.nan],
                              'B' : Series(list('aaabbe')).astype('category',categories=list('cabe')) }).set_index('B')
        assert_frame_equal(result, expected, check_index_type=False)


        # not all labels in the categories
        self.assertRaises(KeyError, lambda : self.df2.loc[['a','d']])

    def test_read_only_source(self):
        # GH 10043
        rw_array = np.eye(10)
        rw_df = DataFrame(rw_array)

        ro_array = np.eye(10)
        ro_array.setflags(write=False)
        ro_df = DataFrame(ro_array)

        assert_frame_equal(rw_df.iloc[[1,2,3]],ro_df.iloc[[1,2,3]])
        assert_frame_equal(rw_df.iloc[[1]],ro_df.iloc[[1]])
        assert_series_equal(rw_df.iloc[1],ro_df.iloc[1])
        assert_frame_equal(rw_df.iloc[1:3],ro_df.iloc[1:3])

        assert_frame_equal(rw_df.loc[[1,2,3]],ro_df.loc[[1,2,3]])
        assert_frame_equal(rw_df.loc[[1]],ro_df.loc[[1]])
        assert_series_equal(rw_df.loc[1],ro_df.loc[1])
        assert_frame_equal(rw_df.loc[1:3],ro_df.loc[1:3])

    def test_reindexing(self):

        # reindexing
        # convert to a regular index
        result = self.df2.reindex(['a','b','e'])
        expected = DataFrame({'A' : [0,1,5,2,3,np.nan],
                              'B' : Series(list('aaabbe')) }).set_index('B')
        assert_frame_equal(result, expected)

        result = self.df2.reindex(['a','b'])
        expected = DataFrame({'A' : [0,1,5,2,3],
                              'B' : Series(list('aaabb')) }).set_index('B')
        assert_frame_equal(result, expected)

        result = self.df2.reindex(['e'])
        expected = DataFrame({'A' : [np.nan],
                              'B' : Series(['e']) }).set_index('B')
        assert_frame_equal(result, expected)

        result = self.df2.reindex(['d'])
        expected = DataFrame({'A' : [np.nan],
                              'B' : Series(['d']) }).set_index('B')
        assert_frame_equal(result, expected)

        # since we are actually reindexing with a Categorical
        # then return a Categorical
        cats = list('cabe')

        result = self.df2.reindex(pd.Categorical(['a','d'],categories=cats))
        expected = DataFrame({'A' : [0,1,5,np.nan],
                              'B' : Series(list('aaad')).astype('category',categories=cats) }).set_index('B')
        assert_frame_equal(result, expected)

        result = self.df2.reindex(pd.Categorical(['a'],categories=cats))
        expected = DataFrame({'A' : [0,1,5],
                              'B' : Series(list('aaa')).astype('category',categories=cats) }).set_index('B')
        assert_frame_equal(result, expected)

        result = self.df2.reindex(['a','b','e'])
        expected = DataFrame({'A' : [0,1,5,2,3,np.nan],
                              'B' : Series(list('aaabbe')) }).set_index('B')
        assert_frame_equal(result, expected)

        result = self.df2.reindex(['a','b'])
        expected = DataFrame({'A' : [0,1,5,2,3],
                              'B' : Series(list('aaabb')) }).set_index('B')
        assert_frame_equal(result, expected)

        result = self.df2.reindex(['e'])
        expected = DataFrame({'A' : [np.nan],
                              'B' : Series(['e']) }).set_index('B')
        assert_frame_equal(result, expected)

        # give back the type of categorical that we received
        result = self.df2.reindex(pd.Categorical(['a','d'],categories=cats,ordered=True))
        expected = DataFrame({'A' : [0,1,5,np.nan],
                              'B' : Series(list('aaad')).astype('category',categories=cats,ordered=True) }).set_index('B')
        assert_frame_equal(result, expected)

        result = self.df2.reindex(pd.Categorical(['a','d'],categories=['a','d']))
        expected = DataFrame({'A' : [0,1,5,np.nan],
                              'B' : Series(list('aaad')).astype('category',categories=['a','d']) }).set_index('B')
        assert_frame_equal(result, expected)

        # passed duplicate indexers are not allowed
        self.assertRaises(ValueError, lambda : self.df2.reindex(['a','a']))

        # args NotImplemented ATM
        self.assertRaises(NotImplementedError, lambda : self.df2.reindex(['a'],method='ffill'))
        self.assertRaises(NotImplementedError, lambda : self.df2.reindex(['a'],level=1))
        self.assertRaises(NotImplementedError, lambda : self.df2.reindex(['a'],limit=2))

    def test_loc_slice(self):

        # slicing
        # not implemented ATM
        # GH9748

        self.assertRaises(TypeError, lambda : self.df.loc[1:5])

        #result = df.loc[1:5]
        #expected = df.iloc[[1,2,3,4]]
        #assert_frame_equal(result, expected)

    def test_boolean_selection(self):

        df3 = self.df3
        df4 = self.df4

        result = df3[df3.index == 'a']
        expected = df3.iloc[[]]
        assert_frame_equal(result,expected)

        result = df4[df4.index == 'a']
        expected = df4.iloc[[]]
        assert_frame_equal(result,expected)

        result = df3[df3.index == 1]
        expected = df3.iloc[[0,1,3]]
        assert_frame_equal(result,expected)

        result = df4[df4.index == 1]
        expected = df4.iloc[[0,1,3]]
        assert_frame_equal(result,expected)

        # since we have an ordered categorical

        # CategoricalIndex([1, 1, 2, 1, 3, 2],
        #         categories=[3, 2, 1],
        #         ordered=True,
        #         name=u'B')
        result = df3[df3.index < 2]
        expected = df3.iloc[[4]]
        assert_frame_equal(result,expected)

        result = df3[df3.index > 1]
        expected = df3.iloc[[]]
        assert_frame_equal(result,expected)

        # unordered
        # cannot be compared

        # CategoricalIndex([1, 1, 2, 1, 3, 2],
        #         categories=[3, 2, 1],
        #         ordered=False,
        #         name=u'B')
        self.assertRaises(TypeError, lambda : df4[df4.index < 2])
        self.assertRaises(TypeError, lambda : df4[df4.index > 1])

class TestSeriesNoneCoercion(tm.TestCase):
    EXPECTED_RESULTS = [
        # For numeric series, we should coerce to NaN.
        ([1, 2, 3], [np.nan, 2, 3]),
        ([1.0, 2.0, 3.0], [np.nan, 2.0, 3.0]),

        # For datetime series, we should coerce to NaT.
        ([datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)],
         [NaT, datetime(2000, 1, 2), datetime(2000, 1, 3)]),

        # For objects, we should preserve the None value.
        (["foo", "bar", "baz"], [None, "bar", "baz"]),
    ]

    def test_coercion_with_setitem(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series[0] = None

            expected_series = Series(expected_result)

            assert_attr_equal('dtype', start_series, expected_series)
            tm.assert_numpy_array_equal(
                start_series.values,
                expected_series.values, strict_nan=True)

    def test_coercion_with_loc_setitem(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series.loc[0] = None

            expected_series = Series(expected_result)

            assert_attr_equal('dtype', start_series, expected_series)
            tm.assert_numpy_array_equal(
                start_series.values,
                expected_series.values, strict_nan=True)

    def test_coercion_with_setitem_and_series(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series[start_series == start_series[0]] = None

            expected_series = Series(expected_result)

            assert_attr_equal('dtype', start_series, expected_series)
            tm.assert_numpy_array_equal(
                start_series.values,
                expected_series.values, strict_nan=True)

    def test_coercion_with_loc_and_series(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series.loc[start_series == start_series[0]] = None

            expected_series = Series(expected_result)

            assert_attr_equal('dtype', start_series, expected_series)
            tm.assert_numpy_array_equal(
                start_series.values,
                expected_series.values, strict_nan=True)


class TestDataframeNoneCoercion(tm.TestCase):
    EXPECTED_SINGLE_ROW_RESULTS = [
        # For numeric series, we should coerce to NaN.
        ([1, 2, 3], [np.nan, 2, 3]),
        ([1.0, 2.0, 3.0], [np.nan, 2.0, 3.0]),

        # For datetime series, we should coerce to NaT.
        ([datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)],
         [NaT, datetime(2000, 1, 2), datetime(2000, 1, 3)]),

        # For objects, we should preserve the None value.
        (["foo", "bar", "baz"], [None, "bar", "baz"]),
    ]

    def test_coercion_with_loc(self):
        for start_data, expected_result, in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({'foo': start_data})
            start_dataframe.loc[0, ['foo']] = None

            expected_dataframe = DataFrame({'foo': expected_result})

            assert_attr_equal('dtype', start_dataframe['foo'], expected_dataframe['foo'])
            tm.assert_numpy_array_equal(
                start_dataframe['foo'].values,
                expected_dataframe['foo'].values, strict_nan=True)

    def test_coercion_with_setitem_and_dataframe(self):
        for start_data, expected_result, in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({'foo': start_data})
            start_dataframe[start_dataframe['foo'] == start_dataframe['foo'][0]] = None

            expected_dataframe = DataFrame({'foo': expected_result})

            assert_attr_equal('dtype', start_dataframe['foo'], expected_dataframe['foo'])
            tm.assert_numpy_array_equal(
                start_dataframe['foo'].values,
                expected_dataframe['foo'].values, strict_nan=True)

    def test_none_coercion_loc_and_dataframe(self):
        for start_data, expected_result, in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({'foo': start_data})
            start_dataframe.loc[start_dataframe['foo'] == start_dataframe['foo'][0]] = None

            expected_dataframe = DataFrame({'foo': expected_result})

            assert_attr_equal('dtype', start_dataframe['foo'], expected_dataframe['foo'])
            tm.assert_numpy_array_equal(
                start_dataframe['foo'].values,
                expected_dataframe['foo'].values, strict_nan=True)

    def test_none_coercion_mixed_dtypes(self):
        start_dataframe = DataFrame({
            'a': [1, 2, 3],
            'b': [1.0, 2.0, 3.0],
            'c': [datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)],
            'd': ['a', 'b', 'c']})
        start_dataframe.iloc[0] = None

        expected_dataframe = DataFrame({
            'a': [np.nan, 2, 3],
            'b': [np.nan, 2.0, 3.0],
            'c': [NaT, datetime(2000, 1, 2), datetime(2000, 1, 3)],
            'd': [None, 'b', 'c']})

        for column in expected_dataframe.columns:
            assert_attr_equal('dtype', start_dataframe[column], expected_dataframe[column])
            tm.assert_numpy_array_equal(
                start_dataframe[column].values,
                expected_dataframe[column].values, strict_nan=True)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
