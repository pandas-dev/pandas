# pylint: disable-msg=W0612,E1101
import unittest
import nose
import itertools
from StringIO import StringIO

from numpy import random, nan
from numpy.random import randn
import numpy as np
from numpy.testing import assert_array_equal

import pandas as pd
import pandas.core.common as com
from pandas.core.api import (DataFrame, Index, Series, Panel, notnull, isnull,
                             MultiIndex, DatetimeIndex, Timestamp)
from pandas.util.testing import (assert_almost_equal, assert_series_equal,
                                 assert_frame_equal, assert_panel_equal)
from pandas.util import py3compat

import pandas.util.testing as tm
import pandas.lib as lib
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
        axes = [ range(len(a)) for a in axes ]

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


class TestIndexing(unittest.TestCase):

    _multiprocess_can_split_ = True

    _objs = set(['series','frame','panel'])
    _typs = set(['ints','labels','mixed','ts','floats','empty'])

    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)

        self.series_ints   = Series(np.random.rand(4), index=range(0,8,2))
        self.frame_ints    = DataFrame(np.random.randn(4, 4), index=range(0, 8, 2), columns=range(0,12,3))
        self.panel_ints    = Panel(np.random.rand(4,4,4), items=range(0,8,2),major_axis=range(0,12,3),minor_axis=range(0,16,4))

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
                    print(v)

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
                        self.assert_(rs == xp)
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

            except (AssertionError):
                raise
            except (TypeError):
                raise AssertionError(_print('type error'))
            except (Exception), detail:

                # if we are in fails, the ok, otherwise raise it
                if fails is not None:
                    if fails == type(detail):
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

    def test_at_timestamp(self):

        # as timestamp is not a tuple!
        dates = date_range('1/1/2000', periods=8)
        df = DataFrame(randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])
        s = df['A']

        result = s.at[dates[5]]
        xp     = s.values[5]
        self.assert_(result == xp)

    def test_iat_invalid_args(self):
        pass

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
        self.check_result('list int', 'iloc', [0,1,2], 'indexer', [0,1,2], typs = ['labels','mixed','ts','floats','empty'], fails = IndexError)

    def test_iloc_getitem_dups(self):

        # no dups in panel (bug?)
        self.check_result('list int (dups)', 'iloc', [0,1,1,3], 'ix', { 0 : [0,2,2,6], 1 : [0,3,3,9] }, objs = ['series','frame'], typs = ['ints'])

    def test_iloc_getitem_array(self):

        # array like
        s = Series(index=range(1,4))
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

    def test_iloc_getitem_out_of_bounds(self):

        # out-of-bounds slice
        self.assertRaises(IndexError, self.frame_ints.iloc.__getitem__, tuple([slice(None),slice(1,5,None)]))
        self.assertRaises(IndexError, self.frame_ints.iloc.__getitem__, tuple([slice(None),slice(-5,3,None)]))
        self.assertRaises(IndexError, self.frame_ints.iloc.__getitem__, tuple([slice(1,5,None)]))
        self.assertRaises(IndexError, self.frame_ints.iloc.__getitem__, tuple([slice(-5,3,None)]))

    def test_iloc_setitem(self):
        df = self.frame_ints

        df.iloc[1,1] = 1
        result = df.iloc[1,1]
        self.assert_(result == 1)

        df.iloc[:,2:3] = 0
        expected = df.iloc[:,2:3]
        result = df.iloc[:,2:3]
        assert_frame_equal(result, expected)

    def test_iloc_multiindex(self):
        df = DataFrame(np.random.randn(3, 3),
                       columns=[[2,2,4],[6,8,10]],
                       index=[[4,4,8],[8,10,12]])

        rs = df.iloc[2]
        xp = df.irow(2)
        assert_series_equal(rs, xp)

        rs = df.iloc[:,2]
        xp = df.icol(2)
        assert_series_equal(rs, xp)

        rs = df.iloc[2,2]
        xp = df.values[2,2]
        self.assert_(rs == xp)

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
        self.check_result('label range', 'loc', 'f', 'ix', 'f', typs = ['ints','labels','mixed','ts','floats'], fails=KeyError)

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

        # fails
        self.check_result('list lbl', 'loc', [0,1,2], 'indexer', [0,1,2], typs = ['empty'], fails = KeyError)
        self.check_result('list lbl', 'loc', [0,2,3], 'ix', [0,2,3], typs = ['ints'], axes=0, fails = KeyError)
        self.check_result('list lbl', 'loc', [3,6,7], 'ix', [3,6,9], typs = ['ints'], axes=1, fails = KeyError)
        self.check_result('list lbl', 'loc', [4,8,10], 'ix', [4,8,12], typs = ['ints'], axes=2, fails = KeyError)

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

        # int slices in int
        self.check_result('int slice1', 'loc', slice(2,4), 'ix', { 0 : [2,4], 1: [3,6], 2: [4,8] }, typs = ['ints'], fails=KeyError)

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

    def test_loc_getitem_label_slice(self):

        # label slices (with ints)
        self.check_result('lab slice', 'loc', slice(1,3), 'ix', slice(1,3), typs = ['labels','mixed','ts','floats','empty'], fails=KeyError)

        # real label slices
        self.check_result('lab slice', 'loc', slice('a','c'), 'ix', slice('a','c'), typs = ['labels'], axes=0)
        self.check_result('lab slice', 'loc', slice('A','C'), 'ix', slice('A','C'), typs = ['labels'], axes=1)
        self.check_result('lab slice', 'loc', slice('W','Z'), 'ix', slice('W','Z'), typs = ['labels'], axes=2)

        self.check_result('ts  slice', 'loc', slice('20130102','20130104'), 'ix', slice('20130102','20130104'), typs = ['ts'], axes=0)
        self.check_result('ts  slice', 'loc', slice('20130102','20130104'), 'ix', slice('20130102','20130104'), typs = ['ts'], axes=1, fails=KeyError)
        self.check_result('ts  slice', 'loc', slice('20130102','20130104'), 'ix', slice('20130102','20130104'), typs = ['ts'], axes=2, fails=KeyError)

        self.check_result('mixed slice', 'loc', slice(2,8), 'ix', slice(2,8), typs = ['mixed'], axes=0, fails=KeyError)
        self.check_result('mixed slice', 'loc', slice(2,8), 'ix', slice(2,8), typs = ['mixed'], axes=1, fails=KeyError)
        self.check_result('mixed slice', 'loc', slice(2,8), 'ix', slice(2,8), typs = ['mixed'], axes=2, fails=KeyError)

        self.check_result('mixed slice', 'loc', slice(2,4,2), 'ix', slice(2,4,2), typs = ['mixed'], axes=0)

    def test_loc_general(self):

        # GH 2922 (these are fails)
        df = DataFrame(np.random.rand(4,4),columns=['A','B','C','D'])
        self.assertRaises(KeyError, df.loc.__getitem__, tuple([slice(0,2),slice(0,2)]))

        df = DataFrame(np.random.rand(4,4),columns=['A','B','C','D'], index=['A','B','C','D'])
        self.assertRaises(KeyError, df.loc.__getitem__, tuple([slice(0,2),df.columns[0:2]]))

        # want this to work
        result = df.loc[:,"A":"B"].iloc[0:2,:]
        self.assert_((result.columns == ['A','B']).all() == True)
        self.assert_((result.index == ['A','B']).all() == True)

    def test_loc_setitem_frame(self):
        df = self.frame_labels

        result = df.iloc[0,0]

        df.loc['a','A'] = 1
        result = df.loc['a','A']
        self.assert_(result == 1)

        result = df.iloc[0,0]
        self.assert_(result == 1)

        df.loc[:,'B':'D'] = 0
        expected = df.loc[:,'B':'D']
        result = df.ix[:,1:]
        assert_frame_equal(result, expected)

    def test_iloc_getitem_frame(self):
        """ originally from test_frame.py"""
        df = DataFrame(np.random.randn(10, 4), index=range(0, 20, 2), columns=range(0,8,2))

        result = df.iloc[2]
        exp = df.ix[4]
        assert_series_equal(result, exp)

        result = df.iloc[2,2]
        exp = df.ix[4,4]
        self.assert_(result == exp)

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
        s = Series(index=range(1,5))
        result = df.iloc[s.index]
        expected = df.ix[[2,4,6,8]]
        assert_frame_equal(result, expected)

        # out-of-bounds slice
        self.assertRaises(IndexError, df.iloc.__getitem__, tuple([slice(None),slice(1,5,None)]))
        self.assertRaises(IndexError, df.iloc.__getitem__, tuple([slice(None),slice(-5,3,None)]))
        self.assertRaises(IndexError, df.iloc.__getitem__, tuple([slice(1,11,None)]))
        self.assertRaises(IndexError, df.iloc.__getitem__, tuple([slice(-11,3,None)]))

        # try with labelled frame
        df = DataFrame(np.random.randn(10, 4), index=list('abcdefghij'), columns=list('ABCD'))

        result = df.iloc[1,1]
        exp = df.ix['b','B']
        self.assert_(result == exp)

        result = df.iloc[:,2:3]
        expected = df.ix[:,['C']]
        assert_frame_equal(result, expected)

        # negative indexing
        result = df.iloc[-1,-1]
        exp = df.ix['j','D']
        self.assert_(result == exp)

        # out-of-bounds exception
        self.assertRaises(IndexError, df.iloc.__getitem__, tuple([10,5]))

        # trying to use a label
        self.assertRaises(ValueError, df.iloc.__getitem__, tuple(['j','D']))

    def test_iloc_setitem_series(self):
        """ originally from test_series.py """
        df = DataFrame(np.random.randn(10, 4), index=list('abcdefghij'), columns=list('ABCD'))

        df.iloc[1,1] = 1
        result = df.iloc[1,1]
        self.assert_(result == 1)

        df.iloc[:,2:3] = 0
        expected = df.iloc[:,2:3]
        result = df.iloc[:,2:3]
        assert_frame_equal(result, expected)

    def test_iloc_setitem_series(self):
        s = Series(np.random.randn(10), index=range(0,20,2))

        s.iloc[1] = 1
        result = s.iloc[1]
        self.assert_(result == 1)

        s.iloc[:4] = 0
        expected = s.iloc[:4]
        result = s.iloc[:4]
        assert_series_equal(result, expected)

    def test_iloc_multiindex(self):
        mi_labels = DataFrame(np.random.randn(4, 3), columns=[['i', 'i', 'j'],
                                                              ['A', 'A', 'B']],
                              index=[['i', 'i', 'j', 'k'], ['X', 'X', 'Y','Y']])

        mi_int    = DataFrame(np.random.randn(3, 3),
                              columns=[[2,2,4],[6,8,10]],
                              index=[[4,4,8],[8,10,12]])


        # the first row
        rs = mi_int.iloc[0]
        xp = mi_int.ix[4].ix[8]
        assert_series_equal(rs, xp)

        # 2nd (last) columns
        rs = mi_int.iloc[:,2]
        xp = mi_int.ix[:,2]
        assert_series_equal(rs, xp)

        # corner column
        rs = mi_int.iloc[2,2]
        xp = mi_int.ix[:,2].ix[2]
        self.assert_(rs == xp)

        # this is basically regular indexing
        rs = mi_labels.iloc[2,2]
        xp = mi_labels.ix['j'].ix[:,'j'].ix[0,0]
        self.assert_(rs == xp)

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

    def test_ix_general(self):

        # ix general issues

        # GH 2817
        data={'amount': {0: 700, 1: 600, 2: 222, 3: 333, 4: 444},
              'col': {0: 3.5, 1: 3.5, 2: 4.0, 3: 4.0, 4: 4.0},
              'year': {0: 2012, 1: 2011, 2: 2012, 3: 2012, 4: 2012}}
        df = DataFrame(data).set_index(keys=['col','year'])

        # this should raise correct error
        self.assertRaises(KeyError, df.ix.__getitem__, tuple([4.0,2012]))

        # this is ok
        df.sortlevel(inplace=True)
        df.ix[(4.0,2012)]

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

    def test_setitem_dtype_upcast(self):

        # GH3216
        df = DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df['c'] = np.nan
        self.assert_(df['c'].dtype == np.float64)

        df.ix[0,'c'] = 'foo'
        expected = DataFrame([{"a": 1, "c" : 'foo'}, {"a": 3, "b": 2, "c" : np.nan}])
        assert_frame_equal(df,expected)

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
        self.assert_(result.equals(expected))

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
        ind = ['A', 'A', 'B', 'C']
        df = DataFrame({'test':range(len(ind))}, index=ind)
        rows = ['C', 'B']
        res = df.ix[rows]
        self.assert_(rows == list(res.index))

        res = df.ix[Index(rows)]
        self.assert_(Index(rows).equals(res.index))

        rows = ['C','B','E']
        res = df.ix[rows]
        self.assert_(rows == list(res.index))

        # inconcistent returns for unique/duplicate indices when values are missing
        df = DataFrame(randn(4,3),index=list('ABCD'))
        expected = df.ix[['E']]

        dfnu = DataFrame(randn(5,3),index=list('AABCD'))
        result = dfnu.ix[['E']]
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
        self.assert_(df.iloc[0,2] == '-----')

        #if I look at df, then element [0,2] equals '_'. If instead I type df.ix[idx,'test'], I get '-----', finally by typing df.iloc[0,2] I get '_'.


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

        self.assert_(p.iloc[:3, :3, :3].shape == (3,3,3))
        self.assert_(p.iloc[1, :3, :3].shape == (3,3))
        self.assert_(p.iloc[:3, 1, :3].shape == (3,3))
        self.assert_(p.iloc[:3, :3, 1].shape == (3,3))
        self.assert_(p.iloc[1, 1, :3].shape == (3,))
        self.assert_(p.iloc[1, :3, 1].shape == (3,))
        self.assert_(p.iloc[:3, 1, 1].shape == (3,))

    def test_multi_assign(self):

        # GH 3626, an assignement of a sub-df to a df
        df = DataFrame({'FC':['a','b','a','b','a','b'],
                        'PF':[0,0,0,0,1,1],
                        'col1':range(6),
                        'col2':range(6,12)})
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

    def test_ix_assign_column_mixed(self):
        # GH #1142
        df = DataFrame(tm.getSeriesData())
        df['foo'] = 'bar'

        orig = df.ix[:, 'B'].copy()
        df.ix[:, 'B'] = df.ix[:, 'B'] + 1
        assert_series_equal(df.B, orig + 1)

        # GH 3668, mixed frame with series value
        df = DataFrame({'x':range(10), 'y':range(10,20),'z' : 'bar'})
        expected = df.copy()
        expected.ix[0, 'y'] = 1000
        expected.ix[2, 'y'] = 1200
        expected.ix[4, 'y'] = 1400
        expected.ix[6, 'y'] = 1600
        expected.ix[8, 'y'] = 1800

        df.ix[df.x % 2 == 0, 'y'] = df.ix[df.x % 2 == 0, 'y'] * 100
        assert_frame_equal(df,expected)

    def test_iloc_mask(self):

        # GH 3631, iloc with a mask (of a series) should raise
        df = DataFrame(range(5), list('ABCDE'), columns=['a'])
        mask = (df.a%2 == 0)
        self.assertRaises(ValueError, df.iloc.__getitem__, tuple([mask]))
        mask.index = range(len(mask))
        self.assertRaises(NotImplementedError, df.iloc.__getitem__, tuple([mask]))

        # ndarray ok
        result = df.iloc[np.array([True] * len(mask),dtype=bool)]
        assert_frame_equal(result,df)

        # the possibilities
        locs = np.arange(4)
        nums = 2**locs
        reps = map(bin, nums)
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

        import warnings
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
                except Exception, e:
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
        df = DataFrame({'A' : [1,2,3,4,5,6], 'B' : [3,4,5,6,7,8]}, index = [0,1,0,1,2,3]).sort(axis=0)
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
        self.assert_(result == 'index_name')

        result = df.ix[[0, 1]].index.name
        self.assert_(result == 'index_name')

        result = df.loc[[0, 1]].index.name
        self.assert_(result == 'index_name')

    def test_iloc_non_unique_indexing(self):

        #GH 4017, non-unique indexing (on the axis)
        df = DataFrame({'A' : [0.1] * 3000, 'B' : [1] * 3000})
        idx = np.array(range(30)) * 99
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
        assert_frame_equal(result, expected)

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
        expected = DataFrame([['a',1,1]],index=['A1'],columns=['h1','h3','h5'],).T
        result = df2['A']['A1']
        assert_frame_equal(result,expected)

        # selecting a non_unique from the 2nd level
        expected = DataFrame([['d',4,4],['e',5,5]],index=Index(['B2','B2'],name='sub'),columns=['h1','h3','h5'],).T
        result = df2['A']['B2']
        assert_frame_equal(result,expected)

    def test_non_unique_loc_memory_error(self):

        # GH 4280
        # non_unique index with a large selection triggers a memory error

        columns = list('ABCDEFG')
        def gen_test(l,l2):
            return pd.concat([ DataFrame(randn(l,len(columns)),index=range(l),columns=columns),
                               DataFrame(np.ones((l2,len(columns))),index=[0]*l2,columns=columns) ])


        def gen_expected(df,mask):
            l = len(mask)
            return pd.concat([
                df.take([0],convert=False),
                DataFrame(np.ones((l,len(columns))),index=[0]*l,columns=columns),
                df.take(mask[1:],convert=False) ])

        df = gen_test(900,100)
        self.assert_(not df.index.is_unique)

        mask = np.arange(100)
        result = df.loc[mask]
        expected = gen_expected(df,mask)
        assert_frame_equal(result,expected)

        df = gen_test(900000,100000)
        self.assert_(not df.index.is_unique)

        mask = np.arange(100000)
        result = df.loc[mask]
        expected = gen_expected(df,mask)
        assert_frame_equal(result,expected)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
