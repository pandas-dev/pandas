# pylint: disable-msg=W0612,E1101
import unittest
import nose
import itertools

from numpy import random, nan
from numpy.random import randn
import numpy as np
from numpy.testing import assert_array_equal

import pandas as pan
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


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
