# -*- coding: utf-8 -*-
# pylint: disable-msg=E1101,W0612

import pytest
import numpy as np
from functools import reduce
from collections import OrderedDict
from operator import __or__, __xor__, __and__, __sub__

from pandas import Series, concat
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm
import pandas.core.sets as sets

ops = {'union': __or__, 'xor': __xor__, 'intersect': __and__, 'diff': __sub__}


class TestSetMethods(object):

    def test_api(self):
        assert Series.set is sets.SetMethods
        assert isinstance(Series([{1, 2}]).set, sets.SetMethods)

    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_self(self, opname):
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        exp = reduce(ops[opname], s.values)
        assert getattr(s.set, opname)() == exp

        # with NaN
        t = Series([{2, 3}, {1, 2}, np.nan])
        exp = reduce(ops[opname], t.dropna().values)
        assert getattr(t.set, opname)() == exp

    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_broadcast(self, opname):
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        x = {3}
        exp = Series(ops[opname](s.values, x), index=s.index)
        assert_series_equal(getattr(s.set, opname)(x), exp)

        # with NaN
        t = Series([{2, 3}, {1, 2}, np.nan])
        exp = Series(index=t.index)
        exp.loc[t.notna()] = ops[opname](t.dropna().values, x)
        assert_series_equal(getattr(t.set, opname)(x), exp)

    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_with_series(self, opname):
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        t = Series([{2, 3}, {1, 2}, np.nan])

        exp = Series(index=s.index)
        na_mask = s.isna() | t.isna()
        exp.loc[~na_mask] = ops[opname](s.loc[~na_mask].values,
                                        t.loc[~na_mask].values)

        assert_series_equal(getattr(s.set, opname)(t), exp)

    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_with_1darray(self, opname):
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        t = Series([{2, 3}, {1, 2}, np.nan])

        exp = getattr(s.set, opname)(t)  # tested in test_set_op_with_series
        assert_series_equal(getattr(s.set, opname)(t.values), exp)

        # errors for incorrect lengths
        rgx = 'If `others` contains arrays, these must all be of the same.*'
        with tm.assert_raises_regex(ValueError, rgx):
            getattr(s.set, opname)(t.iloc[:2].values)

    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_parameter_errors(self, opname):
        s = Series([{1, 2}, {2, 4}, {3, 1}])

        rgx = "The parameter 'fill_value' must be list-like!"
        with tm.assert_raises_regex(ValueError, rgx):
            getattr(s.set, opname)(fill_value=1)

        rgx = "Received illegal value for parameter 'errors'.*"
        with tm.assert_raises_regex(ValueError, rgx):
            getattr(s.set, opname)(errors='abcd')

        rgx = "Must pass Series for validating inputs of set accessor.*"
        with tm.assert_raises_regex(ValueError, rgx):
            s.set._validate(s.values)

    @pytest.mark.parametrize('fill_value', [None, set(), {5}])
    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_fill_value(self, opname, fill_value):
        s = Series([{1, 2}, np.nan, {3, 1}])
        t = Series([{2, 3}, {1, 2}, np.nan])

        # cannot use fillna with sets due to GH21329
        sf = s.copy()
        sf.loc[sf.isna()] = [fill_value] * sf.isna().sum()
        tf = t.copy()
        tf.loc[tf.isna()] = [fill_value] * tf.isna().sum()

        exp = getattr(sf.set, opname)(tf)
        assert_series_equal(getattr(s.set, opname)(t, fill_value=fill_value),
                            exp)

    @pytest.mark.parametrize('errors',
                             ['raise', 'ignore', 'coerce', 'wrap', 'skip'])
    @pytest.mark.parametrize('opname', ops.keys())
    def test_conversion_caller(self, opname, errors):
        u = Series([{1, 2}, [2, 4], 'abcd', 5, np.nan])
        ui = Series([{1, 2}, {2, 4}, np.nan, np.nan, np.nan])  # ignore
        uc = Series([{1, 2}, {2, 4}, set('abcd'), np.nan, np.nan])  # coerce
        uw = Series([{1, 2}, {2, 4}, {'abcd'}, {5}, np.nan])  # wrap
        x = {3}

        if errors == 'raise':
            rgx = '.*can only use .set accessor with values that.*'
            with tm.assert_raises_regex(ValueError, rgx):
                getattr(u.set, opname)(errors=errors)
        elif errors == 'skip':  # raw error from numpy
            rgx = '.*unsupported operand type.*'
            with tm.assert_raises_regex(TypeError, rgx):
                getattr(u.set, opname)(errors=errors)
            # but if series is already converted (like uw/uc here), skip works
            exp = getattr(uc.set, opname)(uw)  # default: errors='raise'
            assert_series_equal(getattr(uc.set, opname)(uw, errors=errors),
                                exp)
        else:  # 'ignore', 'coerce', 'wrap'
            if errors == 'ignore':
                u_exp = ui
            elif errors == 'coerce':
                u_exp = uc
            else:  # 'wrap'
                u_exp = uw

            # apply to self
            exp = getattr(u_exp.set, opname)()
            assert getattr(u.set, opname)(errors=errors) == exp

            # apply to single set with broadcasting
            exp = getattr(u_exp.set, opname)(x)
            assert_series_equal(getattr(u.set, opname)(x, errors=errors), exp)

    @pytest.mark.parametrize('errors',
                             ['raise', 'ignore', 'coerce', 'wrap', 'skip'])
    @pytest.mark.parametrize('opname', ops.keys())
    def test_conversion_others(self, opname, errors):
        s = Series([{1}, {2}, {3}, {4}, {5}])
        u = Series([{1, 2}, [2, 4], 'abcd', 5, np.nan])

        if errors == 'raise':
            rgx = '.*can only use .set accessor with values that.*'
            with tm.assert_raises_regex(ValueError, rgx):
                getattr(s.set, opname)(u, errors=errors)
        elif errors == 'skip':  # raw error from numpy
            rgx = '.*unsupported operand type.*'
            with tm.assert_raises_regex(TypeError, rgx):
                getattr(s.set, opname)(u, errors=errors)
            # but if series is already converted, skip works
            u_conv = u.set.union(set(), errors='wrap')  # explanation below
            exp = getattr(s.set, opname)(u_conv)  # default: errors='raise'
            assert_series_equal(getattr(s.set, opname)(u_conv, errors=errors),
                                exp)
        else:  # 'ignore', 'coerce', 'wrap'
            # union with set() does not change sets, only applies conversion
            # correctness of this behavior tested in test_conversion_caller
            u_exp = u.set.union(set(), errors=errors)

            exp = getattr(s.set, opname)(u_exp)
            assert_series_equal(getattr(s.set, opname)(u, errors=errors), exp)

    @pytest.mark.parametrize('fill_value', [None, set(), {5}])
    @pytest.mark.parametrize('join', ['left', 'outer', 'inner', 'right'])
    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_align(self, opname, join, fill_value):
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        t = Series([{2, 3}, {1, 2}, np.nan], index=[3, 2, 1])

        sa, ta = s.align(t, join=join)
        exp = getattr(sa.set, opname)(ta, fill_value=fill_value)
        assert_series_equal(getattr(s.set, opname)(t, join=join,
                                                   fill_value=fill_value), exp)

    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_mixed_inputs(self, opname):
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        t = Series([None, {1, 2}, np.nan])  # test if NaN/None diff. matters
        d = concat([t, s], axis=1)

        # all methods below are equivalent to sequential application
        # (at least when indexes in rhs are all the same)
        tmp = getattr(s.set, opname)(t)
        exp = getattr(tmp.set, opname)(s)

        # Series with DataFrame
        assert_series_equal(getattr(s.set, opname)(d), exp)

        # Series with two-dimensional ndarray
        assert_series_equal(getattr(s.set, opname)(d.values), exp)

        # Series with list of Series
        assert_series_equal(getattr(s.set, opname)([t, s]), exp)

        # Series with mixed list of Series/ndarray
        assert_series_equal(getattr(s.set, opname)([t, s.values]), exp)

        # Series with iterator of Series
        assert_series_equal(getattr(s.set, opname)(iter([t, s])), exp)

        # Series with dict_view of Series
        dv = d.to_dict('series', into=OrderedDict).values()
        assert_series_equal(getattr(s.set, opname)(dv), exp)

        # errors for incorrect lengths
        rgx = 'If `others` contains arrays, these must all be'

        # two-dimensional ndarray
        with tm.assert_raises_regex(ValueError, rgx):
            getattr(s.set, opname)(d.iloc[:2].values)

        # mixed list with Series/ndarray
        with tm.assert_raises_regex(ValueError, rgx):
            getattr(s.set, opname)([t, s.iloc[:2].values])

        # errors for incorrect arguments in list-like
        rgx = 'others must be Series, DataFrame, np.ndarrary or list-like.*'

        # mix of string and Series
        with tm.assert_raises_regex(TypeError, rgx):
            getattr(s.set, opname)([t, 't'])

        # DataFrame in list
        with tm.assert_raises_regex(TypeError, rgx):
            getattr(s.set, opname)([t, d])

        # 2-dim ndarray in list
        with tm.assert_raises_regex(TypeError, rgx):
            getattr(s.set, opname)([t, d.values])

        # nested lists
        with tm.assert_raises_regex(TypeError, rgx):
            getattr(s.set, opname)([t, [t, t]])

        # forbidden input type, e.g. int
        with tm.assert_raises_regex(TypeError, rgx):
            getattr(s.set, opname)(1)

    @pytest.mark.parametrize('fill_value', [None, set(), {5}])
    @pytest.mark.parametrize('join', ['left', 'outer', 'inner', 'right'])
    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_align_several(self, opname, join, fill_value):
        # no differences in the indexes of the right-hand side yet!
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        t = Series([{2, 3}, {1, 2}, np.nan], index=[3, 2, 1])
        d = concat([t, t], axis=1)

        sa, ta = s.align(t, join=join)
        exp = getattr(sa.set, opname)([ta, ta], fill_value=fill_value)

        # list of Series
        tm.assert_series_equal(getattr(s.set, opname)([t, t], join=join,
                                                      fill_value=fill_value),
                               exp)

        # DataFrame
        tm.assert_series_equal(getattr(s.set, opname)(d, join=join,
                                                      fill_value=fill_value),
                               exp)

    @pytest.mark.parametrize('fill_value', [None, set(), {5}])
    @pytest.mark.parametrize('join', ['left', 'outer', 'inner', 'right'])
    @pytest.mark.parametrize('opname', ops.keys())
    def test_set_op_align_mixed(self, opname, join, fill_value):
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        t = Series([{2, 3}, {1, 2}, np.nan], index=[3, 2, 1])
        u = Series([{5}, {3}, None], index=[2, 4, 1])

        # the index of the right-hand side is the union of the rhs indexes,
        # except for 'inner' - this is only really relevant for 'right', which
        # would not have a well-defined index otherwise.
        ta, ua = t.align(u, join=(join if join == 'inner' else 'outer'))

        # reuse case of same rhs-index; tested in test_set_op_align_several
        exp = getattr(s.set, opname)([ta, ua], join=join,
                                     fill_value=fill_value)

        # list of Series
        tm.assert_series_equal(getattr(s.set, opname)([t, u], join=join,
                                                      fill_value=fill_value),
                               exp)

        # unindexed -> use index of caller
        # reuses test directly above that differently-indexed series work
        tu = Series(t.values, index=s.index)
        exp = getattr(s.set, opname)([tu, u], join=join,
                                     fill_value=fill_value)

        # mixed list of indexed/unindexed
        tm.assert_series_equal(getattr(s.set, opname)([t.values, u], join=join,
                                                      fill_value=fill_value),
                               exp)

    @pytest.mark.parametrize('fill_value', [set(), {5}])
    def test_set_len(self, fill_value):
        s = Series([{1, 2}, {2, 4}, {3, 1}])
        t = Series([{2, 3}, {1, 2}, np.nan])

        rgx = "The parameter 'fill_value' must be list-like!"
        with tm.assert_raises_regex(ValueError, rgx):
            s.set.len(fill_value=1)

        rgx = "Received illegal value for parameter 'errors'.*"
        with tm.assert_raises_regex(ValueError, rgx):
            s.set.len(errors='abcd')

        # no NaN
        exp = s.map(len)
        assert_series_equal(s.set.len(fill_value=fill_value), exp)

        # cannot use fillna with sets due to GH21329
        tf = t.copy()
        tf.loc[tf.isna()] = [fill_value] * tf.isna().sum()
        exp = tf.map(len)
        assert_series_equal(t.set.len(fill_value=fill_value), exp)

    @pytest.mark.parametrize('errors',
                             ['raise', 'ignore', 'coerce', 'wrap', 'skip'])
    def test_set_len_conversion(self, errors):
        u = Series([{1, 2}, [2, 4], 'abcd', 5, np.nan])
        ui = Series([{1, 2}, {2, 4}, np.nan, np.nan, np.nan])  # ignore
        uc = Series([{1, 2}, {2, 4}, set('abcd'), np.nan, np.nan])  # coerce
        uw = Series([{1, 2}, {2, 4}, {'abcd'}, {5}, np.nan])  # wrap

        if errors == 'raise':
            rgx = '.*can only use .set accessor with values that.*'
            with tm.assert_raises_regex(ValueError, rgx):
                u.set.len(errors=errors)
        elif errors == 'skip':  # raw error from trying len(5)
            rgx = "object of type 'int' has no len()"
            with tm.assert_raises_regex(TypeError, rgx):
                u.set.len(errors=errors)
            # but if series is already converted (like uc here), skip works
            exp = uc.set.len()  # default: errors='raise'
            assert_series_equal(uc.set.len(errors=errors), exp)
        else:  # 'ignore', 'coerce', 'wrap'
            if errors == 'ignore':
                u_exp = ui
            elif errors == 'coerce':
                u_exp = uc
            else:  # 'wrap'
                u_exp = uw

            assert_series_equal(u.set.len(errors=errors), u_exp.set.len())
