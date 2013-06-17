#!/usr/bin/env python

import unittest
import itertools
from itertools import product

import nose
from nose.tools import assert_raises, assert_tuple_equal, assert_equal
from nose.tools import assert_true

from numpy.random import randn, rand
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
from numpy.testing.decorators import slow

import pandas as pd
from pandas.core import common as com
from pandas import DataFrame, Series
from pandas.util.testing import makeCustomDataframe as mkdf
from pandas.computation.engines import (_engines, _align_core,
                                        _reconstruct_object)
from pandas.computation.ops import _binary_ops_dict, _unary_ops_dict
import pandas.computation.expr as expr
from pandas.computation.expressions import _USE_NUMEXPR
from pandas.computation.eval import Scope
from pandas.computation.eval import _scope_has_series_and_frame_datetime_index
from pandas.computation.eval import _maybe_convert_engine
from pandas.util.testing import assert_frame_equal, randbool


def skip_numexpr_engine(engine):
    if not _USE_NUMEXPR and engine == 'numexpr':
        raise nose.SkipTest


def engine_has_neg_frac(engine):
    return _engines[engine].has_neg_frac


def fractional(x):
    frac, _ = np.modf(np.asanyarray(x))
    return frac


def hasfractional(x):
    return np.any(fractional(x) != 0.0)


def _eval_from_expr(lhs, cmp1, rhs, binop, cmp2):
    f1 = _binary_ops_dict[cmp1]
    f2 = _binary_ops_dict[cmp2]
    bf = _binary_ops_dict[binop]
    typ, (lhs, rhs), axes = _align_core((lhs, rhs))
    return _reconstruct_object(typ, bf(f1(lhs, rhs), f2(lhs, rhs)), axes)


def _eval_single_bin(lhs, cmp1, rhs, has_neg_frac):
    c = _binary_ops_dict[cmp1]
    if has_neg_frac:
        try:
            result = c(lhs, rhs)
        except ValueError:
            result = np.nan
    else:
        result = c(lhs, rhs)
    return result


def isframe(x):
    return isinstance(x, pd.DataFrame)


def isseries(x):
    return isinstance(x, pd.Series)


def are_compatible_types(op, lhs, rhs):
    if op in ('&', '|'):
        if isframe(lhs) and isseries(rhs) or isframe(rhs) and isseries(lhs):
            return False
    return True


def _eval_bin_and_unary(unary, lhs, arith1, rhs):
    binop = _binary_ops_dict[arith1]
    unop = expr._unary_ops_dict[unary]
    return unop(binop(lhs, rhs))


def _series_and_2d_ndarray(lhs, rhs):
    return (com.is_series(lhs) and isinstance(rhs, np.ndarray) and rhs.ndim > 1
            or com.is_series(rhs) and isinstance(lhs, np.ndarray) and lhs.ndim
            > 1)


# Smoke testing
class TestBasicEval(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cmp_ops = expr._cmp_ops_syms
        self.cmp2_ops = self.cmp_ops[::-1]
        self.bin_ops = expr._bool_ops_syms
        self.arith_ops = tuple(o for o in expr._arith_ops_syms if o != '//')
        self.unary_ops = '+', '-'

    def set_current_engine(self):
        self.engine = 'numexpr'

    def setup_data(self):
        nan_df = DataFrame(rand(10, 5))
        nan_df[nan_df > 0.5] = np.nan
        self.lhses = (DataFrame(randn(10, 5)), Series(randn(5)), randn(),
                      np.float64(randn()), randn(10, 5), randn(5), np.nan,
                      Series([1, 2, np.nan, np.nan, 5]), nan_df)
        self.rhses = (DataFrame(randn(10, 5)), Series(randn(5)), randn(),
                      np.float64(randn()), randn(10, 5), randn(5), np.nan,
                      Series([1, 2, np.nan, np.nan, 5]), nan_df)

    def setUp(self):
        try:
            import numexpr as ne
            self.ne = ne
        except ImportError:
            raise nose.SkipTest
        self.set_current_engine()
        self.setup_data()
        self.current_engines = filter(lambda x: x != self.engine,
                                      _engines.iterkeys())

    @slow
    def test_complex_cmp_ops(self):
        self.setUp()
        lhses, rhses = self.lhses, self.rhses
        args = itertools.product(lhses, self.cmp_ops, rhses, self.bin_ops,
                                 self.cmp2_ops)
        for lhs, cmp1, rhs, binop, cmp2 in args:
            self._create_cmp_op_t(lhs, cmp1, rhs, binop, cmp2)

    def test_simple_cmp_ops(self):
        bool_lhses = (DataFrame(randbool(size=(10, 5))),
                      Series(randbool((5,))), randbool())
        bool_rhses = (DataFrame(randbool(size=(10, 5))),
                      Series(randbool((5,))), randbool())
        args = itertools.product(bool_lhses, bool_rhses, self.cmp_ops)
        for lhs, rhs, cmp_op in args:
            self._create_simple_cmp_op_t(lhs, rhs, cmp_op)

    def test_binary_arith_ops(self):
        self.setUp()
        lhses = DataFrame(randn(10, 5)), Series(randn(5)), randn()
        rhses = DataFrame(randn(10, 5)), Series(randn(5)), randn()
        args = itertools.product(lhses, self.arith_ops, rhses)
        for lhs, op, rhs in args:
            self._create_arith_op_t(lhs, op, rhs)

    def test_unary_arith_ops(self):
        self.setUp()
        lhses = DataFrame(randn(10, 5)), Series(randn(5)), randn()
        rhses = DataFrame(randn(10, 5)), Series(randn(5)), randn()
        aops = tuple(aop for aop in self.arith_ops if aop not in '+-')
        args = itertools.product(self.unary_ops, lhses, aops, rhses)
        for unary_op, lhs, arith_op, rhs in args:
            self._create_unary_arith_op_t(unary_op, lhs, arith_op, rhs)

    def test_invert(self):
        self.setUp()
        lhses = DataFrame(randn(10, 5)), Series(randn(5)), randn()
        rhses = DataFrame(randn(10, 5)), Series(randn(5)), randn()
        args = itertools.product(lhses, self.cmp_ops, rhses)
        for lhs, op, rhs in args:
            self._create_invert_op_t(lhs, op, rhs)

    def _create_cmp_op_t(self, lhs, cmp1, rhs, binop, cmp2):
        ex = '(lhs {cmp1} rhs) {binop} (lhs {cmp2} rhs)'.format(cmp1=cmp1,
                                                                binop=binop,
                                                                cmp2=cmp2)
        if _series_and_2d_ndarray(lhs, rhs):
            self.assertRaises(Exception, _eval_from_expr, lhs, cmp1, rhs,
                              binop, cmp2)
            self.assertRaises(Exception, pd.eval, ex, engine=self.engine)
        else:
            expected = _eval_from_expr(lhs, cmp1, rhs, binop, cmp2)
            result = pd.eval(ex, engine=self.engine)
            assert_array_equal(result, expected)

    def _create_simple_cmp_op_t(self, lhs, rhs, cmp1):
        ex = 'lhs {0} rhs'.format(cmp1)

        if are_compatible_types(cmp1, lhs, rhs):
            expected = _eval_single_bin(lhs, cmp1, rhs,
                                        engine_has_neg_frac(self.engine))
            result = pd.eval(ex, engine=self.engine)
            assert_array_equal(result, expected)
        else:
            assert_raises(TypeError, _eval_single_bin, lhs, cmp1, rhs,
                          engine_has_neg_frac(self.engine))

    def _create_arith_op_t(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)
        nan_frac_neg = (arith1 == '**' and np.any(lhs < 0) and
                        hasfractional(rhs) and np.isscalar(lhs) and
                        np.isscalar(rhs) and
                        not (isinstance(lhs, tuple(np.typeDict.values()))
                             or isinstance(rhs, tuple(np.typeDict.values()))))
        if nan_frac_neg and not engine_has_neg_frac(self.engine):
                assert_raises(ValueError, pd.eval, ex, engine=self.engine,
                              local_dict=locals(), global_dict=globals())
        else:
            result = pd.eval(ex, engine=self.engine)

            if arith1 != '//':
                expected = _eval_single_bin(lhs, arith1, rhs,
                                            engine_has_neg_frac(self.engine))
                # roundoff error with modulus
                if arith1 == '%':
                    assert_allclose(result, expected)
                else:
                    assert_array_equal(result, expected)

            # sanity check on recursive parsing
            try:
                ghs = rhs.copy()
            except AttributeError:
                ghs = rhs

        if nan_frac_neg and not engine_has_neg_frac(self.engine):
            assert_raises(ValueError, pd.eval, ex, engine=self.engine,
                          local_dict=locals(), global_dict=globals())
        else:
            if arith1 == '**':
                ex = '(lhs {0} rhs) {0} ghs'.format(arith1)
            else:
                ex = 'lhs {0} rhs {0} ghs'.format(arith1)
            result = pd.eval(ex, engine=self.engine)

            try:
                nlhs = _eval_single_bin(lhs, arith1, rhs,
                                        engine_has_neg_frac(self.engine))
            except ValueError:
                assert_raises(ValueError, _eval_single_bin, lhs, arith1, rhs,
                              engine_has_neg_frac(self.engine))
            else:
                try:
                    nlhs, ghs = nlhs.align(ghs)
                except:
                    pass
                if arith1 != '//':
                    expected = self.ne.evaluate('nlhs {0} ghs'.format(arith1))

                    # roundoff error with modulus
                    if arith1 == '%':
                        assert_allclose(result, expected)
                    else:
                        assert_array_equal(result, expected)

    def _create_invert_op_t(self, lhs, cmp1, rhs):
        # simple
        for el in (lhs, rhs):
            try:
                elb = el.astype(bool)
            except AttributeError:
                elb = np.array([bool(el)])
            expected = ~elb
            result = pd.eval('~elb', engine=self.engine)
            assert_array_equal(expected, result)

            for engine in self.current_engines:
                assert_array_equal(result, pd.eval('~elb', engine=engine))

        # compound
        ex = '~(lhs {0} rhs)'.format(cmp1)
        if np.isscalar(lhs) and np.isscalar(rhs):
            lhs, rhs = map(lambda x: np.array([x]), (lhs, rhs))
        expected = ~_eval_single_bin(lhs, cmp1, rhs,
                                     engine_has_neg_frac(self.engine))
        result = pd.eval(ex, engine=self.engine)
        assert_array_equal(expected, result)

        # make sure the other engines work
        for engine in self.current_engines:
            ev = pd.eval(ex, engine=self.engine)
            assert_array_equal(ev, result)

    def _create_unary_arith_op_t(self, unary_op, lhs, arith1, rhs):
        # simple
        ex = '{0}lhs'.format(unary_op, arith1)
        f = _unary_ops_dict[unary_op]
        bad_types = tuple(np.typeDict.values())

        nan_frac_neg = (arith1 == '**' and
                        np.any(lhs < 0) and
                        hasfractional(rhs) and
                        np.isscalar(lhs) and np.isscalar(rhs) and
                        not (isinstance(lhs, bad_types) or
                             isinstance(rhs, bad_types))
                        and not engine_has_neg_frac(self.engine))
        try:
            expected = f(lhs.values)
        except AttributeError:
            expected = f(lhs)
        result = pd.eval(ex, engine=self.engine)
        assert_array_equal(result, expected)

        for engine in self.current_engines:
            assert_array_equal(result, pd.eval(ex, engine=engine))

        ex = '{0}(lhs {1} rhs)'.format(unary_op, arith1)

        if nan_frac_neg:
            assert_raises(ValueError, pd.eval, ex, engine=self.engine,
                          local_dict=locals(), global_dict=globals())
        else:
            # compound
            result = pd.eval(ex, engine=self.engine)

            #(lhs, rhs), _ = _align((lhs, rhs))
            #if arith1 != '//':
                #expected = self.ne.evaluate(ex)
                #assert_array_equal(result, expected)
            #else:
                #assert_raises(TypeError, self.ne.evaluate, ex)

            #for engine in self.current_engines:
                #if arith1 != '//':
                    #if engine_has_neg_frac(engine):
                        #assert_array_equal(result, pd.eval(ex, engine=engine))
                #else:
                    #assert_raises(TypeError, pd.eval, ex, engine=engine,
                                  #local_dict=locals(), global_dict=globals())


class TestBasicEvalPython(TestBasicEval):

    @classmethod
    def setUpClass(cls):
        cls.cmp_ops = expr._cmp_ops_syms
        cls.cmp2_ops = cls.cmp_ops[::-1]
        cls.bin_ops = expr._bool_ops_syms
        cls.arith_ops = expr._arith_ops_syms
        cls.unary_ops = '+', '-'

    def set_current_engine(self):
        self.engine = 'python'


def test_syntax_error_exprs():
    for engine in _engines:
        e = 's +'
        assert_raises(SyntaxError, pd.eval, e, engine=engine)


def test_name_error_exprs():
    for engine in _engines:
        e = 's + t'
        assert_raises(NameError, pd.eval, e, engine=engine)


def test_align_nested_unary_op():
    for engine in _engines:
        yield check_align_nested_unary_op, engine


f = lambda *args, **kwargs: np.random.randn()


def check_align_nested_unary_op(engine):
    skip_numexpr_engine(engine)
    s = 'df * ~2'
    df = mkdf(10, 10, data_gen_f=f)
    res = pd.eval(s, engine)
    assert_frame_equal(res, df * ~2)


def check_basic_frame_alignment(engine):
    df = mkdf(10, 10, data_gen_f=f)
    df2 = mkdf(20, 10, data_gen_f=f)
    res = pd.eval('df + df2', engine=engine)
    assert_frame_equal(res, df + df2)


def test_basic_frame_alignment():
    for engine in _engines:
        yield check_basic_frame_alignment, engine


def check_medium_complex_frame_alignment(engine, r1, r2, c1, c2):
    skip_numexpr_engine(engine)
    df = mkdf(5, 2, data_gen_f=f, r_idx_type=r1, c_idx_type=c1)
    df2 = mkdf(10, 2, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
    df3 = mkdf(15, 2, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
    res = pd.eval('df + df2 + df3', engine=engine)
    assert_frame_equal(res, df + df2 + df3)


@slow
def test_medium_complex_frame_alignment():
    args = product(_engines, *([INDEX_TYPES[:4]] * 4))
    for engine, r1, r2, c1, c2 in args:
        check_medium_complex_frame_alignment(engine, r1, r2, c1, c2)


def check_basic_frame_series_alignment(engine, r_idx_type, c_idx_type,
                                       index_name):
    skip_numexpr_engine(engine)
    df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
              c_idx_type=c_idx_type)
    index = getattr(df, index_name)
    s = Series(np.random.randn(5), index[:5])

    if r_idx_type != 'p' and c_idx_type == 'p' and index_name == 'index':
        assert_raises(ValueError, pd.eval, 'df + s', local_dict=locals())
        assert_raises(ValueError, df.add, s, axis=1)
    else:
        res = pd.eval('df + s', engine=engine)
        expected = df + s
        assert_frame_equal(res, expected)


def check_not_both_period_fails_otherwise_succeeds(lhs, rhs, r_idx_type,
                                                   c_idx_type, index_name, s,
                                                   df, *terms):
    if r_idx_type != 'p' and c_idx_type == 'p' and index_name == 'index':
        assert_raises(ValueError, pd.eval, lhs, local_dict=locals())
        assert_raises(ValueError, pd.eval, rhs, local_dict=locals())
    else:
        a, b = pd.eval(lhs), pd.eval(rhs)
        assert_frame_equal(a, b)


def check_basic_series_frame_alignment(engine, r_idx_type, c_idx_type,
                                       index_name):
    skip_numexpr_engine(engine)
    df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
              c_idx_type=c_idx_type)
    index = getattr(df, index_name)
    s = Series(np.random.randn(5), index[:5])

    if r_idx_type != 'p' and c_idx_type == 'p' and index_name == 'index':
        assert_raises(ValueError, pd.eval, 's + df', local_dict=locals())
        assert_raises(ValueError, df.add, s, axis=1)
    else:
        res = pd.eval('s + df', engine=engine)
        expected = s + df
        assert_frame_equal(res, expected)


@slow
def check_basic_series_frame_alignment_datetime(engine, r_idx_type, c_idx_type,
                                                index_name):
    skip_numexpr_engine(engine)
    df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
              c_idx_type=c_idx_type)
    index = getattr(df, index_name)
    s = Series(np.random.randn(5), index[:5])
    if r_idx_type != 'p' and c_idx_type == 'p' and index_name == 'index':
        assert_raises(ValueError, pd.eval, 's + df', local_dict=locals())
        assert_raises(ValueError, df.add, s, axis=1)
    else:
        res = pd.eval('s + df', engine=engine)
        expected = s + df
        assert_frame_equal(res, expected)

    if r_idx_type != 'p' and c_idx_type == 'p' and index_name == 'index':
        assert_raises(ValueError, pd.eval, 'df + s', local_dict=locals())
        assert_raises(ValueError, df.add, s, axis=1)
    else:
        res = pd.eval('df + s', engine=engine)
        expected = df + s
        assert_frame_equal(res, expected)


def check_series_frame_commutativity(engine, r_idx_type, c_idx_type, op,
                                     index_name):
    skip_numexpr_engine(engine)
    df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
              c_idx_type=c_idx_type)
    index = getattr(df, index_name)
    s = Series(np.random.randn(5), index[:5])

    lhs = 's {0} df'.format(op)
    rhs = 'df {0} s'.format(op)
    check_not_both_period_fails_otherwise_succeeds(lhs, rhs, r_idx_type,
                                                   c_idx_type, index_name, s,
                                                   df)


INDEX_TYPES = 'i', 'f', 's', 'u', 'dt',  # 'p'


@slow
def test_series_frame_commutativity():
    args = product(_engines, INDEX_TYPES, INDEX_TYPES, ('+', '*'), ('index',
                                                                    'columns'))
    for engine, r_idx_type, c_idx_type, op, index_name in args:
        check_series_frame_commutativity(engine, r_idx_type, c_idx_type, op,
                                         index_name)


def test_basic_frame_series_alignment():
    args = product(_engines, INDEX_TYPES, INDEX_TYPES, ('index', 'columns'))
    for engine, r_idx_type, c_idx_type, index_name in args:
        check_basic_frame_series_alignment(engine, r_idx_type, c_idx_type,
                                           index_name)


@slow
def test_basic_series_frame_alignment_datetime():
    idx_types = INDEX_TYPES
    args = product(_engines, idx_types, idx_types, ('index', 'columns'))
    for engine, r_idx_type, c_idx_type, index_name in args:
        check_basic_series_frame_alignment_datetime(engine, r_idx_type,
                                                    c_idx_type, index_name)


def test_basic_series_frame_alignment():
    args = product(_engines, INDEX_TYPES, INDEX_TYPES, ('index', 'columns'))
    for engine, r_idx_type, c_idx_type, index_name in args:
        check_basic_series_frame_alignment(engine, r_idx_type, c_idx_type,
                                           index_name)


def check_complex_series_frame_alignment(engine, index_name, obj, r1, r2, c1,
                                         c2):
    skip_numexpr_engine(engine)
    df = mkdf(10, 10, data_gen_f=f, r_idx_type=r1, c_idx_type=c1)
    df2 = mkdf(20, 10, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
    index = getattr(locals()[obj], index_name)
    s = Series(np.random.randn(5), index[:5])
    if engine != 'python':
        expected = df2.add(s, axis=1).add(df)
    else:
        expected = df2 + s + df
    res = pd.eval('df2 + s + df', engine=engine)
    expected = df2 + s + df
    assert_tuple_equal(res.shape, expected.shape)
    assert_frame_equal(res, expected)


@slow
def test_complex_series_frame_alignment():
    args = product(_engines, ('index', 'columns'), ('df', 'df2'),
                   *([INDEX_TYPES[:4]] * 4))
    for engine, index_name, obj, r1, r2, c1, c2 in args:
        check_complex_series_frame_alignment(engine, index_name, obj, r1, r2,
                                             c1, c2)


def check_datetime_index_rows_punts_to_python(engine):
    df = mkdf(10, 10, data_gen_f=f, r_idx_type='dt', c_idx_type='dt')
    index = getattr(df, 'index')
    s = Series(np.random.randn(5), index[:5])
    env = Scope(globals(), locals())
    assert_true(_scope_has_series_and_frame_datetime_index(env))
    assert_equal(_maybe_convert_engine(env, engine), 'python')


def test_datetime_index_rows_punts_to_python():
    for engine in _engines:
        check_datetime_index_rows_punts_to_python(engine)


def test_truediv():
    for engine in _engines:
        check_truediv(engine)


def check_truediv(engine):
    s = randn(10)


__var_s = randn(10)


def check_global_scope(engine):
    e = '__var_s * 2'
    assert_array_equal(__var_s * 2, pd.eval(e, engine=engine))


def test_global_scope():
    for engine in _engines:
        yield check_global_scope, engine


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
