#!/usr/bin/env python

import unittest
import functools
import numbers
from itertools import product
import ast

import nose
from nose.tools import assert_raises, assert_tuple_equal
from nose.tools import assert_true, assert_false, assert_equal

from numpy.random import randn, rand
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
from numpy.testing.decorators import slow

import pandas as pd
from pandas.core import common as com
from pandas import DataFrame, Series, Panel
from pandas.util.testing import makeCustomDataframe as mkdf
from pandas.computation.engines import _engines
from pandas.computation.expr import PythonExprVisitor, PandasExprVisitor
from pandas.computation.ops import (_binary_ops_dict, _unary_ops_dict,
                                    _special_case_arith_ops_syms,
                                    _arith_ops_syms, Constant)
import pandas.computation.expr as expr
from pandas.computation import pytables
from pandas.computation.expressions import _USE_NUMEXPR
from pandas.computation.eval import Scope
from pandas.util.testing import (assert_frame_equal, randbool,
                                 assertRaisesRegexp)
from pandas.util.py3compat import PY3


def skip_numexpr_engine(engine):
    if not _USE_NUMEXPR and engine == 'numexpr':
        raise nose.SkipTest


def engine_has_neg_frac(engine):
    return _engines[engine].has_neg_frac


def _eval_single_bin(lhs, cmp1, rhs, engine):
    c = _binary_ops_dict[cmp1]
    if engine_has_neg_frac(engine):
        try:
            result = c(lhs, rhs)
        except ValueError:
            result = np.nan
    else:
        result = c(lhs, rhs)
    return result


def _series_and_2d_ndarray(lhs, rhs):
    return (com.is_series(lhs) and isinstance(rhs, np.ndarray) and rhs.ndim > 1
            or com.is_series(rhs) and isinstance(lhs, np.ndarray) and lhs.ndim
            > 1)


def skip_incompatible_operand(f):
    @functools.wraps(f)
    def wrapper(self, lhs, arith1, rhs, *args, **kwargs):
        if _series_and_2d_ndarray(lhs, rhs):
            self.assertRaises(Exception, pd.eval, 'lhs {0} rhs'.format(arith1),
                              local_dict={'lhs': lhs, 'rhs': rhs},
                              engine=self.engine)
        else:
            f(self, lhs, arith1, rhs, *args, **kwargs)
    return wrapper


_good_arith_ops = tuple(set(_arith_ops_syms) -
                        set(_special_case_arith_ops_syms))


class TestEvalPandas(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cmp_ops = expr._cmp_ops_syms
        cls.cmp2_ops = cls.cmp_ops[::-1]
        cls.bin_ops = expr._bool_ops_syms
        cls.special_case_ops = _special_case_arith_ops_syms
        cls.arith_ops = _good_arith_ops
        cls.unary_ops = '+', '-'

    def set_current_engine(self):
        self.engine = 'numexpr'

    def setup_data(self):
        nan_df1 = DataFrame(rand(10, 5))
        nan_df1[nan_df1 > 0.5] = np.nan
        nan_df2 = DataFrame(rand(10, 5))
        nan_df2[nan_df2 > 0.5] = np.nan

        self.pandas_lhses = (DataFrame(randn(10, 5)), Series(randn(5)),
                             Series([1, 2, np.nan, np.nan, 5]), nan_df1)
        self.pandas_rhses = (DataFrame(randn(10, 5)), Series(randn(5)),
                             Series([1, 2, np.nan, np.nan, 5]), nan_df2)
        self.scalar_lhses = randn(), np.float64(randn()), np.nan
        self.scalar_rhses = randn(), np.float64(randn()), np.nan

        self.lhses = self.pandas_lhses + self.scalar_lhses
        self.rhses = self.pandas_rhses + self.scalar_rhses

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
        for lhs, cmp1, rhs, binop, cmp2 in product(self.lhses, self.cmp_ops,
                                                   self.rhses, self.bin_ops,
                                                   self.cmp2_ops):
            self.check_complex_cmp_op(lhs, cmp1, rhs, binop, cmp2)

    def test_simple_cmp_ops(self):
        bool_lhses = (DataFrame(randbool(size=(10, 5))),
                      Series(randbool((5,))), randbool())
        bool_rhses = (DataFrame(randbool(size=(10, 5))),
                      Series(randbool((5,))), randbool())
        for lhs, rhs, cmp_op in product(bool_lhses, bool_rhses, self.cmp_ops):
            self.check_simple_cmp_op(lhs, cmp_op, rhs)

    @slow
    def test_binary_arith_ops(self):
        for lhs, op, rhs in product(self.lhses, self.arith_ops, self.rhses):
            self.check_binary_arith_op(lhs, op, rhs)

    def test_modulus(self):
        for lhs, rhs in product(self.lhses, self.rhses):
            self.check_modulus(lhs, '%', rhs)

    def test_floor_division(self):
        for lhs, rhs in product(self.lhses, self.rhses):
            self.check_floor_division(lhs, '//', rhs)

    def test_pow(self):
        for lhs, rhs in product(self.lhses, self.rhses):
            self.check_pow(lhs, '**', rhs)

    @slow
    def test_unary_arith_ops(self):
        for unary_op, lhs, arith_op, rhs in product(self.unary_ops, self.lhses,
                                                    self.arith_ops,
                                                    self.rhses):
            self.check_unary_arith_op(lhs, arith_op, rhs, unary_op)

    @slow
    def test_single_invert_op(self):
        for lhs, op, rhs in product(self.lhses, self.cmp_ops, self.rhses):
            self.check_single_invert_op(lhs, op, rhs)

    @slow
    def test_compound_invert_op(self):
        for lhs, op, rhs in product(self.lhses, self.cmp_ops, self.rhses):
            self.check_compound_invert_op(lhs, op, rhs)

    @skip_incompatible_operand
    def check_complex_cmp_op(self, lhs, cmp1, rhs, binop, cmp2):
        ex = '(lhs {cmp1} rhs) {binop} (lhs {cmp2} rhs)'.format(cmp1=cmp1,
                                                                binop=binop,
                                                                cmp2=cmp2)
        lhs_new = _eval_single_bin(lhs, cmp1, rhs, self.engine)
        rhs_new = _eval_single_bin(lhs, cmp2, rhs, self.engine)
        expected = _eval_single_bin(lhs_new, binop, rhs_new, self.engine)
        result = pd.eval(ex, engine=self.engine)
        assert_array_equal(result, expected)

    @skip_incompatible_operand
    def check_simple_cmp_op(self, lhs, cmp1, rhs):
        ex = 'lhs {0} rhs'.format(cmp1)
        expected = _eval_single_bin(lhs, cmp1, rhs, self.engine)
        result = pd.eval(ex, engine=self.engine)
        assert_array_equal(result, expected)

    @skip_incompatible_operand
    def check_binary_arith_op(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)
        result = pd.eval(ex, engine=self.engine)
        expected = _eval_single_bin(lhs, arith1, rhs, self.engine)
        assert_array_equal(result, expected)
        ex = 'lhs {0} rhs {0} rhs'.format(arith1)
        result = pd.eval(ex, engine=self.engine)
        nlhs = _eval_single_bin(lhs, arith1, rhs,
                                self.engine)
        self.check_alignment(result, nlhs, rhs, arith1)

    def check_alignment(self, result, nlhs, ghs, op):
        try:
            nlhs, ghs = nlhs.align(ghs)
        except (ValueError, TypeError, AttributeError):
            # ValueError: series frame or frame series align
            # TypeError, AttributeError: series or frame with scalar align
            pass
        else:
            expected = self.ne.evaluate('nlhs {0} ghs'.format(op))
            assert_array_equal(result, expected)

    # the following 3 tests require special casing

    @skip_incompatible_operand
    def check_modulus(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)
        result = pd.eval(ex, engine=self.engine)
        expected = lhs % rhs
        assert_allclose(result, expected)
        expected = self.ne.evaluate('expected {0} rhs'.format(arith1))
        assert_allclose(result, expected)

    @skip_incompatible_operand
    def check_floor_division(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)

        if self.engine == 'python':
            res = pd.eval(ex, engine=self.engine)
            expected = lhs // rhs
            assert_array_equal(res, expected)
        else:
            self.assertRaises(TypeError, pd.eval, ex, local_dict={'lhs': lhs,
                                                                  'rhs': rhs},
                              engine=self.engine)

    def get_expected_pow_result(self, lhs, rhs):
        try:
            expected = _eval_single_bin(lhs, '**', rhs, self.engine)
        except ValueError as e:
            msg = 'negative number cannot be raised to a fractional power'
            if e.message == msg:
                if self.engine == 'python':
                    raise nose.SkipTest(e.message)
                else:
                    expected = np.nan
            # raise on other, possibly valid ValueErrors
            else:
                raise
        return expected

    @skip_incompatible_operand
    def check_pow(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)
        expected = self.get_expected_pow_result(lhs, rhs)
        result = pd.eval(ex, engine=self.engine)
        assert_array_equal(result, expected)

        ex = '(lhs {0} rhs) {0} rhs'.format(arith1)
        result = pd.eval(ex, engine=self.engine)
        expected = self.get_expected_pow_result(
            self.get_expected_pow_result(lhs, rhs), rhs)
        assert_array_equal(result, expected)

    @skip_incompatible_operand
    def check_single_invert_op(self, lhs, cmp1, rhs):
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

    @skip_incompatible_operand
    def check_compound_invert_op(self, lhs, cmp1, rhs):
        # compound
        ex = '~(lhs {0} rhs)'.format(cmp1)
        if np.isscalar(lhs) and np.isscalar(rhs):
            lhs, rhs = map(lambda x: np.array([x]), (lhs, rhs))
        expected = ~_eval_single_bin(lhs, cmp1, rhs, self.engine)
        result = pd.eval(ex, engine=self.engine)
        assert_array_equal(expected, result)

        # make sure the other engines work the same as this one
        for engine in self.current_engines:
            ev = pd.eval(ex, engine=self.engine)
            assert_array_equal(ev, result)

    @skip_incompatible_operand
    def check_unary_arith_op(self, lhs, arith1, rhs, unary_op):
        # simple
        ex = '{0}lhs'.format(unary_op, arith1)
        f = _unary_ops_dict[unary_op]
        bad_types = np.floating, float, numbers.Real

        if isinstance(lhs, bad_types):
            raise nose.SkipTest("Incompatiable type for ~ operator")
        if isinstance(rhs, bad_types):
            raise nose.SkipTest("Incompatiable type for ~ operator")

        try:
            expected = f(lhs.values)
        except AttributeError:
            expected = f(lhs)
        result = pd.eval(ex, engine=self.engine)
        assert_array_equal(result, expected)

        for engine in self.current_engines:
            assert_array_equal(result, pd.eval(ex, engine=engine))

        ex = '{0}(lhs {1} rhs)'.format(unary_op, arith1)
        result = pd.eval(ex, engine=self.engine)


class TestEvalPython(TestEvalPandas):

    def set_current_engine(self):
        self.engine = 'python'


class TestEvalPandasWithMixedTypeOperands(TestEvalPandas):
    def setup_data(self):
        super(TestEvalPandasWithMixedTypeOperands, self).setup_data()
        self.lhses += randn(10, 5), randn(5)
        self.rhses += randn(10, 5), randn(5)


class TestEvalPythonWithMixedTypeOperands(TestEvalPandasWithMixedTypeOperands):
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
    res = pd.eval(s, engine=engine)
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


INDEX_TYPES = 'i', 'f', 's', 'u', # 'dt',  # 'p'


@slow
def test_series_frame_commutativity():
    args = product(_engines, INDEX_TYPES, INDEX_TYPES, ('+', '*'), ('index',
                                                                    'columns'))
    for engine, r_idx_type, c_idx_type, op, index_name in args:
        check_series_frame_commutativity(engine, r_idx_type, c_idx_type, op,
                                         index_name)


@slow
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


@slow
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


def test_datetime_index_rows_punts_to_python():
    for engine in _engines:
        check_datetime_index_rows_punts_to_python(engine)


def test_truediv():
    for engine in _engines:
        check_truediv(engine)


def check_truediv(engine):
    s = np.array([1])
    ex = 's / 1'

    if PY3:
        res = pd.eval(ex, truediv=False)
        assert_array_equal(res, np.array([1.0]))

        res = pd.eval(ex, truediv=True)
        assert_array_equal(res, np.array([1.0]))
    else:
        res = pd.eval(ex, truediv=False)
        assert_array_equal(res, np.array([1]))

        res = pd.eval(ex, truediv=True)
        assert_array_equal(res, np.array([1.0]))


__var_s = randn(10)


def check_global_scope(engine):
    e = '__var_s * 2'
    assert_array_equal(__var_s * 2, pd.eval(e, engine=engine))


def test_global_scope():
    for engine in _engines:
        yield check_global_scope, engine


def check_is_expr(engine):
    s = 1
    valid = 's + 1'
    invalid = 's +'
    assert_true(expr.isexpr(valid, check_names=True))
    assert_true(expr.isexpr(valid, check_names=False))
    assert_false(expr.isexpr(invalid, check_names=False))
    assert_false(expr.isexpr(invalid, check_names=True))


def test_is_expr():
    for engine in _engines:
        check_is_expr(engine)


_parsers = {'python': PythonExprVisitor, 'pytables': pytables.ExprVisitor,
             'pandas': PandasExprVisitor}


def check_disallowed_nodes(visitor):
    """make sure the disallowed decorator works"""
    VisitorClass = _parsers[visitor]
    uns_ops = VisitorClass.unsupported_nodes
    inst = VisitorClass('x + 1')

    for ops in uns_ops:
        assert_raises(NotImplementedError, getattr(inst, ops))


def test_disallowed_nodes():
    for visitor in _parsers:
        check_disallowed_nodes(visitor)


def check_simple_arith_ops(engine):
    ops = expr._arith_ops_syms + expr._cmp_ops_syms

    for op in filter(lambda x: x != '//', ops):
        expec = _eval_single_bin(1, op, 1, engine)
        x = pd.eval('1 {0} 1'.format(op), engine=engine)
        assert_equal(x, expec)

        expec = _eval_single_bin(x, op, 1, engine)
        y = pd.eval('x {0} 1'.format(op), engine=engine)
        assert_equal(y, expec)

        expec = _eval_single_bin(1, op, x + 1, engine)
        y = pd.eval('1 {0} (x + 1)'.format(op), engine=engine)
        assert_equal(y, expec)


def check_simple_bool_ops(engine):
    for op, lhs, rhs in product(expr._bool_ops_syms, (True, False), (True,
                                                                     False)):
        expec = _eval_single_bin(lhs, op, rhs, engine)
        x = pd.eval('lhs {0} rhs'.format(op), engine=engine)
        assert_equal(x, expec)


def check_bool_ops_with_constants(engine):
    asteval = ast.literal_eval
    for op, lhs, rhs in product(expr._bool_ops_syms, ('True', 'False'),
                                ('True', 'False')):
        expec = _eval_single_bin(asteval(lhs), op, asteval(rhs), engine)
        x = pd.eval('{0} {1} {2}'.format(lhs, op, rhs), engine=engine)
        assert_equal(x, expec)


def test_simple_arith_ops():
    for engine in _engines:
        check_simple_arith_ops(engine)


def test_simple_bool_ops():
    for engine in _engines:
        check_simple_bool_ops(engine)


def test_bool_ops_with_constants():
    for engine in _engines:
        check_bool_ops_with_constants(engine)


def check_no_new_locals(engine):
    x = 1
    lcls = locals().copy()
    pd.eval('x + 1', local_dict=lcls)
    lcls2 = locals().copy()
    lcls2.pop('lcls')
    assert_equal(lcls, lcls2)


def test_no_new_locals():
    for engine in _engines:
        check_no_new_locals(engine)


def check_no_new_globals(engine):
    x = 1
    gbls = globals().copy()
    pd.eval('x + 1')
    gbls2 = globals().copy()
    assert_equal(gbls, gbls2)


def test_no_new_globals():
    for engine in _engines:
        check_no_new_globals(engine)


def check_panel_fails(engine):
    x = Panel(randn(3, 4, 5))
    y = Series(randn(10))
    assert_raises(NotImplementedError, pd.eval, 'x + y', local_dict={'x': x,
                                                                     'y': y},
                  engine=engine)


def test_panel_fails():
    for engine in _engines:
        check_panel_fails(engine)


def check_4d_ndarray_fails(engine):
    x = randn(3, 4, 5, 6)
    y = Series(randn(10))
    assert_raises(NotImplementedError, pd.eval, 'x + y', local_dict={'x': x,
                                                                     'y': y},
                  engine=engine)


def test_4d_ndarray_fails():
    for engine in _engines:
        check_4d_ndarray_fails(engine)


def check_constant(engine):
    x = pd.eval('1', engine=engine)
    assert_equal(x, 1)


def test_constant():
    for engine in _engines:
        check_constant(engine)


def check_single_variable(engine):
    df = DataFrame(randn(10, 2))
    df2 = pd.eval('df', engine=engine)
    assert_frame_equal(df, df2)


def test_single_variable():
    for engine in _engines:
        check_single_variable(engine)


def test_invalid_engine():
    assertRaisesRegexp(KeyError, 'Invalid engine \'asdf\' passed',
                       pd.eval, 'x + y', local_dict={'x': 1, 'y': 2},
                       engine='asdf')


def test_invalid_parser():
    assertRaisesRegexp(KeyError, 'Invalid parser \'asdf\' passed',
                       pd.eval, 'x + y', local_dict={'x': 1, 'y': 2},
                       parser='asdf')


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
