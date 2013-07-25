#!/usr/bin/env python

import unittest
import functools
import numbers
from itertools import product
import ast

import nose
from nose.tools import assert_raises, assert_true, assert_false, assert_equal

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
                                    _arith_ops_syms)
import pandas.computation.expr as expr
from pandas.computation import pytables
from pandas.computation.expressions import _USE_NUMEXPR
from pandas.util.testing import (assert_frame_equal, randbool,
                                 assertRaisesRegexp,
                                 assert_produces_warning, assert_series_equal)
from pandas.util.py3compat import PY3


def skip_numexpr_engine(engine):
    if not _USE_NUMEXPR and engine == 'numexpr':
        raise nose.SkipTest("not using numexpr")


def engine_has_neg_frac(engine):
    return _engines[engine].has_neg_frac


def _eval_single_bin(lhs, cmp1, rhs, engine):
    c = _binary_ops_dict[cmp1]
    if engine_has_neg_frac(engine):
        try:
            return c(lhs, rhs)
        except ValueError as e:
            if e.message == ('negative number cannot be raised to a '
                             'fractional power'):
                return np.nan
            raise
    return c(lhs, rhs)


def _series_and_2d_ndarray(lhs, rhs):
    return (com.is_series(lhs) and isinstance(rhs, np.ndarray) and rhs.ndim > 1
            or com.is_series(rhs) and isinstance(lhs, np.ndarray) and lhs.ndim
            > 1)


def _bool_and_frame(lhs, rhs):
    return isinstance(lhs, bool) and com.is_frame(rhs)


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


_good_arith_ops = com.difference(_arith_ops_syms, _special_case_arith_ops_syms)

class TestEvalPandas(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cmp_ops = expr._cmp_ops_syms
        cls.cmp2_ops = cls.cmp_ops[::-1]
        cls.bin_ops = expr._bool_ops_syms
        cls.special_case_ops = _special_case_arith_ops_syms
        cls.arith_ops = _good_arith_ops
        cls.unary_ops = '+', '-'

    @classmethod
    def tearDownClass(cls):
        del cls.cmp_ops, cls.cmp2_ops, cls.bin_ops, cls.special_case_ops
        del cls.arith_ops, cls.unary_ops

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

        self.lhses = self.pandas_lhses + self.scalar_lhses + (randn(10, 5),
                                                              randn(5))
        self.rhses = self.pandas_rhses + self.scalar_rhses + (randn(10, 5),
                                                              randn(5))

    def setUp(self):
        try:
            import numexpr as ne
            self.ne = ne
        except ImportError:
            raise nose.SkipTest
        self.set_current_engine()
        self.setup_data()
        self.current_engines = filter(lambda x: x != self.engine, _engines)

    def tearDown(self):
        del self.lhses, self.rhses, self.scalar_rhses, self.scalar_lhses
        del self.pandas_rhses, self.pandas_lhses, self.current_engines, self.ne
        del self.engine

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

    @slow
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

    @slow
    def test_chained_cmp_op(self):
        mids = self.lhses
        cmp_ops = tuple(set(self.cmp_ops) - set(['==', '!=', '<=', '>=']))
        for lhs, cmp1, mid, cmp2, rhs in product(self.lhses, self.cmp_ops,
                                                 mids, cmp_ops, self.rhses):
            self.check_chained_cmp_op(lhs, cmp1, mid, cmp2, rhs)

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

    def check_chained_cmp_op(self, lhs, cmp1, mid, cmp2, rhs):
        # these are not compatible operands
        if _series_and_2d_ndarray(lhs, mid):
            self.assertRaises(ValueError, _eval_single_bin, lhs, cmp2, mid,
                              self.engine)
        else:
            lhs_new = _eval_single_bin(lhs, cmp1, mid, self.engine)

        if _series_and_2d_ndarray(mid, rhs):
            self.assertRaises(ValueError, _eval_single_bin, mid, cmp2, rhs,
                              self.engine)
        else:
            rhs_new = _eval_single_bin(mid, cmp2, rhs, self.engine)

        try:
            lhs_new
            rhs_new
        except NameError:
            pass
        else:
            # these are not compatible operands
            if (com.is_series(lhs_new) and com.is_frame(rhs_new) or
                _bool_and_frame(lhs_new, rhs_new)):
                self.assertRaises(TypeError, _eval_single_bin, lhs_new, '&',
                                  rhs_new, self.engine)
            elif _series_and_2d_ndarray(lhs_new, rhs_new):
                # TODO: once #4319 is fixed add this test back in
                #self.assertRaises(Exception, _eval_single_bin, lhs_new, '&',
                                  #rhs_new, self.engine)
                pass
            else:
                ex1 = 'lhs {0} mid {1} rhs'.format(cmp1, cmp2)
                ex2 = 'lhs {0} mid and mid {1} rhs'.format(cmp1, cmp2)
                ex3 = '(lhs {0} mid) & (mid {1} rhs)'.format(cmp1, cmp2)
                expected = _eval_single_bin(lhs_new, '&', rhs_new, self.engine)

                for ex in (ex1, ex2, ex3):
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


f = lambda *args, **kwargs: np.random.randn()


class TestAlignment(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.index_types = 'i', 'f', 's', 'u', 'dt', # 'p'

    @classmethod
    def tearDownClass(cls):
        del cls.index_types

    def check_align_nested_unary_op(self, engine):
        skip_numexpr_engine(engine)
        s = 'df * ~2'
        df = mkdf(10, 10, data_gen_f=f)
        res = pd.eval(s, engine=engine)
        assert_frame_equal(res, df * ~2)

    def test_align_nested_unary_op(self):
        for engine in _engines:
            self.check_align_nested_unary_op(engine)

    def check_basic_frame_alignment(self, engine, r_idx_type, c_idx_type):
        df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                  c_idx_type=c_idx_type)
        df2 = mkdf(20, 10, data_gen_f=f, r_idx_type=r_idx_type,
                   c_idx_type=c_idx_type)
        res = pd.eval('df + df2', engine=engine)
        assert_frame_equal(res, df + df2)

    @slow
    def test_basic_frame_alignment(self):
        args = product(_engines, self.index_types, self.index_types)
        for engine, r, c in args:
            self.check_basic_frame_alignment(engine, r, c)

    def check_frame_comparison(self, engine, r_idx_type, c_idx_type):
        df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                  c_idx_type=c_idx_type)
        res = pd.eval('df < 2', engine=engine)
        assert_frame_equal(res, df < 2)

        df3 = DataFrame(randn(*df.shape), index=df.index, columns=df.columns)
        res = pd.eval('df < df3', engine=engine)
        assert_frame_equal(res, df < df3)

    @slow
    def test_frame_comparison(self):
        args = product(_engines, self.index_types, self.index_types)
        for engine, r, c in args:
            self.check_frame_comparison(engine, r, c)

    def check_medium_complex_frame_alignment(self, engine, r1, r2, c1, c2):
        skip_numexpr_engine(engine)
        df = mkdf(5, 2, data_gen_f=f, r_idx_type=r1, c_idx_type=c1)
        df2 = mkdf(10, 2, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
        df3 = mkdf(15, 2, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
        res = pd.eval('df + df2 + df3', engine=engine)
        assert_frame_equal(res, df + df2 + df3)

    @slow
    def test_medium_complex_frame_alignment(self):
        args = product(_engines, *([self.index_types] * 4))
        for engine, r1, r2, c1, c2 in args:
            self.check_medium_complex_frame_alignment(engine, r1, r2, c1, c2)

    def check_basic_frame_series_alignment(self, engine, r_idx_type,
                                           c_idx_type, index_name):
        def testit():
            skip_numexpr_engine(engine)
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                    c_idx_type=c_idx_type)
            index = getattr(df, index_name)
            s = Series(np.random.randn(5), index[:5])

            res = pd.eval('df + s', engine=engine)
            if r_idx_type == 'dt' or c_idx_type == 'dt':
                if engine == 'numexpr':
                    expected = df.add(s)
                else:
                    expected = df + s
            else:
                expected = df + s
            assert_frame_equal(res, expected)

        testit()

    @slow
    def test_basic_frame_series_alignment(self):
        args = product(_engines, self.index_types, self.index_types,
                       ('index', 'columns'))
        for engine, r_idx_type, c_idx_type, index_name in args:
            self.check_basic_frame_series_alignment(engine, r_idx_type,
                                                    c_idx_type, index_name)

    def check_basic_series_frame_alignment(self, engine, r_idx_type,
                                           c_idx_type, index_name):
        def testit():
            skip_numexpr_engine(engine)
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                    c_idx_type=c_idx_type)
            index = getattr(df, index_name)
            s = Series(np.random.randn(5), index[:5])

            res = pd.eval('s + df', engine=engine)
            if r_idx_type == 'dt' or c_idx_type == 'dt':
                if engine == 'numexpr':
                    expected = df.add(s)
                else:
                    expected = s + df
            else:
                expected = s + df
            assert_frame_equal(res, expected)

        testit()

    @slow
    def test_basic_series_frame_alignment(self):
        args = product(_engines, self.index_types, self.index_types,
                       ('index', 'columns'))
        for engine, r_idx_type, c_idx_type, index_name in args:
            self.check_basic_series_frame_alignment(engine, r_idx_type,
                                                    c_idx_type, index_name)

    def check_series_frame_commutativity(self, engine, r_idx_type, c_idx_type,
                                         op, index_name):
        skip_numexpr_engine(engine)
        df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                  c_idx_type=c_idx_type)
        index = getattr(df, index_name)
        s = Series(np.random.randn(5), index[:5])

        lhs = 's {0} df'.format(op)
        rhs = 'df {0} s'.format(op)
        a = pd.eval(lhs, engine=engine)
        b = pd.eval(rhs, engine=engine)

        if r_idx_type != 'dt' and c_idx_type != 'dt':
            if engine == 'numexpr':
                assert_frame_equal(a, b)

    @slow
    def test_series_frame_commutativity(self):
        args = product(_engines, self.index_types, self.index_types, ('+',
                                                                      '*'),
                       ('index', 'columns'))
        for engine, r_idx_type, c_idx_type, op, index_name in args:
            self.check_series_frame_commutativity(engine, r_idx_type,
                                                  c_idx_type, op, index_name)

    def check_complex_series_frame_alignment(self, engine, index_name, obj, r1,
                                             r2, c1, c2):
        skip_numexpr_engine(engine)
        df = mkdf(10, 5, data_gen_f=f, r_idx_type=r1, c_idx_type=c1)
        df2 = mkdf(20, 5, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
        index = getattr(locals()[obj], index_name)
        s = Series(np.random.randn(5), index[:5])

        if r2 == 'dt' or c2 == 'dt':
            if engine == 'numexpr':
                expected2 = df2.add(s)
            else:
                expected2 = df2 + s
        else:
            expected2 = df2 + s

        if r1 == 'dt' or c1 == 'dt':
            if engine == 'numexpr':
                expected = expected2.add(df)
            else:
                expected = expected2 + df
        else:
            expected = expected2 + df

        res = pd.eval('df2 + s + df', engine=engine)
        self.assertEqual(res.shape, expected.shape)
        assert_frame_equal(res, expected)

    @slow
    def test_complex_series_frame_alignment(self):
        index_types = [self.index_types] * 4
        args = product(_engines, ('index', 'columns'), ('df', 'df2'),
                       *index_types)
        for engine, index_name, obj, r1, r2, c1, c2 in args:
            self.check_complex_series_frame_alignment(engine, index_name, obj,
                                                      r1, r2, c1, c2)

    @slow
    def test_performance_warning_for_asenine_alignment(self):
        df = DataFrame(randn(1000, 10))
        s = Series(randn(10000))
        with assert_produces_warning(pd.io.common.PerformanceWarning):
            pd.eval('df + s')

        s = Series(randn(1000))
        with assert_produces_warning(False):
            pd.eval('df + s')

        df = DataFrame(randn(10, 10000))
        s = Series(randn(10000))
        with assert_produces_warning(False):
            pd.eval('df + s')


class TestOperations(unittest.TestCase):

    def check_simple_arith_ops(self, engine):
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

    def check_simple_bool_ops(self, engine):
        for op, lhs, rhs in product(expr._bool_ops_syms, (True, False), (True,
                                                                        False)):
            expec = _eval_single_bin(lhs, op, rhs, engine)
            x = pd.eval('lhs {0} rhs'.format(op), engine=engine)
            assert_equal(x, expec)

    def check_bool_ops_with_constants(self, engine):
        asteval = ast.literal_eval
        for op, lhs, rhs in product(expr._bool_ops_syms, ('True', 'False'),
                                    ('True', 'False')):
            expec = _eval_single_bin(asteval(lhs), op, asteval(rhs), engine)
            x = pd.eval('{0} {1} {2}'.format(lhs, op, rhs), engine=engine)
            assert_equal(x, expec)

    def test_simple_arith_ops(self):
        for engine in _engines:
            self.check_simple_arith_ops(engine)

    def test_simple_bool_ops(self):
        for engine in _engines:
            self.check_simple_bool_ops(engine)

    def test_bool_ops_with_constants(self):
        for engine in _engines:
            self.check_bool_ops_with_constants(engine)

    def check_panel_fails(self, engine):
        x = Panel(randn(3, 4, 5))
        y = Series(randn(10))
        assert_raises(NotImplementedError, pd.eval, 'x + y',
                      local_dict={'x': x, 'y': y}, engine=engine)

    def test_panel_fails(self):
        for engine in _engines:
            self.check_panel_fails(engine)

    def check_4d_ndarray_fails(self, engine):
        x = randn(3, 4, 5, 6)
        y = Series(randn(10))
        assert_raises(NotImplementedError, pd.eval, 'x + y', local_dict={'x': x,
                                                                        'y': y},
                    engine=engine)

    def test_4d_ndarray_fails(self):
        for engine in _engines:
            self.check_4d_ndarray_fails(engine)

    def check_constant(self, engine):
        x = pd.eval('1', engine=engine)
        assert_equal(x, 1)

    def test_constant(self):
        for engine in _engines:
            self.check_constant(engine)

    def check_single_variable(self, engine):
        df = DataFrame(randn(10, 2))
        df2 = pd.eval('df', engine=engine)
        assert_frame_equal(df, df2)

    def test_single_variable(self):
        for engine in _engines:
            self.check_single_variable(engine)

    def test_truediv(self):
        for engine in _engines:
            self.check_truediv(engine)

    def check_truediv(self, engine):
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

    def test_python_fails_and(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(NotImplementedError, pd.eval, 'df > 2 and df > 3',
                          local_dict={'df': df}, parser='python')

    def test_python_fails_or(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(NotImplementedError, pd.eval, 'df > 2 or df > 3',
                          local_dict={'df': df}, parser='python')

    def test_python_fails_not(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(NotImplementedError, pd.eval, 'not df > 2',
                          local_dict={'df': df}, parser='python')

    def test_python_fails_ampersand(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(TypeError, pd.eval,
                          '(df + 2)[df > 1] > 0 & (df > 0)',
                          local_dict={'df': df}, parser='python')

    def test_python_fails_pipe(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(TypeError, pd.eval,
                          '(df + 2)[df > 1] > 0 | (df > 0)',
                          local_dict={'df': df}, parser='python')

    def check_failing_subscript_with_name_error(self, engine):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(NameError, pd.eval, 'df[x > 2] > 2',
                          local_dict={'df': df}, engine=engine)

    def test_failing_subscript_with_name_error(self):
        for engine in _engines:
            self.check_failing_subscript_with_name_error(engine)

    def check_lhs_expression_subscript(self, engine):
        df = DataFrame(np.random.randn(5, 3))
        result = pd.eval('(df + 1)[df > 2]', engine=engine)
        expected = (df + 1)[df > 2]
        assert_frame_equal(result, expected)

    def test_lhs_expression_subscript(self):
        for engine in _engines:
            self.check_lhs_expression_subscript(engine)

    def check_attr_expression(self, engine):
        df = DataFrame(np.random.randn(5, 3), columns=list('abc'))
        expr1 = 'df.a < df.b'
        expec1 = df.a < df.b
        expr2 = 'df.a + df.b + df.c'
        expec2 = df.a + df.b + df.c
        expr3 = 'df.a + df.b + df.c[df.b < 0]'
        expec3 = df.a + df.b + df.c[df.b < 0]
        exprs = expr1, expr2, expr3
        expecs = expec1, expec2, expec3
        for e, expec in zip(exprs, expecs):
            assert_series_equal(expec, pd.eval(e, engine=engine))

    def test_attr_expression(self):
        for engine in _engines:
            self.check_attr_expression(engine)

    def check_assignment_fails(self, engine, parser):
        df = DataFrame(np.random.randn(5, 3), columns=list('abc'))
        df2 = DataFrame(np.random.randn(5, 3))
        expr1 = 'df = df2'
        self.assertRaises(NotImplementedError, pd.eval, expr1,
                          local_dict={'df': df, 'df2': df2}, engine=engine,
                          parser=parser)

    def test_assignment_fails(self):
        for engine, parser in product(_engines.iterkeys(), ('pandas',
                                                            'python')):
            self.check_assignment_fails(engine, parser)

    def check_basic_period_index_boolean_expression(self, engine):
        df = mkdf(2, 2, data_gen_f=f, c_idx_type='p', r_idx_type='i')

        e = df < 2
        r = pd.eval('df < 2', engine=engine)
        x = df < 2

        assert_frame_equal(r, e)
        assert_frame_equal(x, e)

    def test_basic_period_index_expression_python(self):
        for engine in _engines:
            self.check_basic_period_index_boolean_expression(engine)

    def check_basic_period_index_subscript_expression(self, engine):
        df = mkdf(2, 2, data_gen_f=f, c_idx_type='p', r_idx_type='i')
        r = pd.eval('df[df < 2 + 3]', engine=engine)
        e = df[df < 2 + 3]
        assert_frame_equal(r, e)

    def test_basic_period_index_subscript_expression(self):
        for engine in _engines:
            self.check_basic_period_index_subscript_expression(engine)

    def check_nested_period_index_subscript_expression(self, engine):
        df = mkdf(2, 2, data_gen_f=f, c_idx_type='p', r_idx_type='i')
        r = pd.eval('df[df[df < 2] < 2] + df * 2', engine=engine)
        e = df[df[df < 2] < 2] + df * 2
        assert_frame_equal(r, e)

    def test_nested_period_index_subscript_expression(self):
        for engine in _engines:
            self.check_nested_period_index_subscript_expression(engine)

    def test_simple_not_expression(self):
        df = DataFrame(randn(10, 3), columns=list('abc'))
        df['bools'] = rand(len(df)) > 0.5
        res = df['not bools']
        res2 = df['~bools']
        expec = df[~df.bools]
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

    def test_complex_boolean_expression(self):
        df = DataFrame(randn(10, 3), columns=list('abc'))
        df['bools'] = rand(len(df)) > 0.5
        res = df['a < b < c and (not bools) or bools > 2']
        expec = df[(df.a < df.b) & (df.b < df.c) & (~df.bools) | (df.bools > 2)]
        assert_frame_equal(res, expec)

_var_s = randn(10)

class TestScope(unittest.TestCase):

    def check_global_scope(self, engine):
        e = '_var_s * 2'
        assert_array_equal(_var_s * 2, pd.eval(e, engine=engine))

    def test_global_scope(self):
        for engine in _engines:
            self.check_global_scope(engine)

    def check_no_new_locals(self, engine):
        x = 1
        lcls = locals().copy()
        pd.eval('x + 1', local_dict=lcls)
        lcls2 = locals().copy()
        lcls2.pop('lcls')
        assert_equal(lcls, lcls2)

    def test_no_new_locals(self):
        for engine in _engines:
            self.check_no_new_locals(engine)

    def check_no_new_globals(self, engine):
        x = 1
        gbls = globals().copy()
        pd.eval('x + 1')
        gbls2 = globals().copy()
        assert_equal(gbls, gbls2)

    def test_no_new_globals(self):
        for engine in _engines:
            self.check_no_new_globals(engine)

    def check_nested_scope(self, engine):
        # smoke test
        x = 1
        result = pd.eval('x + 1', engine=engine)
        self.assertEqual(result, 2)

        df  = DataFrame(np.random.randn(5, 3))
        df2 = DataFrame(np.random.randn(5, 3))
        expected = df[(df>0) & (df2>0)]

        result = df['(df>0) & (df2>0)']
        assert_frame_equal(result, expected)

        result = df.query('(df>0) & (df2>0)', engine=engine)
        assert_frame_equal(result, expected)

        result = pd.eval('df[(df > 0) and (df2 > 0)]', engine=engine)
        assert_frame_equal(result, expected)

        result = pd.eval('df[(df > 0) and (df2 > 0) and df[df > 0] > 0]', engine=engine)
        expected = df[(df > 0) & (df2 > 0) & (df[df > 0] > 0)]
        assert_frame_equal(result, expected)

        result = pd.eval('df[(df>0) & (df2>0)]',engine=engine)
        expected = df.query('(df>0) & (df2>0)', engine=engine)
        assert_frame_equal(result, expected)

    def test_nested_scope(self):
        for engine in _engines:
            self.check_nested_scope(engine)


def test_invalid_engine():
    assertRaisesRegexp(KeyError, 'Invalid engine \'asdf\' passed',
                       pd.eval, 'x + y', local_dict={'x': 1, 'y': 2},
                       engine='asdf')


def test_invalid_parser():
    assertRaisesRegexp(KeyError, 'Invalid parser \'asdf\' passed',
                       pd.eval, 'x + y', local_dict={'x': 1, 'y': 2},
                       parser='asdf')


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


def test_syntax_error_exprs():
    for engine in _engines:
        e = 's +'
        assert_raises(SyntaxError, pd.eval, e, engine=engine)


def test_name_error_exprs():
    for engine in _engines:
        e = 's + t'
        assert_raises(NameError, pd.eval, e, engine=engine)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
