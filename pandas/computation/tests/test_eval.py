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
from pandas import DataFrame, Series, Panel, date_range
from pandas.util.testing import makeCustomDataframe as mkdf

from pandas.computation import pytables
from pandas.computation.expressions import _USE_NUMEXPR
from pandas.computation.engines import _engines
from pandas.computation.expr import PythonExprVisitor, PandasExprVisitor
from pandas.computation.ops import (_binary_ops_dict, _unary_ops_dict,
                                    _special_case_arith_ops_syms,
                                    _arith_ops_syms, _bool_ops_syms)
import pandas.computation.expr as expr
from pandas.util.testing import (assert_frame_equal, randbool,
                                 assertRaisesRegexp,
                                 assert_produces_warning, assert_series_equal)
from pandas.compat import PY3, u

_series_frame_incompatible = _bool_ops_syms
_scalar_skip = 'in', 'not in'

def skip_if_no_ne(engine='numexpr'):
    if not _USE_NUMEXPR and engine == 'numexpr':
        raise nose.SkipTest("numexpr engine not installed or disabled")


def engine_has_neg_frac(engine):
    return _engines[engine].has_neg_frac


def _eval_single_bin(lhs, cmp1, rhs, engine):
    c = _binary_ops_dict[cmp1]
    if engine_has_neg_frac(engine):
        try:
            return c(lhs, rhs)
        except ValueError as e:
            try:
                msg = e.message
            except AttributeError:
                msg = e
            msg = u(msg)
            if msg == u('negative number cannot be raised to a fractional'
                        ' power'):
                return np.nan
            raise
    return c(lhs, rhs)


def _series_and_2d_ndarray(lhs, rhs):
    return ((isinstance(lhs, Series) and
             isinstance(rhs, np.ndarray) and rhs.ndim > 1)
            or (isinstance(rhs, Series) and
                isinstance(lhs, np.ndarray) and lhs.ndim > 1))


def _series_and_frame(lhs, rhs):
    return ((isinstance(lhs, Series) and isinstance(rhs, DataFrame))
            or (isinstance(rhs, Series) and isinstance(lhs, DataFrame)))


def _bool_and_frame(lhs, rhs):
    return isinstance(lhs, bool) and isinstance(rhs, pd.core.generic.NDFrame)


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


def _is_py3_complex_incompat(result, expected):
    return (PY3 and isinstance(expected, (complex, np.complexfloating)) and
            np.isnan(result))


_good_arith_ops = com.difference(_arith_ops_syms, _special_case_arith_ops_syms)


class TestEvalNumexprPandas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        skip_if_no_ne()
        import numexpr as ne
        cls.ne = ne
        cls.engine = 'numexpr'
        cls.parser = 'pandas'

    @classmethod
    def tearDownClass(cls):
        del cls.engine, cls.parser
        if hasattr(cls, 'ne'):
            del cls.ne

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

    def setup_ops(self):
        self.cmp_ops = expr._cmp_ops_syms
        self.cmp2_ops = self.cmp_ops[::-1]
        self.bin_ops = expr._bool_ops_syms
        self.special_case_ops = _special_case_arith_ops_syms
        self.arith_ops = _good_arith_ops
        self.unary_ops = '+', '-', '~', 'not '

    def setUp(self):
        self.setup_ops()
        self.setup_data()
        self.current_engines = filter(lambda x: x != self.engine, _engines)

    def tearDown(self):
        del self.lhses, self.rhses, self.scalar_rhses, self.scalar_lhses
        del self.pandas_rhses, self.pandas_lhses, self.current_engines

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
        for lhs, cmp1, mid, cmp2, rhs in product(self.lhses, cmp_ops,
                                                 mids, cmp_ops, self.rhses):
            self.check_chained_cmp_op(lhs, cmp1, mid, cmp2, rhs)

    def check_complex_cmp_op(self, lhs, cmp1, rhs, binop, cmp2):
        skip_these = 'in', 'not in'
        ex = '(lhs {cmp1} rhs) {binop} (lhs {cmp2} rhs)'.format(cmp1=cmp1,
                                                                binop=binop,
                                                                cmp2=cmp2)
        scalar_with_in_notin = (np.isscalar(rhs) and (cmp1 in skip_these or
                                cmp2 in skip_these))
        if scalar_with_in_notin:
            self.assertRaises(TypeError, pd.eval, ex, engine=self.engine,
                              parser=self.parser, local_dict={'lhs': lhs,
                                                              'rhs': rhs})
        elif (_series_and_frame(lhs, rhs) and (cmp1 in
              _series_frame_incompatible or
              cmp2 in _series_frame_incompatible)):
            self.assertRaises(TypeError, pd.eval, ex,
                              local_dict={'lhs': lhs, 'rhs': rhs},
                              engine=self.engine, parser=self.parser)
        else:
            lhs_new = _eval_single_bin(lhs, cmp1, rhs, self.engine)
            rhs_new = _eval_single_bin(lhs, cmp2, rhs, self.engine)
            if (isinstance(lhs_new, Series) and isinstance(rhs_new, DataFrame)
                    and binop in _series_frame_incompatible):
                pass
                # TODO: the code below should be added back when left and right
                # hand side bool ops are fixed.

                #try:
                    #self.assertRaises(Exception, pd.eval, ex,
                                    #local_dict={'lhs': lhs, 'rhs': rhs},
                                    #engine=self.engine, parser=self.parser)
                #except AssertionError:
                    #import ipdb; ipdb.set_trace()
                    #raise

            else:
                expected = _eval_single_bin(lhs_new, binop, rhs_new, self.engine)
                result = pd.eval(ex, engine=self.engine, parser=self.parser)
                assert_array_equal(result, expected)

    @skip_incompatible_operand
    def check_chained_cmp_op(self, lhs, cmp1, mid, cmp2, rhs):
        skip_these = 'in', 'not in'

        def check_operands(left, right, cmp_op):
            if (np.isscalar(right) and not np.isscalar(left) and cmp_op in
                    skip_these):
                self.assertRaises(Exception, _eval_single_bin, left, cmp_op,
                                  right, self.engine)
            elif _series_and_2d_ndarray(right, left):
                self.assertRaises(Exception, _eval_single_bin, right, cmp_op,
                                  left, self.engine)
            elif (np.isscalar(right) and np.isscalar(left) and cmp_op in
                    skip_these):
                self.assertRaises(Exception, _eval_single_bin, right, cmp_op,
                                  left, self.engine)
            else:
                new = _eval_single_bin(left, cmp_op, right, self.engine)
                return new
            return

        lhs_new = check_operands(lhs, mid, cmp1)
        rhs_new = check_operands(mid, rhs, cmp2)

        if lhs_new is not None and rhs_new is not None:
            # these are not compatible operands
            if isinstance(lhs_new, Series) and isinstance(rhs_new, DataFrame):
                self.assertRaises(TypeError, _eval_single_bin, lhs_new, '&',
                                  rhs_new, self.engine)
            elif (_bool_and_frame(lhs_new, rhs_new)):
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
                try:
                    expected = _eval_single_bin(lhs_new, '&', rhs_new, self.engine)
                except TypeError:
                    import ipdb; ipdb.set_trace()
                    raise

                for ex in (ex1, ex2, ex3):
                    result = pd.eval(ex, engine=self.engine,
                                     parser=self.parser)
                    assert_array_equal(result, expected)

    @skip_incompatible_operand
    def check_simple_cmp_op(self, lhs, cmp1, rhs):
        ex = 'lhs {0} rhs'.format(cmp1)
        if cmp1 in ('in', 'not in') and not com.is_list_like(rhs):
            self.assertRaises(TypeError, pd.eval, ex, engine=self.engine,
                              parser=self.parser, local_dict={'lhs': lhs,
                                                              'rhs': rhs})
        else:
            expected = _eval_single_bin(lhs, cmp1, rhs, self.engine)
            result = pd.eval(ex, engine=self.engine, parser=self.parser)
            assert_array_equal(result, expected)

    @skip_incompatible_operand
    def check_binary_arith_op(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)
        result = pd.eval(ex, engine=self.engine, parser=self.parser)
        expected = _eval_single_bin(lhs, arith1, rhs, self.engine)
        assert_array_equal(result, expected)
        ex = 'lhs {0} rhs {0} rhs'.format(arith1)
        result = pd.eval(ex, engine=self.engine, parser=self.parser)
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
        result = pd.eval(ex, engine=self.engine, parser=self.parser)
        expected = lhs % rhs
        assert_allclose(result, expected)
        expected = self.ne.evaluate('expected {0} rhs'.format(arith1))
        assert_allclose(result, expected)

    @skip_incompatible_operand
    def check_floor_division(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)

        if self.engine == 'python':
            res = pd.eval(ex, engine=self.engine, parser=self.parser)
            expected = lhs // rhs
            assert_array_equal(res, expected)
        else:
            self.assertRaises(TypeError, pd.eval, ex, local_dict={'lhs': lhs,
                                                                  'rhs': rhs},
                              engine=self.engine, parser=self.parser)

    def get_expected_pow_result(self, lhs, rhs):
        try:
            expected = _eval_single_bin(lhs, '**', rhs, self.engine)
        except ValueError as e:
            msg = 'negative number cannot be raised to a fractional power'
            try:
                emsg = e.message
            except AttributeError:
                emsg = e

            emsg = u(emsg)

            if emsg == msg:
                if self.engine == 'python':
                    raise nose.SkipTest(emsg)
                else:
                    expected = np.nan
            else:
                raise
        return expected

    @skip_incompatible_operand
    def check_pow(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)
        expected = self.get_expected_pow_result(lhs, rhs)
        result = pd.eval(ex, engine=self.engine, parser=self.parser)

        if (np.isscalar(lhs) and np.isscalar(rhs) and
            _is_py3_complex_incompat(result, expected)):
            self.assertRaises(AssertionError, assert_array_equal, result,
                              expected)
        else:
            assert_array_equal(result, expected)

            ex = '(lhs {0} rhs) {0} rhs'.format(arith1)
            result = pd.eval(ex, engine=self.engine, parser=self.parser)
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
            result = pd.eval('~elb', engine=self.engine, parser=self.parser)
            assert_array_equal(expected, result)

            for engine in self.current_engines:
                skip_if_no_ne(engine)
                assert_array_equal(result, pd.eval('~elb', engine=engine,
                                                   parser=self.parser))

    @skip_incompatible_operand
    def check_compound_invert_op(self, lhs, cmp1, rhs):
        skip_these = 'in', 'not in'
        ex = '~(lhs {0} rhs)'.format(cmp1)

        if np.isscalar(rhs) and cmp1 in skip_these:
            self.assertRaises(TypeError, pd.eval, ex, engine=self.engine,
                              parser=self.parser, local_dict={'lhs': lhs,
                                                              'rhs': rhs})
        else:
            # compound
            if np.isscalar(lhs) and np.isscalar(rhs):
                lhs, rhs = map(lambda x: np.array([x]), (lhs, rhs))
            expected = ~_eval_single_bin(lhs, cmp1, rhs, self.engine)
            result = pd.eval(ex, engine=self.engine, parser=self.parser)
            assert_array_equal(expected, result)

            # make sure the other engines work the same as this one
            for engine in self.current_engines:
                skip_if_no_ne(engine)
                ev = pd.eval(ex, engine=self.engine, parser=self.parser)
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

        result = pd.eval(ex, engine=self.engine, parser=self.parser)
        assert_array_equal(result, expected)

        for engine in self.current_engines:
            skip_if_no_ne(engine)
            assert_array_equal(result, pd.eval(ex, engine=engine,
                                               parser=self.parser))

        ex = '{0}(lhs {1} rhs)'.format(unary_op, arith1)
        result = pd.eval(ex, engine=self.engine, parser=self.parser)


class TestEvalNumexprPython(TestEvalNumexprPandas):
    @classmethod
    def setUpClass(cls):
        skip_if_no_ne()
        import numexpr as ne
        cls.ne = ne
        cls.engine = 'numexpr'
        cls.parser = 'python'

    def setup_ops(self):
        self.cmp_ops = expr._cmp_ops_syms
        self.cmp2_ops = self.cmp_ops[::-1]
        self.bin_ops = (s for s in expr._bool_ops_syms if s not in ('and', 'or'))
        self.special_case_ops = _special_case_arith_ops_syms
        self.arith_ops = _good_arith_ops
        self.unary_ops = '+', '-', '~'

    def check_chained_cmp_op(self, lhs, cmp1, mid, cmp2, rhs):
        ex1 = 'lhs {0} mid {1} rhs'.format(cmp1, cmp2)
        self.assertRaises(NotImplementedError, pd.eval, ex1,
                          local_dict={'lhs': lhs, 'mid': mid, 'rhs': rhs},
                          engine=self.engine, parser=self.parser)


class TestEvalPythonPython(TestEvalNumexprPython):
    @classmethod
    def setUpClass(cls):
        cls.engine = 'python'
        cls.parser = 'python'

    @skip_incompatible_operand
    def check_modulus(self, lhs, arith1, rhs):
        ex = 'lhs {0} rhs'.format(arith1)
        result = pd.eval(ex, engine=self.engine)
        expected = lhs % rhs
        assert_allclose(result, expected)
        expected = eval('expected {0} rhs'.format(arith1))
        assert_allclose(result, expected)

    def check_alignment(self, result, nlhs, ghs, op):
        try:
            nlhs, ghs = nlhs.align(ghs)
        except (ValueError, TypeError, AttributeError):
            # ValueError: series frame or frame series align
            # TypeError, AttributeError: series or frame with scalar align
            pass
        else:
            expected = eval('nlhs {0} ghs'.format(op))
            assert_array_equal(result, expected)


class TestEvalPythonPandas(TestEvalPythonPython):
    @classmethod
    def setUpClass(cls):
        cls.engine = 'python'
        cls.parser = 'pandas'

    def check_chained_cmp_op(self, lhs, cmp1, mid, cmp2, rhs):
        TestEvalNumexprPandas.check_chained_cmp_op(self, lhs, cmp1, mid, cmp2,
                                                   rhs)


f = lambda *args, **kwargs: np.random.randn()


ENGINES_PARSERS = list(product(_engines, expr._parsers))


#-------------------------------------
# basic and complex alignment

class TestAlignment(object):

    index_types = 'i', 'f', 's', 'u', 'dt', # 'p'

    def check_align_nested_unary_op(self, engine, parser):
        skip_if_no_ne(engine)
        s = 'df * ~2'
        df = mkdf(5, 3, data_gen_f=f)
        res = pd.eval(s, engine=engine, parser=parser)
        assert_frame_equal(res, df * ~2)

    def test_align_nested_unary_op(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_align_nested_unary_op, engine, parser

    def check_basic_frame_alignment(self, engine, parser):
        skip_if_no_ne(engine)
        args = product(self.index_types, repeat=2)
        for r_idx_type, c_idx_type in args:
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                    c_idx_type=c_idx_type)
            df2 = mkdf(20, 10, data_gen_f=f, r_idx_type=r_idx_type,
                    c_idx_type=c_idx_type)
            res = pd.eval('df + df2', engine=engine, parser=parser)
            assert_frame_equal(res, df + df2)

    @slow
    def test_basic_frame_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_basic_frame_alignment, engine, parser

    def check_frame_comparison(self, engine, parser):
        skip_if_no_ne(engine)
        args = product(self.index_types, repeat=2)
        for r_idx_type, c_idx_type in args:
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                    c_idx_type=c_idx_type)
            res = pd.eval('df < 2', engine=engine, parser=parser)
            assert_frame_equal(res, df < 2)

            df3 = DataFrame(randn(*df.shape), index=df.index,
                            columns=df.columns)
            res = pd.eval('df < df3', engine=engine, parser=parser)
            assert_frame_equal(res, df < df3)

    @slow
    def test_frame_comparison(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_frame_comparison, engine, parser

    def check_medium_complex_frame_alignment(self, engine, parser):
        skip_if_no_ne(engine)
        args = product(self.index_types, repeat=4)
        for r1, c1, r2, c2 in args:
            df = mkdf(5, 2, data_gen_f=f, r_idx_type=r1, c_idx_type=c1)
            df2 = mkdf(10, 2, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
            df3 = mkdf(15, 2, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
            res = pd.eval('df + df2 + df3', engine=engine, parser=parser)
            assert_frame_equal(res, df + df2 + df3)

    @slow
    def test_medium_complex_frame_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_medium_complex_frame_alignment, engine, parser

    def check_basic_frame_series_alignment(self, engine, parser):
        skip_if_no_ne(engine)
        def testit(r_idx_type, c_idx_type, index_name):
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                    c_idx_type=c_idx_type)
            index = getattr(df, index_name)
            s = Series(np.random.randn(5), index[:5])

            res = pd.eval('df + s', engine=engine, parser=parser)
            if r_idx_type == 'dt' or c_idx_type == 'dt':
                if engine == 'numexpr':
                    expected = df.add(s)
                else:
                    expected = df + s
            else:
                expected = df + s
            assert_frame_equal(res, expected)

        args = product(self.index_types, self.index_types, ('index',
                                                            'columns'))
        for r_idx_type, c_idx_type, index_name in args:
            testit(r_idx_type, c_idx_type, index_name)

    @slow
    def test_basic_frame_series_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_basic_frame_series_alignment, engine, parser

    def check_basic_series_frame_alignment(self, engine, parser):
        skip_if_no_ne(engine)
        def testit(r_idx_type, c_idx_type, index_name):
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                    c_idx_type=c_idx_type)
            index = getattr(df, index_name)
            s = Series(np.random.randn(5), index[:5])

            res = pd.eval('s + df', engine=engine, parser=parser)
            if r_idx_type == 'dt' or c_idx_type == 'dt':
                if engine == 'numexpr':
                    expected = df.add(s)
                else:
                    expected = s + df
            else:
                expected = s + df
            assert_frame_equal(res, expected)

        args = product(self.index_types, self.index_types, ('index',
                                                            'columns'))
        for r_idx_type, c_idx_type, index_name in args:
            testit(r_idx_type, c_idx_type, index_name)

    @slow
    def test_basic_series_frame_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_basic_series_frame_alignment, engine, parser

    def check_series_frame_commutativity(self, engine, parser):
        skip_if_no_ne(engine)
        args = product(self.index_types, self.index_types, ('+', '*'),
                       ('index', 'columns'))
        for r_idx_type, c_idx_type, op, index_name in args:
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                    c_idx_type=c_idx_type)
            index = getattr(df, index_name)
            s = Series(np.random.randn(5), index[:5])

            lhs = 's {0} df'.format(op)
            rhs = 'df {0} s'.format(op)
            a = pd.eval(lhs, engine=engine, parser=parser)
            b = pd.eval(rhs, engine=engine, parser=parser)

            if r_idx_type != 'dt' and c_idx_type != 'dt':
                if engine == 'numexpr':
                    assert_frame_equal(a, b)

    @slow
    def test_series_frame_commutativity(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_series_frame_commutativity, engine, parser

    def check_complex_series_frame_alignment(self, engine, parser):
        skip_if_no_ne(engine)
        index_types = [self.index_types] * 4
        args = product(('index', 'columns'), ('df', 'df2'), *index_types)
        for index_name, obj, r1, r2, c1, c2 in args:
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

            res = pd.eval('df2 + s + df', engine=engine, parser=parser)
            assert_equal(res.shape, expected.shape)
            assert_frame_equal(res, expected)

    @slow
    def test_complex_series_frame_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_complex_series_frame_alignment, engine, parser

    def check_performance_warning_for_asenine_alignment(self, engine, parser):
        skip_if_no_ne(engine)
        df = DataFrame(randn(1000, 10))
        s = Series(randn(10000))
        if engine == 'numexpr':
            seen = pd.io.common.PerformanceWarning
        else:
            seen = False

        with assert_produces_warning(seen):
            pd.eval('df + s', engine=engine, parser=parser)

        s = Series(randn(1000))
        with assert_produces_warning(False):
            pd.eval('df + s', engine=engine, parser=parser)

        df = DataFrame(randn(10, 10000))
        s = Series(randn(10000))
        with assert_produces_warning(False):
            pd.eval('df + s', engine=engine, parser=parser)

    def test_performance_warning_for_asenine_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_performance_warning_for_asenine_alignment, engine, parser


#------------------------------------
# slightly more complex ops

class TestOperationsNumExprPandas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        skip_if_no_ne()
        cls.engine = 'numexpr'
        cls.parser = 'pandas'

    @classmethod
    def tearDownClass(cls):
        del cls.engine, cls.parser

    def eval(self, *args, **kwargs):
        kwargs['engine'] = self.engine
        kwargs['parser'] = self.parser
        return pd.eval(*args, **kwargs)

    def test_simple_arith_ops(self):
        ops = expr._arith_ops_syms + expr._cmp_ops_syms

        for op in filter(lambda x: x != '//', ops):
            ex = '1 {0} 1'.format(op)
            ex2 = 'x {0} 1'.format(op)
            ex3 = '1 {0} (x + 1)'.format(op)

            if op in ('in', 'not in'):
                self.assertRaises(TypeError, pd.eval, ex,
                                  engine=self.engine, parser=self.parser)
            else:
                expec = _eval_single_bin(1, op, 1, self.engine)
                x = self.eval(ex, engine=self.engine, parser=self.parser)
                assert_equal(x, expec)

                expec = _eval_single_bin(x, op, 1, self.engine)
                y = self.eval(ex2, local_dict={'x': x}, engine=self.engine,
                              parser=self.parser)
                assert_equal(y, expec)

                expec = _eval_single_bin(1, op, x + 1, self.engine)
                y = self.eval(ex3, local_dict={'x': x},
                              engine=self.engine, parser=self.parser)
                assert_equal(y, expec)

    def test_simple_bool_ops(self):
        for op, lhs, rhs in product(expr._bool_ops_syms, (True, False),
                                    (True, False)):
            expec = _eval_single_bin(lhs, op, rhs, self.engine)
            x = self.eval('lhs {0} rhs'.format(op), local_dict={'lhs': lhs,
                                                                'rhs': rhs},
                          engine=self.engine, parser=self.parser)
            assert_equal(x, expec)

    def test_bool_ops_with_constants(self):
        asteval = ast.literal_eval
        for op, lhs, rhs in product(expr._bool_ops_syms, ('True', 'False'),
                                    ('True', 'False')):
            expec = _eval_single_bin(asteval(lhs), op, asteval(rhs),
                                     self.engine)
            x = self.eval('{0} {1} {2}'.format(lhs, op, rhs),
                          local_dict={'lhs': lhs, 'rhs': rhs})
            assert_equal(x, expec)

    def test_panel_fails(self):
        x = Panel(randn(3, 4, 5))
        y = Series(randn(10))
        assert_raises(NotImplementedError, self.eval, 'x + y',
                      local_dict={'x': x, 'y': y})

    def test_4d_ndarray_fails(self):
        x = randn(3, 4, 5, 6)
        y = Series(randn(10))
        assert_raises(NotImplementedError, self.eval, 'x + y',
                      local_dict={'x': x, 'y': y})

    def test_constant(self):
        x = self.eval('1')
        assert_equal(x, 1)

    def test_single_variable(self):
        df = DataFrame(randn(10, 2))
        df2 = self.eval('df', local_dict={'df': df})
        assert_frame_equal(df, df2)

    def test_truediv(self):
        s = np.array([1])
        ex = 's / 1'
        d = {'s': s}

        if PY3:
            res = self.eval(ex, truediv=False, local_dict=d)
            assert_array_equal(res, np.array([1.0]))

            res = self.eval(ex, truediv=True, local_dict=d)
            assert_array_equal(res, np.array([1.0]))

            res = self.eval('1 / 2', truediv=True)
            expec = 0.5
            self.assertEqual(res, expec)

            res = self.eval('1 / 2', truediv=False)
            expec = 0.5
            self.assertEqual(res, expec)

            res = self.eval('s / 2', truediv=False, local_dict={'s': s})
            expec = 0.5
            self.assertEqual(res, expec)

            res = self.eval('s / 2', truediv=True, local_dict={'s': s})
            expec = 0.5
            self.assertEqual(res, expec)
        else:
            res = self.eval(ex, truediv=False, local_dict=d)
            assert_array_equal(res, np.array([1]))

            res = self.eval(ex, truediv=True, local_dict=d)
            assert_array_equal(res, np.array([1.0]))

            res = self.eval('1 / 2', truediv=True)
            expec = 0.5
            self.assertEqual(res, expec)

            res = self.eval('1 / 2', truediv=False)
            expec = 0
            self.assertEqual(res, expec)

            res = self.eval('s / 2', truediv=False, local_dict={'s': s})
            expec = 0
            self.assertEqual(res, expec)

            res = self.eval('s / 2', truediv=True, local_dict={'s': s})
            expec = 0.5
            self.assertEqual(res, expec)

    def test_failing_subscript_with_name_error(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(NameError, self.eval, 'df[x > 2] > 2',
                          local_dict={'df': df})

    def test_lhs_expression_subscript(self):
        df = DataFrame(np.random.randn(5, 3))
        result = self.eval('(df + 1)[df > 2]', local_dict={'df': df})
        expected = (df + 1)[df > 2]
        assert_frame_equal(result, expected)

    def test_attr_expression(self):
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
            assert_series_equal(expec, self.eval(e, local_dict={'df': df}))

    def test_assignment_fails(self):
        df = DataFrame(np.random.randn(5, 3), columns=list('abc'))
        df2 = DataFrame(np.random.randn(5, 3))
        expr1 = 'df = df2'
        self.assertRaises(NotImplementedError, self.eval, expr1,
                          local_dict={'df': df, 'df2': df2})

    def test_basic_period_index_boolean_expression(self):
        df = mkdf(2, 2, data_gen_f=f, c_idx_type='p', r_idx_type='i')

        e = df < 2
        r = self.eval('df < 2', local_dict={'df': df})
        x = df < 2

        assert_frame_equal(r, e)
        assert_frame_equal(x, e)

    def test_basic_period_index_subscript_expression(self):
        df = mkdf(2, 2, data_gen_f=f, c_idx_type='p', r_idx_type='i')
        r = self.eval('df[df < 2 + 3]', local_dict={'df': df})
        e = df[df < 2 + 3]
        assert_frame_equal(r, e)

    def test_nested_period_index_subscript_expression(self):
        df = mkdf(2, 2, data_gen_f=f, c_idx_type='p', r_idx_type='i')
        r = self.eval('df[df[df < 2] < 2] + df * 2', local_dict={'df': df})
        e = df[df[df < 2] < 2] + df * 2
        assert_frame_equal(r, e)

    def test_date_boolean(self):
        df = DataFrame(randn(5, 3))
        df['dates1'] = date_range('1/1/2012', periods=5)
        res = self.eval('df.dates1 < 20130101', local_dict={'df': df})
        expec = df.dates1 < '20130101'
        assert_series_equal(res, expec)


class TestOperationsNumExprPython(TestOperationsNumExprPandas):
    @classmethod
    def setUpClass(cls):
        if not _USE_NUMEXPR:
            raise nose.SkipTest("numexpr engine not installed")
        cls.engine = 'numexpr'
        cls.parser = 'python'

    def test_fails_and(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(NotImplementedError, pd.eval, 'df > 2 and df > 3',
                          local_dict={'df': df}, parser=self.parser,
                          engine=self.engine)

    def test_fails_or(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(NotImplementedError, pd.eval, 'df > 2 or df > 3',
                          local_dict={'df': df}, parser=self.parser,
                          engine=self.engine)

    def test_fails_not(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(NotImplementedError, pd.eval, 'not df > 2',
                          local_dict={'df': df}, parser=self.parser,
                          engine=self.engine)

    def test_fails_ampersand(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(TypeError, pd.eval,
                          '(df + 2)[df > 1] > 0 & (df > 0)',
                          local_dict={'df': df}, parser=self.parser,
                          engine=self.engine)

    def test_fails_pipe(self):
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(TypeError, pd.eval,
                          '(df + 2)[df > 1] > 0 | (df > 0)',
                          local_dict={'df': df}, parser=self.parser,
                          engine=self.engine)

    def test_bool_ops_with_constants(self):
        from ast import literal_eval as asteval
        for op, lhs, rhs in product(expr._bool_ops_syms, ('True', 'False'),
                                    ('True', 'False')):
            if op not in ('and', 'or'):
                expec = _eval_single_bin(asteval(lhs), op, asteval(rhs),
                                        self.engine)
                x = self.eval('{0} {1} {2}'.format(lhs, op, rhs),
                            local_dict={'lhs': lhs, 'rhs': rhs})
                assert_equal(x, expec)
            else:
                self.assertRaises(NotImplementedError,
                                  self.eval,
                                  '{0} {1} {2}'.format(lhs, op, rhs),
                                  local_dict={'lhs': lhs, 'rhs': rhs})

    def test_simple_bool_ops(self):
        for op, lhs, rhs in product(expr._bool_ops_syms, (True, False), (True,
                                                                        False)):
            if op not in ('and', 'or'):
                expec = _eval_single_bin(lhs, op, rhs, self.engine)
                x = self.eval('lhs {0} rhs'.format(op), local_dict={'lhs': lhs,
                                                                    'rhs': rhs})
                assert_equal(x, expec)
            else:
                self.assertRaises(NotImplementedError,
                                  self.eval,
                                  'lhs {0} rhs'.format(op),
                                  local_dict={'lhs': lhs, 'rhs': rhs})


class TestOperationsPythonPython(TestOperationsNumExprPython):
    @classmethod
    def setUpClass(cls):
        cls.engine = cls.parser = 'python'

    def test_fails_ampersand(self):
        raise nose.SkipTest("known failer for now")
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(TypeError, pd.eval,
                          '(df + 2)[df > 1] > 0 & (df > 0)',
                          local_dict={'df': df}, parser=self.parser,
                          engine=self.engine)

    def test_fails_pipe(self):
        raise nose.SkipTest("known failer for now")
        df = DataFrame(np.random.randn(5, 3))
        self.assertRaises(TypeError, pd.eval,
                          '(df + 2)[df > 1] > 0 | (df > 0)',
                          local_dict={'df': df}, parser=self.parser,
                          engine=self.engine)


class TestOperationsPythonPandas(TestOperationsNumExprPandas):
    @classmethod
    def setUpClass(cls):
        cls.engine = 'python'
        cls.parser = 'pandas'


_var_s = randn(10)


class TestScope(object):
    def check_global_scope(self, e, engine, parser):
        skip_if_no_ne(engine)
        assert_array_equal(_var_s * 2, pd.eval(e, engine=engine,
                                               parser=parser))

    def test_global_scope(self):
        e = '_var_s * 2'
        for engine, parser in product(_engines, expr._parsers):
            yield self.check_global_scope, e, engine, parser

    def check_no_new_locals(self, engine, parser):
        skip_if_no_ne(engine)
        x = 1
        lcls = locals().copy()
        pd.eval('x + 1', local_dict=lcls, engine=engine, parser=parser)
        lcls2 = locals().copy()
        lcls2.pop('lcls')
        assert_equal(lcls, lcls2)

    def test_no_new_locals(self):
        for engine, parser in product(_engines, expr._parsers):
            yield self.check_no_new_locals, engine, parser

    def check_no_new_globals(self, engine, parser):
        skip_if_no_ne(engine)
        x = 1
        gbls = globals().copy()
        pd.eval('x + 1', engine=engine, parser=parser)
        gbls2 = globals().copy()
        assert_equal(gbls, gbls2)

    def test_no_new_globals(self):
        for engine, parser in product(_engines, expr._parsers):
            yield self.check_no_new_globals, engine, parser


def test_invalid_engine():
    skip_if_no_ne()
    assertRaisesRegexp(KeyError, 'Invalid engine \'asdf\' passed',
                       pd.eval, 'x + y', local_dict={'x': 1, 'y': 2},
                       engine='asdf')


def test_invalid_parser():
    skip_if_no_ne()
    assertRaisesRegexp(KeyError, 'Invalid parser \'asdf\' passed',
                       pd.eval, 'x + y', local_dict={'x': 1, 'y': 2},
                       parser='asdf')


def check_is_expr_syntax(engine):
    skip_if_no_ne(engine)
    s = 1
    valid1 = 's + 1'
    valid2 = '__y + _xx'
    assert_true(expr.isexpr(valid1, check_names=False))
    assert_true(expr.isexpr(valid2, check_names=False))


def check_is_expr_names(engine):
    skip_if_no_ne(engine)
    r, s = 1, 2
    valid = 's + r'
    invalid = '__y + __x'
    assert_true(expr.isexpr(valid, check_names=True))
    assert_false(expr.isexpr(invalid, check_names=True))


def test_is_expr_syntax():
    for engine in _engines:
        yield check_is_expr_syntax, engine


def test_is_expr_names():
    for engine in _engines:
        yield check_is_expr_names, engine


_parsers = {'python': PythonExprVisitor, 'pytables': pytables.ExprVisitor,
            'pandas': PandasExprVisitor}

def check_disallowed_nodes(engine, parser):
    skip_if_no_ne(engine)
    VisitorClass = _parsers[parser]
    uns_ops = VisitorClass.unsupported_nodes
    inst = VisitorClass('x + 1', engine, parser)

    for ops in uns_ops:
        assert_raises(NotImplementedError, getattr(inst, ops))


def test_disallowed_nodes():
    for engine, visitor in product(_parsers, repeat=2):
        yield check_disallowed_nodes, engine, visitor


def check_syntax_error_exprs(engine, parser):
    skip_if_no_ne(engine)
    e = 's +'
    assert_raises(SyntaxError, pd.eval, e, engine=engine, parser=parser)


def test_syntax_error_exprs():
    for engine, parser in ENGINES_PARSERS:
        yield check_syntax_error_exprs, engine, parser


def check_name_error_exprs(engine, parser):
    skip_if_no_ne(engine)
    e = 's + t'
    assert_raises(NameError, pd.eval, e, engine=engine, parser=parser)


def test_name_error_exprs():
    for engine, parser in ENGINES_PARSERS:
        yield check_name_error_exprs, engine, parser


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
