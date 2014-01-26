#!/usr/bin/env python

import functools
from itertools import product
from distutils.version import LooseVersion

import nose
from nose.tools import assert_raises, assert_true, assert_false, assert_equal

from numpy.random import randn, rand, randint
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
from numpy.testing.decorators import slow

import pandas as pd
from pandas.core import common as com
from pandas import DataFrame, Series, Panel, date_range
from pandas.util.testing import makeCustomDataframe as mkdf

from pandas.computation import pytables
from pandas.computation.engines import _engines
from pandas.computation.expr import PythonExprVisitor, PandasExprVisitor
from pandas.computation.ops import (_binary_ops_dict,
                                    _special_case_arith_ops_syms,
                                    _arith_ops_syms, _bool_ops_syms)
from pandas.computation.common import NameResolutionError

import pandas.computation.expr as expr
import pandas.util.testing as tm
from pandas.util.testing import (assert_frame_equal, randbool,
                                 assertRaisesRegexp,
                                 assert_produces_warning, assert_series_equal)
from pandas.compat import PY3, u

_series_frame_incompatible = _bool_ops_syms
_scalar_skip = 'in', 'not in'


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
                              engine=self.engine, parser=self.parser)
        elif (np.isscalar(lhs) and np.isscalar(rhs) and arith1 in
                _bool_ops_syms):
            with tm.assertRaises(NotImplementedError):
                pd.eval('lhs {0} rhs'.format(arith1), engine=self.engine,
                        parser=self.parser)
        else:
            f(self, lhs, arith1, rhs, *args, **kwargs)
    return wrapper


def _is_py3_complex_incompat(result, expected):
    return (PY3 and isinstance(expected, (complex, np.complexfloating)) and
            np.isnan(result))


_good_arith_ops = com.difference(_arith_ops_syms, _special_case_arith_ops_syms)


class TestEvalNumexprPandas(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestEvalNumexprPandas, cls).setUpClass()
        tm.skip_if_no_ne()
        import numexpr as ne
        cls.ne = ne
        cls.engine = 'numexpr'
        cls.parser = 'pandas'

    @classmethod
    def tearDownClass(cls):
        super(TestEvalNumexprPandas, cls).tearDownClass()
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
        self.scalar_lhses = randn(),
        self.scalar_rhses = randn(),

        self.lhses = self.pandas_lhses + self.scalar_lhses
        self.rhses = self.pandas_rhses + self.scalar_rhses

    def setup_ops(self):
        self.cmp_ops = expr._cmp_ops_syms
        self.cmp2_ops = self.cmp_ops[::-1]
        self.bin_ops = expr._bool_ops_syms
        self.special_case_ops = _special_case_arith_ops_syms
        self.arith_ops = _good_arith_ops
        self.unary_ops = '-', '~', 'not '

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

    def test_pow(self):
        import platform
        if platform.system() == 'Windows':
            raise nose.SkipTest('not testing pow on Windows')

        # odd failure on win32 platform, so skip
        for lhs, rhs in product(self.lhses, self.rhses):
            self.check_pow(lhs, '**', rhs)

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
        # tuple(set(self.cmp_ops) - set(['==', '!=', '<=', '>=']))
        cmp_ops = '<', '>'
        for lhs, cmp1, mid, cmp2, rhs in product(self.lhses, cmp_ops,
                                                 mids, cmp_ops, self.rhses):
            self.check_chained_cmp_op(lhs, cmp1, mid, cmp2, rhs)

    def check_complex_cmp_op(self, lhs, cmp1, rhs, binop, cmp2):
        skip_these = _scalar_skip
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
        elif _bool_and_frame(lhs, rhs):
            self.assertRaises(TypeError, _eval_single_bin, lhs, '&',
                              rhs, self.engine)
            self.assertRaises(TypeError, pd.eval, ex,
                              local_dict={'lhs': lhs, 'rhs': rhs},
                              engine=self.engine, parser=self.parser)
        elif (np.isscalar(lhs) and np.isnan(lhs) and
                not np.isscalar(rhs) and (cmp1 in skip_these or cmp2 in
                                          skip_these)):
            with tm.assertRaises(TypeError):
                _eval_single_bin(lhs, binop, rhs, self.engine)
        else:
            lhs_new = _eval_single_bin(lhs, cmp1, rhs, self.engine)
            rhs_new = _eval_single_bin(lhs, cmp2, rhs, self.engine)
            if (isinstance(lhs_new, Series) and isinstance(rhs_new, DataFrame)
                    and binop in _series_frame_incompatible):
                pass
                # TODO: the code below should be added back when left and right
                # hand side bool ops are fixed.

                # try:
                    # self.assertRaises(Exception, pd.eval, ex,
                                    #local_dict={'lhs': lhs, 'rhs': rhs},
                                    # engine=self.engine, parser=self.parser)
                # except AssertionError:
                    #import ipdb; ipdb.set_trace()
                    # raise
            elif (np.isscalar(lhs_new) and np.isnan(lhs_new) and
                    not np.isscalar(rhs_new) and binop in skip_these):
                with tm.assertRaises(TypeError):
                    _eval_single_bin(lhs_new, binop, rhs_new, self.engine)
            else:
                expected = _eval_single_bin(
                    lhs_new, binop, rhs_new, self.engine)
                result = pd.eval(ex, engine=self.engine, parser=self.parser)
                assert_array_equal(result, expected)

    @skip_incompatible_operand
    def check_chained_cmp_op(self, lhs, cmp1, mid, cmp2, rhs):
        skip_these = _scalar_skip

        def check_operands(left, right, cmp_op):
            if (np.isscalar(left) and np.isnan(left) and not np.isscalar(right)
                    and cmp_op in skip_these):
                ex = 'left {0} right'.format(cmp_op)
                with tm.assertRaises(ValueError):
                    pd.eval(ex, engine=self.engine, parser=self.parser)
                return
            if (np.isscalar(left) and np.isscalar(right) and
                    cmp_op in _bool_ops_syms):
                ex1 = 'lhs {0} mid {1} rhs'.format(cmp1, cmp2)
                ex2 = 'lhs {0} mid and mid {1} rhs'.format(cmp1, cmp2)
                ex3 = '(lhs {0} mid) & (mid {1} rhs)'.format(cmp1, cmp2)
                for ex in (ex1, ex2, ex3):
                    with tm.assertRaises(NotImplementedError):
                        pd.eval(ex, engine=self.engine, parser=self.parser)
                return
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
            ex1 = 'lhs {0} mid {1} rhs'.format(cmp1, cmp2)
            ex2 = 'lhs {0} mid and mid {1} rhs'.format(cmp1, cmp2)
            ex3 = '(lhs {0} mid) & (mid {1} rhs)'.format(cmp1, cmp2)
            expected = _eval_single_bin(lhs_new, '&', rhs_new, self.engine)

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
            assert_allclose(result, expected)

            ex = '(lhs {0} rhs) {0} rhs'.format(arith1)
            result = pd.eval(ex, engine=self.engine, parser=self.parser)
            expected = self.get_expected_pow_result(
                self.get_expected_pow_result(lhs, rhs), rhs)
            assert_allclose(result, expected)

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
                tm.skip_if_no_ne(engine)
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
        elif (np.isscalar(lhs) and np.isnan(lhs) and not np.isscalar(rhs)
                and cmp1 in skip_these):
            with tm.assertRaises(ValueError):
                pd.eval(ex, engine=self.engine, parser=self.parser)
        else:
            # compound
            if np.isscalar(lhs) and np.isscalar(rhs):
                lhs, rhs = map(lambda x: np.array([x]), (lhs, rhs))
            expected = _eval_single_bin(lhs, cmp1, rhs, self.engine)
            if np.isscalar(expected):
                expected = not expected
            else:
                expected = ~expected
            result = pd.eval(ex, engine=self.engine, parser=self.parser)
            assert_array_equal(expected, result)

            # make sure the other engines work the same as this one
            for engine in self.current_engines:
                tm.skip_if_no_ne(engine)
                ev = pd.eval(ex, engine=self.engine, parser=self.parser)
                assert_array_equal(ev, result)

    def ex(self, op, var_name='lhs'):
        return '{0}{1}'.format(op, var_name)

    def test_frame_invert(self):
        expr = self.ex('~')

        # ~ ##
        # frame
        # float always raises
        lhs = DataFrame(randn(5, 2))
        if self.engine == 'numexpr':
            with tm.assertRaises(NotImplementedError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)

        # int raises on numexpr
        lhs = DataFrame(randint(5, size=(5, 2)))
        if self.engine == 'numexpr':
            with tm.assertRaises(NotImplementedError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = ~lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_frame_equal(expect, result)

        # bool always works
        lhs = DataFrame(rand(5, 2) > 0.5)
        expect = ~lhs
        result = pd.eval(expr, engine=self.engine, parser=self.parser)
        assert_frame_equal(expect, result)

        # object raises
        lhs = DataFrame({'b': ['a', 1, 2.0], 'c': rand(3) > 0.5})
        if self.engine == 'numexpr':
            with tm.assertRaises(ValueError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)

    def test_series_invert(self):
        # ~ ####
        expr = self.ex('~')

        # series
        # float raises
        lhs = Series(randn(5))
        if self.engine == 'numexpr':
            with tm.assertRaises(NotImplementedError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)

        # int raises on numexpr
        lhs = Series(randint(5, size=5))
        if self.engine == 'numexpr':
            with tm.assertRaises(NotImplementedError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = ~lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_series_equal(expect, result)

        # bool
        lhs = Series(rand(5) > 0.5)
        expect = ~lhs
        result = pd.eval(expr, engine=self.engine, parser=self.parser)
        assert_series_equal(expect, result)

        # float
        # int
        # bool

        # object
        lhs = Series(['a', 1, 2.0])
        if self.engine == 'numexpr':
            with tm.assertRaises(ValueError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)

    def test_frame_negate(self):
        expr = self.ex('-')

        # float
        lhs = DataFrame(randn(5, 2))
        expect = -lhs
        result = pd.eval(expr, engine=self.engine, parser=self.parser)
        assert_frame_equal(expect, result)

        # int
        lhs = DataFrame(randint(5, size=(5, 2)))
        expect = -lhs
        result = pd.eval(expr, engine=self.engine, parser=self.parser)
        assert_frame_equal(expect, result)

        # bool doesn't work with numexpr but works elsewhere
        lhs = DataFrame(rand(5, 2) > 0.5)
        if self.engine == 'numexpr':
            with tm.assertRaises(NotImplementedError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = -lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_frame_equal(expect, result)

    def test_series_negate(self):
        expr = self.ex('-')

        # float
        lhs = Series(randn(5))
        expect = -lhs
        result = pd.eval(expr, engine=self.engine, parser=self.parser)
        assert_series_equal(expect, result)

        # int
        lhs = Series(randint(5, size=5))
        expect = -lhs
        result = pd.eval(expr, engine=self.engine, parser=self.parser)
        assert_series_equal(expect, result)

        # bool doesn't work with numexpr but works elsewhere
        lhs = Series(rand(5) > 0.5)
        if self.engine == 'numexpr':
            with tm.assertRaises(NotImplementedError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = -lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_series_equal(expect, result)

    def test_frame_pos(self):
        expr = self.ex('+')

        # float
        lhs = DataFrame(randn(5, 2))
        if self.engine == 'python':
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_frame_equal(expect, result)

        # int
        lhs = DataFrame(randint(5, size=(5, 2)))
        if self.engine == 'python':
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_frame_equal(expect, result)

        # bool doesn't work with numexpr but works elsewhere
        lhs = DataFrame(rand(5, 2) > 0.5)
        if self.engine == 'python':
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_frame_equal(expect, result)

    def test_series_pos(self):
        expr = self.ex('+')

        # float
        lhs = Series(randn(5))
        if self.engine == 'python':
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_series_equal(expect, result)

        # int
        lhs = Series(randint(5, size=5))
        if self.engine == 'python':
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_series_equal(expect, result)

        # bool doesn't work with numexpr but works elsewhere
        lhs = Series(rand(5) > 0.5)
        if self.engine == 'python':
            with tm.assertRaises(TypeError):
                result = pd.eval(expr, engine=self.engine, parser=self.parser)
        else:
            expect = lhs
            result = pd.eval(expr, engine=self.engine, parser=self.parser)
            assert_series_equal(expect, result)

    def test_scalar_unary(self):
        with tm.assertRaises(TypeError):
            pd.eval('~1.0', engine=self.engine, parser=self.parser)

        self.assertEqual(
            pd.eval('-1.0', parser=self.parser, engine=self.engine), -1.0)
        self.assertEqual(
            pd.eval('+1.0', parser=self.parser, engine=self.engine), +1.0)

        self.assertEqual(
            pd.eval('~1', parser=self.parser, engine=self.engine), ~1)
        self.assertEqual(
            pd.eval('-1', parser=self.parser, engine=self.engine), -1)
        self.assertEqual(
            pd.eval('+1', parser=self.parser, engine=self.engine), +1)

        self.assertEqual(
            pd.eval('~True', parser=self.parser, engine=self.engine), ~True)
        self.assertEqual(
            pd.eval('~False', parser=self.parser, engine=self.engine), ~False)
        self.assertEqual(
            pd.eval('-True', parser=self.parser, engine=self.engine), -True)
        self.assertEqual(
            pd.eval('-False', parser=self.parser, engine=self.engine), -False)
        self.assertEqual(
            pd.eval('+True', parser=self.parser, engine=self.engine), +True)
        self.assertEqual(
            pd.eval('+False', parser=self.parser, engine=self.engine), +False)

    def test_disallow_scalar_bool_ops(self):
        exprs = '1 or 2', '1 and 2'
        exprs += 'a and b', 'a or b'
        exprs += '1 or 2 and (3 + 2) > 3',
        exprs += '2 * x > 2 or 1 and 2',
        exprs += '2 * df > 3 and 1 or a',

        x, a, b, df = np.random.randn(3), 1, 2, DataFrame(randn(3, 2))
        for ex in exprs:
            with tm.assertRaises(NotImplementedError):
                pd.eval(ex, engine=self.engine, parser=self.parser)


class TestEvalNumexprPython(TestEvalNumexprPandas):

    @classmethod
    def setUpClass(cls):
        super(TestEvalNumexprPython, cls).setUpClass()
        tm.skip_if_no_ne()
        import numexpr as ne
        cls.ne = ne
        cls.engine = 'numexpr'
        cls.parser = 'python'

    def setup_ops(self):
        self.cmp_ops = list(filter(lambda x: x not in ('in', 'not in'),
                                   expr._cmp_ops_syms))
        self.cmp2_ops = self.cmp_ops[::-1]
        self.bin_ops = [s for s in expr._bool_ops_syms
                        if s not in ('and', 'or')]
        self.special_case_ops = _special_case_arith_ops_syms
        self.arith_ops = _good_arith_ops
        self.unary_ops = '+', '-', '~'

    def check_chained_cmp_op(self, lhs, cmp1, mid, cmp2, rhs):
        ex1 = 'lhs {0} mid {1} rhs'.format(cmp1, cmp2)
        with tm.assertRaises(NotImplementedError):
            pd.eval(ex1, engine=self.engine, parser=self.parser)


class TestEvalPythonPython(TestEvalNumexprPython):

    @classmethod
    def setUpClass(cls):
        super(TestEvalPythonPython, cls).setUpClass()
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
        super(TestEvalPythonPandas, cls).setUpClass()
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

    index_types = 'i', 'u', 'dt'
    lhs_index_types = index_types + ('s',)  # 'p'

    def check_align_nested_unary_op(self, engine, parser):
        tm.skip_if_no_ne(engine)
        s = 'df * ~2'
        df = mkdf(5, 3, data_gen_f=f)
        res = pd.eval(s, engine=engine, parser=parser)
        assert_frame_equal(res, df * ~2)

    def test_align_nested_unary_op(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_align_nested_unary_op, engine, parser

    def check_basic_frame_alignment(self, engine, parser):
        tm.skip_if_no_ne(engine)
        args = product(self.lhs_index_types, self.index_types,
                       self.index_types)
        for lr_idx_type, rr_idx_type, c_idx_type in args:
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=lr_idx_type,
                      c_idx_type=c_idx_type)
            df2 = mkdf(20, 10, data_gen_f=f, r_idx_type=rr_idx_type,
                       c_idx_type=c_idx_type)
            res = pd.eval('df + df2', engine=engine, parser=parser)
            assert_frame_equal(res, df + df2)

    def test_basic_frame_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_basic_frame_alignment, engine, parser

    def check_frame_comparison(self, engine, parser):
        tm.skip_if_no_ne(engine)
        args = product(self.lhs_index_types, repeat=2)
        for r_idx_type, c_idx_type in args:
            df = mkdf(10, 10, data_gen_f=f, r_idx_type=r_idx_type,
                      c_idx_type=c_idx_type)
            res = pd.eval('df < 2', engine=engine, parser=parser)
            assert_frame_equal(res, df < 2)

            df3 = DataFrame(randn(*df.shape), index=df.index,
                            columns=df.columns)
            res = pd.eval('df < df3', engine=engine, parser=parser)
            assert_frame_equal(res, df < df3)

    def test_frame_comparison(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_frame_comparison, engine, parser

    def check_medium_complex_frame_alignment(self, engine, parser):
        tm.skip_if_no_ne(engine)
        args = product(self.lhs_index_types, self.index_types,
                       self.index_types, self.index_types)

        for r1, c1, r2, c2 in args:
            df = mkdf(3, 2, data_gen_f=f, r_idx_type=r1, c_idx_type=c1)
            df2 = mkdf(4, 2, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
            df3 = mkdf(5, 2, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
            res = pd.eval('df + df2 + df3', engine=engine, parser=parser)
            assert_frame_equal(res, df + df2 + df3)

    @slow
    def test_medium_complex_frame_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_medium_complex_frame_alignment, engine, parser

    def check_basic_frame_series_alignment(self, engine, parser):
        tm.skip_if_no_ne(engine)

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

        args = product(self.lhs_index_types, self.index_types,
                       ('index', 'columns'))
        for r_idx_type, c_idx_type, index_name in args:
            testit(r_idx_type, c_idx_type, index_name)

    def test_basic_frame_series_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_basic_frame_series_alignment, engine, parser

    def check_basic_series_frame_alignment(self, engine, parser):
        tm.skip_if_no_ne(engine)

        def testit(r_idx_type, c_idx_type, index_name):
            df = mkdf(10, 7, data_gen_f=f, r_idx_type=r_idx_type,
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

        # only test dt with dt, otherwise weird joins result
        args = product(['i', 'u', 's'], ['i', 'u', 's'], ('index', 'columns'))
        for r_idx_type, c_idx_type, index_name in args:
            testit(r_idx_type, c_idx_type, index_name)

        # dt with dt
        args = product(['dt'], ['dt'], ('index', 'columns'))
        for r_idx_type, c_idx_type, index_name in args:
            testit(r_idx_type, c_idx_type, index_name)

    def test_basic_series_frame_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_basic_series_frame_alignment, engine, parser

    def check_series_frame_commutativity(self, engine, parser):
        tm.skip_if_no_ne(engine)
        args = product(self.lhs_index_types, self.index_types, ('+', '*'),
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

    def test_series_frame_commutativity(self):
        for engine, parser in ENGINES_PARSERS:
            yield self.check_series_frame_commutativity, engine, parser

    def check_complex_series_frame_alignment(self, engine, parser):
        tm.skip_if_no_ne(engine)

        import random
        args = product(self.lhs_index_types, self.index_types,
                       self.index_types, self.index_types)
        n = 3
        m1 = 5
        m2 = 2 * m1

        for r1, r2, c1, c2 in args:
            index_name = random.choice(['index', 'columns'])
            obj_name = random.choice(['df', 'df2'])

            df = mkdf(m1, n, data_gen_f=f, r_idx_type=r1, c_idx_type=c1)
            df2 = mkdf(m2, n, data_gen_f=f, r_idx_type=r2, c_idx_type=c2)
            index = getattr(locals().get(obj_name), index_name)
            s = Series(np.random.randn(n), index[:n])

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

    def check_performance_warning_for_poor_alignment(self, engine, parser):
        tm.skip_if_no_ne(engine)
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

        df = DataFrame(randn(10, 10))
        s = Series(randn(10000))

        is_python_engine = engine == 'python'

        if not is_python_engine:
            wrn = pd.io.common.PerformanceWarning
        else:
            wrn = False

        with assert_produces_warning(wrn) as w:
            pd.eval('df + s', engine=engine, parser=parser)

            if not is_python_engine:
                assert_equal(len(w), 1)
                msg = str(w[0].message)
                expected = ("Alignment difference on axis {0} is larger"
                            " than an order of magnitude on term {1!r}, "
                            "by more than {2:.4g}; performance may suffer"
                            "".format(1, 'df', np.log10(s.size - df.shape[1])))
                assert_equal(msg, expected)

    def test_performance_warning_for_poor_alignment(self):
        for engine, parser in ENGINES_PARSERS:
            yield (self.check_performance_warning_for_poor_alignment, engine,
                   parser)


#------------------------------------
# slightly more complex ops

class TestOperationsNumExprPandas(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestOperationsNumExprPandas, cls).setUpClass()
        tm.skip_if_no_ne()
        cls.engine = 'numexpr'
        cls.parser = 'pandas'
        cls.arith_ops = expr._arith_ops_syms + expr._cmp_ops_syms

    @classmethod
    def tearDownClass(cls):
        super(TestOperationsNumExprPandas, cls).tearDownClass()
        del cls.engine, cls.parser

    def eval(self, *args, **kwargs):
        kwargs['engine'] = self.engine
        kwargs['parser'] = self.parser
        return pd.eval(*args, **kwargs)

    def test_simple_arith_ops(self):
        ops = self.arith_ops

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
            ex = '{0} {1} {2}'.format(lhs, op, rhs)
            res = self.eval(ex)
            exp = eval(ex)
            self.assertEqual(res, exp)

    def test_bool_ops_with_constants(self):
        for op, lhs, rhs in product(expr._bool_ops_syms, ('True', 'False'),
                                    ('True', 'False')):
            ex = '{0} {1} {2}'.format(lhs, op, rhs)
            res = self.eval(ex)
            exp = eval(ex)
            self.assertEqual(res, exp)

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
        self.assertRaises(ValueError, self.eval, expr1,
                          local_dict={'df': df, 'df2': df2})

    def test_assignment_column(self):
        tm.skip_if_no_ne('numexpr')
        df = DataFrame(np.random.randn(5, 2), columns=list('ab'))
        orig_df = df.copy()

        # multiple assignees
        self.assertRaises(SyntaxError, df.eval, 'd c = a + b')

        # invalid assignees
        self.assertRaises(SyntaxError, df.eval, 'd,c = a + b')
        self.assertRaises(
            SyntaxError, df.eval, 'Timestamp("20131001") = a + b')

        # single assignment - existing variable
        expected = orig_df.copy()
        expected['a'] = expected['a'] + expected['b']
        df = orig_df.copy()
        df.eval('a = a + b')
        assert_frame_equal(df, expected)

        # single assignment - new variable
        expected = orig_df.copy()
        expected['c'] = expected['a'] + expected['b']
        df = orig_df.copy()
        df.eval('c = a + b')
        assert_frame_equal(df, expected)

        # with a local name overlap
        def f():
            df = orig_df.copy()
            a = 1
            df.eval('a = 1 + b')
            return df

        df = f()
        expected = orig_df.copy()
        expected['a'] = 1 + expected['b']
        assert_frame_equal(df, expected)

        df = orig_df.copy()

        def f():
            a = 1
            df.eval('a=a+b')
        self.assertRaises(NameResolutionError, f)

        # multiple assignment
        df = orig_df.copy()
        df.eval('c = a + b')
        self.assertRaises(SyntaxError, df.eval, 'c = a = b')

        # explicit targets
        df = orig_df.copy()
        self.eval('c = df.a + df.b', local_dict={'df': df}, target=df)
        expected = orig_df.copy()
        expected['c'] = expected['a'] + expected['b']
        assert_frame_equal(df, expected)

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
        res = self.eval('df.dates1 < 20130101', local_dict={'df': df},
                        engine=self.engine, parser=self.parser)
        expec = df.dates1 < '20130101'
        assert_series_equal(res, expec)

    def test_simple_in_ops(self):
        if self.parser != 'python':
            res = pd.eval('1 in [1, 2]', engine=self.engine,
                          parser=self.parser)
            self.assertTrue(res)

            res = pd.eval('2 in (1, 2)', engine=self.engine,
                          parser=self.parser)
            self.assertTrue(res)

            res = pd.eval('3 in (1, 2)', engine=self.engine,
                          parser=self.parser)
            self.assertFalse(res)

            res = pd.eval('3 not in (1, 2)', engine=self.engine,
                          parser=self.parser)
            self.assertTrue(res)

            res = pd.eval('[3] not in (1, 2)', engine=self.engine,
                          parser=self.parser)
            self.assertTrue(res)

            res = pd.eval('[3] in ([3], 2)', engine=self.engine,
                          parser=self.parser)
            self.assertTrue(res)

            res = pd.eval('[[3]] in [[[3]], 2]', engine=self.engine,
                          parser=self.parser)
            self.assertTrue(res)

            res = pd.eval('(3,) in [(3,), 2]', engine=self.engine,
                          parser=self.parser)
            self.assertTrue(res)

            res = pd.eval('(3,) not in [(3,), 2]', engine=self.engine,
                          parser=self.parser)
            self.assertFalse(res)

            res = pd.eval('[(3,)] in [[(3,)], 2]', engine=self.engine,
                          parser=self.parser)
            self.assertTrue(res)
        else:
            with tm.assertRaises(NotImplementedError):
                pd.eval('1 in [1, 2]', engine=self.engine, parser=self.parser)
            with tm.assertRaises(NotImplementedError):
                pd.eval('2 in (1, 2)', engine=self.engine, parser=self.parser)
            with tm.assertRaises(NotImplementedError):
                pd.eval('3 in (1, 2)', engine=self.engine, parser=self.parser)
            with tm.assertRaises(NotImplementedError):
                pd.eval('3 not in (1, 2)', engine=self.engine,
                        parser=self.parser)
            with tm.assertRaises(NotImplementedError):
                pd.eval('[(3,)] in (1, 2, [(3,)])', engine=self.engine,
                        parser=self.parser)
            with tm.assertRaises(NotImplementedError):
                pd.eval('[3] not in (1, 2, [[3]])', engine=self.engine,
                        parser=self.parser)


class TestOperationsNumExprPython(TestOperationsNumExprPandas):

    @classmethod
    def setUpClass(cls):
        super(TestOperationsNumExprPython, cls).setUpClass()
        cls.engine = 'numexpr'
        cls.parser = 'python'
        tm.skip_if_no_ne(cls.engine)
        cls.arith_ops = expr._arith_ops_syms + expr._cmp_ops_syms
        cls.arith_ops = filter(lambda x: x not in ('in', 'not in'),
                               cls.arith_ops)

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
        ex = '(df + 2)[df > 1] > 0 & (df > 0)'
        with tm.assertRaises(NotImplementedError):
            pd.eval(ex, parser=self.parser, engine=self.engine)

    def test_fails_pipe(self):
        df = DataFrame(np.random.randn(5, 3))
        ex = '(df + 2)[df > 1] > 0 | (df > 0)'
        with tm.assertRaises(NotImplementedError):
            pd.eval(ex, parser=self.parser, engine=self.engine)

    def test_bool_ops_with_constants(self):
        for op, lhs, rhs in product(expr._bool_ops_syms, ('True', 'False'),
                                    ('True', 'False')):
            ex = '{0} {1} {2}'.format(lhs, op, rhs)
            if op in ('and', 'or'):
                with tm.assertRaises(NotImplementedError):
                    self.eval(ex)
            else:
                res = self.eval(ex)
                exp = eval(ex)
                self.assertEqual(res, exp)

    def test_simple_bool_ops(self):
        for op, lhs, rhs in product(expr._bool_ops_syms, (True, False),
                                    (True, False)):
            ex = 'lhs {0} rhs'.format(op)
            if op in ('and', 'or'):
                with tm.assertRaises(NotImplementedError):
                    pd.eval(ex, engine=self.engine, parser=self.parser)
            else:
                res = pd.eval(ex, engine=self.engine, parser=self.parser)
                exp = eval(ex)
                self.assertEqual(res, exp)


class TestOperationsPythonPython(TestOperationsNumExprPython):

    @classmethod
    def setUpClass(cls):
        super(TestOperationsPythonPython, cls).setUpClass()
        cls.engine = cls.parser = 'python'
        cls.arith_ops = expr._arith_ops_syms + expr._cmp_ops_syms
        cls.arith_ops = filter(lambda x: x not in ('in', 'not in'),
                               cls.arith_ops)


class TestOperationsPythonPandas(TestOperationsNumExprPandas):

    @classmethod
    def setUpClass(cls):
        super(TestOperationsPythonPandas, cls).setUpClass()
        cls.engine = 'python'
        cls.parser = 'pandas'
        cls.arith_ops = expr._arith_ops_syms + expr._cmp_ops_syms


_var_s = randn(10)


class TestScope(object):

    def check_global_scope(self, e, engine, parser):
        tm.skip_if_no_ne(engine)
        assert_array_equal(_var_s * 2, pd.eval(e, engine=engine,
                                               parser=parser))

    def test_global_scope(self):
        e = '_var_s * 2'
        for engine, parser in product(_engines, expr._parsers):
            yield self.check_global_scope, e, engine, parser

    def check_no_new_locals(self, engine, parser):
        tm.skip_if_no_ne(engine)
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
        tm.skip_if_no_ne(engine)
        x = 1
        gbls = globals().copy()
        pd.eval('x + 1', engine=engine, parser=parser)
        gbls2 = globals().copy()
        assert_equal(gbls, gbls2)

    def test_no_new_globals(self):
        for engine, parser in product(_engines, expr._parsers):
            yield self.check_no_new_globals, engine, parser


def test_invalid_engine():
    tm.skip_if_no_ne()
    assertRaisesRegexp(KeyError, 'Invalid engine \'asdf\' passed',
                       pd.eval, 'x + y', local_dict={'x': 1, 'y': 2},
                       engine='asdf')


def test_invalid_parser():
    tm.skip_if_no_ne()
    assertRaisesRegexp(KeyError, 'Invalid parser \'asdf\' passed',
                       pd.eval, 'x + y', local_dict={'x': 1, 'y': 2},
                       parser='asdf')


def check_is_expr_syntax(engine):
    tm.skip_if_no_ne(engine)
    s = 1
    valid1 = 's + 1'
    valid2 = '__y + _xx'
    assert_true(expr.isexpr(valid1, check_names=False))
    assert_true(expr.isexpr(valid2, check_names=False))


def check_is_expr_names(engine):
    tm.skip_if_no_ne(engine)
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
    tm.skip_if_no_ne(engine)
    VisitorClass = _parsers[parser]
    uns_ops = VisitorClass.unsupported_nodes
    inst = VisitorClass('x + 1', engine, parser)

    for ops in uns_ops:
        assert_raises(NotImplementedError, getattr(inst, ops))


def test_disallowed_nodes():
    for engine, visitor in product(_parsers, repeat=2):
        yield check_disallowed_nodes, engine, visitor


def check_syntax_error_exprs(engine, parser):
    tm.skip_if_no_ne(engine)
    e = 's +'
    assert_raises(SyntaxError, pd.eval, e, engine=engine, parser=parser)


def test_syntax_error_exprs():
    for engine, parser in ENGINES_PARSERS:
        yield check_syntax_error_exprs, engine, parser


def check_name_error_exprs(engine, parser):
    tm.skip_if_no_ne(engine)
    e = 's + t'
    assert_raises(NameError, pd.eval, e, engine=engine, parser=parser)


def test_name_error_exprs():
    for engine, parser in ENGINES_PARSERS:
        yield check_name_error_exprs, engine, parser


def check_invalid_numexpr_version(engine, parser):
    def testit():
        a, b = 1, 2
        res = pd.eval('a + b', engine=engine, parser=parser)
        tm.assert_equal(res, 3)

    if engine == 'numexpr':
        try:
            import numexpr as ne
        except ImportError:
            raise nose.SkipTest("no numexpr")
        else:
            if ne.__version__ < LooseVersion('2.0'):
                with tm.assertRaisesRegexp(ImportError, "'numexpr' version is "
                                           ".+, must be >= 2.0"):
                    testit()
            else:
                testit()
    else:
        testit()


def test_invalid_numexpr_version():
    for engine, parser in ENGINES_PARSERS:
        yield check_invalid_numexpr_version, engine, parser


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
