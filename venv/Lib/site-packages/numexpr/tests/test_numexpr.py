
###################################################################
#  Numexpr - Fast numerical array expression evaluator for NumPy.
#
#      License: MIT
#      Author:  See AUTHORS.txt
#
#  See LICENSE.txt and LICENSES/*.txt for details about copyright and
#  rights to use.
####################################################################


import os
import sys
import platform
import warnings
from contextlib import contextmanager
import subprocess

import numpy as np
from numpy import (
    array, arange, empty, zeros, int32, int64, uint16, cdouble, float64, rec,
    copy, ones_like, where, all as alltrue, linspace,
    sum, prod, sqrt, fmod, floor, ceil,
    sin, cos, tan, arcsin, arccos, arctan, arctan2,
    sinh, cosh, tanh, arcsinh, arccosh, arctanh,
    log, log1p, log10, exp, expm1, conj)
import numpy
from numpy.testing import (assert_equal, assert_array_equal,
                           assert_array_almost_equal, assert_allclose)
from numpy import shape, allclose, array_equal, ravel, isnan, isinf

import numexpr
from numexpr import E, NumExpr, evaluate, re_evaluate, validate, disassemble, use_vml
from numexpr.expressions import ConstantNode
from numexpr.utils import detect_number_of_cores

import unittest

TestCase = unittest.TestCase

double = np.double
long = int
MAX_THREADS = 16


class test_numexpr(TestCase):
    """Testing with 1 thread"""
    nthreads = 1

    def setUp(self):
        numexpr.set_num_threads(self.nthreads)

    def test_simple(self):
        ex = 2.0 * E.a + 3.0 * E.b * E.c
        sig = [('a', double), ('b', double), ('c', double)]
        func = NumExpr(ex, signature=sig)
        x = func(array([1., 2, 3]), array([4., 5, 6]), array([7., 8, 9]))
        assert_array_equal(x, array([86., 124., 168.]))

    def test_simple_expr_small_array(self):
        func = NumExpr(E.a)
        x = arange(100.0)
        y = func(x)
        assert_array_equal(x, y)

    def test_simple_expr(self):
        func = NumExpr(E.a)
        x = arange(1e6)
        y = func(x)
        assert_array_equal(x, y)

    def test_rational_expr(self):
        func = NumExpr((E.a + 2.0 * E.b) / (1 + E.a + 4 * E.b * E.b))
        a = arange(1e6)
        b = arange(1e6) * 0.1
        x = (a + 2 * b) / (1 + a + 4 * b * b)
        y = func(a, b)
        assert_array_almost_equal(x, y)

    def test_reductions(self):
        # Check that they compile OK.
        assert_equal(disassemble(
            NumExpr("sum(x**2+2, axis=None)", [('x', double)])),
                     [(b'mul_ddd', b't3', b'r1[x]', b'r1[x]'),
                      (b'add_ddd', b't3', b't3', b'c2[2.0]'),
                      (b'sum_ddn', b'r0', b't3', None)])
        assert_equal(disassemble(
            NumExpr("sum(x**2+2, axis=1)", [('x', double)])),
                     [(b'mul_ddd', b't3', b'r1[x]', b'r1[x]'),
                      (b'add_ddd', b't3', b't3', b'c2[2.0]'),
                      (b'sum_ddn', b'r0', b't3', 1)])
        assert_equal(disassemble(
            NumExpr("prod(x**2+2, axis=2)", [('x', double)])),
                     [(b'mul_ddd', b't3', b'r1[x]', b'r1[x]'),
                      (b'add_ddd', b't3', b't3', b'c2[2.0]'),
                      (b'prod_ddn', b'r0', b't3', 2)])
        # Check that full reductions work.
        x = zeros(100000) + .01  # checks issue #41
        assert_allclose(evaluate("sum(x+2,axis=None)"), sum(x + 2, axis=None))
        assert_allclose(evaluate("sum(x+2,axis=0)"), sum(x + 2, axis=0))
        assert_allclose(evaluate("prod(x,axis=0)"), prod(x, axis=0))
        assert_allclose(evaluate("min(x)"), np.min(x))
        assert_allclose(evaluate("max(x,axis=0)"), np.max(x, axis=0))

        # Fix for #277, array with leading singleton dimension
        x = np.arange(10).reshape(1,10)
        assert_allclose(evaluate("sum(x,axis=None)"), sum(x, axis=None) )
        assert_allclose(evaluate("sum(x,axis=0)"), sum(x, axis=0) )
        assert_allclose(evaluate("sum(x,axis=1)"), sum(x, axis=1) )

        x = arange(10.0)
        assert_allclose(evaluate("sum(x**2+2,axis=0)"), sum(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("prod(x**2+2,axis=0)"), prod(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("min(x**2+2,axis=0)"), np.min(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("max(x**2+2,axis=0)"), np.max(x ** 2 + 2, axis=0))

        x = arange(100.0)
        assert_allclose(evaluate("sum(x**2+2,axis=0)"), sum(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("prod(x-1,axis=0)"), prod(x - 1, axis=0))
        assert_allclose(evaluate("min(x-1,axis=0)"), np.min(x - 1, axis=0))
        assert_allclose(evaluate("max(x-1,axis=0)"), np.max(x - 1, axis=0))
        x = linspace(0.1, 1.0, 2000)
        assert_allclose(evaluate("sum(x**2+2,axis=0)"), sum(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("prod(x-1,axis=0)"), prod(x - 1, axis=0))
        assert_allclose(evaluate("min(x-1,axis=0)"), np.min(x - 1, axis=0))
        assert_allclose(evaluate("max(x-1,axis=0)"), np.max(x - 1, axis=0))

        # Check that reductions along an axis work
        y = arange(9.0).reshape(3, 3)
        assert_allclose(evaluate("sum(y**2, axis=1)"), sum(y ** 2, axis=1))
        assert_allclose(evaluate("sum(y**2, axis=0)"), sum(y ** 2, axis=0))
        assert_allclose(evaluate("sum(y**2, axis=None)"), sum(y ** 2, axis=None))
        assert_allclose(evaluate("prod(y**2, axis=1)"), prod(y ** 2, axis=1))
        assert_allclose(evaluate("prod(y**2, axis=0)"), prod(y ** 2, axis=0))
        assert_allclose(evaluate("prod(y**2, axis=None)"), prod(y ** 2, axis=None))
        assert_allclose(evaluate("min(y**2, axis=1)"), np.min(y ** 2, axis=1))
        assert_allclose(evaluate("min(y**2, axis=0)"), np.min(y ** 2, axis=0))
        assert_allclose(evaluate("min(y**2, axis=None)"), np.min(y ** 2, axis=None))
        assert_allclose(evaluate("max(y**2, axis=1)"), np.max(y ** 2, axis=1))
        assert_allclose(evaluate("max(y**2, axis=0)"), np.max(y ** 2, axis=0))
        assert_allclose(evaluate("max(y**2, axis=None)"), np.max(y ** 2, axis=None))
        # Check integers
        x = arange(10.)
        x = x.astype(int)
        assert_allclose(evaluate("sum(x**2+2,axis=0)"), sum(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("prod(x**2+2,axis=0)"), prod(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("min(x**2+2,axis=0)"), np.min(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("max(x**2+2,axis=0)"), np.max(x ** 2 + 2, axis=0))
        # Check longs
        x = x.astype(int)
        assert_allclose(evaluate("sum(x**2+2,axis=0)"), sum(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("prod(x**2+2,axis=0)"), prod(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("min(x**2+2,axis=0)"), np.min(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("max(x**2+2,axis=0)"), np.max(x ** 2 + 2, axis=0))
        # Check complex
        x = x + .1j
        assert_allclose(evaluate("sum(x**2+2,axis=0)"), sum(x ** 2 + 2, axis=0))
        assert_allclose(evaluate("prod(x-1,axis=0)"), prod(x - 1, axis=0))

    def test_in_place(self):
        x = arange(10000.).reshape(1000, 10)
        evaluate("x + 3", out=x)
        assert_equal(x, arange(10000.).reshape(1000, 10) + 3)
        y = arange(10)
        evaluate("(x - 3) * y + (x - 3)", out=x)
        assert_equal(x, arange(10000.).reshape(1000, 10) * (arange(10) + 1))

    def test_axis(self):
        y = arange(9.0).reshape(3, 3)
        try:
            evaluate("sum(y, axis=2)")
        except ValueError:
            pass
        else:
            raise ValueError("should raise exception!")
        try:
            evaluate("sum(y, axis=-3)")
        except ValueError:
            pass
        else:
            raise ValueError("should raise exception!")
        try:
            # Negative axis are not supported
            evaluate("sum(y, axis=-1)")
        except ValueError:
            pass
        else:
            raise ValueError("should raise exception!")

    def test_r0_reuse(self):
        assert_equal(disassemble(NumExpr("x * x + 2", [('x', double)])),
                     [(b'mul_ddd', b'r0', b'r1[x]', b'r1[x]'),
                      (b'add_ddd', b'r0', b'r0', b'c2[2.0]')])

    def test_str_contains_basic0(self):
        res = evaluate('contains(b"abc", b"ab")')
        assert_equal(res, True)

    def test_str_contains_basic1(self):
        haystack = array([b'abc', b'def', b'xyz', b'x11', b'za'])
        res = evaluate('contains(haystack, b"ab")')
        assert_equal(res, [True, False, False, False, False])

    def test_str_contains_basic2(self):
        haystack = array([b'abc', b'def', b'xyz', b'x11', b'za'])
        res = evaluate('contains(b"abcd", haystack)')
        assert_equal(res, [True, False, False, False, False])

    def test_str_contains_basic3(self):
        haystacks = array(
            [b'abckkk', b'adef', b'xyz', b'x11abcp', b'za', b'abc'])
        needles = array(
            [b'abc', b'def', b'aterr', b'oot', b'zu', b'ab'])
        res = evaluate('contains(haystacks, needles)')
        assert_equal(res, [True, True, False, False, False, True])

    def test_str_contains_basic4(self):
        needles = array(
            [b'abc', b'def', b'aterr', b'oot', b'zu', b'ab c', b' abc',
             b'abc '])
        res = evaluate('contains(b"test abc here", needles)')
        assert_equal(res, [True, False, False, False, False, False, True, True])

    def test_str_contains_basic5(self):
        needles = array(
            [b'abc', b'ab c', b' abc', b' abc ', b'\tabc', b'c h'])
        res = evaluate('contains(b"test abc here", needles)')
        assert_equal(res, [True, False, True, True, False, True])

        # Compare operation of Python 'in' operator with 'contains' using a
        # product of two lists of strings.

    def test_str_contains_listproduct(self):
        from itertools import product

        small = [
            'It w', 'as th', 'e Whit', 'e Rab', 'bit,', ' tro', 'tting',
            ' sl', 'owly', ' back ', 'again,', ' and', ' lo', 'okin', 'g a',
            'nxious', 'ly a', 'bou', 't a', 's it w', 'ent,', ' as i', 'f it',
            ' had l', 'ost', ' some', 'thi', 'ng; a', 'nd ', 'she ', 'heard ',
            'it mut', 'terin', 'g to ', 'its', 'elf ', "'The",
            ' Duch', 'ess! T', 'he ', 'Duches', 's! Oh ', 'my dea', 'r paws',
            '! Oh ', 'my f', 'ur ', 'and ', 'whiske', 'rs! ', 'She', "'ll g",
            'et me', ' ex', 'ecu', 'ted, ', 'as su', 're a', 's f', 'errets',
            ' are f', 'errets', '! Wh', 'ere ', 'CAN', ' I hav', 'e d',
            'roppe', 'd t', 'hem,', ' I wo', 'nder?', "' A", 'lice',
            ' gu', 'essed', ' in a', ' mom', 'ent ', 'tha', 't it w', 'as ',
            'looki', 'ng f', 'or ', 'the fa', 'n and ', 'the', ' pai',
            'r of w', 'hit', 'e kid', ' glo', 'ves', ', and ', 'she ',
            'very g', 'ood', '-na', 'turedl', 'y be', 'gan h', 'unt', 'ing',
            ' about', ' for t', 'hem', ', but', ' they ', 'wer', 'e nowh',
            'ere to', ' be', ' se', 'en--', 'ever', 'ythin', 'g seem', 'ed ',
            'to ', 'have c', 'hang', 'ed ', 'since', ' he', 'r swim', ' in',
            ' the', ' pool,', ' and', ' the g', 'reat ', 'hal', 'l, w', 'ith',
            ' th', 'e gl', 'ass t', 'abl', 'e and ', 'the', ' li', 'ttle',
            ' doo', 'r, ha', 'd v', 'ani', 'shed c', 'omp', 'lete', 'ly.']
        big = [
            'It wa', 's the', ' W', 'hit', 'e ', 'Ra', 'bb', 'it, t', 'ro',
            'tting s', 'lowly', ' back ', 'agai', 'n, and', ' l', 'ookin',
            'g ', 'an', 'xiously', ' about ', 'as it w', 'ent, as', ' if ',
            'it had', ' los', 't ', 'so', 'mething', '; and', ' she h',
            'eard ', 'it ', 'mutteri', 'ng to', ' itself', " 'The ",
            'Duchess', '! ', 'Th', 'e ', 'Duchess', '! Oh m', 'y de',
            'ar paws', '! ', 'Oh my ', 'fu', 'r and w', 'hiskers', "! She'",
            'll ', 'get', ' me ', 'execute', 'd,', ' a', 's ', 'su', 're as ',
            'fe', 'rrets', ' are f', 'errets!', ' Wher', 'e CAN', ' I ha',
            've dro', 'pped t', 'hem', ', I ', 'won', "der?' A",
            'lice g', 'uess', 'ed ', 'in a m', 'omen', 't that', ' i',
            't was l', 'ook', 'ing f', 'or th', 'e ', 'fan and', ' th', 'e p',
            'air o', 'f whit', 'e ki', 'd glove', 's, and ', 'she v', 'ery ',
            'good-na', 'tu', 'redl', 'y be', 'gan hun', 'ti', 'ng abou',
            't for t', 'he', 'm, bu', 't t', 'hey ', 'were n', 'owhere',
            ' to b', 'e s', 'een-', '-eve', 'rythi', 'ng see', 'me', 'd ',
            'to ha', 've', ' c', 'hanged', ' sinc', 'e her s', 'wim ',
            'in the ', 'pool,', ' an', 'd the g', 'rea', 't h', 'all, wi',
            'th the ', 'glas', 's t', 'able an', 'd th', 'e littl', 'e door,',
            ' had va', 'ni', 'shed co', 'mpletel', 'y.']
        p = list(product(small, big))
        python_in = [x[0] in x[1] for x in p]
        a = [x[0].encode() for x in p]
        b = [x[1].encode() for x in p]
        res = [bool(x) for x in evaluate('contains(b, a)')]
        assert_equal(res, python_in)

    def test_str_contains_withemptystr1(self):
        withemptystr = array([b'abc', b'def', b''])
        res = evaluate('contains(b"abcd", withemptystr)')
        assert_equal(res, [True, False, True])

    def test_str_contains_withemptystr2(self):
        withemptystr = array([b'abc', b'def', b''])
        res = evaluate('contains(withemptystr, b"")')
        assert_equal(res, [True, True, True])

    def test_str_contains_long_needle(self):
        a = b'1' + b'a' * 40
        b = b'a' * 40
        res = evaluate('contains(a, b)')
        assert_equal(res, True)

    def test_where_scalar_bool(self):
        a = True
        b = array([1, 2])
        c = array([3, 4])
        res = evaluate('where(a, b, c)')
        assert_array_equal(res, b)
        a = False
        res = evaluate('where(a, b, c)')
        assert_array_equal(res, c)

    @unittest.skipIf(hasattr(sys, "pypy_version_info"),
                     "PyPy does not have sys.getrefcount()")
    def test_refcount(self):
        # Regression test for issue #310
        a = array([1])
        assert sys.getrefcount(a) == 2
        evaluate('1')
        assert sys.getrefcount(a) == 2

    def test_locals_clears_globals(self):
        # Check for issue #313, whereby clearing f_locals also clear f_globals
        # if in the top-frame. This cannot be done inside `unittest` as it is always
        # executing code in a child frame.
        script = r';'.join([
                r"import numexpr as ne",
                r"a=10",
                r"ne.evaluate('1')",
                r"a += 1",
                r"ne.evaluate('2', local_dict={})",
                r"a += 1",
                r"ne.evaluate('3', global_dict={})",
                r"a += 1",
                r"ne.evaluate('4', local_dict={}, global_dict={})",
                r"a += 1",
            ])
        # Raises CalledProcessError on a non-normal exit
        check = subprocess.check_call([sys.executable, '-c', script])
        # Ideally this test should also be done against ipython but it's not
        # a requirement.



class test_numexpr2(test_numexpr):
    """Testing with 2 threads"""
    nthreads = 2


class test_evaluate(TestCase):
    def test_simple(self):
        a = array([1., 2., 3.])
        b = array([4., 5., 6.])
        c = array([7., 8., 9.])
        x = evaluate("2*a + 3*b*c")
        assert_array_equal(x, array([86., 124., 168.]))

    def test_simple_expr_small_array(self):
        x = arange(100.0)
        y = evaluate("x")
        assert_array_equal(x, y)

    def test_simple_expr(self):
        x = arange(1e6)
        y = evaluate("x")
        assert_array_equal(x, y)

    def test_re_evaluate(self):
        a = array([1., 2., 3.])
        b = array([4., 5., 6.])
        c = array([7., 8., 9.])
        x = evaluate("2*a + 3*b*c")
        x = re_evaluate()
        assert_array_equal(x, array([86., 124., 168.]))

    def test_re_evaluate_dict(self):
        a1 = array([1., 2., 3.])
        b1 = array([4., 5., 6.])
        c1 = array([7., 8., 9.])
        local_dict={'a': a1, 'b': b1, 'c': c1}
        x = evaluate("2*a + 3*b*c", local_dict=local_dict)
        x = re_evaluate(local_dict=local_dict)
        assert_array_equal(x, array([86., 124., 168.]))

    def test_validate(self):
        a = array([1., 2., 3.])
        b = array([4., 5., 6.])
        c = array([7., 8., 9.])
        retval = validate("2*a + 3*b*c")
        assert(retval is None)
        x = re_evaluate()
        assert_array_equal(x, array([86., 124., 168.]))

    def test_validate_missing_var(self):
        a = array([1., 2., 3.])
        b = array([4., 5., 6.])
        retval = validate("2*a + 3*b*c")
        assert(isinstance(retval, KeyError))

    def test_validate_syntax(self):
        retval = validate("2+")
        assert(isinstance(retval, SyntaxError))

    def test_validate_dict(self):
        a1 = array([1., 2., 3.])
        b1 = array([4., 5., 6.])
        c1 = array([7., 8., 9.])
        local_dict={'a': a1, 'b': b1, 'c': c1}
        retval = validate("2*a + 3*b*c", local_dict=local_dict)
        assert(retval is None)
        x = re_evaluate(local_dict=local_dict)
        assert_array_equal(x, array([86., 124., 168.]))

    # Test for issue #22
    def test_true_div(self):
        x = arange(10, dtype='i4')
        assert_array_equal(evaluate("x/2"), x / 2)
        assert_array_equal(evaluate("x/2", truediv=False), x / 2)
        assert_array_equal(evaluate("x/2", truediv='auto'), x / 2)
        assert_array_equal(evaluate("x/2", truediv=True), x / 2.0)

    def test_left_shift(self):
        x = arange(10, dtype='i4')
        assert_array_equal(evaluate("x<<2"), x << 2)

    def test_right_shift(self):
        x = arange(10, dtype='i4')
        assert_array_equal(evaluate("x>>2"), x >> 2)

    # PyTables uses __nonzero__ among ExpressionNode objects internally
    # so this should be commented out for the moment.  See #24.
    def test_boolean_operator(self):
        x = arange(10, dtype='i4')
        try:
            evaluate("(x > 1) and (x < 9)")
        except TypeError:
            pass
        else:
            raise ValueError("should raise exception!")

    def test_rational_expr(self):
        a = arange(1e6)
        b = arange(1e6) * 0.1
        x = (a + 2 * b) / (1 + a + 4 * b * b)
        y = evaluate("(a + 2*b) / (1 + a + 4*b*b)")
        assert_array_almost_equal(x, y)

    def test_complex_expr(self):
        def complex(a, b):
            c = zeros(a.shape, dtype=cdouble)
            c.real = a
            c.imag = b
            return c

        a = arange(1e4)
        b = arange(1e4) ** 1e-5
        z = a + 1j * b
        x = z.imag
        x = sin(complex(a, b)).real + z.imag
        y = evaluate("sin(complex(a, b)).real + z.imag")
        assert_array_almost_equal(x, y)

    def test_complex_strides(self):
        a = arange(100).reshape(10, 10)[::2]
        b = arange(50).reshape(5, 10)
        assert_array_equal(evaluate("a+b"), a + b)
        c = empty([10], dtype=[('c1', int32), ('c2', uint16)])
        c['c1'] = arange(10)
        c['c2'].fill(0xaaaa)
        c1 = c['c1']
        a0 = a[0]
        assert_array_equal(evaluate("c1"), c1)
        assert_array_equal(evaluate("a0+c1"), a0 + c1)

    def test_recarray_strides(self):
        a = arange(100)
        b = arange(100,200)
        recarr = np.rec.array(None, formats='f4,f4', shape=(100,))
        recarr['f0'] = a
        recarr['f1'] = b
        c = recarr['f1']
        assert_array_almost_equal(evaluate("sqrt(c) > 1."), sqrt(c) > 1.)
        assert_array_almost_equal(evaluate("log10(c)"), log10(c))

    def test_broadcasting(self):
        a = arange(100).reshape(10, 10)[::2]
        c = arange(10)
        d = arange(5).reshape(5, 1)
        assert_array_equal(evaluate("a+c"), a + c)
        assert_array_equal(evaluate("a+d"), a + d)
        expr = NumExpr("2.0*a+3.0*c", [('a', double), ('c', double)])
        assert_array_equal(expr(a, c), 2.0 * a + 3.0 * c)

    def test_all_scalar(self):
        a = 3.
        b = 4.
        assert_allclose(evaluate("a+b"), a + b)
        expr = NumExpr("2*a+3*b", [('a', double), ('b', double)])
        assert_equal(expr(a, b), 2 * a + 3 * b)

    def test_run(self):
        a = arange(100).reshape(10, 10)[::2]
        b = arange(10)
        expr = NumExpr("2*a+3*b", [('a', double), ('b', double)])
        assert_array_equal(expr(a, b), expr.run(a, b))

    def test_illegal_value(self):
        a = arange(3)
        try:
            evaluate("a < [0, 0, 0]")
        except (ValueError, TypeError):
            pass
        else:
            self.fail()

    def test_sanitize(self):
        with _environment('NUMEXPR_SANITIZE', '1'):
            # Forbid dunder
            try:
                evaluate('__builtins__')
            except ValueError:
                pass
            else:
                self.fail()

            # Forbid colon for lambda funcs
            try:
                evaluate('lambda x: x')
            except ValueError:
                pass
            else:
                self.fail()

            # Forbid indexing
            try:
                evaluate('locals()["evaluate"]')
            except ValueError:
                pass
            else:
                self.fail()

            # Forbid semicolon
            try:
                evaluate('import os;')
            except ValueError:
                pass
            else:
                self.fail()

            # Attribute access with spaces
            try:
                evaluate('os. cpu_count()')
            except ValueError:
                pass
            else:
                self.fail()

            # Attribute access with funny unicode characters that eval translates
            # into ASCII.
            try:
                evaluate("(3+1).ᵇit_length()")
            except ValueError:
                pass
            else:
                self.fail()

            # Pass decimal points including scientific notation
            a = 3.0
            evaluate('a*2.e-5')
            evaluate('a*2.e+5')
            evaluate('a*2e-5')
            evaluate('a*2e+5')
            evaluate('a*2E-5')
            evaluate('a*2.0e5')
            evaluate('a*2.2e5')
            evaluate('2.+a')

            # pass .real and .imag
            c = 2.5 + 1.5j
            evaluate('c.real')
            evaluate('c.imag')

            # pass imaginary unit j
            evaluate('1.5j')
            evaluate('3.j')

            # pass forbidden characters within quotes
            x = np.array(['a', 'b'], dtype=bytes)
            evaluate("x == 'b:'")


    def test_no_sanitize(self):
        try: # Errors on compile() after eval()
            evaluate('import os;', sanitize=False)
        except SyntaxError:
            pass
        else:
            self.fail()

        with _environment('NUMEXPR_SANITIZE', '0'):
            try: # Errors on compile() after eval()
                evaluate('import os;', sanitize=None)
            except SyntaxError:
                pass
            else:
                self.fail()

    def test_disassemble(self):
        assert_equal(disassemble(NumExpr(
            "where(m, a, -1)", [('m', bool), ('a', float)])),
            [[b'where_fbff', b'r0', b'r1[m]', b'r2[a]', b'c3[-1.0]'],
             [b'noop', None, None, None]])

    def test_constant_deduplication(self):
        assert_equal(NumExpr("(a + 1)*(a - 1)", [('a', np.int32)]).constants, (1,))

    def test_nan_constant(self):
        assert_equal(str(ConstantNode(float("nan")).value), 'nan')

        # check de-duplication works for nan
        _nan = ConstantNode(float("nan"))
        expr = (E.a + _nan)*(E.b + _nan)
        assert_equal(NumExpr(expr, [('a', double), ('b', double)]).constants, (float("nan"),))


    def test_f32_constant(self):
        assert_equal(ConstantNode(numpy.float32(1)).astKind, "float")
        assert_equal(ConstantNode(numpy.float32("nan")).astKind, "float")
        assert_equal(ConstantNode(numpy.float32(3)).value.dtype, numpy.dtype("float32"))
        assert_array_equal(NumExpr(ConstantNode(numpy.float32(1))).run(),
                           numpy.array(1, dtype="float32"))

    def test_unaligned_singleton(self):
        # Test for issue #397 whether singletons outputs assigned to consts must be
        # aligned or not.
        a = np.empty(5, dtype=np.uint8)[1:].view(np.int32)
        evaluate('3', out=a)
        assert_equal(a, 3)

    def test_negative_mod(self):
        # Test for issue #413, modulus of negative integers. C modulus is
        # actually remainder op, and hence different from Python modulus.
        a = np.array([-500, -135, 0, 0, 135, 500], dtype=np.int32)
        n = np.array([-360, -360, -360, 360, 360, 360], dtype=np.int32)
        out_i = evaluate('a % n')
        assert_equal(out_i, np.mod(a, n))

        b = a.astype(np.int64)
        m = n.astype(np.int64)
        out_l = evaluate('b % m')
        assert_equal(out_l, np.mod(b, m))

    def test_negative_power_scalar(self):
        # Test for issue #428, where the power is negative and the base is an
        # integer. This was running afoul in the precomputation in `expressions.py:pow_op()`
        base = np.array([-2, -1, 1, 2, 3], dtype=np.int32)
        out_i = evaluate('base ** -1.0')
        assert_equal(out_i, np.power(base, -1.0))

        base = np.array([-2, -1, 1, 2, 3], dtype=np.int64)
        out_l = evaluate('base ** -1.0')
        assert_equal(out_l, np.power(base, -1.0))


    def test_ex_uses_vml(self):
        vml_funcs = [ "sin", "cos", "tan", "arcsin", "arccos", "arctan",
                      "sinh", "cosh", "tanh", "arcsinh", "arccosh", "arctanh",
                      "log", "log1p","log10", "exp", "expm1", "abs", "conj",
                      "arctan2", "fmod"]
        for func in vml_funcs:
            strexpr = func+'(a)'
            _, ex_uses_vml = numexpr.necompiler.getExprNames(strexpr, {})
            assert_equal(ex_uses_vml, use_vml, strexpr)

    if 'sparc' not in platform.machine():
        # Execution order set here so as to not use too many threads
        # during the rest of the execution.  See #33 for details.
        def test_changing_nthreads_00_inc(self):
            a = linspace(-1, 1, 1000000)
            b = ((.25 * a + .75) * a - 1.5) * a - 2
            for nthreads in range(1, 7):
                numexpr.set_num_threads(nthreads)
                c = evaluate("((.25*a + .75)*a - 1.5)*a - 2")
                assert_array_almost_equal(b, c)

        def test_changing_nthreads_01_dec(self):
            a = linspace(-1, 1, 1000000)
            b = ((.25 * a + .75) * a - 1.5) * a - 2
            for nthreads in range(6, 1, -1):
                numexpr.set_num_threads(nthreads)
                c = evaluate("((.25*a + .75)*a - 1.5)*a - 2")
                assert_array_almost_equal(b, c)


tests = [
    ('MISC', ['b*c+d*e',
              '2*a+3*b',
              '-a',
              'sinh(a)',
              '2*a + (cos(3)+5)*sinh(cos(b))',
              '2*a + arctan2(a, b)',
              'arcsin(0.5)',
              'where(a != 0.0, 2, a)',
              'where(a > 10, b < a, b > a)',
              'where((a-10).real != 0.0, a, 2)',
              '0.25 * (a < 5) + 0.33 * (a >= 5)',
              'cos(1+1)',
              '1+1',
              '1',
              'cos(a2)',
    ])]

optests = []
for op in list('+-*/%') + ['**']:
    optests.append("(a+1) %s (b+3)" % op)
    optests.append("3 %s (b+3)" % op)
    optests.append("(a+1) %s 4" % op)
    optests.append("2 %s (b+3)" % op)
    optests.append("(a+1) %s 2" % op)
    optests.append("(a+1) %s -1" % op)
    optests.append("(a+1) %s 0.5" % op)
    # Check divisions and modulus by zero (see ticket #107)
    optests.append("(a+1) %s 0" % op)
tests.append(('OPERATIONS', optests))

cmptests = []
for op in ['<', '<=', '==', '>=', '>', '!=']:
    cmptests.append("a/2+5 %s b" % op)
    cmptests.append("a/2+5 %s 7" % op)
    cmptests.append("7 %s b" % op)
    cmptests.append("7.0 %s 5" % op)
tests.append(('COMPARISONS', cmptests))

func1tests = []
for func in ['copy', 'ones_like', 'sqrt',
             'sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan',
             'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh',
             'log', 'log1p', 'log10', 'exp', 'expm1', 'abs', 'conj',
             'ceil', 'floor']:
    func1tests.append("a + %s(b+c)" % func)
tests.append(('1_ARG_FUNCS', func1tests))

func2tests = []
for func in ['arctan2', 'fmod']:
    func2tests.append("a + %s(b+c, d+1)" % func)
    func2tests.append("a + %s(b+c, 1)" % func)
    func2tests.append("a + %s(1, d+1)" % func)
tests.append(('2_ARG_FUNCS', func2tests))

powtests = []
# n = -1, 0.5, 2, 4 already handled in section "OPERATIONS"
for n in (-7, -2.5, -1.5, -1.3, -.5, 0, 0.0, 1, 2.3, 2.5, 3):
    powtests.append("(a+1)**%s" % n)
tests.append(('POW_TESTS', powtests))


def equal(a, b, exact):
    if array_equal(a, b):
        return True

    if hasattr(a, 'dtype') and a.dtype in ['f4', 'f8']:
        nnans = isnan(a).sum()
        if nnans > 0:
            # For results containing NaNs, just check that the number
            # of NaNs is the same in both arrays.  This check could be
            # made more exhaustive, but checking element by element in
            # python space is very expensive in general.
            return nnans == isnan(b).sum()
        ninfs = isinf(a).sum()
        if ninfs > 0:
            # Ditto for Inf's
            return ninfs == isinf(b).sum()
    if exact:
        return (shape(a) == shape(b)) and alltrue(ravel(a) == ravel(b), axis=0)
    else:
        if hasattr(a, 'dtype') and a.dtype == 'f4':
            atol = 1e-5  # Relax precision for special opcodes, like fmod
        else:
            atol = 1e-8
        return (shape(a) == shape(b) and
                allclose(ravel(a), ravel(b), atol=atol))


class Skip(Exception): pass


def test_expressions():
    test_no = [0]

    def make_test_method(a, a2, b, c, d, e, x, expr,
                         test_scalar, dtype, optimization, exact, section):
        this_locals = locals()

        def method():
            try:
                # We don't want to listen at RuntimeWarnings like
                # "overflows" or "divide by zero" in plain eval().
                warnings.simplefilter("ignore")
                npval = eval(expr, globals(), this_locals)
                warnings.simplefilter("always")
                npval = eval(expr, globals(), this_locals)
            except Exception as ex:
                # just store the exception in a variable
                # compatibility with numpy v1.12
                # see also https://github.com/pydata/numexpr/issues/239
                np_exception = ex
                npval = None
            else:
                np_exception = None

            try:
                neval = evaluate(expr, local_dict=this_locals,
                                 optimization=optimization)
            except AssertionError:
                raise
            except NotImplementedError:
                print('%r not implemented for %s (scalar=%d, opt=%s)'
                      % (expr, dtype.__name__, test_scalar, optimization))
            except Exception as ne_exception:
                same_exc_type = issubclass(type(ne_exception),
                                           type(np_exception))
                if np_exception is None or not same_exc_type:
                    print('numexpr error for expression %r' % (expr,))
                    raise
            except:
                print('numexpr error for expression %r' % (expr,))
                raise
            else:
                msg = ('expected numexpr error not raised for expression '
                       '%r' % (expr,))
                assert np_exception is None, msg

                assert equal(npval, neval, exact), """%r
(test_scalar=%r, dtype=%r, optimization=%r, exact=%r,
 npval=%r (%r - %r)\n neval=%r (%r - %r))""" % (expr, test_scalar, dtype.__name__,
                                                optimization, exact,
                                                npval, type(npval), shape(npval),
                                                neval, type(neval), shape(neval))

        method.description = ('test_expressions(%s, test_scalar=%r, '
                              'dtype=%r, optimization=%r, exact=%r)') % (expr, test_scalar, dtype.__name__, optimization, exact)
        test_no[0] += 1
        method.__name__ = 'test_scalar%d_%s_%s_%s_%04d' % (test_scalar,
                                                           dtype.__name__,
                                                           optimization.encode('ascii'),
                                                           section.encode('ascii'),
                                                           test_no[0])
        return method

    x = None
    for test_scalar in (0, 1, 2):
        for dtype in (int, int, np.float32, double, complex):
            array_size = 100
            a = arange(2 * array_size, dtype=dtype)[::2]
            a2 = zeros([array_size, array_size], dtype=dtype)
            b = arange(array_size, dtype=dtype) / array_size
            c = arange(array_size, dtype=dtype)
            d = arange(array_size, dtype=dtype)
            e = arange(array_size, dtype=dtype)
            if dtype == complex:
                a = a.real
                for x in [a2, b, c, d, e]:
                    x += 1j
                    x *= 1 + 1j
            if test_scalar == 1:
                a = a[array_size // 2]
            if test_scalar == 2:
                b = b[array_size // 2]
            for optimization, exact in [
                ('none', False), ('moderate', False), ('aggressive', False)]:
                for section_name, section_tests in tests:
                    for expr in section_tests:
                        if (dtype == complex and
                            ('<' in expr or '>' in expr or '%' in expr
                             or "arctan2" in expr or "fmod" in expr
                             or "floor" in expr or "ceil" in expr)):
                            # skip complex comparisons or functions not
                            # defined in complex domain.
                            continue
                        if (dtype in (int, int) and test_scalar and
                                    expr == '(a+1) ** -1'):
                            continue

                        m = make_test_method(a, a2, b, c, d, e, x,
                                             expr, test_scalar, dtype,
                                             optimization, exact,
                                             section_name)
                        yield m


class test_int64(TestCase):
    def test_neg(self):
        a = array([2 ** 31 - 1, 2 ** 31, 2 ** 32, 2 ** 63 - 1], dtype=int64)
        res = evaluate('-a')
        assert_array_equal(res, [1 - 2 ** 31, -(2 ** 31), -(2 ** 32), 1 - 2 ** 63])
        self.assertEqual(res.dtype.name, 'int64')


class test_int32_int64(TestCase):

    def test_small_int(self):
        # Small ints (32-bit ones) should not be promoted to longs.
        res = evaluate('2')
        assert_array_equal(res, 2)
        self.assertEqual(res.dtype.name, 'int32')

    def test_big_int(self):
        # Big ints should be promoted to longs.
        res = evaluate('2**40')
        assert_array_equal(res, 2 ** 40)
        self.assertEqual(res.dtype.name, 'int64')

    def test_long_constant_promotion(self):
        int32array = arange(100, dtype='int32')
        itwo = np.int32(2)
        ltwo = np.int64(2)
        res = int32array * 2
        res32 = evaluate('int32array * itwo')
        res64 = evaluate('int32array * ltwo')
        assert_array_equal(res, res32)
        assert_array_equal(res, res64)
        self.assertEqual(res32.dtype.name, 'int32')
        self.assertEqual(res64.dtype.name, 'int64')

    def test_int64_array_promotion(self):
        int32array = arange(100, dtype='int32')
        int64array = arange(100, dtype='int64')
        respy = int32array * int64array
        resnx = evaluate('int32array * int64array')
        assert_array_equal(respy, resnx)
        self.assertEqual(resnx.dtype.name, 'int64')


class test_uint32_int64(TestCase):
    def test_small_uint32(self):
        # Small uint32 should not be downgraded to ints.
        a = np.uint32(42)
        res = evaluate('a')
        assert_array_equal(res, 42)
        self.assertEqual(res.dtype.name, 'int64')

    def test_uint32_constant_promotion(self):
        int32array = arange(100, dtype='int32')
        stwo = np.int32(2)
        utwo = np.uint32(2)
        res = int32array * utwo
        res32 = evaluate('int32array * stwo')
        res64 = evaluate('int32array * utwo')
        assert_array_equal(res, res32)
        assert_array_equal(res, res64)
        self.assertEqual(res32.dtype.name, 'int32')
        self.assertEqual(res64.dtype.name, 'int64')

    def test_int64_array_promotion(self):
        uint32array = arange(100, dtype='uint32')
        int64array = arange(100, dtype='int64')
        respy = uint32array * int64array
        resnx = evaluate('uint32array * int64array')
        assert_array_equal(respy, resnx)
        self.assertEqual(resnx.dtype.name, 'int64')


class test_strings(TestCase):
    BLOCK_SIZE1 = 128
    BLOCK_SIZE2 = 8
    str_list1 = [b'foo', b'bar', b'', b'  ']
    str_list2 = [b'foo', b'', b'x', b' ']
    str_nloops = len(str_list1) * (BLOCK_SIZE1 + BLOCK_SIZE2 + 1)
    str_array1 = array(str_list1 * str_nloops)
    str_array2 = array(str_list2 * str_nloops)
    str_constant = b'doodoo'

    def test_null_chars(self):
        str_list = [
            b'\0\0\0', b'\0\0foo\0', b'\0\0foo\0b', b'\0\0foo\0b\0',
            b'foo\0', b'foo\0b', b'foo\0b\0', b'foo\0bar\0baz\0\0']
        for s in str_list:
            r = evaluate('s')
            self.assertEqual(s, r.tobytes())  # check *all* stored data

    def test_compare_copy(self):
        sarr = self.str_array1
        expr = 'sarr'
        res1 = eval(expr)
        res2 = evaluate(expr)
        assert_array_equal(res1, res2)

    def test_compare_array(self):
        sarr1 = self.str_array1
        sarr2 = self.str_array2
        expr = 'sarr1 >= sarr2'
        res1 = eval(expr)
        res2 = evaluate(expr)
        assert_array_equal(res1, res2)

    def test_compare_variable(self):
        sarr = self.str_array1
        svar = self.str_constant
        expr = 'sarr >= svar'
        res1 = eval(expr)
        res2 = evaluate(expr)
        assert_array_equal(res1, res2)

    def test_compare_constant(self):
        sarr = self.str_array1
        expr = 'sarr >= %r' % self.str_constant
        res1 = eval(expr)
        res2 = evaluate(expr)
        assert_array_equal(res1, res2)

    def test_add_string_array(self):
        sarr1 = self.str_array1
        sarr2 = self.str_array2
        expr = 'sarr1 + sarr2'
        self.assert_missing_op('add_sss', expr, locals())

    def test_empty_string1(self):
        a = np.array([b"", b"pepe"])
        b = np.array([b"pepe2", b""])
        res = evaluate("(a == b'') & (b == b'pepe2')")
        assert_array_equal(res, np.array([True, False]))
        res2 = evaluate("(a == b'pepe') & (b == b'')")
        assert_array_equal(res2, np.array([False, True]))

    def test_empty_string2(self):
        a = np.array([b"p", b"pepe"])
        b = np.array([b"pepe2", b""])
        res = evaluate("(a == b'') & (b == b'pepe2')")
        assert_array_equal(res, np.array([False, False]))
        res2 = evaluate("(a == b'pepe') & (b == b'')")
        assert_array_equal(res, np.array([False, False]))

    def test_add_numeric_array(self):
        sarr = self.str_array1
        narr = arange(len(sarr), dtype='int32')
        expr = 'sarr >= narr'
        self.assert_missing_op('ge_bsi', expr, locals())

    def assert_missing_op(self, op, expr, local_dict):
        msg = "expected NotImplementedError regarding '%s'" % op
        try:
            evaluate(expr, local_dict)
        except NotImplementedError as nie:
            if "'%s'" % op not in nie.args[0]:
                self.fail(msg)
        else:
            self.fail(msg)

    def test_compare_prefix(self):
        # Check comparing two strings where one is a prefix of the
        # other.
        for s1, s2 in [(b'foo', b'foobar'), (b'foo', b'foo\0bar'),
                       (b'foo\0a', b'foo\0bar')]:
            self.assertTrue(evaluate('s1 < s2'))
            self.assertTrue(evaluate('s1 <= s2'))
            self.assertTrue(evaluate('~(s1 == s2)'))
            self.assertTrue(evaluate('~(s1 >= s2)'))
            self.assertTrue(evaluate('~(s1 > s2)'))

        # Check for NumPy array-style semantics in string equality.
        s1, s2 = b'foo', b'foo\0\0'
        self.assertTrue(evaluate('s1 == s2'))


# Case for testing selections in fields which are aligned but whose
# data length is not an exact multiple of the length of the record.
# The following test exposes the problem only in 32-bit machines,
# because in 64-bit machines 'c2' is unaligned.  However, this should
# check most platforms where, while not unaligned, 'len(datatype) >
# boundary_alignment' is fullfilled.
class test_irregular_stride(TestCase):
    def test_select(self):
        f0 = arange(10, dtype=int32)
        f1 = arange(10, dtype=float64)

        irregular = rec.fromarrays([f0, f1])

        f0 = irregular['f0']
        f1 = irregular['f1']

        i0 = evaluate('f0 < 5')
        i1 = evaluate('f1 < 5')

        assert_array_equal(f0[i0], arange(5, dtype=int32))
        assert_array_equal(f1[i1], arange(5, dtype=float64))


# Cases for testing arrays with dimensions that can be zero.
class test_zerodim(TestCase):
    def test_zerodim1d(self):
        a0 = array([], dtype=int32)
        a1 = array([], dtype=float64)

        r0 = evaluate('a0 + a1')
        r1 = evaluate('a0 * a1')

        assert_array_equal(r0, a1)
        assert_array_equal(r1, a1)

    def test_zerodim3d(self):
        a0 = array([], dtype=int32).reshape(0, 2, 4)
        a1 = array([], dtype=float64).reshape(0, 2, 4)

        r0 = evaluate('a0 + a1')
        r1 = evaluate('a0 * a1')

        assert_array_equal(r0, a1)
        assert_array_equal(r1, a1)


@contextmanager
def _environment(key, value):
    old = os.environ.get(key)
    os.environ[key] = value
    try:
        yield
    finally:
        if old:
            os.environ[key] = old
        else:
            del os.environ[key]

# Test cases for the threading configuration
class test_threading_config(TestCase):
    def test_max_threads_unset(self):
        # Has to be done in a subprocess as `importlib.reload` doesn't let us
        # re-initialize the threadpool
        script = '\n'.join([
                "import os",
                "if 'NUMEXPR_MAX_THREADS' in os.environ: os.environ.pop('NUMEXPR_MAX_THREADS')",
                "if 'OMP_NUM_THREADS' in os.environ: os.environ.pop('OMP_NUM_THREADS')",
                "import numexpr",
                f"assert(numexpr.nthreads <= {MAX_THREADS})",
                "exit(0)"])
        subprocess.check_call([sys.executable, '-c', script])

    def test_max_threads_set(self):
        # Has to be done in a subprocess as `importlib.reload` doesn't let us
        # re-initialize the threadpool
        script = '\n'.join([
                "import os",
                "os.environ['NUMEXPR_MAX_THREADS'] = '4'",
                "import numexpr",
                "assert(numexpr.MAX_THREADS == 4)",
                "exit(0)"])
        subprocess.check_call([sys.executable, '-c', script])

    def test_numexpr_num_threads(self):
        with _environment('OMP_NUM_THREADS', '5'):
            # NUMEXPR_NUM_THREADS has priority
            with _environment('NUMEXPR_NUM_THREADS', '3'):
                 if 'sparc' in platform.machine():
                     self.assertEqual(1, numexpr._init_num_threads())
                 else:
                     self.assertEqual(3, numexpr._init_num_threads())

    def test_omp_num_threads(self):
        with _environment('OMP_NUM_THREADS', '5'):
            if 'sparc' in platform.machine():
                self.assertEqual(1, numexpr._init_num_threads())
            else:
                self.assertEqual(5, numexpr._init_num_threads())

    def test_omp_num_threads_empty_string(self):
        with _environment('OMP_NUM_THREADS', ''):
            if 'sparc' in platform.machine():
                self.assertEqual(1, numexpr._init_num_threads())
            else:
                self.assertEqual(min(detect_number_of_cores(), MAX_THREADS), numexpr._init_num_threads())

    def test_numexpr_max_threads_empty_string(self):
        with _environment('NUMEXPR_MAX_THREADS', ''):
            if 'sparc' in platform.machine():
                self.assertEqual(1, numexpr._init_num_threads())
            else:
                self.assertEqual(min(detect_number_of_cores(), MAX_THREADS), numexpr._init_num_threads())

    def test_vml_threads_round_trip(self):
        n_threads = 3
        if use_vml:
            numexpr.utils.set_vml_num_threads(n_threads)
            set_threads = numexpr.utils.get_vml_num_threads()
            self.assertEqual(n_threads, set_threads)
        else:
            self.assertIsNone(numexpr.utils.set_vml_num_threads(n_threads))
            self.assertIsNone(numexpr.utils.get_vml_num_threads())


# Case test for threads
class test_threading(TestCase):

    def test_thread(self):
        import threading

        class ThreadTest(threading.Thread):
            def run(self):
                a = arange(3)
                assert_array_equal(evaluate('a**3'), array([0, 1, 8]))

        test = ThreadTest()
        test.start()
        test.join()

    def test_multithread(self):

        import threading

        # Running evaluate() from multiple threads shouldn't crash
        def work(n):
            a = arange(n)
            evaluate('a+a')

        work(10)  # warm compilation cache

        nthreads = 30
        threads = [threading.Thread(target=work, args=(1e5,))
                   for i in range(nthreads)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

    def test_thread_safety(self):
        """
        Expected output

        When not safe (before the pr this test is commited)
        AssertionError: Thread-0 failed: result does not match expected

        When safe (after the pr this test is commited)
        Should pass without failure
        """
        import threading
        import time

        barrier = threading.Barrier(4)

        # Function that each thread will run with different expressions
        def thread_function(a_value, b_value, expression, expected_result, results, index):
            validate(expression, local_dict={"a": a_value, "b": b_value})
            # Wait for all threads to reach this point
            # such that they all set _numexpr_last
            barrier.wait()

            # Simulate some work or a context switch delay
            time.sleep(0.1)

            result = re_evaluate(local_dict={"a": a_value, "b": b_value})
            results[index] = np.array_equal(result, expected_result)

        def test_thread_safety_with_numexpr():
            num_threads = 4
            array_size = 1000000

            expressions = [
                "a + b",
                "a - b",
                "a * b",
                "a / b"
            ]

            a_value = [np.full(array_size, i + 1) for i in range(num_threads)]
            b_value = [np.full(array_size, (i + 1) * 2) for i in range(num_threads)]

            expected_results = [
                a_value[i] + b_value[i] if expr == "a + b" else
                a_value[i] - b_value[i] if expr == "a - b" else
                a_value[i] * b_value[i] if expr == "a * b" else
                a_value[i] / b_value[i] if expr == "a / b" else None
                for i, expr in enumerate(expressions)
            ]

            results = [None] * num_threads
            threads = []

            # Create and start threads with different expressions
            for i in range(num_threads):
                thread = threading.Thread(
                    target=thread_function,
                    args=(a_value[i], b_value[i], expressions[i], expected_results[i], results, i)
                )
                threads.append(thread)
                thread.start()

            for thread in threads:
                thread.join()

            for i in range(num_threads):
                if not results[i]:
                    self.fail(f"Thread-{i} failed: result does not match expected")

        test_thread_safety_with_numexpr()


# The worker function for the subprocess (needs to be here because Windows
# has problems pickling nested functions with the multiprocess module :-/)
def _worker(qout=None):
    ra = np.arange(1e3)
    rows = evaluate('ra > 0')
    #print "Succeeded in evaluation!\n"
    if qout is not None:
        qout.put("Done")


# Case test for subprocesses (via multiprocessing module)
class test_subprocess(TestCase):
    def test_multiprocess(self):
        try:
            import multiprocessing as mp
        except ImportError:
            return
        # Check for two threads at least
        numexpr.set_num_threads(2)
        #print "**** Running from main process:"
        _worker()
        #print "**** Running from subprocess:"
        qout = mp.Queue()
        ps = mp.Process(target=_worker, args=(qout,))
        ps.daemon = True
        ps.start()

        result = qout.get()
        #print result


def print_versions():
    """Print the versions of software that numexpr relies on."""
    # from pkg_resources import parse_version
    from numexpr.cpuinfo import cpu
    import platform

    print('-=' * 38)
    print('Numexpr version:   %s' % numexpr.__version__)
    print('NumPy version:     %s' % np.__version__)
    print('Python version:    %s' % sys.version)
    (sysname, nodename, release, os_version, machine, processor) = platform.uname()
    print('Platform:          %s-%s-%s' % (sys.platform, machine, os_version))
    try:
        # cpuinfo doesn't work on OSX well it seems, so protect these outputs
        # with a try block
        cpu_info = cpu.info[0]
        print('CPU vendor:        %s' % cpu_info.get('VendorIdentifier', ''))
        print('CPU model:         %s' % cpu_info.get('ProcessorNameString', ''))
        print('CPU clock speed:   %s MHz' % cpu_info.get('~MHz',''))
    except KeyError:
        pass
    print('VML available?     %s' % use_vml)
    if use_vml:
        print('VML/MKL version:   %s' % numexpr.get_vml_version())
    print('Number of threads used by default: %d '
          '(out of %d detected cores)' % (numexpr.nthreads, numexpr.ncores))
    print('Maximum number of threads: %s' % numexpr.MAX_THREADS)
    print('-=' * 38)


def test(verbosity=1):
    """
    Run all the tests in the test suite.
    """
    print_versions()
    # For some reason, NumPy issues all kinds of warnings when using Python3.
    # Ignoring them in tests should be ok, as all results are checked out.
    # See https://github.com/pydata/numexpr/issues/183 for details.
    np.seterr(divide='ignore', invalid='ignore', over='ignore', under='ignore')
    return unittest.TextTestRunner(verbosity=verbosity).run(suite())


test.__test__ = False


def suite():
    import unittest
    import platform as pl

    theSuite = unittest.TestSuite()
    niter = 1

    class TestExpressions(TestCase):
        pass

    def add_method(func):
        def method(self):
            return func()

        setattr(TestExpressions, func.__name__,
                method.__get__(None, TestExpressions))

    for func in test_expressions():
        add_method(func)

    for n in range(niter):
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_numexpr))
        if 'sparc' not in platform.machine():
            theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_numexpr2))
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_evaluate))
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestExpressions))
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_int32_int64))
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_uint32_int64))
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_strings))
        theSuite.addTest(
            unittest.defaultTestLoader.loadTestsFromTestCase(test_irregular_stride))
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_zerodim))
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_threading_config))

        # multiprocessing module is not supported on Hurd/kFreeBSD
        if (pl.system().lower() not in ('gnu', 'gnu/kfreebsd')):
            theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_subprocess))

        # I need to put this test after test_subprocess because
        # if not, the test suite locks immediately before test_subproces.
        # This only happens with Windows, so I suspect of a subtle bad
        # interaction with threads and subprocess :-/
        theSuite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_threading))

    return theSuite


if __name__ == '__main__':
    print_versions()
    unittest.main(defaultTest='suite')
#    suite = suite()
#    unittest.TextTestRunner(verbosity=2).run(suite)
