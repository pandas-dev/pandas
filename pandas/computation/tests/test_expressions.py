# pylint: disable-msg=W0612,E1101

import unittest
import operator

import nose


import numpy as np
from numpy.testing import assert_array_equal

from pandas.core.api import DataFrame
from pandas.computation import expressions as expr

if not expr._USE_NUMEXPR:
    raise nose.SkipTest

import numexpr as ne


_frame = DataFrame(np.random.randn(10000, 4), columns=list('ABCD'),
                   dtype='float64')
_frame2 = DataFrame(np.random.randn(100, 4), columns=list('ABCD'),
                    dtype='float64')
_mixed = DataFrame({'A': _frame['A'].copy(),
                    'B': _frame['B'].astype('float32'),
                    'C': _frame['C'].astype('int64'),
                    'D': _frame['D'].astype('int32')})
_mixed2 = DataFrame({'A': _frame2['A'].copy(),
                     'B': _frame2['B'].astype('float32'),
                     'C': _frame2['C'].astype('int64'),
                     'D': _frame2['D'].astype('int32')})


class TestExpressions(unittest.TestCase):

    _multiprocess_can_split_ = False

    def setUp(self):
        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.mixed = _mixed.copy()
        self.mixed2 = _mixed2.copy()

    def test_invalid(self):
        # no op
        result = expr._can_use_numexpr(operator.add, None, self.frame,
                                       self.frame, 'evaluate')
        self.assertFalse(result)

        # mixed
        result = expr._can_use_numexpr(
            operator.add, '+', self.mixed, self.frame, 'evaluate')
        self.assertFalse(result)

        # min elements
        result = expr._can_use_numexpr(
            operator.add, '+', self.frame2, self.frame2, 'evaluate')
        self.assertFalse(result)

        # ok, we only check on first part of expression
        result = expr._can_use_numexpr(
            operator.add, '+', self.frame, self.frame2, 'evaluate')
        self.assert_(result)

    def test_binary_ops(self):
        def testit():

            for f, f2 in [(self.frame, self.frame2),
                          (self.mixed, self.mixed2)]:

                for op, op_str in [('add', '+'), ('sub', '-'), ('mul', '*'),
                                   ('div', '/'), ('pow', '**')]:

                    op = getattr(operator, op, None)
                    if op is not None:
                        result = expr._can_use_numexpr(
                            op, op_str, f, f, 'evaluate')
                        self.assert_(result == (not f._is_mixed_type))

                        result = expr.evaluate(
                            op, op_str, f, f, use_numexpr=True)
                        expected = expr.evaluate(
                            op, op_str, f, f, use_numexpr=False)
                        assert_array_equal(result, expected.values)

                        result = expr._can_use_numexpr(
                            op, op_str, f2, f2, 'evaluate')
                        self.assertFalse(result)

        expr.set_use_numexpr(False)
        testit()
        expr.set_use_numexpr(True)
        expr.set_numexpr_threads(1)
        testit()
        expr.set_numexpr_threads()
        testit()

    def test_boolean_ops(self):
        def testit():
            for f, f2 in [(self.frame, self.frame2),
                          (self.mixed, self.mixed2)]:

                f11 = f
                f12 = f + 1

                f21 = f2
                f22 = f2 + 1

                for op, op_str in [('gt', '>'), ('lt', '<'), ('ge', '>='),
                                   ('le', '<='), ('eq', '=='), ('ne', '!=')]:

                    op = getattr(operator, op)

                    result = expr._can_use_numexpr(
                        op, op_str, f11, f12, 'evaluate')
                    self.assert_(result == (not f11._is_mixed_type))

                    result = expr.evaluate(
                        op, op_str, f11, f12, use_numexpr=True)
                    expected = expr.evaluate(
                        op, op_str, f11, f12, use_numexpr=False)
                    assert_array_equal(result, expected.values)

                    result = expr._can_use_numexpr(
                        op, op_str, f21, f22, 'evaluate')
                    self.assertFalse(result)

        expr.set_use_numexpr(False)
        testit()
        expr.set_use_numexpr(True)
        expr.set_numexpr_threads(1)
        testit()
        expr.set_numexpr_threads()
        testit()

    def test_where(self):
        def testit():
            for f in [self.frame, self.frame2, self.mixed, self.mixed2]:

                for cond in [True, False]:

                    c = np.empty(f.shape, dtype=np.bool_)
                    c.fill(cond)
                    result = expr.where(c, f.values, f.values + 1)
                    expected = np.where(c, f.values, f.values + 1)
                    assert_array_equal(result, expected)

        expr.set_use_numexpr(False)
        testit()
        expr.set_use_numexpr(True)
        expr.set_numexpr_threads(1)
        testit()
        expr.set_numexpr_threads()
        testit()


