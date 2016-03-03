# -*- coding: utf-8 -*-
"""
Testing that functions from compat work as expected
"""

from pandas.compat import (range, zip, map, filter, lrange, lzip, lmap,
                           lfilter, builtins, iterkeys, itervalues, iteritems,
                           next)
import pandas.util.testing as tm


class TestBuiltinIterators(tm.TestCase):

    def check_result(self, actual, expected, lengths):
        for (iter_res, list_res), exp, length in zip(actual, expected,
                                                     lengths):
            self.assertNotIsInstance(iter_res, list)
            tm.assertIsInstance(list_res, list)
            iter_res = list(iter_res)
            self.assertEqual(len(list_res), length)
            self.assertEqual(len(iter_res), length)
            self.assertEqual(iter_res, exp)
            self.assertEqual(list_res, exp)

    def test_range(self):
        actual1 = range(10)
        actual2 = lrange(10)
        actual = [actual1, actual2],
        expected = list(builtins.range(10)),
        lengths = 10,

        actual1 = range(1, 10, 2)
        actual2 = lrange(1, 10, 2)
        actual += [actual1, actual2],
        lengths += 5,
        expected += list(builtins.range(1, 10, 2)),
        self.check_result(actual, expected, lengths)

    def test_map(self):
        func = lambda x, y, z: x + y + z
        lst = [builtins.range(10), builtins.range(10), builtins.range(10)]
        actual1 = map(func, *lst)
        actual2 = lmap(func, *lst)
        actual = [actual1, actual2],
        expected = list(builtins.map(func, *lst)),
        lengths = 10,
        self.check_result(actual, expected, lengths)

    def test_filter(self):
        func = lambda x: x
        lst = list(builtins.range(10))
        actual1 = filter(func, lst)
        actual2 = lfilter(func, lst)
        actual = [actual1, actual2],
        lengths = 9,
        expected = list(builtins.filter(func, lst)),
        self.check_result(actual, expected, lengths)

    def test_zip(self):
        lst = [builtins.range(10), builtins.range(10), builtins.range(10)]
        actual = [zip(*lst), lzip(*lst)],
        expected = list(builtins.zip(*lst)),
        lengths = 10,
        self.check_result(actual, expected, lengths)

    def test_dict_iterators(self):
        self.assertEqual(next(itervalues({1: 2})), 2)
        self.assertEqual(next(iterkeys({1: 2})), 1)
        self.assertEqual(next(iteritems({1: 2})), (1, 2))
