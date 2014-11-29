# -*- coding: utf-8 -*-
# pylint: disable=E1101,E1103,W0232

from datetime import datetime
from pandas.compat import range, lrange, u
import os
import pickle
import re
from distutils.version import LooseVersion

import numpy as np
import pandas as pd

from pandas import Categorical, Index, Series, DataFrame, PeriodIndex, Timestamp

import pandas.core.common as com
import pandas.compat as compat
import pandas.util.testing as tm

class TestCategorical(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.factor = Categorical.from_array(['a', 'b', 'b', 'a',
                                              'a', 'c', 'c', 'c'])

    def test_getitem(self):
        self.assertEqual(self.factor[0], 'a')
        self.assertEqual(self.factor[-1], 'c')

        subf = self.factor[[0, 1, 2]]
        tm.assert_almost_equal(subf._codes, [0, 1, 1])

        subf = self.factor[np.asarray(self.factor) == 'c']
        tm.assert_almost_equal(subf._codes, [2, 2, 2])

    def test_constructor_unsortable(self):

        # it works!
        arr = np.array([1, 2, 3, datetime.now()], dtype='O')
        factor = Categorical.from_array(arr)
        self.assertFalse(factor.ordered)

    def test_constructor(self):

        exp_arr = np.array(["a", "b", "c", "a", "b", "c"])
        c1 = Categorical(exp_arr)
        self.assert_numpy_array_equal(c1.__array__(), exp_arr)
        c2 = Categorical(exp_arr, categories=["a","b","c"])
        self.assert_numpy_array_equal(c2.__array__(), exp_arr)
        c2 = Categorical(exp_arr, categories=["c","b","a"])
        self.assert_numpy_array_equal(c2.__array__(), exp_arr)

        # categories must be unique
        def f():
            Categorical([1,2], [1,2,2])
        self.assertRaises(ValueError, f)
        def f():
            Categorical(["a","b"], ["a","b","b"])
        self.assertRaises(ValueError, f)
        def f():
            Categorical([1,2], [1,2,np.nan, np.nan])
        self.assertRaises(ValueError, f)


        # Categorical as input
        c1 = Categorical(["a", "b", "c", "a"])
        c2 = Categorical(c1)
        self.assertTrue(c1.equals(c2))

        c1 = Categorical(["a", "b", "c", "a"], categories=["a","b","c","d"])
        c2 = Categorical(c1)
        self.assertTrue(c1.equals(c2))

        c1 = Categorical(["a", "b", "c", "a"], categories=["a","c","b"])
        c2 = Categorical(c1)
        self.assertTrue(c1.equals(c2))

        c1 = Categorical(["a", "b", "c", "a"], categories=["a","c","b"])
        c2 = Categorical(c1, categories=["a","b","c"])
        self.assert_numpy_array_equal(c1.__array__(), c2.__array__())
        self.assert_numpy_array_equal(c2.categories, np.array(["a","b","c"]))

        # Series of dtype category
        c1 = Categorical(["a", "b", "c", "a"], categories=["a","b","c","d"])
        c2 = Categorical(Series(c1))
        self.assertTrue(c1.equals(c2))

        c1 = Categorical(["a", "b", "c", "a"], categories=["a","c","b"])
        c2 = Categorical(Series(c1))
        self.assertTrue(c1.equals(c2))

        # Series
        c1 = Categorical(["a", "b", "c", "a"])
        c2 = Categorical(Series(["a", "b", "c", "a"]))
        self.assertTrue(c1.equals(c2))

        c1 = Categorical(["a", "b", "c", "a"], categories=["a","b","c","d"])
        c2 = Categorical(Series(["a", "b", "c", "a"]), categories=["a","b","c","d"])
        self.assertTrue(c1.equals(c2))

        # This should result in integer categories, not float!
        cat = pd.Categorical([1,2,3,np.nan], categories=[1,2,3])
        self.assertTrue(com.is_integer_dtype(cat.categories))

        # https://github.com/pydata/pandas/issues/3678
        cat = pd.Categorical([np.nan,1, 2, 3])
        self.assertTrue(com.is_integer_dtype(cat.categories))

        # this should result in floats
        cat = pd.Categorical([np.nan, 1, 2., 3 ])
        self.assertTrue(com.is_float_dtype(cat.categories))

        cat = pd.Categorical([np.nan, 1., 2., 3. ])
        self.assertTrue(com.is_float_dtype(cat.categories))

        # preserve int as far as possible by converting to object if NaN is in categories
        cat = pd.Categorical([np.nan, 1, 2, 3], categories=[np.nan, 1, 2, 3])
        self.assertTrue(com.is_object_dtype(cat.categories))
        # This doesn't work -> this would probably need some kind of "remember the original type"
        # feature to try to cast the array interface result to...
        #vals = np.asarray(cat[cat.notnull()])
        #self.assertTrue(com.is_integer_dtype(vals))
        cat = pd.Categorical([np.nan,"a", "b", "c"], categories=[np.nan,"a", "b", "c"])
        self.assertTrue(com.is_object_dtype(cat.categories))
        # but don't do it for floats
        cat = pd.Categorical([np.nan, 1., 2., 3.], categories=[np.nan, 1., 2., 3.])
        self.assertTrue(com.is_float_dtype(cat.categories))


        # corner cases
        cat = pd.Categorical([1])
        self.assertTrue(len(cat.categories) == 1)
        self.assertTrue(cat.categories[0] == 1)
        self.assertTrue(len(cat.codes) == 1)
        self.assertTrue(cat.codes[0] == 0)

        cat = pd.Categorical(["a"])
        self.assertTrue(len(cat.categories) == 1)
        self.assertTrue(cat.categories[0] == "a")
        self.assertTrue(len(cat.codes) == 1)
        self.assertTrue(cat.codes[0] == 0)

        # Scalars should be converted to lists
        cat = pd.Categorical(1)
        self.assertTrue(len(cat.categories) == 1)
        self.assertTrue(cat.categories[0] == 1)
        self.assertTrue(len(cat.codes) == 1)
        self.assertTrue(cat.codes[0] == 0)

        cat = pd.Categorical([1], categories=1)
        self.assertTrue(len(cat.categories) == 1)
        self.assertTrue(cat.categories[0] == 1)
        self.assertTrue(len(cat.codes) == 1)
        self.assertTrue(cat.codes[0] == 0)

        # Catch old style constructor useage: two arrays, codes + categories
        # We can only catch two cases:
        #  - when the first is an integer dtype and the second is not
        #  - when the resulting codes are all -1/NaN
        with tm.assert_produces_warning(RuntimeWarning):
            c_old = Categorical([0,1,2,0,1,2], categories=["a","b","c"])

        with tm.assert_produces_warning(RuntimeWarning):
            c_old = Categorical([0,1,2,0,1,2], categories=[3,4,5])

        # the next one are from the old docs, but unfortunately these don't trigger :-(
        with tm.assert_produces_warning(None):
            c_old2 = Categorical([0, 1, 2, 0, 1, 2], [1, 2, 3])
            cat = Categorical([1,2], categories=[1,2,3])

    def test_constructor_with_generator(self):
        # This was raising an Error in isnull(single_val).any() because isnull returned a scalar
        # for a generator
        from pandas.compat import range as xrange

        exp = Categorical([0,1,2])
        cat = Categorical((x for x in [0,1,2]))
        self.assertTrue(cat.equals(exp))
        cat = Categorical(xrange(3))
        self.assertTrue(cat.equals(exp))

        # This uses xrange internally
        from pandas.core.index import MultiIndex
        MultiIndex.from_product([range(5), ['a', 'b', 'c']])

        # check that categories accept generators and sequences
        cat = pd.Categorical([0,1,2], categories=(x for x in [0,1,2]))
        self.assertTrue(cat.equals(exp))
        cat = pd.Categorical([0,1,2], categories=xrange(3))
        self.assertTrue(cat.equals(exp))


    def test_from_codes(self):

        # too few categories
        def f():
            Categorical.from_codes([1,2], [1,2])
        self.assertRaises(ValueError, f)

        # no int codes
        def f():
            Categorical.from_codes(["a"], [1,2])
        self.assertRaises(ValueError, f)

        # no unique categories
        def f():
            Categorical.from_codes([0,1,2], ["a","a","b"])
        self.assertRaises(ValueError, f)

        # too negative
        def f():
            Categorical.from_codes([-2,1,2], ["a","b","c"])
        self.assertRaises(ValueError, f)


        exp = Categorical(["a","b","c"], ordered=False)
        res = Categorical.from_codes([0,1,2], ["a","b","c"])
        self.assertTrue(exp.equals(res))

        # Not available in earlier numpy versions
        if hasattr(np.random, "choice"):
            codes = np.random.choice([0,1], 5, p=[0.9,0.1])
            pd.Categorical.from_codes(codes, categories=["train", "test"])

    def test_comparisons(self):
        result = self.factor[self.factor == 'a']
        expected = self.factor[np.asarray(self.factor) == 'a']
        self.assertTrue(result.equals(expected))

        result = self.factor[self.factor != 'a']
        expected = self.factor[np.asarray(self.factor) != 'a']
        self.assertTrue(result.equals(expected))

        result = self.factor[self.factor < 'c']
        expected = self.factor[np.asarray(self.factor) < 'c']
        self.assertTrue(result.equals(expected))

        result = self.factor[self.factor > 'a']
        expected = self.factor[np.asarray(self.factor) > 'a']
        self.assertTrue(result.equals(expected))

        result = self.factor[self.factor >= 'b']
        expected = self.factor[np.asarray(self.factor) >= 'b']
        self.assertTrue(result.equals(expected))

        result = self.factor[self.factor <= 'b']
        expected = self.factor[np.asarray(self.factor) <= 'b']
        self.assertTrue(result.equals(expected))

        n = len(self.factor)

        other = self.factor[np.random.permutation(n)]
        result = self.factor == other
        expected = np.asarray(self.factor) == np.asarray(other)
        self.assert_numpy_array_equal(result, expected)

        result = self.factor == 'd'
        expected = np.repeat(False, len(self.factor))
        self.assert_numpy_array_equal(result, expected)

        # comparisons with categoricals
        cat_rev = pd.Categorical(["a","b","c"], categories=["c","b","a"])
        cat_rev_base = pd.Categorical(["b","b","b"], categories=["c","b","a"])
        cat = pd.Categorical(["a","b","c"])
        cat_base = pd.Categorical(["b","b","b"], categories=cat.categories)

        # comparisons need to take categories ordering into account
        res_rev = cat_rev > cat_rev_base
        exp_rev = np.array([True, False, False])
        self.assert_numpy_array_equal(res_rev, exp_rev)

        res_rev = cat_rev < cat_rev_base
        exp_rev = np.array([False, False, True])
        self.assert_numpy_array_equal(res_rev, exp_rev)

        res = cat > cat_base
        exp = np.array([False, False, True])
        self.assert_numpy_array_equal(res, exp)

        # Only categories with same categories can be compared
        def f():
            cat > cat_rev
        self.assertRaises(TypeError, f)

        cat_rev_base2 = pd.Categorical(["b","b","b"], categories=["c","b","a","d"])
        def f():
            cat_rev > cat_rev_base2
        self.assertRaises(TypeError, f)

        # Only categories with same ordering information can be compared
        cat_unorderd = cat.copy()
        cat_unorderd.ordered = False
        self.assertFalse((cat > cat).any())
        def f():
            cat > cat_unorderd
        self.assertRaises(TypeError, f)

        # comparison (in both directions) with Series will raise
        s = Series(["b","b","b"])
        self.assertRaises(TypeError, lambda: cat > s)
        self.assertRaises(TypeError, lambda: cat_rev > s)
        self.assertRaises(TypeError, lambda: s < cat)
        self.assertRaises(TypeError, lambda: s < cat_rev)

        # comparison with numpy.array will raise in both direction, but only on newer
        # numpy versions
        a = np.array(["b","b","b"])
        self.assertRaises(TypeError, lambda: cat > a)
        self.assertRaises(TypeError, lambda: cat_rev > a)

        # The following work via '__array_priority__ = 1000'
        # works only on numpy >= 1.7.1 and not on PY3.2
        if LooseVersion(np.__version__) > "1.7.1" and not compat.PY3_2:
            self.assertRaises(TypeError, lambda: a < cat)
            self.assertRaises(TypeError, lambda: a < cat_rev)

    def test_na_flags_int_categories(self):
        # #1457

        categories = lrange(10)
        labels = np.random.randint(0, 10, 20)
        labels[::5] = -1

        cat = Categorical(labels, categories, fastpath=True)
        repr(cat)

        self.assert_numpy_array_equal(com.isnull(cat), labels == -1)

    def test_categories_none(self):
        factor = Categorical(['a', 'b', 'b', 'a',
                              'a', 'c', 'c', 'c'])
        self.assertTrue(factor.equals(self.factor))

    def test_describe(self):
        # string type
        desc = self.factor.describe()
        expected = DataFrame.from_dict(dict(counts=[3, 2, 3],
                                            freqs=[3/8., 2/8., 3/8.],
                                            categories=['a', 'b', 'c'])
                                            ).set_index('categories')
        tm.assert_frame_equal(desc, expected)

        # check unused categories
        cat = self.factor.copy()
        cat.set_categories(["a","b","c","d"], inplace=True)
        desc = cat.describe()
        expected = DataFrame.from_dict(dict(counts=[3, 2, 3, np.nan],
                                            freqs=[3/8., 2/8., 3/8., np.nan],
                                            categories=['a', 'b', 'c', 'd'])
                                            ).set_index('categories')
        tm.assert_frame_equal(desc, expected)

        # check an integer one
        desc = Categorical([1,2,3,1,2,3,3,2,1,1,1]).describe()
        expected = DataFrame.from_dict(dict(counts=[5, 3, 3],
                                            freqs=[5/11., 3/11., 3/11.],
                                            categories=[1,2,3]
                                            )
                                            ).set_index('categories')
        tm.assert_frame_equal(desc, expected)

        # https://github.com/pydata/pandas/issues/3678
        # describe should work with NaN
        cat = pd.Categorical([np.nan,1, 2, 2])
        desc = cat.describe()
        expected = DataFrame.from_dict(dict(counts=[1, 2, 1],
                                            freqs=[1/4., 2/4., 1/4.],
                                            categories=[1,2,np.nan]
                                            )
                                            ).set_index('categories')
        tm.assert_frame_equal(desc, expected)

        # having NaN as category and as "not available" should also print two NaNs in describe!
        cat = pd.Categorical([np.nan,1, 2, 2])
        cat.set_categories([1,2,np.nan], rename=True, inplace=True)
        desc = cat.describe()
        expected = DataFrame.from_dict(dict(counts=[1, 2, np.nan, 1],
                                            freqs=[1/4., 2/4., np.nan, 1/4.],
                                            categories=[1,2,np.nan,np.nan]
                                            )
                                            ).set_index('categories')
        tm.assert_frame_equal(desc, expected)

        # empty categories show up as NA
        cat = Categorical(["a","b","b","b"], categories=['a','b','c'], ordered=True)
        result = cat.describe()

        expected = DataFrame([[1,0.25],[3,0.75],[np.nan,np.nan]],
                             columns=['counts','freqs'],
                             index=Index(['a','b','c'],name='categories'))
        tm.assert_frame_equal(result,expected)

        # NA as a category
        cat = pd.Categorical(["a","c","c",np.nan], categories=["b","a","c",np.nan] )
        result = cat.describe()

        expected = DataFrame([[np.nan, np.nan],[1,0.25],[2,0.5], [1,0.25]],
                             columns=['counts','freqs'],
                             index=Index(['b','a','c',np.nan],name='categories'))
        tm.assert_frame_equal(result,expected)


    def test_print(self):
        expected = ["[a, b, b, a, a, c, c, c]",
                    "Categories (3, object): [a < b < c]"]
        expected = "\n".join(expected)
        actual = repr(self.factor)
        self.assertEqual(actual, expected)

    def test_big_print(self):
        factor = Categorical([0,1,2,0,1,2]*100, ['a', 'b', 'c'], name='cat', fastpath=True)
        expected = ["[a, b, c, a, b, ..., b, c, a, b, c]",
                    "Name: cat, Length: 600",
                    "Categories (3, object): [a, b, c]"]
        expected = "\n".join(expected)

        actual = repr(factor)

        self.assertEqual(actual, expected)

    def test_empty_print(self):
        factor = Categorical([], ["a","b","c"], name="cat")
        expected = ("[], Name: cat, Categories (3, object): [a < b < c]")
        # hack because array_repr changed in numpy > 1.6.x
        actual = repr(factor)

        self.assertEqual(actual, expected)

        factor = Categorical([], ["a","b","c"])
        expected = ("[], Categories (3, object): [a < b < c]")
        actual = repr(factor)

        self.assertEqual(expected, actual)

        factor = Categorical([], [])
        expected = ("[], Categories (0, object): []")
        self.assertEqual(expected, repr(factor))

    def test_periodindex(self):
        idx1 = PeriodIndex(['2014-01', '2014-01', '2014-02', '2014-02',
                            '2014-03', '2014-03'], freq='M')

        cat1 = Categorical.from_array(idx1)
        str(cat1)
        exp_arr = np.array([0, 0, 1, 1, 2, 2],dtype='int64')
        exp_idx = PeriodIndex(['2014-01', '2014-02', '2014-03'], freq='M')
        self.assert_numpy_array_equal(cat1._codes, exp_arr)
        self.assertTrue(cat1.categories.equals(exp_idx))

        idx2 = PeriodIndex(['2014-03', '2014-03', '2014-02', '2014-01',
                            '2014-03', '2014-01'], freq='M')
        cat2 = Categorical.from_array(idx2)
        str(cat2)
        exp_arr = np.array([2, 2, 1, 0, 2, 0],dtype='int64')
        exp_idx2 = PeriodIndex(['2014-01', '2014-02', '2014-03'], freq='M')
        self.assert_numpy_array_equal(cat2._codes, exp_arr)
        self.assertTrue(cat2.categories.equals(exp_idx2))

        idx3 = PeriodIndex(['2013-12', '2013-11', '2013-10', '2013-09',
                            '2013-08', '2013-07', '2013-05'], freq='M')
        cat3 = Categorical.from_array(idx3)
        exp_arr = np.array([6, 5, 4, 3, 2, 1, 0],dtype='int64')
        exp_idx = PeriodIndex(['2013-05', '2013-07', '2013-08', '2013-09',
                               '2013-10', '2013-11', '2013-12'], freq='M')
        self.assert_numpy_array_equal(cat3._codes, exp_arr)
        self.assertTrue(cat3.categories.equals(exp_idx))

    def test_categories_assigments(self):
        s = pd.Categorical(["a","b","c","a"])
        exp = np.array([1,2,3,1])
        s.categories = [1,2,3]
        self.assert_numpy_array_equal(s.__array__(), exp)
        self.assert_numpy_array_equal(s.categories, np.array([1,2,3]))
        # lengthen
        def f():
            s.categories = [1,2,3,4]
        self.assertRaises(ValueError, f)
        # shorten
        def f():
            s.categories = [1,2]
        self.assertRaises(ValueError, f)

    def test_set_categories(self):
        cat = Categorical(["a","b","c","a"], ordered=True)
        exp_categories = np.array(["c","b","a"])
        exp_values = np.array(["a","b","c","a"])

        res = cat.set_categories(["c","b","a"], inplace=True)
        self.assert_numpy_array_equal(cat.categories, exp_categories)
        self.assert_numpy_array_equal(cat.__array__(), exp_values)
        self.assertIsNone(res)

        res = cat.set_categories(["a","b","c"])
        # cat must be the same as before
        self.assert_numpy_array_equal(cat.categories, exp_categories)
        self.assert_numpy_array_equal(cat.__array__(), exp_values)
        # only res is changed
        exp_categories_back = np.array(["a","b","c"])
        self.assert_numpy_array_equal(res.categories, exp_categories_back)
        self.assert_numpy_array_equal(res.__array__(), exp_values)

        # not all "old" included in "new" -> all not included ones are now np.nan
        cat = Categorical(["a","b","c","a"], ordered=True)
        res = cat.set_categories(["a"])
        self.assert_numpy_array_equal(res.codes, np.array([0,-1,-1,0]))

        # still not all "old" in "new"
        res = cat.set_categories(["a","b","d"])
        self.assert_numpy_array_equal(res.codes, np.array([0,1,-1,0]))
        self.assert_numpy_array_equal(res.categories, np.array(["a","b","d"]))

        # all "old" included in "new"
        cat = cat.set_categories(["a","b","c","d"])
        exp_categories = np.array(["a","b","c","d"])
        self.assert_numpy_array_equal(cat.categories, exp_categories)

        # internals...
        c = Categorical([1,2,3,4,1], categories=[1,2,3,4])
        self.assert_numpy_array_equal(c._codes, np.array([0,1,2,3,0]))
        self.assert_numpy_array_equal(c.categories , np.array([1,2,3,4] ))
        self.assert_numpy_array_equal(c.get_values(), np.array([1,2,3,4,1] ))
        c = c.set_categories([4,3,2,1]) # all "pointers" to '4' must be changed from 3 to 0,...
        self.assert_numpy_array_equal(c._codes, np.array([3,2,1,0,3])) # positions are changed
        self.assert_numpy_array_equal(c.categories, np.array([4,3,2,1])) # categories are now in new order
        self.assert_numpy_array_equal(c.get_values(), np.array([1,2,3,4,1])) # output is the same
        self.assertTrue(c.min(), 4)
        self.assertTrue(c.max(), 1)

    def test_rename_categories(self):
        cat = pd.Categorical(["a","b","c","a"])

        # inplace=False: the old one must not be changed
        res = cat.rename_categories([1,2,3])
        self.assert_numpy_array_equal(res.__array__(), np.array([1,2,3,1]))
        self.assert_numpy_array_equal(res.categories, np.array([1,2,3]))
        self.assert_numpy_array_equal(cat.__array__(), np.array(["a","b","c","a"]))
        self.assert_numpy_array_equal(cat.categories, np.array(["a","b","c"]))
        res = cat.rename_categories([1,2,3], inplace=True)

        # and now inplace
        self.assertIsNone(res)
        self.assert_numpy_array_equal(cat.__array__(), np.array([1,2,3,1]))
        self.assert_numpy_array_equal(cat.categories, np.array([1,2,3]))

        # lengthen
        def f():
            cat.rename_categories([1,2,3,4])
        self.assertRaises(ValueError, f)
        # shorten
        def f():
            cat.rename_categories([1,2])
        self.assertRaises(ValueError, f)

    def test_reorder_categories(self):
        cat = Categorical(["a","b","c","a"], ordered=True)
        old = cat.copy()
        new = Categorical(["a","b","c","a"], categories=["c","b","a"], ordered=True)

        # first inplace == False
        res = cat.reorder_categories(["c","b","a"])
        # cat must be the same as before
        self.assert_categorical_equal(cat, old)
        # only res is changed
        self.assert_categorical_equal(res, new)

        # inplace == True
        res = cat.reorder_categories(["c","b","a"], inplace=True)
        self.assertIsNone(res)
        self.assert_categorical_equal(cat, new)

        # not all "old" included in "new"
        cat = Categorical(["a","b","c","a"], ordered=True)
        def f():
            cat.reorder_categories(["a"])
        self.assertRaises(ValueError, f)

        # still not all "old" in "new"
        def f():
            cat.reorder_categories(["a","b","d"])
        self.assertRaises(ValueError, f)

        # all "old" included in "new", but too long
        def f():
            cat.reorder_categories(["a","b","c","d"])
        self.assertRaises(ValueError, f)

    def test_add_categories(self):
        cat = Categorical(["a","b","c","a"], ordered=True)
        old = cat.copy()
        new = Categorical(["a","b","c","a"], categories=["a","b","c","d"], ordered=True)

        # first inplace == False
        res = cat.add_categories("d")
        self.assert_categorical_equal(cat, old)
        self.assert_categorical_equal(res, new)

        res = cat.add_categories(["d"])
        self.assert_categorical_equal(cat, old)
        self.assert_categorical_equal(res, new)

        # inplace == True
        res = cat.add_categories("d", inplace=True)
        self.assert_categorical_equal(cat, new)
        self.assertIsNone(res)

        # new is in old categories
        def f():
            cat.add_categories(["d"])
        self.assertRaises(ValueError, f)

    def test_remove_categories(self):
        cat = Categorical(["a","b","c","a"], ordered=True)
        old = cat.copy()
        new = Categorical(["a","b",np.nan,"a"], categories=["a","b"], ordered=True)

        # first inplace == False
        res = cat.remove_categories("c")
        self.assert_categorical_equal(cat, old)
        self.assert_categorical_equal(res, new)

        res = cat.remove_categories(["c"])
        self.assert_categorical_equal(cat, old)
        self.assert_categorical_equal(res, new)

        # inplace == True
        res = cat.remove_categories("c", inplace=True)
        self.assert_categorical_equal(cat, new)
        self.assertIsNone(res)

        # removal is not in categories
        def f():
            cat.remove_categories(["c"])
        self.assertRaises(ValueError, f)

    def test_remove_unused_categories(self):
        c = Categorical(["a","b","c","d","a"], categories=["a","b","c","d","e"])
        exp_categories_all = np.array(["a","b","c","d","e"])
        exp_categories_dropped = np.array(["a","b","c","d"])

        self.assert_numpy_array_equal(c.categories, exp_categories_all)

        res = c.remove_unused_categories()
        self.assert_numpy_array_equal(res.categories, exp_categories_dropped)
        self.assert_numpy_array_equal(c.categories, exp_categories_all)

        res = c.remove_unused_categories(inplace=True)
        self.assert_numpy_array_equal(c.categories, exp_categories_dropped)
        self.assertIsNone(res)


    def test_nan_handling(self):

        # Nans are represented as -1 in codes
        c = Categorical(["a","b",np.nan,"a"])
        self.assert_numpy_array_equal(c.categories , np.array(["a","b"]))
        self.assert_numpy_array_equal(c._codes , np.array([0,1,-1,0]))
        c[1] = np.nan
        self.assert_numpy_array_equal(c.categories , np.array(["a","b"]))
        self.assert_numpy_array_equal(c._codes , np.array([0,-1,-1,0]))

        # If categories have nan included, the code should point to that instead
        c = Categorical(["a","b",np.nan,"a"], categories=["a","b",np.nan])
        self.assert_numpy_array_equal(c.categories , np.array(["a","b",np.nan],dtype=np.object_))
        self.assert_numpy_array_equal(c._codes , np.array([0,1,2,0]))
        c[1] = np.nan
        self.assert_numpy_array_equal(c.categories , np.array(["a","b",np.nan],dtype=np.object_))
        self.assert_numpy_array_equal(c._codes , np.array([0,2,2,0]))

        # Changing categories should also make the replaced category np.nan
        c = Categorical(["a","b","c","a"])
        c.categories = ["a","b",np.nan]
        self.assert_numpy_array_equal(c.categories , np.array(["a","b",np.nan],dtype=np.object_))
        self.assert_numpy_array_equal(c._codes , np.array([0,1,2,0]))

        # Adding nan to categories should make assigned nan point to the category!
        c = Categorical(["a","b",np.nan,"a"])
        self.assert_numpy_array_equal(c.categories , np.array(["a","b"]))
        self.assert_numpy_array_equal(c._codes , np.array([0,1,-1,0]))
        c.set_categories(["a","b",np.nan], rename=True, inplace=True)
        self.assert_numpy_array_equal(c.categories , np.array(["a","b",np.nan],dtype=np.object_))
        self.assert_numpy_array_equal(c._codes , np.array([0,1,-1,0]))
        c[1] = np.nan
        self.assert_numpy_array_equal(c.categories , np.array(["a","b",np.nan],dtype=np.object_))
        self.assert_numpy_array_equal(c._codes , np.array([0,2,-1,0]))


    def test_isnull(self):
        exp = np.array([False, False, True])
        c = Categorical(["a","b",np.nan])
        res = c.isnull()
        self.assert_numpy_array_equal(res, exp)

        c = Categorical(["a","b",np.nan], categories=["a","b",np.nan])
        res = c.isnull()
        self.assert_numpy_array_equal(res, exp)

        # test both nan in categories and as -1
        exp = np.array([True, False, True])
        c = Categorical(["a","b",np.nan])
        c.set_categories(["a","b",np.nan], rename=True, inplace=True)
        c[0] = np.nan
        res = c.isnull()
        self.assert_numpy_array_equal(res, exp)

    def test_codes_immutable(self):

        # Codes should be read only
        c = Categorical(["a","b","c","a", np.nan])
        exp = np.array([0,1,2,0,-1],dtype='int8')
        self.assert_numpy_array_equal(c.codes, exp)

        # Assignments to codes should raise
        def f():
            c.codes = np.array([0,1,2,0,1],dtype='int8')
        self.assertRaises(ValueError, f)

        # changes in the codes array should raise
        # np 1.6.1 raises RuntimeError rather than ValueError
        codes= c.codes
        def f():
            codes[4] = 1
        self.assertRaises(ValueError, f)

        # But even after getting the codes, the original array should still be writeable!
        c[4] = "a"
        exp = np.array([0,1,2,0,0],dtype='int8')
        self.assert_numpy_array_equal(c.codes, exp)
        c._codes[4] = 2
        exp = np.array([0,1,2,0, 2],dtype='int8')
        self.assert_numpy_array_equal(c.codes, exp)


    def test_min_max(self):

        # unordered cats have no min/max
        cat = Categorical(["a","b","c","d"], ordered=False)
        self.assertRaises(TypeError, lambda : cat.min())
        self.assertRaises(TypeError, lambda : cat.max())
        cat = Categorical(["a","b","c","d"], ordered=True)
        _min = cat.min()
        _max = cat.max()
        self.assertEqual(_min, "a")
        self.assertEqual(_max, "d")
        cat = Categorical(["a","b","c","d"], categories=['d','c','b','a'], ordered=True)
        _min = cat.min()
        _max = cat.max()
        self.assertEqual(_min, "d")
        self.assertEqual(_max, "a")
        cat = Categorical([np.nan,"b","c",np.nan], categories=['d','c','b','a'], ordered=True)
        _min = cat.min()
        _max = cat.max()
        self.assertTrue(np.isnan(_min))
        self.assertEqual(_max, "b")

        _min = cat.min(numeric_only=True)
        self.assertEqual(_min, "c")
        _max = cat.max(numeric_only=True)
        self.assertEqual(_max, "b")

        cat = Categorical([np.nan,1,2,np.nan], categories=[5,4,3,2,1], ordered=True)
        _min = cat.min()
        _max = cat.max()
        self.assertTrue(np.isnan(_min))
        self.assertEqual(_max, 1)

        _min = cat.min(numeric_only=True)
        self.assertEqual(_min, 2)
        _max = cat.max(numeric_only=True)
        self.assertEqual(_max, 1)

    def test_unique(self):
        cat = Categorical(["a","b"])
        exp = np.asarray(["a","b"])
        res = cat.unique()
        self.assert_numpy_array_equal(res, exp)
        cat = Categorical(["a","b","a","a"], categories=["a","b","c"])
        res = cat.unique()
        self.assert_numpy_array_equal(res, exp)
        cat = Categorical(["a","b","a", np.nan], categories=["a","b","c"])
        res = cat.unique()
        exp = np.asarray(["a","b", np.nan], dtype=object)
        self.assert_numpy_array_equal(res, exp)

    def test_mode(self):
        s = Categorical([1,1,2,4,5,5,5], categories=[5,4,3,2,1], ordered=True)
        res = s.mode()
        exp = Categorical([5], categories=[5,4,3,2,1], ordered=True)
        self.assertTrue(res.equals(exp))
        s = Categorical([1,1,1,4,5,5,5], categories=[5,4,3,2,1], ordered=True)
        res = s.mode()
        exp = Categorical([5,1], categories=[5,4,3,2,1], ordered=True)
        self.assertTrue(res.equals(exp))
        s = Categorical([1,2,3,4,5], categories=[5,4,3,2,1], ordered=True)
        res = s.mode()
        exp = Categorical([], categories=[5,4,3,2,1], ordered=True)
        self.assertTrue(res.equals(exp))
        # NaN should not become the mode!
        s = Categorical([np.nan,np.nan,np.nan,4,5], categories=[5,4,3,2,1], ordered=True)
        res = s.mode()
        exp = Categorical([], categories=[5,4,3,2,1], ordered=True)
        self.assertTrue(res.equals(exp))
        s = Categorical([np.nan,np.nan,np.nan,4,5,4], categories=[5,4,3,2,1], ordered=True)
        res = s.mode()
        exp = Categorical([4], categories=[5,4,3,2,1], ordered=True)
        self.assertTrue(res.equals(exp))
        s = Categorical([np.nan,np.nan,4,5,4], categories=[5,4,3,2,1], ordered=True)
        res = s.mode()
        exp = Categorical([4], categories=[5,4,3,2,1], ordered=True)
        self.assertTrue(res.equals(exp))


    def test_sort(self):

        # unordered cats are not sortable
        cat = Categorical(["a","b","b","a"], ordered=False)
        self.assertRaises(TypeError, lambda : cat.sort())
        cat = Categorical(["a","c","b","d"], ordered=True)

        # order
        res = cat.order()
        exp = np.array(["a","b","c","d"],dtype=object)
        self.assert_numpy_array_equal(res.__array__(), exp)

        cat = Categorical(["a","c","b","d"], categories=["a","b","c","d"], ordered=True)
        res = cat.order()
        exp = np.array(["a","b","c","d"],dtype=object)
        self.assert_numpy_array_equal(res.__array__(), exp)

        res = cat.order(ascending=False)
        exp = np.array(["d","c","b","a"],dtype=object)
        self.assert_numpy_array_equal(res.__array__(), exp)

        # sort (inplace order)
        cat1 = cat.copy()
        cat1.sort()
        exp = np.array(["a","b","c","d"],dtype=object)
        self.assert_numpy_array_equal(cat1.__array__(), exp)

    def test_slicing_directly(self):
        cat = Categorical(["a","b","c","d","a","b","c"])
        sliced = cat[3]
        tm.assert_equal(sliced, "d")
        sliced = cat[3:5]
        expected = Categorical(["d","a"], categories=['a', 'b', 'c', 'd'])
        self.assert_numpy_array_equal(sliced._codes, expected._codes)
        tm.assert_index_equal(sliced.categories, expected.categories)

    def test_set_item_nan(self):
        cat = pd.Categorical([1,2,3])
        exp = pd.Categorical([1,np.nan,3], categories=[1,2,3])
        cat[1] = np.nan
        self.assertTrue(cat.equals(exp))

        # if nan in categories, the proper code should be set!
        cat = pd.Categorical([1,2,3, np.nan], categories=[1,2,3])
        cat.set_categories([1,2,3, np.nan], rename=True, inplace=True)
        cat[1] = np.nan
        exp = np.array([0,3,2,-1])
        self.assert_numpy_array_equal(cat.codes, exp)

        cat = pd.Categorical([1,2,3, np.nan], categories=[1,2,3])
        cat.set_categories([1,2,3, np.nan], rename=True, inplace=True)
        cat[1:3] = np.nan
        exp = np.array([0,3,3,-1])
        self.assert_numpy_array_equal(cat.codes, exp)

        cat = pd.Categorical([1,2,3, np.nan], categories=[1,2,3])
        cat.set_categories([1,2,3, np.nan], rename=True, inplace=True)
        cat[1:3] = [np.nan, 1]
        exp = np.array([0,3,0,-1])
        self.assert_numpy_array_equal(cat.codes, exp)

        cat = pd.Categorical([1,2,3, np.nan], categories=[1,2,3])
        cat.set_categories([1,2,3, np.nan], rename=True, inplace=True)
        cat[1:3] = [np.nan, np.nan]
        exp = np.array([0,3,3,-1])
        self.assert_numpy_array_equal(cat.codes, exp)

        cat = pd.Categorical([1,2, np.nan, 3], categories=[1,2,3])
        cat.set_categories([1,2,3, np.nan], rename=True, inplace=True)
        cat[pd.isnull(cat)] = np.nan
        exp = np.array([0,1,3,2])
        self.assert_numpy_array_equal(cat.codes, exp)

    def test_nbytes(self):
        cat = pd.Categorical([1,2,3])
        exp = cat._codes.nbytes + cat._categories.values.nbytes
        self.assertEqual(cat.nbytes, exp)

    def test_searchsorted(self):
        # https://github.com/pydata/pandas/issues/8420
        s1 = pd.Series(['apple', 'bread', 'bread', 'cheese', 'milk' ])
        s2 = pd.Series(['apple', 'bread', 'bread', 'cheese', 'milk', 'donuts' ])
        c1 = pd.Categorical(s1)
        c2 = pd.Categorical(s2)

        # Single item array
        res = c1.searchsorted(['bread'])
        chk = s1.searchsorted(['bread'])
        exp = np.array([1])
        self.assert_numpy_array_equal(res, exp)
        self.assert_numpy_array_equal(res, chk)

        # Scalar version of single item array
        # Categorical return np.array like pd.Series, but different from np.array.searchsorted()
        res = c1.searchsorted('bread')
        chk = s1.searchsorted('bread')
        exp = np.array([1])
        self.assert_numpy_array_equal(res, exp)
        self.assert_numpy_array_equal(res, chk)
       
        # Searching for a value that is not present in the Categorical 
        res = c1.searchsorted(['bread', 'eggs'])
        chk = s1.searchsorted(['bread', 'eggs'])
        exp = np.array([1, 4])
        self.assert_numpy_array_equal(res, exp)
        self.assert_numpy_array_equal(res, chk)

        # Searching for a value that is not present, to the right
        res = c1.searchsorted(['bread', 'eggs'], side='right')
        chk = s1.searchsorted(['bread', 'eggs'], side='right')
        exp = np.array([3, 4])	    # eggs before milk
        self.assert_numpy_array_equal(res, exp)
        self.assert_numpy_array_equal(res, chk)

        # As above, but with a sorter array to reorder an unsorted array
        res = c2.searchsorted(['bread', 'eggs'], side='right', sorter=[0, 1, 2, 3, 5, 4])
        chk = s2.searchsorted(['bread', 'eggs'], side='right', sorter=[0, 1, 2, 3, 5, 4])
        exp = np.array([3, 5])       # eggs after donuts, after switching milk and donuts 
        self.assert_numpy_array_equal(res, exp)
        self.assert_numpy_array_equal(res, chk)

    def test_deprecated_labels(self):
        # TODO: labels is deprecated and should be removed in 0.18 or 2017, whatever is earlier
        cat = pd.Categorical([1,2,3, np.nan], categories=[1,2,3])
        exp = cat.codes
        with tm.assert_produces_warning(FutureWarning):
            res = cat.labels
        self.assert_numpy_array_equal(res, exp)
        self.assertFalse(LooseVersion(pd.__version__) >= '0.18')

    def test_deprecated_levels(self):
        # TODO: levels is deprecated and should be removed in 0.18 or 2017, whatever is earlier
        cat = pd.Categorical([1,2,3, np.nan], categories=[1,2,3])
        exp = cat.categories
        with tm.assert_produces_warning(FutureWarning):
            res = cat.levels
        self.assert_numpy_array_equal(res, exp)
        with tm.assert_produces_warning(FutureWarning):
            res = pd.Categorical([1,2,3, np.nan], levels=[1,2,3])
        self.assert_numpy_array_equal(res.categories, exp)

        self.assertFalse(LooseVersion(pd.__version__) >= '0.18')

    def test_datetime_categorical_comparison(self):
        dt_cat = pd.Categorical(pd.date_range('2014-01-01', periods=3))
        self.assert_numpy_array_equal(dt_cat > dt_cat[0], [False, True, True])
        self.assert_numpy_array_equal(dt_cat[0] < dt_cat, [False, True, True])

    def test_reflected_comparison_with_scalars(self):
        # GH8658
        cat = pd.Categorical([1, 2, 3])
        self.assert_numpy_array_equal(cat > cat[0], [False, True, True])
        self.assert_numpy_array_equal(cat[0] < cat, [False, True, True])


class TestCategoricalAsBlock(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.factor = Categorical.from_array(['a', 'b', 'b', 'a',
                                              'a', 'c', 'c', 'c'])

        df = DataFrame({'value': np.random.randint(0, 10000, 100)})
        labels = [ "{0} - {1}".format(i, i + 499) for i in range(0, 10000, 500) ]

        df = df.sort(columns=['value'], ascending=True)
        df['value_group'] = pd.cut(df.value, range(0, 10500, 500), right=False, labels=labels)
        self.cat = df

    def test_dtypes(self):

        dtype = com.CategoricalDtype()
        hash(dtype)
        self.assertTrue(com.is_categorical_dtype(dtype))

        s = Series(self.factor,name='A')

        # dtypes
        self.assertTrue(com.is_categorical_dtype(s.dtype))
        self.assertTrue(com.is_categorical_dtype(s))
        self.assertFalse(com.is_categorical_dtype(np.dtype('float64')))

        # np.dtype doesn't know about our new dtype
        def f():
            np.dtype(dtype)
        self.assertRaises(TypeError, f)

        self.assertFalse(dtype == np.str_)
        self.assertFalse(np.str_ == dtype)

        # GH8143
        index = ['cat','obj','num']
        cat = pd.Categorical(['a', 'b', 'c'])
        obj = pd.Series(['a', 'b', 'c'])
        num = pd.Series([1, 2, 3])
        df = pd.concat([pd.Series(cat), obj, num], axis=1, keys=index)

        result = df.dtypes == 'object'
        expected = Series([False,True,False],index=index)
        tm.assert_series_equal(result, expected)

        result = df.dtypes == 'int64'
        expected = Series([False,False,True],index=index)
        tm.assert_series_equal(result, expected)

        result = df.dtypes == 'category'
        expected = Series([True,False,False],index=index)
        tm.assert_series_equal(result, expected)

    def test_codes_dtypes(self):

        # GH 8453
        result = Categorical(['foo','bar','baz'])
        self.assertTrue(result.codes.dtype == 'int8')

        result = Categorical(['foo%05d' % i for i in range(400) ])
        self.assertTrue(result.codes.dtype == 'int16')

        result = Categorical(['foo%05d' % i for i in range(40000) ])
        self.assertTrue(result.codes.dtype == 'int32')

        # adding cats
        result = Categorical(['foo','bar','baz'])
        self.assertTrue(result.codes.dtype == 'int8')
        result = result.add_categories(['foo%05d' % i for i in range(400) ])
        self.assertTrue(result.codes.dtype == 'int16')

        # removing cats
        result = result.remove_categories(['foo%05d' % i for i in range(300) ])
        self.assertTrue(result.codes.dtype == 'int8')

    def test_basic(self):

        # test basic creation / coercion of categoricals
        s = Series(self.factor,name='A')
        self.assertEqual(s.dtype,'category')
        self.assertEqual(len(s),len(self.factor))
        str(s.values)
        str(s)

        # in a frame
        df = DataFrame({'A' : self.factor })
        result = df['A']
        tm.assert_series_equal(result,s)
        result = df.iloc[:,0]
        tm.assert_series_equal(result,s)
        self.assertEqual(len(df),len(self.factor))
        str(df.values)
        str(df)

        df = DataFrame({'A' : s })
        result = df['A']
        tm.assert_series_equal(result,s)
        self.assertEqual(len(df),len(self.factor))
        str(df.values)
        str(df)

        # multiples
        df = DataFrame({'A' : s, 'B' : s, 'C' : 1})
        result1 = df['A']
        result2 = df['B']
        tm.assert_series_equal(result1,s)
        tm.assert_series_equal(result2,s)
        self.assertEqual(len(df),len(self.factor))
        str(df.values)
        str(df)

        # GH8623
        x = pd.DataFrame([[1,'John P. Doe'],[2,'Jane Dove'],[1,'John P. Doe']],
                         columns=['person_id','person_name'])
        x['person_name'] = pd.Categorical(x.person_name) # doing this breaks transform

        expected = x.iloc[0].person_name
        result = x.person_name.iloc[0]
        self.assertEqual(result,expected)

        result = x.person_name[0]
        self.assertEqual(result,expected)

        result = x.person_name.loc[0]
        self.assertEqual(result,expected)

    def test_creation_astype(self):
        l = ["a","b","c","a"]
        s = pd.Series(l)
        exp = pd.Series(Categorical(l))
        res = s.astype('category')
        tm.assert_series_equal(res, exp)

        l = [1,2,3,1]
        s = pd.Series(l)
        exp = pd.Series(Categorical(l))
        res = s.astype('category')
        tm.assert_series_equal(res, exp)

        df = pd.DataFrame({"cats":[1,2,3,4,5,6], "vals":[1,2,3,4,5,6]})
        cats = Categorical([1,2,3,4,5,6])
        exp_df = pd.DataFrame({"cats":cats, "vals":[1,2,3,4,5,6]})
        df["cats"] =  df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

        df = pd.DataFrame({"cats":['a', 'b', 'b', 'a', 'a', 'd'], "vals":[1,2,3,4,5,6]})
        cats = Categorical(['a', 'b', 'b', 'a', 'a', 'd'])
        exp_df = pd.DataFrame({"cats":cats, "vals":[1,2,3,4,5,6]})
        df["cats"] =  df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

    def test_construction_series(self):

        l = [1,2,3,1]
        exp = Series(l).astype('category')
        res = Series(l,dtype='category')
        tm.assert_series_equal(res, exp)

        l = ["a","b","c","a"]
        exp = Series(l).astype('category')
        res = Series(l,dtype='category')
        tm.assert_series_equal(res, exp)

        # insert into frame with different index
        # GH 8076
        index = pd.date_range('20000101', periods=3)
        expected = Series(Categorical(values=[np.nan,np.nan,np.nan],categories=['a', 'b', 'c']))
        expected.index = index

        expected = DataFrame({'x': expected})
        df = DataFrame({'x': Series(['a', 'b', 'c'],dtype='category')}, index=index)
        tm.assert_frame_equal(df, expected)

    def test_construction_frame(self):

        # GH8626

        # dict creation
        df = DataFrame({ 'A' : list('abc') },dtype='category')
        expected = Series(list('abc'),dtype='category')
        tm.assert_series_equal(df['A'],expected)

        # to_frame
        s = Series(list('abc'),dtype='category')
        result = s.to_frame()
        expected = Series(list('abc'),dtype='category')
        tm.assert_series_equal(result[0],expected)
        result = s.to_frame(name='foo')
        expected = Series(list('abc'),dtype='category')
        tm.assert_series_equal(result['foo'],expected)

        # list-like creation
        df = DataFrame(list('abc'),dtype='category')
        expected = Series(list('abc'),dtype='category')
        tm.assert_series_equal(df[0],expected)

        # ndim != 1
        df = DataFrame([pd.Categorical(list('abc'))])
        expected = DataFrame({ 0 : Series(list('abc'),dtype='category')})
        tm.assert_frame_equal(df,expected)

        df = DataFrame([pd.Categorical(list('abc')),pd.Categorical(list('abd'))])
        expected = DataFrame({ 0 : Series(list('abc'),dtype='category'),
                               1 : Series(list('abd'),dtype='category')},columns=[0,1])
        tm.assert_frame_equal(df,expected)

        # mixed
        df = DataFrame([pd.Categorical(list('abc')),list('def')])
        expected = DataFrame({ 0 : Series(list('abc'),dtype='category'),
                               1 : list('def')},columns=[0,1])
        tm.assert_frame_equal(df,expected)

        # invalid (shape)
        self.assertRaises(ValueError, lambda : DataFrame([pd.Categorical(list('abc')),pd.Categorical(list('abdefg'))]))

        # ndim > 1
        self.assertRaises(NotImplementedError, lambda : pd.Categorical(np.array([list('abcd')])))

    def test_reshaping(self):

        p = tm.makePanel()
        p['str'] = 'foo'
        df = p.to_frame()
        df['category'] = df['str'].astype('category')
        result = df['category'].unstack()

        c = Categorical(['foo']*len(p.major_axis))
        expected = DataFrame({'A' : c.copy(),
                              'B' : c.copy(),
                              'C' : c.copy(),
                              'D' : c.copy()},
                             columns=Index(list('ABCD'),name='minor'),
                             index=p.major_axis.set_names('major'))
        tm.assert_frame_equal(result, expected)

    def test_reindex(self):

        index = pd.date_range('20000101', periods=3)

        # reindexing to an invalid Categorical
        s = Series(['a', 'b', 'c'],dtype='category')
        result = s.reindex(index)
        expected = Series(Categorical(values=[np.nan,np.nan,np.nan],categories=['a', 'b', 'c']))
        expected.index = index
        tm.assert_series_equal(result, expected)

        # partial reindexing
        expected = Series(Categorical(values=['b','c'],categories=['a', 'b', 'c']))
        expected.index = [1,2]
        result = s.reindex([1,2])
        tm.assert_series_equal(result, expected)

        expected = Series(Categorical(values=['c',np.nan],categories=['a', 'b', 'c']))
        expected.index = [2,3]
        result = s.reindex([2,3])
        tm.assert_series_equal(result, expected)



    def test_sideeffects_free(self):

        # Passing a categorical to a Series and then changing values in either the series or the
        # categorical should not change the values in the other one, IF you specify copy!
        cat = Categorical(["a","b","c","a"])
        s =  pd.Series(cat, copy=True)
        self.assertFalse(s.cat is cat)
        s.cat.categories = [1,2,3]
        exp_s = np.array([1,2,3,1])
        exp_cat = np.array(["a","b","c","a"])
        self.assert_numpy_array_equal(s.__array__(), exp_s)
        self.assert_numpy_array_equal(cat.__array__(), exp_cat)

        # setting
        s[0] = 2
        exp_s2 = np.array([2,2,3,1])
        self.assert_numpy_array_equal(s.__array__(), exp_s2)
        self.assert_numpy_array_equal(cat.__array__(), exp_cat)

        # however, copy is False by default
        # so this WILL change values
        cat = Categorical(["a","b","c","a"])
        s =  pd.Series(cat)
        self.assertTrue(s.values is cat)
        s.cat.categories = [1,2,3]
        exp_s = np.array([1,2,3,1])
        self.assert_numpy_array_equal(s.__array__(), exp_s)
        self.assert_numpy_array_equal(cat.__array__(), exp_s)

        s[0] = 2
        exp_s2 = np.array([2,2,3,1])
        self.assert_numpy_array_equal(s.__array__(), exp_s2)
        self.assert_numpy_array_equal(cat.__array__(), exp_s2)

    def test_nan_handling(self):

        # Nans are represented as -1 in labels
        s = Series(Categorical(["a","b",np.nan,"a"]))
        self.assert_numpy_array_equal(s.cat.categories, np.array(["a","b"]))
        self.assert_numpy_array_equal(s.values.codes, np.array([0,1,-1,0]))

        # If categories have nan included, the label should point to that instead
        s2 = Series(Categorical(["a","b",np.nan,"a"], categories=["a","b",np.nan]))
        self.assert_numpy_array_equal(s2.cat.categories,
                                      np.array(["a","b",np.nan], dtype=np.object_))
        self.assert_numpy_array_equal(s2.values.codes, np.array([0,1,2,0]))

        # Changing categories should also make the replaced category np.nan
        s3 = Series(Categorical(["a","b","c","a"]))
        s3.cat.categories = ["a","b",np.nan]
        self.assert_numpy_array_equal(s3.cat.categories,
                                      np.array(["a","b",np.nan], dtype=np.object_))
        self.assert_numpy_array_equal(s3.values.codes, np.array([0,1,2,0]))

    def test_cat_accessor(self):
        s = Series(Categorical(["a","b",np.nan,"a"]))
        self.assert_numpy_array_equal(s.cat.categories, np.array(["a","b"]))
        self.assertEqual(s.cat.ordered, True)
        exp = Categorical(["a","b",np.nan,"a"], categories=["b","a"])
        s.cat.set_categories(["b", "a"], inplace=True)
        self.assertTrue(s.values.equals(exp))
        res = s.cat.set_categories(["b", "a"])
        self.assertTrue(res.values.equals(exp))
        exp = Categorical(["a","b",np.nan,"a"], categories=["b","a"])
        s[:] = "a"
        s = s.cat.remove_unused_categories()
        self.assert_numpy_array_equal(s.cat.categories, np.array(["a"]))

    def test_sequence_like(self):

        # GH 7839
        # make sure can iterate
        df = DataFrame({"id":[1,2,3,4,5,6], "raw_grade":['a', 'b', 'b', 'a', 'a', 'e']})
        df['grade'] = Categorical(df['raw_grade'])

        # basic sequencing testing
        result = list(df.grade.values)
        expected = np.array(df.grade.values).tolist()
        tm.assert_almost_equal(result,expected)

        # iteration
        for t in df.itertuples(index=False):
            str(t)

        for row, s in df.iterrows():
            str(s)

        for c, col in df.iteritems():
            str(s)

    def test_series_delegations(self):

        # invalid accessor
        self.assertRaises(TypeError, lambda : Series([1,2,3]).cat)
        tm.assertRaisesRegexp(TypeError,
                              r"Can only use .cat accessor with a 'category' dtype",
                              lambda : Series([1,2,3]).cat)
        self.assertRaises(TypeError, lambda : Series(['a','b','c']).cat)
        self.assertRaises(TypeError, lambda : Series(np.arange(5.)).cat)
        self.assertRaises(TypeError, lambda : Series([Timestamp('20130101')]).cat)

        # Series should delegate calls to '.categories', '.codes', '.ordered' and the
        # methods '.set_categories()' 'drop_unused_categories()' to the categorical
        s = Series(Categorical(["a","b","c","a"], ordered=True))
        exp_categories = np.array(["a","b","c"])
        self.assert_numpy_array_equal(s.cat.categories, exp_categories)
        s.cat.categories = [1,2,3]
        exp_categories = np.array([1,2,3])
        self.assert_numpy_array_equal(s.cat.categories, exp_categories)

        exp_codes = Series([0,1,2,0],dtype='int8')
        tm.assert_series_equal(s.cat.codes, exp_codes)

        self.assertEqual(s.cat.ordered, True)
        s.cat.ordered = False
        self.assertEqual(s.cat.ordered, False)

        # reorder
        s = Series(Categorical(["a","b","c","a"], ordered=True))
        exp_categories = np.array(["c","b","a"])
        exp_values = np.array(["a","b","c","a"])
        s = s.cat.set_categories(["c","b","a"])
        self.assert_numpy_array_equal(s.cat.categories, exp_categories)
        self.assert_numpy_array_equal(s.values.__array__(), exp_values)
        self.assert_numpy_array_equal(s.__array__(), exp_values)

        # remove unused categories
        s = Series(Categorical(["a","b","b","a"], categories=["a","b","c"]))
        exp_categories = np.array(["a","b"])
        exp_values = np.array(["a","b","b","a"])
        s = s.cat.remove_unused_categories()
        self.assert_numpy_array_equal(s.cat.categories, exp_categories)
        self.assert_numpy_array_equal(s.values.__array__(), exp_values)
        self.assert_numpy_array_equal(s.__array__(), exp_values)

        # This method is likely to be confused, so test that it raises an error on wrong inputs:
        def f():
            s.set_categories([4,3,2,1])
        self.assertRaises(Exception, f)
        # right: s.cat.set_categories([4,3,2,1])

    def test_series_functions_no_warnings(self):
        df = pd.DataFrame({'value': np.random.randint(0, 100, 20)})
        labels = [ "{0} - {1}".format(i, i + 9) for i in range(0, 100, 10)]
        with tm.assert_produces_warning(False):
            df['group'] = pd.cut(df.value, range(0, 105, 10), right=False, labels=labels)

    def test_assignment_to_dataframe(self):
        # assignment
        df = DataFrame({'value': np.array(np.random.randint(0, 10000, 100),dtype='int32')})
        labels = [ "{0} - {1}".format(i, i + 499) for i in range(0, 10000, 500) ]

        df = df.sort(columns=['value'], ascending=True)
        s = pd.cut(df.value, range(0, 10500, 500), right=False, labels=labels)
        d = s.values
        df['D'] = d
        str(df)

        result = df.dtypes
        expected = Series([np.dtype('int32'), com.CategoricalDtype()],index=['value','D'])
        tm.assert_series_equal(result,expected)

        df['E'] = s
        str(df)

        result = df.dtypes
        expected = Series([np.dtype('int32'), com.CategoricalDtype(), com.CategoricalDtype()],
                          index=['value','D','E'])
        tm.assert_series_equal(result,expected)

        result1 = df['D']
        result2 = df['E']
        self.assertTrue(result1._data._block.values.equals(d))

        # sorting
        s.name = 'E'
        self.assertTrue(result2.sort_index().equals(s.sort_index()))

        cat = pd.Categorical([1,2,3,10], categories=[1,2,3,4,10])
        df = pd.DataFrame(pd.Series(cat))

    def test_describe(self):

        # Categoricals should not show up together with numerical columns
        result = self.cat.describe()
        self.assertEquals(len(result.columns),1)


        # In a frame, describe() for the cat should be the same as for string arrays (count, unique,
        # top, freq)

        cat = Categorical(["a","b","b","b"], categories=['a','b','c'], ordered=True)
        s = Series(cat)
        result = s.describe()
        expected = Series([4,2,"b",3],index=['count','unique','top', 'freq'])
        tm.assert_series_equal(result,expected)

        cat = pd.Series(pd.Categorical(["a","b","c","c"]))
        df3 = pd.DataFrame({"cat":cat, "s":["a","b","c","c"]})
        res = df3.describe()
        self.assert_numpy_array_equal(res["cat"].values, res["s"].values)

    def test_repr(self):
        a = pd.Series(pd.Categorical([1,2,3,4], name="a"))
        exp = u("0    1\n1    2\n2    3\n3    4\n" +
              "Name: a, dtype: category\nCategories (4, int64): [1 < 2 < 3 < 4]")

        self.assertEqual(exp, a.__unicode__())

        a = pd.Series(pd.Categorical(["a","b"] *25, name="a"))
        exp = u("".join(["%s    a\n%s    b\n"%(i,i+1) for i in range(0,10,2)]) + "...\n" +
                "".join(["%s    a\n%s    b\n"%(i,i+1) for i in range(40,50,2)]) +
                "Name: a, Length: 50, dtype: category\n" +
                "Categories (2, object): [a < b]")
        self.assertEqual(exp,a._tidy_repr())

        levs = list("abcdefghijklmnopqrstuvwxyz")
        a = pd.Series(pd.Categorical(["a","b"], name="a", categories=levs))
        exp = u("0    a\n1    b\n" +
                "Name: a, dtype: category\n"
                "Categories (26, object): [a < b < c < d ... w < x < y < z]")
        self.assertEqual(exp,a.__unicode__())


    def test_info(self):

        # make sure it works
        n = 2500
        df = DataFrame({ 'int64' : np.random.randint(100,size=n) })
        df['category'] = Series(np.array(list('abcdefghij')).take(np.random.randint(0,10,size=n))).astype('category')
        df.isnull()
        df.info()

        df2 = df[df['category']=='d']
        df2.info()

    def test_groupby_sort(self):

        # http://stackoverflow.com/questions/23814368/sorting-pandas-categorical-labels-after-groupby
        # This should result in a properly sorted Series so that the plot
        # has a sorted x axis
        #self.cat.groupby(['value_group'])['value_group'].count().plot(kind='bar')

        res = self.cat.groupby(['value_group'])['value_group'].count()
        exp = res[sorted(res.index, key=lambda x: float(x.split()[0]))]
        tm.assert_series_equal(res, exp)

    def test_min_max(self):
        # unordered cats have no min/max
        cat = Series(Categorical(["a","b","c","d"], ordered=False))
        self.assertRaises(TypeError, lambda : cat.min())
        self.assertRaises(TypeError, lambda : cat.max())

        cat = Series(Categorical(["a","b","c","d"], ordered=True))
        _min = cat.min()
        _max = cat.max()
        self.assertEqual(_min, "a")
        self.assertEqual(_max, "d")

        cat = Series(Categorical(["a","b","c","d"], categories=['d','c','b','a'], ordered=True))
        _min = cat.min()
        _max = cat.max()
        self.assertEqual(_min, "d")
        self.assertEqual(_max, "a")

        cat = Series(Categorical([np.nan,"b","c",np.nan], categories=['d','c','b','a'], ordered=True))
        _min = cat.min()
        _max = cat.max()
        self.assertTrue(np.isnan(_min))
        self.assertEqual(_max, "b")

        cat = Series(Categorical([np.nan,1,2,np.nan], categories=[5,4,3,2,1], ordered=True))
        _min = cat.min()
        _max = cat.max()
        self.assertTrue(np.isnan(_min))
        self.assertEqual(_max, 1)

    def test_mode(self):
        s = Series(Categorical([1,1,2,4,5,5,5], categories=[5,4,3,2,1], ordered=True))
        res = s.mode()
        exp = Series(Categorical([5], categories=[5,4,3,2,1], ordered=True))
        tm.assert_series_equal(res, exp)
        s = Series(Categorical([1,1,1,4,5,5,5], categories=[5,4,3,2,1], ordered=True))
        res = s.mode()
        exp = Series(Categorical([5,1], categories=[5,4,3,2,1], ordered=True))
        tm.assert_series_equal(res, exp)
        s = Series(Categorical([1,2,3,4,5], categories=[5,4,3,2,1], ordered=True))
        res = s.mode()
        exp = Series(Categorical([], categories=[5,4,3,2,1], ordered=True))
        tm.assert_series_equal(res, exp)

    def test_value_counts(self):

        s = pd.Series(pd.Categorical(["a","b","c","c","c","b"], categories=["c","a","b","d"]))
        res = s.value_counts(sort=False)
        exp = Series([3,1,2,0], index=["c","a","b","d"])
        tm.assert_series_equal(res, exp)
        res = s.value_counts(sort=True)
        exp = Series([3,2,1,0], index=["c","b","a","d"])
        tm.assert_series_equal(res, exp)

    def test_groupby(self):

        cats = Categorical(["a", "a", "a", "b", "b", "b", "c", "c", "c"], categories=["a","b","c","d"])
        data = DataFrame({"a":[1,1,1,2,2,2,3,4,5], "b":cats})

        expected = DataFrame({ 'a' : Series([1,2,4,np.nan],index=Index(['a','b','c','d'],name='b')) })
        result = data.groupby("b").mean()
        tm.assert_frame_equal(result, expected)

        raw_cat1 = Categorical(["a","a","b","b"], categories=["a","b","z"])
        raw_cat2 = Categorical(["c","d","c","d"], categories=["c","d","y"])
        df = DataFrame({"A":raw_cat1,"B":raw_cat2, "values":[1,2,3,4]})

        # single grouper
        gb = df.groupby("A")
        expected = DataFrame({ 'values' : Series([3,7,np.nan],index=Index(['a','b','z'],name='A')) })
        result = gb.sum()
        tm.assert_frame_equal(result, expected)

        # multiple groupers
        gb = df.groupby(['A','B'])
        expected = DataFrame({ 'values' : Series([1,2,np.nan,3,4,np.nan,np.nan,np.nan,np.nan],
                                                 index=pd.MultiIndex.from_product([['a','b','z'],['c','d','y']],names=['A','B'])) })
        result = gb.sum()
        tm.assert_frame_equal(result, expected)

        # multiple groupers with a non-cat
        df = df.copy()
        df['C'] = ['foo','bar']*2
        gb = df.groupby(['A','B','C'])
        expected = DataFrame({ 'values' :
                               Series(np.nan,index=pd.MultiIndex.from_product([['a','b','z'],
                                                                               ['c','d','y'],
                                                                               ['foo','bar']],
                                                                              names=['A','B','C']))
                               }).sortlevel()
        expected.iloc[[1,2,7,8],0] = [1,2,3,4]
        result = gb.sum()
        tm.assert_frame_equal(result, expected)

        # GH 8623
        x=pd.DataFrame([[1,'John P. Doe'],[2,'Jane Dove'],[1,'John P. Doe']],
                       columns=['person_id','person_name'])
        x['person_name'] = pd.Categorical(x.person_name)

        g = x.groupby(['person_id'])
        result = g.transform(lambda x:x)
        tm.assert_frame_equal(result, x[['person_name']])

        result = x.drop_duplicates('person_name')
        expected = x.iloc[[0,1]]
        tm.assert_frame_equal(result, expected)

        def f(x):
            return x.drop_duplicates('person_name').iloc[0]

        result = g.apply(f)
        expected = x.iloc[[0,1]].copy()
        expected.index = Index([1,2],name='person_id')
        expected['person_name'] = expected['person_name'].astype('object')
        tm.assert_frame_equal(result, expected)

    def test_pivot_table(self):

        raw_cat1 = Categorical(["a","a","b","b"], categories=["a","b","z"])
        raw_cat2 = Categorical(["c","d","c","d"], categories=["c","d","y"])
        df = DataFrame({"A":raw_cat1,"B":raw_cat2, "values":[1,2,3,4]})
        result = pd.pivot_table(df, values='values', index=['A', 'B'])

        expected = Series([1,2,np.nan,3,4,np.nan,np.nan,np.nan,np.nan],
                          index=pd.MultiIndex.from_product([['a','b','z'],['c','d','y']],names=['A','B']),
                          name='values')
        tm.assert_series_equal(result, expected)

    def test_count(self):

        s = Series(Categorical([np.nan,1,2,np.nan], categories=[5,4,3,2,1], ordered=True))
        result = s.count()
        self.assertEqual(result, 2)

    def test_sort(self):

        # unordered cats are not sortable
        cat = Series(Categorical(["a","b","b","a"], ordered=False))
        self.assertRaises(TypeError, lambda : cat.sort())

        cat = Series(Categorical(["a","c","b","d"], ordered=True))

        res = cat.order()
        exp = np.array(["a","b","c","d"])
        self.assert_numpy_array_equal(res.__array__(), exp)

        cat = Series(Categorical(["a","c","b","d"], categories=["a","b","c","d"], ordered=True))
        res = cat.order()
        exp = np.array(["a","b","c","d"])
        self.assert_numpy_array_equal(res.__array__(), exp)

        res = cat.order(ascending=False)
        exp = np.array(["d","c","b","a"])
        self.assert_numpy_array_equal(res.__array__(), exp)

        raw_cat1 = Categorical(["a","b","c","d"], categories=["a","b","c","d"], ordered=False)
        raw_cat2 = Categorical(["a","b","c","d"], categories=["d","c","b","a"])
        s = ["a","b","c","d"]
        df = DataFrame({"unsort":raw_cat1,"sort":raw_cat2, "string":s, "values":[1,2,3,4]})

        # Cats must be sorted in a dataframe
        res = df.sort(columns=["string"], ascending=False)
        exp = np.array(["d", "c", "b", "a"])
        self.assert_numpy_array_equal(res["sort"].values.__array__(), exp)
        self.assertEqual(res["sort"].dtype, "category")

        res = df.sort(columns=["sort"], ascending=False)
        exp = df.sort(columns=["string"], ascending=True)
        self.assert_numpy_array_equal(res["values"], exp["values"])
        self.assertEqual(res["sort"].dtype, "category")
        self.assertEqual(res["unsort"].dtype, "category")

        def f():
            df.sort(columns=["unsort"], ascending=False)
        self.assertRaises(TypeError, f)

        # multi-columns sort
        # GH 7848
        df = DataFrame({"id":[6,5,4,3,2,1], "raw_grade":['a', 'b', 'b', 'a', 'a', 'e']})
        df["grade"] = pd.Categorical(df["raw_grade"])
        df['grade'] = df['grade'].cat.set_categories(['b', 'e', 'a'])

        # sorts 'grade' according to the order of the categories
        result = df.sort(columns=['grade'])
        expected = df.iloc[[1,2,5,0,3,4]]
        tm.assert_frame_equal(result,expected)

        # multi
        result = df.sort(columns=['grade', 'id'])
        expected = df.iloc[[2,1,5,4,3,0]]
        tm.assert_frame_equal(result,expected)

        # reverse
        cat = Categorical(["a","c","c","b","d"], ordered=True)
        res = cat.order(ascending=False)
        exp_val = np.array(["d","c", "c", "b","a"],dtype=object)
        exp_categories = np.array(["a","b","c","d"],dtype=object)
        self.assert_numpy_array_equal(res.__array__(), exp_val)
        self.assert_numpy_array_equal(res.categories, exp_categories)

        # some NaN positions

        cat = Categorical(["a","c","b","d", np.nan], ordered=True)
        res = cat.order(ascending=False, na_position='last')
        exp_val = np.array(["d","c","b","a", np.nan],dtype=object)
        exp_categories = np.array(["a","b","c","d"],dtype=object)
        self.assert_numpy_array_equal(res.__array__(), exp_val)
        self.assert_numpy_array_equal(res.categories, exp_categories)

        cat = Categorical(["a","c","b","d", np.nan], ordered=True)
        res = cat.order(ascending=False, na_position='first')
        exp_val = np.array([np.nan, "d","c","b","a"],dtype=object)
        exp_categories = np.array(["a","b","c","d"],dtype=object)
        self.assert_numpy_array_equal(res.__array__(), exp_val)
        self.assert_numpy_array_equal(res.categories, exp_categories)

        cat = Categorical(["a","c","b","d", np.nan], ordered=True)
        res = cat.order(ascending=False, na_position='first')
        exp_val = np.array([np.nan, "d","c","b","a"],dtype=object)
        exp_categories = np.array(["a","b","c","d"],dtype=object)
        self.assert_numpy_array_equal(res.__array__(), exp_val)
        self.assert_numpy_array_equal(res.categories, exp_categories)

        cat = Categorical(["a","c","b","d", np.nan], ordered=True)
        res = cat.order(ascending=False, na_position='last')
        exp_val = np.array(["d","c","b","a",np.nan],dtype=object)
        exp_categories = np.array(["a","b","c","d"],dtype=object)
        self.assert_numpy_array_equal(res.__array__(), exp_val)
        self.assert_numpy_array_equal(res.categories, exp_categories)

    def test_slicing(self):
        cat = Series(Categorical([1,2,3,4]))
        reversed = cat[::-1]
        exp = np.array([4,3,2,1])
        self.assert_numpy_array_equal(reversed.__array__(), exp)

        df = DataFrame({'value': (np.arange(100)+1).astype('int64')})
        df['D'] = pd.cut(df.value, bins=[0,25,50,75,100])

        expected = Series([11,'(0, 25]'],index=['value','D'])
        result = df.iloc[10]
        tm.assert_series_equal(result,expected)

        expected = DataFrame({'value': np.arange(11,21).astype('int64')},
                             index=np.arange(10,20).astype('int64'))
        expected['D'] = pd.cut(expected.value, bins=[0,25,50,75,100])
        result = df.iloc[10:20]
        tm.assert_frame_equal(result,expected)

        expected = Series([9,'(0, 25]'],index=['value','D'])
        result = df.loc[8]
        tm.assert_series_equal(result,expected)

    def test_slicing_and_getting_ops(self):

        # systematically test the slicing operations:
        #  for all slicing ops:
        #   - returning a dataframe
        #   - returning a column
        #   - returning a row
        #   - returning a single value

        cats = pd.Categorical(["a","c","b","c","c","c","c"], categories=["a","b","c"])
        idx = pd.Index(["h","i","j","k","l","m","n"])
        values= [1,2,3,4,5,6,7]
        df = pd.DataFrame({"cats":cats,"values":values}, index=idx)

        # the expected values
        cats2 = pd.Categorical(["b","c"], categories=["a","b","c"])
        idx2 = pd.Index(["j","k"])
        values2= [3,4]

        # 2:4,: | "j":"k",:
        exp_df = pd.DataFrame({"cats":cats2,"values":values2}, index=idx2)

        # :,"cats" | :,0
        exp_col = pd.Series(cats,index=idx,name='cats')

        # "j",: | 2,:
        exp_row = pd.Series(["b",3], index=["cats","values"], dtype="object", name="j")

        # "j","cats | 2,0
        exp_val = "b"

        # iloc
        # frame
        res_df = df.iloc[2:4,:]
        tm.assert_frame_equal(res_df, exp_df)
        self.assertTrue(com.is_categorical_dtype(res_df["cats"]))

        # row
        res_row = df.iloc[2,:]
        tm.assert_series_equal(res_row, exp_row)
        tm.assert_isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.iloc[:,0]
        tm.assert_series_equal(res_col, exp_col)
        self.assertTrue(com.is_categorical_dtype(res_col))

        # single value
        res_val = df.iloc[2,0]
        self.assertEqual(res_val, exp_val)

        # loc
        # frame
        res_df = df.loc["j":"k",:]
        tm.assert_frame_equal(res_df, exp_df)
        self.assertTrue(com.is_categorical_dtype(res_df["cats"]))

        # row
        res_row = df.loc["j",:]
        tm.assert_series_equal(res_row, exp_row)
        tm.assert_isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.loc[:,"cats"]
        tm.assert_series_equal(res_col, exp_col)
        self.assertTrue(com.is_categorical_dtype(res_col))

        # single value
        res_val = df.loc["j","cats"]
        self.assertEqual(res_val, exp_val)

        # ix
        # frame
        #res_df = df.ix["j":"k",[0,1]] # doesn't work?
        res_df = df.ix["j":"k",:]
        tm.assert_frame_equal(res_df, exp_df)
        self.assertTrue(com.is_categorical_dtype(res_df["cats"]))

        # row
        res_row = df.ix["j",:]
        tm.assert_series_equal(res_row, exp_row)
        tm.assert_isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.ix[:,"cats"]
        tm.assert_series_equal(res_col, exp_col)
        self.assertTrue(com.is_categorical_dtype(res_col))

        # single value
        res_val = df.ix["j",0]
        self.assertEqual(res_val, exp_val)

        # iat
        res_val = df.iat[2,0]
        self.assertEqual(res_val, exp_val)

        # at
        res_val = df.at["j","cats"]
        self.assertEqual(res_val, exp_val)

        # fancy indexing
        exp_fancy = df.iloc[[2]]

        res_fancy = df[df["cats"] == "b"]
        tm.assert_frame_equal(res_fancy,exp_fancy)
        res_fancy = df[df["values"] == 3]
        tm.assert_frame_equal(res_fancy,exp_fancy)

        # get_value
        res_val = df.get_value("j","cats")
        self.assertEqual(res_val, exp_val)

        # i : int, slice, or sequence of integers
        res_row = df.irow(2)
        tm.assert_series_equal(res_row, exp_row)
        tm.assert_isinstance(res_row["cats"], compat.string_types)

        res_df = df.irow(slice(2,4))
        tm.assert_frame_equal(res_df, exp_df)
        self.assertTrue(com.is_categorical_dtype(res_df["cats"]))

        res_df = df.irow([2,3])
        tm.assert_frame_equal(res_df, exp_df)
        self.assertTrue(com.is_categorical_dtype(res_df["cats"]))

        res_col = df.icol(0)
        tm.assert_series_equal(res_col, exp_col)
        self.assertTrue(com.is_categorical_dtype(res_col))

        res_df = df.icol(slice(0,2))
        tm.assert_frame_equal(res_df, df)
        self.assertTrue(com.is_categorical_dtype(res_df["cats"]))

        res_df = df.icol([0,1])
        tm.assert_frame_equal(res_df, df)
        self.assertTrue(com.is_categorical_dtype(res_df["cats"]))

    def test_slicing_doc_examples(self):

        #GH 7918
        cats = Categorical(["a","b","b","b","c","c","c"], categories=["a","b","c"])
        idx = Index(["h","i","j","k","l","m","n",])
        values= [1,2,2,2,3,4,5]
        df = DataFrame({"cats":cats,"values":values}, index=idx)

        result = df.iloc[2:4,:]
        expected = DataFrame({"cats":Categorical(['b','b'],categories=['a','b','c']),"values":[2,2]}, index=['j','k'])
        tm.assert_frame_equal(result, expected)

        result = df.iloc[2:4,:].dtypes
        expected = Series(['category','int64'],['cats','values'])
        tm.assert_series_equal(result, expected)

        result = df.loc["h":"j","cats"]
        expected = Series(Categorical(['a','b','b'],categories=['a','b','c']),index=['h','i','j'])
        tm.assert_series_equal(result, expected)

        result = df.ix["h":"j",0:1]
        expected = DataFrame({'cats' : Series(Categorical(['a','b','b'],categories=['a','b','c']),index=['h','i','j']) })
        tm.assert_frame_equal(result, expected)

    def test_assigning_ops(self):

        # systematically test the assigning operations:
        # for all slicing ops:
        #  for value in categories and value not in categories:
        #   - assign a single value -> exp_single_cats_value
        #   - assign a complete row (mixed values) -> exp_single_row
        #   - assign multiple rows (mixed values) (-> array) -> exp_multi_row
        #   - assign a part of a column with dtype == categorical -> exp_parts_cats_col
        #   - assign a part of a column with dtype != categorical -> exp_parts_cats_col

        cats = pd.Categorical(["a","a","a","a","a","a","a"], categories=["a","b"])
        idx = pd.Index(["h","i","j","k","l","m","n"])
        values = [1,1,1,1,1,1,1]
        orig = pd.DataFrame({"cats":cats,"values":values}, index=idx)

        ### the expected values
        # changed single row
        cats1 = pd.Categorical(["a","a","b","a","a","a","a"], categories=["a","b"])
        idx1 = pd.Index(["h","i","j","k","l","m","n"])
        values1 = [1,1,2,1,1,1,1]
        exp_single_row = pd.DataFrame({"cats":cats1,"values":values1}, index=idx1)

        #changed multiple rows
        cats2 = pd.Categorical(["a","a","b","b","a","a","a"], categories=["a","b"])
        idx2 = pd.Index(["h","i","j","k","l","m","n"])
        values2 = [1,1,2,2,1,1,1]
        exp_multi_row = pd.DataFrame({"cats":cats2,"values":values2}, index=idx2)

        # changed part of the cats column
        cats3 = pd.Categorical(["a","a","b","b","a","a","a"], categories=["a","b"])
        idx3 = pd.Index(["h","i","j","k","l","m","n"])
        values3 = [1,1,1,1,1,1,1]
        exp_parts_cats_col = pd.DataFrame({"cats":cats3,"values":values3}, index=idx3)

        # changed single value in cats col
        cats4 = pd.Categorical(["a","a","b","a","a","a","a"], categories=["a","b"])
        idx4 = pd.Index(["h","i","j","k","l","m","n"])
        values4 = [1,1,1,1,1,1,1]
        exp_single_cats_value = pd.DataFrame({"cats":cats4,"values":values4}, index=idx4)

        ####  iloc #####
        ################
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.iloc[2,0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)


        df = orig.copy()
        df.iloc[df.index == "j",0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)


        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.iloc[2,0] = "c"
        self.assertRaises(ValueError, f)

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.iloc[2,:] = ["b",2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        def f():
            df = orig.copy()
            df.iloc[2,:] = ["c",2]
        self.assertRaises(ValueError, f)

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.iloc[2:4,:] = [["b",2],["b",2]]
        tm.assert_frame_equal(df, exp_multi_row)

        def f():
            df = orig.copy()
            df.iloc[2:4,:] = [["c",2],["c",2]]
        self.assertRaises(ValueError, f)

        #   - assign a part of a column with dtype == categorical -> exp_parts_cats_col
        df = orig.copy()
        df.iloc[2:4,0] = pd.Categorical(["b","b"], categories=["a","b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with tm.assertRaises(ValueError):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.iloc[2:4,0] = pd.Categorical(["b","b"], categories=["a","b","c"])

        with tm.assertRaises(ValueError):
            # different values
            df = orig.copy()
            df.iloc[2:4,0] = pd.Categorical(["c","c"], categories=["a","b","c"])

        #   - assign a part of a column with dtype != categorical -> exp_parts_cats_col
        df = orig.copy()
        df.iloc[2:4,0] = ["b","b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with tm.assertRaises(ValueError):
            df.iloc[2:4,0] = ["c","c"]

        ####  loc  #####
        ################
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.loc["j","cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        df = orig.copy()
        df.loc[df.index == "j","cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.loc["j","cats"] = "c"
        self.assertRaises(ValueError, f)

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.loc["j",:] = ["b",2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        def f():
            df = orig.copy()
            df.loc["j",:] = ["c",2]
        self.assertRaises(ValueError, f)

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.loc["j":"k",:] = [["b",2],["b",2]]
        tm.assert_frame_equal(df, exp_multi_row)

        def f():
            df = orig.copy()
            df.loc["j":"k",:] = [["c",2],["c",2]]
        self.assertRaises(ValueError, f)

        #   - assign a part of a column with dtype == categorical -> exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k","cats"] = pd.Categorical(["b","b"], categories=["a","b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with tm.assertRaises(ValueError):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.loc["j":"k","cats"] = pd.Categorical(["b","b"], categories=["a","b","c"])

        with tm.assertRaises(ValueError):
            # different values
            df = orig.copy()
            df.loc["j":"k","cats"] = pd.Categorical(["c","c"], categories=["a","b","c"])

        #   - assign a part of a column with dtype != categorical -> exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k","cats"] = ["b","b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with tm.assertRaises(ValueError):
            df.loc["j":"k","cats"] = ["c","c"]

        ####  ix   #####
        ################
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.ix["j",0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        df = orig.copy()
        df.ix[df.index == "j",0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.ix["j",0] = "c"
        self.assertRaises(ValueError, f)

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.ix["j",:] = ["b",2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        def f():
            df = orig.copy()
            df.ix["j",:] = ["c",2]
        self.assertRaises(ValueError, f)

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.ix["j":"k",:] = [["b",2],["b",2]]
        tm.assert_frame_equal(df, exp_multi_row)

        def f():
            df = orig.copy()
            df.ix["j":"k",:] = [["c",2],["c",2]]
        self.assertRaises(ValueError, f)

        #   - assign a part of a column with dtype == categorical -> exp_parts_cats_col
        df = orig.copy()
        df.ix["j":"k",0] = pd.Categorical(["b","b"], categories=["a","b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with tm.assertRaises(ValueError):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.ix["j":"k",0] = pd.Categorical(["b","b"], categories=["a","b","c"])

        with tm.assertRaises(ValueError):
            # different values
            df = orig.copy()
            df.ix["j":"k",0] = pd.Categorical(["c","c"], categories=["a","b","c"])

        #   - assign a part of a column with dtype != categorical -> exp_parts_cats_col
        df = orig.copy()
        df.ix["j":"k",0] = ["b","b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with tm.assertRaises(ValueError):
            df.ix["j":"k",0] = ["c","c"]

        # iat
        df = orig.copy()
        df.iat[2,0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.iat[2,0] = "c"
        self.assertRaises(ValueError, f)

        # at
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.at["j","cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.at["j","cats"] = "c"
        self.assertRaises(ValueError, f)

        # fancy indexing
        catsf = pd.Categorical(["a","a","c","c","a","a","a"], categories=["a","b","c"])
        idxf = pd.Index(["h","i","j","k","l","m","n"])
        valuesf = [1,1,3,3,1,1,1]
        df = pd.DataFrame({"cats":catsf,"values":valuesf}, index=idxf)

        exp_fancy = exp_multi_row.copy()
        exp_fancy["cats"].cat.set_categories(["a","b","c"], inplace=True)

        df[df["cats"] == "c"] = ["b",2]
        tm.assert_frame_equal(df, exp_multi_row)

        # set_value
        df = orig.copy()
        df.set_value("j","cats", "b")
        tm.assert_frame_equal(df, exp_single_cats_value)

        def f():
            df = orig.copy()
            df.set_value("j","cats", "c")
        self.assertRaises(ValueError, f)

        # Assigning a Category to parts of a int/... column uses the values of the Catgorical
        df = pd.DataFrame({"a":[1,1,1,1,1], "b":["a","a","a","a","a"]})
        exp = pd.DataFrame({"a":[1,"b","b",1,1], "b":["a","a","b","b","a"]})
        df.loc[1:2,"a"] = pd.Categorical(["b","b"], categories=["a","b"])
        df.loc[2:3,"b"] = pd.Categorical(["b","b"], categories=["a","b"])
        tm.assert_frame_equal(df, exp)

        ######### Series ##########
        orig = Series(pd.Categorical(["b","b"], categories=["a","b"]))
        s = orig.copy()
        s[:] = "a"
        exp = Series(pd.Categorical(["a","a"], categories=["a","b"]))
        tm.assert_series_equal(s, exp)

        s = orig.copy()
        s[1] = "a"
        exp = Series(pd.Categorical(["b","a"], categories=["a","b"]))
        tm.assert_series_equal(s, exp)

        s = orig.copy()
        s[s.index > 0] = "a"
        exp = Series(pd.Categorical(["b","a"], categories=["a","b"]))
        tm.assert_series_equal(s, exp)

        s = orig.copy()
        s[[False, True]] = "a"
        exp = Series(pd.Categorical(["b","a"], categories=["a","b"]))
        tm.assert_series_equal(s, exp)

        s = orig.copy()
        s.index = ["x", "y"]
        s["y"] = "a"
        exp = Series(pd.Categorical(["b","a"], categories=["a","b"]), index=["x", "y"])
        tm.assert_series_equal(s, exp)

        # ensure that one can set something to np.nan
        s = Series(Categorical([1,2,3]))
        exp = Series(Categorical([1,np.nan,3]))
        s[1] = np.nan
        tm.assert_series_equal(s, exp)


    def test_comparisons(self):
        tests_data = [(list("abc"), list("cba"), list("bbb")),
                      ([1,2,3], [3,2,1], [2,2,2])]
        for data , reverse, base in tests_data:
            cat_rev = pd.Series(pd.Categorical(data, categories=reverse))
            cat_rev_base = pd.Series(pd.Categorical(base, categories=reverse))
            cat = pd.Series(pd.Categorical(data))
            cat_base = pd.Series(pd.Categorical(base, categories=cat.cat.categories))
            s = Series(base)
            a = np.array(base)

            # comparisons need to take categories ordering into account
            res_rev = cat_rev > cat_rev_base
            exp_rev = Series([True, False, False])
            tm.assert_series_equal(res_rev, exp_rev)

            res_rev = cat_rev < cat_rev_base
            exp_rev = Series([False, False, True])
            tm.assert_series_equal(res_rev, exp_rev)

            res = cat > cat_base
            exp = Series([False, False, True])
            tm.assert_series_equal(res, exp)

            # Only categories with same categories can be compared
            def f():
                cat > cat_rev
            self.assertRaises(TypeError, f)

            # categorical cannot be compared to Series or numpy array, and also not the other way
            # around
            self.assertRaises(TypeError, lambda: cat > s)
            self.assertRaises(TypeError, lambda: cat_rev > s)
            self.assertRaises(TypeError, lambda: cat > a)
            self.assertRaises(TypeError, lambda: cat_rev > a)

            self.assertRaises(TypeError, lambda: s < cat)
            self.assertRaises(TypeError, lambda: s < cat_rev)

            self.assertRaises(TypeError, lambda: a < cat)
            self.assertRaises(TypeError, lambda: a < cat_rev)

            # Categoricals can be compared to scalar values
            res = cat_rev > base[0]
            tm.assert_series_equal(res, exp)

        # And test NaN handling...
        cat = Series(Categorical(["a","b","c", np.nan]))
        exp = Series([True, True, True, False])
        res = (cat == cat)
        tm.assert_series_equal(res, exp)

    def test_cat_equality(self):

        # GH 8938
        # allow equality comparisons
        a = Series(list('abc'),dtype="category")
        b = Series(list('abc'),dtype="object")
        c = Series(['a','b','cc'],dtype="object")
        d = Series(list('acb'),dtype="object")
        e = Categorical(list('abc'))
        f = Categorical(list('acb'))

        # vs scalar
        self.assertFalse((a=='a').all())
        self.assertTrue(((a!='a') == ~(a=='a')).all())

        self.assertFalse(('a'==a).all())
        self.assertTrue((a=='a')[0])
        self.assertTrue(('a'==a)[0])
        self.assertFalse(('a'!=a)[0])

        # vs list-like
        self.assertTrue((a==a).all())
        self.assertFalse((a!=a).all())

        self.assertTrue((a==list(a)).all())
        self.assertTrue((a==b).all())
        self.assertTrue((b==a).all())
        self.assertTrue(((~(a==b))==(a!=b)).all())
        self.assertTrue(((~(b==a))==(b!=a)).all())

        self.assertFalse((a==c).all())
        self.assertFalse((c==a).all())
        self.assertFalse((a==d).all())
        self.assertFalse((d==a).all())

        # vs a cat-like
        self.assertTrue((a==e).all())
        self.assertTrue((e==a).all())
        self.assertFalse((a==f).all())
        self.assertFalse((f==a).all())

        self.assertTrue(((~(a==e)==(a!=e)).all()))
        self.assertTrue(((~(e==a)==(e!=a)).all()))
        self.assertTrue(((~(a==f)==(a!=f)).all()))
        self.assertTrue(((~(f==a)==(f!=a)).all()))

        # non-equality is not comparable
        self.assertRaises(TypeError, lambda: a < b)
        self.assertRaises(TypeError, lambda: b < a)
        self.assertRaises(TypeError, lambda: a > b)
        self.assertRaises(TypeError, lambda: b > a)

    def test_concat(self):
        cat = pd.Categorical(["a","b"], categories=["a","b"])
        vals = [1,2]
        df = pd.DataFrame({"cats":cat, "vals":vals})
        cat2 = pd.Categorical(["a","b","a","b"], categories=["a","b"])
        vals2 = [1,2,1,2]
        exp = pd.DataFrame({"cats":cat2, "vals":vals2}, index=pd.Index([0, 1, 0, 1]))

        res = pd.concat([df,df])
        tm.assert_frame_equal(exp, res)

        # Concat should raise if the two categoricals do not have the same categories
        cat3 = pd.Categorical(["a","b"], categories=["a","b","c"])
        vals3 = [1,2]
        df_wrong_categories = pd.DataFrame({"cats":cat3, "vals":vals3})

        def f():
            pd.concat([df,df_wrong_categories])
        self.assertRaises(ValueError, f)

        # GH 7864
        # make sure ordering is preserverd
        df = pd.DataFrame({"id":[1,2,3,4,5,6], "raw_grade":['a', 'b', 'b', 'a', 'a', 'e']})
        df["grade"] = pd.Categorical(df["raw_grade"])
        df['grade'].cat.set_categories(['e', 'a', 'b'])

        df1 = df[0:3]
        df2 = df[3:]

        self.assert_numpy_array_equal(df['grade'].cat.categories, df1['grade'].cat.categories)
        self.assert_numpy_array_equal(df['grade'].cat.categories, df2['grade'].cat.categories)

        dfx = pd.concat([df1, df2])
        dfx['grade'].cat.categories
        self.assert_numpy_array_equal(df['grade'].cat.categories, dfx['grade'].cat.categories)

        # GH 8641
        # series concat not preserving category dtype
        s = Series(list('abc'),dtype='category')
        s2 = Series(list('abd'),dtype='category')

        def f():
            pd.concat([s,s2])
        self.assertRaises(ValueError, f)

        result = pd.concat([s,s],ignore_index=True)
        expected = Series(list('abcabc')).astype('category')
        tm.assert_series_equal(result, expected)

        result = pd.concat([s,s])
        expected = Series(list('abcabc'),index=[0,1,2,0,1,2]).astype('category')
        tm.assert_series_equal(result, expected)

    def test_append(self):
        cat = pd.Categorical(["a","b"], categories=["a","b"])
        vals = [1,2]
        df = pd.DataFrame({"cats":cat, "vals":vals})
        cat2 = pd.Categorical(["a","b","a","b"], categories=["a","b"])
        vals2 = [1,2,1,2]
        exp = pd.DataFrame({"cats":cat2, "vals":vals2}, index=pd.Index([0, 1, 0, 1]))

        res = df.append(df)
        tm.assert_frame_equal(exp, res)

        # Concat should raise if the two categoricals do not have the same categories
        cat3 = pd.Categorical(["a","b"], categories=["a","b","c"])
        vals3 = [1,2]
        df_wrong_categories = pd.DataFrame({"cats":cat3, "vals":vals3})

        def f():
            df.append(df_wrong_categories)
        self.assertRaises(ValueError, f)

    def test_na_actions(self):

        cat = pd.Categorical([1,2,3,np.nan], categories=[1,2,3])
        vals = ["a","b",np.nan,"d"]
        df = pd.DataFrame({"cats":cat, "vals":vals})
        cat2 = pd.Categorical([1,2,3,3], categories=[1,2,3])
        vals2 = ["a","b","b","d"]
        df_exp_fill = pd.DataFrame({"cats":cat2, "vals":vals2})
        cat3 = pd.Categorical([1,2,3], categories=[1,2,3])
        vals3 = ["a","b",np.nan]
        df_exp_drop_cats = pd.DataFrame({"cats":cat3, "vals":vals3})
        cat4 = pd.Categorical([1,2], categories=[1,2,3])
        vals4 = ["a","b"]
        df_exp_drop_all = pd.DataFrame({"cats":cat4, "vals":vals4})

        # fillna
        res = df.fillna(value={"cats":3, "vals":"b"})
        tm.assert_frame_equal(res, df_exp_fill)

        def f():
            df.fillna(value={"cats":4, "vals":"c"})
        self.assertRaises(ValueError, f)

        res = df.fillna(method='pad')
        tm.assert_frame_equal(res, df_exp_fill)

        res = df.dropna(subset=["cats"])
        tm.assert_frame_equal(res, df_exp_drop_cats)

        res = df.dropna()
        tm.assert_frame_equal(res, df_exp_drop_all)

        # make sure that fillna takes both missing values and NA categories into account
        c = Categorical(["a","b",np.nan])
        c.set_categories(["a","b",np.nan], rename=True, inplace=True)
        c[0] = np.nan
        df = pd.DataFrame({"cats":c, "vals":[1,2,3]})
        df_exp = pd.DataFrame({"cats": Categorical(["a","b","a"]), "vals": [1,2,3]})
        res = df.fillna("a")
        tm.assert_frame_equal(res, df_exp)


    def test_astype_to_other(self):

        s = self.cat['value_group']
        expected = s
        tm.assert_series_equal(s.astype('category'),expected)
        tm.assert_series_equal(s.astype(com.CategoricalDtype()),expected)
        self.assertRaises(ValueError, lambda : s.astype('float64'))

        cat = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c']))
        exp = Series(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])
        tm.assert_series_equal(cat.astype('str'), exp)
        s2 = Series(Categorical.from_array(['1', '2', '3', '4']))
        exp2 = Series([1,2,3,4]).astype(int)
        tm.assert_series_equal(s2.astype('int') , exp2)

        # object don't sort correctly, so just compare that we have the same values
        def cmp(a,b):
            tm.assert_almost_equal(np.sort(np.unique(a)),np.sort(np.unique(b)))
        expected = Series(np.array(s.values),name='value_group')
        cmp(s.astype('object'),expected)
        cmp(s.astype(np.object_),expected)

        # array conversion
        tm.assert_almost_equal(np.array(s),np.array(s.values))

        # valid conversion
        for valid in [lambda x: x.astype('category'),
                      lambda x: x.astype(com.CategoricalDtype()),
                      lambda x: x.astype('object').astype('category'),
                      lambda x: x.astype('object').astype(com.CategoricalDtype())]:

            result = valid(s)
            tm.assert_series_equal(result,s)

        # invalid conversion (these are NOT a dtype)
        for invalid in [lambda x: x.astype(pd.Categorical),
                        lambda x: x.astype('object').astype(pd.Categorical)]:
            self.assertRaises(TypeError, lambda : invalid(s))


    def test_to_records(self):

        # GH8626

        # dict creation
        df = DataFrame({ 'A' : list('abc') },dtype='category')
        expected = Series(list('abc'),dtype='category')
        tm.assert_series_equal(df['A'],expected)

        # list-like creation
        df = DataFrame(list('abc'),dtype='category')
        expected = Series(list('abc'),dtype='category')
        tm.assert_series_equal(df[0],expected)

        # to record array
        # this coerces
        result = df.to_records()
        expected = np.rec.array([(0, 'a'), (1, 'b'), (2, 'c')],
                                dtype=[('index', '<i8'), ('0', 'O')])
        tm.assert_almost_equal(result,expected)

    def test_numeric_like_ops(self):

        # numeric ops should not succeed
        for op in ['__add__','__sub__','__mul__','__truediv__']:
            self.assertRaises(TypeError, lambda : getattr(self.cat,op)(self.cat))

        # reduction ops should not succeed (unless specifically defined, e.g. min/max)
        s = self.cat['value_group']
        for op in ['kurt','skew','var','std','mean','sum','median']:
            self.assertRaises(TypeError, lambda : getattr(s,op)(numeric_only=False))

        # mad technically works because it takes always the numeric data

        # numpy ops
        s = pd.Series(pd.Categorical([1,2,3,4]))
        self.assertRaises(TypeError, lambda : np.sum(s))

        # numeric ops on a Series
        for op in ['__add__','__sub__','__mul__','__truediv__']:
            self.assertRaises(TypeError, lambda : getattr(s,op)(2))

        # invalid ufunc
        self.assertRaises(TypeError, lambda : np.log(s))

    def test_cat_tab_completition(self):
         # test the tab completion display
        ok_for_cat = ['categories','codes','ordered','set_categories',
                      'add_categories', 'remove_categories', 'rename_categories',
                      'reorder_categories', 'remove_unused_categories']
        def get_dir(s):
            results = [ r for r in s.cat.__dir__() if not r.startswith('_') ]
            return list(sorted(set(results)))

        s = Series(list('aabbcde')).astype('category')
        results = get_dir(s)
        tm.assert_almost_equal(results,list(sorted(set(ok_for_cat))))

    def test_pickle_v0_14_1(self):
        cat = pd.Categorical(values=['a', 'b', 'c'],
                             categories=['a', 'b', 'c', 'd'],
                             name='foobar', ordered=False)
        pickle_path = os.path.join(tm.get_data_path(),
                                   'categorical_0_14_1.pickle')
        # This code was executed once on v0.14.1 to generate the pickle:
        #
        # cat = Categorical(labels=np.arange(3), levels=['a', 'b', 'c', 'd'],
        #                   name='foobar')
        # with open(pickle_path, 'wb') as f: pickle.dump(cat, f)
        #
        self.assert_categorical_equal(cat, pd.read_pickle(pickle_path))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core']
                    exit=False)
