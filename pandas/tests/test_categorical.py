# pylint: disable=E1101,E1103,W0232

from datetime import datetime
from pandas.compat import range, lrange, u
import nose
import re

import numpy as np

from pandas.core.categorical import Categorical
from pandas.core.index import Index, Int64Index, MultiIndex
from pandas.core.frame import DataFrame
from pandas.util.testing import assert_almost_equal
import pandas.core.common as com

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
        tm.assert_almost_equal(subf.labels, [0, 1, 1])

        subf = self.factor[np.asarray(self.factor) == 'c']
        tm.assert_almost_equal(subf.labels, [2, 2, 2])

    def test_constructor_unsortable(self):
        raise nose.SkipTest('skipping for now')

        arr = np.array([1, 2, 3, datetime.now()], dtype='O')

        # it works!
        factor = Categorical.from_array(arr)

    def test_factor_agg(self):
        import pandas.core.frame as frame

        arr = np.arange(len(self.factor))

        f = np.sum
        agged = frame.factor_agg(self.factor, arr, f)
        labels = self.factor.labels
        for i, idx in enumerate(self.factor.levels):
            self.assertEqual(f(arr[labels == i]), agged[i])

    def test_comparisons(self):
        result = self.factor[self.factor == 'a']
        expected = self.factor[np.asarray(self.factor) == 'a']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor != 'a']
        expected = self.factor[np.asarray(self.factor) != 'a']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor < 'c']
        expected = self.factor[np.asarray(self.factor) < 'c']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor > 'a']
        expected = self.factor[np.asarray(self.factor) > 'a']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor >= 'b']
        expected = self.factor[np.asarray(self.factor) >= 'b']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor <= 'b']
        expected = self.factor[np.asarray(self.factor) <= 'b']
        self.assert_(result.equals(expected))

        n = len(self.factor)

        other = self.factor[np.random.permutation(n)]
        result = self.factor == other
        expected = np.asarray(self.factor) == np.asarray(other)
        self.assert_(np.array_equal(result, expected))

        result = self.factor == 'd'
        expected = np.repeat(False, len(self.factor))
        self.assert_(np.array_equal(result, expected))

    def test_na_flags_int_levels(self):
        # #1457

        levels = lrange(10)
        labels = np.random.randint(0, 10, 20)
        labels[::5] = -1

        cat = Categorical(labels, levels)
        repr(cat)

        self.assert_(np.array_equal(com.isnull(cat), labels == -1))

    def test_levels_none(self):
        factor = Categorical(['a', 'b', 'b', 'a',
                              'a', 'c', 'c', 'c'])
        self.assert_(factor.equals(self.factor))

    def test_describe(self):
        # string type
        desc = self.factor.describe()
        expected = DataFrame.from_dict(dict(counts=[3, 2, 3],
                                            freqs=[3/8., 2/8., 3/8.],
                                            levels=['a', 'b', 'c'])
                                            ).set_index('levels')
        tm.assert_frame_equal(desc, expected)

        # check an integer one
        desc = Categorical([1,2,3,1,2,3,3,2,1,1,1]).describe()
        expected = DataFrame.from_dict(dict(counts=[5, 3, 3],
                                            freqs=[5/11., 3/11., 3/11.],
                                            levels=[1,2,3]
                                            )
                                            ).set_index('levels')
        tm.assert_frame_equal(desc, expected)

    def test_print(self):
        expected = [" a", " b", " b", " a", " a", " c", " c", " c",
                    "Levels (3): Index([a, b, c], dtype=object)"]
        expected = "\n".join(expected)
        # hack because array_repr changed in numpy > 1.6.x
        actual = repr(self.factor)
        pat = "Index\(\['a', 'b', 'c']"
        sub = "Index([a, b, c]"
        actual = re.sub(pat, sub, actual)

        self.assertEquals(actual, expected)

    def test_big_print(self):
        factor = Categorical([0,1,2,0,1,2]*100, ['a', 'b', 'c'], name='cat')
        expected = [" a", " b", " c", " a", " b", " c", " a", " b", " c",
                    " a", " b", " c", " a", "...", " c", " a", " b", " c",
                    " a", " b", " c", " a", " b", " c", " a", " b", " c",
                    "Levels (3): Index([a, b, c], dtype=object)",
                    "Name: cat, Length: 600" ]
        expected = "\n".join(expected)

        # hack because array_repr changed in numpy > 1.6.x
        actual = repr(factor)
        pat = "Index\(\['a', 'b', 'c']"
        sub = "Index([a, b, c]"
        actual = re.sub(pat, sub, actual)

        self.assertEquals(actual, expected)

    def test_empty_print(self):
        factor = Categorical([], ["a","b","c"], name="cat")
        expected = ("Categorical([], Name: cat, Levels (3): "
                    "Index([a, b, c], dtype=object)")
        # hack because array_repr changed in numpy > 1.6.x
        actual = repr(factor)
        pat = "Index\(\['a', 'b', 'c']"
        sub = "Index([a, b, c]"
        actual = re.sub(pat, sub, actual)

        self.assertEqual(actual, expected)

        factor = Categorical([], ["a","b","c"])
        expected = ("Categorical([], Levels (3): "
                    "Index([a, b, c], dtype=object)")
        # hack because array_repr changed in numpy > 1.6.x
        actual = repr(factor)
        pat = "Index\(\['a', 'b', 'c']"
        sub = "Index([a, b, c]"
        actual = re.sub(pat, sub, actual)

        self.assertEqual(actual, expected)

        factor = Categorical([], [])
        expected = ("Categorical([], Levels (0): "
                    "Index([], dtype=object)")
        self.assertEqual(repr(factor), expected)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
