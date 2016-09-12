# -*- coding: utf-8 -*-

import numpy as np
from pandas import Series, DataFrame, Index, Float64Index
from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm


class TestFloatIndexers(tm.TestCase):

    def check(self, result, original, indexer, getitem):
        """
        comparator for results
        we need to take care if we are indexing on a
        Series or a frame
        """
        if isinstance(original, Series):
            expected = original.iloc[indexer]
        else:
            if getitem:
                expected = original.iloc[:, indexer]
            else:
                expected = original.iloc[indexer]

        assert_almost_equal(result, expected)

    def test_scalar_error(self):

        # GH 4892
        # float_indexers should raise exceptions
        # on appropriate Index types & accessors
        # this duplicates the code below
        # but is spefically testing for the error
        # message

        for index in [tm.makeStringIndex, tm.makeUnicodeIndex,
                      tm.makeCategoricalIndex,
                      tm.makeDateIndex, tm.makeTimedeltaIndex,
                      tm.makePeriodIndex, tm.makeIntIndex,
                      tm.makeRangeIndex]:

            i = index(5)

            s = Series(np.arange(len(i)), index=i)

            def f():
                s.iloc[3.0]
            self.assertRaisesRegexp(TypeError,
                                    'cannot do positional indexing',
                                    f)

            def f():
                s.iloc[3.0] = 0
            self.assertRaises(TypeError, f)

    def test_scalar_non_numeric(self):

        # GH 4892
        # float_indexers should raise exceptions
        # on appropriate Index types & accessors

        for index in [tm.makeStringIndex, tm.makeUnicodeIndex,
                      tm.makeCategoricalIndex,
                      tm.makeDateIndex, tm.makeTimedeltaIndex,
                      tm.makePeriodIndex]:

            i = index(5)

            for s in [Series(
                    np.arange(len(i)), index=i), DataFrame(
                        np.random.randn(
                            len(i), len(i)), index=i, columns=i)]:

                # getting
                for idxr, getitem in [(lambda x: x.ix, False),
                                      (lambda x: x.iloc, False),
                                      (lambda x: x, True)]:

                    def f():
                        idxr(s)[3.0]

                    # gettitem on a DataFrame is a KeyError as it is indexing
                    # via labels on the columns
                    if getitem and isinstance(s, DataFrame):
                        error = KeyError
                    else:
                        error = TypeError
                    self.assertRaises(error, f)

                # label based can be a TypeError or KeyError
                def f():
                    s.loc[3.0]

                if s.index.inferred_type in ['string', 'unicode', 'mixed']:
                    error = KeyError
                else:
                    error = TypeError
                self.assertRaises(error, f)

                # contains
                self.assertFalse(3.0 in s)

                # setting with a float fails with iloc
                def f():
                    s.iloc[3.0] = 0
                self.assertRaises(TypeError, f)

                # setting with an indexer
                if s.index.inferred_type in ['categorical']:
                    # Value or Type Error
                    pass
                elif s.index.inferred_type in ['datetime64', 'timedelta64',
                                               'period']:

                    # these should prob work
                    # and are inconsisten between series/dataframe ATM
                    # for idxr in [lambda x: x.ix,
                    #             lambda x: x]:
                    #    s2 = s.copy()
                    #    def f():
                    #        idxr(s2)[3.0] = 0
                    #    self.assertRaises(TypeError, f)
                    pass

                else:

                    s2 = s.copy()
                    s2.loc[3.0] = 10
                    self.assertTrue(s2.index.is_object())

                    for idxr in [lambda x: x.ix,
                                 lambda x: x]:
                        s2 = s.copy()
                        idxr(s2)[3.0] = 0
                        self.assertTrue(s2.index.is_object())

            # fallsback to position selection, series only
            s = Series(np.arange(len(i)), index=i)
            s[3]
            self.assertRaises(TypeError, lambda: s[3.0])

    def test_scalar_with_mixed(self):

        s2 = Series([1, 2, 3], index=['a', 'b', 'c'])
        s3 = Series([1, 2, 3], index=['a', 'b', 1.5])

        # lookup in a pure string index
        # with an invalid indexer
        for idxr in [lambda x: x.ix,
                     lambda x: x,
                     lambda x: x.iloc]:

            def f():
                idxr(s2)[1.0]

            self.assertRaises(TypeError, f)

        self.assertRaises(KeyError, lambda: s2.loc[1.0])

        result = s2.loc['b']
        expected = 2
        self.assertEqual(result, expected)

        # mixed index so we have label
        # indexing
        for idxr in [lambda x: x.ix,
                     lambda x: x]:

            def f():
                idxr(s3)[1.0]

            self.assertRaises(TypeError, f)

            result = idxr(s3)[1]
            expected = 2
            self.assertEqual(result, expected)

        self.assertRaises(TypeError, lambda: s3.iloc[1.0])
        self.assertRaises(KeyError, lambda: s3.loc[1.0])

        result = s3.loc[1.5]
        expected = 3
        self.assertEqual(result, expected)

    def test_scalar_integer(self):

        # test how scalar float indexers work on int indexes

        # integer index
        for index in [tm.makeIntIndex, tm.makeRangeIndex]:

            i = index(5)
            for s in [Series(np.arange(len(i))),
                      DataFrame(np.random.randn(len(i), len(i)),
                                index=i, columns=i)]:

                # coerce to equal int
                for idxr, getitem in [(lambda x: x.ix, False),
                                      (lambda x: x.loc, False),
                                      (lambda x: x, True)]:

                    result = idxr(s)[3.0]
                    self.check(result, s, 3, getitem)

                # coerce to equal int
                for idxr, getitem in [(lambda x: x.ix, False),
                                      (lambda x: x.loc, False),
                                      (lambda x: x, True)]:

                    if isinstance(s, Series):
                        compare = self.assertEqual
                        expected = 100
                    else:
                        compare = tm.assert_series_equal
                        if getitem:
                            expected = Series(100,
                                              index=range(len(s)), name=3)
                        else:
                            expected = Series(100.,
                                              index=range(len(s)), name=3)

                    s2 = s.copy()
                    idxr(s2)[3.0] = 100

                    result = idxr(s2)[3.0]
                    compare(result, expected)

                    result = idxr(s2)[3]
                    compare(result, expected)

                # contains
                # coerce to equal int
                self.assertTrue(3.0 in s)

    def test_scalar_float(self):

        # scalar float indexers work on a float index
        index = Index(np.arange(5.))
        for s in [Series(np.arange(len(index)), index=index),
                  DataFrame(np.random.randn(len(index), len(index)),
                            index=index, columns=index)]:

            # assert all operations except for iloc are ok
            indexer = index[3]
            for idxr, getitem in [(lambda x: x.ix, False),
                                  (lambda x: x.loc, False),
                                  (lambda x: x, True)]:

                # getting
                result = idxr(s)[indexer]
                self.check(result, s, 3, getitem)

                # setting
                s2 = s.copy()

                def f():
                    idxr(s2)[indexer] = expected
                result = idxr(s2)[indexer]
                self.check(result, s, 3, getitem)

                # random integer is a KeyError
                self.assertRaises(KeyError, lambda: idxr(s)[3.5])

            # contains
            self.assertTrue(3.0 in s)

            # iloc succeeds with an integer
            expected = s.iloc[3]
            s2 = s.copy()

            s2.iloc[3] = expected
            result = s2.iloc[3]
            self.check(result, s, 3, False)

            # iloc raises with a float
            self.assertRaises(TypeError, lambda: s.iloc[3.0])

            def g():
                s2.iloc[3.0] = 0
            self.assertRaises(TypeError, g)

    def test_slice_non_numeric(self):

        # GH 4892
        # float_indexers should raise exceptions
        # on appropriate Index types & accessors

        for index in [tm.makeStringIndex, tm.makeUnicodeIndex,
                      tm.makeDateIndex, tm.makeTimedeltaIndex,
                      tm.makePeriodIndex]:

            index = index(5)
            for s in [Series(range(5), index=index),
                      DataFrame(np.random.randn(5, 2), index=index)]:

                # getitem
                for l in [slice(3.0, 4),
                          slice(3, 4.0),
                          slice(3.0, 4.0)]:

                    def f():
                        s.iloc[l]
                    self.assertRaises(TypeError, f)

                    for idxr in [lambda x: x.ix,
                                 lambda x: x.loc,
                                 lambda x: x.iloc,
                                 lambda x: x]:

                        def f():
                            idxr(s)[l]
                        self.assertRaises(TypeError, f)

                # setitem
                for l in [slice(3.0, 4),
                          slice(3, 4.0),
                          slice(3.0, 4.0)]:

                    def f():
                        s.iloc[l] = 0
                    self.assertRaises(TypeError, f)

                    for idxr in [lambda x: x.ix,
                                 lambda x: x.loc,
                                 lambda x: x.iloc,
                                 lambda x: x]:
                        def f():
                            idxr(s)[l] = 0
                        self.assertRaises(TypeError, f)

    def test_slice_integer(self):

        # same as above, but for Integer based indexes
        # these coerce to a like integer
        # oob indiciates if we are out of bounds
        # of positional indexing
        for index, oob in [(tm.makeIntIndex(5), False),
                           (tm.makeRangeIndex(5), False),
                           (tm.makeIntIndex(5) + 10, True)]:

            # s is an in-range index
            s = Series(range(5), index=index)

            # getitem
            for l in [slice(3.0, 4),
                      slice(3, 4.0),
                      slice(3.0, 4.0)]:

                for idxr in [lambda x: x.loc,
                             lambda x: x.ix]:

                    result = idxr(s)[l]

                    # these are all label indexing
                    # except getitem which is positional
                    # empty
                    if oob:
                        indexer = slice(0, 0)
                    else:
                        indexer = slice(3, 5)
                    self.check(result, s, indexer, False)

                # positional indexing
                def f():
                    s[l]

                self.assertRaises(TypeError, f)

            # getitem out-of-bounds
            for l in [slice(-6, 6),
                      slice(-6.0, 6.0)]:

                for idxr in [lambda x: x.loc,
                             lambda x: x.ix]:
                    result = idxr(s)[l]

                    # these are all label indexing
                    # except getitem which is positional
                    # empty
                    if oob:
                        indexer = slice(0, 0)
                    else:
                        indexer = slice(-6, 6)
                    self.check(result, s, indexer, False)

            # positional indexing
            def f():
                s[slice(-6.0, 6.0)]

            self.assertRaises(TypeError, f)

            # getitem odd floats
            for l, res1 in [(slice(2.5, 4), slice(3, 5)),
                            (slice(2, 3.5), slice(2, 4)),
                            (slice(2.5, 3.5), slice(3, 4))]:

                for idxr in [lambda x: x.loc,
                             lambda x: x.ix]:

                    result = idxr(s)[l]
                    if oob:
                        res = slice(0, 0)
                    else:
                        res = res1

                    self.check(result, s, res, False)

                # positional indexing
                def f():
                    s[l]

                self.assertRaises(TypeError, f)

            # setitem
            for l in [slice(3.0, 4),
                      slice(3, 4.0),
                      slice(3.0, 4.0)]:

                for idxr in [lambda x: x.loc,
                             lambda x: x.ix]:
                    sc = s.copy()
                    idxr(sc)[l] = 0
                    result = idxr(sc)[l].values.ravel()
                    self.assertTrue((result == 0).all())

                # positional indexing
                def f():
                    s[l] = 0

                self.assertRaises(TypeError, f)

    def test_integer_positional_indexing(self):
        """ make sure that we are raising on positional indexing
        w.r.t. an integer index """

        s = Series(range(2, 6), index=range(2, 6))

        result = s[2:4]
        expected = s.iloc[2:4]
        assert_series_equal(result, expected)

        for idxr in [lambda x: x,
                     lambda x: x.iloc]:

            for l in [slice(2, 4.0),
                      slice(2.0, 4),
                      slice(2.0, 4.0)]:

                def f():
                    idxr(s)[l]

                self.assertRaises(TypeError, f)

    def test_slice_integer_frame_getitem(self):

        # similar to above, but on the getitem dim (of a DataFrame)
        for index in [tm.makeIntIndex, tm.makeRangeIndex]:

            index = index(5)
            s = DataFrame(np.random.randn(5, 2), index=index)

            for idxr in [lambda x: x.loc,
                         lambda x: x.ix]:

                # getitem
                for l in [slice(0.0, 1),
                          slice(0, 1.0),
                          slice(0.0, 1.0)]:

                    result = idxr(s)[l]
                    indexer = slice(0, 2)
                    self.check(result, s, indexer, False)

                    # positional indexing
                    def f():
                        s[l]

                    self.assertRaises(TypeError, f)

                # getitem out-of-bounds
                for l in [slice(-10, 10),
                          slice(-10.0, 10.0)]:

                    result = idxr(s)[l]
                    self.check(result, s, slice(-10, 10), True)

                # positional indexing
                def f():
                    s[slice(-10.0, 10.0)]

                self.assertRaises(TypeError, f)

                # getitem odd floats
                for l, res in [(slice(0.5, 1), slice(1, 2)),
                               (slice(0, 0.5), slice(0, 1)),
                               (slice(0.5, 1.5), slice(1, 2))]:

                    result = idxr(s)[l]
                    self.check(result, s, res, False)

                    # positional indexing
                    def f():
                        s[l]

                    self.assertRaises(TypeError, f)

                # setitem
                for l in [slice(3.0, 4),
                          slice(3, 4.0),
                          slice(3.0, 4.0)]:

                    sc = s.copy()
                    idxr(sc)[l] = 0
                    result = idxr(sc)[l].values.ravel()
                    self.assertTrue((result == 0).all())

                    # positional indexing
                    def f():
                        s[l] = 0

                    self.assertRaises(TypeError, f)

    def test_slice_float(self):

        # same as above, but for floats
        index = Index(np.arange(5.)) + 0.1
        for s in [Series(range(5), index=index),
                  DataFrame(np.random.randn(5, 2), index=index)]:

            for l in [slice(3.0, 4),
                      slice(3, 4.0),
                      slice(3.0, 4.0)]:

                expected = s.iloc[3:4]
                for idxr in [lambda x: x.ix,
                             lambda x: x.loc,
                             lambda x: x]:

                    # getitem
                    result = idxr(s)[l]
                    if isinstance(s, Series):
                        self.assert_series_equal(result, expected)
                    else:
                        self.assert_frame_equal(result, expected)
                    # setitem
                    s2 = s.copy()
                    idxr(s2)[l] = 0
                    result = idxr(s2)[l].values.ravel()
                    self.assertTrue((result == 0).all())

    def test_floating_index_doc_example(self):

        index = Index([1.5, 2, 3, 4.5, 5])
        s = Series(range(5), index=index)
        self.assertEqual(s[3], 2)
        self.assertEqual(s.ix[3], 2)
        self.assertEqual(s.loc[3], 2)
        self.assertEqual(s.iloc[3], 3)

    def test_floating_misc(self):

        # related 236
        # scalar/slicing of a float index
        s = Series(np.arange(5), index=np.arange(5) * 2.5, dtype=np.int64)

        # label based slicing
        result1 = s[1.0:3.0]
        result2 = s.ix[1.0:3.0]
        result3 = s.loc[1.0:3.0]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)

        # exact indexing when found
        result1 = s[5.0]
        result2 = s.loc[5.0]
        result3 = s.ix[5.0]
        self.assertEqual(result1, result2)
        self.assertEqual(result1, result3)

        result1 = s[5]
        result2 = s.loc[5]
        result3 = s.ix[5]
        self.assertEqual(result1, result2)
        self.assertEqual(result1, result3)

        self.assertEqual(s[5.0], s[5])

        # value not found (and no fallbacking at all)

        # scalar integers
        self.assertRaises(KeyError, lambda: s.loc[4])
        self.assertRaises(KeyError, lambda: s.ix[4])
        self.assertRaises(KeyError, lambda: s[4])

        # fancy floats/integers create the correct entry (as nan)
        # fancy tests
        expected = Series([2, 0], index=Float64Index([5.0, 0.0]))
        for fancy_idx in [[5.0, 0.0], np.array([5.0, 0.0])]:  # float
            assert_series_equal(s[fancy_idx], expected)
            assert_series_equal(s.loc[fancy_idx], expected)
            assert_series_equal(s.ix[fancy_idx], expected)

        expected = Series([2, 0], index=Index([5, 0], dtype='int64'))
        for fancy_idx in [[5, 0], np.array([5, 0])]:  # int
            assert_series_equal(s[fancy_idx], expected)
            assert_series_equal(s.loc[fancy_idx], expected)
            assert_series_equal(s.ix[fancy_idx], expected)

        # all should return the same as we are slicing 'the same'
        result1 = s.loc[2:5]
        result2 = s.loc[2.0:5.0]
        result3 = s.loc[2.0:5]
        result4 = s.loc[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        # previously this did fallback indexing
        result1 = s[2:5]
        result2 = s[2.0:5.0]
        result3 = s[2.0:5]
        result4 = s[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        result1 = s.ix[2:5]
        result2 = s.ix[2.0:5.0]
        result3 = s.ix[2.0:5]
        result4 = s.ix[2.1:5]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        # combined test
        result1 = s.loc[2:5]
        result2 = s.ix[2:5]
        result3 = s[2:5]

        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)

        # list selection
        result1 = s[[0.0, 5, 10]]
        result2 = s.loc[[0.0, 5, 10]]
        result3 = s.ix[[0.0, 5, 10]]
        result4 = s.iloc[[0, 2, 4]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, result4)

        result1 = s[[1.6, 5, 10]]
        result2 = s.loc[[1.6, 5, 10]]
        result3 = s.ix[[1.6, 5, 10]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series(
            [np.nan, 2, 4], index=[1.6, 5, 10]))

        result1 = s[[0, 1, 2]]
        result2 = s.ix[[0, 1, 2]]
        result3 = s.loc[[0, 1, 2]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series(
            [0.0, np.nan, np.nan], index=[0, 1, 2]))

        result1 = s.loc[[2.5, 5]]
        result2 = s.ix[[2.5, 5]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, Series([1, 2], index=[2.5, 5.0]))

        result1 = s[[2.5]]
        result2 = s.ix[[2.5]]
        result3 = s.loc[[2.5]]
        assert_series_equal(result1, result2)
        assert_series_equal(result1, result3)
        assert_series_equal(result1, Series([1], index=[2.5]))

    def test_floating_tuples(self):
        # GH13509
        s = Series([(1, 1), (2, 2), (3, 3)], index=[0.0, 0.1, 0.2], name='foo')
        result = s[0.0]
        self.assertEqual(result, (1, 1))

        s = Series([(1, 1), (2, 2), (3, 3)], index=[0.0, 0.0, 0.2], name='foo')
        result = s[0.0]
        expected = Series([(1, 1), (2, 2)], index=[0.0, 0.0], name='foo')
        assert_series_equal(result, expected)
