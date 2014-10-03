"""
Tests related to deprecation warnings. Also a convenient place
to document how deprecations should eventually be turned into errors.

"""
from __future__ import division, absolute_import, print_function

import sys
import operator
import warnings
from nose.plugins.skip import SkipTest

import numpy as np
from numpy.testing import (dec, run_module_suite, assert_raises,
                           assert_warns, assert_array_equal, assert_)


class _DeprecationTestCase(object):
    # Just as warning: warnings uses re.match, so the start of this message
    # must match.
    message = ''

    def setUp(self):
        self.warn_ctx = warnings.catch_warnings(record=True)
        self.log = self.warn_ctx.__enter__()

        # Do *not* ignore other DeprecationWarnings. Ignoring warnings
        # can give very confusing results because of
        # http://bugs.python.org/issue4180 and it is probably simplest to
        # try to keep the tests cleanly giving only the right warning type.
        # (While checking them set to "error" those are ignored anyway)
        # We still have them show up, because otherwise they would be raised
        warnings.filterwarnings("always", category=DeprecationWarning)
        warnings.filterwarnings("always", message=self.message,
                                    category=DeprecationWarning)


    def tearDown(self):
        self.warn_ctx.__exit__()


    def assert_deprecated(self, function, num=1, ignore_others=False,
                        function_fails=False,
                        exceptions=(DeprecationWarning,), args=(), kwargs={}):
        """Test if DeprecationWarnings are given and raised.

        This first checks if the function when called gives `num`
        DeprecationWarnings, after that it tries to raise these
        DeprecationWarnings and compares them with `exceptions`.
        The exceptions can be different for cases where this code path
        is simply not anticipated and the exception is replaced.

        Parameters
        ----------
        f : callable
            The function to test
        num : int
            Number of DeprecationWarnings to expect. This should normally be 1.
        ignore_other : bool
            Whether warnings of the wrong type should be ignored (note that
            the message is not checked)
        function_fails : bool
            If the function would normally fail, setting this will check for
            warnings inside a try/except block.
        exceptions : Exception or tuple of Exceptions
            Exception to expect when turning the warnings into an error.
            The default checks for DeprecationWarnings. If exceptions is
            empty the function is expected to run successfull.
        args : tuple
            Arguments for `f`
        kwargs : dict
            Keyword arguments for `f`
        """
        # reset the log
        self.log[:] = []

        try:
            function(*args, **kwargs)
        except (Exception if function_fails else tuple()):
            pass
        # just in case, clear the registry
        num_found = 0
        for warning in self.log:
            if warning.category is DeprecationWarning:
                num_found += 1
            elif not ignore_others:
                raise AssertionError("expected DeprecationWarning but %s given"
                                                            % warning.category)
        if num is not None and num_found != num:
            raise AssertionError("%i warnings found but %i expected"
                                                        % (len(self.log), num))

        with warnings.catch_warnings():
            warnings.filterwarnings("error", message=self.message,
                                        category=DeprecationWarning)

            try:
                function(*args, **kwargs)
                if exceptions != tuple():
                    raise AssertionError("No error raised during function call")
            except exceptions:
                if exceptions == tuple():
                    raise AssertionError("Error raised during function call")


    def assert_not_deprecated(self, function, args=(), kwargs={}):
        """Test if DeprecationWarnings are given and raised.

        This is just a shorthand for:

        self.assert_deprecated(function, num=0, ignore_others=True,
                        exceptions=tuple(), args=args, kwargs=kwargs)
        """
        self.assert_deprecated(function, num=0, ignore_others=True,
                        exceptions=tuple(), args=args, kwargs=kwargs)


class TestFloatNonIntegerArgumentDeprecation(_DeprecationTestCase):
    """
    These test that ``DeprecationWarning`` is given when you try to use
    non-integers as arguments to for indexing and slicing e.g. ``a[0.0:5]``
    and ``a[0.5]``, or other functions like ``array.reshape(1., -1)``.

    After deprecation, changes need to be done inside conversion_utils.c
    in PyArray_PyIntAsIntp and possibly PyArray_IntpConverter.
    In iterators.c the function slice_GetIndices could be removed in favor
    of its python equivalent and in mapping.c the function _tuple_of_integers
    can be simplified (if ``np.array([1]).__index__()`` is also deprecated).

    As for the deprecation time-frame: via Ralf Gommers,

    "Hard to put that as a version number, since we don't know if the
    version after 1.8 will be 6 months or 2 years after. I'd say 2
    years is reasonable."

    I interpret this to mean 2 years after the 1.8 release. Possibly
    giving a PendingDeprecationWarning before that (which is visible
    by default)

    """
    message = "using a non-integer number instead of an integer " \
              "will result in an error in the future"

    def test_indexing(self):
        a = np.array([[[5]]])
        def assert_deprecated(*args, **kwargs):
            self.assert_deprecated(*args, exceptions=(IndexError,), **kwargs)

        assert_deprecated(lambda: a[0.0])
        assert_deprecated(lambda: a[0, 0.0])
        assert_deprecated(lambda: a[0.0, 0])
        assert_deprecated(lambda: a[0.0,:])
        assert_deprecated(lambda: a[:, 0.0])
        assert_deprecated(lambda: a[:, 0.0,:])
        assert_deprecated(lambda: a[0.0,:,:])
        assert_deprecated(lambda: a[0, 0, 0.0])
        assert_deprecated(lambda: a[0.0, 0, 0])
        assert_deprecated(lambda: a[0, 0.0, 0])
        assert_deprecated(lambda: a[-1.4])
        assert_deprecated(lambda: a[0, -1.4])
        assert_deprecated(lambda: a[-1.4, 0])
        assert_deprecated(lambda: a[-1.4,:])
        assert_deprecated(lambda: a[:, -1.4])
        assert_deprecated(lambda: a[:, -1.4,:])
        assert_deprecated(lambda: a[-1.4,:,:])
        assert_deprecated(lambda: a[0, 0, -1.4])
        assert_deprecated(lambda: a[-1.4, 0, 0])
        assert_deprecated(lambda: a[0, -1.4, 0])

        # Test that the slice parameter deprecation warning doesn't mask
        # the scalar index warning.
        assert_deprecated(lambda: a[0.0:, 0.0], num=2)
        assert_deprecated(lambda: a[0.0:, 0.0,:], num=2)


    def test_valid_indexing(self):
        a = np.array([[[5]]])
        assert_not_deprecated = self.assert_not_deprecated

        assert_not_deprecated(lambda: a[np.array([0])])
        assert_not_deprecated(lambda: a[[0, 0]])
        assert_not_deprecated(lambda: a[:, [0, 0]])
        assert_not_deprecated(lambda: a[:, 0,:])
        assert_not_deprecated(lambda: a[:,:,:])


    def test_slicing(self):
        a = np.array([[5]])
        def assert_deprecated(*args, **kwargs):
            self.assert_deprecated(*args, exceptions=(IndexError,), **kwargs)

        # start as float.
        assert_deprecated(lambda: a[0.0:])
        assert_deprecated(lambda: a[0:, 0.0:2])
        assert_deprecated(lambda: a[0.0::2, :0])
        assert_deprecated(lambda: a[0.0:1:2,:])
        assert_deprecated(lambda: a[:, 0.0:])
        # stop as float.
        assert_deprecated(lambda: a[:0.0])
        assert_deprecated(lambda: a[:0, 1:2.0])
        assert_deprecated(lambda: a[:0.0:2, :0])
        assert_deprecated(lambda: a[:0.0,:])
        assert_deprecated(lambda: a[:, 0:4.0:2])
        # step as float.
        assert_deprecated(lambda: a[::1.0])
        assert_deprecated(lambda: a[0:, :2:2.0])
        assert_deprecated(lambda: a[1::4.0, :0])
        assert_deprecated(lambda: a[::5.0,:])
        assert_deprecated(lambda: a[:, 0:4:2.0])
        # mixed.
        assert_deprecated(lambda: a[1.0:2:2.0], num=2)
        assert_deprecated(lambda: a[1.0::2.0], num=2)
        assert_deprecated(lambda: a[0:, :2.0:2.0], num=2)
        assert_deprecated(lambda: a[1.0:1:4.0, :0], num=2)
        assert_deprecated(lambda: a[1.0:5.0:5.0,:], num=3)
        assert_deprecated(lambda: a[:, 0.4:4.0:2.0], num=3)
        # should still get the DeprecationWarning if step = 0.
        assert_deprecated(lambda: a[::0.0], function_fails=True)


    def test_valid_slicing(self):
        a = np.array([[[5]]])
        assert_not_deprecated = self.assert_not_deprecated

        assert_not_deprecated(lambda: a[::])
        assert_not_deprecated(lambda: a[0:])
        assert_not_deprecated(lambda: a[:2])
        assert_not_deprecated(lambda: a[0:2])
        assert_not_deprecated(lambda: a[::2])
        assert_not_deprecated(lambda: a[1::2])
        assert_not_deprecated(lambda: a[:2:2])
        assert_not_deprecated(lambda: a[1:2:2])


    def test_non_integer_argument_deprecations(self):
        a = np.array([[5]])

        self.assert_deprecated(np.reshape, args=(a, (1., 1., -1)), num=2)
        self.assert_deprecated(np.reshape, args=(a, (np.array(1.), -1)))
        self.assert_deprecated(np.take, args=(a, [0], 1.))
        self.assert_deprecated(np.take, args=(a, [0], np.float64(1.)))


    def test_non_integer_sequence_multiplication(self):
        # Numpy scalar sequence multiply should not work with non-integers
        def mult(a, b):
            return a * b
        self.assert_deprecated(mult, args=([1], np.float_(3)))
        self.assert_not_deprecated(mult, args=([1], np.int_(3)))


    def test_reduce_axis_float_index(self):
        d = np.zeros((3,3,3))
        self.assert_deprecated(np.min, args=(d, 0.5))
        self.assert_deprecated(np.min, num=1, args=(d, (0.5, 1)))
        self.assert_deprecated(np.min, num=1, args=(d, (1, 2.2)))
        self.assert_deprecated(np.min, num=2, args=(d, (.2, 1.2)))


class TestBooleanArgumentDeprecation(_DeprecationTestCase):
    """This tests that using a boolean as integer argument/indexing is
    deprecated.

    This should be kept in sync with TestFloatNonIntegerArgumentDeprecation
    and like it is handled in PyArray_PyIntAsIntp.
    """
    message = "using a boolean instead of an integer " \
              "will result in an error in the future"

    def test_bool_as_int_argument(self):
        a = np.array([[[1]]])

        self.assert_deprecated(np.reshape, args=(a, (True, -1)))
        self.assert_deprecated(np.reshape, args=(a, (np.bool_(True), -1)))
        # Note that operator.index(np.array(True)) does not work, a boolean
        # array is thus also deprecated, but not with the same message:
        assert_raises(TypeError, operator.index, np.array(True))
        self.assert_deprecated(np.take, args=(a, [0], False))
        self.assert_deprecated(lambda: a[False:True:True], exceptions=IndexError, num=3)
        self.assert_deprecated(lambda: a[False, 0], exceptions=IndexError)
        self.assert_deprecated(lambda: a[False, 0, 0], exceptions=IndexError)


class TestArrayToIndexDeprecation(_DeprecationTestCase):
    """This tests that creating an an index from an array is deprecated
    if the array is not 0d.

    This can probably be deprecated somewhat faster then the integer
    deprecations. The deprecation period started with NumPy 1.8.
    For deprecation this needs changing of array_index in number.c
    """
    message = "converting an array with ndim \> 0 to an index will result " \
              "in an error in the future"

    def test_array_to_index_deprecation(self):
        # This drops into the non-integer deprecation, which is ignored here,
        # so no exception is expected. The raising is effectively tested above.
        a = np.array([[[1]]])

        self.assert_deprecated(operator.index, args=(np.array([1]),))
        self.assert_deprecated(np.reshape, args=(a, (a, -1)), exceptions=())
        self.assert_deprecated(np.take, args=(a, [0], a), exceptions=())
        # Check slicing. Normal indexing checks arrays specifically.
        self.assert_deprecated(lambda: a[a:a:a], exceptions=(), num=3)


class TestNonIntegerArrayLike(_DeprecationTestCase):
    """Tests that array likes, i.e. lists give a deprecation warning
    when they cannot be safely cast to an integer.
    """
    message = "non integer \(and non boolean\) array-likes will not be " \
              "accepted as indices in the future"

    def test_basic(self):
        a = np.arange(10)
        self.assert_deprecated(a.__getitem__, args=([0.5, 1.5],),
                               exceptions=IndexError)
        self.assert_deprecated(a.__getitem__, args=((['1', '2'],),),
                               exceptions=IndexError)

        self.assert_not_deprecated(a.__getitem__, ([],))


    def test_boolean_futurewarning(self):
        a = np.arange(10)
        with warnings.catch_warnings():
            warnings.filterwarnings('always')
            assert_warns(FutureWarning, a.__getitem__, [True])
            # Unfortunatly, the deprecation warning takes precedence:
            #assert_warns(FutureWarning, a.__getitem__, True)

        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            assert_raises(FutureWarning, a.__getitem__, [True])
            #assert_raises(FutureWarning, a.__getitem__, True)


class TestMultipleEllipsisDeprecation(_DeprecationTestCase):
    message = "an index can only have a single Ellipsis \(`...`\); replace " \
              "all but one with slices \(`:`\)."

    def test_basic(self):
        a = np.arange(10)
        self.assert_deprecated(a.__getitem__, args=((Ellipsis, Ellipsis),))

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', '', DeprecationWarning)
            # Just check that this works:
            b = a[...,...]
            assert_array_equal(a, b)
            assert_raises(IndexError, a.__getitem__, ((Ellipsis, ) * 3,))


class TestBooleanSubtractDeprecations(_DeprecationTestCase):
    """Test deprecation of boolean `-`. While + and * are well
    defined, - is not and even a corrected form seems to have
    no real uses.

    The deprecation process was started in NumPy 1.9.
    """
    message = r"numpy boolean .* \(the .* `-` operator\) is deprecated, " \
              "use the bitwise"

    def test_operator_deprecation(self):
        array = np.array([True])
        generic = np.bool_(True)

        # Minus operator/subtract ufunc:
        self.assert_deprecated(operator.sub, args=(array, array))
        self.assert_deprecated(operator.sub, args=(generic, generic))

        # Unary minus/negative ufunc:
        self.assert_deprecated(operator.neg, args=(array,))
        self.assert_deprecated(operator.neg, args=(generic,))


class TestRankDeprecation(_DeprecationTestCase):
    """Test that np.rank is deprecated. The function should simply be
    removed. The VisibleDeprecationWarning may become unnecessary.
    """
    def test(self):
        a = np.arange(10)
        assert_warns(np.VisibleDeprecationWarning, np.rank, a)


class TestComparisonDepreactions(_DeprecationTestCase):
    """This tests the deprecation, for non-elementwise comparison logic.
    This used to mean that when an error occured during element-wise comparison
    (i.e. broadcasting) NotImplemented was returned, but also in the comparison
    itself, False was given instead of the error.

    Also test FutureWarning for the None comparison.
    """

    message = "elementwise comparison failed; " \
              "this will raise the error in the future."

    def test_normal_types(self):
        for op in (operator.eq, operator.ne):
            # Broadcasting errors:
            self.assert_deprecated(op, args=(np.zeros(3), []))
            a = np.zeros(3, dtype='i,i')
            # (warning is issued a couple of times here)
            self.assert_deprecated(op, args=(a, a[:-1]), num=None)

            # Element comparison error (numpy array can't be compared).
            a = np.array([1, np.array([1,2,3])], dtype=object)
            b = np.array([1, np.array([1,2,3])], dtype=object)
            self.assert_deprecated(op, args=(a, b), num=None)


    def test_string(self):
        # For two string arrays, strings always raised the broadcasting error:
        a = np.array(['a', 'b'])
        b = np.array(['a', 'b', 'c'])
        assert_raises(ValueError, lambda x, y: x == y, a, b)

        # The empty list is not cast to string, this is only to document
        # that fact (it likely should be changed). This means that the
        # following works (and returns False) due to dtype mismatch:
        a == []


    def test_none_comparison(self):
        # Test comparison of None, which should result in elementwise
        # comparison in the future. [1, 2] == None should be [False, False].
        with warnings.catch_warnings():
            warnings.filterwarnings('always', '', FutureWarning)
            assert_warns(FutureWarning, operator.eq, np.arange(3), None)
            assert_warns(FutureWarning, operator.ne, np.arange(3), None)

        with warnings.catch_warnings():
            warnings.filterwarnings('error', '', FutureWarning)
            assert_raises(FutureWarning, operator.eq, np.arange(3), None)
            assert_raises(FutureWarning, operator.ne, np.arange(3), None)

    def test_scalar_none_comparison(self):
        # Scalars should still just return false and not give a warnings.
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', '', FutureWarning)
            assert_(not np.float32(1) == None)
            assert_(not np.str_('test') == None)
            # This is dubious (see below):
            assert_(not np.datetime64('NaT') == None)

            assert_(np.float32(1) != None)
            assert_(np.str_('test') != None)
            # This is dubious (see below):
            assert_(np.datetime64('NaT') != None)
        assert_(len(w) == 0)

        # For documentaiton purpose, this is why the datetime is dubious.
        # At the time of deprecation this was no behaviour change, but
        # it has to be considered when the deprecations is done.
        assert_(np.equal(np.datetime64('NaT'), None))


class TestIdentityComparisonDepreactions(_DeprecationTestCase):
    """This tests the equal and not_equal object ufuncs identity check
    deprecation. This was due to the usage of PyObject_RichCompareBool.

    This tests that for example for `a = np.array([np.nan], dtype=object)`
    `a == a` it is warned that False and not `np.nan is np.nan` is returned.

    Should be kept in sync with TestComparisonDepreactions and new tests
    added when the deprecation is over. Requires only removing of @identity@
    (and blocks) from the ufunc loops.c.src of the OBJECT comparisons.
    """

    message = "numpy .* will not check object identity in the future."

    def test_identity_equality_mismatch(self):
        a = np.array([np.nan], dtype=object)

        with warnings.catch_warnings():
            warnings.filterwarnings('always', '', FutureWarning)
            assert_warns(FutureWarning, np.equal, a, a)
            assert_warns(FutureWarning, np.not_equal, a, a)

        with warnings.catch_warnings():
            warnings.filterwarnings('error', '', FutureWarning)
            assert_raises(FutureWarning, np.equal, a, a)
            assert_raises(FutureWarning, np.not_equal, a, a)
            # And the other do not warn:
            with np.errstate(invalid='ignore'):
                np.less(a, a)
                np.greater(a, a)
                np.less_equal(a, a)
                np.greater_equal(a, a)


    def test_comparison_error(self):
        class FunkyType(object):
            def __eq__(self, other):
                raise TypeError("I won't compare")
            def __ne__(self, other):
                raise TypeError("I won't compare")

        a = np.array([FunkyType()])
        self.assert_deprecated(np.equal, args=(a, a))
        self.assert_deprecated(np.not_equal, args=(a, a))


    def test_bool_error(self):
        # The comparison result cannot be interpreted as a bool
        a = np.array([np.array([1, 2, 3]), None], dtype=object)
        self.assert_deprecated(np.equal, args=(a, a))
        self.assert_deprecated(np.not_equal, args=(a, a))


if __name__ == "__main__":
    run_module_suite()
