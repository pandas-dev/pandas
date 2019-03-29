import numpy as np

from pandas.core.dtypes.missing import isna, array_equivalent
from pandas.core.dtypes.common import is_dtype_equal

cdef NUMERIC_TYPES = (
    bool,
    int,
    float,
    np.bool,
    np.int8,
    np.int16,
    np.int32,
    np.int64,
    np.uint8,
    np.uint16,
    np.uint32,
    np.uint64,
    np.float16,
    np.float32,
    np.float64,
)


cdef bint is_comparable_as_number(obj):
    return isinstance(obj, NUMERIC_TYPES)


cdef bint isiterable(obj):
    return hasattr(obj, '__iter__')


cdef bint has_length(obj):
    return hasattr(obj, '__len__')


cdef bint is_dictlike(obj):
    return hasattr(obj, 'keys') and hasattr(obj, '__getitem__')


cdef bint decimal_almost_equal(double desired, double actual, int decimal):
    # Code from
    # http://docs.scipy.org/doc/numpy/reference/generated
    # /numpy.testing.assert_almost_equal.html
    return abs(desired - actual) < (0.5 * 10.0 ** -decimal)


cpdef assert_dict_equal(a, b, bint compare_keys=True):
    assert is_dictlike(a) and is_dictlike(b), (
        "Cannot compare dict objects, one or both is not dict-like"
    )

    a_keys = frozenset(a.keys())
    b_keys = frozenset(b.keys())

    if compare_keys:
        assert a_keys == b_keys

    for k in a_keys:
        assert_almost_equal(a[k], b[k])

    return True


cpdef assert_almost_equal(a, b,
                          check_less_precise=False,
                          bint check_dtype=True,
                          obj=None, lobj=None, robj=None):
    """Check that left and right objects are almost equal.

    Parameters
    ----------
    a : object
    b : object
    check_less_precise : bool or int, default False
        Specify comparison precision.
        5 digits (False) or 3 digits (True) after decimal points are
        compared. If an integer, then this will be the number of decimal
        points to compare
    check_dtype: bool, default True
        check dtype if both a and b are np.ndarray
    obj : str, default None
        Specify object name being compared, internally used to show
        appropriate assertion message
    lobj : str, default None
        Specify left object name being compared, internally used to show
        appropriate assertion message
    robj : str, default None
        Specify right object name being compared, internally used to show
        appropriate assertion message
    """

    cdef:
        int decimal
        double diff = 0.0
        Py_ssize_t i, na, nb
        double fa, fb
        bint is_unequal = False, a_is_ndarray, b_is_ndarray

    if lobj is None:
        lobj = a
    if robj is None:
        robj = b

    assert isinstance(check_less_precise, (int, bool))

    if isinstance(a, dict) or isinstance(b, dict):
        return assert_dict_equal(a, b)

    if isinstance(a, str) or isinstance(b, str):
        assert a == b, "%r != %r" % (a, b)
        return True

    a_is_ndarray = isinstance(a, np.ndarray)
    b_is_ndarray = isinstance(b, np.ndarray)

    if obj is None:
        if a_is_ndarray or b_is_ndarray:
            obj = 'numpy array'
        else:
            obj = 'Iterable'

    if isiterable(a):

        if not isiterable(b):
            from pandas.util.testing import assert_class_equal
            # classes can't be the same, to raise error
            assert_class_equal(a, b, obj=obj)

        assert has_length(a) and has_length(b), (
            "Can't compare objects without length, one or both is invalid: "
            "(%r, %r)" % (a, b))

        if a_is_ndarray and b_is_ndarray:
            na, nb = a.size, b.size
            if a.shape != b.shape:
                from pandas.util.testing import raise_assert_detail
                raise_assert_detail(
                    obj, '{0} shapes are different'.format(obj),
                    a.shape, b.shape)

            if check_dtype and not is_dtype_equal(a, b):
                from pandas.util.testing import assert_attr_equal
                assert_attr_equal('dtype', a, b, obj=obj)

            try:
                if array_equivalent(a, b, strict_nan=True):
                    return True
            except:
                pass
        else:
            na, nb = len(a), len(b)

        if na != nb:
            from pandas.util.testing import raise_assert_detail

            # if we have a small diff set, print it
            if abs(na - nb) < 10:
                r = list(set(a) ^ set(b))
            else:
                r = None

            raise_assert_detail(obj, '{0} length are different'.format(obj),
                                na, nb, r)

        for i in xrange(len(a)):
            try:
                assert_almost_equal(a[i], b[i],
                                    check_less_precise=check_less_precise)
            except AssertionError:
                is_unequal = True
                diff += 1

        if is_unequal:
            from pandas.util.testing import raise_assert_detail
            msg = '{0} values are different ({1} %)'.format(
                obj, np.round(diff * 100.0 / na, 5))
            raise_assert_detail(obj, msg, lobj, robj)

        return True

    elif isiterable(b):
        from pandas.util.testing import assert_class_equal
        # classes can't be the same, to raise error
        assert_class_equal(a, b, obj=obj)

    if a == b:
        # object comparison
        return True
    if isna(a) and isna(b):
        # nan / None comparison
        return True
    if is_comparable_as_number(a) and is_comparable_as_number(b):
        if array_equivalent(a, b, strict_nan=True):
            # inf comparison
            return True

        if check_less_precise is True:
            decimal = 3
        elif check_less_precise is False:
            decimal = 5
        else:
            decimal = check_less_precise

        fa, fb = a, b

        # case for zero
        if abs(fa) < 1e-5:
            if not decimal_almost_equal(fa, fb, decimal):
                assert False, ('(very low values) expected %.5f but '
                               'got %.5f, with decimal %d' % (fb, fa, decimal))
        else:
            if not decimal_almost_equal(1, fb / fa, decimal):
                assert False, ('expected %.5f but got %.5f, '
                               'with decimal %d' % (fb, fa, decimal))
        return True

    raise AssertionError("{0} != {1}".format(a, b))
