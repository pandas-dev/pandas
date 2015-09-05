import numpy as np

from pandas import compat
from pandas.core.common import isnull, array_equivalent

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

cpdef assert_almost_equal(a, b, bint check_less_precise=False,
                          obj=None, lobj=None, robj=None):
    """Check that left and right objects are almost equal.

    Parameters
    ----------
    a : object
    b : object
    check_less_precise : bool, default False
        Specify comparison precision.
        5 digits (False) or 3 digits (True) after decimal points are compared.
    obj : str, default None
        Specify object name being compared, internally used to show appropriate
        assertion message
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
        bint is_unequal = False

    if lobj is None:
        lobj = a
    if robj is None:
        robj = b

    if isinstance(a, dict) or isinstance(b, dict):
        return assert_dict_equal(a, b)

    if (isinstance(a, compat.string_types) or
            isinstance(b, compat.string_types)):
        assert a == b, "%r != %r" % (a, b)
        return True

    if isiterable(a):

        if not isiterable(b):
            from pandas.util.testing import raise_assert_detail
            if obj is None:
                obj = 'Iterable'
            msg = "First object is iterable, second isn't"
            raise_assert_detail(obj, msg, a, b)

        assert has_length(a) and has_length(b), (
            "Can't compare objects without length, one or both is invalid: "
            "(%r, %r)" % (a, b)
        )

        if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
            if obj is None:
                obj = 'numpy array'
            na, nb = a.size, b.size
            if a.shape != b.shape:
                from pandas.util.testing import raise_assert_detail
                raise_assert_detail(obj, '{0} shapes are different'.format(obj),
                                    a.shape, b.shape)
            try:
                if array_equivalent(a, b, strict_nan=True):
                    return True
            except:
                pass
        else:
            if obj is None:
                obj = 'Iterable'
            na, nb = len(a), len(b)

        if na != nb:
            from pandas.util.testing import raise_assert_detail
            raise_assert_detail(obj, '{0} length are different'.format(obj),
                                na, nb)

        for i in xrange(len(a)):
            try:
                assert_almost_equal(a[i], b[i], check_less_precise)
            except AssertionError:
                is_unequal = True
                diff += 1

        if is_unequal:
            from pandas.util.testing import raise_assert_detail
            msg = '{0} values are different ({1} %)'.format(obj, np.round(diff * 100.0 / na, 5))
            raise_assert_detail(obj, msg, lobj, robj)

        return True

    elif isiterable(b):
        from pandas.util.testing import raise_assert_detail
        if obj is None:
            obj = 'Iterable'
        msg = "Second object is iterable, first isn't"
        raise_assert_detail(obj, msg, a, b)

    if isnull(a):
        assert isnull(b), (
            "First object is null, second isn't: %r != %r" % (a, b)
        )
        return True
    elif isnull(b):
        assert isnull(a), (
            "First object is not null, second is null: %r != %r" % (a, b)
        )
        return True

    if is_comparable_as_number(a):
        assert is_comparable_as_number(b), (
            "First object is numeric, second is not: %r != %r" % (a, b)
        )

        decimal = 5

        # deal with differing dtypes
        if check_less_precise:
            decimal = 3

        if np.isinf(a):
            assert np.isinf(b), "First object is inf, second isn't"
            if np.isposinf(a):
                assert np.isposinf(b), "First object is positive inf, second is negative inf"
            else:
                assert np.isneginf(b), "First object is negative inf, second is positive inf"
        else:
            fa, fb = a, b

            # case for zero
            if abs(fa) < 1e-5:
                if not decimal_almost_equal(fa, fb, decimal):
                    assert False, (
                        '(very low values) expected %.5f but got %.5f, with decimal %d' % (fb, fa, decimal)
                    )
            else:
                if not decimal_almost_equal(1, fb / fa, decimal):
                    assert False, 'expected %.5f but got %.5f, with decimal %d' % (fb, fa, decimal)

    else:
        assert a == b, "%r != %r" % (a, b)

    return True
