import numpy as np

from pandas import compat
from pandas.core.common import isnull

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

cpdef assert_almost_equal(a, b, bint check_less_precise=False):
    cdef:
        int decimal
        Py_ssize_t i, na, nb
        double fa, fb

    if isinstance(a, dict) or isinstance(b, dict):
        return assert_dict_equal(a, b)

    if (isinstance(a, compat.string_types) or
            isinstance(b, compat.string_types)):
        assert a == b, "%r != %r" % (a, b)
        return True

    if isiterable(a):
        assert isiterable(b), (
            "First object is iterable, second isn't: %r != %r" % (a, b)
        )
        assert has_length(a) and has_length(b), (
            "Can't compare objects without length, one or both is invalid: "
            "(%r, %r)" % (a, b)
        )

        na, nb = len(a), len(b)
        assert na == nb, (
            "Length of two iterators not the same: %r != %r" % (na, nb)
        )
        if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
            try:
                if np.array_equal(a, b):
                    return True
            except:
                pass

        for i in xrange(na):
            assert_almost_equal(a[i], b[i], check_less_precise)

        return True
    elif isiterable(b):
        assert False, (
            "Second object is iterable, first isn't: %r != %r" % (a, b)
        )

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
