import numpy as np

from pandas import compat
from pandas.core.common import isnull

cdef bint isiterable(obj):
    return hasattr(obj, '__iter__')

cdef bint decimal_almost_equal(double desired, double actual, int decimal):
    # Code from
    # http://docs.scipy.org/doc/numpy/reference/generated
    # /numpy.testing.assert_almost_equal.html
    return abs(desired - actual) < (0.5 * 10.0 ** -decimal)

cpdef assert_dict_equal(a, b, bint compare_keys=True):
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

    if isinstance(a, compat.string_types):
        assert a == b, "%r != %r" % (a, b)
        return True

    if isiterable(a):
        assert isiterable(b), "First object is iterable, second isn't"
        na, nb = len(a), len(b)
        assert na == nb, "%s != %s" % (na, nb)
        if (isinstance(a, np.ndarray) and
                isinstance(b, np.ndarray) and
                np.array_equal(a, b)):
            return True
        else:
            for i in xrange(na):
                assert_almost_equal(a[i], b[i], check_less_precise)
        return True

    if isnull(a):
        assert isnull(b), "First object is null, second isn't"
        return True

    if isinstance(a, (bool, float, int, np.float32)):
        decimal = 5

        # deal with differing dtypes
        if check_less_precise:
            dtype_a = np.dtype(type(a))
            dtype_b = np.dtype(type(b))
            if dtype_a.kind == 'f' and dtype_b == 'f':
                if dtype_a.itemsize <= 4 and dtype_b.itemsize <= 4:
                    decimal = 3

        if np.isinf(a):
            assert np.isinf(b), "First object is inf, second isn't"
        else:
            fa, fb = a, b

            # case for zero
            if abs(fa) < 1e-5:
                if not decimal_almost_equal(fa, fb, decimal):
                    assert False, (
                        '(very low values) expected %.5f but got %.5f' % (b, a)
                    )
            else:
                if not decimal_almost_equal(1, fb / fa, decimal):
                    assert False, 'expected %.5f but got %.5f' % (b, a)

    else:
        assert a == b, "%s != %s" % (a, b)

    return True
