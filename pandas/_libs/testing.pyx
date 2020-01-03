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
    """
    Check that left and right objects are almost equal.

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
        assert a == b, f"{a} != {b}"
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
            from pandas._testing import assert_class_equal
            # classes can't be the same, to raise error
            assert_class_equal(a, b, obj=obj)

        assert has_length(a) and has_length(b), ("Can't compare objects without "
                                                 "length, one or both is invalid: "
                                                 f"({a}, {b})")

        if a_is_ndarray and b_is_ndarray:
            na, nb = a.size, b.size
            if a.shape != b.shape:
                from pandas._testing import raise_assert_detail
                raise_assert_detail(
                    obj, f'{obj} shapes are different', a.shape, b.shape)

            if check_dtype and not is_dtype_equal(a.dtype, b.dtype):
                from pandas._testing import assert_attr_equal
                assert_attr_equal('dtype', a, b, obj=obj)

            if array_equivalent(a, b, strict_nan=True):
                return True

        else:
            na, nb = len(a), len(b)

        if na != nb:
            from pandas._testing import raise_assert_detail

            # if we have a small diff set, print it
            if abs(na - nb) < 10:
                r = list(set(a) ^ set(b))
            else:
                r = None

            raise_assert_detail(obj, f"{obj} length are different", na, nb, r)

        for i in range(len(a)):
            try:
                assert_almost_equal(a[i], b[i],
                                    check_less_precise=check_less_precise)
            except AssertionError:
                is_unequal = True
                diff += 1

        if is_unequal:
            from pandas._testing import raise_assert_detail
            msg = (f"{obj} values are different "
                   f"({np.round(diff * 100.0 / na, 5)} %)")
            raise_assert_detail(obj, msg, lobj, robj)

        return True

    elif isiterable(b):
        from pandas._testing import assert_class_equal
        # classes can't be the same, to raise error
        assert_class_equal(a, b, obj=obj)

    if isna(a) and isna(b):
        # TODO: Should require same-dtype NA?
        # nan / None comparison
        return True

    if a == b:
        # object comparison
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
                assert False, (f'(very low values) expected {fb:.5f} '
                               f'but got {fa:.5f}, with decimal {decimal}')
        else:
            if not decimal_almost_equal(1, fb / fa, decimal):
                assert False, (f'expected {fb:.5f} but got {fa:.5f}, '
                               f'with decimal {decimal}')
        return True

    raise AssertionError(f"{a} != {b}")
