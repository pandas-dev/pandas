import string

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm

UNARY_UFUNCS = [np.positive, np.floor, np.exp]
BINARY_UFUNCS = [
    np.add,  # dunder op
    np.logaddexp,
]
SPARSE = [
    pytest.param(True,
                 marks=pytest.mark.xfail(reason="Series.__array_ufunc__")),
    False,
]
SPARSE_IDS = ['sparse', 'dense']
SHUFFLE = [
    pytest.param(True, marks=pytest.mark.xfail(reason="GH-26945",
                                               strict=False)),
    False
]


@pytest.fixture
def arrays_for_binary_ufunc():
    """
    A pair of random, length-100 integer-dtype arrays, that are mostly 0.
    """
    a1 = np.random.randint(0, 10, 100, dtype='int64')
    a2 = np.random.randint(0, 10, 100, dtype='int64')
    a1[::3] = 0
    a2[::4] = 0
    return a1, a2


@pytest.mark.parametrize("ufunc", UNARY_UFUNCS)
@pytest.mark.parametrize("sparse", SPARSE, ids=SPARSE_IDS)
def test_unary_ufunc(ufunc, sparse):
    # Test that ufunc(Series) == Series(ufunc)
    array = np.random.randint(0, 10, 10, dtype='int64')
    array[::2] = 0
    if sparse:
        array = pd.SparseArray(array, dtype=pd.SparseDtype('int', 0))

    index = list(string.ascii_letters[:10])
    name = "name"
    series = pd.Series(array, index=index, name=name)

    result = ufunc(series)
    expected = pd.Series(ufunc(array), index=index, name=name)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ufunc", BINARY_UFUNCS)
@pytest.mark.parametrize("sparse", SPARSE, ids=SPARSE_IDS)
@pytest.mark.parametrize("flip", [True, False], ids=['flipped', 'straight'])
def test_binary_ufunc_with_array(flip, sparse, ufunc, arrays_for_binary_ufunc):
    # Test that ufunc(Series(a), array) == Series(ufunc(a, b))
    a1, a2 = arrays_for_binary_ufunc
    if sparse:
        a1 = pd.SparseArray(a1, dtype=pd.SparseDtype('int', 0))
        a2 = pd.SparseArray(a2, dtype=pd.SparseDtype('int', 0))

    name = "name"  # op(Series, array) preserves the name.
    series = pd.Series(a1, name=name)
    other = a2

    array_args = (a1, a2)
    series_args = (series, other)            # ufunc(series, array)

    if flip:
        array_args = reversed(array_args)
        series_args = reversed(series_args)  # ufunc(array, series)

    expected = pd.Series(ufunc(*array_args), name=name)
    result = ufunc(*series_args)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ufunc", BINARY_UFUNCS)
@pytest.mark.parametrize("sparse", SPARSE, ids=SPARSE_IDS)
@pytest.mark.parametrize("flip", [
    pytest.param(True, marks=pytest.mark.xfail(reason="Index should defer")),
    False
], ids=['flipped', 'straight'])
def test_binary_ufunc_with_index(flip, sparse, ufunc, arrays_for_binary_ufunc):
    # Test that
    #   * func(Series(a), Series(b)) == Series(ufunc(a, b))
    #   * ufunc(Index, Series) dispatches to Series (returns a Series)
    a1, a2 = arrays_for_binary_ufunc
    if sparse:
        a1 = pd.SparseArray(a1, dtype=pd.SparseDtype('int', 0))
        a2 = pd.SparseArray(a2, dtype=pd.SparseDtype('int', 0))

    name = "name"  # op(Series, array) preserves the name.
    series = pd.Series(a1, name=name)
    other = pd.Index(a2, name=name).astype("int64")

    array_args = (a1, a2)
    series_args = (series, other)            # ufunc(series, array)

    if flip:
        array_args = reversed(array_args)
        series_args = reversed(series_args)  # ufunc(array, series)

    expected = pd.Series(ufunc(*array_args), name=name)
    result = ufunc(*series_args)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ufunc", BINARY_UFUNCS)
@pytest.mark.parametrize("sparse", SPARSE, ids=SPARSE_IDS)
@pytest.mark.parametrize("shuffle", [True, False], ids=['unaligned',
                                                        'aligned'])
@pytest.mark.parametrize("flip", [True, False], ids=['flipped', 'straight'])
def test_binary_ufunc_with_series(flip, shuffle, sparse, ufunc,
                                  arrays_for_binary_ufunc):
    # Test that
    #   * func(Series(a), Series(b)) == Series(ufunc(a, b))
    #   with alignment between the indices

    if flip and shuffle:
        pytest.xfail(reason="Fix with Series.__array_ufunc__")

    a1, a2 = arrays_for_binary_ufunc
    if sparse:
        a1 = pd.SparseArray(a1, dtype=pd.SparseDtype('int', 0))
        a2 = pd.SparseArray(a2, dtype=pd.SparseDtype('int', 0))

    name = "name"  # op(Series, array) preserves the name.
    series = pd.Series(a1, name=name)
    other = pd.Series(a2, name=name)

    idx = np.random.permutation(len(a1))

    if shuffle:
        other = other.take(idx)
        a2 = a2.take(idx)
        # alignment, so the expected index is the first index in the op.
        if flip:
            index = other.align(series)[0].index
        else:
            index = series.align(other)[0].index
    else:
        index = series.index

    array_args = (a1, a2)
    series_args = (series, other)  # ufunc(series, array)

    if flip:
        array_args = tuple(reversed(array_args))
        series_args = tuple(reversed(series_args))  # ufunc(array, series)

    expected = pd.Series(ufunc(*array_args), index=index, name=name)
    result = ufunc(*series_args)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ufunc", BINARY_UFUNCS)
@pytest.mark.parametrize("sparse", SPARSE, ids=SPARSE_IDS)
@pytest.mark.parametrize("flip", [True, False])
def test_binary_ufunc_scalar(ufunc, sparse, flip, arrays_for_binary_ufunc):
    # Test that
    #   * ufunc(Series, scalar) == Series(ufunc(array, scalar))
    #   * ufunc(Series, scalar) == ufunc(scalar, Series)
    array, _ = arrays_for_binary_ufunc
    if sparse:
        array = pd.SparseArray(array)
    other = 2
    series = pd.Series(array, name="name")

    series_args = (series, other)
    array_args = (array, other)

    if flip:
        series_args = tuple(reversed(series_args))
        array_args = tuple(reversed(array_args))

    expected = pd.Series(ufunc(*array_args), name="name")
    result = ufunc(*series_args)

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ufunc", [np.divmod])  # any others?
@pytest.mark.parametrize("sparse", SPARSE, ids=SPARSE_IDS)
@pytest.mark.parametrize("shuffle", SHUFFLE)
@pytest.mark.filterwarnings("ignore:divide by zero:RuntimeWarning")
def test_multiple_ouput_binary_ufuncs(ufunc, sparse, shuffle,
                                      arrays_for_binary_ufunc):
    # Test that
    #  the same conditions from binary_ufunc_scalar apply to
    #  ufuncs with multiple outputs.
    if sparse and ufunc is np.divmod:
        pytest.skip("sparse divmod not implemented.")

    a1, a2 = arrays_for_binary_ufunc

    if sparse:
        a1 = pd.SparseArray(a1, dtype=pd.SparseDtype('int', 0))
        a2 = pd.SparseArray(a2, dtype=pd.SparseDtype('int', 0))

    s1 = pd.Series(a1)
    s2 = pd.Series(a2)

    if shuffle:
        # ensure we align before applying the ufunc
        s2 = s2.sample(frac=1)

    expected = ufunc(a1, a2)
    assert isinstance(expected, tuple)

    result = ufunc(s1, s2)
    assert isinstance(result, tuple)
    tm.assert_series_equal(result[0], pd.Series(expected[0]))
    tm.assert_series_equal(result[1], pd.Series(expected[1]))


@pytest.mark.parametrize("sparse", SPARSE, ids=SPARSE_IDS)
def test_multiple_ouput_ufunc(sparse, arrays_for_binary_ufunc):
    # Test that the same conditions from unary input apply to multi-output
    # ufuncs
    array, _ = arrays_for_binary_ufunc

    if sparse:
        array = pd.SparseArray(array)

    series = pd.Series(array, name="name")
    result = np.modf(series)
    expected = np.modf(array)

    assert isinstance(result, tuple)
    assert isinstance(expected, tuple)

    tm.assert_series_equal(result[0], pd.Series(expected[0], name="name"))
    tm.assert_series_equal(result[1], pd.Series(expected[1], name="name"))


@pytest.mark.parametrize("sparse", SPARSE, ids=SPARSE_IDS)
@pytest.mark.parametrize("ufunc", BINARY_UFUNCS)
@pytest.mark.xfail(reason="Series.__array_ufunc__")
def test_binary_ufunc_drops_series_name(ufunc, sparse,
                                        arrays_for_binary_ufunc):
    # Drop the names when they differ.
    a1, a2 = arrays_for_binary_ufunc
    s1 = pd.Series(a1, name='a')
    s2 = pd.Series(a2, name='b')

    result = ufunc(s1, s2)
    assert result.name is None
