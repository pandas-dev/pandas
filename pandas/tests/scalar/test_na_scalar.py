import numpy as np
import pytest

from pandas._libs.missing import NA

from pandas.core.dtypes.common import is_scalar

import pandas as pd
import pandas.util.testing as tm


def test_singleton():
    assert NA is NA
    new_NA = type(NA)()
    assert new_NA is NA


def test_repr():
    assert repr(NA) == "NA"
    assert str(NA) == "NA"


def test_truthiness():
    with pytest.raises(TypeError):
        bool(NA)

    with pytest.raises(TypeError):
        not NA


def test_hashable():
    assert hash(NA) == hash(NA)
    d = {NA: "test"}
    assert d[NA] == "test"


def test_arithmetic_ops(all_arithmetic_functions):
    op = all_arithmetic_functions

    for other in [NA, 1, 1.0, "a", np.int64(1), np.nan]:
        if op.__name__ in ("pow", "rpow", "rmod") and isinstance(other, str):
            continue
        if op.__name__ in ("divmod", "rdivmod"):
            assert op(NA, other) is (NA, NA)
        else:
            if op.__name__ == "rpow":
                # avoid special case
                other += 1
            assert op(NA, other) is NA


def test_comparison_ops():

    for other in [NA, 1, 1.0, "a", np.int64(1), np.nan, np.bool_(True)]:
        assert (NA == other) is NA
        assert (NA != other) is NA
        assert (NA > other) is NA
        assert (NA >= other) is NA
        assert (NA < other) is NA
        assert (NA <= other) is NA
        assert (other == NA) is NA
        assert (other != NA) is NA
        assert (other > NA) is NA
        assert (other >= NA) is NA
        assert (other < NA) is NA
        assert (other <= NA) is NA


@pytest.mark.parametrize(
    "value",
    [
        0,
        0.0,
        -0,
        -0.0,
        False,
        np.bool_(False),
        np.int_(0),
        np.float_(0),
        np.int_(-0),
        np.float_(-0),
    ],
)
def test_pow_special(value):
    result = pd.NA ** value
    assert isinstance(result, type(value))
    assert result == 1


@pytest.mark.parametrize(
    "value",
    [
        1,
        1.0,
        -1,
        -1.0,
        True,
        np.bool_(True),
        np.int_(1),
        np.float_(1),
        np.int_(-1),
        np.float_(-1),
    ],
)
def test_rpow_special(value):
    result = value ** pd.NA
    assert result == value
    if not isinstance(value, (np.float_, np.bool_, np.int_)):
        assert isinstance(result, type(value))


def test_unary_ops():
    assert +NA is NA
    assert -NA is NA
    assert abs(NA) is NA
    assert ~NA is NA


def test_logical_and():

    assert NA & True is NA
    assert True & NA is NA
    assert NA & False is False
    assert False & NA is False
    assert NA & NA is NA

    with pytest.raises(TypeError):
        NA & 5


def test_logical_or():

    assert NA | True is True
    assert True | NA is True
    assert NA | False is NA
    assert False | NA is NA
    assert NA | NA is NA

    with pytest.raises(TypeError):
        NA | 5


def test_logical_xor():

    assert NA ^ True is NA
    assert True ^ NA is NA
    assert NA ^ False is NA
    assert False ^ NA is NA
    assert NA ^ NA is NA

    with pytest.raises(TypeError):
        NA ^ 5


def test_logical_not():
    assert ~NA is NA


@pytest.mark.parametrize(
    "shape",
    [
        # (),  #  TODO: this fails, since we return a scalar. What to do?
        (3,),
        (3, 3),
        (1, 2, 3),
    ],
)
def test_arithmetic_ndarray(shape, all_arithmetic_functions):
    op = all_arithmetic_functions
    a = np.zeros(shape)
    if op.__name__ == "pow":
        a += 5
    result = op(pd.NA, a)
    expected = np.full(a.shape, pd.NA, dtype=object)
    tm.assert_numpy_array_equal(result, expected)


# TODO: test pow special with ndarray


def test_is_scalar():
    assert is_scalar(NA) is True


def test_isna():
    assert pd.isna(NA) is True
    assert pd.notna(NA) is False


def test_series_isna():
    s = pd.Series([1, NA], dtype=object)
    expected = pd.Series([False, True])
    tm.assert_series_equal(s.isna(), expected)


def test_ufunc():
    assert np.log(pd.NA) is pd.NA
    assert np.add(pd.NA, 1) is pd.NA
    result = np.divmod(pd.NA, 1)
    assert result[0] is pd.NA and result[1] is pd.NA

    result = np.frexp(pd.NA)
    assert result[0] is pd.NA and result[1] is pd.NA


def test_ufunc_raises():
    with pytest.raises(ValueError, match="ufunc method 'at'"):
        np.log.at(pd.NA, 0)


def test_binary_input_not_dunder():
    a = np.array([1, 2, 3])
    expected = np.array([pd.NA, pd.NA, pd.NA], dtype=object)
    result = np.logaddexp(a, pd.NA)
    tm.assert_numpy_array_equal(result, expected)

    result = np.logaddexp(pd.NA, a)
    tm.assert_numpy_array_equal(result, expected)
