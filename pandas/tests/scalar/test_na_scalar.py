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

        if isinstance(other, (np.int64, np.bool_)):
            # for numpy scalars we get a deprecation warning and False as result
            # for equality or error for larger/lesser than
            continue

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


def test_is_scalar():
    assert is_scalar(NA) is True


def test_isna():
    assert pd.isna(NA) is True
    assert pd.notna(NA) is False


def test_series_isna():
    s = pd.Series([1, NA], dtype=object)
    expected = pd.Series([False, True])
    tm.assert_series_equal(s.isna(), expected)


def test_integer_hash_collision_dict():
    # GH 30013
    result = {NA: "foo", hash(NA): "bar"}

    assert result[NA] == "foo"
    assert result[hash(NA)] == "bar"


def test_integer_hash_collision_set():
    # GH 30013
    result = {NA, hash(NA)}

    assert len(result) == 2
    assert NA in result
    assert hash(NA) in result
