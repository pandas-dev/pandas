from types import SimpleNamespace

import pytest

from pandas.core.dtypes.common import is_float

import pandas._testing as tm


def test_assert_attr_equal(nulls_fixture):
    obj = SimpleNamespace()
    obj.na_value = nulls_fixture
    tm.assert_attr_equal("na_value", obj, obj)


def test_assert_attr_equal_different_nulls(nulls_fixture, nulls_fixture2):
    obj = SimpleNamespace()
    obj.na_value = nulls_fixture

    obj2 = SimpleNamespace()
    obj2.na_value = nulls_fixture2

    if nulls_fixture is nulls_fixture2:
        tm.assert_attr_equal("na_value", obj, obj2)
    elif is_float(nulls_fixture) and is_float(nulls_fixture2):
        # we consider float("nan") and np.float64("nan") to be equivalent
        tm.assert_attr_equal("na_value", obj, obj2)
    elif type(nulls_fixture) is type(nulls_fixture2):
        # e.g. Decimal("NaN")
        tm.assert_attr_equal("na_value", obj, obj2)
    else:
        with pytest.raises(AssertionError, match='"na_value" are different'):
            tm.assert_attr_equal("na_value", obj, obj2)


def test_assert_attr_equal_tuple_raising_type_error():
    class TypeErrorOnComparison:
        def __eq__(self, other):
            raise TypeError(f"Cannot compare with {other!r}")

    left = SimpleNamespace(name=(TypeErrorOnComparison(),))
    right = SimpleNamespace(name=(TypeErrorOnComparison(),))

    with pytest.raises(AssertionError, match='"name" are different'):
        tm.assert_attr_equal("name", left, right)


@pytest.mark.parametrize(
    "left_name,right_name",
    [
        # ragged, so the nested objects survive conversion to an object array
        (("a", (1, 2)), ("a", range(1, 3))),
        # equal-length, so a 2D conversion would discard the nested types
        (((1, 2), (3, 4)), ([1, 2], [3, 4])),
        (((1, 2), (3, 4)), (range(1, 3), range(3, 5))),
    ],
)
def test_assert_attr_equal_tuple_nested_type_mismatch(left_name, right_name):
    # GH#54521 comparing tuple attributes NA-aware must not treat nested
    # sequences of differing type as equal, the way == does not.
    left = SimpleNamespace(name=left_name)
    right = SimpleNamespace(name=right_name)

    with pytest.raises(AssertionError, match='"name" are different'):
        tm.assert_attr_equal("name", left, right)
