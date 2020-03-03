import pytest

from pandas.core.base import PandasObject

pandas_object = PandasObject()


class SubclassPandasObject(PandasObject):
    pass


subclass_pandas_object = SubclassPandasObject()


@pytest.mark.parametrize("other_object", [pandas_object, subclass_pandas_object])
def test_pandas_object_ensure_type(other_object):
    pandas_object = PandasObject()
    assert pandas_object._ensure_type(other_object)


def test_pandas_object_ensure_type_for_same_object():
    pandas_object_a = PandasObject()
    pandas_object_b = pandas_object_a
    assert pandas_object_a._ensure_type(pandas_object_b)


class OtherClass:
    pass


other_class = OtherClass()


@pytest.mark.parametrize("other_object", [other_class])
def test_pandas_object_ensure_type_for_false(other_object):
    pandas_object = PandasObject()
    with pytest.raises(AssertionError):
        assert pandas_object._ensure_type(other_object)
