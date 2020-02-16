import pytest

from pandas import Index


@pytest.mark.parametrize(
    "from_type, to_type", [("float64", "int64"), ("int64", "float64")]
)
def test_astype_preserves_name(from_type, to_type):
    idx = Index([1, 2], name="abc", dtype=from_type).astype(to_type)

    assert idx.name == "abc"


@pytest.mark.parametrize(
    "from_type, to_type", [("float64", "int64"), ("int64", "float64")]
)
def test_copy_with_astype_preserves_name(from_type, to_type):
    idx = Index([1, 2], name="abc", dtype=from_type).copy(dtype=to_type)

    assert idx.name == "abc"
