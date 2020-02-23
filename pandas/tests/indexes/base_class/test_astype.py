import pytest


@pytest.mark.parametrize("dtype", ["int64", "uint64", "float64", "category"])
def test_astype_preserves_name(indices, dtype):
    indices.name = "idx"
    try:
        result = indices.astype(dtype)
    except (ValueError, TypeError, NotImplementedError):
        return

    assert result.name == indices.name


@pytest.mark.parametrize("dtype", ["int64", "uint64", "float64", "category"])
def test_astype_with_copy_preserves_name(indices, dtype):
    indices.name = "idx"
    try:
        result = indices.copy(dtype=dtype)
    except (ValueError, TypeError, NotImplementedError):
        return

    assert result.name == indices.name
