import pytest


@pytest.mark.parametrize("dtype", ["int64", "uint64", "float64", "category"])
def test_astype_preserves_name(indices, dtype):
    try:
        result = indices.astype(dtype)
    except:
        return

    assert result.name == indices.name


@pytest.mark.parametrize("dtype", ["int64", "uint64", "float64", "category"])
def test_astype_with_copy_preserves_name(indices, dtype):
    try:
        result = indices.copy(dtype=dtype)
    except:
        return

    assert result.name == indices.name
