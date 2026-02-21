"""
Tests that apply specifically to the PyArrow parser.
"""

from io import StringIO

import pytest

pa = pytest.importorskip("pyarrow")


@pytest.mark.parametrize("split_blocks", [True, False])
def test_to_pandas_kwargs_split_blocks(pyarrow_parser_only, split_blocks):
    # GH#34823
    # split_blocks=True prevents consolidation of same-dtype columns
    data = "a,b\n1,2\n3,4"

    result = pyarrow_parser_only.read_csv(
        StringIO(data),
        to_pandas_kwargs={"split_blocks": split_blocks},
    )
    assert list(result.columns) == ["a", "b"]
    assert len(result) == 2
    # With split_blocks=True, each column should be in its own block
    assert len(result._mgr.blocks) == len(result.columns) if split_blocks else 1


def test_to_pandas_kwargs_zero_copy_only_raises(pyarrow_parser_only):
    # zero_copy_only=True raises if zero-copy conversion not possible
    data = "a,b\n1,2\n3,4"

    # zero_copy_only=True raises without split_blocks for multi-column data
    with pytest.raises(pa.lib.ArrowInvalid, match="zero copy"):
        pyarrow_parser_only.read_csv(
            StringIO(data),
            to_pandas_kwargs={"zero_copy_only": True},
        )


@pytest.mark.parametrize("zero_copy_only", [True, False])
def test_to_pandas_kwargs_zero_copy_only_success(pyarrow_parser_only, zero_copy_only):
    # GH#34823
    # zero_copy_only with split_blocks=True enables zero-copy conversion
    # No exception means pyarrow confirmed zero-copy conversion is possible
    data = "a,b\n1,2\n3,4"

    result = pyarrow_parser_only.read_csv(
        StringIO(data),
        to_pandas_kwargs={"zero_copy_only": zero_copy_only, "split_blocks": True},
    )
    assert list(result.columns) == ["a", "b"]
    assert len(result) == 2
    # Zero-copy arrays share memory with pyarrow and are not writeable
    assert not result["a"].values.flags.writeable
    assert not result["b"].values.flags.writeable
