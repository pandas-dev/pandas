from typing import Any
from unittest.mock import PropertyMock, patch

import pytest

from pandas.core.computation.pytables import BinOp, TermValue
from pandas.core.series import Series


@patch(
    "pandas.core.computation.pytables.BinOp.kind",
    new_callable=PropertyMock,
    return_value="integer",
)
@patch(
    "pandas.core.computation.pytables.BinOp.meta",
    new_callable=PropertyMock,
    return_value="category",
)
@patch(
    "pandas.core.computation.pytables.BinOp.metadata",
    new_callable=PropertyMock,
    return_value=Series(data=["a", "b", "s"]),
)
@pytest.mark.parametrize(
    "value, expected_results",
    [("q", TermValue(-1, -1, "integer")), ("a", TermValue(0, 0, "integer"))],
)
def test_convert_value(
    mock_kind, mock_meta, mock_metadata, value: Any, expected_results: TermValue
):

    with patch.object(BinOp, "__init__", lambda p1, p2, p3, p4, p5, p6: None):
        bin_op = BinOp(None, None, None, None, None)
        bin_op.encoding = "UTF-8"

        result = bin_op.convert_value(value)

        assert (
            result.kind == expected_results.kind
            and result.value == expected_results.value
            and result.converted == expected_results.converted
        )
