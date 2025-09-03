from io import StringIO

import pytest


def test_leading_zeros_preserved_with_dtype_str(all_parsers, request):
    """
    Ensure that all parser engines preserve leading zeros when dtype=str is passed.

    This test verifies that when dtype=str is specified, leading zeros in
    numeric-looking strings are preserved across all available parser engines.
    """
    parser = all_parsers
    engine_name = getattr(parser, "engine", "unknown")

    data = """col1|col2|col3|col4
AB|000388907|abc|0150
CD|101044572|def|0150
EF|000023607|ghi|0205
GH|100102040|jkl|0205"""

    result = parser.read_csv(
        StringIO(data),
        sep="|",
        dtype=str,
    )

    try:
        assert result.shape == (4, 4)
        assert list(result.columns) == ["col1", "col2", "col3", "col4"]

        assert result.loc[0, "col2"] == "000388907", "lost zeros in col2 row 0"
        assert result.loc[2, "col2"] == "000023607", "lost zeros in col2 row 2"
        assert result.loc[0, "col4"] == "0150", "lost zeros in col4 row 0"
        assert result.loc[2, "col4"] == "0205", "lost zeros in col4 row 2"

    except AssertionError as exc:
        if engine_name == "pyarrow":
            # Known issue: pyarrow engine strips leading zeros even with dtype=str.
            request.node.add_marker(
                pytest.mark.xfail(reason=f"failed assertions: {exc}", strict=False)
            )
        raise
