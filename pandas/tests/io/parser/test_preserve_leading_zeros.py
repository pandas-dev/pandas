from io import StringIO

import pytest


def test_leading_zeros_preserved_with_dtype_str(all_parsers, request):
    # GH#57666: pyarrow engine strips leading zeros when dtype=str is passed
    # GH#61618: further discussion on ensuring string dtype preservation across engines

    parser = all_parsers
    engine_name = getattr(parser, "engine", "unknown")

    data = """col1,col2,col3,col4
AB,000388907,abc,0150
CD,101044572,def,0150
EF,000023607,ghi,0205
GH,100102040,jkl,0205"""

    result = parser.read_csv(
        StringIO(data),
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
            # Temporary workaround for GH#57666
            # Remove once type preservation is fixed in pyarrow engine.
            request.node.add_marker(
                pytest.mark.xfail(reason=f"failed assertions: {exc}", strict=False)
            )
        raise


def test_leading_zeros_preserved_with_dtype_dict(all_parsers):
    # GH#57666: pyarrow engine strips leading zeros when dtype=str is passed
    # GH#61618: further discussion on ensuring string dtype preservation across engines

    parser = all_parsers

    data = """col1,col2,col3,col4
AB,000388907,199,0150
CD,101044572,200,0150
EF,000023607,201,0205
GH,100102040,202,0205"""

    result = parser.read_csv(
        StringIO(data),
        dtype={"col2": str, "col4": str},
    )

    assert result.shape == (4, 4)
    assert list(result.columns) == ["col1", "col2", "col3", "col4"]

    assert result.loc[0, "col2"] == "000388907", "lost zeros in col2 row 0"
    assert result.loc[2, "col2"] == "000023607", "lost zeros in col2 row 2"
    assert result.loc[0, "col4"] == "0150", "lost zeros in col4 row 0"
    assert result.loc[2, "col4"] == "0205", "lost zeros in col4 row 2"

    assert result.loc[0, "col3"] == 199
    assert result.loc[1, "col3"] == 200
    assert result.loc[2, "col3"] == 201
    assert result.loc[3, "col3"] == 202
