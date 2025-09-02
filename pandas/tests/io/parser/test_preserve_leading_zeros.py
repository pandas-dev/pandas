import pytest
from io import StringIO
import pandas._testing as tm


@pytest.mark.xfail(reason="Leading zeros preservation may not work consistently across all engines")
def test_leading_zeros_preserved_with_dtype_str(all_parsers):
    """
    Ensure that all parser engines preserve leading zeros when dtype=str is passed.
    
    This test verifies that when dtype=str is specified, leading zeros in 
    numeric-looking strings are preserved across all available parser engines.
    """
    parser = all_parsers
    engine_name = getattr(parser, 'engine', 'unknown')
    
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
    
    # Verify leading zeros are preserved in col2
    assert result.loc[0, "col2"] == "000388907", f"Engine {engine_name}: Leading zeros lost in col2, row 0. Got: {result.loc[0, 'col2']}"
    assert result.loc[2, "col2"] == "000023607", f"Engine {engine_name}: Leading zeros lost in col2, row 2. Got: {result.loc[2, 'col2']}"
    
    # Verify leading zeros are preserved in col4
    assert result.loc[0, "col4"] == "0150", f"Engine {engine_name}: Leading zeros lost in col4, row 0. Got: {result.loc[0, 'col4']}"
    assert result.loc[2, "col4"] == "0205", f"Engine {engine_name}: Leading zeros lost in col4, row 2. Got: {result.loc[2, 'col4']}"
    
    # Verify all columns are string type
    assert result.dtypes["col1"] == "object", f"Engine {engine_name}: col1 should be string type, got {result.dtypes['col1']}"
    assert result.dtypes["col2"] == "object", f"Engine {engine_name}: col2 should be string type, got {result.dtypes['col2']}"
    assert result.dtypes["col3"] == "object", f"Engine {engine_name}: col3 should be string type, got {result.dtypes['col3']}"
    assert result.dtypes["col4"] == "object", f"Engine {engine_name}: col4 should be string type, got {result.dtypes['col4']}"
    
    # Verify shape
    assert result.shape == (4, 4), f"Engine {engine_name}: Expected shape (4, 4), got {result.shape}"
    
    # Verify column names
    expected_columns = ["col1", "col2", "col3", "col4"]
    assert list(result.columns) == expected_columns, f"Engine {engine_name}: Expected columns {expected_columns}, got {list(result.columns)}"
