import decimal
import pytest
import pandas as pd
import pandas._testing as tm

try:
    import pyarrow as pa
except ImportError:
    pytest.skip("pyarrow not installed", allow_module_level=True)

def test_groupby_var_arrow_dtype_GH54627():
    # GH#54627
    df = pd.DataFrame({
        "A": pd.Series([True, True], dtype="bool[pyarrow]"),
        "B": pd.Series([decimal.Decimal(123), decimal.Decimal(12)], 
                       dtype=pd.ArrowDtype(pa.decimal128(6,3)))
    })
    result = df.groupby("A").var()
    
    # Check if the dtype is Arrow-backed float64
    assert isinstance(result["B"].dtype, pd.ArrowDtype)
    assert result["B"].dtype.pyarrow_dtype == pa.float64()
