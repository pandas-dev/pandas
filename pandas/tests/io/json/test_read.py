import pytest
import pandas as pd
import io
import json
import io
import json
import pandas as pd
import pytest


def test_read_json_type_inference_cases():
    d = [{"col1": 1, "col2": 1.0, "col3": "1"}, {"col1": 2, "col2": 2.0, "col3": "2"}]

    # Default engine (ujson): float column should stay float64
    df = pd.read_json(io.StringIO(json.dumps(d)))
    assert df["col1"].dtype == "int64"
    assert df["col2"].dtype == "float64"

    # Default engine (ujson) with dtype=False: prevents string conversion
    df_no_dtype = pd.read_json(io.StringIO(json.dumps(d)), dtype=False)
    assert df_no_dtype["col1"].dtype == "int64"
    assert df_no_dtype["col2"].dtype == "float64"
    assert df_no_dtype["col3"].dtype == "str"

    # PyArrow engine (json_lines): preserves floats and strings
    pytest.importorskip("pyarrow")
    json_lines = b'{"col1": 1, "col2": 1.0, "col3": "1"}\n{"col1": 2, "col2": 2.0, "col3": "2"}'
    df_pa = pd.read_json(io.BytesIO(json_lines), lines=True, engine="pyarrow")
    assert df_pa["col1"].dtype == "int64"
    assert df_pa["col2"].dtype == "float64"
    assert df_pa["col3"].dtype == "str"

    
