import pytest
import numpy as np
import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize('test_dtype', [object, 'int8', 'int16', 'int32', 'int64'])
def test_dtypes(test_dtype):
    df = pd.DataFrame({'A': pd.Series([1, 2, 3], dtype=test_dtype), 'B': [1, 2, 3]})
    expected = df.dtypes.values[0].type

    df.set_index('A')
    result = df.index.dtype.type
    assert result == expected
