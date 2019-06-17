import numpy as np
import pytest

import pandas as pd
from pandas import Series
from pandas.util import testing as tm


@pytest.mark.filterwarnings("ignore:Sparse:FutureWarning")
def test_deprecated_to_sparse():
    # GH 26557
    # Deprecated 0.25.0

    ser = Series([1, np.nan, 3])
    sparse_ser = pd.SparseSeries([1, np.nan, 3])

    # Deprecated 0.25.0
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = ser.to_sparse()
    tm.assert_series_equal(result, sparse_ser)
