import numpy as np
from numpy import nan
import pytest

import pandas as pd
from pandas import Series
from pandas.core.sparse import series as sps
from pandas.util import testing as tm


@pytest.mark.filterwarnings("ignore:Sparse:FutureWarning")
def test_deprecated_to_sparse():
    ser = Series([1, np.nan, 3])
    sparse_ser = pd.SparseSeries([1, np.nan, 3])

    # Deprecated 0.25.0
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = ser.to_sparse()
    tm.assert_series_equal(result, sparse_ser)
