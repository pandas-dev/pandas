import numpy as np
import pytest

import pandas as pd
from pandas.util import testing as tm


@pytest.mark.filterwarnings("ignore:Sparse:FutureWarning")
def test_deprecated_to_sparse():
    # GH 26557
    # Deprecated 0.25.0

    df = pd.DataFrame({"A": [1, np.nan, 3]})
    sparse_df = pd.SparseDataFrame({"A": [1, np.nan, 3]})

    # Deprecated 0.25.0
    with tm.assert_produces_warning(FutureWarning,
                                    check_stacklevel=False):
        result = df.to_sparse()
    tm.assert_frame_equal(result, sparse_df)
