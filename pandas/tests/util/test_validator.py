import numpy as np
import pytest

import pandas as pd


def test_series_quantile_invalid_percentile():
    s = pd.Series(np.random.randn(100))
    percentile_array = [-0.5, 0.25, 1.5]
    msg = "percentiles should all be in the interval \\[0, 1\\]"
    with pytest.raises(ValueError, match=msg):
        s.quantile(percentile_array)
