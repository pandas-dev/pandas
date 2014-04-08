import numpy as np
import pandas as pd
from pandas.util.testing import assert_series_equal
import warnings

def test_describe_with_warnings_raised():
    with warnings.catch_warnings():
        # escalate warnings
        warnings.simplefilter("error")
        df = pd.DataFrame({"A": [1, 2, 3], "B": [1.2, 4.2, 5.2]})
        act = df.groupby('A')['B'].describe().unstack(0)
        exp = pd.Series(
            [1.0, 1.2, np.NaN, 1.2, 1.2, 1.2, 1.2, 1.2],
            ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'],
        )
        assert_series_equal(act[1], exp)

