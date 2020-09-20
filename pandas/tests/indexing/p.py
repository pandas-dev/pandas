import pandas as pd
import pytest
import numpy as np
from pandas.core.common import SettingWithCopyError

def test_loc_replace_subset_with_subset():
    # GH#36424 Should raise a warning
    df1 = pd.DataFrame(
        data=np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
        columns=['A', 'B', 'C'],
        index=[0, 0, 1]
    )
    df2 = df1.copy()
    df2[:] = np.nan
    # it fail if a warning is not raised
    with pytest.raises(SettingWithCopyError):
        df2.loc[0]['A'] = df1.loc[0]['A']

if __name__ == "__main__":
    test_loc_replace_subset_with_subset()