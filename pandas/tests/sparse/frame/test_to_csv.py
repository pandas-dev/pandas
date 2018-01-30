import numpy as np
from pandas import SparseDataFrame, read_csv
from pandas.util import testing as tm


def test_to_csv_sparse_dataframe():
    fill_values = [0, 1, np.nan, None]

    for fill_value in fill_values:
        sdf = SparseDataFrame({'a': fill_values},
                              default_fill_value=fill_value)

        with tm.ensure_clean('sparse_df.csv') as path:
            sdf.to_csv(path, index=False)
            df = read_csv(path)

            tm.assert_sp_frame_equal(df.to_sparse(fill_value=fill_value), sdf)
