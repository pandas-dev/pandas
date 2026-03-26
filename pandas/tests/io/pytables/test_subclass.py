import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm

from pandas.io.pytables import (
    HDFStore,
    read_hdf,
)

pytest.importorskip("tables")


class TestHDFStoreSubclass:
    # GH 33748
    def test_supported_for_subclass_dataframe(self, temp_h5_path):
        data = {"a": [1, 2], "b": [3, 4]}
        sdf = tm.SubclassedDataFrame(data, dtype=np.intp)

        expected = DataFrame(data, dtype=np.intp)

        sdf.to_hdf(temp_h5_path, key="df")
        result = read_hdf(temp_h5_path, "df")
        tm.assert_frame_equal(result, expected)

        with HDFStore(temp_h5_path) as store:
            store.put("df", sdf)
        result = read_hdf(temp_h5_path, "df")
        tm.assert_frame_equal(result, expected)

    def test_supported_for_subclass_series(self, temp_h5_path):
        data = [1, 2, 3]
        sser = tm.SubclassedSeries(data, dtype=np.intp)

        expected = Series(data, dtype=np.intp)

        sser.to_hdf(temp_h5_path, key="ser")
        result = read_hdf(temp_h5_path, "ser")
        tm.assert_series_equal(result, expected)

        with HDFStore(temp_h5_path) as store:
            store.put("ser", sser)
        result = read_hdf(temp_h5_path, "ser")
        tm.assert_series_equal(result, expected)
