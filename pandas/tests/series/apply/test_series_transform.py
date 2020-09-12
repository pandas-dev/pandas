import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_transform(string_series):
    # transforming functions

    with np.errstate(all="ignore"):
        f_sqrt = np.sqrt(string_series)
        f_abs = np.abs(string_series)

        # ufunc
        result = string_series.transform(np.sqrt)
        expected = f_sqrt.copy()
        tm.assert_series_equal(result, expected)

        # list-like
        result = string_series.transform([np.sqrt])
        expected = f_sqrt.to_frame().copy()
        expected.columns = ["sqrt"]
        tm.assert_frame_equal(result, expected)

        result = string_series.transform([np.sqrt])
        tm.assert_frame_equal(result, expected)

        result = string_series.transform(["sqrt"])
        tm.assert_frame_equal(result, expected)

        # multiple items in list
        # these are in the order as if we are applying both functions per
        # series and then concatting
        expected = pd.concat([f_sqrt, f_abs], axis=1)
        result = string_series.transform(["sqrt", "abs"])
        expected.columns = ["sqrt", "abs"]
        tm.assert_frame_equal(result, expected)


def test_transform_and_agg_error(string_series):
    # we are trying to transform with an aggregator
    msg = "transforms cannot produce aggregated results"
    with pytest.raises(ValueError, match=msg):
        string_series.transform(["min", "max"])

    msg = "cannot combine transform and aggregation operations"
    with pytest.raises(ValueError, match=msg):
        with np.errstate(all="ignore"):
            string_series.transform(["sqrt", "max"])


def test_transform_none_to_type():
    # GH34377
    df = pd.DataFrame({"a": [None]})

    msg = "DataFrame constructor called with incompatible data and dtype"
    with pytest.raises(TypeError, match=msg):
        df.transform({"a": int})
