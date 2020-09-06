import operator

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.frame.common import zip_frames


def test_agg_transform(axis, float_frame):
    other_axis = 1 if axis in {0, "index"} else 0

    with np.errstate(all="ignore"):

        f_abs = np.abs(float_frame)
        f_sqrt = np.sqrt(float_frame)

        # ufunc
        result = float_frame.transform(np.sqrt, axis=axis)
        expected = f_sqrt.copy()
        tm.assert_frame_equal(result, expected)

        result = float_frame.transform(np.sqrt, axis=axis)
        tm.assert_frame_equal(result, expected)

        # list-like
        expected = f_sqrt.copy()
        if axis in {0, "index"}:
            expected.columns = pd.MultiIndex.from_product(
                [float_frame.columns, ["sqrt"]]
            )
        else:
            expected.index = pd.MultiIndex.from_product([float_frame.index, ["sqrt"]])
        result = float_frame.transform([np.sqrt], axis=axis)
        tm.assert_frame_equal(result, expected)

        # multiple items in list
        # these are in the order as if we are applying both
        # functions per series and then concatting
        expected = zip_frames([f_abs, f_sqrt], axis=other_axis)
        if axis in {0, "index"}:
            expected.columns = pd.MultiIndex.from_product(
                [float_frame.columns, ["absolute", "sqrt"]]
            )
        else:
            expected.index = pd.MultiIndex.from_product(
                [float_frame.index, ["absolute", "sqrt"]]
            )
        result = float_frame.transform([np.abs, "sqrt"], axis=axis)
        tm.assert_frame_equal(result, expected)


def test_transform_and_agg_err(axis, float_frame):
    # cannot both transform and agg
    msg = "transforms cannot produce aggregated results"
    with pytest.raises(ValueError, match=msg):
        float_frame.transform(["max", "min"], axis=axis)

    msg = "cannot combine transform and aggregation operations"
    with pytest.raises(ValueError, match=msg):
        with np.errstate(all="ignore"):
            float_frame.transform(["max", "sqrt"], axis=axis)


@pytest.mark.parametrize("method", ["abs", "shift", "pct_change", "cumsum", "rank"])
def test_transform_method_name(method):
    # GH 19760
    df = pd.DataFrame({"A": [-1, 2]})
    result = df.transform(method)
    expected = operator.methodcaller(method)(df)
    tm.assert_frame_equal(result, expected)
