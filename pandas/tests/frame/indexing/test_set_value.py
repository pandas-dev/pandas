import numpy as np
import pytest

from pandas.core.dtypes.common import is_float_dtype

from pandas import isna


class TestSetValue:
    def test_set_value(self, float_frame):
        for idx in float_frame.index:
            for col in float_frame.columns:
                float_frame._set_value(idx, col, 1)
                assert float_frame[col][idx] == 1

    def test_set_value_resize(self, float_frame):

        res = float_frame._set_value("foobar", "B", 0)
        assert res is None
        assert float_frame.index[-1] == "foobar"
        assert float_frame._get_value("foobar", "B") == 0

        float_frame.loc["foobar", "qux"] = 0
        assert float_frame._get_value("foobar", "qux") == 0

        res = float_frame.copy()
        res._set_value("foobar", "baz", "sam")
        assert res["baz"].dtype == np.object_

        res = float_frame.copy()
        res._set_value("foobar", "baz", True)
        assert res["baz"].dtype == np.object_

        res = float_frame.copy()
        res._set_value("foobar", "baz", 5)
        assert is_float_dtype(res["baz"])
        assert isna(res["baz"].drop(["foobar"])).all()
        msg = "could not convert string to float: 'sam'"
        with pytest.raises(ValueError, match=msg):
            res._set_value("foobar", "baz", "sam")
