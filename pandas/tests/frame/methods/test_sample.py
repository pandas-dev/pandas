import numpy as np
import pytest

from pandas.compat.numpy import np_version_under1p17

from pandas import DataFrame
import pandas._testing as tm
import pandas.core.common as com


class TestSample:
    @pytest.mark.parametrize(
        "func_str,arg",
        [
            ("np.array", [2, 3, 1, 0]),
            pytest.param(
                "np.random.MT19937",
                3,
                marks=pytest.mark.skipif(np_version_under1p17, reason="NumPy<1.17"),
            ),
            pytest.param(
                "np.random.PCG64",
                11,
                marks=pytest.mark.skipif(np_version_under1p17, reason="NumPy<1.17"),
            ),
        ],
    )
    def test_sample_random_state(self, func_str, arg):
        # GH#32503
        df = DataFrame({"col1": range(10, 20), "col2": range(20, 30)})
        result = df.sample(n=3, random_state=eval(func_str)(arg))
        expected = df.sample(n=3, random_state=com.random_state(eval(func_str)(arg)))
        tm.assert_frame_equal(result, expected)

    def test_sample_upsampling_without_replacement(self):
        # GH#27451

        df = DataFrame({"A": list("abc")})
        msg = (
            "Replace has to be set to `True` when "
            "upsampling the population `frac` > 1."
        )
        with pytest.raises(ValueError, match=msg):
            df.sample(frac=2, replace=False)

    def test_sample_is_copy(self):
        # GH#27357, GH#30784: ensure the result of sample is an actual copy and
        # doesn't track the parent dataframe / doesn't give SettingWithCopy warnings
        df = DataFrame(np.random.randn(10, 3), columns=["a", "b", "c"])
        df2 = df.sample(3)

        with tm.assert_produces_warning(None):
            df2["d"] = 1
