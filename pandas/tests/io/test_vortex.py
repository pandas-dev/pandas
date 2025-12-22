import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm

pytestmark = pytest.mark.single_cpu


class TestVortexIO:
    def test_roundtrip(self, tmp_path):
        pytest.importorskip("vortex")
        from pandas import read_vortex

        df = DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]})
        path = tmp_path / "test.vortex"
        df.to_vortex(path)
        result = read_vortex(path)
        tm.assert_frame_equal(df, result)

    def test_columns(self, tmp_path):
        pytest.importorskip("vortex")
        from pandas import read_vortex

        df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
        path = tmp_path / "test.vortex"
        df.to_vortex(path)
        result = read_vortex(path, columns=["a", "c"])
        tm.assert_frame_equal(result, df[["a", "c"]])
