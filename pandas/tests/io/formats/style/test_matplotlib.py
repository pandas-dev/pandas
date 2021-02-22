import numpy as np
import pytest

from pandas import (
    DataFrame,
    IndexSlice,
    Series,
)

pytest.importorskip("matplotlib")
pytest.importorskip("jinja2")


class TestStylerMatplotlibDep:
    def test_background_gradient(self):
        df = DataFrame([[1, 2], [2, 4]], columns=["A", "B"])

        for c_map in [None, "YlOrRd"]:
            result = df.style.background_gradient(cmap=c_map)._compute().ctx
            assert all("#" in x[0][1] for x in result.values())
            assert result[(0, 0)] == result[(0, 1)]
            assert result[(1, 0)] == result[(1, 1)]

        result = df.style.background_gradient(subset=IndexSlice[1, "A"])._compute().ctx

        assert result[(1, 0)] == [("background-color", "#fff7fb"), ("color", "#000000")]

    @pytest.mark.parametrize(
        "cmap, expected",
        [
            (
                "PuBu",
                {
                    (4, 5): [("background-color", "#86b0d3"), ("color", "#000000")],
                    (4, 6): [("background-color", "#83afd3"), ("color", "#f1f1f1")],
                },
            ),
            (
                "YlOrRd",
                {
                    (4, 8): [("background-color", "#fd913e"), ("color", "#000000")],
                    (4, 9): [("background-color", "#fd8f3d"), ("color", "#f1f1f1")],
                },
            ),
            (
                None,
                {
                    (7, 0): [("background-color", "#48c16e"), ("color", "#f1f1f1")],
                    (7, 1): [("background-color", "#4cc26c"), ("color", "#000000")],
                },
            ),
        ],
    )
    def test_text_color_threshold(self, cmap, expected):
        df = DataFrame(np.arange(100).reshape(10, 10))
        result = df.style.background_gradient(cmap=cmap, axis=None)._compute().ctx
        for k in expected.keys():
            assert result[k] == expected[k]

    @pytest.mark.parametrize("text_color_threshold", [1.1, "1", -1, [2, 2]])
    def test_text_color_threshold_raises(self, text_color_threshold):
        df = DataFrame([[1, 2], [2, 4]], columns=["A", "B"])
        msg = "`text_color_threshold` must be a value from 0 to 1."
        with pytest.raises(ValueError, match=msg):
            df.style.background_gradient(
                text_color_threshold=text_color_threshold
            )._compute()

    def test_background_gradient_axis(self):
        df = DataFrame([[1, 2], [2, 4]], columns=["A", "B"])

        low = [("background-color", "#f7fbff"), ("color", "#000000")]
        high = [("background-color", "#08306b"), ("color", "#f1f1f1")]
        mid = [("background-color", "#abd0e6"), ("color", "#000000")]
        result = df.style.background_gradient(cmap="Blues", axis=0)._compute().ctx
        assert result[(0, 0)] == low
        assert result[(0, 1)] == low
        assert result[(1, 0)] == high
        assert result[(1, 1)] == high

        result = df.style.background_gradient(cmap="Blues", axis=1)._compute().ctx
        assert result[(0, 0)] == low
        assert result[(0, 1)] == high
        assert result[(1, 0)] == low
        assert result[(1, 1)] == high

        result = df.style.background_gradient(cmap="Blues", axis=None)._compute().ctx
        assert result[(0, 0)] == low
        assert result[(0, 1)] == mid
        assert result[(1, 0)] == mid
        assert result[(1, 1)] == high

    def test_background_gradient_vmin_vmax(self):
        # GH 12145
        df = DataFrame(range(5))
        ctx = df.style.background_gradient(vmin=1, vmax=3)._compute().ctx
        assert ctx[(0, 0)] == ctx[(1, 0)]
        assert ctx[(4, 0)] == ctx[(3, 0)]

    def test_background_gradient_int64(self):
        # GH 28869
        df1 = Series(range(3)).to_frame()
        df2 = Series(range(3), dtype="Int64").to_frame()
        ctx1 = df1.style.background_gradient()._compute().ctx
        ctx2 = df2.style.background_gradient()._compute().ctx
        assert ctx2[(0, 0)] == ctx1[(0, 0)]
        assert ctx2[(1, 0)] == ctx1[(1, 0)]
        assert ctx2[(2, 0)] == ctx1[(2, 0)]
