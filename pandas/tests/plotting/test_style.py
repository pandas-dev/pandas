import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    IndexSlice,
    Series,
)

pytest.importorskip("matplotlib")
from pandas.plotting._matplotlib.style import get_standard_colors

pytestmark = pytest.mark.slow


class TestGetStandardColors:
    @pytest.mark.parametrize(
        "num_colors, expected",
        [
            (3, ["red", "green", "blue"]),
            (5, ["red", "green", "blue", "red", "green"]),
            (7, ["red", "green", "blue", "red", "green", "blue", "red"]),
            (2, ["red", "green"]),
            (1, ["red"]),
        ],
    )
    def test_default_colors_named_from_prop_cycle(self, num_colors, expected):
        import matplotlib as mpl
        from matplotlib.pyplot import cycler

        mpl_params = {
            "axes.prop_cycle": cycler(color=["red", "green", "blue"]),
        }
        with mpl.rc_context(rc=mpl_params):
            result = get_standard_colors(num_colors=num_colors)
            assert result == expected

    @pytest.mark.parametrize(
        "num_colors, expected",
        [
            (1, ["b"]),
            (3, ["b", "g", "r"]),
            (4, ["b", "g", "r", "y"]),
            (5, ["b", "g", "r", "y", "b"]),
            (7, ["b", "g", "r", "y", "b", "g", "r"]),
        ],
    )
    def test_default_colors_named_from_prop_cycle_string(self, num_colors, expected):
        import matplotlib as mpl
        from matplotlib.pyplot import cycler

        mpl_params = {
            "axes.prop_cycle": cycler(color="bgry"),
        }
        with mpl.rc_context(rc=mpl_params):
            result = get_standard_colors(num_colors=num_colors)
            assert result == expected

    @pytest.mark.parametrize(
        "num_colors, expected_name",
        [
            (1, ["C0"]),
            (3, ["C0", "C1", "C2"]),
            (
                12,
                [
                    "C0",
                    "C1",
                    "C2",
                    "C3",
                    "C4",
                    "C5",
                    "C6",
                    "C7",
                    "C8",
                    "C9",
                    "C0",
                    "C1",
                ],
            ),
        ],
    )
    def test_default_colors_named_undefined_prop_cycle(self, num_colors, expected_name):
        import matplotlib as mpl
        import matplotlib.colors as mcolors

        with mpl.rc_context(rc={}):
            expected = [mcolors.to_hex(x) for x in expected_name]
            result = get_standard_colors(num_colors=num_colors)
            assert result == expected

    @pytest.mark.parametrize(
        "num_colors, expected",
        [
            (1, ["red", "green", (0.1, 0.2, 0.3)]),
            (2, ["red", "green", (0.1, 0.2, 0.3)]),
            (3, ["red", "green", (0.1, 0.2, 0.3)]),
            (4, ["red", "green", (0.1, 0.2, 0.3), "red"]),
        ],
    )
    def test_user_input_color_sequence(self, num_colors, expected):
        color = ["red", "green", (0.1, 0.2, 0.3)]
        result = get_standard_colors(color=color, num_colors=num_colors)
        assert result == expected

    @pytest.mark.parametrize(
        "num_colors, expected",
        [
            (1, ["r", "g", "b", "k"]),
            (2, ["r", "g", "b", "k"]),
            (3, ["r", "g", "b", "k"]),
            (4, ["r", "g", "b", "k"]),
            (5, ["r", "g", "b", "k", "r"]),
            (6, ["r", "g", "b", "k", "r", "g"]),
        ],
    )
    def test_user_input_color_string(self, num_colors, expected):
        color = "rgbk"
        result = get_standard_colors(color=color, num_colors=num_colors)
        assert result == expected

    @pytest.mark.parametrize(
        "num_colors, expected",
        [
            (1, [(0.1, 0.2, 0.3)]),
            (2, [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3)]),
            (3, [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.1, 0.2, 0.3)]),
        ],
    )
    def test_user_input_color_floats(self, num_colors, expected):
        color = (0.1, 0.2, 0.3)
        result = get_standard_colors(color=color, num_colors=num_colors)
        assert result == expected

    @pytest.mark.parametrize(
        "color, num_colors, expected",
        [
            ("Crimson", 1, ["Crimson"]),
            ("DodgerBlue", 2, ["DodgerBlue", "DodgerBlue"]),
            ("firebrick", 3, ["firebrick", "firebrick", "firebrick"]),
        ],
    )
    def test_user_input_named_color_string(self, color, num_colors, expected):
        result = get_standard_colors(color=color, num_colors=num_colors)
        assert result == expected

    @pytest.mark.parametrize("color", ["", [], (), Series([], dtype="object")])
    def test_empty_color_raises(self, color):
        with pytest.raises(ValueError, match="Invalid color argument"):
            get_standard_colors(color=color, num_colors=1)

    @pytest.mark.parametrize(
        "color",
        [
            "bad_color",
            ("red", "green", "bad_color"),
            (0.1,),
            (0.1, 0.2),
            (0.1, 0.2, 0.3, 0.4, 0.5),  # must be either 3 or 4 floats
        ],
    )
    def test_bad_color_raises(self, color):
        with pytest.raises(ValueError, match="Invalid color"):
            get_standard_colors(color=color, num_colors=5)


@td.skip_if_no_mpl
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
        "c_map,expected",
        [
            (
                None,
                {
                    (0, 0): [("background-color", "#440154"), ("color", "#f1f1f1")],
                    (1, 0): [("background-color", "#fde725"), ("color", "#000000")],
                },
            ),
            (
                "YlOrRd",
                {
                    (0, 0): [("background-color", "#ffffcc"), ("color", "#000000")],
                    (1, 0): [("background-color", "#800026"), ("color", "#f1f1f1")],
                },
            ),
        ],
    )
    def test_text_color_threshold(self, c_map, expected):
        df = DataFrame([1, 2], columns=["A"])
        result = df.style.background_gradient(cmap=c_map)._compute().ctx
        assert result == expected

    @pytest.mark.parametrize("text_color_threshold", [1.1, "1", -1, [2, 2]])
    def test_text_color_threshold_raises(self, text_color_threshold):
        df = DataFrame([[1, 2], [2, 4]], columns=["A", "B"])
        msg = "`text_color_threshold` must be a value from 0 to 1."
        with pytest.raises(ValueError, match=msg):
            df.style.background_gradient(
                text_color_threshold=text_color_threshold
            )._compute()

    @td.skip_if_no_mpl
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
