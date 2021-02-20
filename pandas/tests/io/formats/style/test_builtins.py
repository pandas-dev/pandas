import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame

jinja2 = pytest.importorskip("jinja2")
from pandas.io.formats.style import Styler


def bar_grad(a=None, b=None, c=None, d=None):
    """Used in multiple tests to simplify formatting of expected result"""
    ret = [("width", "10em"), ("height", "80%")]
    if all(x is None for x in [a, b, c, d]):
        return ret
    return ret + [
        (
            "background",
            f"linear-gradient(90deg,{','.join(x for x in [a, b, c, d] if x)})",
        )
    ]


class TestStyler:
    def test_bar_align_left(self):
        df = DataFrame({"A": [0, 1, 2]})
        result = df.style.bar()._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad("#d65f5f 50.0%", " transparent 50.0%"),
            (2, 0): bar_grad("#d65f5f 100.0%", " transparent 100.0%"),
        }
        assert result == expected

        result = df.style.bar(color="red", width=50)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad("red 25.0%", " transparent 25.0%"),
            (2, 0): bar_grad("red 50.0%", " transparent 50.0%"),
        }
        assert result == expected

        df["C"] = ["a"] * len(df)
        result = df.style.bar(color="red", width=50)._compute().ctx
        assert result == expected
        df["C"] = df["C"].astype("category")
        result = df.style.bar(color="red", width=50)._compute().ctx
        assert result == expected

    def test_bar_align_left_0points(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = df.style.bar()._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (0, 1): bar_grad(),
            (0, 2): bar_grad(),
            (1, 0): bar_grad("#d65f5f 50.0%", " transparent 50.0%"),
            (1, 1): bar_grad("#d65f5f 50.0%", " transparent 50.0%"),
            (1, 2): bar_grad("#d65f5f 50.0%", " transparent 50.0%"),
            (2, 0): bar_grad("#d65f5f 100.0%", " transparent 100.0%"),
            (2, 1): bar_grad("#d65f5f 100.0%", " transparent 100.0%"),
            (2, 2): bar_grad("#d65f5f 100.0%", " transparent 100.0%"),
        }
        assert result == expected

        result = df.style.bar(axis=1)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (0, 1): bar_grad("#d65f5f 50.0%", " transparent 50.0%"),
            (0, 2): bar_grad("#d65f5f 100.0%", " transparent 100.0%"),
            (1, 0): bar_grad(),
            (1, 1): bar_grad("#d65f5f 50.0%", " transparent 50.0%"),
            (1, 2): bar_grad("#d65f5f 100.0%", " transparent 100.0%"),
            (2, 0): bar_grad(),
            (2, 1): bar_grad("#d65f5f 50.0%", " transparent 50.0%"),
            (2, 2): bar_grad("#d65f5f 100.0%", " transparent 100.0%"),
        }
        assert result == expected

    def test_bar_align_mid_pos_and_neg(self):
        df = DataFrame({"A": [-10, 0, 20, 90]})
        result = df.style.bar(align="mid", color=["#d65f5f", "#5fba7d"])._compute().ctx
        expected = {
            (0, 0): bar_grad(
                "#d65f5f 10.0%",
                " transparent 10.0%",
            ),
            (1, 0): bar_grad(),
            (2, 0): bar_grad(
                " transparent 10.0%",
                " #5fba7d 10.0%",
                " #5fba7d 30.0%",
                " transparent 30.0%",
            ),
            (3, 0): bar_grad(
                " transparent 10.0%",
                " #5fba7d 10.0%",
                " #5fba7d 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_all_pos(self):
        df = DataFrame({"A": [10, 20, 50, 100]})

        result = df.style.bar(align="mid", color=["#d65f5f", "#5fba7d"])._compute().ctx

        expected = {
            (0, 0): bar_grad(
                "#5fba7d 10.0%",
                " transparent 10.0%",
            ),
            (1, 0): bar_grad(
                "#5fba7d 20.0%",
                " transparent 20.0%",
            ),
            (2, 0): bar_grad(
                "#5fba7d 50.0%",
                " transparent 50.0%",
            ),
            (3, 0): bar_grad(
                "#5fba7d 100.0%",
                " transparent 100.0%",
            ),
        }

        assert result == expected

    def test_bar_align_mid_all_neg(self):
        df = DataFrame({"A": [-100, -60, -30, -20]})

        result = df.style.bar(align="mid", color=["#d65f5f", "#5fba7d"])._compute().ctx

        expected = {
            (0, 0): bar_grad(
                "#d65f5f 100.0%",
                " transparent 100.0%",
            ),
            (1, 0): bar_grad(
                " transparent 40.0%",
                " #d65f5f 40.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
            (2, 0): bar_grad(
                " transparent 70.0%",
                " #d65f5f 70.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
            (3, 0): bar_grad(
                " transparent 80.0%",
                " #d65f5f 80.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_zero_pos_and_neg(self):
        # See https://github.com/pandas-dev/pandas/pull/14757
        df = DataFrame({"A": [-10, 0, 20, 90]})

        result = (
            df.style.bar(align="zero", color=["#d65f5f", "#5fba7d"], width=90)
            ._compute()
            .ctx
        )
        expected = {
            (0, 0): bar_grad(
                " transparent 40.0%",
                " #d65f5f 40.0%",
                " #d65f5f 45.0%",
                " transparent 45.0%",
            ),
            (1, 0): bar_grad(),
            (2, 0): bar_grad(
                " transparent 45.0%",
                " #5fba7d 45.0%",
                " #5fba7d 55.0%",
                " transparent 55.0%",
            ),
            (3, 0): bar_grad(
                " transparent 45.0%",
                " #5fba7d 45.0%",
                " #5fba7d 90.0%",
                " transparent 90.0%",
            ),
        }
        assert result == expected

    def test_bar_align_left_axis_none(self):
        df = DataFrame({"A": [0, 1], "B": [2, 4]})
        result = df.style.bar(axis=None)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                "#d65f5f 25.0%",
                " transparent 25.0%",
            ),
            (0, 1): bar_grad(
                "#d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (1, 1): bar_grad(
                "#d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_zero_axis_none(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="zero", axis=None)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 50.0%",
                " #d65f5f 50.0%",
                " #d65f5f 62.5%",
                " transparent 62.5%",
            ),
            (0, 1): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (1, 1): bar_grad(
                " transparent 50.0%",
                " #d65f5f 50.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_axis_none(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 33.3%",
                " #d65f5f 33.3%",
                " #d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (0, 1): bar_grad(
                "#d65f5f 33.3%",
                " transparent 33.3%",
            ),
            (1, 1): bar_grad(
                " transparent 33.3%",
                " #d65f5f 33.3%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_vmin(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None, vmin=-6)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 60.0%",
                " #d65f5f 60.0%",
                " #d65f5f 70.0%",
                " transparent 70.0%",
            ),
            (0, 1): bar_grad(
                " transparent 40.0%",
                " #d65f5f 40.0%",
                " #d65f5f 60.0%",
                " transparent 60.0%",
            ),
            (1, 1): bar_grad(
                " transparent 60.0%",
                " #d65f5f 60.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_vmax(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None, vmax=8)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 20.0%",
                " #d65f5f 20.0%",
                " #d65f5f 30.0%",
                " transparent 30.0%",
            ),
            (0, 1): bar_grad(
                "#d65f5f 20.0%",
                " transparent 20.0%",
            ),
            (1, 1): bar_grad(
                " transparent 20.0%",
                " #d65f5f 20.0%",
                " #d65f5f 60.0%",
                " transparent 60.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_vmin_vmax_wide(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None, vmin=-3, vmax=7)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 30.0%",
                " #d65f5f 30.0%",
                " #d65f5f 40.0%",
                " transparent 40.0%",
            ),
            (0, 1): bar_grad(
                " transparent 10.0%",
                " #d65f5f 10.0%",
                " #d65f5f 30.0%",
                " transparent 30.0%",
            ),
            (1, 1): bar_grad(
                " transparent 30.0%",
                " #d65f5f 30.0%",
                " #d65f5f 70.0%",
                " transparent 70.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_vmin_vmax_clipping(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None, vmin=-1, vmax=3)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (0, 1): bar_grad("#d65f5f 25.0%", " transparent 25.0%"),
            (1, 1): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_nans(self):
        df = DataFrame({"A": [1, None], "B": [-1, 3]})
        result = df.style.bar(align="mid", axis=None)._compute().ctx
        expected = {
            (0, 0): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (0, 1): bar_grad("#d65f5f 25.0%", " transparent 25.0%"),
            (1, 1): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_zero_nans(self):
        df = DataFrame({"A": [1, None], "B": [-1, 2]})
        result = df.style.bar(align="zero", axis=None)._compute().ctx
        expected = {
            (0, 0): bar_grad(
                " transparent 50.0%",
                " #d65f5f 50.0%",
                " #d65f5f 75.0%",
                " transparent 75.0%",
            ),
            (0, 1): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (1, 1): bar_grad(
                " transparent 50.0%",
                " #d65f5f 50.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_bad_align_raises(self):
        df = DataFrame({"A": [-100, -60, -30, -20]})
        msg = "`align` must be one of {'left', 'zero',' mid'}"
        with pytest.raises(ValueError, match=msg):
            df.style.bar(align="poorly", color=["#d65f5f", "#5fba7d"])

    def test_highlight_null(self):
        df = DataFrame({"A": [0, np.nan]})
        result = df.style.highlight_null()._compute().ctx
        expected = {(1, 0): [("background-color", "red")]}
        assert result == expected

    def test_highlight_null_subset(self):
        # GH 31345
        df = DataFrame({"A": [0, np.nan], "B": [0, np.nan]})
        result = (
            df.style.highlight_null(null_color="red", subset=["A"])
            .highlight_null(null_color="green", subset=["B"])
            ._compute()
            .ctx
        )
        expected = {
            (1, 0): [("background-color", "red")],
            (1, 1): [("background-color", "green")],
        }
        assert result == expected

    def test_highlight_max(self):
        df = DataFrame([[1, 2], [3, 4]], columns=["A", "B"])
        css_seq = [("background-color", "yellow")]
        # max(df) = min(-df)
        for max_ in [True, False]:
            if max_:
                attr = "highlight_max"
            else:
                df = -df
                attr = "highlight_min"
            result = getattr(df.style, attr)()._compute().ctx
            assert result[(1, 1)] == css_seq

            result = getattr(df.style, attr)(color="green")._compute().ctx
            assert result[(1, 1)] == [("background-color", "green")]

            result = getattr(df.style, attr)(subset="A")._compute().ctx
            assert result[(1, 0)] == css_seq

            result = getattr(df.style, attr)(axis=0)._compute().ctx
            expected = {
                (1, 0): css_seq,
                (1, 1): css_seq,
            }
            assert result == expected

            result = getattr(df.style, attr)(axis=1)._compute().ctx
            expected = {
                (0, 1): css_seq,
                (1, 1): css_seq,
            }
            assert result == expected

        # separate since we can't negate the strs
        df["C"] = ["a", "b"]
        result = df.style.highlight_max()._compute().ctx
        expected = {(1, 1): css_seq}

        result = df.style.highlight_min()._compute().ctx
        expected = {(0, 0): css_seq}

    @pytest.mark.parametrize(
        "ttips",
        [
            DataFrame(
                data=[["Min", "Max"], [np.nan, ""]],
                columns=["A", "B"],
                index=["a", "b"],
            ),
            DataFrame(data=[["Max", "Min"]], columns=["B", "A"], index=["a"]),
            DataFrame(
                data=[["Min", "Max", None]], columns=["A", "B", "C"], index=["a"]
            ),
        ],
    )
    def test_tooltip_render(self, ttips):
        # GH 21266
        df = DataFrame(data=[[0, 3], [1, 2]], columns=["A", "B"], index=["a", "b"])
        s = Styler(df, uuid_len=0).set_tooltips(ttips).render()

        # test tooltip table level class
        assert "#T__ .pd-t {\n  visibility: hidden;\n" in s

        # test 'Min' tooltip added
        assert (
            "#T__ #T__row0_col0:hover .pd-t {\n  visibility: visible;\n}\n"
            + '#T__ #T__row0_col0 .pd-t::after {\n  content: "Min";\n}'
            in s
        )
        assert (
            '<td id="T__row0_col0" class="data row0 col0" >0<span class="pd-t">'
            + "</span></td>"
            in s
        )

        # test 'Max' tooltip added
        assert (
            "#T__ #T__row0_col1:hover .pd-t {\n  visibility: visible;\n}\n"
            + '#T__ #T__row0_col1 .pd-t::after {\n  content: "Max";\n}'
            in s
        )
        assert (
            '<td id="T__row0_col1" class="data row0 col1" >3<span class="pd-t">'
            + "</span></td>"
            in s
        )

    def test_tooltip_reindex(self):
        # GH 39317
        df = DataFrame(
            data=[[0, 1, 2], [3, 4, 5], [6, 7, 8]], columns=[0, 1, 2], index=[0, 1, 2]
        )
        ttips = DataFrame(
            data=[["Mi", "Ma"], ["Mu", "Mo"]],
            columns=[0, 2],
            index=[0, 2],
        )
        s = Styler(df, uuid_len=0).set_tooltips(DataFrame(ttips)).render()
        assert '#T__ #T__row0_col0 .pd-t::after {\n  content: "Mi";\n}' in s
        assert '#T__ #T__row0_col2 .pd-t::after {\n  content: "Ma";\n}' in s
        assert '#T__ #T__row2_col0 .pd-t::after {\n  content: "Mu";\n}' in s
        assert '#T__ #T__row2_col2 .pd-t::after {\n  content: "Mo";\n}' in s

    def test_tooltip_ignored(self):
        # GH 21266
        df = DataFrame(data=[[0, 1], [2, 3]])
        s = Styler(df).set_tooltips_class("pd-t").render()  # no set_tooltips()
        assert '<style type="text/css">\n</style>' in s
        assert '<span class="pd-t"></span>' not in s

    def test_tooltip_class(self):
        # GH 21266
        df = DataFrame(data=[[0, 1], [2, 3]])
        s = (
            Styler(df, uuid_len=0)
            .set_tooltips(DataFrame([["tooltip"]]))
            .set_tooltips_class(name="other-class", properties=[("color", "green")])
            .render()
        )
        assert "#T__ .other-class {\n  color: green;\n" in s
        assert '#T__ #T__row0_col0 .other-class::after {\n  content: "tooltip";\n' in s

        # GH 39563
        s = (
            Styler(df, uuid_len=0)
            .set_tooltips(DataFrame([["tooltip"]]))
            .set_tooltips_class(name="other-class", properties="color:green;color:red;")
            .render()
        )
        assert "#T__ .other-class {\n  color: green;\n  color: red;\n}" in s


@td.skip_if_no_mpl
class TestStylerMatplotlibDep:
    def test_background_gradient(self):
        df = DataFrame([[1, 2], [2, 4]], columns=["A", "B"])

        for c_map in [None, "YlOrRd"]:
            result = df.style.background_gradient(cmap=c_map)._compute().ctx
            assert all("#" in x[0][1] for x in result.values())
            assert result[(0, 0)] == result[(0, 1)]
            assert result[(1, 0)] == result[(1, 1)]

        result = (
            df.style.background_gradient(subset=pd.IndexSlice[1, "A"])._compute().ctx
        )

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
        df1 = pd.Series(range(3)).to_frame()
        df2 = pd.Series(range(3), dtype="Int64").to_frame()
        ctx1 = df1.style.background_gradient()._compute().ctx
        ctx2 = df2.style.background_gradient()._compute().ctx
        assert ctx2[(0, 0)] == ctx1[(0, 0)]
        assert ctx2[(1, 0)] == ctx1[(1, 0)]
        assert ctx2[(2, 0)] == ctx1[(2, 0)]
