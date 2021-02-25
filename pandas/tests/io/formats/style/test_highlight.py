import numpy as np
import pytest

from pandas import DataFrame

pytest.importorskip("jinja2")


class TestStylerHighlight:
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
