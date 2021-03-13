from textwrap import dedent

import pytest

from pandas import DataFrame

from pandas.io.formats.style import (
    _parse_latex_cell_styles,
    _parse_latex_table_styles,
)

pytest.importorskip("jinja2")


class TestStylerLatex:
    def setup_method(self, method):
        self.df = DataFrame({"A": [0, 1], "B": [-0.61, -1.22], "C": ["ab", "cd"]})

    def test_parse_latex_table_styles(self):
        s = self.df.style.set_table_styles(
            [
                {"selector": "foo", "props": [("attr", "value")]},
                {"selector": "bar", "props": [("attr", "overwritten")]},
                {"selector": "bar", "props": [("attr", "baz"), ("attr2", "ignored")]},
                {"selector": "label", "props": [("", "fig§item")]},
            ]
        )
        assert _parse_latex_table_styles(s.table_styles, "bar") == "baz"

        # test '§' replaced by ':' [for CSS compatibility]
        assert _parse_latex_table_styles(s.table_styles, "label") == "fig:item"

    def test_parse_latex_cell_styles_basic(self):
        cell_style = [("emph", ""), ("cellcolor", "[rgb]{0,1,1}")]
        expected = "\\emph{\\cellcolor[rgb]{0,1,1}{text}}"
        assert _parse_latex_cell_styles(cell_style, "text") == expected

    def test_parse_latex_cell_styles_braces_wrap(self):
        cell_style = [("Huge", "-wrap-"), ("cellcolor", "-wrap-[rgb]{0,1,1}")]
        expected = "{\\Huge {\\cellcolor[rgb]{0,1,1} text}}"
        assert _parse_latex_cell_styles(cell_style, "text") == expected

    def test_minimal_latex_tabular(self):
        s = self.df.style.format(precision=2)
        expected = dedent(
            """\
            \\begin{tabular}{lrrl}
             & A & B & C \\\\
            0 & 0 & -0.61 & ab \\\\
            1 & 1 & -1.22 & cd \\\\
            \\end{tabular}
            """
        )
        assert s.to_latex() == expected

    def test_latex_tabular_hrules(self):
        s = self.df.style.format(precision=2)
        expected = dedent(
            """\
            \\begin{tabular}{lrrl}
            \\toprule
             & A & B & C \\\\
            \\midrule
            0 & 0 & -0.61 & ab \\\\
            1 & 1 & -1.22 & cd \\\\
            \\bottomrule
            \\end{tabular}
            """
        )
        assert s.to_latex(hrules=True) == expected

    def test_latex_tabular_custom_hrules(self):
        s = self.df.style.format(precision=2)
        s.set_table_styles(
            [
                {"selector": "toprule", "props": ":hline"},
                {"selector": "bottomrule", "props": ":otherline"},
            ]
        )  # no midrule
        expected = dedent(
            """\
            \\begin{tabular}{lrrl}
            \\hline
             & A & B & C \\\\
            0 & 0 & -0.61 & ab \\\\
            1 & 1 & -1.22 & cd \\\\
            \\otherline
            \\end{tabular}
            """
        )
        assert s.to_latex() == expected
