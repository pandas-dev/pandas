from textwrap import dedent

import pytest

from pandas import (
    DataFrame,
    MultiIndex,
)

from pandas.io.formats.style import (
    _parse_latex_cell_styles,
    _parse_latex_table_styles,
)

pytest.importorskip("jinja2")


class TestStylerLatex:
    def setup_method(self, method):
        self.df = DataFrame({"A": [0, 1], "B": [-0.61, -1.22], "C": ["ab", "cd"]})
        self.s = self.df.style.format(precision=2)

    def test_parse_latex_table_styles(self):
        s = self.df.style.set_table_styles(
            [
                {"selector": "foo", "props": [("attr", "value")]},
                {"selector": "bar", "props": [("attr", "overwritten")]},
                {"selector": "bar", "props": [("attr", "baz"), ("attr2", "ignored")]},
                {"selector": "label", "props": [("", "{fig§item}")]},
            ]
        )
        assert _parse_latex_table_styles(s.table_styles, "bar") == "baz"

        # test '§' replaced by ':' [for CSS compatibility]
        assert _parse_latex_table_styles(s.table_styles, "label") == "{fig:item}"

    def test_parse_latex_cell_styles_basic(self):
        cell_style = [("emph", ""), ("cellcolor", "[rgb]{0,1,1}")]
        expected = "\\emph{\\cellcolor[rgb]{0,1,1}{text}}"
        assert _parse_latex_cell_styles(cell_style, "text") == expected

    def test_parse_latex_cell_styles_braces_wrap(self):
        cell_style = [("emph", "-wrap-"), ("cellcolor", "-wrap-[rgb]{0,1,1}")]
        expected = "{\\emph {\\cellcolor[rgb]{0,1,1} text}}"
        assert _parse_latex_cell_styles(cell_style, "text") == expected

    def test_minimal_latex_tabular(self):
        expected = dedent(
            """\
            \\begin{tabular}{lrrl}
             & A & B & C \\\\
            0 & 0 & -0.61 & ab \\\\
            1 & 1 & -1.22 & cd \\\\
            \\end{tabular}
            """
        )
        assert self.s.to_latex() == expected

    def test_tabular_hrules(self):
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
        assert self.s.to_latex(hrules=True) == expected

    def test_tabular_custom_hrules(self):
        self.s.set_table_styles(
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
        assert self.s.to_latex() == expected

    def test_column_format(self):
        # default setting is already tested in `test_latex_minimal_tabular`
        self.s.set_table_styles([{"selector": "column_format", "props": ":cccc"}])

        assert "\\begin{tabular}{rrrr}" in self.s.to_latex(column_format="rrrr")
        self.s.set_table_styles([{"selector": "column_format", "props": ":rrrr"}])
        assert "\\begin{tabular}{rrrr}" in self.s.to_latex()

    def test_position(self):
        assert "\\begin{table}[h!]" in self.s.to_latex(position="h!")
        assert "\\end{table}" in self.s.to_latex(position="h!")
        self.s.set_table_styles([{"selector": "position", "props": ":h!"}])
        assert "\\begin{table}[h!]" in self.s.to_latex()
        assert "\\end{table}" in self.s.to_latex()

    def test_label(self):
        assert "\\label{text}" in self.s.to_latex(label="text")
        self.s.set_table_styles([{"selector": "label", "props": ":{text}"}])
        assert "\\label{text}" in self.s.to_latex()

    @pytest.mark.parametrize("label", [(None, ""), ("text", "\\label{text}")])
    @pytest.mark.parametrize("position", [(None, ""), ("h!", "{table}[h!]")])
    @pytest.mark.parametrize("caption", [(None, ""), ("text", "\\caption{text}")])
    @pytest.mark.parametrize("column_format", [(None, ""), ("rcrl", "{tabular}{rcrl}")])
    def test_kwargs_combinations(self, label, position, caption, column_format):
        result = self.s.to_latex(
            label=label[0],
            position=position[0],
            caption=caption[0],
            column_format=column_format[0],
        )
        assert label[1] in result
        assert position[1] in result
        assert caption[1] in result
        assert column_format[1] in result

    def test_custom_table_styles(self):
        self.s.set_table_styles(
            [
                {"selector": "mycommand", "props": ":{myoptions}"},
                {"selector": "mycommand2", "props": ":{myoptions2}"},
            ]
        )
        expected = dedent(
            """\
            \\begin{table}
            \\mycommand{myoptions}
            \\mycommand2{myoptions2}
            """
        )
        assert expected in self.s.to_latex()

    def test_cell_styling(self):
        self.s.highlight_max(props="emph:;Huge:-wrap-;")
        expected = dedent(
            """\
            \\begin{tabular}{lrrl}
             & A & B & C \\\\
            0 & 0 & \\emph{{\\Huge -0.61}} & ab \\\\
            1 & \\emph{{\\Huge 1}} & -1.22 & \\emph{{\\Huge cd}} \\\\
            \\end{tabular}
            """
        )
        assert expected == self.s.to_latex()

    def test_multiindex_columns(self):
        cidx = MultiIndex.from_tuples([("A", "a"), ("A", "b"), ("B", "c")])
        self.df.columns = cidx
        expected = dedent(
            """\
            \\begin{tabular}{lrrl}
             & \\multicolumn{2}{r}{A} & B \\\\
             & a & b & c \\\\
            0 & 0 & -0.61 & ab \\\\
            1 & 1 & -1.22 & cd \\\\
            \\end{tabular}
            """
        )
        assert expected == self.df.style.format(precision=2).to_latex()

    def test_multiindex_row(self):
        ridx = MultiIndex.from_tuples([("A", "a"), ("A", "b"), ("B", "c")])
        self.df.loc[2, :] = [2, -2.22, "de"]
        self.df = self.df.astype({"A": int})
        self.df.index = ridx
        expected = dedent(
            """\
            \\begin{tabular}{llrrl}
             &  & A & B & C \\\\
            \\multirow{2}{*}{A} & a & 0 & -0.61 & ab \\\\
             & b & 1 & -1.22 & cd \\\\
            B & c & 2 & -2.22 & de \\\\
            \\end{tabular}
            """
        )
        assert expected == self.df.style.format(precision=2).to_latex()

    def test_multiindex_row_and_col(self):
        cidx = MultiIndex.from_tuples([("Z", "a"), ("Z", "b"), ("Y", "c")])
        ridx = MultiIndex.from_tuples([("A", "a"), ("A", "b"), ("B", "c")])
        self.df.loc[2, :] = [2, -2.22, "de"]
        self.df = self.df.astype({"A": int})
        self.df.index, self.df.columns = ridx, cidx
        expected = dedent(
            """\
            \\begin{tabular}{llrrl}
             &  & \\multicolumn{2}{r}{Z} & Y \\\\
             &  & a & b & c \\\\
            \\multirow{2}{*}{A} & a & 0 & -0.61 & ab \\\\
             & b & 1 & -1.22 & cd \\\\
            B & c & 2 & -2.22 & de \\\\
            \\end{tabular}
            """
        )
        assert expected == self.df.style.format(precision=2).to_latex()
