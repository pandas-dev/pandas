import numpy as np
import pytest

from pandas import DataFrame

pytest.importorskip("jinja2")
from pandas.io.formats.style import Styler


class TestStylerTooltip:
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
        s = Styler(df).render()  # no set_tooltips() creates no <span>
        assert '<style type="text/css">\n</style>' in s
        assert '<span class="pd-t"></span>' not in s

    def test_tooltip_css_class(self):
        # GH 21266
        df = DataFrame(data=[[0, 1], [2, 3]])
        s = (
            Styler(df, uuid_len=0)
            .set_tooltips(
                DataFrame([["tooltip"]]),
                css_class="other-class",
                props=[("color", "green")],
            )
            .render()
        )
        assert "#T__ .other-class {\n  color: green;\n" in s
        assert '#T__ #T__row0_col0 .other-class::after {\n  content: "tooltip";\n' in s

        # GH 39563
        s = (
            Styler(df, uuid_len=0)
            .set_tooltips(
                DataFrame([["tooltip"]]),
                css_class="other-class",
                props="color:green;color:red;",
            )
            .render()
        )
        assert "#T__ .other-class {\n  color: green;\n  color: red;\n}" in s
