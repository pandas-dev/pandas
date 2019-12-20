import pandas.util._test_decorators as td

import pandas as pd


@td.skip_if_no_tabulate
class TestToMarkdown:
    def test_simple(self):
        df = pd.DataFrame([1, 2, 3])
        result = df.to_markdown()
        assert (
            result
            == "|    |   0 |\n|---:|----:|\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
        )

    def test_other_tablefmt(self):
        df = pd.DataFrame([1, 2, 3])
        result = df.to_markdown(tablefmt="jira")
        assert result == "||    ||   0 ||\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"

    def test_other_headers(self):
        df = pd.DataFrame([1, 2, 3])
        result = df.to_markdown(headers=["foo", "bar"])
        assert (
            result
            == "|   foo |   bar |\n|------:|------:|\n|     0 |     1 |\n|     1 |     2 |\n|     2 |     3 |"
        )

    def test_series(self):
        s = pd.Series([1, 2, 3], name="foo")
        result = s.to_markdown()
        assert (
            result
            == "|    |   foo |\n|---:|------:|\n|  0 |     1 |\n|  1 |     2 |\n|  2 |     3 |"
        )
