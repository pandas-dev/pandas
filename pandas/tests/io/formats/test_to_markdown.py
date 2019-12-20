import pandas.util._test_decorators as td

import pandas as pd


@td.skip_if_no_tabulate
class TestToMarkdown:
    def test_to_markdown(self):
        df = pd.DataFrame([1, 2, 3])
        result = df.to_markdown()
        assert (
            result
            == "|    |   0 |\n|---:|----:|\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
        )

    def test_to_markdown_other_tablefmt(self):
        df = pd.DataFrame([1, 2, 3])
        result = df.to_markdown(tablefmt="jira")
        assert result == "||    ||   0 ||\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
