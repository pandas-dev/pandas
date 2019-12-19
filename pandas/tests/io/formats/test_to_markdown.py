import pandas.util._test_decorators as td

import pandas as pd


@td.skip_if_no_tabulate
def test_to_markdown():
    df = pd.DataFrame([1, 2, 3])
    result = df.to_markdown()
    assert (
        result == "|    |   0 |\n|---:|----:|\n|  0 |   1 |\n|  1 |   2 |\n|  2 |   3 |"
    )
