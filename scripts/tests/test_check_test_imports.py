import pytest

from scripts.check_test_imports import main


@pytest.mark.parametrize(
    "src, expected_out, expected_ret",
    [
        (
            "from pandas import Series\n",
            "t.py:1:0 use 'pd.Series' rather than importing 'Series' "
            "from the pandas namespace\n",
            1,
        ),
        (
            "from pandas import (\n    DataFrame,\n    Series,\n)\n",
            "t.py:1:0 use 'pd.DataFrame' rather than importing 'DataFrame' "
            "from the pandas namespace\n"
            "t.py:1:0 use 'pd.Series' rather than importing 'Series' "
            "from the pandas namespace\n",
            1,
        ),
        (
            "from pandas import compat\n",
            "t.py:1:0 import 'compat' from the module which defines it rather "
            "than from the pandas namespace\n",
            1,
        ),
        (
            "from pandas import _testing as tm\n",
            "t.py:1:0 import '_testing' from the module which defines it rather "
            "than from the pandas namespace\n",
            1,
        ),
        (
            "import pandas\n",
            "t.py:1:0 'pandas' should be imported as 'pd'\n",
            1,
        ),
        ("import pandas as pd\n", "", 0),
        ("import pandas._testing as tm\n", "", 0),
        ("import pandas.util._test_decorators as td\n", "", 0),
        ("from pandas.util import _test_decorators as td\n", "", 0),
        ("from pandas.core.dtypes.common import is_integer_dtype\n", "", 0),
        ("from pandas.tseries import offsets\n", "", 0),
        ("from .common import fixture\n", "", 0),
    ],
)
def test_main(capsys, src, expected_out, expected_ret):
    ret = main(src, "t.py")
    out, _ = capsys.readouterr()
    assert out == expected_out
    assert ret == expected_ret
