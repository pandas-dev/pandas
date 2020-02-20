import re

import pytest

import pandas as pd


@pytest.mark.filterwarnings(
    # openpyxl
    "ignore:defusedxml.lxml is no longer supported:DeprecationWarning"
)
@pytest.mark.filterwarnings(
    # html5lib
    "ignore:Using or importing the ABCs from:DeprecationWarning"
)
@pytest.mark.filterwarnings(
    # fastparquet
    "ignore:pandas.core.index is deprecated:FutureWarning"
)
@pytest.mark.filterwarnings(
    # pandas_datareader
    "ignore:pandas.util.testing is deprecated:FutureWarning"
)
def test_show_versions(capsys):
    # gh-32041
    pd.show_versions()
    captured = capsys.readouterr()
    result = captured.out

    # check header
    assert "INSTALLED VERSIONS" in result

    # check full commit hash
    if not re.search(r"commit\s*:\s[0-9a-f]{40}\n", result):
        # GH#32120  If test is being run in a branch that has uncommited
        #  changes, then we will not see the full commit hash, but this
        #  should show up in the pandas version number.
        assert re.search(r"pandas\s*: .*\.dirty\n", result)

    # check required dependency
    assert re.search(r"numpy\s*:\s([0-9\.\+a-f]|dev)+\n", result)

    # check optional dependency
    assert re.search(r"pyarrow\s*:\s([0-9\.]+|None)\n", result)
