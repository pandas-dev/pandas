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
@pytest.mark.filterwarnings(
    # https://github.com/pandas-dev/pandas/issues/35252
    "ignore:Distutils:UserWarning"
)
@pytest.mark.filterwarnings("ignore:Setuptools is replacing distutils:UserWarning")
def test_show_versions(capsys):
    # gh-32041
    pd.show_versions()
    captured = capsys.readouterr()
    result = captured.out

    # check header
    assert "INSTALLED VERSIONS" in result

    # check full commit hash
    assert re.search(r"commit\s*:\s[0-9a-f]{40}\n", result)

    # check required dependency
    assert re.search(r"numpy\s*:\s([0-9\.\+a-f]|dev)+\n", result)

    # check optional dependency
    assert re.search(r"pyarrow\s*:\s([0-9\.]+|None)\n", result)
