import json
import os
import re

import pytest

from pandas.util._print_versions import (
    _get_dependency_info,
    _get_sys_info,
)

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
def test_show_versions(tmpdir):
    as_json = os.path.join(tmpdir, "test_output.json")

    pd.show_versions(as_json=as_json)

    # make sure that the file was created
    assert os.path.exists(as_json)

    with open(as_json) as fd:
        contents = fd.readlines()
        str_contents = "".join(contents)

        # make sure that there was output to the file
        assert str_contents

        # check if file output is valid JSON, will raise an exception if not
        dict_check = json.loads(str_contents)

    # Basic check that each version element is found in output
    version_elements = {
        "system": _get_sys_info(),
        "dependencies": _get_dependency_info(),
    }

    assert version_elements == dict_check


def test_show_versions_console_json(capsys):
    pd.show_versions(as_json=True)
    captured = capsys.readouterr()
    result = captured.out

    # check valid json is printed to the console if as_json is True
    dict_check = json.loads(result)

    # Basic check that each version element is found in output
    version_elements = {
        "system": _get_sys_info(),
        "dependencies": _get_dependency_info(),
    }

    assert version_elements == dict_check


def test_show_versions_console(capsys):
    # gh-32041
    pd.show_versions(as_json=False)
    captured = capsys.readouterr()
    result = captured.out

    assert "INSTALLED VERSIONS" in result

    # check full commit hash
    assert re.search(r"commit\s*:\s[0-9a-f]{40}\n", result)

    # check required dependency
    # 2020-12-09 npdev has "dirty" in the tag
    assert re.search(r"numpy\s*:\s([0-9\.\+a-g\_]|dev)+(dirty)?\n", result)

    # check optional dependency
    assert re.search(r"pyarrow\s*:\s([0-9\.]+|None)\n", result)


def test_json_output_match(capsys, tmpdir):
    pd.show_versions(as_json=True)
    result_console = capsys.readouterr().out

    out_path = os.path.join(tmpdir, "test_json.json")
    pd.show_versions(as_json=out_path)
    with open(out_path) as out_fd:
        result_file = "".join(out_fd.readlines())

    assert result_console == result_file
