import re

import pandas as pd


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
