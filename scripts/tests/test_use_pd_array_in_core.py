import pytest

from scripts.use_pd_array_in_core import use_pd_array

BAD_FILE_0 = "import pandas as pd\npd.array"
BAD_FILE_1 = "\nfrom pandas import array"
GOOD_FILE_0 = "from pandas import array as pd_array"
GOOD_FILE_1 = "from pandas.core.construction import array as pd_array"
PATH = "t.py"


@pytest.mark.parametrize("content", [BAD_FILE_0, BAD_FILE_1])
def test_inconsistent_usage(content, capsys):
    result_msg = (
        "t.py:2:0: Don't use pd.array in core, import array as pd_array instead\n"
    )
    with pytest.raises(SystemExit, match=None):
        use_pd_array(content, PATH)
    expected_msg, _ = capsys.readouterr()
    assert result_msg == expected_msg


@pytest.mark.parametrize("content", [GOOD_FILE_0, GOOD_FILE_1])
def test_consistent_usage(content):
    # should not raise
    use_pd_array(content, PATH)
