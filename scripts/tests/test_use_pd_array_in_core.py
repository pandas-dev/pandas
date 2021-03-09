import pytest

from scripts.use_pd_array_in_core import use_pd_array

BAD_FILE_0 = "\nfrom pandas import array as pd_array"
BAD_FILE_1 = "import pandas as pd\npd.array"
GOOD_FILE = "from pandas.core.construction import pd_array"
PATH = "t.py"


@pytest.mark.parametrize("content", [BAD_FILE_0, BAD_FILE_1])
def test_inconsistent_usage(content, capsys):
    result_msg = (
        "t.py:2:0: Don't use pd.array in core, "
        "instead use 'from pandas.core.construction import pd_array'\n"
    )
    with pytest.raises(SystemExit):
        use_pd_array(content, PATH)
    expected_msg, _ = capsys.readouterr()
    assert result_msg == expected_msg


def test_consistent_usage():
    # should not raise
    use_pd_array(GOOD_FILE, PATH)
