import pytest

from scripts.use_io_common_urlopen import use_io_common_urlopen

PATH = "t.py"


def test_inconsistent_usage(capsys):
    content = "from urllib.request import urlopen"
    result_msg = (
        "t.py:1:0: Don't use urllib.request.urlopen, "
        "use pandas.io.common.urlopen instead\n"
    )
    with pytest.raises(SystemExit, match=None):
        use_io_common_urlopen(content, PATH)
    expected_msg, _ = capsys.readouterr()
    assert result_msg == expected_msg


def test_consistent_usage():
    # should not raise
    content = "from pandas.io.common import urlopen"
    use_io_common_urlopen(content, PATH)
