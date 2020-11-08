from pathlib import Path

import pytest

from scripts.check_for_inconsistent_pandas_namespace import main

BAD_FILE_0 = "cat_0 = Categorical()\ncat_1 = pd.Categorical()"
BAD_FILE_1 = "cat_0 = pd.Categorical()\ncat_1 = Categorical()"
GOOD_FILE_0 = "cat_0 = Categorical()\ncat_1 = Categorical()"
GOOD_FILE_1 = "cat_0 = pd.Categorical()\ncat_1 = pd.Categorical()"


@pytest.mark.parametrize("content", [BAD_FILE_0, BAD_FILE_1])
def test_inconsistent_usage(tmpdir, content):
    tmpfile = Path(tmpdir / "tmpfile.py")
    tmpfile.touch()
    tmpfile.write_text(content)
    msg = fr"Found both `pd\.Categorical` and `Categorical` in {str(tmpfile)}"
    with pytest.raises(AssertionError, match=msg):
        main((str(tmpfile),))


@pytest.mark.parametrize("content", [GOOD_FILE_0, GOOD_FILE_1])
def test_consistent_usage(tmpdir, content):
    tmpfile = Path(tmpdir / "tmpfile.py")
    tmpfile.touch()
    tmpfile.write_text(content)
    main((str(tmpfile),))  # Should not raise.
