import pytest

from scripts.check_for_inconsistent_pandas_namespace import (
    check_for_inconsistent_pandas_namespace,
)

BAD_FILE_0 = "cat_0 = Categorical()\ncat_1 = pd.Categorical()"
BAD_FILE_1 = "cat_0 = pd.Categorical()\ncat_1 = Categorical()"
GOOD_FILE_0 = "cat_0 = Categorical()\ncat_1 = Categorical()"
GOOD_FILE_1 = "cat_0 = pd.Categorical()\ncat_1 = pd.Categorical()"
PATH = "t.py"


@pytest.mark.parametrize("content", [BAD_FILE_0, BAD_FILE_1])
def test_inconsistent_usage(content):
    msg = r"Found both `pd\.Categorical` and `Categorical` in t\.py"
    with pytest.raises(RuntimeError, match=msg):
        check_for_inconsistent_pandas_namespace(content, PATH, replace=False)


@pytest.mark.parametrize("content", [GOOD_FILE_0, GOOD_FILE_1])
def test_consistent_usage(content):
    # should not raise
    check_for_inconsistent_pandas_namespace(content, PATH, replace=False)


@pytest.mark.parametrize("content", [BAD_FILE_0, BAD_FILE_1])
def test_inconsistent_usage_with_replace(content):
    result = check_for_inconsistent_pandas_namespace(content, PATH, replace=True)
    expected = "cat_0 = Categorical()\ncat_1 = Categorical()"
    assert result == expected


@pytest.mark.parametrize("content", [GOOD_FILE_0, GOOD_FILE_1])
def test_consistent_usage_with_replace(content):
    result = check_for_inconsistent_pandas_namespace(content, PATH, replace=True)
    expected = content
    assert result == expected
