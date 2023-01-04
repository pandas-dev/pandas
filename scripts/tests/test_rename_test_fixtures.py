import pytest

from scripts.rename_test_fixtures import main


@pytest.mark.parametrize(
    "src, expected",
    [
        pytest.param(
            "@pytest.fixture\ndef foo():\n    pass",
            '@pytest.fixture(name="foo")\ndef fixture_foo():\n    pass',
            id="attribute",
        ),
        pytest.param(
            "@pytest.fixture()\ndef foo():\n    pass",
            '@pytest.fixture(name="foo")\ndef fixture_foo():\n    pass',
            id="call with no arguments",
        ),
        pytest.param(
            "@pytest.fixture(autouse=True)\ndef foo():\n    pass",
            '@pytest.fixture(name="foo", autouse=True)\n'
            "def fixture_foo():\n"
            "    pass",
            id="call with arguments",
        ),
        pytest.param(
            '@pytest.fixture(name="foo")\ndef foo():\n    pass',
            '@pytest.fixture(name="foo")\ndef fixture_foo():\n    pass',
            id="already renamed",
        ),
        pytest.param(
            '@pytest.fixture(name="foo")\ndef fixture_foo():\n    pass',
            None,
            id="already renamed",
        ),
    ],
)
def test_rename_test_fixtures(src, expected):
    result = main(src, "t.py")
    assert result == expected
