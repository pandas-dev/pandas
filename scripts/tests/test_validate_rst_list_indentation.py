import pytest

from scripts.validate_rst_list_indentation import find_indented_lists


@pytest.mark.parametrize(
    "content",
    [
        "Text:\n\n* item\n* item\n",
        "* item\n\n  * nested\n  * nested\n\n* item\n",
        "#. item\n\n   #. nested\n",
        # a quote that happens to contain a list alongside other content
        "Text:\n\n    A quote.\n\n    * item\n",
        # lists inside directives unknown to docutils are not parsed
        ".. ipython:: python\n\n    * not a list\n",
    ],
)
def test_find_indented_lists_valid(content):
    assert find_indented_lists(content) == []


@pytest.mark.parametrize(
    "content, expected",
    [
        ("Text:\n\n  * item\n  * item\n", [3]),
        ("* item\n\n    * nested\n    * nested\n\n* item\n", [3]),
        ("Text:\n\n  1. item\n  2. item\n", [3]),
        (".. note::\n\n   Text:\n\n     - item\n", [5]),
    ],
)
def test_find_indented_lists_invalid(content, expected):
    assert find_indented_lists(content) == expected
