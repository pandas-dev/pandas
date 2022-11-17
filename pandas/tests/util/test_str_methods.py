import sys

import pytest

if sys.version_info < (3, 9):
    from pandas.util._str_methods import (
        removeprefix,
        removesuffix,
    )

    @pytest.mark.parametrize(
        "string, prefix, expected",
        (
            ("wildcat", "wild", "cat"),
            ("blackbird", "black", "bird"),
            ("housefly", "house", "fly"),
            ("ladybug", "lady", "bug"),
            ("rattlesnake", "rattle", "snake"),
            ("baboon", "badger", "baboon"),
            ("quetzal", "elk", "quetzal"),
        ),
    )
    def test_remove_prefix(string, prefix, expected):
        result = removeprefix(string, prefix)
        assert result == expected

    @pytest.mark.parametrize(
        "string, suffix, expected",
        (
            ("wildcat", "cat", "wild"),
            ("blackbird", "bird", "black"),
            ("housefly", "fly", "house"),
            ("ladybug", "bug", "lady"),
            ("rattlesnake", "snake", "rattle"),
            ("seahorse", "horse", "sea"),
            ("baboon", "badger", "baboon"),
            ("quetzal", "elk", "quetzal"),
        ),
    )
    def test_remove_suffix(string, suffix, expected):
        result = removesuffix(string, suffix)
        assert result == expected

else:
    # NOTE: remove this file when pyupgrade --py39-plus removes
    # the above block
    pass
