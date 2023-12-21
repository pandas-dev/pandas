import io

import pytest

from scripts import validate_unwanted_patterns


class TestBarePytestRaises:
    @pytest.mark.parametrize(
        "data",
        [
            (
                """
    with pytest.raises(ValueError, match="foo"):
        pass
    """
            ),
            (
                """
    # with pytest.raises(ValueError, match="foo"):
    #    pass
    """
            ),
            (
                """
    # with pytest.raises(ValueError):
    #    pass
    """
            ),
            (
                """
    with pytest.raises(
        ValueError,
        match="foo"
    ):
        pass
    """
            ),
        ],
    )
    def test_pytest_raises(self, data) -> None:
        fd = io.StringIO(data.strip())
        result = list(validate_unwanted_patterns.bare_pytest_raises(fd))
        assert result == []

    @pytest.mark.parametrize(
        "data, expected",
        [
            (
                (
                    """
    with pytest.raises(ValueError):
        pass
    """
                ),
                [
                    (
                        1,
                        (
                            "Bare pytests raise have been found. "
                            "Please pass in the argument 'match' "
                            "as well the exception."
                        ),
                    ),
                ],
            ),
            (
                (
                    """
    with pytest.raises(ValueError, match="foo"):
        with pytest.raises(ValueError):
            pass
        pass
    """
                ),
                [
                    (
                        2,
                        (
                            "Bare pytests raise have been found. "
                            "Please pass in the argument 'match' "
                            "as well the exception."
                        ),
                    ),
                ],
            ),
            (
                (
                    """
    with pytest.raises(ValueError):
        with pytest.raises(ValueError, match="foo"):
            pass
        pass
    """
                ),
                [
                    (
                        1,
                        (
                            "Bare pytests raise have been found. "
                            "Please pass in the argument 'match' "
                            "as well the exception."
                        ),
                    ),
                ],
            ),
            (
                (
                    """
    with pytest.raises(
        ValueError
    ):
        pass
    """
                ),
                [
                    (
                        1,
                        (
                            "Bare pytests raise have been found. "
                            "Please pass in the argument 'match' "
                            "as well the exception."
                        ),
                    ),
                ],
            ),
            (
                (
                    """
    with pytest.raises(
        ValueError,
        # match = "foo"
    ):
        pass
    """
                ),
                [
                    (
                        1,
                        (
                            "Bare pytests raise have been found. "
                            "Please pass in the argument 'match' "
                            "as well the exception."
                        ),
                    ),
                ],
            ),
        ],
    )
    def test_pytest_raises_raises(self, data, expected) -> None:
        fd = io.StringIO(data.strip())
        result = list(validate_unwanted_patterns.bare_pytest_raises(fd))
        assert result == expected


class TestStringsWithWrongPlacedWhitespace:
    @pytest.mark.parametrize(
        "data",
        [
            (
                """
    msg = (
        "foo\n"
        " bar"
    )
    """
            ),
            (
                """
    msg = (
        "foo"
        "  bar"
        "baz"
    )
    """
            ),
            (
                """
    msg = (
        f"foo"
        "  bar"
    )
    """
            ),
            (
                """
    msg = (
        "foo"
        f"  bar"
    )
    """
            ),
            (
                """
    msg = (
        "foo"
        rf"  bar"
    )
    """
            ),
        ],
    )
    def test_strings_with_wrong_placed_whitespace(self, data) -> None:
        fd = io.StringIO(data.strip())
        result = list(
            validate_unwanted_patterns.strings_with_wrong_placed_whitespace(fd)
        )
        assert result == []

    @pytest.mark.parametrize(
        "data, expected",
        [
            (
                (
                    """
    msg = (
        "foo"
        " bar"
    )
    """
                ),
                [
                    (
                        3,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    )
                ],
            ),
            (
                (
                    """
    msg = (
        f"foo"
        " bar"
    )
    """
                ),
                [
                    (
                        3,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    )
                ],
            ),
            (
                (
                    """
    msg = (
        "foo"
        f" bar"
    )
    """
                ),
                [
                    (
                        3,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    )
                ],
            ),
            (
                (
                    """
    msg = (
        f"foo"
        f" bar"
    )
    """
                ),
                [
                    (
                        3,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    )
                ],
            ),
            (
                (
                    """
    msg = (
        "foo"
        rf" bar"
        " baz"
    )
    """
                ),
                [
                    (
                        3,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    ),
                    (
                        4,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    ),
                ],
            ),
            (
                (
                    """
    msg = (
        "foo"
        " bar"
        rf" baz"
    )
    """
                ),
                [
                    (
                        3,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    ),
                    (
                        4,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    ),
                ],
            ),
            (
                (
                    """
    msg = (
        "foo"
        rf" bar"
        rf" baz"
    )
    """
                ),
                [
                    (
                        3,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    ),
                    (
                        4,
                        (
                            "String has a space at the beginning instead "
                            "of the end of the previous string."
                        ),
                    ),
                ],
            ),
        ],
    )
    def test_strings_with_wrong_placed_whitespace_raises(self, data, expected) -> None:
        fd = io.StringIO(data.strip())
        result = list(
            validate_unwanted_patterns.strings_with_wrong_placed_whitespace(fd)
        )
        assert result == expected


class TestNoDefaultUsedNotOnlyForTyping:
    @pytest.mark.parametrize(
        "data",
        [
            (
                """
def f(
    a: int | NoDefault,
    b: float | lib.NoDefault = 0.1,
    c: pandas._libs.lib.NoDefault = lib.no_default,
) -> lib.NoDefault | None:
    pass
"""
            ),
            (
                """
# var = lib.NoDefault
# the above is incorrect
a: NoDefault | int
b: lib.NoDefault = lib.no_default
"""
            ),
        ],
    )
    def test_nodefault_used_not_only_for_typing(self, data) -> None:
        fd = io.StringIO(data.strip())
        result = list(validate_unwanted_patterns.nodefault_used_not_only_for_typing(fd))
        assert result == []

    @pytest.mark.parametrize(
        "data, expected",
        [
            (
                (
                    """
def f(
    a = lib.NoDefault,
    b: Any
        = pandas._libs.lib.NoDefault,
):
    pass
"""
                ),
                [
                    (2, "NoDefault is used not only for typing"),
                    (4, "NoDefault is used not only for typing"),
                ],
            ),
            (
                (
                    """
a: Any = lib.NoDefault
if a is NoDefault:
    pass
"""
                ),
                [
                    (1, "NoDefault is used not only for typing"),
                    (2, "NoDefault is used not only for typing"),
                ],
            ),
        ],
    )
    def test_nodefault_used_not_only_for_typing_raises(self, data, expected) -> None:
        fd = io.StringIO(data.strip())
        result = list(validate_unwanted_patterns.nodefault_used_not_only_for_typing(fd))
        assert result == expected
