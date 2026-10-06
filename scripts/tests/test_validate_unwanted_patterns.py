import io

import pytest

from scripts import validate_unwanted_patterns


class TestStringsWithWrongPlacedWhitespace:
    @pytest.mark.parametrize(
        "data",
        [
            (
                r"""
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
            pytest.param(
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
                marks=pytest.mark.xfail(reason="TODO: Changed with PY3.13+"),
            ),
            pytest.param(
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
                marks=pytest.mark.xfail(reason="TODO: Changed with PY3.13+"),
            ),
            pytest.param(
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
                marks=pytest.mark.xfail(reason="TODO: Changed with PY3.13+"),
            ),
            pytest.param(
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
                marks=pytest.mark.xfail(reason="TODO: Changed with PY3.13+"),
            ),
            pytest.param(
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
                marks=pytest.mark.xfail(reason="TODO: Changed with PY3.13+"),
            ),
            pytest.param(
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
                marks=pytest.mark.xfail(reason="TODO: Changed with PY3.13+"),
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


@pytest.mark.parametrize("function", ["warnings.warn", "warn"])
@pytest.mark.parametrize("positional", [True, False])
@pytest.mark.parametrize(
    "category",
    [
        "FutureWarning",
        "DeprecationWarning",
        "PendingDeprecationWarning",
        "Pandas4Warning",
        "RuntimeWarning",
    ],
)
@pytest.mark.parametrize("pdlint_ignore", [True, False])
def test_doesnt_use_pandas_warnings(function, positional, category, pdlint_ignore):
    code = (
        f"{function}({'  # pdlint: ignore[warning_class]' if pdlint_ignore else ''}\n"
        f'   "message",\n'
        f"   {'' if positional else 'category='}{category},\n"
        f")\n"
    )
    flag_issue = (
        category in ["FutureWarning", "DeprecationWarning", "PendingDeprecationWarning"]
        and not pdlint_ignore
    )
    fd = io.StringIO(code)
    result = list(validate_unwanted_patterns.doesnt_use_pandas_warnings(fd))
    if flag_issue:
        assert len(result) == 1
        assert result[0][0] == 1
        assert result[0][1].startswith(f"Don't use {category}")
    else:
        assert len(result) == 0


@pytest.mark.parametrize(
    "data",
    [
        # already using the explicit join form
        'msg = "|".join(["first message", "second message"])',
        'pytest.raises(ValueError, match="|".join(["foo", "bar"]))',
        # inline alternation correctly wrapped in a group
        'msg = "cannot (add|subtract)"',
        'msg = r"index (0|1) out of bounds"',
        # f-string with a grouped inline alternation
        'msg = f"cannot cast (object|str) to {dtype}"',
        # an escaped pipe is a literal "|", not an alternation
        r'msg = r"value must be one of None\|12"',
        # pipe inside a character class
        'msg = "[<|>] not supported"',
        # match= referring to a name rather than a literal
        "pytest.raises(ValueError, match=msg)",
        # no alternation at all
        'msg = "a single message"',
    ],
)
def test_bare_pipe_alternation_in_message(data) -> None:
    fd = io.StringIO(data.strip())
    result = list(validate_unwanted_patterns.bare_pipe_alternation_in_message(fd))
    assert result == []


@pytest.mark.parametrize(
    "data, expected",
    [
        (
            'msg = "first message|second message"',
            [(1, validate_unwanted_patterns.BARE_PIPE_MESSAGE)],
        ),
        (
            'pat = "foo|bar|baz"',
            [(1, validate_unwanted_patterns.BARE_PIPE_MESSAGE)],
        ),
        (
            'pytest.raises(ValueError, match="foo|bar")',
            [(1, validate_unwanted_patterns.BARE_PIPE_MESSAGE)],
        ),
        (
            # a raw string with escaped parens still has a bare top-level pipe
            r'msg = r"foo\(s\)|bar"',
            [(1, validate_unwanted_patterns.BARE_PIPE_MESSAGE)],
        ),
        (
            # substitutions on each side of a bare pipe in an f-string
            'msg = f"{first}|{second}"',
            [(1, validate_unwanted_patterns.BARE_PIPE_MESSAGE)],
        ),
        (
            # implicit concatenation across lines
            """
msg = (
    "first message|"
    "second message"
)
""",
            [(2, validate_unwanted_patterns.BARE_PIPE_MESSAGE)],
        ),
    ],
)
def test_bare_pipe_alternation_in_message_raises(data, expected) -> None:
    fd = io.StringIO(data.strip())
    result = list(validate_unwanted_patterns.bare_pipe_alternation_in_message(fd))
    assert result == expected
