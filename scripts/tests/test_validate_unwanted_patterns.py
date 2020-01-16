import io

# the "as" is just so I won't have to change this in many places later.
import validate_string_concatenation as validate_unwanted_patterns


class TestBarePytestRaises:
    def test_bare_pytest_raises(self):
        fd = io.StringIO(
            """
with pytest.raises(ValueError):
    pass
""".strip()
        )
        result = list(validate_unwanted_patterns.bare_pytest_raises(fd))
        expected = [
            (
                1,
                (
                    "Bare pytests raise have been found. "
                    "Please pass in the argument 'match' as well the exception."
                ),
            )
        ]
        assert result == expected

    def test_pytest_raises(self):
        fd = io.StringIO(
            """
with pytest.raises(ValueError, match="foo"):
    pass
""".strip()
        )
        result = list(validate_unwanted_patterns.bare_pytest_raises(fd))
        expected = []
        assert result == expected

    def test_bare_pytest_raises_as_comment(self):
        fd = io.StringIO(
            """
# with pytest.raises(ValueError):
#    pass
""".strip()
        )
        result = list(validate_unwanted_patterns.bare_pytest_raises(fd))
        expected = []
        assert result == expected

    def test_pytest_raises_as_comment(self):
        fd = io.StringIO(
            """
# with pytest.raises(ValueError, match="foo"):
#    pass
""".strip()
        )
        result = list(validate_unwanted_patterns.bare_pytest_raises(fd))
        expected = []
        assert result == expected


class TestStringsToConcatenate:
    def test_strings_to_concatenate(self):
        fd = io.StringIO(
            """
msg = ("bar " "baz")
""".strip()
        )
        result = list(validate_unwanted_patterns.strings_to_concatenate(fd))
        expected = [
            (
                1,
                (
                    "String unnecessarily split in two by black. "
                    "Please merge them manually."
                ),
            )
        ]
        assert result == expected

    def test_strings_to_concatenate_two_occurrences_same_line(self):
        fd = io.StringIO(
            """
msg = ("foo " "bar " "baz")
""".strip()
        )
        result = list(validate_unwanted_patterns.strings_to_concatenate(fd))
        expected = [
            (
                1,
                (
                    "String unnecessarily split in two by black. "
                    "Please merge them manually."
                ),
            ),
            (
                1,
                (
                    "String unnecessarily split in two by black. "
                    "Please merge them manually."
                ),
            ),
        ]
        assert result == expected


class TestStringsWithWrongPlacedWhitespace:
    def test_strings_with_wrong_placed_whitespace(self):
        fd = io.StringIO(
            """
msg = (
    "foo"
    " bar"
)
""".strip()
        )
        result = list(
            validate_unwanted_patterns.strings_with_wrong_placed_whitespace(fd)
        )
        expected = [
            (
                3,
                (
                    "String has a space at the beginning instead "
                    "of the end of the previous string."
                ),
            ),
        ]
        assert result == expected

    def test_strings_with_wrong_placed_whitespace_ends_new_line(self):
        fd = io.StringIO(
            """
msg = (
    "foo\n"
    " bar"
)
""".strip()
        )
        result = list(
            validate_unwanted_patterns.strings_with_wrong_placed_whitespace(fd)
        )
        expected = []
        assert result == expected

    def test_strings_with_wrong_placed_whitespace_first_line_starts_double_sapce(self):
        fd = io.StringIO(
            """
msg = (
    "  foo"
    "bar"
)
""".strip()
        )
        result = list(
            validate_unwanted_patterns.strings_with_wrong_placed_whitespace(fd)
        )
        expected = []
        assert result == expected

    def test_strings_with_wrong_placed_whitespace_second_line_starts_double_sapce(self):
        fd = io.StringIO(
            """
msg = (
    "foo"
    "  bar"
)
""".strip()
        )
        result = list(
            validate_unwanted_patterns.strings_with_wrong_placed_whitespace(fd)
        )
        expected = []
        assert result == expected

    def test_strings_with_wrong_placed_whitespace_first_line_string_literal_prefix(
        self,
    ):
        fd = io.StringIO(
            """
msg = (
    f"foo"
    "  bar"
)
""".strip()
        )
        result = list(
            validate_unwanted_patterns.strings_with_wrong_placed_whitespace(fd)
        )
        expected = []
        assert result == expected

    def test_strings_with_wrong_placed_whitespace_second_line_string_literal_prefix(
        self,
    ):
        fd = io.StringIO(
            """
msg = (
    "foo"
    f" bar"
)
""".strip()
        )
        result = list(
            validate_unwanted_patterns.strings_with_wrong_placed_whitespace(fd)
        )
        expected = [
            (
                3,
                (
                    "String has a space at the beginning instead "
                    "of the end of the previous string."
                ),
            )
        ]
        assert result == expected
