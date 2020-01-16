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
    pass


class TestStringsWithWrongPlacedWhitespace:
    pass
