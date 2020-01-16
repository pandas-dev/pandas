import io

# the "as" is just so I won't have to change this in many places later.
import validate_string_concatenation as validate_unwanted_patterns


class GoodBarePytestRaises:
    def normal_pytest_raises(self):
        return io.StringIO(
            """
with pytest.raises(ValueError, match="foo"):
    pass
"""
        )

    def pytest_raises_as_comment(self):
        return io.StringIO(
            """
# with pytest.raises(ValueError, match="foo"):
#     pass
"""
        )

    def bare_pytest_raises_as_comment(self):
        """
        We do not care if the bare pytest raises is commented out.
        """
        return io.StringIO(
            """
# with pytest.raises(ValueError):
#     pass
"""
        )


class BadBarePytestRaises:
    def normal_bare_pytest_raises(self):
        return io.StringIO(
            """
with pytest.raises(ValueError):
    raise ValueError("foo")
"""
        )


class TestGoodBarePytestRaises:
    def test_good_normal_pytest_raises(self):
        result = validate_unwanted_patterns.bare_pytest_raises(
            GoodBarePytestRaises.normal_pytest_raises(self)
        )

        for _ in result:
            assert False  # Not sure about this

    def test_pytest_raises_as_comment(self):
        result = validate_unwanted_patterns.bare_pytest_raises(
            GoodBarePytestRaises.pytest_raises_as_comment(self)
        )
        for _ in result:
            assert False  # Not sure about this

    def test_bare_pytest_raises_as_comment(self):
        result = validate_unwanted_patterns.bare_pytest_raises(
            GoodBarePytestRaises.bare_pytest_raises_as_comment(self)
        )
        for _ in result:
            assert False  # Not sure about this


class TestBadBarePytestRaises:
    def test_bad_normal_pytest_raises(self):
        result = validate_unwanted_patterns.bare_pytest_raises(
            BadBarePytestRaises.normal_bare_pytest_raises(self)
        )

        for bare_pytest_raise in result:
            assert bare_pytest_raise[0] == 2
