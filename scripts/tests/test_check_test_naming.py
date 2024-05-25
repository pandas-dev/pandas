import pytest
from scripts.check_test_naming import main

@pytest.mark.parametrize(
    "src, expected_out, expected_ret",
    [
        # Function without `test` prefix should be flagged
        (
            "def foo(): pass\n",
            "t.py:1:0 found test function which does not start with 'test'\n",
            1,
        ),
        # Class without `Test` prefix should be flagged even if it has correctly named test methods
        (
            "class Foo:\n    def test_foo(): pass\n",
            "t.py:1:0 found test class which does not start with 'Test'\n",
            1,
        ),
        # Correctly named test function should not be flagged
        ("def test_foo(): pass\n", "", 0),
        # Correctly named test class with test function should not be flagged
        ("class TestFoo:\n    def test_foo(): pass\n", "", 0),
        # Helper function used in a test function should not be flagged
        (
            "def foo():\n    pass\ndef test_foo():\n    foo()\n",
            "",
            0,
        ),
        # Class with pragma `# not a test` should not be flagged
        (
            "class Foo:  # not a test\n"
            "    pass\n"
            "def test_foo():\n"
            "    Class.foo()\n",
            "",
            0,
        ),
        # Pytest fixture without parentheses should not be flagged
        ("@pytest.fixture\ndef foo(): pass\n", "", 0),
        # Pytest fixture with parentheses should not be flagged
        ("@pytest.fixture()\ndef foo(): pass\n", "", 0),
        # Registered extension dtype class should not be flagged
        ("@register_extension_dtype\nclass Foo: pass\n", "", 0),
        # Nested classes where inner class is not a test should be flagged
        (
            "class TestOuter:\n    class Inner: pass\n",
            "t.py:2:4 found test class which does not start with 'Test'\n",
            1,
        ),
        # Function with multiple decorators including `@pytest.fixture` should not be flagged
        (
            "@pytest.fixture\n@another_decorator\ndef foo(): pass\n", "", 0
        ),
        # Mixed content with both correctly and incorrectly named elements
        (
            "def test_foo(): pass\nclass Foo: pass\n",
            "t.py:2:0 found test class which does not start with 'Test'\n",
            1,
        ),
    ],
)
def test_main(capsys, src, expected_out, expected_ret):
    """
    Test the main function from check_test_naming.py script.
    
    Parameters:
    - src: Source code to test.
    - expected_out: Expected output message.
    - expected_ret: Expected return code (0 for no errors, 1 for errors).
    
    This test uses pytest's capsys to capture printed output and checks
    both the output and return code against expected values.
    """
    ret = main(src, "t.py")
    out, _ = capsys.readouterr()
    assert out == expected_out
    assert ret == expected_ret
