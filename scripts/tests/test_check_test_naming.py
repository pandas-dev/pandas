import pytest

from scripts.check_test_naming import main


@pytest.mark.parametrize(
    "src, expected_out, expected_ret",
    [
        (
            "def foo(): pass\n",
            "t.py:1:0 found test function which does not start with 'test'\n",
            1,
        ),
        (
            "class Foo:\n    def test_foo(): pass\n",
            "t.py:1:0 found test class which does not start with 'Test'\n",
            1,
        ),
        ("def test_foo(): pass\n", "", 0),
        ("class TestFoo:\n    def test_foo(): pass\n", "", 0),
        (
            "def foo():\n    pass\ndef test_foo():\n    foo()\n",
            "",
            0,
        ),
        (
            "class Foo:  # not a test\n"
            "    pass\n"
            "def test_foo():\n"
            "    Class.foo()\n",
            "",
            0,
        ),
        ("@pytest.fixture\ndef foo(): pass\n", "", 0),
        ("@pytest.fixture()\ndef foo(): pass\n", "", 0),
        ("@register_extension_dtype\nclass Foo: pass\n", "", 0),
    ],
)
def test_main(capsys, src, expected_out, expected_ret):
    ret = main(src, "t.py")
    out, _ = capsys.readouterr()
    assert out == expected_out
    assert ret == expected_ret
