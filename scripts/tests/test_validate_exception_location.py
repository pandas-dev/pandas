import pytest

from scripts.validate_exception_location import validate_exception_and_warning_placement

PATH = "t.py"

# ERRORS_IN_TESTING_RST is the set returned when parsing testing.rst for all the
# exceptions and warnings.
CUSTOM_EXCEPTION_NOT_IN_TESTING_RST = "MyException"
CUSTOM_EXCEPTION__IN_TESTING_RST = "MyOldException"
ERRORS_IN_TESTING_RST = {CUSTOM_EXCEPTION__IN_TESTING_RST}

TEST_CODE = """
import numpy as np
import sys

def my_func():
  pass

class {custom_name}({error_type}):
  pass

"""

# Test with various python-defined exceptions to ensure they are all flagged.
testdata = [
    "Exception",
    "ValueError",
    "Warning",
    "UserWarning",
]


@pytest.mark.parametrize("error_type", testdata)
def test_class_that_inherits_an_exception_and_is_not_in_the_testing_rst_is_flagged(
    capsys, error_type
):
    content = TEST_CODE.format(
        custom_name=CUSTOM_EXCEPTION_NOT_IN_TESTING_RST, error_type=error_type
    )
    result_msg = (
        "t.py:8:0: {exception_name}: Please don't place exceptions or "
        "warnings outside of pandas/errors/__init__.py or "
        "pandas/_libs\n".format(exception_name=CUSTOM_EXCEPTION_NOT_IN_TESTING_RST)
    )
    with pytest.raises(SystemExit, match=None):
        validate_exception_and_warning_placement(PATH, content, ERRORS_IN_TESTING_RST)
    expected_msg, _ = capsys.readouterr()
    assert result_msg == expected_msg


@pytest.mark.parametrize("error_type", testdata)
def test_class_that_inherits_an_exception_but_is_in_the_testing_rst_is_not_flagged(
    capsys, error_type
):
    content = TEST_CODE.format(
        custom_name=CUSTOM_EXCEPTION__IN_TESTING_RST, error_type=error_type
    )
    validate_exception_and_warning_placement(PATH, content, ERRORS_IN_TESTING_RST)


def test_class_that_does_not_inherit_an_exception_is_not_flagged(capsys):
    content = "class MyClass(NonExceptionClass): pass"
    validate_exception_and_warning_placement(PATH, content, ERRORS_IN_TESTING_RST)
