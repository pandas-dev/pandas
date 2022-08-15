import pytest

from scripts.validate_exception_location import validate_exception_and_warning_placement

PATH = "t.py"
CUSTOM_EXCEPTION = "MyException"

TEST_CODE = """
import numpy as np
import sys

def my_func():
  pass

class {custom_name}({error_type}):
  pass

"""

testdata = [
    "Exception",
    "ValueError",
    "Warning",
    "UserWarning",
]


@pytest.mark.parametrize("error_type", testdata)
def test_class_that_inherits_an_exception_is_flagged(capsys, error_type):
    content = TEST_CODE.format(custom_name=CUSTOM_EXCEPTION, error_type=error_type)
    result_msg = (
        "t.py:8:0: {exception_name}: Please don't place exceptions or "
        "warnings outside of pandas/errors/__init__.py or "
        "pandas/_libs\n".format(exception_name=CUSTOM_EXCEPTION)
    )
    with pytest.raises(SystemExit, match=None):
        validate_exception_and_warning_placement(PATH, content)
    expected_msg, _ = capsys.readouterr()
    assert result_msg == expected_msg


def test_class_that_does_not_inherit_an_exception_is_flagged(capsys):
    content = "class MyClass(NonExceptionClass): pass"
    validate_exception_and_warning_placement(PATH, content)
