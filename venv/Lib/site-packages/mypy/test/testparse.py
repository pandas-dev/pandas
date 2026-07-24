"""Tests for the mypy parser."""

from __future__ import annotations

import sys

from pytest import skip

from mypy import defaults
from mypy.config_parser import parse_mypy_comments
from mypy.errors import CompileError, Errors
from mypy.options import Options
from mypy.parse import parse
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.test.helpers import assert_string_arrays_equal, find_test_files, parse_options
from mypy.util import get_mypy_comments


class ParserSuite(DataSuite):
    required_out_section = True
    base_path = "."
    files = find_test_files(pattern="parse*.test", exclude=["parse-errors.test"])

    if sys.version_info < (3, 12):
        files.remove("parse-python312.test")
    if sys.version_info < (3, 13):
        files.remove("parse-python313.test")
    if sys.version_info < (3, 14):
        files.remove("parse-python314.test")

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        test_parser(testcase)


def test_parser(testcase: DataDrivenTestCase) -> None:
    """Perform a single parser test case.

    The argument contains the description of the test case.
    """
    options = Options()
    options.hide_error_codes = True

    if testcase.file.endswith("python310.test"):
        options.python_version = (3, 10)
    elif testcase.file.endswith("python312.test"):
        options.python_version = (3, 12)
    elif testcase.file.endswith("python313.test"):
        options.python_version = (3, 13)
    elif testcase.file.endswith("python314.test"):
        options.python_version = (3, 14)
    else:
        options.python_version = defaults.PYTHON3_VERSION

    source = "\n".join(testcase.input)

    # Apply mypy: comments to options.
    comments = get_mypy_comments(source)
    changes, _ = parse_mypy_comments(comments, options)
    options = options.apply_changes(changes)

    try:
        errors = Errors(options)
        n = parse(
            bytes(source, "ascii"),
            fnam="main",
            module="__main__",
            errors=errors,
            options=options,
            file_exists=False,
            eager=True,
        )
        if errors.is_errors():
            errors.raise_error()
        a = n.str_with_options(options).split("\n")
    except CompileError as e:
        a = e.messages
    assert_string_arrays_equal(
        testcase.output, a, f"Invalid parser output ({testcase.file}, line {testcase.line})"
    )


# The file name shown in test case output. This is displayed in error
# messages, and must match the file name in the test case descriptions.
INPUT_FILE_NAME = "file"


class ParseErrorSuite(DataSuite):
    required_out_section = True
    base_path = "."
    files = ["parse-errors.test"]

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        test_parse_error(testcase)


def test_parse_error(testcase: DataDrivenTestCase) -> None:
    try:
        options = parse_options("\n".join(testcase.input), testcase, 0)
        if options.python_version < defaults.PYTHON3_VERSION:
            options.python_version = defaults.PYTHON3_VERSION
        if options.python_version != sys.version_info[:2]:
            skip()
        # Compile temporary file. The test file contains non-ASCII characters.
        errors = Errors(options)
        parse(
            bytes("\n".join(testcase.input), "utf-8"),
            INPUT_FILE_NAME,
            "__main__",
            errors=errors,
            options=options,
            file_exists=False,
            eager=True,
        )
        if errors.is_errors():
            errors.raise_error()
        raise AssertionError("No errors reported")
    except CompileError as e:
        if e.module_with_blocker is not None:
            assert e.module_with_blocker == "__main__"
        # Verify that there was a compile error and that the error messages
        # are equivalent.
        assert_string_arrays_equal(
            testcase.output,
            e.messages,
            f"Invalid compiler output ({testcase.file}, line {testcase.line})",
        )
