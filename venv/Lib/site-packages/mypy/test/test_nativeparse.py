"""Tests for the experimental native mypy parser.

To run these, you will need to manually install ast_serialize from
https://github.com/mypyc/ast_serialize first (see the README for the details).
"""

from __future__ import annotations

import contextlib
import os
import tempfile
import unittest
from collections.abc import Iterator

from librt.internal import ReadBuffer

from mypy import defaults, nodes
from mypy.cache import (
    END_TAG,
    LIST_GEN,
    LIST_INT,
    LITERAL_INT,
    LITERAL_NONE,
    LITERAL_STR,
    LOCATION,
    read_int,
)
from mypy.config_parser import parse_mypy_comments
from mypy.errors import CompileError
from mypy.nodes import MypyFile, ParseError
from mypy.options import Options
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.test.helpers import assert_string_arrays_equal
from mypy.util import get_mypy_comments

# If the experimental ast_serialize module isn't installed, the following import will fail
# and we won't run any native parser tests.
try:
    from mypy.nativeparse import (
        State,
        deserialize_imports,
        native_parse,
        parse_to_binary_ast,
        read_statements,
    )

    has_nativeparse = True
except ImportError:
    has_nativeparse = False


class NativeParserSuite(DataSuite):
    required_out_section = True
    base_path = "."
    files = ["native-parser.test"] if has_nativeparse else []

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        test_parser(testcase)


class NativeParserImportsSuite(DataSuite):
    required_out_section = True
    base_path = "."
    files = ["native-parser-imports.test"] if has_nativeparse else []

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        test_parser_imports(testcase)


def test_parser(testcase: DataDrivenTestCase) -> None:
    """Perform a single native parser test case.

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

    # Check if we should skip function bodies (when ignoring errors)
    skip_function_bodies = "# mypy: ignore-errors=True" in source

    try:
        with temp_source(source) as fnam:
            node, errors, type_ignores = native_parse(fnam, options, skip_function_bodies)
            errors += load_tree(node, options)
            node.path = "main"
            a = node.str_with_options(options).split("\n")
            a = [format_error(err) for err in errors] + a
            a = [format_ignore(ignore) for ignore in type_ignores] + a
    except CompileError as e:
        a = e.messages
    assert_string_arrays_equal(
        testcase.output, a, f"Invalid parser output ({testcase.file}, line {testcase.line})"
    )


def format_error(err: ParseError) -> str:
    return f"{err['line']}:{err['column']}: error: {err['message']}"


def format_ignore(ignore: tuple[int, list[str]]) -> str:
    line, codes = ignore
    if not codes:
        return f"ignore: {line}"
    else:
        return f"ignore: {line} [{', '.join(codes)}]"


def load_tree(node: MypyFile, options: Options) -> list[ParseError]:
    """Deserialize full AST from serialized raw data."""
    assert node.raw_data is not None
    state = State(options)
    data = ReadBuffer(node.raw_data.defs)
    n = read_int(data)
    node.defs = read_statements(state, data, n)
    node.imports = deserialize_imports(node.raw_data.imports)
    node.raw_data = None
    return state.errors


def test_parser_imports(testcase: DataDrivenTestCase) -> None:
    """Perform a single native parser imports test case.

    The argument contains the description of the test case.
    This test outputs only reachable import information.
    """
    options = Options()
    options.hide_error_codes = True
    options.python_version = (3, 10)

    source = "\n".join(testcase.input)

    try:
        with temp_source(source) as fnam:
            node, errors, type_ignores = native_parse(fnam, options)
            errors += load_tree(node, options)
            # Extract and format reachable imports
            a = format_reachable_imports(node)
            a = [format_error(err) for err in errors] + a
    except CompileError as e:
        a = e.messages

    assert_string_arrays_equal(
        testcase.output, a, f"Invalid parser output ({testcase.file}, line {testcase.line})"
    )


def format_reachable_imports(node: MypyFile) -> list[str]:
    """Format reachable imports from a MypyFile node.

    Returns a list of strings representing reachable imports with line numbers and flags.
    """
    from mypy.nodes import Import, ImportAll, ImportFrom

    output: list[str] = []

    # Filter for reachable imports (is_unreachable == False)
    reachable_imports = [imp for imp in node.imports if not imp.is_unreachable]

    for imp in reachable_imports:
        line_num = imp.line

        # Collect flags (only show when flag is False/not set)
        flags = []
        if not imp.is_top_level:
            flags.append("not top_level")
        if imp.is_mypy_only:
            flags.append("mypy_only")

        flags_str = " [" + ", ".join(flags) + "]" if flags else ""

        if isinstance(imp, Import):
            # Format: line: import foo [as bar] [flags]
            for module_id, as_id in imp.ids:
                if as_id:
                    output.append(f"{line_num}: import {module_id} as {as_id}{flags_str}")
                else:
                    output.append(f"{line_num}: import {module_id}{flags_str}")
        elif isinstance(imp, ImportFrom):
            # Format: line: from foo import bar, baz [as b] [flags]
            # Handle relative imports
            if imp.relative > 0:
                prefix = "." * imp.relative
                if imp.id:
                    module = f"{prefix}{imp.id}"
                else:
                    module = prefix
            else:
                module = imp.id

            # Group all names together
            name_parts = []
            for name, as_name in imp.names:
                if as_name:
                    name_parts.append(f"{name} as {as_name}")
                else:
                    name_parts.append(name)

            names_str = ", ".join(name_parts)
            output.append(f"{line_num}: from {module} import {names_str}{flags_str}")
        elif isinstance(imp, ImportAll):
            # Format: line: from foo import * [flags]
            # Handle relative imports
            if imp.relative > 0:
                prefix = "." * imp.relative
                if imp.id:
                    module = f"{prefix}{imp.id}"
                else:
                    module = prefix
            else:
                module = imp.id

            output.append(f"{line_num}: from {module} import *{flags_str}")

    return output


@unittest.skipUnless(has_nativeparse, "nativeparse not available")
class TestNativeParserBinaryFormat(unittest.TestCase):
    def test_trivial_binary_data(self) -> None:
        # A quick sanity check to ensure the serialized data looks as expected. Only covers
        # a few AST nodes.

        def int_enc(n: int) -> int:
            return (n + 10) << 1

        def locs(start_line: int, start_column: int, end_line: int, end_column: int) -> list[int]:
            return [
                LOCATION,
                int_enc(start_line),
                int_enc(start_column),
                int_enc(end_line - start_line),
                int_enc(end_column - start_column),
            ]

        with temp_source("print('hello')") as fnam:
            b, _, _, _, _, _, _, _ = parse_to_binary_ast(fnam, Options())
            assert list(b) == (
                [LITERAL_INT, 22, nodes.EXPR_STMT, nodes.CALL_EXPR]
                + [nodes.NAME_EXPR, LITERAL_STR]
                + [int_enc(5)]
                + list(b"print")
                + locs(1, 0, 1, 5)
                + [END_TAG, LIST_GEN, 22, nodes.STR_EXPR]
                + [LITERAL_STR, int_enc(5)]
                + list(b"hello")
                + locs(1, 6, 1, 13)
                + [END_TAG]
                # arg_kinds: [ARG_POS]
                + [LIST_INT, 22, int_enc(0)]
                # arg_names: [None]
                + [LIST_GEN, 22, LITERAL_NONE]
                + locs(1, 0, 1, 14)
                + [END_TAG, END_TAG]
            )


@contextlib.contextmanager
def temp_source(text: str) -> Iterator[str]:
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = os.path.join(temp_dir, "t.py")
        with open(temp_path, "w") as f:
            f.write(text)
        yield temp_path
