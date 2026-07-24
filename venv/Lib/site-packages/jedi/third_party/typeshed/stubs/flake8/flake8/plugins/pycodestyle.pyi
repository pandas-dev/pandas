from collections.abc import Generator
from typing import Any

def pycodestyle_logical(
    blank_before: Any,
    blank_lines: Any,
    checker_state: Any,
    hang_closing: Any,
    indent_char: Any,
    indent_level: Any,
    indent_size: Any,
    line_number: Any,
    lines: Any,
    logical_line: Any,
    max_doc_length: Any,
    noqa: Any,
    previous_indent_level: Any,
    previous_logical: Any,
    previous_unindented_logical_line: Any,
    tokens: Any,
    verbose: Any,
) -> Generator[tuple[int, str]]: ...
def pycodestyle_physical(
    indent_char: Any,
    line_number: Any,
    lines: Any,
    max_line_length: Any,
    multiline: Any,
    noqa: Any,
    physical_line: Any,
    total_lines: Any,
) -> Generator[tuple[int, str]]: ...
