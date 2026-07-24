from __future__ import annotations

import re

from librt.internal import ReadBuffer

from mypy import errorcodes as codes
from mypy.cache import read_int
from mypy.errors import Errors
from mypy.nodes import FileRawData, MypyFile, ParseError
from mypy.options import Options


def parse(
    source: str | bytes,
    fnam: str,
    module: str | None,
    errors: Errors,
    options: Options,
    file_exists: bool,
    eager: bool = False,
) -> MypyFile:
    """Parse a source file, without doing any semantic analysis.

    Return the parse tree, use the errors object to report parse errors.
    The python_version (major, minor) option determines the Python syntax variant.

    New parser returns empty tree with serialized data. To get the full tree and
    the parse errors, use eager=True.
    """
    if options.native_parser:
        # Native parser only works with actual files on disk
        # Fall back to fastparse for in-memory source or non-existent files
        if file_exists:
            import mypy.nativeparse

            ignore_errors = options.ignore_errors or fnam in errors.ignored_files
            # If errors are ignored, we can drop many function bodies to speed up type checking.
            strip_function_bodies = ignore_errors and not options.preserve_asts
            tree, _, _ = mypy.nativeparse.native_parse(
                fnam, options, skip_function_bodies=strip_function_bodies
            )
            # Set is_stub based on file extension
            tree.is_stub = fnam.endswith(".pyi")
            # Note: tree.imports is populated directly by load_from_raw() with deserialized
            # import metadata, so we don't need to collect imports via AST traversal
            if eager and tree.raw_data is not None:
                tree = load_from_raw(fnam, module, tree.raw_data, errors, options)
            return tree
        # Fall through to fastparse for non-existent files

    if options.transform_source is not None:
        source = options.transform_source(source)
    import mypy.fastparse

    return mypy.fastparse.parse(source, fnam=fnam, module=module, errors=errors, options=options)


def load_from_raw(
    fnam: str,
    module: str | None,
    raw_data: FileRawData,
    errors: Errors,
    options: Options,
    imports_only: bool = False,
) -> MypyFile:
    """Load AST from parsed binary data and report stored errors.

    If imports_only is true, only deserialize imports and return a mostly
    empty AST.
    """
    from mypy.nativeparse import State, deserialize_imports, read_statements

    state = State(options)
    if imports_only:
        defs = []
    else:
        data = ReadBuffer(raw_data.defs)
        n = read_int(data)
        defs = read_statements(state, data, n)
    imports = deserialize_imports(raw_data.imports)

    tree = MypyFile(defs, imports)
    tree.path = fnam
    tree.ignored_lines = raw_data.ignored_lines
    tree.is_partial_stub_package = raw_data.is_partial_stub_package
    tree.uses_template_strings = raw_data.uses_template_strings
    tree.is_stub = fnam.endswith(".pyi")
    if module is not None:
        tree._fullname = module

    # Report parse errors, this replicates the logic in parse().
    all_errors = raw_data.raw_errors + state.errors
    errors.set_file(fnam, module, options=options)
    for error in all_errors:
        # Note we never raise in this function, so it should not be called in coordinator.
        report_parse_error(error, errors)
    if imports_only:
        # Preserve raw data when only de-serializing imports, it will be sent to
        # the parallel workers.
        tree.raw_data = raw_data
    return tree


def report_parse_error(error: ParseError, errors: Errors) -> None:
    message = error["message"]
    # Standardize error message by capitalizing the first word
    message = re.sub(r"^(\s*\w)", lambda m: m.group(1).upper(), message)
    # Respect blocker status from error, default to True for syntax errors
    is_blocker = error.get("blocker", True)
    error_code = error.get("code")
    if error_code is None:
        error_code = codes.SYNTAX
    else:
        # Fallback to [syntax] for backwards compatibility.
        error_code = codes.error_codes.get(error_code) or codes.SYNTAX
    errors.report(error["line"], error["column"], message, blocker=is_blocker, code=error_code)
