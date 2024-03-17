#!/usr/bin/env python3
"""
Analyze docstrings to detect errors.

If no argument is provided, it does a quick check of docstrings and returns
a csv with all API functions and results of basic checks.

If a function or method is provided in the form "pandas.function",
"pandas.module.class.method", etc. a list of all errors in the docstring for
the specified function or method.

Usage::
    $ ./validate_docstrings.py
    $ ./validate_docstrings.py pandas.DataFrame.head
"""
from __future__ import annotations

import argparse
import doctest
import importlib
import json
import os
import pathlib
import subprocess
import sys
import tempfile

import matplotlib
import matplotlib.pyplot as plt
from numpydoc.docscrape import get_doc_object
from numpydoc.validate import (
    ERROR_MSGS as NUMPYDOC_ERROR_MSGS,
    Validator,
    validate,
)

# With template backend, matplotlib plots nothing
matplotlib.use("template")

# Styler methods are Jinja2 objects who's docstrings we don't own.
IGNORE_VALIDATION = {
    "Styler.env",
    "Styler.template_html",
    "Styler.template_html_style",
    "Styler.template_html_table",
    "Styler.template_latex",
    "Styler.template_string",
    "Styler.loader",
    "errors.InvalidComparison",
    "errors.LossySetitemError",
    "errors.NoBufferPresent",
    "errors.IncompatibilityWarning",
    "errors.PyperclipException",
    "errors.PyperclipWindowsException",
}
PRIVATE_CLASSES = ["NDFrame", "IndexOpsMixin"]
ERROR_MSGS = {
    "GL04": "Private classes ({mentioned_private_classes}) should not be "
            "mentioned in public docstrings",
    "PD01": "Use 'array-like' rather than 'array_like' in docstrings.",
    "SA05": "{reference_name} in `See Also` section does not need `pandas` "
            "prefix, use {right_reference} instead.",
    "EX03": "flake8 error: line {line_number}, col {col_number}: {error_code} "
            "{error_message}",
    "EX04": "Do not import {imported_library}, as it is imported "
            "automatically for the examples (numpy as np, pandas as pd)",
}


def pandas_error(code, **kwargs):
    """
    Copy of the numpydoc error function, since ERROR_MSGS can't be updated
    with our custom errors yet.
    """
    return code, ERROR_MSGS[code].format(**kwargs)


def get_api_items(api_doc_fd):
    """
    Yield information about all public API items.

    Parse api.rst file from the documentation, and extract all the functions,
    methods, classes, attributes... This should include all pandas public API.

    Parameters
    ----------
    api_doc_fd : file descriptor
        A file descriptor of the API documentation page, containing the table
        of contents with all the public API.

    Yields
    ------
    name : str
        The name of the object (e.g. 'pandas.Series.str.upper').
    func : function
        The object itself. In most cases this will be a function or method,
        but it can also be classes, properties, cython objects...
    section : str
        The name of the section in the API page where the object item is
        located.
    subsection : str
        The name of the subsection in the API page where the object item is
        located.
    """
    current_module = "pandas"
    previous_line = current_section = current_subsection = ""
    position = None
    for line in api_doc_fd:
        line_stripped = line.strip()
        if len(line_stripped) == len(previous_line):
            if set(line_stripped) == set("-"):
                current_section = previous_line
                continue
            if set(line_stripped) == set("~"):
                current_subsection = previous_line
                continue

        if line_stripped.startswith(".. currentmodule::"):
            current_module = line_stripped.replace(".. currentmodule::", "").strip()
            continue

        if line_stripped == ".. autosummary::":
            position = "autosummary"
            continue

        if position == "autosummary":
            if line_stripped == "":
                position = "items"
                continue

        if position == "items":
            if line_stripped == "":
                position = None
                continue
            if line_stripped in IGNORE_VALIDATION:
                continue
            func = importlib.import_module(current_module)
            for part in line_stripped.split("."):
                func = getattr(func, part)

            yield (
                f"{current_module}.{line_stripped}",
                func,
                current_section,
                current_subsection,
            )

        previous_line = line_stripped


def validate_pep8_for_examples(docs: list[PandasDocstring] | PandasDocstring
                               ) -> dict[PandasDocstring, list[tuple]]:
    """
    Call the pep8 validation for docstrings with examples and add the found errors.

    Parameters
    ----------
    docs : list[PandasDocString]
        List of docstrings to validate.

    Returns
    -------
    dict[PandasDocstring, list]
       Dict of function names and the pep8 error messages found in their docstrings.
       The errors messages are of the form
       (error_code, message, line_number, col_number).
    """
    if isinstance(docs, PandasDocstring):
        docs = [docs]

    with tempfile.TemporaryDirectory() as temp_dir:
        doc_to_filename = {}
        for doc in docs:
            if not doc.examples:
                continue

            # F401 is needed to not generate flake8 errors in examples
            # that do not use numpy or pandas
            content = "".join(
                (
                    "import numpy as np  # noqa: F401\n",
                    "import pandas as pd  # noqa: F401\n",
                    *doc.examples_source_code,
                )
            )

            temp_file = tempfile.NamedTemporaryFile(mode="w",
                                                    dir=temp_dir,
                                                    encoding="utf-8",
                                                    delete=False)
            temp_file.write(content)
            temp_file.flush()
            doc_to_filename[doc] = temp_file.name

        # No docs with examples to process
        if not doc_to_filename:
            return {}

        cmd = [
            sys.executable,
            "-m",
            "flake8",
            "--format=%(row)d\t%(col)d\t%(code)s\t%(text)s",
            "--max-line-length=88",
            "--ignore=E203,E3,W503,W504,E402,E731,E128,E124,E704",
        ]
        cmd.extend(doc_to_filename.values())
        response = subprocess.run(cmd, capture_output=True, check=False,
                                  text=True)

    all_docs_error_messages = {doc: [] for doc in docs}
    for doc, temp_file_name in doc_to_filename.items():
        # one output for each error, each error must be mapped to the func_name
        for output in ("stdout", "stderr"):
            out = getattr(response, output)
            out = out.replace(temp_file_name, "").strip("\n").splitlines()
            if out:
                all_docs_error_messages[doc].extend(out)

    for doc, raw_error_messages in all_docs_error_messages.items():
        doc_error_messages = []
        for raw_error_message in raw_error_messages:
            line_num, col_num, err_code, msg = raw_error_message.split("\t", maxsplit=3)
            # Note: we subtract 2 from the line number because
            # 'import numpy as np\nimport pandas as pd\n'
            # is prepended to the docstrings.
            doc_error_messages.append(
                (
                    err_code,
                    msg,
                    int(line_num) - 2,
                    int(col_num)
                )
            )
        all_docs_error_messages[doc] = doc_error_messages

    for doc in docs:
        if doc.examples and doc not in all_docs_error_messages.keys():
            raise KeyError(f"Docstring\n###\n{doc}\n###\nhas examples but "
                           f"no pep8 validation results.")

    return all_docs_error_messages


class PandasDocstring(Validator):
    def __init__(self, func_name: str, doc_obj=None) -> None:
        self.func_name = func_name
        if doc_obj is None:
            doc_obj = get_doc_object(Validator._load_obj(func_name))
        super().__init__(doc_obj)

    @property
    def name(self):
        return self.func_name

    @property
    def mentioned_private_classes(self):
        return [klass for klass in PRIVATE_CLASSES if klass in self.raw_doc]

    @property
    def examples_source_code(self):
        lines = doctest.DocTestParser().get_examples(self.raw_doc)
        return [line.source for line in lines]

    def non_hyphenated_array_like(self):
        return "array_like" in self.raw_doc


def pandas_validate_single_func(func_name: str) -> dict:
    func_names = [func_name]
    results = pandas_validate(func_names)
    return results[func_name]


def pandas_validate(func_names: str | list[str]) -> dict[str, dict]:
    """
    Call the numpydoc validation, and add the errors specific to pandas.

    Parameters
    ----------
    func_names : list[str]
        List of name of the object of the docstrings to validate.

    Returns
    -------
    dict[str, dict]
        For each function, information about the docstring and the errors found.
    """
    if isinstance(func_names, str):
        func_names = [func_names]

    docs_to_results = {}
    for func_name in func_names:
        func_obj = Validator._load_obj(func_name)
        # Some objects are instances, e.g. IndexSlice, which numpydoc can't validate
        doc_obj = get_doc_object(func_obj, doc=func_obj.__doc__)
        doc = PandasDocstring(func_name, doc_obj)
        result = validate(doc_obj)
        docs_to_results[doc] = result

    # add errors not from examples to the result
    for doc, result in docs_to_results.items():
        mentioned_errs = doc.mentioned_private_classes
        if mentioned_errs:
            result["errors"].append(
                pandas_error(
                    "GL04",
                    mentioned_private_classes=", ".join(mentioned_errs))
            )

        if doc.see_also:
            see_also_prefix_errors = [
                pandas_error("SA05",
                             reference_name=rel_name,
                             right_reference=rel_name[len("pandas."):],
                             )
                for rel_name in doc.see_also
                if rel_name.startswith("pandas.")
            ]
            result["errors"].extend(see_also_prefix_errors)

        if doc.non_hyphenated_array_like():
            result["errors"].append(pandas_error("PD01"))

    pep8_results = validate_pep8_for_examples(list(docs_to_results.keys()))

    for doc, pep8_errors in pep8_results.items():
        result = docs_to_results[doc]
        pep8_pandas_errors = [
            pandas_error(
                "EX03",
                error_code=err_code,
                error_message=err_msg,
                line_number=line_number,
                col_number=col_number,
            ) for err_code, err_msg, line_number, col_number in pep8_errors
        ]
        result["errors"].extend(pep8_pandas_errors)
        examples_source_code = "".join(doc.examples_source_code)
        import_errors = [pandas_error("EX04", imported_library=wrong_import)
                         for wrong_import in ("numpy", "pandas")
                         if f"import {wrong_import}" in examples_source_code]
        result["errors"].extend(import_errors)

    plt.close("all")
    validation_results = {doc.func_name: result
                          for doc, result
                          in docs_to_results.items()}
    return validation_results


def validate_all(prefix, ignore_deprecated=False):
    """
    Execute the validation of all docstrings, and return a dict with the
    results.

    Parameters
    ----------
    prefix : str or None
        If provided, only the docstrings that start with this pattern will be
        validated. If None, all docstrings will be validated.
    ignore_deprecated: bool, default False
        If True, deprecated objects are ignored when validating docstrings.

    Returns
    -------
    dict
        A dictionary with an item for every function/method... containing
        all the validation information.
    """
    result = {}
    seen = {}

    def matches_prefix(func_name):
        return func_name.startswith(prefix) if prefix else True

    func_names = [func_name for func_name, _, _, _
                  in get_all_api_items()
                  if matches_prefix(func_name)]
    doc_infos = pandas_validate(func_names)

    for func_name, _, section, subsection in get_all_api_items():
        doc_info = doc_infos[func_name]
        if ignore_deprecated and doc_info["deprecated"]:
            continue
        result[func_name] = doc_info

        shared_code_key = doc_info["file"], doc_info["file_line"]
        shared_code = seen.get(shared_code_key, "")
        result[func_name].update(
            {
                "in_api": True,
                "section": section,
                "subsection": subsection,
                "shared_code_with": shared_code,
            }
        )

        seen[shared_code_key] = func_name

    return result


def get_all_api_items():
    base_path = pathlib.Path(__file__).parent.parent
    api_doc_fnames = pathlib.Path(base_path, "doc", "source", "reference")
    for api_doc_fname in api_doc_fnames.glob("*.rst"):
        with open(api_doc_fname, encoding="utf-8") as f:
            yield from get_api_items(f)


def print_validate_all_results(
        output_format: str,
        prefix: str | None,
        errors: list[str] | None,
        ignore_deprecated: bool,
        ignore_errors: dict[str, list[str]] | None,
):
    if output_format not in ("default", "json", "actions"):
        raise ValueError(f'Unknown output_format "{output_format}"')
    if ignore_errors is None:
        ignore_errors = {}

    result = validate_all(prefix, ignore_deprecated)

    if output_format == "json":
        sys.stdout.write(json.dumps(result))
        return 0

    prefix = "##[error]" if output_format == "actions" else ""
    exit_status = 0
    for func_name, res in result.items():
        for err_code, err_desc in res["errors"]:
            is_not_requested_error = errors and err_code not in errors
            is_ignored_error = err_code in ignore_errors.get(func_name, [])
            if is_not_requested_error or is_ignored_error:
                continue

            sys.stdout.write(
                f'{prefix}{res["file"]}:{res["file_line"]}:'
                f"{err_code}:{func_name}:{err_desc}\n"
            )
            exit_status += 1

    return exit_status


def print_validate_one_results(func_name: str) -> None:
    def header(title, width=80, char="#") -> str:
        full_line = char * width
        side_len = (width - len(title) - 2) // 2
        adj = "" if len(title) % 2 == 0 else " "
        title_line = f"{char * side_len} {title}{adj} {char * side_len}"

        return f"\n{full_line}\n{title_line}\n{full_line}\n\n"

    result = pandas_validate(func_name)

    sys.stderr.write(header(f"Docstring ({func_name})"))
    sys.stderr.write(f"{result['docstring']}\n")

    sys.stderr.write(header("Validation"))
    if result["errors"]:
        sys.stderr.write(f'{len(result["errors"])} Errors found for `{func_name}`:\n')
        for err_code, err_desc in result["errors"]:
            sys.stderr.write(f"\t{err_code}\t{err_desc}\n")
    else:
        sys.stderr.write(f'Docstring for "{func_name}" correct. :)\n')


def validate_error_codes(errors):
    overlapped_errors = set(NUMPYDOC_ERROR_MSGS).intersection(set(ERROR_MSGS))
    assert not overlapped_errors, f"{overlapped_errors} is overlapped."
    all_errors = set(NUMPYDOC_ERROR_MSGS).union(set(ERROR_MSGS))
    nonexistent_errors = set(errors) - all_errors
    assert not nonexistent_errors, f"{nonexistent_errors} don't exist."


def main(
        func_name,
        output_format,
        prefix,
        errors,
        ignore_deprecated,
        ignore_errors
):
    """
    Main entry point. Call the validation for one or for all docstrings.
    """
    if func_name is None:
        error_str = ", ".join(errors)
        msg = f"Validate docstrings ({error_str})\n"
    else:
        msg = f"Validate docstring in function {func_name}\n"
    sys.stdout.write(msg)

    validate_error_codes(errors)
    if ignore_errors is not None:
        for error_codes in ignore_errors.values():
            validate_error_codes(error_codes)

    if func_name is None:
        exit_status = print_validate_all_results(
            output_format,
            prefix,
            errors,
            ignore_deprecated,
            ignore_errors
        )
    else:
        print_validate_one_results(func_name)
        exit_status = 0
    sys.stdout.write(msg + "DONE" + os.linesep)

    return exit_status


if __name__ == "__main__":
    format_opts = "default", "json", "actions"
    func_help = (
        "function or method to validate (e.g. pandas.DataFrame.head) "
        "if not provided, all docstrings are validated and returned "
        "as JSON"
    )
    argparser = argparse.ArgumentParser(description="validate pandas docstrings")
    argparser.add_argument(
        "function",
        nargs="?",
        default=None,
        help=func_help)
    argparser.add_argument(
        "--format",
        default="default",
        choices=format_opts,
        help="format of the output when validating "
             "multiple docstrings (ignored when validating one). "
             "It can be {str(format_opts)[1:-1]}",
    )
    argparser.add_argument(
        "--prefix",
        default=None,
        help="pattern for the "
             "docstring names, in order to decide which ones "
             'will be validated. A prefix "pandas.Series.str."'
             "will make the script validate all the docstrings "
             "of methods starting by this pattern. It is "
             "ignored if parameter function is provided",
    )
    argparser.add_argument(
        "--errors",
        default=None,
        help="comma separated "
             "list of error codes to validate. By default it "
             "validates all errors (ignored when validating "
             "a single docstring)",
    )
    argparser.add_argument(
        "--ignore_deprecated",
        default=False,
        action="store_true",
        help="if this flag is set, "
             "deprecated objects are ignored when validating "
             "all docstrings",
    )
    argparser.add_argument(
        "--ignore_errors",
        default=None,
        action="append",
        nargs=2,
        metavar=("function", "error_codes"),
        help="function for which comma separated list "
             "of error codes should not be validated"
             "(e.g. pandas.DataFrame.head PR01,SA01). "
             "Partial validation for more than one function"
             "can be achieved by repeating this parameter.",
    )
    args = argparser.parse_args(sys.argv[1:])

    args.errors = args.errors.split(",") if args.errors else None
    if args.ignore_errors:
        args.ignore_errors = {function: error_codes.split(",")
                              for function, error_codes
                              in args.ignore_errors}

    sys.exit(
        main(args.function,
             args.format,
             args.prefix,
             args.errors,
             args.ignore_deprecated,
             args.ignore_errors
             )
    )
