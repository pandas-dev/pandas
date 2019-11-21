#!/usr/bin/env python
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
import argparse
import doctest
import glob
import importlib
import inspect
import json
import os
import sys
import tempfile

import flake8.main.application

try:
    from io import StringIO
except ImportError:
    from cStringIO import StringIO

# Template backend makes matplotlib to not plot anything. This is useful
# to avoid that plot windows are open from the doctests while running the
# script. Setting here before matplotlib is loaded.
# We don't warn for the number of open plots, as none is actually being opened
os.environ["MPLBACKEND"] = "Template"
import matplotlib  # noqa: E402 isort:skip

matplotlib.rc("figure", max_open_warning=10000)

import numpy  # noqa: E402 isort:skip

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, os.path.join(BASE_PATH))
import pandas  # noqa: E402 isort:skip

sys.path.insert(1, os.path.join(BASE_PATH, "doc", "sphinxext"))
from numpydoc.docscrape import NumpyDocString  # noqa: E402 isort:skip
from numpydoc.validate import validate, Docstring, error  # noqa: E402 isort:skip
from pandas.io.formats.printing import pprint_thing  # noqa: E402 isort:skip


PRIVATE_CLASSES = ["NDFrame", "IndexOpsMixin"]
ERROR_MSGS = {
    "GL04": "Private classes ({mentioned_private_classes}) should not be "
    "mentioned in public docstrings",
    "SA05": "{reference_name} in `See Also` section does not need `pandas` "
    "prefix, use {right_reference} instead.",
    "EX02": "Examples do not pass tests:\n{doctest_log}",
    "EX03": "flake8 error: {error_code} {error_message}{times_happening}",
    "EX04": "Do not import {imported_library}, as it is imported "
    "automatically for the examples (numpy as np, pandas as pd)",
}


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
        The name of the object (e.g. 'pandas.Series.str.upper).
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
        line = line.strip()
        if len(line) == len(previous_line):
            if set(line) == set("-"):
                current_section = previous_line
                continue
            if set(line) == set("~"):
                current_subsection = previous_line
                continue

        if line.startswith(".. currentmodule::"):
            current_module = line.replace(".. currentmodule::", "").strip()
            continue

        if line == ".. autosummary::":
            position = "autosummary"
            continue

        if position == "autosummary":
            if line == "":
                position = "items"
                continue

        if position == "items":
            if line == "":
                position = None
                continue
            item = line.strip()
            func = importlib.import_module(current_module)
            for part in item.split("."):
                func = getattr(func, part)

            yield (
                ".".join([current_module, item]),
                func,
                current_section,
                current_subsection,
            )

        previous_line = line


class PandasDocstring(Docstring):
    @property
    def mentioned_private_classes(self):
        return [klass for klass in PRIVATE_CLASSES if klass in self.raw_doc]

    @property
    def examples_errors(self):
        flags = doctest.NORMALIZE_WHITESPACE | doctest.IGNORE_EXCEPTION_DETAIL
        finder = doctest.DocTestFinder()
        runner = doctest.DocTestRunner(optionflags=flags)
        context = {"np": numpy, "pd": pandas}
        error_msgs = ""
        for test in finder.find(self.raw_doc, self.name, globs=context):
            f = StringIO()
            runner.run(test, out=f.write)
            error_msgs += f.getvalue()
        return error_msgs

    @property
    def examples_source_code(self):
        lines = doctest.DocTestParser().get_examples(self.raw_doc)
        return [line.source for line in lines]

    def validate_pep8(self):
        if not self.examples:
            return

        # F401 is needed to not generate flake8 errors in examples
        # that do not user numpy or pandas
        content = "".join(
            (
                "import numpy as np  # noqa: F401\n",
                "import pandas as pd  # noqa: F401\n",
                *self.examples_source_code,
            )
        )

        application = flake8.main.application.Application()
        application.initialize(["--quiet"])

        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8") as file:
            file.write(content)
            file.flush()
            application.run_checks([file.name])

        # We need this to avoid flake8 printing the names of the files to
        # the standard output
        application.formatter.write = lambda line, source: None
        application.report()

        yield from application.guide.stats.statistics_for("")


def pandas_validation(doc):
    """
    Validation of errors specific to pandas.

    Parameters
    ----------
    doc : PandasDocString
        Instance of the PandasDocString corresponding to the docstring to validate.

    Returns
    -------
    errs : list of error (namedtuple)
        List of errors detected in the docstring.
    example_errs : str
        Error messages captured from running the examples.
    """
    errs = []

    mentioned_errs = doc.mentioned_private_classes
    if mentioned_errs:
        errs.append(error("GL04", mentioned_private_classes=", ".join(mentioned_errs)))

    if doc.see_also:
        for rel_name, rel_desc in doc.see_also.items():
            if rel_name.startswith("pandas."):
                errs.append(
                    error(
                        "SA05",
                        reference_name=rel_name,
                        right_reference=rel_name[len("pandas.") :],
                    )
                )

    examples_errs = ""
    if doc.examples:
        examples_errs = doc.examples_errors
        if examples_errs:
            errs.append(error("EX02", doctest_log=examples_errs))
        for err in doc.validate_pep8():
            errs.append(
                error(
                    "EX03",
                    error_code=err.error_code,
                    error_message=err.message,
                    times_happening=" ({} times)".format(err.count)
                    if err.count > 1
                    else "",
                )
            )
        examples_source_code = "".join(doc.examples_source_code)
        for wrong_import in ("numpy", "pandas"):
            if "import {}".format(wrong_import) in examples_source_code:
                errs.append(error("EX04", imported_library=wrong_import))
    return errs, examples_errs


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

    # functions from the API docs
    api_doc_fnames = os.path.join(BASE_PATH, "doc", "source", "reference", "*.rst")
    api_items = []
    for api_doc_fname in glob.glob(api_doc_fnames):
        with open(api_doc_fname) as f:
            api_items += list(get_api_items(f))
    for func_name, func_obj, section, subsection in api_items:
        if prefix and not func_name.startswith(prefix):
            continue
        doc_info = validate_one(func_name)
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

    # functions from introspecting Series and DataFrame
    api_item_names = set(list(zip(*api_items))[0])
    for class_ in (pandas.Series, pandas.DataFrame):
        for member in inspect.getmembers(class_):
            func_name = "pandas.{}.{}".format(class_.__name__, member[0])
            if not member[0].startswith("_") and func_name not in api_item_names:
                if prefix and not func_name.startswith(prefix):
                    continue
                doc_info = validate_one(func_name)
                if ignore_deprecated and doc_info["deprecated"]:
                    continue
                result[func_name] = doc_info
                result[func_name]["in_api"] = False

    return result


def main(func_name, prefix, errors, output_format, ignore_deprecated):
    def header(title, width=80, char="#"):
        full_line = char * width
        side_len = (width - len(title) - 2) // 2
        adj = "" if len(title) % 2 == 0 else " "
        title_line = "{side} {title}{adj} {side}".format(
            side=char * side_len, title=title, adj=adj
        )

        return "\n{full_line}\n{title_line}\n{full_line}\n\n".format(
            full_line=full_line, title_line=title_line
        )

    exit_status = 0
    if func_name is None:
        result = validate_all(prefix, ignore_deprecated)

        if output_format == "json":
            output = json.dumps(result)
        else:
            if output_format == "default":
                output_format = "{text}\n"
            elif output_format == "azure":
                output_format = (
                    "##vso[task.logissue type=error;"
                    "sourcepath={path};"
                    "linenumber={row};"
                    "code={code};"
                    "]{text}\n"
                )
            else:
                raise ValueError('Unknown output_format "{}"'.format(output_format))

            output = ""
            for name, res in result.items():
                for err_code, err_desc in res["errors"]:
                    # The script would be faster if instead of filtering the
                    # errors after validating them, it didn't validate them
                    # initially. But that would complicate the code too much
                    if errors and err_code not in errors:
                        continue
                    exit_status += 1
                    output += output_format.format(
                        name=name,
                        path=res["file"],
                        row=res["file_line"],
                        code=err_code,
                        text="{}: {}".format(name, err_desc),
                    )

        sys.stdout.write(output)

    else:
        result = validate_one(func_name)
        sys.stderr.write(header("Docstring ({})".format(func_name)))
        sys.stderr.write("{}\n".format(result["docstring"]))
        sys.stderr.write(header("Validation"))
        if result["errors"]:
            sys.stderr.write("{} Errors found:\n".format(len(result["errors"])))
            for err_code, err_desc in result["errors"]:
                # Failing examples are printed at the end
                if err_code == "EX02":
                    sys.stderr.write("\tExamples do not pass tests\n")
                    continue
                sys.stderr.write("\t{}\n".format(err_desc))
        if result["warnings"]:
            sys.stderr.write("{} Warnings found:\n".format(len(result["warnings"])))
            for wrn_code, wrn_desc in result["warnings"]:
                sys.stderr.write("\t{}\n".format(wrn_desc))

        if not result["errors"]:
            sys.stderr.write('Docstring for "{}" correct. :)\n'.format(func_name))

        if result["examples_errors"]:
            sys.stderr.write(header("Doctests"))
            sys.stderr.write(result["examples_errors"])

    return exit_status


if __name__ == "__main__":
    format_opts = "default", "json", "azure"
    func_help = (
        "function or method to validate (e.g. pandas.DataFrame.head) "
        "if not provided, all docstrings are validated and returned "
        "as JSON"
    )
    argparser = argparse.ArgumentParser(description="validate pandas docstrings")
    argparser.add_argument("function", nargs="?", default=None, help=func_help)
    argparser.add_argument(
        "--format",
        default="default",
        choices=format_opts,
        help="format of the output when validating "
        "multiple docstrings (ignored when validating one)."
        "It can be {}".format(str(format_opts)[1:-1]),
    )
    argparser.add_argument(
        "--prefix",
        default=None,
        help="pattern for the "
        "docstring names, in order to decide which ones "
        'will be validated. A prefix "pandas.Series.str.'
        "will make the script validate all the docstrings"
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

    args = argparser.parse_args()
    sys.exit(
        main(
            args.function,
            args.prefix,
            args.errors.split(",") if args.errors else None,
            args.format,
            args.ignore_deprecated,
        )
    )
