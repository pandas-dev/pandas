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
import ast
import doctest
import functools
import glob
import importlib
import inspect
import json
import os
import pydoc
import re
import sys
import tempfile
import textwrap

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
from pandas.io.formats.printing import pprint_thing  # noqa: E402 isort:skip


PRIVATE_CLASSES = ["NDFrame", "IndexOpsMixin"]
DIRECTIVES = ["versionadded", "versionchanged", "deprecated"]
DIRECTIVE_PATTERN = re.compile(rf"^\s*\.\. ({'|'.join(DIRECTIVES)})(?!::)", re.I | re.M)
ALLOWED_SECTIONS = [
    "Parameters",
    "Attributes",
    "Methods",
    "Returns",
    "Yields",
    "Other Parameters",
    "Raises",
    "Warns",
    "See Also",
    "Notes",
    "References",
    "Examples",
]
ERROR_MSGS = {
    "GL01": "Docstring text (summary) should start in the line immediately "
    "after the opening quotes (not in the same line, or leaving a "
    "blank line in between)",
    "GL02": "Closing quotes should be placed in the line after the last text "
    "in the docstring (do not close the quotes in the same line as "
    "the text, or leave a blank line between the last text and the "
    "quotes)",
    "GL03": "Double line break found; please use only one blank line to "
    "separate sections or paragraphs, and do not leave blank lines "
    "at the end of docstrings",
    "GL04": "Private classes ({mentioned_private_classes}) should not be "
    "mentioned in public docstrings",
    "GL05": 'Tabs found at the start of line "{line_with_tabs}", please use '
    "whitespace only",
    "GL06": 'Found unknown section "{section}". Allowed sections are: '
    "{allowed_sections}",
    "GL07": "Sections are in the wrong order. Correct order is: {correct_sections}",
    "GL08": "The object does not have a docstring",
    "GL09": "Deprecation warning should precede extended summary",
    "GL10": "reST directives {directives} must be followed by two colons",
    "SS01": "No summary found (a short summary in a single line should be "
    "present at the beginning of the docstring)",
    "SS02": "Summary does not start with a capital letter",
    "SS03": "Summary does not end with a period",
    "SS04": "Summary contains heading whitespaces",
    "SS05": "Summary must start with infinitive verb, not third person "
    '(e.g. use "Generate" instead of "Generates")',
    "SS06": "Summary should fit in a single line",
    "ES01": "No extended summary found",
    "PR01": "Parameters {missing_params} not documented",
    "PR02": "Unknown parameters {unknown_params}",
    "PR03": "Wrong parameters order. Actual: {actual_params}. "
    "Documented: {documented_params}",
    "PR04": 'Parameter "{param_name}" has no type',
    "PR05": 'Parameter "{param_name}" type should not finish with "."',
    "PR06": 'Parameter "{param_name}" type should use "{right_type}" instead '
    'of "{wrong_type}"',
    "PR07": 'Parameter "{param_name}" has no description',
    "PR08": 'Parameter "{param_name}" description should start with a '
    "capital letter",
    "PR09": 'Parameter "{param_name}" description should finish with "."',
    "PR10": 'Parameter "{param_name}" requires a space before the colon '
    "separating the parameter name and type",
    "RT01": "No Returns section found",
    "RT02": "The first line of the Returns section should contain only the "
    "type, unless multiple values are being returned",
    "RT03": "Return value has no description",
    "RT04": "Return value description should start with a capital letter",
    "RT05": 'Return value description should finish with "."',
    "YD01": "No Yields section found",
    "SA01": "See Also section not found",
    "SA02": "Missing period at end of description for See Also "
    '"{reference_name}" reference',
    "SA03": "Description should be capitalized for See Also "
    '"{reference_name}" reference',
    "SA04": 'Missing description for See Also "{reference_name}" reference',
    "SA05": "{reference_name} in `See Also` section does not need `pandas` "
    "prefix, use {right_reference} instead.",
    "EX01": "No examples section found",
    "EX02": "Examples do not pass tests:\n{doctest_log}",
    "EX03": "flake8 error: {error_code} {error_message}{times_happening}",
    "EX04": "Do not import {imported_library}, as it is imported "
    "automatically for the examples (numpy as np, pandas as pd)",
}


def error(code, **kwargs):
    """
    Return a tuple with the error code and the message with variables replaced.

    This is syntactic sugar so instead of:
    - `('EX02', ERROR_MSGS['EX02'].format(doctest_log=log))`

    We can simply use:
    - `error('EX02', doctest_log=log)`

    Parameters
    ----------
    code : str
        Error code.
    **kwargs
        Values for the variables in the error messages

    Returns
    -------
    code : str
        Error code.
    message : str
        Error message with variables replaced.
    """
    return (code, ERROR_MSGS[code].format(**kwargs))


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


class Docstring:
    def __init__(self, name):
        self.name = name
        obj = self._load_obj(name)
        self.obj = obj
        self.code_obj = self._to_original_callable(obj)
        self.raw_doc = obj.__doc__ or ""
        self.clean_doc = pydoc.getdoc(obj)
        self.doc = NumpyDocString(self.clean_doc)

    def __len__(self) -> int:
        return len(self.raw_doc)

    @staticmethod
    def _load_obj(name):
        """
        Import Python object from its name as string.

        Parameters
        ----------
        name : str
            Object name to import (e.g. pandas.Series.str.upper)

        Returns
        -------
        object
            Python object that can be a class, method, function...

        Examples
        --------
        >>> Docstring._load_obj('pandas.Series')
        <class 'pandas.core.series.Series'>
        """
        for maxsplit in range(1, name.count(".") + 1):
            # TODO when py3 only replace by: module, *func_parts = ...
            func_name_split = name.rsplit(".", maxsplit)
            module = func_name_split[0]
            func_parts = func_name_split[1:]
            try:
                obj = importlib.import_module(module)
            except ImportError:
                pass
            else:
                continue

        if "obj" not in locals():
            raise ImportError(f'No module can be imported from "{name}"')

        for part in func_parts:
            obj = getattr(obj, part)
        return obj

    @staticmethod
    def _to_original_callable(obj):
        """
        Find the Python object that contains the source code of the object.

        This is useful to find the place in the source code (file and line
        number) where a docstring is defined. It does not currently work for
        all cases, but it should help find some (properties...).
        """
        while True:
            if inspect.isfunction(obj) or inspect.isclass(obj):
                f = inspect.getfile(obj)
                if f.startswith("<") and f.endswith(">"):
                    return None
                return obj
            if inspect.ismethod(obj):
                obj = obj.__func__
            elif isinstance(obj, functools.partial):
                obj = obj.func
            elif isinstance(obj, property):
                obj = obj.fget
            else:
                return None

    @property
    def type(self):
        return type(self.obj).__name__

    @property
    def is_function_or_method(self):
        # TODO(py27): remove ismethod
        return inspect.isfunction(self.obj) or inspect.ismethod(self.obj)

    @property
    def source_file_name(self):
        """
        File name where the object is implemented (e.g. pandas/core/frame.py).
        """
        try:
            fname = inspect.getsourcefile(self.code_obj)
        except TypeError:
            # In some cases the object is something complex like a cython
            # object that can't be easily introspected. An it's better to
            # return the source code file of the object as None, than crash
            pass
        else:
            if fname:
                fname = os.path.relpath(fname, BASE_PATH)
                return fname

    @property
    def source_file_def_line(self):
        """
        Number of line where the object is defined in its file.
        """
        try:
            return inspect.getsourcelines(self.code_obj)[-1]
        except (OSError, TypeError):
            # In some cases the object is something complex like a cython
            # object that can't be easily introspected. An it's better to
            # return the line number as None, than crash
            pass

    @property
    def github_url(self):
        url = "https://github.com/pandas-dev/pandas/blob/master/"
        url += "{}#L{}".format(self.source_file_name, self.source_file_def_line)
        return url

    @property
    def start_blank_lines(self):
        i = None
        if self.raw_doc:
            for i, row in enumerate(self.raw_doc.split("\n")):
                if row.strip():
                    break
        return i

    @property
    def end_blank_lines(self):
        i = None
        if self.raw_doc:
            for i, row in enumerate(reversed(self.raw_doc.split("\n"))):
                if row.strip():
                    break
        return i

    @property
    def double_blank_lines(self):
        prev = True
        for row in self.raw_doc.split("\n"):
            if not prev and not row.strip():
                return True
            prev = row.strip()
        return False

    @property
    def section_titles(self):
        sections = []
        self.doc._doc.reset()
        while not self.doc._doc.eof():
            content = self.doc._read_to_next_section()
            if (
                len(content) > 1
                and len(content[0]) == len(content[1])
                and set(content[1]) == {"-"}
            ):
                sections.append(content[0])
        return sections

    @property
    def summary(self):
        return " ".join(self.doc["Summary"])

    @property
    def num_summary_lines(self):
        return len(self.doc["Summary"])

    @property
    def extended_summary(self):
        if not self.doc["Extended Summary"] and len(self.doc["Summary"]) > 1:
            return " ".join(self.doc["Summary"])
        return " ".join(self.doc["Extended Summary"])

    @property
    def needs_summary(self):
        return not (bool(self.summary) and bool(self.extended_summary))

    @property
    def doc_parameters(self):
        parameters = {}
        for names, type_, desc in self.doc["Parameters"]:
            for name in names.split(", "):
                parameters[name] = (type_, "".join(desc))
        return parameters

    @property
    def signature_parameters(self):
        def add_stars(param_name: str, info: inspect.Parameter):
            """
            Add stars to *args and **kwargs parameters
            """
            if info.kind == inspect.Parameter.VAR_POSITIONAL:
                return f"*{param_name}"
            elif info.kind == inspect.Parameter.VAR_KEYWORD:
                return f"**{param_name}"
            else:
                return param_name

        if inspect.isclass(self.obj):
            if hasattr(self.obj, "_accessors") and (
                self.name.split(".")[-1] in self.obj._accessors
            ):
                # accessor classes have a signature but don't want to show this
                return tuple()
        try:
            sig = inspect.signature(self.obj)
        except (TypeError, ValueError):
            # Some objects, mainly in C extensions do not support introspection
            # of the signature
            return tuple()

        params = tuple(
            add_stars(parameter, sig.parameters[parameter])
            for parameter in sig.parameters
        )
        if params and params[0] in ("self", "cls"):
            return params[1:]
        return params

    @property
    def parameter_mismatches(self):
        errs = []
        signature_params = self.signature_parameters
        doc_params = tuple(self.doc_parameters)
        missing = set(signature_params) - set(doc_params)
        if missing:
            errs.append(error("PR01", missing_params=pprint_thing(missing)))
        extra = set(doc_params) - set(signature_params)
        if extra:
            errs.append(error("PR02", unknown_params=pprint_thing(extra)))
        if (
            not missing
            and not extra
            and signature_params != doc_params
            and not (not signature_params and not doc_params)
        ):
            errs.append(
                error(
                    "PR03", actual_params=signature_params, documented_params=doc_params
                )
            )

        return errs

    @property
    def correct_parameters(self):
        return not bool(self.parameter_mismatches)

    @property
    def directives_without_two_colons(self):
        return DIRECTIVE_PATTERN.findall(self.raw_doc)

    def parameter_type(self, param):
        return self.doc_parameters[param][0]

    def parameter_desc(self, param):
        desc = self.doc_parameters[param][1]
        # Find and strip out any sphinx directives
        for directive in DIRECTIVES:
            full_directive = ".. {}".format(directive)
            if full_directive in desc:
                # Only retain any description before the directive
                desc = desc[: desc.index(full_directive)]
        return desc

    @property
    def see_also(self):
        result = {}
        for funcs, desc in self.doc["See Also"]:
            for func, _ in funcs:
                result[func] = "".join(desc)

        return result

    @property
    def examples(self):
        return self.doc["Examples"]

    @property
    def returns(self):
        return self.doc["Returns"]

    @property
    def yields(self):
        return self.doc["Yields"]

    @property
    def method_source(self):
        try:
            source = inspect.getsource(self.obj)
        except TypeError:
            return ""
        return textwrap.dedent(source)

    @property
    def method_returns_something(self):
        """
        Check if the docstrings method can return something.

        Bare returns, returns valued None and returns from nested functions are
        disconsidered.

        Returns
        -------
        bool
            Whether the docstrings method can return something.
        """

        def get_returns_not_on_nested_functions(node):
            returns = [node] if isinstance(node, ast.Return) else []
            for child in ast.iter_child_nodes(node):
                # Ignore nested functions and its subtrees.
                if not isinstance(child, ast.FunctionDef):
                    child_returns = get_returns_not_on_nested_functions(child)
                    returns.extend(child_returns)
            return returns

        tree = ast.parse(self.method_source).body
        if tree:
            returns = get_returns_not_on_nested_functions(tree[0])
            return_values = [r.value for r in returns]
            # Replace NameConstant nodes valued None for None.
            for i, v in enumerate(return_values):
                if isinstance(v, ast.NameConstant) and v.value is None:
                    return_values[i] = None
            return any(return_values)
        else:
            return False

    @property
    def first_line_ends_in_dot(self):
        if self.doc:
            return self.doc.split("\n")[0][-1] == "."

    @property
    def deprecated(self):
        return ".. deprecated:: " in (self.summary + self.extended_summary)

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


def get_validation_data(doc):
    """
    Validate the docstring.

    Parameters
    ----------
    doc : Docstring
        A Docstring object with the given function name.

    Returns
    -------
    tuple
        errors : list of tuple
            Errors occurred during validation.
        warnings : list of tuple
            Warnings occurred during validation.
        examples_errs : str
            Examples usage displayed along the error, otherwise empty string.

    Notes
    -----
    The errors codes are defined as:
    - First two characters: Section where the error happens:
       * GL: Global (no section, like section ordering errors)
       * SS: Short summary
       * ES: Extended summary
       * PR: Parameters
       * RT: Returns
       * YD: Yields
       * RS: Raises
       * WN: Warns
       * SA: See Also
       * NT: Notes
       * RF: References
       * EX: Examples
    - Last two characters: Numeric error code inside the section

    For example, EX02 is the second codified error in the Examples section
    (which in this case is assigned to examples that do not pass the tests).

    The error codes, their corresponding error messages, and the details on how
    they are validated, are not documented more than in the source code of this
    function.
    """

    errs = []
    wrns = []
    if not doc.raw_doc:
        errs.append(error("GL08"))
        return errs, wrns, ""

    if doc.start_blank_lines != 1:
        errs.append(error("GL01"))
    if doc.end_blank_lines != 1:
        errs.append(error("GL02"))
    if doc.double_blank_lines:
        errs.append(error("GL03"))
    mentioned_errs = doc.mentioned_private_classes
    if mentioned_errs:
        errs.append(error("GL04", mentioned_private_classes=", ".join(mentioned_errs)))
    for line in doc.raw_doc.splitlines():
        if re.match("^ *\t", line):
            errs.append(error("GL05", line_with_tabs=line.lstrip()))

    unexpected_sections = [
        section for section in doc.section_titles if section not in ALLOWED_SECTIONS
    ]
    for section in unexpected_sections:
        errs.append(
            error("GL06", section=section, allowed_sections=", ".join(ALLOWED_SECTIONS))
        )

    correct_order = [
        section for section in ALLOWED_SECTIONS if section in doc.section_titles
    ]
    if correct_order != doc.section_titles:
        errs.append(error("GL07", correct_sections=", ".join(correct_order)))

    if doc.deprecated and not doc.extended_summary.startswith(".. deprecated:: "):
        errs.append(error("GL09"))

    directives_without_two_colons = doc.directives_without_two_colons
    if directives_without_two_colons:
        errs.append(error("GL10", directives=directives_without_two_colons))

    if not doc.summary:
        errs.append(error("SS01"))
    else:
        if not doc.summary[0].isupper():
            errs.append(error("SS02"))
        if doc.summary[-1] != ".":
            errs.append(error("SS03"))
        if doc.summary != doc.summary.lstrip():
            errs.append(error("SS04"))
        elif doc.is_function_or_method and doc.summary.split(" ")[0][-1] == "s":
            errs.append(error("SS05"))
        if doc.num_summary_lines > 1:
            errs.append(error("SS06"))

    if not doc.extended_summary:
        wrns.append(("ES01", "No extended summary found"))

    # PR01: Parameters not documented
    # PR02: Unknown parameters
    # PR03: Wrong parameters order
    errs += doc.parameter_mismatches

    for param in doc.doc_parameters:
        if not param.startswith("*"):  # Check can ignore var / kwargs
            if not doc.parameter_type(param):
                if ":" in param:
                    errs.append(error("PR10", param_name=param.split(":")[0]))
                else:
                    errs.append(error("PR04", param_name=param))
            else:
                if doc.parameter_type(param)[-1] == ".":
                    errs.append(error("PR05", param_name=param))
                common_type_errors = [
                    ("integer", "int"),
                    ("boolean", "bool"),
                    ("string", "str"),
                ]
                for wrong_type, right_type in common_type_errors:
                    if wrong_type in doc.parameter_type(param):
                        errs.append(
                            error(
                                "PR06",
                                param_name=param,
                                right_type=right_type,
                                wrong_type=wrong_type,
                            )
                        )
        if not doc.parameter_desc(param):
            errs.append(error("PR07", param_name=param))
        else:
            if not doc.parameter_desc(param)[0].isupper():
                errs.append(error("PR08", param_name=param))
            if doc.parameter_desc(param)[-1] != ".":
                errs.append(error("PR09", param_name=param))

    if doc.is_function_or_method:
        if not doc.returns:
            if doc.method_returns_something:
                errs.append(error("RT01"))
        else:
            if len(doc.returns) == 1 and doc.returns[0].name:
                errs.append(error("RT02"))
            for name_or_type, type_, desc in doc.returns:
                if not desc:
                    errs.append(error("RT03"))
                else:
                    desc = " ".join(desc)
                    if not desc[0].isupper():
                        errs.append(error("RT04"))
                    if not desc.endswith("."):
                        errs.append(error("RT05"))

        if not doc.yields and "yield" in doc.method_source:
            errs.append(error("YD01"))

    if not doc.see_also:
        wrns.append(error("SA01"))
    else:
        for rel_name, rel_desc in doc.see_also.items():
            if rel_desc:
                if not rel_desc.endswith("."):
                    errs.append(error("SA02", reference_name=rel_name))
                if not rel_desc[0].isupper():
                    errs.append(error("SA03", reference_name=rel_name))
            else:
                errs.append(error("SA04", reference_name=rel_name))
            if rel_name.startswith("pandas."):
                errs.append(
                    error(
                        "SA05",
                        reference_name=rel_name,
                        right_reference=rel_name[len("pandas.") :],
                    )
                )

    examples_errs = ""
    if not doc.examples:
        wrns.append(error("EX01"))
    else:
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
    return errs, wrns, examples_errs


def validate_one(func_name):
    """
    Validate the docstring for the given func_name

    Parameters
    ----------
    func_name : function
        Function whose docstring will be evaluated (e.g. pandas.read_csv).

    Returns
    -------
    dict
        A dictionary containing all the information obtained from validating
        the docstring.
    """
    doc = Docstring(func_name)
    errs, wrns, examples_errs = get_validation_data(doc)
    return {
        "type": doc.type,
        "docstring": doc.clean_doc,
        "deprecated": doc.deprecated,
        "file": doc.source_file_name,
        "file_line": doc.source_file_def_line,
        "github_link": doc.github_url,
        "errors": errs,
        "warnings": wrns,
        "examples_errors": examples_errs,
    }


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
                raise ValueError(f'Unknown output_format "{output_format}"')

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
                        path=res["file"],
                        row=res["file_line"],
                        code=err_code,
                        text=f"{name}: {err_desc}",
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
