import io
import textwrap

import pytest

from scripts import validate_docstrings


class BadDocstrings:
    """Everything here has a bad docstring"""

    def private_classes(self) -> None:
        """
        This mentions NDFrame, which is not correct.
        """

    def prefix_pandas(self) -> None:
        """
        Have `pandas` prefix in See Also section.

        See Also
        --------
        pandas.Series.rename : Alter Series index labels or name.
        DataFrame.head : The first `n` rows of the caller object.
        """

    def redundant_import(self, paramx=None, paramy=None) -> None:
        """
        A sample DataFrame method.

        Should not import numpy and pandas.

        Examples
        --------
        >>> import numpy as np
        >>> import pandas as pd
        >>> df = pd.DataFrame(np.ones((3, 3)),
        ...                   columns=('a', 'b', 'c'))
        >>> df.all(1)
        0    True
        1    True
        2    True
        dtype: bool
        >>> df.all(bool_only=True)
        Series([], dtype: bool)
        """

    def unused_import(self) -> None:
        """
        Examples
        --------
        >>> import pandas as pdf
        >>> df = pd.DataFrame(np.ones((3, 3)), columns=('a', 'b', 'c'))
        """

    def missing_whitespace_around_arithmetic_operator(self) -> None:
        """
        Examples
        --------
        >>> 2+5
        7
        """

    def indentation_is_not_a_multiple_of_four(self) -> None:
        """
        Examples
        --------
        >>> if 2 + 5:
        ...   pass
        """

    def missing_whitespace_after_comma(self) -> None:
        """
        Examples
        --------
        >>> df = pd.DataFrame(np.ones((3,3)),columns=('a','b', 'c'))
        """

    def write_array_like_with_hyphen_not_underscore(self) -> None:
        """
        In docstrings, use array-like over array_like
        """

    def leftover_files(self) -> None:
        """
        Examples
        --------
        >>> import pathlib
        >>> pathlib.Path("foo.txt").touch()
        """


class TestValidator:
    def _import_path(self, klass=None, func=None):
        """
        Build the required import path for tests in this module.

        Parameters
        ----------
        klass : str
            Class name of object in module.
        func : str
            Function name of object in module.

        Returns
        -------
        str
            Import path of specified object in this module
        """
        base_path = "scripts.tests.test_validate_docstrings"

        if klass:
            base_path = f"{base_path}.{klass}"

        if func:
            base_path = f"{base_path}.{func}"

        return base_path

    def test_bad_class(self, capsys) -> None:
        errors = validate_docstrings.pandas_validate(
            self._import_path(klass="BadDocstrings")
        )["errors"]
        assert isinstance(errors, list)
        assert errors

    @pytest.mark.parametrize(
        "klass,func,msgs",
        [
            (
                "BadDocstrings",
                "private_classes",
                (
                    "Private classes (NDFrame) should not be mentioned in public "
                    "docstrings",
                ),
            ),
            (
                "BadDocstrings",
                "prefix_pandas",
                (
                    "pandas.Series.rename in `See Also` section "
                    "does not need `pandas` prefix",
                ),
            ),
            # Examples tests
            (
                "BadDocstrings",
                "redundant_import",
                ("Do not import numpy, as it is imported automatically",),
            ),
            (
                "BadDocstrings",
                "redundant_import",
                ("Do not import pandas, as it is imported automatically",),
            ),
            (
                "BadDocstrings",
                "unused_import",
                (
                    "flake8 error: line 1, col 1: F401 'pandas as pdf' "
                    "imported but unused",
                ),
            ),
            (
                "BadDocstrings",
                "missing_whitespace_around_arithmetic_operator",
                (
                    "flake8 error: line 1, col 2: "
                    "E226 missing whitespace around arithmetic operator",
                ),
            ),
            (
                "BadDocstrings",
                "indentation_is_not_a_multiple_of_four",
                # with flake8 3.9.0, the message ends with four spaces,
                #  whereas in earlier versions, it ended with "four"
                (
                    "flake8 error: line 2, col 3: E111 indentation is not a "
                    "multiple of 4",
                ),
            ),
            (
                "BadDocstrings",
                "missing_whitespace_after_comma",
                ("flake8 error: line 1, col 33: E231 missing whitespace after ','",),
            ),
            (
                "BadDocstrings",
                "write_array_like_with_hyphen_not_underscore",
                ("Use 'array-like' rather than 'array_like' in docstrings",),
            ),
        ],
    )
    def test_bad_docstrings(self, capsys, klass, func, msgs) -> None:
        result = validate_docstrings.pandas_validate(
            self._import_path(klass=klass, func=func)
        )
        for msg in msgs:
            assert msg in " ".join([err[1] for err in result["errors"]])

    def test_validate_all_ignore_deprecated(self, monkeypatch) -> None:
        monkeypatch.setattr(
            validate_docstrings,
            "pandas_validate",
            lambda func_name: {
                "docstring": "docstring1",
                "errors": [
                    ("ER01", "err desc"),
                    ("ER02", "err desc"),
                    ("ER03", "err desc"),
                ],
                "warnings": [],
                "examples_errors": "",
                "deprecated": True,
            },
        )
        result = validate_docstrings.validate_all(prefix=None, ignore_deprecated=True)
        assert len(result) == 0

    def test_validate_all_ignore_errors(self, monkeypatch):
        monkeypatch.setattr(
            validate_docstrings,
            "pandas_validate",
            lambda func_name: {
                "docstring": "docstring1",
                "errors": [
                    ("ER01", "err desc"),
                    ("ER02", "err desc"),
                    ("ER03", "err desc")
                ],
                "warnings": [],
                "examples_errors": "",
                "deprecated": True,
                "file": "file1",
                "file_line": "file_line1"
            },
        )
        monkeypatch.setattr(
            validate_docstrings,
            "get_all_api_items",
            lambda: [
                (
                    "pandas.DataFrame.align",
                    "func",
                    "current_section",
                    "current_subsection",
                ),
                (
                    "pandas.Index.all",
                    "func",
                    "current_section",
                    "current_subsection",
                ),
            ],
        )

        exit_status = validate_docstrings.print_validate_all_results(
            output_format="default",
            prefix=None,
            ignore_deprecated=False,
            ignore_errors={None: {"ER03"}},
        )
        # two functions * two not ignored errors
        assert exit_status == 2 * 2

        exit_status = validate_docstrings.print_validate_all_results(
            output_format="default",
            prefix=None,
            ignore_deprecated=False,
            ignore_errors={
                None: {"ER03"},
                "pandas.DataFrame.align": {"ER01"},
                # ignoring an error that is not requested should be of no effect
                "pandas.Index.all": {"ER03"}
            }
        )
        # two functions * two not global ignored errors - one function ignored error
        assert exit_status == 2 * 2 - 1



class TestApiItems:
    @property
    def api_doc(self):
        return io.StringIO(
            textwrap.dedent(
                """
            .. currentmodule:: itertools

            Itertools
            ---------

            Infinite
            ~~~~~~~~

            .. autosummary::

                cycle
                count

            Finite
            ~~~~~~

            .. autosummary::

                chain

            .. currentmodule:: random

            Random
            ------

            All
            ~~~

            .. autosummary::

                seed
                randint
            """
            )
        )

    @pytest.mark.parametrize(
        "idx,name",
        [
            (0, "itertools.cycle"),
            (1, "itertools.count"),
            (2, "itertools.chain"),
            (3, "random.seed"),
            (4, "random.randint"),
        ],
    )
    def test_item_name(self, idx, name) -> None:
        result = list(validate_docstrings.get_api_items(self.api_doc))
        assert result[idx][0] == name

    @pytest.mark.parametrize(
        "idx,func",
        [(0, "cycle"), (1, "count"), (2, "chain"), (3, "seed"), (4, "randint")],
    )
    def test_item_function(self, idx, func) -> None:
        result = list(validate_docstrings.get_api_items(self.api_doc))
        assert callable(result[idx][1])
        assert result[idx][1].__name__ == func

    @pytest.mark.parametrize(
        "idx,section",
        [
            (0, "Itertools"),
            (1, "Itertools"),
            (2, "Itertools"),
            (3, "Random"),
            (4, "Random"),
        ],
    )
    def test_item_section(self, idx, section) -> None:
        result = list(validate_docstrings.get_api_items(self.api_doc))
        assert result[idx][2] == section

    @pytest.mark.parametrize(
        "idx,subsection",
        [(0, "Infinite"), (1, "Infinite"), (2, "Finite"), (3, "All"), (4, "All")],
    )
    def test_item_subsection(self, idx, subsection) -> None:
        result = list(validate_docstrings.get_api_items(self.api_doc))
        assert result[idx][3] == subsection


class TestPandasDocstringClass:
    @pytest.mark.parametrize(
        "name", ["pandas.Series.str.isdecimal", "pandas.Series.str.islower"]
    )
    def test_encode_content_write_to_file(self, name) -> None:
        # GH25466
        docstr = validate_docstrings.PandasDocstring(name).validate_pep8()
        # the list of pep8 errors should be empty
        assert not list(docstr)


class TestMainFunction:
    def test_exit_status_for_main(self, monkeypatch) -> None:
        monkeypatch.setattr(
            validate_docstrings,
            "pandas_validate",
            lambda func_name: {
                "docstring": "docstring1",
                "errors": [
                    ("ER01", "err desc"),
                    ("ER02", "err desc"),
                    ("ER03", "err desc"),
                ],
                "examples_errs": "",
            },
        )
        exit_status = validate_docstrings.main(
            func_name="docstring1",
            prefix=None,
            output_format="default",
            ignore_deprecated=False,
            ignore_errors={},
        )
        assert exit_status == 3

    def test_exit_status_errors_for_validate_all(self, monkeypatch) -> None:
        monkeypatch.setattr(
            validate_docstrings,
            "validate_all",
            lambda prefix, ignore_deprecated=False, ignore_functions=None: {
                "docstring1": {
                    "errors": [
                        ("ER01", "err desc"),
                        ("ER02", "err desc"),
                        ("ER03", "err desc"),
                    ],
                    "file": "module1.py",
                    "file_line": 23,
                },
                "docstring2": {
                    "errors": [("ER04", "err desc"), ("ER05", "err desc")],
                    "file": "module2.py",
                    "file_line": 925,
                },
            },
        )
        exit_status = validate_docstrings.main(
            func_name=None,
            prefix=None,
            output_format="default",
            ignore_deprecated=False,
            ignore_errors={},
        )
        assert exit_status == 5

    def test_no_exit_status_noerrors_for_validate_all(self, monkeypatch) -> None:
        monkeypatch.setattr(
            validate_docstrings,
            "validate_all",
            lambda prefix, ignore_deprecated=False, ignore_functions=None: {
                "docstring1": {"errors": [], "warnings": [("WN01", "warn desc")]},
                "docstring2": {"errors": []},
            },
        )
        exit_status = validate_docstrings.main(
            func_name=None,
            output_format="default",
            prefix=None,
            ignore_deprecated=False,
            ignore_errors={},
        )
        assert exit_status == 0

    def test_exit_status_for_validate_all_json(self, monkeypatch) -> None:
        monkeypatch.setattr(
            validate_docstrings,
            "validate_all",
            lambda prefix, ignore_deprecated=False, ignore_functions=None: {
                "docstring1": {
                    "errors": [
                        ("ER01", "err desc"),
                        ("ER02", "err desc"),
                        ("ER03", "err desc"),
                    ]
                },
                "docstring2": {"errors": [("ER04", "err desc"), ("ER05", "err desc")]},
            },
        )
        exit_status = validate_docstrings.main(
            func_name=None,
            output_format="json",
            prefix=None,
            ignore_deprecated=False,
            ignore_errors={},
        )
        assert exit_status == 0

    def test_errors_param_filters_errors(self, monkeypatch) -> None:
        monkeypatch.setattr(
            validate_docstrings,
            "validate_all",
            lambda prefix, ignore_deprecated=False, ignore_functions=None: {
                "Series.foo": {
                    "errors": [
                        ("ER01", "err desc"),
                        ("ER02", "err desc"),
                        ("ER03", "err desc"),
                    ],
                    "file": "series.py",
                    "file_line": 142,
                },
                "DataFrame.bar": {
                    "errors": [("ER01", "err desc"), ("ER02", "err desc")],
                    "file": "frame.py",
                    "file_line": 598,
                },
                "Series.foobar": {
                    "errors": [("ER01", "err desc")],
                    "file": "series.py",
                    "file_line": 279,
                },
            },
        )
        monkeypatch.setattr(
            validate_docstrings,
            "ERROR_MSGS",
            {
                "ER01": "err desc",
                "ER02": "err desc",
                "ER03": "err desc",
            },
        )
        exit_status = validate_docstrings.main(
            func_name=None,
            output_format="default",
            prefix=None,
            ignore_deprecated=False,
            ignore_errors={None: {"ER02", "ER03"}},
        )
        assert exit_status == 3

        exit_status = validate_docstrings.main(
            func_name=None,
            output_format="default",
            prefix=None,
            ignore_deprecated=False,
            ignore_errors={None: {"ER01", "ER02"}},
        )
        assert exit_status == 1
