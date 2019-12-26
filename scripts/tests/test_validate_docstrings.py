import functools
import io
import random
import string
import textwrap

import numpy as np
import pytest
import validate_docstrings

import pandas as pd

validate_one = validate_docstrings.validate_one


class GoodDocStrings:
    """
    Collection of good doc strings.

    This class contains a lot of docstrings that should pass the validation
    script without any errors.
    """

    def plot(self, kind, color="blue", **kwargs):
        """
        Generate a plot.

        Render the data in the Series as a matplotlib plot of the
        specified kind.

        Parameters
        ----------
        kind : str
            Kind of matplotlib plot.
        color : str, default 'blue'
            Color name or rgb code.
        **kwargs
            These parameters will be passed to the matplotlib plotting
            function.
        """
        pass

    def swap(self, arr, i, j, *args, **kwargs):
        """
        Swap two indicies on an array.

        Parameters
        ----------
        arr : list
            The list having indexes swapped.
        i, j : int
            The indexes being swapped.
        *args, **kwargs
            Extraneous parameters are being permitted.
        """
        pass

    def sample(self):
        """
        Generate and return a random number.

        The value is sampled from a continuous uniform distribution between
        0 and 1.

        Returns
        -------
        float
            Random number generated.
        """
        return random.random()

    @functools.lru_cache(None)
    def decorated_sample(self, max):
        """
        Generate and return a random integer between 0 and max.

        Parameters
        ----------
        max : int
            The maximum value of the random number.

        Returns
        -------
        int
            Random number generated.
        """
        return random.randint(0, max)

    def random_letters(self):
        """
        Generate and return a sequence of random letters.

        The length of the returned string is also random, and is also
        returned.

        Returns
        -------
        length : int
            Length of the returned string.
        letters : str
            String of random letters.
        """
        length = random.randint(1, 10)
        letters = "".join(random.sample(string.ascii_lowercase, length))
        return length, letters

    def sample_values(self):
        """
        Generate an infinite sequence of random numbers.

        The values are sampled from a continuous uniform distribution between
        0 and 1.

        Yields
        ------
        float
            Random number generated.
        """
        while True:
            yield random.random()

    def head(self):
        """
        Return the first 5 elements of the Series.

        This function is mainly useful to preview the values of the
        Series without displaying the whole of it.

        Returns
        -------
        Series
            Subset of the original series with the 5 first values.

        See Also
        --------
        Series.tail : Return the last 5 elements of the Series.
        Series.iloc : Return a slice of the elements in the Series,
            which can also be used to return the first or last n.
        """
        return self.iloc[:5]

    def head1(self, n=5):
        """
        Return the first elements of the Series.

        This function is mainly useful to preview the values of the
        Series without displaying the whole of it.

        Parameters
        ----------
        n : int
            Number of values to return.

        Returns
        -------
        Series
            Subset of the original series with the n first values.

        See Also
        --------
        tail : Return the last n elements of the Series.

        Examples
        --------
        >>> s = pd.Series(['Ant', 'Bear', 'Cow', 'Dog', 'Falcon'])
        >>> s.head()
        0   Ant
        1   Bear
        2   Cow
        3   Dog
        4   Falcon
        dtype: object

        With the `n` parameter, we can change the number of returned rows:

        >>> s.head(n=3)
        0   Ant
        1   Bear
        2   Cow
        dtype: object
        """
        return self.iloc[:n]

    def contains(self, pat, case=True, na=np.nan):
        """
        Return whether each value contains `pat`.

        In this case, we are illustrating how to use sections, even
        if the example is simple enough and does not require them.

        Parameters
        ----------
        pat : str
            Pattern to check for within each element.
        case : bool, default True
            Whether check should be done with case sensitivity.
        na : object, default np.nan
            Fill value for missing data.

        Examples
        --------
        >>> s = pd.Series(['Antelope', 'Lion', 'Zebra', np.nan])
        >>> s.str.contains(pat='a')
        0    False
        1    False
        2     True
        3      NaN
        dtype: object

        **Case sensitivity**

        With `case_sensitive` set to `False` we can match `a` with both
        `a` and `A`:

        >>> s.str.contains(pat='a', case=False)
        0     True
        1    False
        2     True
        3      NaN
        dtype: object

        **Missing values**

        We can fill missing values in the output using the `na` parameter:

        >>> s.str.contains(pat='a', na=False)
        0    False
        1    False
        2     True
        3    False
        dtype: bool
        """
        pass

    def mode(self, axis, numeric_only):
        """
        Ensure reST directives don't affect checks for leading periods.

        Parameters
        ----------
        axis : str
            Sentence ending in period, followed by single directive.

            .. versionchanged:: 0.1.2

        numeric_only : bool
            Sentence ending in period, followed by multiple directives.

            .. versionadded:: 0.1.2
            .. deprecated:: 0.00.0
                A multiline description,
                which spans another line.
        """
        pass

    def good_imports(self):
        """
        Ensure import other than numpy and pandas are fine.

        Examples
        --------
        This example does not import pandas or import numpy.
        >>> import datetime
        >>> datetime.MAXYEAR
        9999
        """
        pass

    def no_returns(self):
        """
        Say hello and have no returns.
        """
        pass

    def empty_returns(self):
        """
        Say hello and always return None.

        Since this function never returns a value, this
        docstring doesn't need a return section.
        """

        def say_hello():
            return "Hello World!"

        say_hello()
        if True:
            return
        else:
            return None

    def multiple_variables_on_one_line(self, matrix, a, b, i, j):
        """
        Swap two values in a matrix.

        Parameters
        ----------
        matrix : list of list
            A double list that represents a matrix.
        a, b : int
            The indicies of the first value.
        i, j : int
            The indicies of the second value.
        """
        pass


class BadGenericDocStrings:
    """Everything here has a bad docstring
    """

    def func(self):

        """Some function.

        With several mistakes in the docstring.

        It has a blank like after the signature `def func():`.

        The text 'Some function' should go in the line after the
        opening quotes of the docstring, not in the same line.

        There is a blank line between the docstring and the first line
        of code `foo = 1`.

        The closing quotes should be in the next line, not in this one."""

        foo = 1
        bar = 2
        return foo + bar

    def astype(self, dtype):
        """
        Casts Series type.

        Verb in third-person of the present simple, should be infinitive.
        """
        pass

    def astype1(self, dtype):
        """
        Method to cast Series type.

        Does not start with verb.
        """
        pass

    def astype2(self, dtype):
        """
        Cast Series type

        Missing dot at the end.
        """
        pass

    def astype3(self, dtype):
        """
        Cast Series type from its current type to the new type defined in
        the parameter dtype.

        Summary is too verbose and doesn't fit in a single line.
        """
        pass

    def two_linebreaks_between_sections(self, foo):
        """
        Test linebreaks message GL03.

        Note 2 blank lines before parameters section.


        Parameters
        ----------
        foo : str
            Description of foo parameter.
        """
        pass

    def linebreak_at_end_of_docstring(self, foo):
        """
        Test linebreaks message GL03.

        Note extra blank line at end of docstring.

        Parameters
        ----------
        foo : str
            Description of foo parameter.

        """
        pass

    def plot(self, kind, **kwargs):
        """
        Generate a plot.

        Render the data in the Series as a matplotlib plot of the
        specified kind.

        Note the blank line between the parameters title and the first
        parameter. Also, note that after the name of the parameter `kind`
        and before the colon, a space is missing.

        Also, note that the parameter descriptions do not start with a
        capital letter, and do not finish with a dot.

        Finally, the `**kwargs` parameter is missing.

        Parameters
        ----------

        kind: str
            kind of matplotlib plot
        """
        pass

    def method(self, foo=None, bar=None):
        """
        A sample DataFrame method.

        Do not import numpy and pandas.

        Try to use meaningful data, when it makes the example easier
        to understand.

        Try to avoid positional arguments like in `df.method(1)`. They
        can be alright if previously defined with a meaningful name,
        like in `present_value(interest_rate)`, but avoid them otherwise.

        When presenting the behavior with different parameters, do not place
        all the calls one next to the other. Instead, add a short sentence
        explaining what the example shows.

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
        pass

    def private_classes(self):
        """
        This mentions NDFrame, which is not correct.
        """

    def unknown_section(self):
        """
        This section has an unknown section title.

        Unknown Section
        ---------------
        This should raise an error in the validation.
        """

    def sections_in_wrong_order(self):
        """
        This docstring has the sections in the wrong order.

        Parameters
        ----------
        name : str
            This section is in the right position.

        Examples
        --------
        >>> print('So far Examples is good, as it goes before Parameters')
        So far Examples is good, as it goes before Parameters

        See Also
        --------
        function : This should generate an error, as See Also needs to go
            before Examples.
        """

    def deprecation_in_wrong_order(self):
        """
        This docstring has the deprecation warning in the wrong order.

        This is the extended summary. The correct order should be
        summary, deprecation warning, extended summary.

        .. deprecated:: 1.0
            This should generate an error as it needs to go before
            extended summary.
        """

    def method_wo_docstrings(self):
        pass

    def directives_without_two_colons(self, first, second):
        """
        Ensure reST directives have trailing colons.

        Parameters
        ----------
        first : str
            Sentence ending in period, followed by single directive w/o colons.

            .. versionchanged 0.1.2

        second : bool
            Sentence ending in period, followed by multiple directives w/o
            colons.

            .. versionadded 0.1.2
            .. deprecated 0.00.0

        """
        pass


class BadSummaries:
    def wrong_line(self):
        """Exists on the wrong line"""
        pass

    def no_punctuation(self):
        """
        Has the right line but forgets punctuation
        """
        pass

    def no_capitalization(self):
        """
        provides a lowercase summary.
        """
        pass

    def no_infinitive(self):
        """
        Started with a verb that is not infinitive.
        """

    def multi_line(self):
        """
        Extends beyond one line
        which is not correct.
        """

    def two_paragraph_multi_line(self):
        """
        Extends beyond one line
        which is not correct.

        Extends beyond one line, which in itself is correct but the
        previous short summary should still be an issue.
        """


class BadParameters:
    """
    Everything here has a problem with its Parameters section.
    """

    def missing_params(self, kind, **kwargs):
        """
        Lacks kwargs in Parameters.

        Parameters
        ----------
        kind : str
            Foo bar baz.
        """

    def bad_colon_spacing(self, kind):
        """
        Has bad spacing in the type line.

        Parameters
        ----------
        kind: str
            Needs a space after kind.
        """

    def no_description_period(self, kind):
        """
        Forgets to add a period to the description.

        Parameters
        ----------
        kind : str
           Doesn't end with a dot
        """

    def no_description_period_with_directive(self, kind):
        """
        Forgets to add a period, and also includes a directive.

        Parameters
        ----------
        kind : str
           Doesn't end with a dot

           .. versionadded:: 0.00.0
        """

    def no_description_period_with_directives(self, kind):
        """
        Forgets to add a period, and also includes multiple directives.

        Parameters
        ----------
        kind : str
           Doesn't end with a dot

           .. versionchanged:: 0.00.0
           .. deprecated:: 0.00.0
        """

    def parameter_capitalization(self, kind):
        """
        Forgets to capitalize the description.

        Parameters
        ----------
        kind : str
           this is not capitalized.
        """

    def blank_lines(self, kind):
        """
        Adds a blank line after the section header.

        Parameters
        ----------

        kind : str
            Foo bar baz.
        """
        pass

    def integer_parameter(self, kind):
        """
        Uses integer instead of int.

        Parameters
        ----------
        kind : integer
            Foo bar baz.
        """
        pass

    def string_parameter(self, kind):
        """
        Uses string instead of str.

        Parameters
        ----------
        kind : string
            Foo bar baz.
        """
        pass

    def boolean_parameter(self, kind):
        """
        Uses boolean instead of bool.

        Parameters
        ----------
        kind : boolean
            Foo bar baz.
        """
        pass

    def list_incorrect_parameter_type(self, kind):
        """
        Uses list of boolean instead of list of bool.

        Parameters
        ----------
        kind : list of boolean, integer, float or string
            Foo bar baz.
        """
        pass

    def bad_parameter_spacing(self, a, b):
        """
        The parameters on the same line have an extra space between them.

        Parameters
        ----------
        a,  b : int
            Foo bar baz.
        """
        pass


class BadReturns:
    def return_not_documented(self):
        """
        Lacks section for Returns
        """
        return "Hello world!"

    def yield_not_documented(self):
        """
        Lacks section for Yields
        """
        yield "Hello world!"

    def no_type(self):
        """
        Returns documented but without type.

        Returns
        -------
        Some value.
        """
        return "Hello world!"

    def no_description(self):
        """
        Provides type but no description.

        Returns
        -------
        str
        """
        return "Hello world!"

    def no_punctuation(self):
        """
        Provides type and description but no period.

        Returns
        -------
        str
           A nice greeting
        """
        return "Hello world!"

    def named_single_return(self):
        """
        Provides name but returns only one value.

        Returns
        -------
        s : str
           A nice greeting.
        """
        return "Hello world!"

    def no_capitalization(self):
        """
        Forgets capitalization in return values description.

        Returns
        -------
        foo : str
           The first returned string.
        bar : str
           the second returned string.
        """
        return "Hello", "World!"

    def no_period_multi(self):
        """
        Forgets period in return values description.

        Returns
        -------
        foo : str
           The first returned string
        bar : str
           The second returned string.
        """
        return "Hello", "World!"


class BadSeeAlso:
    def desc_no_period(self):
        """
        Return the first 5 elements of the Series.

        See Also
        --------
        Series.tail : Return the last 5 elements of the Series.
        Series.iloc : Return a slice of the elements in the Series,
            which can also be used to return the first or last n
        """
        pass

    def desc_first_letter_lowercase(self):
        """
        Return the first 5 elements of the Series.

        See Also
        --------
        Series.tail : return the last 5 elements of the Series.
        Series.iloc : Return a slice of the elements in the Series,
            which can also be used to return the first or last n.
        """
        pass

    def prefix_pandas(self):
        """
        Have `pandas` prefix in See Also section.

        See Also
        --------
        pandas.Series.rename : Alter Series index labels or name.
        DataFrame.head : The first `n` rows of the caller object.
        """
        pass


class BadExamples:
    def unused_import(self):
        """
        Examples
        --------
        >>> import pandas as pdf
        >>> df = pd.DataFrame(np.ones((3, 3)), columns=('a', 'b', 'c'))
        """
        pass

    def missing_whitespace_around_arithmetic_operator(self):
        """
        Examples
        --------
        >>> 2+5
        7
        """
        pass

    def indentation_is_not_a_multiple_of_four(self):
        """
        Examples
        --------
        >>> if 2 + 5:
        ...   pass
        """
        pass

    def missing_whitespace_after_comma(self):
        """
        Examples
        --------
        >>> df = pd.DataFrame(np.ones((3,3)),columns=('a','b', 'c'))
        """
        pass


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
            base_path = ".".join([base_path, klass])

        if func:
            base_path = ".".join([base_path, func])

        return base_path

    def test_good_class(self, capsys):
        errors = validate_one(self._import_path(klass="GoodDocStrings"))["errors"]
        assert isinstance(errors, list)
        assert not errors

    @pytest.mark.parametrize(
        "func",
        [
            "plot",
            "swap",
            "sample",
            "decorated_sample",
            "random_letters",
            "sample_values",
            "head",
            "head1",
            "contains",
            "mode",
            "good_imports",
            "no_returns",
            "empty_returns",
            "multiple_variables_on_one_line",
        ],
    )
    def test_good_functions(self, capsys, func):
        errors = validate_one(self._import_path(klass="GoodDocStrings", func=func))[
            "errors"
        ]
        assert isinstance(errors, list)
        assert not errors

    def test_bad_class(self, capsys):
        errors = validate_one(self._import_path(klass="BadGenericDocStrings"))["errors"]
        assert isinstance(errors, list)
        assert errors

    @pytest.mark.parametrize(
        "func",
        [
            "func",
            "astype",
            "astype1",
            "astype2",
            "astype3",
            "plot",
            "method",
            "private_classes",
            "directives_without_two_colons",
        ],
    )
    def test_bad_generic_functions(self, capsys, func):
        errors = validate_one(
            self._import_path(klass="BadGenericDocStrings", func=func)  # noqa:F821
        )["errors"]
        assert isinstance(errors, list)
        assert errors

    @pytest.mark.parametrize(
        "klass,func,msgs",
        [
            # See Also tests
            (
                "BadGenericDocStrings",
                "private_classes",
                (
                    "Private classes (NDFrame) should not be mentioned in public "
                    "docstrings",
                ),
            ),
            (
                "BadGenericDocStrings",
                "unknown_section",
                ('Found unknown section "Unknown Section".',),
            ),
            (
                "BadGenericDocStrings",
                "sections_in_wrong_order",
                (
                    "Sections are in the wrong order. Correct order is: Parameters, "
                    "See Also, Examples",
                ),
            ),
            (
                "BadGenericDocStrings",
                "deprecation_in_wrong_order",
                ("Deprecation warning should precede extended summary",),
            ),
            (
                "BadGenericDocStrings",
                "directives_without_two_colons",
                (
                    "reST directives ['versionchanged', 'versionadded', "
                    "'deprecated'] must be followed by two colons",
                ),
            ),
            (
                "BadSeeAlso",
                "desc_no_period",
                ('Missing period at end of description for See Also "Series.iloc"',),
            ),
            (
                "BadSeeAlso",
                "desc_first_letter_lowercase",
                ('should be capitalized for See Also "Series.tail"',),
            ),
            # Summary tests
            (
                "BadSummaries",
                "wrong_line",
                ("should start in the line immediately after the opening quotes",),
            ),
            ("BadSummaries", "no_punctuation", ("Summary does not end with a period",)),
            (
                "BadSummaries",
                "no_capitalization",
                ("Summary does not start with a capital letter",),
            ),
            (
                "BadSummaries",
                "no_capitalization",
                ("Summary must start with infinitive verb",),
            ),
            ("BadSummaries", "multi_line", ("Summary should fit in a single line",)),
            (
                "BadSummaries",
                "two_paragraph_multi_line",
                ("Summary should fit in a single line",),
            ),
            # Parameters tests
            (
                "BadParameters",
                "missing_params",
                ("Parameters {**kwargs} not documented",),
            ),
            (
                "BadParameters",
                "bad_colon_spacing",
                (
                    'Parameter "kind" requires a space before the colon '
                    "separating the parameter name and type",
                ),
            ),
            (
                "BadParameters",
                "no_description_period",
                ('Parameter "kind" description should finish with "."',),
            ),
            (
                "BadParameters",
                "no_description_period_with_directive",
                ('Parameter "kind" description should finish with "."',),
            ),
            (
                "BadParameters",
                "parameter_capitalization",
                ('Parameter "kind" description should start with a capital letter',),
            ),
            (
                "BadParameters",
                "integer_parameter",
                ('Parameter "kind" type should use "int" instead of "integer"',),
            ),
            (
                "BadParameters",
                "string_parameter",
                ('Parameter "kind" type should use "str" instead of "string"',),
            ),
            (
                "BadParameters",
                "boolean_parameter",
                ('Parameter "kind" type should use "bool" instead of "boolean"',),
            ),
            (
                "BadParameters",
                "list_incorrect_parameter_type",
                ('Parameter "kind" type should use "bool" instead of "boolean"',),
            ),
            (
                "BadParameters",
                "list_incorrect_parameter_type",
                ('Parameter "kind" type should use "int" instead of "integer"',),
            ),
            (
                "BadParameters",
                "list_incorrect_parameter_type",
                ('Parameter "kind" type should use "str" instead of "string"',),
            ),
            (
                "BadParameters",
                "bad_parameter_spacing",
                ("Parameters {b} not documented", "Unknown parameters { b}"),
            ),
            pytest.param(
                "BadParameters",
                "blank_lines",
                ("No error yet?",),
                marks=pytest.mark.xfail,
            ),
            # Returns tests
            ("BadReturns", "return_not_documented", ("No Returns section found",)),
            ("BadReturns", "yield_not_documented", ("No Yields section found",)),
            pytest.param("BadReturns", "no_type", ("foo",), marks=pytest.mark.xfail),
            ("BadReturns", "no_description", ("Return value has no description",)),
            (
                "BadReturns",
                "no_punctuation",
                ('Return value description should finish with "."',),
            ),
            (
                "BadReturns",
                "named_single_return",
                (
                    "The first line of the Returns section should contain only the "
                    "type, unless multiple values are being returned",
                ),
            ),
            (
                "BadReturns",
                "no_capitalization",
                ("Return value description should start with a capital letter",),
            ),
            (
                "BadReturns",
                "no_period_multi",
                ('Return value description should finish with "."',),
            ),
            # Examples tests
            (
                "BadGenericDocStrings",
                "method",
                ("Do not import numpy, as it is imported automatically",),
            ),
            (
                "BadGenericDocStrings",
                "method",
                ("Do not import pandas, as it is imported automatically",),
            ),
            (
                "BadGenericDocStrings",
                "method_wo_docstrings",
                ("The object does not have a docstring",),
            ),
            # See Also tests
            (
                "BadSeeAlso",
                "prefix_pandas",
                (
                    "pandas.Series.rename in `See Also` section "
                    "does not need `pandas` prefix",
                ),
            ),
            # Examples tests
            (
                "BadExamples",
                "unused_import",
                ("flake8 error: F401 'pandas as pdf' imported but unused",),
            ),
            (
                "BadExamples",
                "indentation_is_not_a_multiple_of_four",
                ("flake8 error: E111 indentation is not a multiple of four",),
            ),
            (
                "BadExamples",
                "missing_whitespace_around_arithmetic_operator",
                (
                    "flake8 error: "
                    "E226 missing whitespace around arithmetic operator",
                ),
            ),
            (
                "BadExamples",
                "missing_whitespace_after_comma",
                ("flake8 error: E231 missing whitespace after ',' (3 times)",),
            ),
            (
                "BadGenericDocStrings",
                "two_linebreaks_between_sections",
                (
                    "Double line break found; please use only one blank line to "
                    "separate sections or paragraphs, and do not leave blank lines "
                    "at the end of docstrings",
                ),
            ),
            (
                "BadGenericDocStrings",
                "linebreak_at_end_of_docstring",
                (
                    "Double line break found; please use only one blank line to "
                    "separate sections or paragraphs, and do not leave blank lines "
                    "at the end of docstrings",
                ),
            ),
        ],
    )
    def test_bad_docstrings(self, capsys, klass, func, msgs):
        result = validate_one(self._import_path(klass=klass, func=func))
        for msg in msgs:
            assert msg in " ".join(err[1] for err in result["errors"])

    def test_validate_all_ignore_deprecated(self, monkeypatch):
        monkeypatch.setattr(
            validate_docstrings,
            "validate_one",
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
    def test_item_name(self, idx, name):
        result = list(validate_docstrings.get_api_items(self.api_doc))
        assert result[idx][0] == name

    @pytest.mark.parametrize(
        "idx,func",
        [(0, "cycle"), (1, "count"), (2, "chain"), (3, "seed"), (4, "randint")],
    )
    def test_item_function(self, idx, func):
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
    def test_item_section(self, idx, section):
        result = list(validate_docstrings.get_api_items(self.api_doc))
        assert result[idx][2] == section

    @pytest.mark.parametrize(
        "idx,subsection",
        [(0, "Infinite"), (1, "Infinite"), (2, "Finite"), (3, "All"), (4, "All")],
    )
    def test_item_subsection(self, idx, subsection):
        result = list(validate_docstrings.get_api_items(self.api_doc))
        assert result[idx][3] == subsection


class TestDocstringClass:
    @pytest.mark.parametrize(
        "name, expected_obj",
        [
            ("pandas.isnull", pd.isnull),
            ("pandas.DataFrame", pd.DataFrame),
            ("pandas.Series.sum", pd.Series.sum),
        ],
    )
    def test_resolves_class_name(self, name, expected_obj):
        d = validate_docstrings.Docstring(name)
        assert d.obj is expected_obj

    @pytest.mark.parametrize("invalid_name", ["panda", "panda.DataFrame"])
    def test_raises_for_invalid_module_name(self, invalid_name):
        msg = 'No module can be imported from "{}"'.format(invalid_name)
        with pytest.raises(ImportError, match=msg):
            validate_docstrings.Docstring(invalid_name)

    @pytest.mark.parametrize(
        "invalid_name", ["pandas.BadClassName", "pandas.Series.bad_method_name"]
    )
    def test_raises_for_invalid_attribute_name(self, invalid_name):
        name_components = invalid_name.split(".")
        obj_name, invalid_attr_name = name_components[-2], name_components[-1]
        msg = "'{}' has no attribute '{}'".format(obj_name, invalid_attr_name)
        with pytest.raises(AttributeError, match=msg):
            validate_docstrings.Docstring(invalid_name)

    @pytest.mark.parametrize(
        "name", ["pandas.Series.str.isdecimal", "pandas.Series.str.islower"]
    )
    def test_encode_content_write_to_file(self, name):
        # GH25466
        docstr = validate_docstrings.Docstring(name).validate_pep8()
        # the list of pep8 errors should be empty
        assert not list(docstr)


class TestMainFunction:
    def test_exit_status_for_validate_one(self, monkeypatch):
        monkeypatch.setattr(
            validate_docstrings,
            "validate_one",
            lambda func_name: {
                "docstring": "docstring1",
                "errors": [
                    ("ER01", "err desc"),
                    ("ER02", "err desc"),
                    ("ER03", "err desc"),
                ],
                "warnings": [],
                "examples_errors": "",
            },
        )
        exit_status = validate_docstrings.main(
            func_name="docstring1",
            prefix=None,
            errors=[],
            output_format="default",
            ignore_deprecated=False,
        )
        assert exit_status == 0

    def test_exit_status_errors_for_validate_all(self, monkeypatch):
        monkeypatch.setattr(
            validate_docstrings,
            "validate_all",
            lambda prefix, ignore_deprecated=False: {
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
            errors=[],
            output_format="default",
            ignore_deprecated=False,
        )
        assert exit_status == 5

    def test_no_exit_status_noerrors_for_validate_all(self, monkeypatch):
        monkeypatch.setattr(
            validate_docstrings,
            "validate_all",
            lambda prefix, ignore_deprecated=False: {
                "docstring1": {"errors": [], "warnings": [("WN01", "warn desc")]},
                "docstring2": {"errors": []},
            },
        )
        exit_status = validate_docstrings.main(
            func_name=None,
            prefix=None,
            errors=[],
            output_format="default",
            ignore_deprecated=False,
        )
        assert exit_status == 0

    def test_exit_status_for_validate_all_json(self, monkeypatch):
        print("EXECUTED")
        monkeypatch.setattr(
            validate_docstrings,
            "validate_all",
            lambda prefix, ignore_deprecated=False: {
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
            prefix=None,
            errors=[],
            output_format="json",
            ignore_deprecated=False,
        )
        assert exit_status == 0

    def test_errors_param_filters_errors(self, monkeypatch):
        monkeypatch.setattr(
            validate_docstrings,
            "validate_all",
            lambda prefix, ignore_deprecated=False: {
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
        exit_status = validate_docstrings.main(
            func_name=None,
            prefix=None,
            errors=["ER01"],
            output_format="default",
            ignore_deprecated=False,
        )
        assert exit_status == 3

        exit_status = validate_docstrings.main(
            func_name=None,
            prefix=None,
            errors=["ER03"],
            output_format="default",
            ignore_deprecated=False,
        )
        assert exit_status == 1
