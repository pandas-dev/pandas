"""
Expose public exceptions & warnings
"""

from __future__ import annotations

import abc
import ctypes

from pandas._config.config import OptionError

from pandas._libs.tslibs import (
    IncompatibleFrequency,
    OutOfBoundsDatetime,
    OutOfBoundsTimedelta,
)

from pandas.util.version import InvalidVersion


class IntCastingNaNError(ValueError):
    """
    Exception raised when converting (``astype``) an array with NaN to an integer type.

    This error occurs when attempting to cast a data structure containing non-finite
    values (such as NaN or infinity) to an integer data type. Integer types do not
    support non-finite values, so such conversions are explicitly disallowed to
    prevent silent data corruption or unexpected behavior.

    See Also
    --------
    DataFrame.astype : Method to cast a pandas DataFrame object to a specified dtype.
    Series.astype : Method to cast a pandas Series object to a specified dtype.

    Examples
    --------
    >>> pd.DataFrame(np.array([[1, np.nan], [2, 3]]), dtype="i8")
    Traceback (most recent call last):
    IntCastingNaNError: Cannot convert non-finite values (NA or inf) to integer
    """


class NullFrequencyError(ValueError):
    """
    Exception raised when a ``freq`` cannot be null.

    Particularly ``DatetimeIndex.shift``, ``TimedeltaIndex.shift``,
    ``PeriodIndex.shift``.

    See Also
    --------
    Index.shift : Shift values of Index.
    Series.shift : Shift values of Series.

    Examples
    --------
    >>> df = pd.DatetimeIndex(["2011-01-01 10:00", "2011-01-01"], freq=None)
    >>> df.shift(2)
    Traceback (most recent call last):
    NullFrequencyError: Cannot shift with no freq
    """


class PerformanceWarning(Warning):
    """
    Warning raised when there is a possible performance impact.

    See Also
    --------
    DataFrame.set_index : Set the DataFrame index using existing columns.
    DataFrame.loc : Access a group of rows and columns by label(s) \
    or a boolean array.

    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {"jim": [0, 0, 1, 1], "joe": ["x", "x", "z", "y"], "jolie": [1, 2, 3, 4]}
    ... )
    >>> df = df.set_index(["jim", "joe"])
    >>> df
              jolie
    jim  joe
    0    x    1
         x    2
    1    z    3
         y    4
    >>> df.loc[(1, "z")]  # doctest: +SKIP
    # PerformanceWarning: indexing past lexsort depth may impact performance.
    df.loc[(1, 'z')]
              jolie
    jim  joe
    1    z        3
    """


class PandasChangeWarning(Warning):
    """
    Warning raised for any upcoming change.

    See Also
    --------
    errors.PandasPendingDeprecationWarning : Class for deprecations that will raise a
        PendingDeprecationWarning.
    errors.PandasDeprecationWarning : Class for deprecations that will raise a
        DeprecationWarning.
    errors.PandasFutureWarning : Class for deprecations that will raise a FutureWarning.

    Examples
    --------
    >>> pd.errors.PandasChangeWarning
    <class 'pandas.errors.PandasChangeWarning'>
    """

    @classmethod
    @abc.abstractmethod
    def version(cls) -> str:
        """Version where change will be enforced."""


class PandasPendingDeprecationWarning(PandasChangeWarning, PendingDeprecationWarning):
    """
    Warning raised for an upcoming change that is a PendingDeprecationWarning.

    See Also
    --------
    errors.PandasChangeWarning: Class for deprecations that will raise any warning.
    errors.PandasDeprecationWarning : Class for deprecations that will raise a
        DeprecationWarning.
    errors.PandasFutureWarning : Class for deprecations that will raise a FutureWarning.

    Examples
    --------
    >>> pd.errors.PandasPendingDeprecationWarning
    <class 'pandas.errors.PandasPendingDeprecationWarning'>
    """


class PandasDeprecationWarning(PandasChangeWarning, DeprecationWarning):
    """
    Warning raised for an upcoming change that is a DeprecationWarning.

    See Also
    --------
    errors.PandasChangeWarning: Class for deprecations that will raise any warning.
    errors.PandasPendingDeprecationWarning : Class for deprecations that will raise a
        PendingDeprecationWarning.
    errors.PandasFutureWarning : Class for deprecations that will raise a FutureWarning.

    Examples
    --------
    >>> pd.errors.PandasDeprecationWarning
    <class 'pandas.errors.PandasDeprecationWarning'>
    """


class PandasFutureWarning(PandasChangeWarning, FutureWarning):
    """
    Warning raised for an upcoming change that is a FutureWarning.

    See Also
    --------
    errors.PandasChangeWarning: Class for deprecations that will raise any warning.
    errors.PandasPendingDeprecationWarning : Class for deprecations that will raise a
        PendingDeprecationWarning.
    errors.PandasDeprecationWarning : Class for deprecations that will raise a
        DeprecationWarning.

    Examples
    --------
    >>> pd.errors.PandasFutureWarning
    <class 'pandas.errors.PandasFutureWarning'>
    """


class Pandas4Warning(PandasDeprecationWarning):
    """
    Warning raised for an upcoming change that will be enforced in pandas 4.0.

    See Also
    --------
    errors.PandasChangeWarning: Class for deprecations that will raise any warning.
    errors.PandasPendingDeprecationWarning : Class for deprecations that will raise a
        PendingDeprecationWarning.
    errors.PandasDeprecationWarning : Class for deprecations that will raise a
        DeprecationWarning.
    errors.PandasFutureWarning : Class for deprecations that will raise a FutureWarning.

    Examples
    --------
    >>> pd.errors.Pandas4Warning
    <class 'pandas.errors.Pandas4Warning'>
    """

    @classmethod
    def version(cls) -> str:
        """Version where change will be enforced."""
        return "4.0"


class Pandas5Warning(PandasPendingDeprecationWarning):
    """
    Warning raised for an upcoming change that will be enforced in pandas 5.0.

    See Also
    --------
    errors.PandasChangeWarning: Class for deprecations that will raise any warning.
    errors.PandasPendingDeprecationWarning : Class for deprecations that will raise a
        PendingDeprecationWarning.
    errors.PandasDeprecationWarning : Class for deprecations that will raise a
        DeprecationWarning.
    errors.PandasFutureWarning : Class for deprecations that will raise a FutureWarning.

    Examples
    --------
    >>> pd.errors.Pandas5Warning
    <class 'pandas.errors.Pandas5Warning'>
    """

    @classmethod
    def version(cls) -> str:
        """Version where change will be enforced."""
        return "5.0"


_CurrentDeprecationWarning = Pandas4Warning


class UnsupportedFunctionCall(ValueError):
    """
    Exception raised when attempting to call a unsupported numpy function.

    For example, ``np.cumsum(groupby_object)``.

    See Also
    --------
    DataFrame.groupby : Group DataFrame using a mapper or by a Series of columns.
    Series.groupby : Group Series using a mapper or by a Series of columns.
    core.groupby.GroupBy.cumsum : Compute cumulative sum for each group.

    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {"A": [0, 0, 1, 1], "B": ["x", "x", "z", "y"], "C": [1, 2, 3, 4]}
    ... )
    >>> np.cumsum(df.groupby(["A"]))
    Traceback (most recent call last):
    UnsupportedFunctionCall: numpy operations are not valid with groupby.
    Use .groupby(...).cumsum() instead
    """


class UnsortedIndexError(KeyError):
    """
    Error raised when slicing a MultiIndex which has not been lexsorted.

    Subclass of `KeyError`.

    See Also
    --------
    DataFrame.sort_index : Sort a DataFrame by its index.
    DataFrame.set_index : Set the DataFrame index using existing columns.

    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {
    ...         "cat": [0, 0, 1, 1],
    ...         "color": ["white", "white", "brown", "black"],
    ...         "lives": [4, 4, 3, 7],
    ...     },
    ... )
    >>> df = df.set_index(["cat", "color"])
    >>> df
                lives
    cat  color
    0    white    4
         white    4
    1    brown    3
         black    7
    >>> df.loc[(0, "black") : (1, "white")]
    Traceback (most recent call last):
    UnsortedIndexError: 'Key length (2) was greater
    than MultiIndex lexsort depth (1)'
    """


class ParserError(ValueError):
    """
    Exception that is raised by an error encountered in parsing file contents.

    This is a generic error raised for errors encountered when functions like
    `read_csv` or `read_html` are parsing contents of a file.

    See Also
    --------
    read_csv : Read CSV (comma-separated) file into a DataFrame.
    read_html : Read HTML table into a DataFrame.

    Examples
    --------
    >>> data = '''a,b,c
    ... cat,foo,bar
    ... dog,foo,"baz'''
    >>> from io import StringIO
    >>> pd.read_csv(StringIO(data), skipfooter=1, engine="python")
    Traceback (most recent call last):
    ParserError: ',' expected after '"'. Error could possibly be due
    to parsing errors in the skipped footer rows
    """


class DtypeWarning(Warning):
    """
    Warning raised when reading different dtypes in a column from a file.

    Raised for a dtype incompatibility. This can happen whenever `read_csv`
    or `read_table` encounter non-uniform dtypes in a column(s) of a given
    CSV file.

    See Also
    --------
    read_csv : Read CSV (comma-separated) file into a DataFrame.
    read_table : Read general delimited file into a DataFrame.

    Notes
    -----
    This warning is issued when dealing with larger files because the dtype
    checking happens per chunk read.

    Despite the warning, the CSV file is read with mixed types in a single
    column which will be an object type. See the examples below to better
    understand this issue.

    Examples
    --------
    This example creates and reads a large CSV file with a column that contains
    `int` and `str`.

    >>> df = pd.DataFrame(
    ...     {
    ...         "a": (["1"] * 100000 + ["X"] * 100000 + ["1"] * 100000),
    ...         "b": ["b"] * 300000,
    ...     }
    ... )  # doctest: +SKIP
    >>> df.to_csv("test.csv", index=False)  # doctest: +SKIP
    >>> df2 = pd.read_csv("test.csv")  # doctest: +SKIP
    ... # DtypeWarning: Columns (0: a) have mixed types

    Important to notice that ``df2`` will contain both `str` and `int` for the
    same input, '1'.

    >>> df2.iloc[262140, 0]  # doctest: +SKIP
    '1'
    >>> type(df2.iloc[262140, 0])  # doctest: +SKIP
    <class 'str'>
    >>> df2.iloc[262150, 0]  # doctest: +SKIP
    1
    >>> type(df2.iloc[262150, 0])  # doctest: +SKIP
    <class 'int'>

    One way to solve this issue is using the `dtype` parameter in the
    `read_csv` and `read_table` functions to explicit the conversion:

    >>> df2 = pd.read_csv("test.csv", sep=",", dtype={"a": str})  # doctest: +SKIP

    No warning was issued.
    """


class EmptyDataError(ValueError):
    """
    Exception raised in ``pd.read_csv`` when empty data or header is encountered.

    This error is typically encountered when attempting to read an empty file or
    an invalid file where no data or headers are present.

    See Also
    --------
    read_csv : Read a comma-separated values (CSV) file into DataFrame.
    errors.ParserError : Exception that is raised by an error encountered in parsing
        file contents.
    errors.DtypeWarning : Warning raised when reading different dtypes in a column
        from a file.

    Examples
    --------
    >>> from io import StringIO
    >>> empty = StringIO()
    >>> pd.read_csv(empty)
    Traceback (most recent call last):
    EmptyDataError: No columns to parse from file
    """


class ParserWarning(Warning):
    """
    Warning raised when reading a file that doesn't use the default 'c' parser.

    Raised by `pd.read_csv` and `pd.read_table` when it is necessary to change
    parsers, generally from the default 'c' parser to 'python'.

    It happens due to a lack of support or functionality for parsing a
    particular attribute of a CSV file with the requested engine.

    Currently, 'c' unsupported options include the following parameters:

    1. `sep` other than a single character (e.g. regex separators)
    2. `skipfooter` higher than 0

    The warning can be avoided by adding `engine='python'` as a parameter in
    `pd.read_csv` and `pd.read_table` methods.

    See Also
    --------
    pd.read_csv : Read CSV (comma-separated) file into DataFrame.
    pd.read_table : Read general delimited file into DataFrame.

    Examples
    --------
    Using a `sep` in `pd.read_csv` other than a single character:

    >>> import io
    >>> csv = '''a;b;c
    ...           1;1,8
    ...           1;2,1'''
    >>> df = pd.read_csv(io.StringIO(csv), sep="[;,]")  # doctest: +SKIP
    ... # ParserWarning: Falling back to the 'python' engine...

    Adding `engine='python'` to `pd.read_csv` removes the Warning:

    >>> df = pd.read_csv(io.StringIO(csv), sep="[;,]", engine="python")
    """


class MergeError(ValueError):
    """
    Exception raised when merging data.

    Subclass of ``ValueError``.

    See Also
    --------
    DataFrame.join : For joining DataFrames on their indexes.
    merge : For merging two DataFrames on a common set of keys.

    Examples
    --------
    >>> left = pd.DataFrame(
    ...     {"a": ["a", "b", "b", "d"], "b": ["cat", "dog", "weasel", "horse"]},
    ...     index=range(4),
    ... )
    >>> right = pd.DataFrame(
    ...     {"a": ["a", "b", "c", "d"], "c": ["meow", "bark", "chirp", "nay"]},
    ...     index=range(4),
    ... ).set_index("a")
    >>> left.join(
    ...     right,
    ...     on="a",
    ...     validate="one_to_one",
    ... )
    Traceback (most recent call last):
    MergeError: Merge keys are not unique in left dataset; not a one-to-one merge
    """


class AbstractMethodError(NotImplementedError):
    """
    Raise this error instead of NotImplementedError for abstract methods.

    The `AbstractMethodError` is designed for use in classes that follow an abstract
    base class pattern. By raising this error in the method, it ensures that a subclass
    must implement the method to provide specific functionality. This is useful in a
    framework or library where certain methods must be implemented by the user to
    ensure correct behavior.

    Parameters
    ----------
    class_instance : object
        The instance of the class where the abstract method is being called.
    methodtype : str, default "method"
        A string indicating the type of method that is abstract.
        Must be one of {"method", "classmethod", "staticmethod", "property"}.

    See Also
    --------
    api.extensions.ExtensionArray
        An example of a pandas extension mechanism that requires implementing
        specific abstract methods.
    NotImplementedError
        A built-in exception that can also be used for abstract methods but lacks
        the specificity of `AbstractMethodError` in indicating the need for subclass
        implementation.

    Examples
    --------
    >>> class Foo:
    ...     @classmethod
    ...     def classmethod(cls):
    ...         raise pd.errors.AbstractMethodError(cls, methodtype="classmethod")
    ...
    ...     def method(self):
    ...         raise pd.errors.AbstractMethodError(self)
    >>> test = Foo.classmethod()
    Traceback (most recent call last):
    AbstractMethodError: This classmethod must be defined in the concrete class Foo

    >>> test2 = Foo().method()
    Traceback (most recent call last):
    AbstractMethodError: This classmethod must be defined in the concrete class Foo
    """

    def __init__(self, class_instance, methodtype: str = "method") -> None:
        types = {"method", "classmethod", "staticmethod", "property"}
        if methodtype not in types:
            raise ValueError(
                f"methodtype must be one of {types}, got {methodtype} instead."
            )
        self.methodtype = methodtype
        self.class_instance = class_instance

    def __str__(self) -> str:
        if self.methodtype == "classmethod":
            name = self.class_instance.__name__
        else:
            name = type(self.class_instance).__name__
        return f"This {self.methodtype} must be defined in the concrete class {name}"


class NumbaUtilError(Exception):
    """
    Error raised for unsupported Numba engine routines.

    See Also
    --------
    DataFrame.groupby : Group DataFrame using a mapper or by a Series of columns.
    Series.groupby : Group Series using a mapper or by a Series of columns.
    DataFrame.agg : Aggregate using one or more operations over the specified axis.
    Series.agg : Aggregate using one or more operations over the specified axis.

    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {"key": ["a", "a", "b", "b"], "data": [1, 2, 3, 4]}, columns=["key", "data"]
    ... )
    >>> def incorrect_function(x):
    ...     return sum(x) * 2.7
    >>> df.groupby("key").agg(incorrect_function, engine="numba")
    Traceback (most recent call last):
    NumbaUtilError: The first 2 arguments to incorrect_function
    must be ['values', 'index']
    """


class DuplicateLabelError(ValueError):
    """
    Error raised when an operation would introduce duplicate labels.

    This error is typically encountered when performing operations on objects
    with `allows_duplicate_labels=False` and the operation would result in
    duplicate labels in the index. Duplicate labels can lead to ambiguities
    in indexing and reduce data integrity.

    See Also
    --------
    Series.set_flags : Return a new ``Series`` object with updated flags.
    DataFrame.set_flags : Return a new ``DataFrame`` object with updated flags.
    Series.reindex : Conform ``Series`` object to new index with optional filling logic.
    DataFrame.reindex : Conform ``DataFrame`` object to new index with optional filling
        logic.

    Examples
    --------
    >>> s = pd.Series([0, 1, 2], index=["a", "b", "c"]).set_flags(
    ...     allows_duplicate_labels=False
    ... )
    >>> s.reindex(["a", "a", "b"])
    Traceback (most recent call last):
       ...
    DuplicateLabelError: Index has duplicates.
          positions
    label
    a        [0, 1]
    """


class InvalidIndexError(Exception):
    """
    Exception raised when attempting to use an invalid index key.

    This exception is triggered when a user attempts to access or manipulate
    data in a pandas DataFrame or Series using an index key that is not valid
    for the given object. This may occur in cases such as using a malformed
    slice, a mismatched key for a ``MultiIndex``, or attempting to access an index
    element that does not exist.

    See Also
    --------
    MultiIndex : A multi-level, or hierarchical, index object for pandas objects.

    Examples
    --------
    >>> idx = pd.MultiIndex.from_product([["x", "y"], [0, 1]])
    >>> df = pd.DataFrame([[1, 1, 2, 2], [3, 3, 4, 4]], columns=idx)
    >>> df
        x       y
        0   1   0   1
    0   1   1   2   2
    1   3   3   4   4
    >>> df[:, 0]
    Traceback (most recent call last):
    InvalidIndexError: (slice(None, None, None), 0)
    """


class DataError(Exception):
    """
    Exception raised when performing an operation on non-numerical data.

    For example, calling ``ohlc`` on a non-numerical column or a function
    on a rolling window.

    See Also
    --------
    Series.rolling : Provide rolling window calculations on Series object.
    DataFrame.rolling : Provide rolling window calculations on DataFrame object.

    Examples
    --------
    >>> ser = pd.Series(["a", "b", "c"])
    >>> ser.rolling(2).sum()
    Traceback (most recent call last):
    DataError: No numeric types to aggregate
    """


class SpecificationError(Exception):
    """
    Exception raised by ``agg`` when the functions are ill-specified.

    The exception raised in two scenarios.

    The first way is calling ``agg`` on a
    Dataframe or Series using a nested renamer (dict-of-dict).

    The second way is calling ``agg`` on a Dataframe with duplicated functions
    names without assigning column name.

    See Also
    --------
    DataFrame.agg : Aggregate using one or more operations over the specified axis.
    Series.agg : Aggregate using one or more operations over the specified axis.

    Examples
    --------
    >>> df = pd.DataFrame({"A": [1, 1, 1, 2, 2], "B": range(5), "C": range(5)})
    >>> df.groupby("A").B.agg({"foo": "count"})  # doctest: +SKIP
    ... # SpecificationError: nested renamer is not supported

    >>> df.groupby("A").agg({"B": {"foo": ["sum", "max"]}})  # doctest: +SKIP
    ... # SpecificationError: nested renamer is not supported

    >>> df.groupby("A").agg(["min", "min"])  # doctest: +SKIP
    ... # SpecificationError: nested renamer is not supported
    """


class ChainedAssignmentError(Warning):
    """
    Warning raised when trying to set using chained assignment.

    With Copy-on-Write now always enabled, chained assignment can
    never work. In such a situation, we are always setting into a temporary
    object that is the result of an indexing operation (getitem), which under
    Copy-on-Write always behaves as a copy. Thus, assigning through a chain
    can never update the original Series or DataFrame.

    For more information on Copy-on-Write,
    see :ref:`the user guide<copy_on_write>`.

    See Also
    --------
    DataFrame.loc : Access a group of rows and columns by label(s) or a boolean array.
    DataFrame.iloc : Purely integer-location based indexing for selection by position.
    Series.loc : Access a group of rows by label(s) or a boolean array.

    Examples
    --------
    >>> df = pd.DataFrame({"A": [1, 1, 1, 2, 2]}, columns=["A"])
    >>> df["A"][0:3] = 10  # doctest: +SKIP
    ... # ChainedAssignmentError: ...
    """


class NumExprClobberingError(NameError):
    """
    Exception raised when trying to use a built-in numexpr name as a variable name.

    ``eval`` or ``query`` will throw the error if the engine is set
    to 'numexpr'. 'numexpr' is the default engine value for these methods if the
    numexpr package is installed.

    See Also
    --------
    eval : Evaluate a Python expression as a string using various backends.
    DataFrame.query : Query the columns of a DataFrame with a boolean expression.

    Examples
    --------
    >>> df = pd.DataFrame({"abs": [1, 1, 1]})
    >>> df.query("abs > 2")  # doctest: +SKIP
    ... # NumExprClobberingError: Variables in expression "(abs) > (2)" overlap...
    >>> sin, a = 1, 2
    >>> pd.eval("sin + a", engine="numexpr")  # doctest: +SKIP
    ... # NumExprClobberingError: Variables in expression "(sin) + (a)" overlap...
    """


class UndefinedVariableError(NameError):
    """
    Exception raised by ``query`` or ``eval`` when using an undefined variable name.

    It will also specify whether the undefined variable is local or not.

    Parameters
    ----------
    name : str
        The name of the undefined variable.
    is_local : bool or None, optional
        Indicates whether the undefined variable is considered a local variable.
        If ``True``, the error message specifies it as a local variable.
        If ``False`` or ``None``, the variable is treated as a non-local name.

    See Also
    --------
    DataFrame.query : Query the columns of a DataFrame with a boolean expression.
    DataFrame.eval : Evaluate a string describing operations on DataFrame columns.

    Examples
    --------
    >>> df = pd.DataFrame({"A": [1, 1, 1]})
    >>> df.query("A > x")  # doctest: +SKIP
    ... # UndefinedVariableError: name 'x' is not defined
    >>> df.query("A > @y")  # doctest: +SKIP
    ... # UndefinedVariableError: local variable 'y' is not defined
    >>> pd.eval("x + 1")  # doctest: +SKIP
    ... # UndefinedVariableError: name 'x' is not defined
    """

    def __init__(self, name: str, is_local: bool | None = None) -> None:
        base_msg = f"{name!r} is not defined"
        if is_local:
            msg = f"local variable {base_msg}"
        else:
            msg = f"name {base_msg}"
        super().__init__(msg)


class IndexingError(Exception):
    """
    Exception is raised when trying to index and there is a mismatch in dimensions.

    Raised by properties like :attr:`.pandas.DataFrame.iloc` when
    an indexer is out of bounds or :attr:`.pandas.DataFrame.loc` when its index is
    unalignable to the frame index.

    See Also
    --------
    DataFrame.iloc : Purely integer-location based indexing for \
    selection by position.
    DataFrame.loc : Access a group of rows and columns by label(s) \
    or a boolean array.

    Examples
    --------
    >>> df = pd.DataFrame({"A": [1, 1, 1]})
    >>> df.loc[..., ..., "A"]  # doctest: +SKIP
    ... # IndexingError: indexer may only contain one '...' entry
    >>> df = pd.DataFrame({"A": [1, 1, 1]})
    >>> df.loc[1, ..., ...]  # doctest: +SKIP
    ... # IndexingError: Too many indexers
    >>> df[pd.Series([True], dtype=bool)]  # doctest: +SKIP
    ... # IndexingError: Unalignable boolean Series provided as indexer...
    >>> s = pd.Series(range(2), index=pd.MultiIndex.from_product([["a", "b"], ["c"]]))
    >>> s.loc["a", "c", "d"]  # doctest: +SKIP
    ... # IndexingError: Too many indexers
    """


class PyperclipException(RuntimeError):
    """
    Exception raised when clipboard functionality is unsupported.

    Raised by ``to_clipboard()`` and ``read_clipboard()``.
    """


class PyperclipWindowsException(PyperclipException):
    """
    Exception raised when clipboard functionality is unsupported by Windows.

    Access to the clipboard handle would be denied due to some other
    window process is accessing it.
    """

    def __init__(self, message: str) -> None:
        # attr only exists on Windows, so typing fails on other platforms
        message += f" ({ctypes.WinError()})"  # type: ignore[attr-defined]
        super().__init__(message)


class CSSWarning(UserWarning):
    """
    Warning is raised when converting css styling fails.

    This can be due to the styling not having an equivalent value or because the
    styling isn't properly formatted.

    See Also
    --------
    DataFrame.style : Returns a Styler object for applying CSS-like styles.
    io.formats.style.Styler : Helps style a DataFrame or Series according to the
        data with HTML and CSS.
    io.formats.style.Styler.to_excel : Export styled DataFrame to Excel.
    io.formats.style.Styler.to_html : Export styled DataFrame to HTML.

    Examples
    --------
    >>> df = pd.DataFrame({"A": [1, 1, 1]})
    >>> df.style.map(lambda x: "background-color: blueGreenRed;").to_excel(
    ...     "styled.xlsx"
    ... )  # doctest: +SKIP
    CSSWarning: Unhandled color format: 'blueGreenRed'
    >>> df.style.map(lambda x: "border: 1px solid red red;").to_excel(
    ...     "styled.xlsx"
    ... )  # doctest: +SKIP
    CSSWarning: Unhandled color format: 'blueGreenRed'
    """


class PossibleDataLossError(Exception):
    """
    Exception raised when trying to open an HDFStore file when already opened.

    This error is triggered when there is a potential risk of data loss due to
    conflicting operations on an HDFStore file. It serves to prevent unintended
    overwrites or data corruption by enforcing exclusive access to the file.

    See Also
    --------
    HDFStore : Dict-like IO interface for storing pandas objects in PyTables.
    HDFStore.open : Open an HDFStore file in the specified mode.

    Examples
    --------
    >>> store = pd.HDFStore("my-store", "a")  # doctest: +SKIP
    >>> store.open("w")  # doctest: +SKIP
    """


class ClosedFileError(Exception):
    """
    Exception is raised when trying to perform an operation on a closed HDFStore file.

    ``ClosedFileError`` is specific to operations on ``HDFStore`` objects. Once an
    HDFStore is closed, its resources are no longer available, and any further attempt
    to access data or perform file operations will raise this exception.

    See Also
    --------
    HDFStore.close : Closes the PyTables file handle.
    HDFStore.open : Opens the file in the specified mode.
    HDFStore.is_open : Returns a boolean indicating whether the file is open.

    Examples
    --------
    >>> store = pd.HDFStore("my-store", "a")  # doctest: +SKIP
    >>> store.close()  # doctest: +SKIP
    >>> store.keys()  # doctest: +SKIP
    ... # ClosedFileError: my-store file is not open!
    """


class IncompatibilityWarning(Warning):
    """
    Warning raised when trying to use where criteria on an incompatible HDF5 file.
    """


class AttributeConflictWarning(Warning):
    """
    Warning raised when index attributes conflict when using HDFStore.

    Occurs when attempting to append an index with a different
    name than the existing index on an HDFStore or attempting to append an index with a
    different frequency than the existing index on an HDFStore.

    See Also
    --------
    HDFStore : Dict-like IO interface for storing pandas objects in PyTables.
    DataFrame.to_hdf : Write the contained data to an HDF5 file using HDFStore.
    read_hdf : Read from an HDF5 file into a DataFrame.

    Examples
    --------
    >>> idx1 = pd.Index(["a", "b"], name="name1")
    >>> df1 = pd.DataFrame([[1, 2], [3, 4]], index=idx1)
    >>> df1.to_hdf("file", "data", "w", append=True)  # doctest: +SKIP
    >>> idx2 = pd.Index(["c", "d"], name="name2")
    >>> df2 = pd.DataFrame([[5, 6], [7, 8]], index=idx2)
    >>> df2.to_hdf("file", "data", "a", append=True)  # doctest: +SKIP
    AttributeConflictWarning: the [index_name] attribute of the existing index is
    [name1] which conflicts with the new [name2]...
    """


class DatabaseError(OSError):
    """
    Error is raised when executing SQL with bad syntax or SQL that throws an error.

    Raised by :func:`.pandas.read_sql` when a bad SQL statement is passed in.

    See Also
    --------
    read_sql : Read SQL query or database table into a DataFrame.

    Examples
    --------
    >>> from sqlite3 import connect
    >>> conn = connect(":memory:")
    >>> pd.read_sql("select * test", conn)  # doctest: +SKIP
    """


class PossiblePrecisionLoss(Warning):
    """
    Warning raised by to_stata on a column with a value outside or equal to int64.

    When the column value is outside or equal to the int64 value the column is
    converted to a float64 dtype.

    See Also
    --------
    DataFrame.to_stata : Export DataFrame object to Stata dta format.

    Examples
    --------
    >>> df = pd.DataFrame({"s": pd.Series([1, 2**53], dtype=np.int64)})
    >>> df.to_stata("test")  # doctest: +SKIP
    """


class ValueLabelTypeMismatch(Warning):
    """
    Warning raised by to_stata on a category column that contains non-string values.

    When exporting data to Stata format using the `to_stata` method, category columns
    must have string values as labels. If a category column contains non-string values
    (e.g., integers, floats, or other types), this warning is raised to indicate that
    the Stata file may not correctly represent the data.

    See Also
    --------
    DataFrame.to_stata : Export DataFrame object to Stata dta format.
    Series.cat : Accessor for categorical properties of the Series values.

    Examples
    --------
    >>> df = pd.DataFrame({"categories": pd.Series(["a", 2], dtype="category")})
    >>> df.to_stata("test")  # doctest: +SKIP
    """


class InvalidColumnName(Warning):
    """
    Warning raised by to_stata the column contains a non-valid stata name.

    Because the column name is an invalid Stata variable, the name needs to be
    converted.

    See Also
    --------
    DataFrame.to_stata : Export DataFrame object to Stata dta format.

    Examples
    --------
    >>> df = pd.DataFrame({"0categories": pd.Series([2, 2])})
    >>> df.to_stata("test")  # doctest: +SKIP
    """


class CategoricalConversionWarning(Warning):
    """
    Warning is raised when reading a partial labeled Stata file using an iterator.

    This warning helps ensure data integrity and alerts users to potential issues
    during the incremental reading of Stata files with labeled data, allowing for
    additional checks and adjustments as necessary.

    See Also
    --------
    read_stata : Read a Stata file into a DataFrame.
    Categorical : Represents a categorical variable in pandas.

    Examples
    --------
    >>> from pandas.io.stata import StataReader
    >>> with StataReader("dta_file", chunksize=2) as reader:  # doctest: +SKIP
    ...     for i, block in enumerate(reader):
    ...         print(i, block)
    ... # CategoricalConversionWarning: One or more series with value labels...
    """


class LossySetitemError(Exception):
    """
    Raised when trying to do a __setitem__ on an np.ndarray that is not lossless.

    Notes
    -----
    This is an internal error.
    """


class NoBufferPresent(Exception):
    """
    Exception is raised in _get_data_buffer to signal that there is no requested buffer.
    """


class InvalidComparison(Exception):
    """
    Exception is raised by _validate_comparison_value to indicate an invalid comparison.

    Notes
    -----
    This is an internal error.
    """


__all__ = [
    "AbstractMethodError",
    "AttributeConflictWarning",
    "CSSWarning",
    "CategoricalConversionWarning",
    "ChainedAssignmentError",
    "ClosedFileError",
    "DataError",
    "DatabaseError",
    "DtypeWarning",
    "DuplicateLabelError",
    "EmptyDataError",
    "IncompatibilityWarning",
    "IncompatibleFrequency",
    "IndexingError",
    "IntCastingNaNError",
    "InvalidColumnName",
    "InvalidComparison",
    "InvalidIndexError",
    "InvalidVersion",
    "LossySetitemError",
    "MergeError",
    "NoBufferPresent",
    "NullFrequencyError",
    "NumExprClobberingError",
    "NumbaUtilError",
    "OptionError",
    "OutOfBoundsDatetime",
    "OutOfBoundsTimedelta",
    "Pandas4Warning",
    "Pandas5Warning",
    "PandasChangeWarning",
    "PandasDeprecationWarning",
    "PandasFutureWarning",
    "PandasPendingDeprecationWarning",
    "ParserError",
    "ParserWarning",
    "PerformanceWarning",
    "PossibleDataLossError",
    "PossiblePrecisionLoss",
    "PyperclipException",
    "PyperclipWindowsException",
    "SpecificationError",
    "UndefinedVariableError",
    "UnsortedIndexError",
    "UnsupportedFunctionCall",
    "ValueLabelTypeMismatch",
]
