"""
Expose public exceptions & warnings
"""
from __future__ import annotations

import ctypes

from pandas._config.config import OptionError

from pandas._libs.tslibs import (
    OutOfBoundsDatetime,
    OutOfBoundsTimedelta,
)


class IntCastingNaNError(ValueError):
    """
    Raised when attempting an astype operation on an array with NaN to an integer
    dtype.
    """

    pass


class NullFrequencyError(ValueError):
    """
    Error raised when a null `freq` attribute is used in an operation
    that needs a non-null frequency, particularly `DatetimeIndex.shift`,
    `TimedeltaIndex.shift`, `PeriodIndex.shift`.
    """

    pass


class PerformanceWarning(Warning):
    """
    Warning raised when there is a possible performance impact.
    """


class UnsupportedFunctionCall(ValueError):
    """
    Exception raised when attempting to call a numpy function
    on a pandas object, but that function is not supported by
    the object e.g. ``np.cumsum(groupby_object)``.
    """


class UnsortedIndexError(KeyError):
    """
    Error raised when attempting to get a slice of a MultiIndex,
    and the index has not been lexsorted. Subclass of `KeyError`.
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

    >>> df = pd.DataFrame({'a': (['1'] * 100000 + ['X'] * 100000 +
    ...                          ['1'] * 100000),
    ...                    'b': ['b'] * 300000})  # doctest: +SKIP
    >>> df.to_csv('test.csv', index=False)  # doctest: +SKIP
    >>> df2 = pd.read_csv('test.csv')  # doctest: +SKIP
    ... # DtypeWarning: Columns (0) have mixed types

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

    >>> df2 = pd.read_csv('test.csv', sep=',', dtype={'a': str})  # doctest: +SKIP

    No warning was issued.
    """


class EmptyDataError(ValueError):
    """
    Exception that is thrown in `pd.read_csv` (by both the C and
    Python engines) when empty data or header is encountered.
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
    3. `sep=None` with `delim_whitespace=False`

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
    >>> df = pd.read_csv(io.StringIO(csv), sep='[;,]')  # doctest: +SKIP
    ... # ParserWarning: Falling back to the 'python' engine...

    Adding `engine='python'` to `pd.read_csv` removes the Warning:

    >>> df = pd.read_csv(io.StringIO(csv), sep='[;,]', engine='python')
    """


class MergeError(ValueError):
    """
    Error raised when problems arise during merging due to problems
    with input data. Subclass of `ValueError`.
    """


class AccessorRegistrationWarning(Warning):
    """
    Warning for attribute conflicts in accessor registration.
    """


class AbstractMethodError(NotImplementedError):
    """
    Raise this error instead of NotImplementedError for abstract methods
    while keeping compatibility with Python 2 and Python 3.
    """

    def __init__(self, class_instance, methodtype="method") -> None:
        types = {"method", "classmethod", "staticmethod", "property"}
        if methodtype not in types:
            raise ValueError(
                f"methodtype must be one of {methodtype}, got {types} instead."
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
    """


class DuplicateLabelError(ValueError):
    """
    Error raised when an operation would introduce duplicate labels.

    .. versionadded:: 1.2.0

    Examples
    --------
    >>> s = pd.Series([0, 1, 2], index=['a', 'b', 'c']).set_flags(
    ...     allows_duplicate_labels=False
    ... )
    >>> s.reindex(['a', 'a', 'b'])
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

    .. versionadded:: 1.1.0
    """


class DataError(Exception):
    """
    Exception raised when trying to perform a ohlc on a non-numnerical column.
    Or, it can be raised when trying to apply a function to a non-numerical
    column on a rolling window.
    """


class SpecificationError(Exception):
    """
    Exception raised in two scenarios. The first way is calling agg on a
    Dataframe or Series using a nested renamer (dict-of-dict).
    The second way is calling agg on a Dataframe with duplicated functions
    names without assigning column name.

    Examples
    --------
    >>> df = pd.DataFrame({'A': [1, 1, 1, 2, 2],
    ...                    'B': range(5),
    ...                    'C': range(5)})
    >>> df.groupby('A').B.agg({'foo': 'count'}) # doctest: +SKIP
    ... # SpecificationError: nested renamer is not supported

    >>> df.groupby('A').agg({'B': {'foo': ['sum', 'max']}}) # doctest: +SKIP
    ... # SpecificationError: nested renamer is not supported

    >>> df.groupby('A').agg(['min', 'min']) # doctest: +SKIP
    ... # SpecificationError: nested renamer is not supported
    """


class SettingWithCopyError(ValueError):
    """
    Exception is raised when trying to set on a copied slice from a dataframe and
    the mode.chained_assignment is set to 'raise.' This can happen unintentionally
    when chained indexing.

    For more information on eveluation order,
    see :ref:`the user guide<indexing.evaluation_order>`.

    For more information on view vs. copy,
    see :ref:`the user guide<indexing.view_versus_copy>`.

    Examples
    --------
    >>> pd.options.mode.chained_assignment = 'raise'
    >>> df = pd.DataFrame({'A': [1, 1, 1, 2, 2]}, columns=['A'])
    >>> df.loc[0:3]['A'] = 'a' # doctest: +SKIP
    ... # SettingWithCopyError: A value is trying to be set on a copy of a...
    """


class SettingWithCopyWarning(Warning):
    """
    Warning is raised when trying to set on a copied slice from a dataframe and
    the mode.chained_assignment is set to 'warn.' 'Warn' is the default option.
    This can happen unintentionally when chained indexing.

    For more information on eveluation order,
    see :ref:`the user guide<indexing.evaluation_order>`.

    For more information on view vs. copy,
    see :ref:`the user guide<indexing.view_versus_copy>`.

    Examples
    --------
    >>> df = pd.DataFrame({'A': [1, 1, 1, 2, 2]}, columns=['A'])
    >>> df.loc[0:3]['A'] = 'a' # doctest: +SKIP
    ... # SettingWithCopyWarning: A value is trying to be set on a copy of a...
    """


class NumExprClobberingError(NameError):
    """
    Exception is raised when trying to use a built-in numexpr name as a variable name
    in a method like query or eval. Eval will throw the error if the engine is set
    to 'numexpr'. 'numexpr' is the default engine value for eval if the numexpr package
    is installed.

    Examples
    --------
    >>> df = pd.DataFrame({'abs': [1, 1, 1]})
    >>> df.query("abs > 2") # doctest: +SKIP
    ... # NumExprClobberingError: Variables in expression "(abs) > (2)" overlap...
    >>> sin, a = 1, 2
    >>> pd.eval("sin + a", engine='numexpr') # doctest: +SKIP
    ... # NumExprClobberingError: Variables in expression "(sin) + (a)" overlap...
    """


class UndefinedVariableError(NameError):
    """
    Exception is raised when trying to use an undefined variable name in a method
    like query or eval. It will also specific whether the undefined variable is
    local or not.

    Examples
    --------
    >>> df = pd.DataFrame({'A': [1, 1, 1]})
    >>> df.query("A > x") # doctest: +SKIP
    ... # UndefinedVariableError: name 'x' is not defined
    >>> df.query("A > @y") # doctest: +SKIP
    ... # UndefinedVariableError: local variable 'y' is not defined
    >>> pd.eval('x + 1') # doctest: +SKIP
    ... # UndefinedVariableError: name 'x' is not defined
    """

    def __init__(self, name: str, is_local: bool | None = None) -> None:
        base_msg = f"{repr(name)} is not defined"
        if is_local:
            msg = f"local variable {base_msg}"
        else:
            msg = f"name {base_msg}"
        super().__init__(msg)


class IndexingError(Exception):
    """
    Exception is raised when trying to index and there is a mismatch in dimensions.

    Examples
    --------
    >>> df = pd.DataFrame({'A': [1, 1, 1]})
    >>> df.loc[..., ..., 'A'] # doctest: +SKIP
    ... # IndexingError: indexer may only contain one '...' entry
    >>> df = pd.DataFrame({'A': [1, 1, 1]})
    >>> df.loc[1, ..., ...] # doctest: +SKIP
    ... # IndexingError: Too many indexers
    >>> df[pd.Series([True], dtype=bool)] # doctest: +SKIP
    ... # IndexingError: Unalignable boolean Series provided as indexer...
    >>> s = pd.Series(range(2),
    ...               index = pd.MultiIndex.from_product([["a", "b"], ["c"]]))
    >>> s.loc["a", "c", "d"] # doctest: +SKIP
    ... # IndexingError: Too many indexers
    """


class PyperclipException(RuntimeError):
    """
    Exception is raised when trying to use methods like to_clipboard() and
    read_clipboard() on an unsupported OS/platform.
    """


class PyperclipWindowsException(PyperclipException):
    """
    Exception is raised when pandas is unable to get access to the clipboard handle
    due to some other window process is accessing it.
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

    Examples
    --------
    >>> df = pd.DataFrame({'A': [1, 1, 1]})
    >>> df.style.applymap(lambda x: 'background-color: blueGreenRed;')
    ...         .to_excel('styled.xlsx') # doctest: +SKIP
    ... # CSSWarning: Unhandled color format: 'blueGreenRed'
    >>> df.style.applymap(lambda x: 'border: 1px solid red red;')
    ...         .to_excel('styled.xlsx') # doctest: +SKIP
    ... # CSSWarning: Too many tokens provided to "border" (expected 1-3)
    """


class PossibleDataLossError(Exception):
    """
    Exception is raised when trying to open a HDFStore file when the file is already
    opened.

    Examples
    --------
    >>> store = pd.HDFStore('my-store', 'a') # doctest: +SKIP
    >>> store.open("w") # doctest: +SKIP
    ... # PossibleDataLossError: Re-opening the file [my-store] with mode [a]...
    """


class ClosedFileError(Exception):
    """
    Exception is raised when trying to perform an operation on a closed HDFStore file.

    Examples
    --------
    >>> store = pd.HDFStore('my-store', 'a') # doctest: +SKIP
    >>> store.close() # doctest: +SKIP
    >>> store.keys() # doctest: +SKIP
    ... # ClosedFileError: my-store file is not open!
    """


class IncompatibilityWarning(Warning):
    """
    Warning is raised when trying to use where criteria on an incompatible
    HDF5 file.
    """


class AttributeConflictWarning(Warning):
    """
    Warning is raised when attempting to append an index with a different
    name than the existing index on an HDFStore or attempting to append an index with a
    different frequency than the existing index on an HDFStore.
    """


class DatabaseError(OSError):
    """
    Error is raised when executing sql with bad syntax or sql that throws an error.

    Examples
    --------
    >>> from sqlite3 import connect
    >>> conn = connect(':memory:')
    >>> pd.read_sql('select * test', conn) # doctest: +SKIP
    ... # DatabaseError: Execution failed on sql 'test': near "test": syntax error
    """


__all__ = [
    "AbstractMethodError",
    "AccessorRegistrationWarning",
    "AttributeConflictWarning",
    "ClosedFileError",
    "CSSWarning",
    "DatabaseError",
    "DataError",
    "DtypeWarning",
    "DuplicateLabelError",
    "EmptyDataError",
    "IncompatibilityWarning",
    "IntCastingNaNError",
    "InvalidIndexError",
    "IndexingError",
    "MergeError",
    "NullFrequencyError",
    "NumbaUtilError",
    "NumExprClobberingError",
    "OptionError",
    "OutOfBoundsDatetime",
    "OutOfBoundsTimedelta",
    "ParserError",
    "ParserWarning",
    "PerformanceWarning",
    "PossibleDataLossError",
    "PyperclipException",
    "PyperclipWindowsException",
    "SettingWithCopyError",
    "SettingWithCopyWarning",
    "SpecificationError",
    "UndefinedVariableError",
    "UnsortedIndexError",
    "UnsupportedFunctionCall",
]
