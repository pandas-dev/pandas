# flake8: noqa

""" expose public exceptions & warnings """

from pandas._libs.tslib import OutOfBoundsDatetime


class PandasError(Exception):
    pass


class PerformanceWarning(Warning):
    pass

class AmbiguousIndexError(PandasError, KeyError):
    pass


class UnsupportedFunctionCall(ValueError):
    pass


class UnsortedIndexError(KeyError):
    """ Error raised when attempting to get a slice of a MultiIndex
    and the index has not been lexsorted. Subclass of `KeyError`.

    .. versionadded:: 0.20.0

    """
    pass


class ParserError(ValueError):
    """
    Exception that is thrown by an error is encountered in `pd.read_csv`
    """
    pass


class DtypeWarning(Warning):
    """
    Warning that is raised whenever `pd.read_csv` encounters non-
    uniform dtypes in a column(s) of a given CSV file
    """
    pass


class EmptyDataError(ValueError):
    """
    Exception that is thrown in `pd.read_csv` (by both the C and
    Python engines) when empty data or header is encountered
    """
    pass


class ParserWarning(Warning):
    """
    Warning that is raised in `pd.read_csv` whenever it is necessary
    to change parsers (generally from 'c' to 'python') contrary to the
    one specified by the user due to lack of support or functionality for
    parsing particular attributes of a CSV file with the requsted engine
    """
    pass
