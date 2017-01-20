
class PandasError(Exception):
    pass


class PerformanceWarning(Warning):
    pass


class SettingWithCopyError(ValueError):
    pass


class SettingWithCopyWarning(Warning):
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


class AbstractMethodError(NotImplementedError):
    """Raise this error instead of NotImplementedError for abstract methods
    while keeping compatibility with Python 2 and Python 3.
    """

    def __init__(self, class_instance):
        self.class_instance = class_instance

    def __str__(self):
        return ("This method must be defined in the concrete class of %s" %
                self.class_instance.__class__.__name__)
