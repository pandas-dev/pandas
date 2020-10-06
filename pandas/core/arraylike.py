"""
Methods that can be shared by many array-like classes or subclasses:
    Series
    Index
    ExtensionArray
"""
import operator

from pandas.errors import AbstractMethodError


class OpsMixin:
    # -------------------------------------------------------------
    # Comparisons

    def _cmp_method(self, other, op):
        raise AbstractMethodError(self)

    def __eq__(self, other):
        return self._cmp_method(other, operator.eq)

    def __ne__(self, other):
        return self._cmp_method(other, operator.ne)

    def __lt__(self, other):
        return self._cmp_method(other, operator.lt)

    def __le__(self, other):
        return self._cmp_method(other, operator.le)

    def __gt__(self, other):
        return self._cmp_method(other, operator.gt)

    def __ge__(self, other):
        return self._cmp_method(other, operator.ge)
