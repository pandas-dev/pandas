import abc
from typing import Pattern, Union

import numpy as np

from pandas._typing import Scalar

from pandas.core.arrays.base import ExtensionArray


class BaseStringArrayMethods(abc.ABC):
    """
    Base class for array _str accessor.

    This is where ExtensionArrays can override the implementation of
    Series.str.<method>. The rough layout is

    * User calls Series.str.<method>
    * pandas extracts the extension array from the Series
    * pandas calls ``extension_array._str.<method>(*args, **kwargs)``
    * pandas wraps the result, to return to the user.

    See :ref:`Series.str` for the docstring of each method.
    """

    def __init__(self, array: ExtensionArray):
        self._array = array

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.slice(start=key.start, stop=key.stop, step=key.step)
        else:
            return self.get(key)

    @abc.abstractmethod
    def count(self, pat, flags=0):
        pass

    @abc.abstractmethod
    def pad(self, width, side="left", fillchar=" "):
        pass

    @abc.abstractmethod
    def contains(self, pat, case=True, flags=0, na=np.nan, regex=True):
        pass

    @abc.abstractmethod
    def startswith(self, pat, na=np.nan):
        pass

    @abc.abstractmethod
    def endswith(self, pat, na=np.nan):
        pass

    @abc.abstractmethod
    def replace(self, pat, repl, n=-1, case=None, flags=0, regex=True):
        pass

    @abc.abstractmethod
    def repeat(self, repeats):
        pass

    @abc.abstractmethod
    def match(
        self,
        pat: Union[str, Pattern],
        case: bool = True,
        flags: int = 0,
        na: Scalar = np.nan,
    ):
        pass

    @abc.abstractmethod
    def fullmatch(
        self,
        pat: Union[str, Pattern],
        case: bool = True,
        flags: int = 0,
        na: Scalar = np.nan,
    ):
        pass

    @abc.abstractmethod
    def encode(self, encoding, errors="strict"):
        pass

    @abc.abstractmethod
    def find(self, sub, start=0, end=None):
        pass

    @abc.abstractmethod
    def rfind(self, sub, start=0, end=None):
        pass

    @abc.abstractmethod
    def findall(self, pat, flags=0):
        pass

    @abc.abstractmethod
    def get(self, i):
        pass

    @abc.abstractmethod
    def index(self, sub, start=0, end=None):
        pass

    @abc.abstractmethod
    def rindex(self, sub, start=0, end=None):
        pass

    @abc.abstractmethod
    def join(self, sep):
        pass

    @abc.abstractmethod
    def partition(self, sep, expand):
        pass

    @abc.abstractmethod
    def rpartition(self, sep, expand):
        pass

    @abc.abstractmethod
    def len(self):
        pass

    @abc.abstractmethod
    def slice(self, start=None, stop=None, step=None):
        pass

    @abc.abstractmethod
    def slice_replace(self, start=None, stop=None, repl=None):
        pass

    @abc.abstractmethod
    def translate(self, table):
        pass

    @abc.abstractmethod
    def wrap(self, width, **kwargs):
        pass

    @abc.abstractmethod
    def get_dummies(self, sep="|"):
        pass

    @abc.abstractmethod
    def isalnum(self):
        pass

    @abc.abstractmethod
    def isalpha(self):
        pass

    @abc.abstractmethod
    def isdecimal(self):
        pass

    @abc.abstractmethod
    def isdigit(self):
        pass

    @abc.abstractmethod
    def islower(self):
        pass

    @abc.abstractmethod
    def isnumeric(self):
        pass

    @abc.abstractmethod
    def isspace(self):
        pass

    @abc.abstractmethod
    def istitle(self):
        pass

    @abc.abstractmethod
    def isupper(self):
        pass

    @abc.abstractmethod
    def capitalize(self):
        pass

    @abc.abstractmethod
    def casefold(self):
        pass

    @abc.abstractmethod
    def title(self):
        pass

    @abc.abstractmethod
    def swapcase(self):
        pass

    @abc.abstractmethod
    def lower(self):
        pass

    @abc.abstractmethod
    def upper(self):
        pass

    @abc.abstractmethod
    def normalize(self, form):
        pass

    @abc.abstractmethod
    def strip(self, to_strip=None):
        pass

    @abc.abstractmethod
    def lstrip(self, to_strip=None):
        pass

    @abc.abstractmethod
    def rstrip(self, to_strip=None):
        pass

    @abc.abstractmethod
    def split(self, pat=None, n=-1, expand=False):
        pass

    @abc.abstractmethod
    def rsplit(self, pat=None, n=-1):
        pass
