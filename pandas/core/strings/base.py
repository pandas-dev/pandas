from __future__ import annotations

import abc
from collections.abc import Callable  # noqa: PDF001
import re
from typing import TYPE_CHECKING

import numpy as np

from pandas._typing import Scalar

if TYPE_CHECKING:
    from pandas import Series


class BaseStringArrayMethods(abc.ABC):
    """
    Base class for extension arrays implementing string methods.

    This is where our ExtensionArrays can override the implementation of
    Series.str.<method>. We don't expect this to work with
    3rd-party extension arrays.

    * User calls Series.str.<method>
    * pandas extracts the extension array from the Series
    * pandas calls ``extension_array._str_<method>(*args, **kwargs)``
    * pandas wraps the result, to return to the user.

    See :ref:`Series.str` for the docstring of each method.
    """

    def _str_getitem(self, key):
        if isinstance(key, slice):
            return self._str_slice(start=key.start, stop=key.stop, step=key.step)
        else:
            return self._str_get(key)

    @abc.abstractmethod
    def _str_count(self, pat, flags=0) -> None:
        pass

    @abc.abstractmethod
    def _str_pad(self, width, side="left", fillchar=" ") -> None:
        pass

    @abc.abstractmethod
    def _str_contains(self, pat, case=True, flags=0, na=None, regex=True) -> None:
        pass

    @abc.abstractmethod
    def _str_startswith(self, pat, na=None) -> None:
        pass

    @abc.abstractmethod
    def _str_endswith(self, pat, na=None) -> None:
        pass

    @abc.abstractmethod
    def _str_replace(
        self,
        pat: str | re.Pattern,
        repl: str | Callable,
        n: int = -1,
        case: bool = True,
        flags: int = 0,
        regex: bool = True,
    ) -> None:
        pass

    @abc.abstractmethod
    def _str_repeat(self, repeats) -> None:
        pass

    @abc.abstractmethod
    def _str_match(
        self, pat: str, case: bool = True, flags: int = 0, na: Scalar = np.nan
    ) -> None:
        pass

    @abc.abstractmethod
    def _str_fullmatch(
        self,
        pat: str | re.Pattern,
        case: bool = True,
        flags: int = 0,
        na: Scalar = np.nan,
    ) -> None:
        pass

    @abc.abstractmethod
    def _str_encode(self, encoding, errors="strict") -> None:
        pass

    @abc.abstractmethod
    def _str_find(self, sub, start=0, end=None) -> None:
        pass

    @abc.abstractmethod
    def _str_rfind(self, sub, start=0, end=None) -> None:
        pass

    @abc.abstractmethod
    def _str_findall(self, pat, flags=0) -> None:
        pass

    @abc.abstractmethod
    def _str_get(self, i) -> None:
        pass

    @abc.abstractmethod
    def _str_index(self, sub, start=0, end=None) -> None:
        pass

    @abc.abstractmethod
    def _str_rindex(self, sub, start=0, end=None) -> None:
        pass

    @abc.abstractmethod
    def _str_join(self, sep) -> None:
        pass

    @abc.abstractmethod
    def _str_partition(self, sep, expand) -> None:
        pass

    @abc.abstractmethod
    def _str_rpartition(self, sep, expand) -> None:
        pass

    @abc.abstractmethod
    def _str_len(self) -> None:
        pass

    @abc.abstractmethod
    def _str_slice(self, start=None, stop=None, step=None) -> None:
        pass

    @abc.abstractmethod
    def _str_slice_replace(self, start=None, stop=None, repl=None) -> None:
        pass

    @abc.abstractmethod
    def _str_translate(self, table) -> None:
        pass

    @abc.abstractmethod
    def _str_wrap(self, width, **kwargs) -> None:
        pass

    @abc.abstractmethod
    def _str_get_dummies(self, sep="|") -> None:
        pass

    @abc.abstractmethod
    def _str_isalnum(self) -> None:
        pass

    @abc.abstractmethod
    def _str_isalpha(self) -> None:
        pass

    @abc.abstractmethod
    def _str_isdecimal(self) -> None:
        pass

    @abc.abstractmethod
    def _str_isdigit(self) -> None:
        pass

    @abc.abstractmethod
    def _str_islower(self) -> None:
        pass

    @abc.abstractmethod
    def _str_isnumeric(self) -> None:
        pass

    @abc.abstractmethod
    def _str_isspace(self) -> None:
        pass

    @abc.abstractmethod
    def _str_istitle(self) -> None:
        pass

    @abc.abstractmethod
    def _str_isupper(self) -> None:
        pass

    @abc.abstractmethod
    def _str_capitalize(self) -> None:
        pass

    @abc.abstractmethod
    def _str_casefold(self) -> None:
        pass

    @abc.abstractmethod
    def _str_title(self) -> None:
        pass

    @abc.abstractmethod
    def _str_swapcase(self) -> None:
        pass

    @abc.abstractmethod
    def _str_lower(self) -> None:
        pass

    @abc.abstractmethod
    def _str_upper(self) -> None:
        pass

    @abc.abstractmethod
    def _str_normalize(self, form) -> None:
        pass

    @abc.abstractmethod
    def _str_strip(self, to_strip=None) -> None:
        pass

    @abc.abstractmethod
    def _str_lstrip(self, to_strip=None) -> None:
        pass

    @abc.abstractmethod
    def _str_rstrip(self, to_strip=None) -> None:
        pass

    @abc.abstractmethod
    def _str_removeprefix(self, prefix: str) -> Series:
        pass

    @abc.abstractmethod
    def _str_removesuffix(self, suffix: str) -> Series:
        pass

    @abc.abstractmethod
    def _str_split(self, pat=None, n=-1, expand=False) -> None:
        pass

    @abc.abstractmethod
    def _str_rsplit(self, pat=None, n=-1) -> None:
        pass

    @abc.abstractmethod
    def _str_extract(self, pat: str, flags: int = 0, expand: bool = True) -> None:
        pass
