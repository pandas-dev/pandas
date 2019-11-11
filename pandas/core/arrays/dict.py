from typing import Dict, Optional, Sequence, Type

import numpy as np

from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import register_extension_dtype

from pandas._typing import Dtype
from pandas.core.arrays import PandasArray
from pandas.core.construction import extract_array


@register_extension_dtype
class DictDtype(ExtensionDtype):
    """
    Extension dtype for nested dictionaries.

    .. versionadded:: 1.0.0

    .. warning::

    DictDtype is considered experimental. The implementation and
    parts of the API may change without warning.

    """

    @property
    def na_value(self) -> float:
        return np.nan

    @property
    def type(self) -> Type:
        return dict

    @property
    def name(self) -> str:
        """
        The alias for DictDtype is ``'dict'``.
        """
        return "dict"

    @classmethod
    def construct_from_string(cls, string: str) -> ExtensionDtype:
        if string == "dict":
            return cls()

        return super().construct_from_string(string)

    @classmethod
    def construct_array_type(cls) -> Type["DictArray"]:
        return DictArray

    def __repr__(self) -> str:
        return "DictDtype"


class DictArray(PandasArray):
    """
    Extension array for nested dictionaries.

    .. versionadded:: 1.0.0

    .. warning::

    DictArray is considered experimental. The implementation and
    parts of the API may change without warning.
    """

    # undo the PandasArray hack
    _typ = "extension"

    def __init__(self, values: Sequence[Dict], copy: bool = False):
        np_values = extract_array(values)
        super().__init__(np_values, copy=copy)
        self._dtype = DictDtype()
        self._validate()

    def _validate(self):
        """Validate that we only store dicts."""
        if self._ndarray.dtype != "object":
            raise ValueError(
                "DictArray requires a sequence of dicts. Got "
                "'{}' dtype instead.".format(self._ndarray.dtype)
            )

    @classmethod
    def _from_sequence(
        cls, dicts: Sequence[Dict], dtype: Optional[Dtype] = None, copy: bool = False
    ) -> "DictArray":
        if dtype:
            assert dtype == "dict"

        result = super()._from_sequence(dicts, dtype=object, copy=copy)
        # convert None to np.nan
        # TODO: it would be nice to do this in _validate / lib.is_string_array
        # We are already doing a scan over the values there.
        result[result.isna()] = np.nan
        return result

    def _reduce(self, name, axis=0, **kwargs):
        raise NotImplementedError


"""
@register_series_accessor("dict")
class DictAccessor:

    def __init__(self, obj: DictArray):
        if not isinstance(obj, DictArray):
            raise AttributeError("Can only use .dict accessor with a DictArray")
        self._obj = obj

    def __getitem__(self, key: Hashable):
        return pd.Series(x[key] for x in self._obj)

    def get(self, key: Hashable, default=np.nan):
        # TODO: justify np.nan - maybe because will coerce that way in
        # resulting Series construction?
        return pd.Series(x.get(key, default=None) for x in self._obj)
"""
