import numpy as np

from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import registry
from pandas import compat


class SparseDtype(ExtensionDtype):
    """
    Dtype for data stored in :class:`SparseArray`.

    This dtype implements the pandas ExtensionDtype interface.

    .. versionadded:: 0.24.0

    Parameters
    ----------
    dtype : numpy.dtype, default numpy.float64
        The dtype of the underlying array storing the non-fill value values.
    fill_value : scalar, optional.
        The scalar value not stored in the SparseArray. By default, this
        depends on `dtype`.

        ========== ==========
        dtype      na_value
        ========== ==========
        float      ``np.nan``
        int        ``0``
        bool       False
        datetime64 ``pd.NaT``
        ========== ==========

        The default value may be overridden by specifying a `fill_value`.
    """

    def __init__(self, dtype=np.float64, fill_value=None):
        from pandas.core.dtypes.missing import na_value_for_dtype

        if isinstance(dtype, type(self)):
            dtype = dtype.subtype
        else:
            dtype = np.dtype(dtype)

        if fill_value is None:
            fill_value = na_value_for_dtype(dtype)

        self._dtype = dtype
        self._fill_value = fill_value

    def __hash__(self):
        # XXX: this needs to be part of the interface.
        return hash(str(self))

    def __eq__(self, other):
        # TODO: test
        if isinstance(other, type(self)):
            return (self.subtype == other.subtype and
                    self._is_na_fill_value is other._is_na_fill_value)
        else:
            return super(SparseDtype, self).__eq__(other)

    @property
    def fill_value(self):
        return self._fill_value

    @property
    def _is_na_fill_value(self):
        from pandas.core.dtypes.missing import isna
        return isna(self.fill_value)

    @property
    def _is_numeric(self):
        from pandas.core.dtypes.common import is_object_dtype
        return not is_object_dtype(self.subtype)

    @property
    def kind(self):
        return self.subtype.kind

    @property
    def type(self):
        return self.subtype.type

    @property
    def subtype(self):
        return self._dtype

    @property
    def name(self):
        return 'Sparse[{}]'.format(self.subtype.name)

    def __repr__(self):
        return 'Sparse[{},{}]'.format(self.subtype.name, self.fill_value)

    @classmethod
    def construct_array_type(cls):
        from .array import SparseArray
        return SparseArray

    @classmethod
    def construct_from_string(cls, string):
        if string.startswith("Sparse"):
            sub_type = cls._parse_subtype(string)
            try:
                return SparseDtype(sub_type)
            except Exception:
                raise TypeError
        else:
            raise TypeError

    @staticmethod
    def _parse_subtype(dtype):
        if dtype.startswith("Sparse["):
            sub_type = dtype[7:-1]
        elif dtype == "Sparse":
            sub_type = 'float64'
        else:
            raise ValueError
        return sub_type

    @classmethod
    def is_dtype(cls, dtype):
        dtype = getattr(dtype, 'dtype', dtype)
        if (isinstance(dtype, compat.string_types) and
                dtype.startswith("Sparse")):
            dtype = np.dtype(cls._parse_subtype(dtype))
        elif isinstance(dtype, cls):
            return True
        return isinstance(dtype, np.dtype) or dtype == 'Sparse'


registry.register(SparseDtype)
