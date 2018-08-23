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
    dtype : str, ExtensionDtype, numpy.dtype, type, default numpy.float64
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
        timedelta64 ``pd.NaT``
        ========== ==========

        The default value may be overridden by specifying a `fill_value`.
    """

    def __init__(self, dtype=np.float64, fill_value=None):
        # type: (Union[str, np.dtype, 'ExtensionDtype', type], Any) -> None
        from pandas.core.dtypes.missing import na_value_for_dtype
        from pandas.core.dtypes.common import pandas_dtype, is_string_dtype
        from pandas.core.dtypes.common import is_scalar

        if isinstance(dtype, type(self)):
            if fill_value is None:
                fill_value = dtype.fill_value
            dtype = dtype.subtype

        dtype = pandas_dtype(dtype)
        if is_string_dtype(dtype):
            dtype = np.dtype('object')

        if fill_value is None:
            fill_value = na_value_for_dtype(dtype)

        if not is_scalar(fill_value):
            raise ValueError("fill_value must be a scalar. Got {} "
                             "instead".format(fill_value))
        self._dtype = dtype
        self._fill_value = fill_value

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, type(self)):
            subtype = self.subtype == other.subtype
            if self._is_na_fill_value:
                # this case is complicated by two things:
                # SparseDtype(float, float(nan)) == SparseDtype(float, np.nan)
                # SparseDtype(float, np.nan)     != SparseDtype(float, pd.NaT)
                # i.e. we want to treat any floating-point NaN as equal, but
                # not a floating-point NaN and a datetime NaT.
                fill_value = (
                    other._is_na_fill_value and
                    isinstance(self.fill_value, type(other.fill_value)) or
                    isinstance(other.fill_value, type(self.fill_value))
                )
            else:
                fill_value = self.fill_value == other.fill_value

            return subtype and fill_value
        else:
            return super(SparseDtype, self).__eq__(other)

    @property
    def fill_value(self):
        """
        The fill value of the array.

        Converting the SparseArray to a dense ndarray will fill the
        array with this value.

        .. warning::

           It's possible to end up with a SparseArray that has ``fill_value``
           values in ``sp_values``. This can occur, for example, when setting
           ``SparseArray.fill_value`` directly.
        """
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
        msg = "Could not construct SparseDtype from '{}'".format(string)
        if string.startswith("Sparse"):
            sub_type = cls._parse_subtype(string)
            try:
                return SparseDtype(sub_type)
            except Exception:
                raise TypeError(msg)
        else:
            raise TypeError(msg)

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
