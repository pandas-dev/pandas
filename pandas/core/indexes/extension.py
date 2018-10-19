import numpy as np
from pandas._libs import index as libindex

# from pandas._libs import (lib, index as libindex, tslibs,
#                           algos as libalgos, join as libjoin,
#                           Timedelta)

from pandas.compat.numpy import function as nv

from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.common import (
    pandas_dtype,
    ensure_platform_int,
    is_dtype_equal,
    is_integer_dtype,
    is_float_dtype,
    is_extension_array_dtype)
from pandas.core.dtypes.generic import (
    ABCSeries, ABCIndex
)
from pandas.util._decorators import (
    Appender, cache_readonly)

from .base import Index


# _index_doc_kwargs = dict(ibase._index_doc_kwargs)
# _index_doc_kwargs.update(
#     dict(klass='IntervalIndex',
#          target_klass='IntervalIndex or list of Intervals',
#          name=textwrap.dedent("""\
#          name : object, optional
#               to be stored in the index.
#          """),
#          ))


class ExtensionIndex(Index):
    """
    Index class that holds an ExtensionArray.

    """
    _typ = 'extensionindex'
    _comparables = ['name']
    _attributes = ['name']

    _can_hold_na = True

    @property
    def _is_numeric_dtype(self):
        return self.dtype._is_numeric

    # TODO
    # # would we like our indexing holder to defer to us
    # _defer_to_indexing = False

    # # prioritize current class for _shallow_copy_with_infer,
    # # used to infer integers as datetime-likes
    # _infer_as_myclass = False

    def __new__(cls, *args, **kwargs):
        return object.__new__(cls)

    def __init__(self, data, dtype=None, name=None, copy=False, **kwargs):
        # needs to accept and ignore kwargs eg for freq passed in
        # Index._shallow_copy_with_infer

        # unbox containers that can contain ExtensionArray
        if isinstance(data, (ABCSeries, ABCIndex)):
            data = data._values

        # check dtype and coerce data to dtype if needed
        if dtype is not None:
            dtype = pandas_dtype(dtype)
            if not is_extension_array_dtype(dtype):
                raise ValueError(
                    "The passed dtype should be an ExtensionDtype")
            if not is_dtype_equal(getattr(data, 'dtype', None), dtype):
                data = dtype.construct_array_type()._from_sequence(
                    data, dtype=dtype, copy=False)

        if not isinstance(data, ExtensionArray):
            raise ValueError("passed data should be an ExtensionArray, or the "
                             "passed dtype should be an ExtensionDtype")

        if copy:
            data = data.copy()

        self._data = data
        self.name = name

    def __len__(self):
        """
        return the length of the Index
        """
        return len(self._data)

    @property
    def size(self):
        # EA does not have .size
        return len(self._data)

    def __array__(self, dtype=None):
        """ the array interface, return my values """
        return np.array(self._data)

    @cache_readonly
    def dtype(self):
        """ return the dtype object of the underlying data """
        return self._values.dtype

    @cache_readonly
    def dtype_str(self):
        """ return the dtype str of the underlying data """
        return str(self.dtype)

    @property
    def _values(self):
        return self._data

    @property
    def values(self):
        """ return the underlying data as an ndarray """
        return self._values

    @cache_readonly
    def _isnan(self):
        """ return if each value is nan"""
        return self._values.isna()

    @cache_readonly
    def _engine_type(self):
        values, na_value = self._values._values_for_factorize()
        if is_integer_dtype(values):
            return libindex.Int64Engine
        elif is_float_dtype(values):
            return libindex.Float64Engine
        # TODO add more
        else:
            return libindex.ObjectEngine

    @cache_readonly
    def _engine(self):
        # property, for now, slow to look up
        values, na_value = self._values._values_for_factorize()
        return self._engine_type(lambda: values, len(self))

    def _format_with_header(self, header, **kwargs):
        return header + list(self._format_native_types(**kwargs))

    @Appender(Index.take.__doc__)
    def take(self, indices, axis=0, allow_fill=True, fill_value=None,
             **kwargs):
        if kwargs:
            nv.validate_take(tuple(), kwargs)
        indices = ensure_platform_int(indices)

        result = self._data.take(indices, allow_fill=allow_fill,
                                 fill_value=fill_value)
        attributes = self._get_attributes_dict()
        return self._simple_new(result, **attributes)

    def __getitem__(self, value):
        result = self._data[value]
        if isinstance(result, self._data.__class__):
            return self._shallow_copy(result)
        else:
            # scalar
            return result
