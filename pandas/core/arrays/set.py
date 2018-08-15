import sys
import warnings
import copy
import numpy as np

import operator

from pandas import Series

from pandas._libs.lib import infer_dtype
from pandas.util._decorators import cache_readonly
from pandas.compat import u, range
from pandas.compat import set_function_name

from pandas.core.dtypes.generic import ABCSeries, ABCIndexClass
from pandas.core.dtypes.common import (
    is_integer, is_scalar, is_float,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_list_like)
from pandas.core.arrays import ExtensionArray, ExtensionOpsMixin, ExtensionScalarOpsMixin
from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import registry
from pandas.core.dtypes.missing import isna, notna

from pandas.io.formats.printing import (
    format_object_summary, format_object_attrs, default_pprint)


class SetDtype(ExtensionDtype):
    """
    An ExtensionDtype to hold sets.
    """
    name = 'set'
    type = object
    na_value = np.nan

    def __hash__(self):
        # XXX: this needs to be part of the interface.
        return hash(str(self))

    def __eq__(self, other):
        # TODO: test
        if isinstance(other, type(self)):
            return True
        else:
            return super(SetDtype, self).__eq__(other)

    @property
    def _is_numeric(self):
        return False

    def __repr__(self):
        return self.name

    @classmethod
    def construct_array_type(cls):
        """Return the array type associated with this dtype

        Returns
        -------
        type
        """
        return SetArray

    @classmethod
    def construct_from_string(cls, string):
        """
        Construction from a string, raise a TypeError if not
        possible
        """
        if string == cls.name or string is set:
            return cls()
        raise TypeError("Cannot construct a '{}' from "
                        "'{}'".format(cls, string))

#     @classmethod
#     def is_dtype(cls, dtype):
#         dtype = getattr(dtype, 'dtype', dtype)
#         if (isinstance(dtype, compat.string_types) and
#                 dtype == 'set'):
#             return True
#         elif isinstance(dtype, cls):
#             return True
#         return isinstance(dtype, np.dtype) or dtype == 'set'


def to_set_array(values):
    """
    Infer and return a set array of the values.

    Parameters
    ----------
    values : 1D list-like of list-likes

    Returns
    -------
    SetArray

    Raises
    ------
    TypeError if incompatible types
    """
    return SetArray(values, copy=False)


def coerce_to_array(values, mask=None, copy=False):
    """
    Coerce the input values array to numpy arrays with a mask

    Parameters
    ----------
    values : 1D list-like
    mask : boolean 1D array, optional
    copy : boolean, default False
        if True, copy the input

    Returns
    -------
    tuple of (values, mask)
    """

    if isinstance(values, SetArray):
        values, mask = values._data, values._mask

        if copy:
            values = values.copy()
            mask = mask.copy()
        return values, mask

    values = np.array(values, copy=copy)
    if not (is_object_dtype(values) or isna(values).all()):
        raise TypeError("{} cannot be converted to a SetDtype".format(
            values.dtype))

    if mask is None:
        mask = isna(values)
    else:
        assert len(mask) == len(values)

    if not values.ndim == 1:
        raise TypeError("values must be a 1D list-like")
    if not mask.ndim == 1:
        raise TypeError("mask must be a 1D list-like")

    if mask.any():
        values = values.copy()
        values[mask] = np.nan

    return values, mask


class SetArray(ExtensionArray, ExtensionOpsMixin):
    """
    We represent a SetArray with 2 numpy arrays
    - data: contains a numpy set array of object dtype
    - mask: a boolean array holding a mask on the data, False is missing
    """

    @cache_readonly
    def dtype(self):
        return SetDtype()

    def __init__(self, values, mask=None, dtype=None, copy=False):
        """
        Parameters
        ----------
        values : 1D list-like / SetArray
        mask : 1D list-like, optional
        copy : bool, default False

        Returns
        -------
        SetArray
        """
        self._data, self._mask = coerce_to_array(
            values, mask=mask, copy=copy)

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        # dtype is ignored
        return cls(scalars, copy=copy)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls(values)

    def __getitem__(self, item):
        if is_integer(item):
            if self._mask[item]:
                return self.dtype.na_value
            return self._data[item]
        return type(self)(self._data[item], mask=self._mask[item])

    def _coerce_to_ndarray(self):
        """
        coerce to an ndarray of object dtype
        """
        data = self._data
        data[self._mask] = self._na_value
        return data

    def __array__(self):
        """
        the array interface, return values
        """
        return self._coerce_to_ndarray()

    def __iter__(self):
        """Iterate over elements of the array.

        """
        # This needs to be implemented so that pandas recognizes extension
        # arrays as list-like. The default implementation makes successive
        # calls to ``__getitem__``, which may be slower than necessary.
        for i in range(len(self)):
            if self._mask[i]:
                yield self.dtype.na_value
            else:
                yield self._data[i]

    def _formatting_values(self):
        # type: () -> np.ndarray
        return self._coerce_to_ndarray()

    def take(self, indices, allow_fill=False, fill_value=None):
        from pandas.core.algorithms import take

        if allow_fill and fill_value is None:
            fill_value = self.dtype.na_value

        result = take(self._data, indices, fill_value=fill_value,
                      allow_fill=allow_fill)
        return self._from_sequence(result)

    def copy(self, deep=False):
        data, mask = self._data, self._mask
        if deep:
            data = copy.deepcopy(data)
            mask = copy.deepcopy(mask)
        else:
            data = data.copy()
            mask = mask.copy()
        return type(self)(data, mask, copy=False)

    def __setitem__(self, key, value):
        _is_scalar = is_scalar(value)
        if _is_scalar:
            value = [value]
        value, mask = coerce_to_array(value)

        if _is_scalar:
            value = value[0]
            mask = mask[0]

        self._data[key] = value
        self._mask[key] = mask

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        """
        Return a string representation for this object.

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        klass = self.__class__.__name__
        data = format_object_summary(self, default_pprint, False)
        attrs = format_object_attrs(self)
        space = " "

        prepr = (u(",%s") %
                 space).join(u("%s=%s") % (k, v) for k, v in attrs)

        res = u("%s(%s%s)") % (klass, data, prepr)

        return res

    @property
    def nbytes(self):
        return self._data.nbytes + self._mask.nbytes

    def isna(self):
        return self._mask

    @property
    def _na_value(self):
        return np.nan

    @classmethod
    def _concat_same_type(cls, to_concat):
        data = np.concatenate([x._data for x in to_concat])
        mask = np.concatenate([x._mask for x in to_concat])
        return cls(data, mask=mask)

    def astype(self, dtype, copy=True, errors='raise', fill_value=None):
        """Cast to a NumPy array or SetArray with 'dtype'.

        Parameters
        ----------
        dtype : str or dtype
            Typecode or data-type to which the array is cast.
        copy : bool, default True
            Whether to copy the data, even if not necessary. If False,
            a copy is made only if the old dtype does not match the
            new dtype.

        Returns
        -------
        array : ndarray or SetArray
            NumPy ndarray or SetArray with 'dtype' for its dtype.

        Raises
        ------
        TypeError
            if incompatible type with a SetDtype, equivalent of same_kind
            casting
        """

        # if we are astyping to an existing IntegerDtype we can fastpath
        if isinstance(dtype, SetDtype):
            result = self._data.astype(dtype.type,
                                       casting='same_kind', copy=False)
            return type(self)(result, mask=self._mask, copy=False)

        # coerce
        data = self._coerce_to_ndarray()
        return data.astype(dtype, copy=False)

    @property
    def _ndarray_values(self):
        # type: () -> np.ndarray
        """Internal pandas method for lossy conversion to a NumPy ndarray.

        This method is not part of the pandas interface.

        The expectation is that this is cheap to compute, and is primarily
        used for interacting with our indexers.
        """
        return self._data

    def fillna(self, value=None, method=None, limit=None):
        res = self._data.copy()
        res[self._mask] = [value] * self._mask.sum()
        return type(self)(res,
                          mask=np.full_like(res, fill_value=False, dtype=bool),
                          copy=False)

    def dropna(self):
        res = self._data[~self._mask]
        return type(self)(res,
                          mask=np.full_like(res, fill_value=False, dtype=bool),
                          copy=False)

    def unique(self):
        raise NotImplementedError

    def factorize(self):
        raise NotImplementedError

    def argsort(self):
        raise NotImplementedError

#     def value_counts(self, dropna=True):
#         """
#         Returns a Series containing counts of each category.
# 
#         Every category will have an entry, even those with a count of 0.
# 
#         Parameters
#         ----------
#         dropna : boolean, default True
#             Don't include counts of NaN.
# 
#         Returns
#         -------
#         counts : Series
# 
#         See Also
#         --------
#         Series.value_counts
# 
#         """
# 
#         from pandas import Index, Series
# 
#         # compute counts on the data with no nans
#         data = self._data[~self._mask]
#         value_counts = Index(data).value_counts()
#         array = value_counts.values
# 
#         # TODO(extension)
#         # if we have allow Index to hold an ExtensionArray
#         # this is easier
#         index = value_counts.index.astype(object)
# 
#         # if we want nans, count the mask
#         if not dropna:
# 
#             # TODO(extension)
#             # appending to an Index *always* infers
#             # w/o passing the dtype
#             array = np.append(array, [self._mask.sum()])
#             index = Index(np.concatenate(
#                 [index.values,
#                  np.array([np.nan], dtype=object)]), dtype=object)
# 
#         return Series(array, index=index)

#     def _values_for_argsort(self):
#         # type: () -> ndarray
#         """Return values for sorting.
# 
#         Returns
#         -------
#         ndarray
#             The transformed values should maintain the ordering between values
#             within the array.
# 
#         See Also
#         --------
#         ExtensionArray.argsort
#         """
#         data = self._data.copy()
#         data[self._mask] = data.min() - 1
#         return data

    @classmethod
    def _create_comparison_method(cls, op):
        def cmp_method(self, other):

            op_name = op.__name__
            mask = None
            if isinstance(other, SetArray):
                other, mask = other._data, other._mask
            elif (isinstance(other, Series)
                  and isinstance(other.values, SetArray)):
                other, mask = other.values._data, other.values._mask
            elif isinstance(other, set) or (is_scalar(other) and isna(other)):
                other = np.array([other] * len(self))
            elif is_list_like(other):
                other = np.asarray(other)
                if other.ndim > 0 and len(self) != len(other):
                    raise ValueError('Lengths must match to compare')

            mask = self._mask | mask if mask is not None else self._mask
            result = np.full_like(self._data, fill_value=np.nan, dtype='O')

            # numpy will show a DeprecationWarning on invalid elementwise
            # comparisons, this will raise in the future
            with warnings.catch_warnings(record=True):
                with np.errstate(all='ignore'):
                    result[~mask] = op(self._data[~mask], other[~mask])

            result[mask] = True if op_name == 'ne' else False
            return result.astype('bool')

        name = '__{name}__'.format(name=op.__name__)
        return set_function_name(cmp_method, name, cls)

    @classmethod
    def _create_arithmetic_method(cls, op):
        def arithmetic_method(self, other):

            op_name = op.__name__
            mask = None
            #print(other)
            if isinstance(other, SetArray):
                other, mask = other._data, other._mask
            elif isinstance(other, set) or (is_scalar(other) and isna(other)):
                other = np.array([other] * len(self))
            elif is_list_like(other):
                other = np.asarray(other)
                #print(other)
                # cannot use isnan due to numpy/numpy#9009
                mask = np.array([x is np.nan for x in other])
                if other.ndim > 0 and len(self) != len(other):
                    raise ValueError('Lengths must match to compare')

            mask = self._mask | mask if mask is not None else self._mask
            result = np.full_like(self._data, fill_value=np.nan, dtype='O')
            #print(result[~mask], self._data[~mask], other[~mask])
            #print(type(result), type(self._data), type(other))

            with np.errstate(all='ignore'):
                result[~mask] = op(self._data[~mask], other[~mask])

            return type(self)(result, mask=mask, copy=False)
        
        name = '__{name}__'.format(name=op.__name__)
        def raiser(self, other):
            raise NotImplementedError
        if name != '__sub__':
            return raiser
        return set_function_name(arithmetic_method, name, cls)


SetArray._add_comparison_ops()
SetArray._add_arithmetic_ops()
SetArray.__or__ = SetArray._create_arithmetic_method(operator.__or__)
SetArray.__xor__ = SetArray._create_arithmetic_method(operator.__xor__)
SetArray.__and__ = SetArray._create_arithmetic_method(operator.__and__)

module = sys.modules[__name__]
setattr(module, 'SetDtype', SetDtype)
registry.register(SetDtype)
