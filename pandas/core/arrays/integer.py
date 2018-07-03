import sys
import warnings
import copy
import numpy as np

from pandas.compat import u
from pandas.core.dtypes.generic import ABCSeries, ABCIndexClass
from pandas.util._decorators import cache_readonly
from pandas.compat import set_function_name
from pandas.api.types import (is_integer, is_scalar, is_float,
                              is_float_dtype, is_integer_dtype,
                              is_object_dtype,
                              is_list_like,
                              infer_dtype)
from pandas.core.arrays import ExtensionArray, ExtensionOpsMixin
from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import registry
from pandas.core.dtypes.missing import isna, notna
from pandas.core.dtypes.cast import maybe_downcast_to_dtype

from pandas.io.formats.printing import (
    format_object_summary, format_object_attrs, default_pprint)


class IntegerDtype(ExtensionDtype):
    type = None
    na_value = np.nan

    @cache_readonly
    def is_signed_integer(self):
        return self.kind == 'i'

    @cache_readonly
    def is_unsigned_integer(self):
        return self.kind == 'u'

    @cache_readonly
    def numpy_dtype(self):
        """ Return an instance of our numpy dtype """
        return np.dtype(self.type)

    @cache_readonly
    def kind(self):
        return self.numpy_dtype.kind

    @classmethod
    def construct_array_type(cls):
        """Return the array type associated with this dtype

        Returns
        -------
        type
        """
        return IntegerArray

    @classmethod
    def construct_from_string(cls, string):
        """
        Construction from a string, raise a TypeError if not
        possible
        """
        if string == cls.name:
            return cls()
        raise TypeError("Cannot construct a '{}' from "
                        "'{}'".format(cls, string))


def to_integer_array(values):
    """
    Parameters
    ----------
    values : 1D list-like

    Returns
    -------
    infer and return an integer array

    Raises
    ------
    TypeError if incompatible types
    """
    values = np.array(values, copy=False)
    try:
        dtype = _dtypes[str(values.dtype)]
    except KeyError:
        if is_float_dtype(values):
            return IntegerArray(values)

        raise TypeError("Incompatible dtype for {}".format(values.dtype))
    return IntegerArray(values, dtype=dtype, copy=False)


def coerce_to_array(values, dtype, mask=None, copy=False):
    """
    Coerce the input values array to numpy arrays with a mask

    Parameters
    ----------
    values : 1D list-like
    dtype : integer dtype
    mask : boolean 1D array, optional
    copy : boolean, default False
        if True, copy the input

    Returns
    -------
    tuple of (values, mask)
    """

    if isinstance(values, IntegerArray):
        values, mask = values.data, values.mask
        if copy:
            values = values.copy()
            mask = mask.copy()
        return values, mask

    values = np.array(values, copy=copy)
    if is_object_dtype(values):
        inferred_type = infer_dtype(values)
        if inferred_type not in ['floating', 'integer',
                                 'mixed-integer', 'mixed-integer-float']:
            raise TypeError("{} cannot be converted to an IntegerDtype".format(
                values.dtype))

    elif not (is_integer_dtype(values) or is_float_dtype(values)):
        raise TypeError("{} cannot be converted to an IntegerDtype".format(
            values.dtype))

    if mask is None:
        mask = isna(values)
    else:
        assert len(mask) == len(values)

    if not values.ndim == 1:
        raise TypeError("values must be a 1D list-like")
    if not mask.ndim == 1:
        raise TypeError("mask must be a 1D list-like")

    # avoid float->int numpy conversion issues
    if is_object_dtype(values):
        mask |= isna(values)

    # infer dtype if needed
    if dtype is None:
        if is_integer_dtype(values):
            dtype = values.dtype
        else:
            dtype = np.dtype('int64')
    else:
        dtype = dtype.type

    # we copy as need to coerce here
    if mask.any():
        values = values.copy()
        values[mask] = 1

        values = values.astype(dtype)
    else:
        values = values.astype(dtype, copy=False)

    return values, mask


class IntegerArray(ExtensionArray, ExtensionOpsMixin):
    """
    We represent an IntegerArray with 2 numpy arrays
    - data: contains a numpy integer array of the appropriate dtype
    - mask: a boolean array holding a mask on the data, False is missing
    """

    @cache_readonly
    def dtype(self):
        return _dtypes[str(self.data.dtype)]

    def __init__(self, values, mask=None, dtype=None, copy=False):
        self.data, self.mask = coerce_to_array(
            values, dtype=dtype, mask=mask, copy=copy)

    @classmethod
    def _from_sequence(cls, scalars, mask=None, dtype=None, copy=False):
        return cls(scalars, mask=mask, dtype=dtype, copy=copy)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls(values, dtype=original.dtype)

    def __getitem__(self, item):
        if is_integer(item):
            if self.mask[item]:
                return self.dtype.na_value
            return self.data[item]
        return type(self)(self.data[item],
                          mask=self.mask[item],
                          dtype=self.dtype)

    def _coerce_to_ndarray(self):
        """ coerce to an ndarary, preserving my scalar types """

        # TODO(jreback) make this better
        data = self.data.astype(object)
        data[self.mask] = self._na_value
        return data

    def __array__(self, dtype=None):
        """
        the array interface, return my values
        We return an object array here to preserve our scalar values
        """
        return self._coerce_to_ndarray()

    def __iter__(self):
        """Iterate over elements of the array.

        """
        # This needs to be implemented so that pandas recognizes extension
        # arrays as list-like. The default implementation makes successive
        # calls to ``__getitem__``, which may be slower than necessary.
        for i in range(len(self)):
            if self.mask[i]:
                yield self.dtype.na_value
            else:
                yield self.data[i]

    def _formatting_values(self):
        # type: () -> np.ndarray
        return self._coerce_to_ndarray()

    def take(self, indexer, allow_fill=False, fill_value=None):
        from pandas.api.extensions import take

        # we always fill with 1 internally
        # to avoid upcasting
        data_fill_value = 1 if isna(fill_value) else fill_value
        result = take(self.data, indexer, fill_value=data_fill_value,
                      allow_fill=allow_fill)

        mask = take(self.mask, indexer, fill_value=True,
                    allow_fill=allow_fill)

        # if we are filling
        # we only fill where the indexer is null
        # not existing missing values
        # TODO(jreback) what if we have a non-na float as a fill value?
        if allow_fill and notna(fill_value):
            fill_mask = np.asarray(indexer) == -1
            result[fill_mask] = fill_value
            mask = mask ^ fill_mask

        return type(self)(result, mask=mask, dtype=self.dtype)

    def copy(self, deep=False):
        data, mask = self.data, self.mask
        if deep:
            data = copy.deepcopy(data)
            mask = copy.deepcopy(mask)
        else:
            data = data.copy()
            mask = mask.copy()
        return type(self)(data, mask, dtype=self.dtype, copy=False)

    def __setitem__(self, key, value):
        _is_scalar = is_scalar(value)
        if _is_scalar:
            value = [value]
        value, mask = coerce_to_array(value, dtype=self.dtype)

        if _is_scalar:
            value = value[0]
            mask = mask[0]

        self.data[key] = value
        self.mask[key] = mask

    def __len__(self):
        return len(self.data)

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
        return self.data.nbytes + self.mask.nbytes

    def isna(self):
        return self.mask

    @property
    def _na_value(self):
        return np.nan

    @classmethod
    def _concat_same_type(cls, to_concat):
        data = np.concatenate([x.data for x in to_concat])
        mask = np.concatenate([x.mask for x in to_concat])
        return cls(data, mask=mask, dtype=to_concat[0].dtype)

    def astype(self, dtype, copy=True):
        """Cast to a NumPy array with 'dtype'.

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
        array : ndarray
            NumPy ndarray with 'dtype' for its dtype.

        Raises
        ------
        TypeError
            if incompatible type with an IntegerDtype, equivalent of same_kind
            casting
        """

        # if we are astyping to an existing IntegerDtype we can fastpath
        if isinstance(dtype, IntegerDtype):
            result = self.data.astype(dtype.numpy_dtype,
                                      casting='same_kind', copy=False)
            return type(self)(result, mask=self.mask,
                              dtype=dtype, copy=False)

        # coerce
        data = self._coerce_to_ndarray()
        return data.astype(dtype=dtype, copy=False)

    @property
    def _ndarray_values(self):
        # type: () -> np.ndarray
        """Internal pandas method for lossy conversion to a NumPy ndarray.

        This method is not part of the pandas interface.

        The expectation is that this is cheap to compute, and is primarily
        used for interacting with our indexers.
        """
        return self.data

    def value_counts(self, dropna=True):
        """
        Returns a Series containing counts of each category.

        Every category will have an entry, even those with a count of 0.

        Parameters
        ----------
        dropna : boolean, default True
            Don't include counts of NaN.

        Returns
        -------
        counts : Series

        See Also
        --------
        Series.value_counts

        """

        from pandas import Index, Series

        # compute counts on the data with no nans
        data = self.data[~self.mask]
        value_counts = Index(data).value_counts()
        array = value_counts.values

        # TODO(extension)
        # if we have allow Index to hold an ExtensionArray
        # this is easier
        index = value_counts.index.astype(object)

        # if we want nans, count the mask
        if not dropna:

            # TODO(extension)
            # appending to an Index *always* infers
            # w/o passing the dtype
            array = np.append(array, [self.mask.sum()])
            index = Index(np.concatenate(
                [index.values,
                 np.array([np.nan], dtype=object)]), dtype=object)

        return Series(array, index=index)

    def _values_for_argsort(self):
        # type: () -> ndarray
        """Return values for sorting.

        Returns
        -------
        ndarray
            The transformed values should maintain the ordering between values
            within the array.

        See Also
        --------
        ExtensionArray.argsort
        """
        data = self.data.copy()
        data[self.mask] = data.min() - 1
        return data

    @classmethod
    def _create_comparison_method(cls, op):
        def cmp_method(self, other):

            op_name = op.__name__
            mask = None
            if isinstance(other, IntegerArray):
                other, mask = other.data, other.mask
            elif is_list_like(other):
                other = np.asarray(other)
                if other.ndim > 0 and len(self) != len(other):
                    raise ValueError('Lengths must match to compare')

            # numpy will show a DeprecationWarning on invalid elementwise
            # comparisons, this will raise in the future
            with warnings.catch_warnings(record=True):
                with np.errstate(all='ignore'):
                    result = op(self.data, other)

            # nans propagate
            if mask is None:
                mask = self.mask
            else:
                mask = self.mask | mask

            result[mask] = True if op_name == 'ne' else False
            return result

        name = '__{name}__'.format(name=op.__name__)
        return set_function_name(cmp_method, name, cls)

    def _maybe_mask_result(self, result, mask, other, op_name):
        """
        Parameters
        ----------
        result : array-like
        mask : array-like bool
        other : scalar or array-like
        op_name : str
        """

        # may need to fill infs
        # and mask wraparound
        if is_float_dtype(result):
            mask |= (result == np.inf) | (result == -np.inf)

        # floor div can be a float or an integer dependending
        # on the operands
        if (op_name in ['rfloordiv', 'floordiv'] and
                (is_float_dtype(other) or is_float(other))):
            result[mask] = np.nan
            return result

        # by definition a float result
        elif op_name in ['rtruediv', 'truediv', 'rdiv', 'div']:
            result[mask] = np.nan
            return result

        elif is_float_dtype(result):
            # if our float result, try to downcast if possible
            # if remains float, then mask and return as float
            nonans = result[notna(result)]
            maybe = maybe_downcast_to_dtype(nonans, self.dtype.numpy_dtype)
            if not is_integer_dtype(maybe):
                result[mask] = np.nan
                return result

        return type(self)(result, mask=mask, dtype=self.dtype, copy=False)

    @classmethod
    def _create_arithmetic_method(cls, op):
        def integer_arithmetic_method(self, other):

            op_name = op.__name__
            mask = None
            if isinstance(other, (ABCSeries, ABCIndexClass)):
                other = getattr(other, 'values', other)

            if isinstance(other, IntegerArray):
                other, mask = other.data, other.mask
            elif getattr(other, 'ndim', 0) > 1:
                raise TypeError("can only perform ops with 1-d structures")
            elif is_list_like(other):
                other = np.asarray(other)
                if not other.ndim:
                    other = other.item()
                elif other.ndim == 1:
                    if not (is_float_dtype(other) or is_integer_dtype(other)):
                        raise TypeError(
                            "can only perform ops with numeric values")
            else:
                if not (is_float(other) or is_integer(other)):
                    raise TypeError("can only perform ops with numeric values")

            # nans propagate
            if mask is None:
                mask = self.mask
            else:
                mask = self.mask | mask

            with np.errstate(all='ignore'):
                result = op(self.data, other)

            # divmod returns a tuple
            if op_name == 'divmod':
                div, mod = result
                return (self._maybe_mask_result(div, mask, other, 'floordiv'),
                        self._maybe_mask_result(mod, mask, other, 'mod'))

            return self._maybe_mask_result(result, mask, other, op_name)

        name = '__{name}__'.format(name=op.__name__)
        return set_function_name(integer_arithmetic_method, name, cls)


IntegerArray._add_arithmetic_ops()
IntegerArray._add_comparison_ops()


module = sys.modules[__name__]


# create the Dtype
_dtypes = {}
for dtype in ['int8', 'int16', 'int32', 'int64',
              'uint8', 'uint16', 'uint32', 'uint64']:

    if dtype.startswith('u'):
        name = "U{}".format(dtype[1:].capitalize())
    else:
        name = dtype.capitalize()
    classname = "{}Dtype".format(name)
    attributes_dict = {'type': getattr(np, dtype),
                       'name': name}
    dtype_type = type(classname, (IntegerDtype, ), attributes_dict)
    setattr(module, classname, dtype_type)

    # register
    registry.register(dtype_type)
    _dtypes[dtype] = dtype_type()
