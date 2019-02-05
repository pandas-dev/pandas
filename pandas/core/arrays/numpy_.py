import numbers

import numpy as np

from pandas._libs import lib
from pandas.compat.numpy import function as nv
from pandas.util._validators import validate_fillna_kwargs

from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import ABCIndexClass, ABCSeries
from pandas.core.dtypes.inference import is_array_like, is_list_like

from pandas import compat
from pandas.core import nanops
from pandas.core.missing import backfill_1d, pad_1d

from .base import ExtensionArray, ExtensionOpsMixin


class PandasDtype(ExtensionDtype):
    """
    A Pandas ExtensionDtype for NumPy dtypes.

    .. versionadded:: 0.24.0

    This is mostly for internal compatibility, and is not especially
    useful on its own.

    Parameters
    ----------
    dtype : numpy.dtype
    """
    _metadata = ('_dtype',)

    def __init__(self, dtype):
        dtype = np.dtype(dtype)
        self._dtype = dtype
        self._name = dtype.name
        self._type = dtype.type

    def __repr__(self):
        return "PandasDtype({!r})".format(self.name)

    @property
    def numpy_dtype(self):
        """The NumPy dtype this PandasDtype wraps."""
        return self._dtype

    @property
    def name(self):
        return self._name

    @property
    def type(self):
        return self._type

    @property
    def _is_numeric(self):
        # exclude object, str, unicode, void.
        return self.kind in set('biufc')

    @property
    def _is_boolean(self):
        return self.kind == 'b'

    @classmethod
    def construct_from_string(cls, string):
        return cls(np.dtype(string))

    def construct_array_type(cls):
        return PandasArray

    @property
    def kind(self):
        return self._dtype.kind

    @property
    def itemsize(self):
        """The element size of this data-type object."""
        return self._dtype.itemsize


# TODO(NumPy1.13): remove this
# Compat for NumPy 1.12, which doesn't provide NDArrayOperatorsMixin
# or __array_ufunc__, so those operations won't be available to people
# on older NumPys.
#
# We would normally write this as bases=(...), then "class Foo(*bases):
# but Python2 doesn't allow unpacking tuples in the class statement.
# So, we fall back to "object", to avoid writing a metaclass.
try:
    from numpy.lib.mixins import NDArrayOperatorsMixin
except ImportError:
    NDArrayOperatorsMixin = object


class PandasArray(ExtensionArray, ExtensionOpsMixin, NDArrayOperatorsMixin):
    """
    A pandas ExtensionArray for NumPy data.

    .. versionadded :: 0.24.0

    This is mostly for internal compatibility, and is not especially
    useful on its own.

    Parameters
    ----------
    values : ndarray
        The NumPy ndarray to wrap. Must be 1-dimensional.
    copy : bool, default False
        Whether to copy `values`.

    Notes
    -----
    Operations like ``+`` and applying ufuncs requires NumPy>=1.13.
    """
    # If you're wondering why pd.Series(cls) doesn't put the array in an
    # ExtensionBlock, search for `ABCPandasArray`. We check for
    # that _typ to ensure that that users don't unnecessarily use EAs inside
    # pandas internals, which turns off things like block consolidation.
    _typ = "npy_extension"
    __array_priority__ = 1000

    # ------------------------------------------------------------------------
    # Constructors

    def __init__(self, values, copy=False):
        if isinstance(values, type(self)):
            values = values._ndarray
        if not isinstance(values, np.ndarray):
            raise ValueError("'values' must be a NumPy array.")

        if values.ndim != 1:
            raise ValueError("PandasArray must be 1-dimensional.")

        if copy:
            values = values.copy()

        self._ndarray = values
        self._dtype = PandasDtype(values.dtype)

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        if isinstance(dtype, PandasDtype):
            dtype = dtype._dtype

        result = np.asarray(scalars, dtype=dtype)
        if copy and result is scalars:
            result = result.copy()
        return cls(result)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls(values)

    @classmethod
    def _concat_same_type(cls, to_concat):
        return cls(np.concatenate(to_concat))

    # ------------------------------------------------------------------------
    # Data

    @property
    def dtype(self):
        return self._dtype

    # ------------------------------------------------------------------------
    # NumPy Array Interface

    def __array__(self, dtype=None):
        return np.asarray(self._ndarray, dtype=dtype)

    _HANDLED_TYPES = (np.ndarray, numbers.Number)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # Lightly modified version of
        # https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/\
        # numpy.lib.mixins.NDArrayOperatorsMixin.html
        # The primary modification is not boxing scalar return values
        # in PandasArray, since pandas' ExtensionArrays are 1-d.
        out = kwargs.get('out', ())
        for x in inputs + out:
            # Only support operations with instances of _HANDLED_TYPES.
            # Use PandasArray instead of type(self) for isinstance to
            # allow subclasses that don't override __array_ufunc__ to
            # handle PandasArray objects.
            if not isinstance(x, self._HANDLED_TYPES + (PandasArray,)):
                return NotImplemented

        # Defer to the implementation of the ufunc on unwrapped values.
        inputs = tuple(x._ndarray if isinstance(x, PandasArray) else x
                       for x in inputs)
        if out:
            kwargs['out'] = tuple(
                x._ndarray if isinstance(x, PandasArray) else x
                for x in out)
        result = getattr(ufunc, method)(*inputs, **kwargs)

        if type(result) is tuple and len(result):
            # multiple return values
            if not lib.is_scalar(result[0]):
                # re-box array-like results
                return tuple(type(self)(x) for x in result)
            else:
                # but not scalar reductions
                return result
        elif method == 'at':
            # no return value
            return None
        else:
            # one return value
            if not lib.is_scalar(result):
                # re-box array-like results, but not scalar reductions
                result = type(self)(result)
            return result

    # ------------------------------------------------------------------------
    # Pandas ExtensionArray Interface

    def __getitem__(self, item):
        if isinstance(item, type(self)):
            item = item._ndarray

        result = self._ndarray[item]
        if not lib.is_scalar(item):
            result = type(self)(result)
        return result

    def __setitem__(self, key, value):
        from pandas.core.internals.arrays import extract_array

        value = extract_array(value, extract_numpy=True)

        if not lib.is_scalar(key) and is_list_like(key):
            key = np.asarray(key)

        if not lib.is_scalar(value):
            value = np.asarray(value)

        values = self._ndarray
        t = np.result_type(value, values)
        if t != self._ndarray.dtype:
            values = values.astype(t, casting='safe')
            values[key] = value
            self._dtype = PandasDtype(t)
            self._ndarray = values
        else:
            self._ndarray[key] = value

    def __len__(self):
        return len(self._ndarray)

    @property
    def nbytes(self):
        return self._ndarray.nbytes

    def isna(self):
        from pandas import isna

        return isna(self._ndarray)

    def fillna(self, value=None, method=None, limit=None):
        # TODO(_values_for_fillna): remove this
        value, method = validate_fillna_kwargs(value, method)

        mask = self.isna()

        if is_array_like(value):
            if len(value) != len(self):
                raise ValueError("Length of 'value' does not match. Got ({}) "
                                 " expected {}".format(len(value), len(self)))
            value = value[mask]

        if mask.any():
            if method is not None:
                func = pad_1d if method == 'pad' else backfill_1d
                new_values = func(self._ndarray, limit=limit,
                                  mask=mask)
                new_values = self._from_sequence(new_values, dtype=self.dtype)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            new_values = self.copy()
        return new_values

    def take(self, indices, allow_fill=False, fill_value=None):
        from pandas.core.algorithms import take

        result = take(self._ndarray, indices, allow_fill=allow_fill,
                      fill_value=fill_value)
        return type(self)(result)

    def copy(self, deep=False):
        return type(self)(self._ndarray.copy())

    def _values_for_argsort(self):
        return self._ndarray

    def _values_for_factorize(self):
        return self._ndarray, -1

    def unique(self):
        from pandas import unique

        return type(self)(unique(self._ndarray))

    # ------------------------------------------------------------------------
    # Reductions

    def _reduce(self, name, skipna=True, **kwargs):
        meth = getattr(self, name, None)
        if meth:
            return meth(skipna=skipna, **kwargs)
        else:
            msg = (
                "'{}' does not implement reduction '{}'"
            )
            raise TypeError(msg.format(type(self).__name__, name))

    def any(self, axis=None, out=None, keepdims=False, skipna=True):
        nv.validate_any((), dict(out=out, keepdims=keepdims))
        return nanops.nanany(self._ndarray, axis=axis, skipna=skipna)

    def all(self, axis=None, out=None, keepdims=False, skipna=True):
        nv.validate_all((), dict(out=out, keepdims=keepdims))
        return nanops.nanall(self._ndarray, axis=axis, skipna=skipna)

    def min(self, axis=None, out=None, keepdims=False, skipna=True):
        nv.validate_min((), dict(out=out, keepdims=keepdims))
        return nanops.nanmin(self._ndarray, axis=axis, skipna=skipna)

    def max(self, axis=None, out=None, keepdims=False, skipna=True):
        nv.validate_max((), dict(out=out, keepdims=keepdims))
        return nanops.nanmax(self._ndarray, axis=axis, skipna=skipna)

    def sum(self, axis=None, dtype=None, out=None, keepdims=False,
            initial=None, skipna=True, min_count=0):
        nv.validate_sum((), dict(dtype=dtype, out=out, keepdims=keepdims,
                                 initial=initial))
        return nanops.nansum(self._ndarray, axis=axis, skipna=skipna,
                             min_count=min_count)

    def prod(self, axis=None, dtype=None, out=None, keepdims=False,
             initial=None, skipna=True, min_count=0):
        nv.validate_prod((), dict(dtype=dtype, out=out, keepdims=keepdims,
                                  initial=initial))
        return nanops.nanprod(self._ndarray, axis=axis, skipna=skipna,
                              min_count=min_count)

    def mean(self, axis=None, dtype=None, out=None, keepdims=False,
             skipna=True):
        nv.validate_mean((), dict(dtype=dtype, out=out, keepdims=keepdims))
        return nanops.nanmean(self._ndarray, axis=axis, skipna=skipna)

    def median(self, axis=None, out=None, overwrite_input=False,
               keepdims=False, skipna=True):
        nv.validate_median((), dict(out=out, overwrite_input=overwrite_input,
                                    keepdims=keepdims))
        return nanops.nanmedian(self._ndarray, axis=axis, skipna=skipna)

    def std(self, axis=None, dtype=None, out=None, ddof=1, keepdims=False,
            skipna=True):
        nv.validate_stat_ddof_func((), dict(dtype=dtype, out=out,
                                            keepdims=keepdims),
                                   fname='std')
        return nanops.nanstd(self._ndarray, axis=axis, skipna=skipna,
                             ddof=ddof)

    def var(self, axis=None, dtype=None, out=None, ddof=1, keepdims=False,
            skipna=True):
        nv.validate_stat_ddof_func((), dict(dtype=dtype, out=out,
                                            keepdims=keepdims),
                                   fname='var')
        return nanops.nanvar(self._ndarray, axis=axis, skipna=skipna,
                             ddof=ddof)

    def sem(self, axis=None, dtype=None, out=None, ddof=1, keepdims=False,
            skipna=True):
        nv.validate_stat_ddof_func((), dict(dtype=dtype, out=out,
                                            keepdims=keepdims),
                                   fname='sem')
        return nanops.nansem(self._ndarray, axis=axis, skipna=skipna,
                             ddof=ddof)

    def kurt(self, axis=None, dtype=None, out=None, keepdims=False,
             skipna=True):
        nv.validate_stat_ddof_func((), dict(dtype=dtype, out=out,
                                            keepdims=keepdims),
                                   fname='kurt')
        return nanops.nankurt(self._ndarray, axis=axis, skipna=skipna)

    def skew(self, axis=None, dtype=None, out=None, keepdims=False,
             skipna=True):
        nv.validate_stat_ddof_func((), dict(dtype=dtype, out=out,
                                            keepdims=keepdims),
                                   fname='skew')
        return nanops.nanskew(self._ndarray, axis=axis, skipna=skipna)

    # ------------------------------------------------------------------------
    # Additional Methods
    def to_numpy(self, dtype=None, copy=False):
        """
        Convert the PandasArray to a :class:`numpy.ndarray`.

        By default, this requires no coercion or copying of data.

        Parameters
        ----------
        dtype : numpy.dtype
            The NumPy dtype to pass to :func:`numpy.asarray`.
        copy : bool, default False
            Whether to copy the underlying data.

        Returns
        -------
        ndarray
        """
        result = np.asarray(self._ndarray, dtype=dtype)
        if copy and result is self._ndarray:
            result = result.copy()

        return result

    # ------------------------------------------------------------------------
    # Ops

    def __invert__(self):
        return type(self)(~self._ndarray)

    @classmethod
    def _create_arithmetic_method(cls, op):
        def arithmetic_method(self, other):
            if isinstance(other, (ABCIndexClass, ABCSeries)):
                return NotImplemented

            elif isinstance(other, cls):
                other = other._ndarray

            with np.errstate(all="ignore"):
                result = op(self._ndarray, other)

            if op is divmod:
                a, b = result
                return cls(a), cls(b)

            return cls(result)

        return compat.set_function_name(arithmetic_method,
                                        "__{}__".format(op.__name__),
                                        cls)

    _create_comparison_method = _create_arithmetic_method


PandasArray._add_arithmetic_ops()
PandasArray._add_comparison_ops()
