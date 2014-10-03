"""
numpy.ma : a package to handle missing or invalid values.

This package was initially written for numarray by Paul F. Dubois
at Lawrence Livermore National Laboratory.
In 2006, the package was completely rewritten by Pierre Gerard-Marchant
(University of Georgia) to make the MaskedArray class a subclass of ndarray,
and to improve support of structured arrays.


Copyright 1999, 2000, 2001 Regents of the University of California.
Released for unlimited redistribution.

* Adapted for numpy_core 2005 by Travis Oliphant and (mainly) Paul Dubois.
* Subclassing of the base `ndarray` 2006 by Pierre Gerard-Marchant
  (pgmdevlist_AT_gmail_DOT_com)
* Improvements suggested by Reggie Dugard (reggie_AT_merfinllc_DOT_com)

.. moduleauthor:: Pierre Gerard-Marchant

"""
# pylint: disable-msg=E1002
from __future__ import division, absolute_import, print_function

import sys
import warnings
from functools import reduce

import numpy as np
import numpy.core.umath as umath
import numpy.core.numerictypes as ntypes
from numpy import ndarray, amax, amin, iscomplexobj, bool_
from numpy import array as narray
from numpy.lib.function_base import angle
from numpy.compat import getargspec, formatargspec, long, basestring
from numpy import expand_dims as n_expand_dims

if sys.version_info[0] >= 3:
    import pickle
else:
    import cPickle as pickle

__author__ = "Pierre GF Gerard-Marchant"
__docformat__ = "restructuredtext en"

__all__ = ['MAError', 'MaskError', 'MaskType', 'MaskedArray',
           'bool_',
           'abs', 'absolute', 'add', 'all', 'allclose', 'allequal', 'alltrue',
           'amax', 'amin', 'angle', 'anom', 'anomalies', 'any', 'append', 'arange',
           'arccos', 'arccosh', 'arcsin', 'arcsinh', 'arctan', 'arctan2',
           'arctanh', 'argmax', 'argmin', 'argsort', 'around',
           'array', 'asarray', 'asanyarray',
           'bitwise_and', 'bitwise_or', 'bitwise_xor',
           'ceil', 'choose', 'clip', 'common_fill_value', 'compress',
           'compressed', 'concatenate', 'conjugate', 'copy', 'cos', 'cosh',
           'count', 'cumprod', 'cumsum',
           'default_fill_value', 'diag', 'diagonal', 'diff', 'divide', 'dump',
           'dumps',
           'empty', 'empty_like', 'equal', 'exp', 'expand_dims',
           'fabs', 'flatten_mask', 'fmod', 'filled', 'floor', 'floor_divide',
           'fix_invalid', 'flatten_structured_array', 'frombuffer', 'fromflex',
           'fromfunction',
           'getdata', 'getmask', 'getmaskarray', 'greater', 'greater_equal',
           'harden_mask', 'hypot',
           'identity', 'ids', 'indices', 'inner', 'innerproduct',
           'isMA', 'isMaskedArray', 'is_mask', 'is_masked', 'isarray',
           'left_shift', 'less', 'less_equal', 'load', 'loads', 'log', 'log2',
           'log10', 'logical_and', 'logical_not', 'logical_or', 'logical_xor',
           'make_mask', 'make_mask_descr', 'make_mask_none', 'mask_or',
           'masked', 'masked_array', 'masked_equal', 'masked_greater',
           'masked_greater_equal', 'masked_inside', 'masked_invalid',
           'masked_less', 'masked_less_equal', 'masked_not_equal',
           'masked_object', 'masked_outside', 'masked_print_option',
           'masked_singleton', 'masked_values', 'masked_where', 'max', 'maximum',
           'maximum_fill_value', 'mean', 'min', 'minimum', 'minimum_fill_value',
           'mod', 'multiply', 'mvoid',
           'negative', 'nomask', 'nonzero', 'not_equal',
           'ones', 'outer', 'outerproduct',
           'power', 'prod', 'product', 'ptp', 'put', 'putmask',
           'rank', 'ravel', 'remainder', 'repeat', 'reshape', 'resize',
           'right_shift', 'round_', 'round',
           'set_fill_value', 'shape', 'sin', 'sinh', 'size', 'sometrue',
           'sort', 'soften_mask', 'sqrt', 'squeeze', 'std', 'subtract', 'sum',
           'swapaxes',
           'take', 'tan', 'tanh', 'trace', 'transpose', 'true_divide',
           'var', 'where',
           'zeros']

MaskType = np.bool_
nomask = MaskType(0)

def doc_note(initialdoc, note):
    """
    Adds a Notes section to an existing docstring.
    """
    if initialdoc is None:
        return
    if note is None:
        return initialdoc
    newdoc = """
    %s

    Notes
    -----
    %s
    """
    return newdoc % (initialdoc, note)

def get_object_signature(obj):
    """
    Get the signature from obj
    """
    try:
        sig = formatargspec(*getargspec(obj))
    except TypeError as errmsg:
        sig = ''
#        msg = "Unable to retrieve the signature of %s '%s'\n"\
#              "(Initial error message: %s)"
#        warnings.warn(msg % (type(obj),
#                             getattr(obj, '__name__', '???'),
#                             errmsg))
    return sig


#####--------------------------------------------------------------------------
#---- --- Exceptions ---
#####--------------------------------------------------------------------------
class MAError(Exception):
    """Class for masked array related errors."""
    pass
class MaskError(MAError):
    "Class for mask related errors."
    pass


#####--------------------------------------------------------------------------
#---- --- Filling options ---
#####--------------------------------------------------------------------------
# b: boolean - c: complex - f: floats - i: integer - O: object - S: string
default_filler = {'b': True,
                  'c' : 1.e20 + 0.0j,
                  'f' : 1.e20,
                  'i' : 999999,
                  'O' : '?',
                  'S' : 'N/A',
                  'u' : 999999,
                  'V' : '???',
                  'U' : 'N/A',
                  'M8[D]' : np.datetime64('NaT', 'D'),
                  'M8[us]' : np.datetime64('NaT', 'us')
                  }
max_filler = ntypes._minvals
max_filler.update([(k, -np.inf) for k in [np.float32, np.float64]])
min_filler = ntypes._maxvals
min_filler.update([(k, +np.inf) for k in [np.float32, np.float64]])
if 'float128' in ntypes.typeDict:
    max_filler.update([(np.float128, -np.inf)])
    min_filler.update([(np.float128, +np.inf)])


def default_fill_value(obj):
    """
    Return the default fill value for the argument object.

    The default filling value depends on the datatype of the input
    array or the type of the input scalar:

       ========  ========
       datatype  default
       ========  ========
       bool      True
       int       999999
       float     1.e20
       complex   1.e20+0j
       object    '?'
       string    'N/A'
       ========  ========


    Parameters
    ----------
    obj : ndarray, dtype or scalar
        The array data-type or scalar for which the default fill value
        is returned.

    Returns
    -------
    fill_value : scalar
        The default fill value.

    Examples
    --------
    >>> np.ma.default_fill_value(1)
    999999
    >>> np.ma.default_fill_value(np.array([1.1, 2., np.pi]))
    1e+20
    >>>  np.ma.default_fill_value(np.dtype(complex))
    (1e+20+0j)

    """
    if hasattr(obj, 'dtype'):
        defval = _check_fill_value(None, obj.dtype)
    elif isinstance(obj, np.dtype):
        if obj.subdtype:
            defval = default_filler.get(obj.subdtype[0].kind, '?')
        elif obj.kind == 'M':
            defval = default_filler.get(obj.str[1:], '?')
        else:
            defval = default_filler.get(obj.kind, '?')
    elif isinstance(obj, float):
        defval = default_filler['f']
    elif isinstance(obj, int) or isinstance(obj, long):
        defval = default_filler['i']
    elif isinstance(obj, str):
        defval = default_filler['S']
    elif isinstance(obj, unicode):
        defval = default_filler['U']
    elif isinstance(obj, complex):
        defval = default_filler['c']
    else:
        defval = default_filler['O']
    return defval


def _recursive_extremum_fill_value(ndtype, extremum):
    names = ndtype.names
    if names:
        deflist = []
        for name in names:
            fval = _recursive_extremum_fill_value(ndtype[name], extremum)
            deflist.append(fval)
        return tuple(deflist)
    return extremum[ndtype]


def minimum_fill_value(obj):
    """
    Return the maximum value that can be represented by the dtype of an object.

    This function is useful for calculating a fill value suitable for
    taking the minimum of an array with a given dtype.

    Parameters
    ----------
    obj : ndarray or dtype
        An object that can be queried for it's numeric type.

    Returns
    -------
    val : scalar
        The maximum representable value.

    Raises
    ------
    TypeError
        If `obj` isn't a suitable numeric type.

    See Also
    --------
    maximum_fill_value : The inverse function.
    set_fill_value : Set the filling value of a masked array.
    MaskedArray.fill_value : Return current fill value.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.int8()
    >>> ma.minimum_fill_value(a)
    127
    >>> a = np.int32()
    >>> ma.minimum_fill_value(a)
    2147483647

    An array of numeric data can also be passed.

    >>> a = np.array([1, 2, 3], dtype=np.int8)
    >>> ma.minimum_fill_value(a)
    127
    >>> a = np.array([1, 2, 3], dtype=np.float32)
    >>> ma.minimum_fill_value(a)
    inf

    """
    errmsg = "Unsuitable type for calculating minimum."
    if hasattr(obj, 'dtype'):
        return _recursive_extremum_fill_value(obj.dtype, min_filler)
    elif isinstance(obj, float):
        return min_filler[ntypes.typeDict['float_']]
    elif isinstance(obj, int):
        return min_filler[ntypes.typeDict['int_']]
    elif isinstance(obj, long):
        return min_filler[ntypes.typeDict['uint']]
    elif isinstance(obj, np.dtype):
        return min_filler[obj]
    else:
        raise TypeError(errmsg)


def maximum_fill_value(obj):
    """
    Return the minimum value that can be represented by the dtype of an object.

    This function is useful for calculating a fill value suitable for
    taking the maximum of an array with a given dtype.

    Parameters
    ----------
    obj : {ndarray, dtype}
        An object that can be queried for it's numeric type.

    Returns
    -------
    val : scalar
        The minimum representable value.

    Raises
    ------
    TypeError
        If `obj` isn't a suitable numeric type.

    See Also
    --------
    minimum_fill_value : The inverse function.
    set_fill_value : Set the filling value of a masked array.
    MaskedArray.fill_value : Return current fill value.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.int8()
    >>> ma.maximum_fill_value(a)
    -128
    >>> a = np.int32()
    >>> ma.maximum_fill_value(a)
    -2147483648

    An array of numeric data can also be passed.

    >>> a = np.array([1, 2, 3], dtype=np.int8)
    >>> ma.maximum_fill_value(a)
    -128
    >>> a = np.array([1, 2, 3], dtype=np.float32)
    >>> ma.maximum_fill_value(a)
    -inf

    """
    errmsg = "Unsuitable type for calculating maximum."
    if hasattr(obj, 'dtype'):
        return _recursive_extremum_fill_value(obj.dtype, max_filler)
    elif isinstance(obj, float):
        return max_filler[ntypes.typeDict['float_']]
    elif isinstance(obj, int):
        return max_filler[ntypes.typeDict['int_']]
    elif isinstance(obj, long):
        return max_filler[ntypes.typeDict['uint']]
    elif isinstance(obj, np.dtype):
        return max_filler[obj]
    else:
        raise TypeError(errmsg)


def _recursive_set_default_fill_value(dtypedescr):
    deflist = []
    for currentdescr in dtypedescr:
        currenttype = currentdescr[1]
        if isinstance(currenttype, list):
            deflist.append(tuple(_recursive_set_default_fill_value(currenttype)))
        else:
            deflist.append(default_fill_value(np.dtype(currenttype)))
    return tuple(deflist)

def _recursive_set_fill_value(fillvalue, dtypedescr):
    fillvalue = np.resize(fillvalue, len(dtypedescr))
    output_value = []
    for (fval, descr) in zip(fillvalue, dtypedescr):
        cdtype = descr[1]
        if isinstance(cdtype, list):
            output_value.append(tuple(_recursive_set_fill_value(fval, cdtype)))
        else:
            output_value.append(np.array(fval, dtype=cdtype).item())
    return tuple(output_value)


def _check_fill_value(fill_value, ndtype):
    """
    Private function validating the given `fill_value` for the given dtype.

    If fill_value is None, it is set to the default corresponding to the dtype
    if this latter is standard (no fields). If the datatype is flexible (named
    fields), fill_value is set to a tuple whose elements are the default fill
    values corresponding to each field.

    If fill_value is not None, its value is forced to the given dtype.

    """
    ndtype = np.dtype(ndtype)
    fields = ndtype.fields
    if fill_value is None:
        if fields:
            descr = ndtype.descr
            fill_value = np.array(_recursive_set_default_fill_value(descr),
                                  dtype=ndtype,)
        else:
            fill_value = default_fill_value(ndtype)
    elif fields:
        fdtype = [(_[0], _[1]) for _ in ndtype.descr]
        if isinstance(fill_value, (ndarray, np.void)):
            try:
                fill_value = np.array(fill_value, copy=False, dtype=fdtype)
            except ValueError:
                err_msg = "Unable to transform %s to dtype %s"
                raise ValueError(err_msg % (fill_value, fdtype))
        else:
            descr = ndtype.descr
            fill_value = np.asarray(fill_value, dtype=object)
            fill_value = np.array(_recursive_set_fill_value(fill_value, descr),
                                  dtype=ndtype)
    else:
        if isinstance(fill_value, basestring) and (ndtype.char not in 'OSVU'):
            err_msg = "Cannot set fill value of string with array of dtype %s"
            raise TypeError(err_msg % ndtype)
        else:
            # In case we want to convert 1e20 to int...
            try:
                fill_value = np.array(fill_value, copy=False, dtype=ndtype)
            except OverflowError:
                # Raise TypeError instead of OverflowError. OverflowError
                # is seldom used, and the real problem here is that the
                # passed fill_value is not compatible with the ndtype.
                err_msg = "Fill value %s overflows dtype %s"
                raise TypeError(err_msg % (fill_value, ndtype))
    return np.array(fill_value)


def set_fill_value(a, fill_value):
    """
    Set the filling value of a, if a is a masked array.

    This function changes the fill value of the masked array `a` in place.
    If `a` is not a masked array, the function returns silently, without
    doing anything.

    Parameters
    ----------
    a : array_like
        Input array.
    fill_value : dtype
        Filling value. A consistency test is performed to make sure
        the value is compatible with the dtype of `a`.

    Returns
    -------
    None
        Nothing returned by this function.

    See Also
    --------
    maximum_fill_value : Return the default fill value for a dtype.
    MaskedArray.fill_value : Return current fill value.
    MaskedArray.set_fill_value : Equivalent method.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(5)
    >>> a
    array([0, 1, 2, 3, 4])
    >>> a = ma.masked_where(a < 3, a)
    >>> a
    masked_array(data = [-- -- -- 3 4],
          mask = [ True  True  True False False],
          fill_value=999999)
    >>> ma.set_fill_value(a, -999)
    >>> a
    masked_array(data = [-- -- -- 3 4],
          mask = [ True  True  True False False],
          fill_value=-999)

    Nothing happens if `a` is not a masked array.

    >>> a = range(5)
    >>> a
    [0, 1, 2, 3, 4]
    >>> ma.set_fill_value(a, 100)
    >>> a
    [0, 1, 2, 3, 4]
    >>> a = np.arange(5)
    >>> a
    array([0, 1, 2, 3, 4])
    >>> ma.set_fill_value(a, 100)
    >>> a
    array([0, 1, 2, 3, 4])

    """
    if isinstance(a, MaskedArray):
        a.set_fill_value(fill_value)
    return

def get_fill_value(a):
    """
    Return the filling value of a, if any.  Otherwise, returns the
    default filling value for that type.

    """
    if isinstance(a, MaskedArray):
        result = a.fill_value
    else:
        result = default_fill_value(a)
    return result

def common_fill_value(a, b):
    """
    Return the common filling value of two masked arrays, if any.

    If ``a.fill_value == b.fill_value``, return the fill value,
    otherwise return None.

    Parameters
    ----------
    a, b : MaskedArray
        The masked arrays for which to compare fill values.

    Returns
    -------
    fill_value : scalar or None
        The common fill value, or None.

    Examples
    --------
    >>> x = np.ma.array([0, 1.], fill_value=3)
    >>> y = np.ma.array([0, 1.], fill_value=3)
    >>> np.ma.common_fill_value(x, y)
    3.0

    """
    t1 = get_fill_value(a)
    t2 = get_fill_value(b)
    if t1 == t2:
        return t1
    return None


#####--------------------------------------------------------------------------
def filled(a, fill_value=None):
    """
    Return input as an array with masked data replaced by a fill value.

    If `a` is not a `MaskedArray`, `a` itself is returned.
    If `a` is a `MaskedArray` and `fill_value` is None, `fill_value` is set to
    ``a.fill_value``.

    Parameters
    ----------
    a : MaskedArray or array_like
        An input object.
    fill_value : scalar, optional
        Filling value. Default is None.

    Returns
    -------
    a : ndarray
        The filled array.

    See Also
    --------
    compressed

    Examples
    --------
    >>> x = np.ma.array(np.arange(9).reshape(3, 3), mask=[[1, 0, 0],
    ...                                                   [1, 0, 0],
    ...                                                   [0, 0, 0]])
    >>> x.filled()
    array([[999999,      1,      2],
           [999999,      4,      5],
           [     6,      7,      8]])

    """
    if hasattr(a, 'filled'):
        return a.filled(fill_value)
    elif isinstance(a, ndarray):
        # Should we check for contiguity ? and a.flags['CONTIGUOUS']:
        return a
    elif isinstance(a, dict):
        return np.array(a, 'O')
    else:
        return np.array(a)

#####--------------------------------------------------------------------------
def get_masked_subclass(*arrays):
    """
    Return the youngest subclass of MaskedArray from a list of (masked) arrays.
    In case of siblings, the first listed takes over.

    """
    if len(arrays) == 1:
        arr = arrays[0]
        if isinstance(arr, MaskedArray):
            rcls = type(arr)
        else:
            rcls = MaskedArray
    else:
        arrcls = [type(a) for a in arrays]
        rcls = arrcls[0]
        if not issubclass(rcls, MaskedArray):
            rcls = MaskedArray
        for cls in arrcls[1:]:
            if issubclass(cls, rcls):
                rcls = cls
    # Don't return MaskedConstant as result: revert to MaskedArray
    if rcls.__name__ == 'MaskedConstant':
        return MaskedArray
    return rcls

#####--------------------------------------------------------------------------
def getdata(a, subok=True):
    """
    Return the data of a masked array as an ndarray.

    Return the data of `a` (if any) as an ndarray if `a` is a ``MaskedArray``,
    else return `a` as a ndarray or subclass (depending on `subok`) if not.

    Parameters
    ----------
    a : array_like
        Input ``MaskedArray``, alternatively a ndarray or a subclass thereof.
    subok : bool
        Whether to force the output to be a `pure` ndarray (False) or to
        return a subclass of ndarray if appropriate (True, default).

    See Also
    --------
    getmask : Return the mask of a masked array, or nomask.
    getmaskarray : Return the mask of a masked array, or full array of False.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = ma.masked_equal([[1,2],[3,4]], 2)
    >>> a
    masked_array(data =
     [[1 --]
     [3 4]],
          mask =
     [[False  True]
     [False False]],
          fill_value=999999)
    >>> ma.getdata(a)
    array([[1, 2],
           [3, 4]])

    Equivalently use the ``MaskedArray`` `data` attribute.

    >>> a.data
    array([[1, 2],
           [3, 4]])

    """
    try:
        data = a._data
    except AttributeError:
        data = np.array(a, copy=False, subok=subok)
    if not subok:
        return data.view(ndarray)
    return data
get_data = getdata


def fix_invalid(a, mask=nomask, copy=True, fill_value=None):
    """
    Return input with invalid data masked and replaced by a fill value.

    Invalid data means values of `nan`, `inf`, etc.

    Parameters
    ----------
    a : array_like
        Input array, a (subclass of) ndarray.
    copy : bool, optional
        Whether to use a copy of `a` (True) or to fix `a` in place (False).
        Default is True.
    fill_value : scalar, optional
        Value used for fixing invalid data. Default is None, in which case
        the ``a.fill_value`` is used.

    Returns
    -------
    b : MaskedArray
        The input array with invalid entries fixed.

    Notes
    -----
    A copy is performed by default.

    Examples
    --------
    >>> x = np.ma.array([1., -1, np.nan, np.inf], mask=[1] + [0]*3)
    >>> x
    masked_array(data = [-- -1.0 nan inf],
                 mask = [ True False False False],
           fill_value = 1e+20)
    >>> np.ma.fix_invalid(x)
    masked_array(data = [-- -1.0 -- --],
                 mask = [ True False  True  True],
           fill_value = 1e+20)

    >>> fixed = np.ma.fix_invalid(x)
    >>> fixed.data
    array([  1.00000000e+00,  -1.00000000e+00,   1.00000000e+20,
             1.00000000e+20])
    >>> x.data
    array([  1.,  -1.,  NaN,  Inf])

    """
    a = masked_array(a, copy=copy, mask=mask, subok=True)
    #invalid = (numpy.isnan(a._data) | numpy.isinf(a._data))
    invalid = np.logical_not(np.isfinite(a._data))
    if not invalid.any():
        return a
    a._mask |= invalid
    if fill_value is None:
        fill_value = a.fill_value
    a._data[invalid] = fill_value
    return a



#####--------------------------------------------------------------------------
#---- --- Ufuncs ---
#####--------------------------------------------------------------------------
ufunc_domain = {}
ufunc_fills = {}

class _DomainCheckInterval:
    """
    Define a valid interval, so that :

    ``domain_check_interval(a,b)(x) == True`` where
    ``x < a`` or ``x > b``.

    """
    def __init__(self, a, b):
        "domain_check_interval(a,b)(x) = true where x < a or y > b"
        if (a > b):
            (a, b) = (b, a)
        self.a = a
        self.b = b

    def __call__ (self, x):
        "Execute the call behavior."
        return umath.logical_or(umath.greater (x, self.b),
                                umath.less(x, self.a))



class _DomainTan:
    """Define a valid interval for the `tan` function, so that:

    ``domain_tan(eps) = True`` where ``abs(cos(x)) < eps``

    """
    def __init__(self, eps):
        "domain_tan(eps) = true where abs(cos(x)) < eps)"
        self.eps = eps

    def __call__ (self, x):
        "Executes the call behavior."
        return umath.less(umath.absolute(umath.cos(x)), self.eps)



class _DomainSafeDivide:
    """Define a domain for safe division."""
    def __init__ (self, tolerance=None):
        self.tolerance = tolerance

    def __call__ (self, a, b):
        # Delay the selection of the tolerance to here in order to reduce numpy
        # import times. The calculation of these parameters is a substantial
        # component of numpy's import time.
        if self.tolerance is None:
            self.tolerance = np.finfo(float).tiny
        return umath.absolute(a) * self.tolerance >= umath.absolute(b)



class _DomainGreater:
    """DomainGreater(v)(x) is True where x <= v."""
    def __init__(self, critical_value):
        "DomainGreater(v)(x) = true where x <= v"
        self.critical_value = critical_value

    def __call__ (self, x):
        "Executes the call behavior."
        return umath.less_equal(x, self.critical_value)



class _DomainGreaterEqual:
    """DomainGreaterEqual(v)(x) is True where x < v."""
    def __init__(self, critical_value):
        "DomainGreaterEqual(v)(x) = true where x < v"
        self.critical_value = critical_value

    def __call__ (self, x):
        "Executes the call behavior."
        return umath.less(x, self.critical_value)

#..............................................................................
class _MaskedUnaryOperation:
    """
    Defines masked version of unary operations, where invalid values are
    pre-masked.

    Parameters
    ----------
    mufunc : callable
        The function for which to define a masked version. Made available
        as ``_MaskedUnaryOperation.f``.
    fill : scalar, optional
        Filling value, default is 0.
    domain : class instance
        Domain for the function. Should be one of the ``_Domain*``
        classes. Default is None.

    """
    def __init__ (self, mufunc, fill=0, domain=None):
        """ _MaskedUnaryOperation(aufunc, fill=0, domain=None)
            aufunc(fill) must be defined
            self(x) returns aufunc(x)
            with masked values where domain(x) is true or getmask(x) is true.
        """
        self.f = mufunc
        self.fill = fill
        self.domain = domain
        self.__doc__ = getattr(mufunc, "__doc__", str(mufunc))
        self.__name__ = getattr(mufunc, "__name__", str(mufunc))
        ufunc_domain[mufunc] = domain
        ufunc_fills[mufunc] = fill
    #
    def __call__ (self, a, *args, **kwargs):
        "Execute the call behavior."
        d = getdata(a)
        # Case 1.1. : Domained function
        if self.domain is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                result = self.f(d, *args, **kwargs)
            # Make a mask
            m = ~umath.isfinite(result)
            m |= self.domain(d)
            m |= getmask(a)
        # Case 1.2. : Function without a domain
        else:
            # Get the result and the mask
            result = self.f(d, *args, **kwargs)
            m = getmask(a)
        # Case 2.1. : The result is scalarscalar
        if not result.ndim:
            if m:
                return masked
            return result
        # Case 2.2. The result is an array
        # We need to fill the invalid data back w/ the input
        # Now, that's plain silly: in C, we would just skip the element and keep
        # the original, but we do have to do it that way in Python
        if m is not nomask:
            # In case result has a lower dtype than the inputs (as in equal)
            try:
                np.copyto(result, d, where=m)
            except TypeError:
                pass
        # Transform to
        if isinstance(a, MaskedArray):
            subtype = type(a)
        else:
            subtype = MaskedArray
        result = result.view(subtype)
        result._mask = m
        result._update_from(a)
        return result
    #
    def __str__ (self):
        return "Masked version of %s. [Invalid values are masked]" % str(self.f)



class _MaskedBinaryOperation:
    """
    Define masked version of binary operations, where invalid
    values are pre-masked.

    Parameters
    ----------
    mbfunc : function
        The function for which to define a masked version. Made available
        as ``_MaskedBinaryOperation.f``.
    domain : class instance
        Default domain for the function. Should be one of the ``_Domain*``
        classes. Default is None.
    fillx : scalar, optional
        Filling value for the first argument, default is 0.
    filly : scalar, optional
        Filling value for the second argument, default is 0.

    """
    def __init__ (self, mbfunc, fillx=0, filly=0):
        """abfunc(fillx, filly) must be defined.
           abfunc(x, filly) = x for all x to enable reduce.
        """
        self.f = mbfunc
        self.fillx = fillx
        self.filly = filly
        self.__doc__ = getattr(mbfunc, "__doc__", str(mbfunc))
        self.__name__ = getattr(mbfunc, "__name__", str(mbfunc))
        ufunc_domain[mbfunc] = None
        ufunc_fills[mbfunc] = (fillx, filly)

    def __call__ (self, a, b, *args, **kwargs):
        "Execute the call behavior."
        # Get the data, as ndarray
        (da, db) = (getdata(a, subok=False), getdata(b, subok=False))
        # Get the mask
        (ma, mb) = (getmask(a), getmask(b))
        if ma is nomask:
            if mb is nomask:
                m = nomask
            else:
                m = umath.logical_or(getmaskarray(a), mb)
        elif mb is nomask:
            m = umath.logical_or(ma, getmaskarray(b))
        else:
            m = umath.logical_or(ma, mb)
        # Get the result
        with np.errstate(divide='ignore', invalid='ignore'):
            result = self.f(da, db, *args, **kwargs)
        # check it worked
        if result is NotImplemented:
            return NotImplemented
        # Case 1. : scalar
        if not result.ndim:
            if m:
                return masked
            return result
        # Case 2. : array
        # Revert result to da where masked
        if m is not nomask:
            np.copyto(result, da, casting='unsafe', where=m)
        # Transforms to a (subclass of) MaskedArray
        result = result.view(get_masked_subclass(a, b))
        result._mask = m
        # Update the optional info from the inputs
        if isinstance(b, MaskedArray):
            if isinstance(a, MaskedArray):
                result._update_from(a)
            else:
                result._update_from(b)
        elif isinstance(a, MaskedArray):
            result._update_from(a)
        return result


    def reduce(self, target, axis=0, dtype=None):
        """Reduce `target` along the given `axis`."""
        if isinstance(target, MaskedArray):
            tclass = type(target)
        else:
            tclass = MaskedArray
        m = getmask(target)
        t = filled(target, self.filly)
        if t.shape == ():
            t = t.reshape(1)
            if m is not nomask:
                m = make_mask(m, copy=1)
                m.shape = (1,)
        if m is nomask:
            return self.f.reduce(t, axis).view(tclass)
        t = t.view(tclass)
        t._mask = m
        tr = self.f.reduce(getdata(t), axis, dtype=dtype or t.dtype)
        mr = umath.logical_and.reduce(m, axis)
        tr = tr.view(tclass)
        if mr.ndim > 0:
            tr._mask = mr
            return tr
        elif mr:
            return masked
        return tr

    def outer (self, a, b):
        """Return the function applied to the outer product of a and b.

        """
        ma = getmask(a)
        mb = getmask(b)
        if ma is nomask and mb is nomask:
            m = nomask
        else:
            ma = getmaskarray(a)
            mb = getmaskarray(b)
            m = umath.logical_or.outer(ma, mb)
        if (not m.ndim) and m:
            return masked
        (da, db) = (getdata(a), getdata(b))
        d = self.f.outer(da, db)
        # check it worked
        if d is NotImplemented:
            return NotImplemented
        if m is not nomask:
            np.copyto(d, da, where=m)
        if d.shape:
            d = d.view(get_masked_subclass(a, b))
            d._mask = m
        return d

    def accumulate (self, target, axis=0):
        """Accumulate `target` along `axis` after filling with y fill
        value.

        """
        if isinstance(target, MaskedArray):
            tclass = type(target)
        else:
            tclass = MaskedArray
        t = filled(target, self.filly)
        return self.f.accumulate(t, axis).view(tclass)

    def __str__ (self):
        return "Masked version of " + str(self.f)



class _DomainedBinaryOperation:
    """
    Define binary operations that have a domain, like divide.

    They have no reduce, outer or accumulate.

    Parameters
    ----------
    mbfunc : function
        The function for which to define a masked version. Made available
        as ``_DomainedBinaryOperation.f``.
    domain : class instance
        Default domain for the function. Should be one of the ``_Domain*``
        classes.
    fillx : scalar, optional
        Filling value for the first argument, default is 0.
    filly : scalar, optional
        Filling value for the second argument, default is 0.

    """
    def __init__ (self, dbfunc, domain, fillx=0, filly=0):
        """abfunc(fillx, filly) must be defined.
           abfunc(x, filly) = x for all x to enable reduce.
        """
        self.f = dbfunc
        self.domain = domain
        self.fillx = fillx
        self.filly = filly
        self.__doc__ = getattr(dbfunc, "__doc__", str(dbfunc))
        self.__name__ = getattr(dbfunc, "__name__", str(dbfunc))
        ufunc_domain[dbfunc] = domain
        ufunc_fills[dbfunc] = (fillx, filly)

    def __call__(self, a, b, *args, **kwargs):
        "Execute the call behavior."
        # Get the data and the mask
        (da, db) = (getdata(a, subok=False), getdata(b, subok=False))
        (ma, mb) = (getmask(a), getmask(b))
        # Get the result
        with np.errstate(divide='ignore', invalid='ignore'):
            result = self.f(da, db, *args, **kwargs)
        # check it worked
        if result is NotImplemented:
            return NotImplemented
        # Get the mask as a combination of ma, mb and invalid
        m = ~umath.isfinite(result)
        m |= ma
        m |= mb
        # Apply the domain
        domain = ufunc_domain.get(self.f, None)
        if domain is not None:
            m |= filled(domain(da, db), True)
        # Take care of the scalar case first
        if (not m.ndim):
            if m:
                return masked
            else:
                return result
        # When the mask is True, put back da
        np.copyto(result, da, casting='unsafe', where=m)
        result = result.view(get_masked_subclass(a, b))
        result._mask = m
        if isinstance(b, MaskedArray):
            if isinstance(a, MaskedArray):
                result._update_from(a)
            else:
                result._update_from(b)
        elif isinstance(a, MaskedArray):
            result._update_from(a)
        return result

    def __str__ (self):
        return "Masked version of " + str(self.f)

#..............................................................................
# Unary ufuncs
exp = _MaskedUnaryOperation(umath.exp)
conjugate = _MaskedUnaryOperation(umath.conjugate)
sin = _MaskedUnaryOperation(umath.sin)
cos = _MaskedUnaryOperation(umath.cos)
tan = _MaskedUnaryOperation(umath.tan)
arctan = _MaskedUnaryOperation(umath.arctan)
arcsinh = _MaskedUnaryOperation(umath.arcsinh)
sinh = _MaskedUnaryOperation(umath.sinh)
cosh = _MaskedUnaryOperation(umath.cosh)
tanh = _MaskedUnaryOperation(umath.tanh)
abs = absolute = _MaskedUnaryOperation(umath.absolute)
angle = _MaskedUnaryOperation(angle) # from numpy.lib.function_base
fabs = _MaskedUnaryOperation(umath.fabs)
negative = _MaskedUnaryOperation(umath.negative)
floor = _MaskedUnaryOperation(umath.floor)
ceil = _MaskedUnaryOperation(umath.ceil)
around = _MaskedUnaryOperation(np.round_)
logical_not = _MaskedUnaryOperation(umath.logical_not)
# Domained unary ufuncs .......................................................
sqrt = _MaskedUnaryOperation(umath.sqrt, 0.0,
                             _DomainGreaterEqual(0.0))
log = _MaskedUnaryOperation(umath.log, 1.0,
                            _DomainGreater(0.0))
log2 = _MaskedUnaryOperation(umath.log2, 1.0,
                            _DomainGreater(0.0))
log10 = _MaskedUnaryOperation(umath.log10, 1.0,
                              _DomainGreater(0.0))
tan = _MaskedUnaryOperation(umath.tan, 0.0,
                            _DomainTan(1e-35))
arcsin = _MaskedUnaryOperation(umath.arcsin, 0.0,
                               _DomainCheckInterval(-1.0, 1.0))
arccos = _MaskedUnaryOperation(umath.arccos, 0.0,
                               _DomainCheckInterval(-1.0, 1.0))
arccosh = _MaskedUnaryOperation(umath.arccosh, 1.0,
                                _DomainGreaterEqual(1.0))
arctanh = _MaskedUnaryOperation(umath.arctanh, 0.0,
                                _DomainCheckInterval(-1.0 + 1e-15, 1.0 - 1e-15))
# Binary ufuncs ...............................................................
add = _MaskedBinaryOperation(umath.add)
subtract = _MaskedBinaryOperation(umath.subtract)
multiply = _MaskedBinaryOperation(umath.multiply, 1, 1)
arctan2 = _MaskedBinaryOperation(umath.arctan2, 0.0, 1.0)
equal = _MaskedBinaryOperation(umath.equal)
equal.reduce = None
not_equal = _MaskedBinaryOperation(umath.not_equal)
not_equal.reduce = None
less_equal = _MaskedBinaryOperation(umath.less_equal)
less_equal.reduce = None
greater_equal = _MaskedBinaryOperation(umath.greater_equal)
greater_equal.reduce = None
less = _MaskedBinaryOperation(umath.less)
less.reduce = None
greater = _MaskedBinaryOperation(umath.greater)
greater.reduce = None
logical_and = _MaskedBinaryOperation(umath.logical_and)
alltrue = _MaskedBinaryOperation(umath.logical_and, 1, 1).reduce
logical_or = _MaskedBinaryOperation(umath.logical_or)
sometrue = logical_or.reduce
logical_xor = _MaskedBinaryOperation(umath.logical_xor)
bitwise_and = _MaskedBinaryOperation(umath.bitwise_and)
bitwise_or = _MaskedBinaryOperation(umath.bitwise_or)
bitwise_xor = _MaskedBinaryOperation(umath.bitwise_xor)
hypot = _MaskedBinaryOperation(umath.hypot)
# Domained binary ufuncs ......................................................
divide = _DomainedBinaryOperation(umath.divide, _DomainSafeDivide(), 0, 1)
true_divide = _DomainedBinaryOperation(umath.true_divide,
                                        _DomainSafeDivide(), 0, 1)
floor_divide = _DomainedBinaryOperation(umath.floor_divide,
                                         _DomainSafeDivide(), 0, 1)
remainder = _DomainedBinaryOperation(umath.remainder,
                                      _DomainSafeDivide(), 0, 1)
fmod = _DomainedBinaryOperation(umath.fmod, _DomainSafeDivide(), 0, 1)
mod = _DomainedBinaryOperation(umath.mod, _DomainSafeDivide(), 0, 1)


#####--------------------------------------------------------------------------
#---- --- Mask creation functions ---
#####--------------------------------------------------------------------------

def _recursive_make_descr(datatype, newtype=bool_):
    "Private function allowing recursion in make_descr."
    # Do we have some name fields ?
    if datatype.names:
        descr = []
        for name in datatype.names:
            field = datatype.fields[name]
            if len(field) == 3:
                # Prepend the title to the name
                name = (field[-1], name)
            descr.append((name, _recursive_make_descr(field[0], newtype)))
        return descr
    # Is this some kind of composite a la (np.float,2)
    elif datatype.subdtype:
        mdescr = list(datatype.subdtype)
        mdescr[0] = newtype
        return tuple(mdescr)
    else:
        return newtype

def make_mask_descr(ndtype):
    """
    Construct a dtype description list from a given dtype.

    Returns a new dtype object, with the type of all fields in `ndtype` to a
    boolean type. Field names are not altered.

    Parameters
    ----------
    ndtype : dtype
        The dtype to convert.

    Returns
    -------
    result : dtype
        A dtype that looks like `ndtype`, the type of all fields is boolean.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> dtype = np.dtype({'names':['foo', 'bar'],
                          'formats':[np.float32, np.int]})
    >>> dtype
    dtype([('foo', '<f4'), ('bar', '<i4')])
    >>> ma.make_mask_descr(dtype)
    dtype([('foo', '|b1'), ('bar', '|b1')])
    >>> ma.make_mask_descr(np.float32)
    <type 'numpy.bool_'>

    """
    # Make sure we do have a dtype
    if not isinstance(ndtype, np.dtype):
        ndtype = np.dtype(ndtype)
    return np.dtype(_recursive_make_descr(ndtype, np.bool))

def getmask(a):
    """
    Return the mask of a masked array, or nomask.

    Return the mask of `a` as an ndarray if `a` is a `MaskedArray` and the
    mask is not `nomask`, else return `nomask`. To guarantee a full array
    of booleans of the same shape as a, use `getmaskarray`.

    Parameters
    ----------
    a : array_like
        Input `MaskedArray` for which the mask is required.

    See Also
    --------
    getdata : Return the data of a masked array as an ndarray.
    getmaskarray : Return the mask of a masked array, or full array of False.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = ma.masked_equal([[1,2],[3,4]], 2)
    >>> a
    masked_array(data =
     [[1 --]
     [3 4]],
          mask =
     [[False  True]
     [False False]],
          fill_value=999999)
    >>> ma.getmask(a)
    array([[False,  True],
           [False, False]], dtype=bool)

    Equivalently use the `MaskedArray` `mask` attribute.

    >>> a.mask
    array([[False,  True],
           [False, False]], dtype=bool)

    Result when mask == `nomask`

    >>> b = ma.masked_array([[1,2],[3,4]])
    >>> b
    masked_array(data =
     [[1 2]
     [3 4]],
          mask =
     False,
          fill_value=999999)
    >>> ma.nomask
    False
    >>> ma.getmask(b) == ma.nomask
    True
    >>> b.mask == ma.nomask
    True

    """
    return getattr(a, '_mask', nomask)
get_mask = getmask

def getmaskarray(arr):
    """
    Return the mask of a masked array, or full boolean array of False.

    Return the mask of `arr` as an ndarray if `arr` is a `MaskedArray` and
    the mask is not `nomask`, else return a full boolean array of False of
    the same shape as `arr`.

    Parameters
    ----------
    arr : array_like
        Input `MaskedArray` for which the mask is required.

    See Also
    --------
    getmask : Return the mask of a masked array, or nomask.
    getdata : Return the data of a masked array as an ndarray.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = ma.masked_equal([[1,2],[3,4]], 2)
    >>> a
    masked_array(data =
     [[1 --]
     [3 4]],
          mask =
     [[False  True]
     [False False]],
          fill_value=999999)
    >>> ma.getmaskarray(a)
    array([[False,  True],
           [False, False]], dtype=bool)

    Result when mask == ``nomask``

    >>> b = ma.masked_array([[1,2],[3,4]])
    >>> b
    masked_array(data =
     [[1 2]
     [3 4]],
          mask =
     False,
          fill_value=999999)
    >>> >ma.getmaskarray(b)
    array([[False, False],
           [False, False]], dtype=bool)

    """
    mask = getmask(arr)
    if mask is nomask:
        mask = make_mask_none(np.shape(arr), getdata(arr).dtype)
    return mask

def is_mask(m):
    """
    Return True if m is a valid, standard mask.

    This function does not check the contents of the input, only that the
    type is MaskType. In particular, this function returns False if the
    mask has a flexible dtype.

    Parameters
    ----------
    m : array_like
        Array to test.

    Returns
    -------
    result : bool
        True if `m.dtype.type` is MaskType, False otherwise.

    See Also
    --------
    isMaskedArray : Test whether input is an instance of MaskedArray.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> m = ma.masked_equal([0, 1, 0, 2, 3], 0)
    >>> m
    masked_array(data = [-- 1 -- 2 3],
          mask = [ True False  True False False],
          fill_value=999999)
    >>> ma.is_mask(m)
    False
    >>> ma.is_mask(m.mask)
    True

    Input must be an ndarray (or have similar attributes)
    for it to be considered a valid mask.

    >>> m = [False, True, False]
    >>> ma.is_mask(m)
    False
    >>> m = np.array([False, True, False])
    >>> m
    array([False,  True, False], dtype=bool)
    >>> ma.is_mask(m)
    True

    Arrays with complex dtypes don't return True.

    >>> dtype = np.dtype({'names':['monty', 'pithon'],
                          'formats':[np.bool, np.bool]})
    >>> dtype
    dtype([('monty', '|b1'), ('pithon', '|b1')])
    >>> m = np.array([(True, False), (False, True), (True, False)],
                     dtype=dtype)
    >>> m
    array([(True, False), (False, True), (True, False)],
          dtype=[('monty', '|b1'), ('pithon', '|b1')])
    >>> ma.is_mask(m)
    False

    """
    try:
        return m.dtype.type is MaskType
    except AttributeError:
        return False

def make_mask(m, copy=False, shrink=True, dtype=MaskType):
    """
    Create a boolean mask from an array.

    Return `m` as a boolean mask, creating a copy if necessary or requested.
    The function can accept any sequence that is convertible to integers,
    or ``nomask``.  Does not require that contents must be 0s and 1s, values
    of 0 are interepreted as False, everything else as True.

    Parameters
    ----------
    m : array_like
        Potential mask.
    copy : bool, optional
        Whether to return a copy of `m` (True) or `m` itself (False).
    shrink : bool, optional
        Whether to shrink `m` to ``nomask`` if all its values are False.
    dtype : dtype, optional
        Data-type of the output mask. By default, the output mask has
        a dtype of MaskType (bool). If the dtype is flexible, each field
        has a boolean dtype.

    Returns
    -------
    result : ndarray
        A boolean mask derived from `m`.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> m = [True, False, True, True]
    >>> ma.make_mask(m)
    array([ True, False,  True,  True], dtype=bool)
    >>> m = [1, 0, 1, 1]
    >>> ma.make_mask(m)
    array([ True, False,  True,  True], dtype=bool)
    >>> m = [1, 0, 2, -3]
    >>> ma.make_mask(m)
    array([ True, False,  True,  True], dtype=bool)

    Effect of the `shrink` parameter.

    >>> m = np.zeros(4)
    >>> m
    array([ 0.,  0.,  0.,  0.])
    >>> ma.make_mask(m)
    False
    >>> ma.make_mask(m, shrink=False)
    array([False, False, False, False], dtype=bool)

    Using a flexible `dtype`.

    >>> m = [1, 0, 1, 1]
    >>> n = [0, 1, 0, 0]
    >>> arr = []
    >>> for man, mouse in zip(m, n):
    ...     arr.append((man, mouse))
    >>> arr
    [(1, 0), (0, 1), (1, 0), (1, 0)]
    >>> dtype = np.dtype({'names':['man', 'mouse'],
                          'formats':[np.int, np.int]})
    >>> arr = np.array(arr, dtype=dtype)
    >>> arr
    array([(1, 0), (0, 1), (1, 0), (1, 0)],
          dtype=[('man', '<i4'), ('mouse', '<i4')])
    >>> ma.make_mask(arr, dtype=dtype)
    array([(True, False), (False, True), (True, False), (True, False)],
          dtype=[('man', '|b1'), ('mouse', '|b1')])

    """
    if m is nomask:
        return nomask
    elif isinstance(m, ndarray):
        # We won't return after this point to make sure we can shrink the mask
        # Fill the mask in case there are missing data
        m = filled(m, True)
        # Make sure the input dtype is valid
        dtype = make_mask_descr(dtype)
        if m.dtype == dtype:
            if copy:
                result = m.copy()
            else:
                result = m
        else:
            result = np.array(m, dtype=dtype, copy=copy)
    else:
        result = np.array(filled(m, True), dtype=MaskType)
    # Bas les masques !
    if shrink and (not result.dtype.names) and (not result.any()):
        return nomask
    else:
        return result


def make_mask_none(newshape, dtype=None):
    """
    Return a boolean mask of the given shape, filled with False.

    This function returns a boolean ndarray with all entries False, that can
    be used in common mask manipulations. If a complex dtype is specified, the
    type of each field is converted to a boolean type.

    Parameters
    ----------
    newshape : tuple
        A tuple indicating the shape of the mask.
    dtype : {None, dtype}, optional
        If None, use a MaskType instance. Otherwise, use a new datatype with
        the same fields as `dtype`, converted to boolean types.

    Returns
    -------
    result : ndarray
        An ndarray of appropriate shape and dtype, filled with False.

    See Also
    --------
    make_mask : Create a boolean mask from an array.
    make_mask_descr : Construct a dtype description list from a given dtype.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> ma.make_mask_none((3,))
    array([False, False, False], dtype=bool)

    Defining a more complex dtype.

    >>> dtype = np.dtype({'names':['foo', 'bar'],
                          'formats':[np.float32, np.int]})
    >>> dtype
    dtype([('foo', '<f4'), ('bar', '<i4')])
    >>> ma.make_mask_none((3,), dtype=dtype)
    array([(False, False), (False, False), (False, False)],
          dtype=[('foo', '|b1'), ('bar', '|b1')])

    """
    if dtype is None:
        result = np.zeros(newshape, dtype=MaskType)
    else:
        result = np.zeros(newshape, dtype=make_mask_descr(dtype))
    return result

def mask_or (m1, m2, copy=False, shrink=True):
    """
    Combine two masks with the ``logical_or`` operator.

    The result may be a view on `m1` or `m2` if the other is `nomask`
    (i.e. False).

    Parameters
    ----------
    m1, m2 : array_like
        Input masks.
    copy : bool, optional
        If copy is False and one of the inputs is `nomask`, return a view
        of the other input mask. Defaults to False.
    shrink : bool, optional
        Whether to shrink the output to `nomask` if all its values are
        False. Defaults to True.

    Returns
    -------
    mask : output mask
        The result masks values that are masked in either `m1` or `m2`.

    Raises
    ------
    ValueError
        If `m1` and `m2` have different flexible dtypes.

    Examples
    --------
    >>> m1 = np.ma.make_mask([0, 1, 1, 0])
    >>> m2 = np.ma.make_mask([1, 0, 0, 0])
    >>> np.ma.mask_or(m1, m2)
    array([ True,  True,  True, False], dtype=bool)

    """
    def _recursive_mask_or(m1, m2, newmask):
        names = m1.dtype.names
        for name in names:
            current1 = m1[name]
            if current1.dtype.names:
                _recursive_mask_or(current1, m2[name], newmask[name])
            else:
                umath.logical_or(current1, m2[name], newmask[name])
        return
    #
    if (m1 is nomask) or (m1 is False):
        dtype = getattr(m2, 'dtype', MaskType)
        return make_mask(m2, copy=copy, shrink=shrink, dtype=dtype)
    if (m2 is nomask) or (m2 is False):
        dtype = getattr(m1, 'dtype', MaskType)
        return make_mask(m1, copy=copy, shrink=shrink, dtype=dtype)
    if m1 is m2 and is_mask(m1):
        return m1
    (dtype1, dtype2) = (getattr(m1, 'dtype', None), getattr(m2, 'dtype', None))
    if (dtype1 != dtype2):
        raise ValueError("Incompatible dtypes '%s'<>'%s'" % (dtype1, dtype2))
    if dtype1.names:
        newmask = np.empty_like(m1)
        _recursive_mask_or(m1, m2, newmask)
        return newmask
    return make_mask(umath.logical_or(m1, m2), copy=copy, shrink=shrink)


def flatten_mask(mask):
    """
    Returns a completely flattened version of the mask, where nested fields
    are collapsed.

    Parameters
    ----------
    mask : array_like
        Input array, which will be interpreted as booleans.

    Returns
    -------
    flattened_mask : ndarray of bools
        The flattened input.

    Examples
    --------
    >>> mask = np.array([0, 0, 1], dtype=np.bool)
    >>> flatten_mask(mask)
    array([False, False,  True], dtype=bool)

    >>> mask = np.array([(0, 0), (0, 1)], dtype=[('a', bool), ('b', bool)])
    >>> flatten_mask(mask)
    array([False, False, False,  True], dtype=bool)

    >>> mdtype = [('a', bool), ('b', [('ba', bool), ('bb', bool)])]
    >>> mask = np.array([(0, (0, 0)), (0, (0, 1))], dtype=mdtype)
    >>> flatten_mask(mask)
    array([False, False, False, False, False,  True], dtype=bool)

    """
    #
    def _flatmask(mask):
        "Flatten the mask and returns a (maybe nested) sequence of booleans."
        mnames = mask.dtype.names
        if mnames:
            return [flatten_mask(mask[name]) for name in mnames]
        else:
            return mask
    #
    def _flatsequence(sequence):
        "Generates a flattened version of the sequence."
        try:
            for element in sequence:
                if hasattr(element, '__iter__'):
                    for f in _flatsequence(element):
                        yield f
                else:
                    yield element
        except TypeError:
            yield sequence
    #
    mask = np.asarray(mask)
    flattened = _flatsequence(_flatmask(mask))
    return np.array([_ for _ in flattened], dtype=bool)


def _check_mask_axis(mask, axis):
    "Check whether there are masked values along the given axis"
    if mask is not nomask:
        return mask.all(axis=axis)
    return nomask


#####--------------------------------------------------------------------------
#--- --- Masking functions ---
#####--------------------------------------------------------------------------

def masked_where(condition, a, copy=True):
    """
    Mask an array where a condition is met.

    Return `a` as an array masked where `condition` is True.
    Any masked values of `a` or `condition` are also masked in the output.

    Parameters
    ----------
    condition : array_like
        Masking condition.  When `condition` tests floating point values for
        equality, consider using ``masked_values`` instead.
    a : array_like
        Array to mask.
    copy : bool
        If True (default) make a copy of `a` in the result.  If False modify
        `a` in place and return a view.

    Returns
    -------
    result : MaskedArray
        The result of masking `a` where `condition` is True.

    See Also
    --------
    masked_values : Mask using floating point equality.
    masked_equal : Mask where equal to a given value.
    masked_not_equal : Mask where `not` equal to a given value.
    masked_less_equal : Mask where less than or equal to a given value.
    masked_greater_equal : Mask where greater than or equal to a given value.
    masked_less : Mask where less than a given value.
    masked_greater : Mask where greater than a given value.
    masked_inside : Mask inside a given interval.
    masked_outside : Mask outside a given interval.
    masked_invalid : Mask invalid values (NaNs or infs).

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(4)
    >>> a
    array([0, 1, 2, 3])
    >>> ma.masked_where(a <= 2, a)
    masked_array(data = [-- -- -- 3],
          mask = [ True  True  True False],
          fill_value=999999)

    Mask array `b` conditional on `a`.

    >>> b = ['a', 'b', 'c', 'd']
    >>> ma.masked_where(a == 2, b)
    masked_array(data = [a b -- d],
          mask = [False False  True False],
          fill_value=N/A)

    Effect of the `copy` argument.

    >>> c = ma.masked_where(a <= 2, a)
    >>> c
    masked_array(data = [-- -- -- 3],
          mask = [ True  True  True False],
          fill_value=999999)
    >>> c[0] = 99
    >>> c
    masked_array(data = [99 -- -- 3],
          mask = [False  True  True False],
          fill_value=999999)
    >>> a
    array([0, 1, 2, 3])
    >>> c = ma.masked_where(a <= 2, a, copy=False)
    >>> c[0] = 99
    >>> c
    masked_array(data = [99 -- -- 3],
          mask = [False  True  True False],
          fill_value=999999)
    >>> a
    array([99,  1,  2,  3])

    When `condition` or `a` contain masked values.

    >>> a = np.arange(4)
    >>> a = ma.masked_where(a == 2, a)
    >>> a
    masked_array(data = [0 1 -- 3],
          mask = [False False  True False],
          fill_value=999999)
    >>> b = np.arange(4)
    >>> b = ma.masked_where(b == 0, b)
    >>> b
    masked_array(data = [-- 1 2 3],
          mask = [ True False False False],
          fill_value=999999)
    >>> ma.masked_where(a == 3, b)
    masked_array(data = [-- 1 -- --],
          mask = [ True False  True  True],
          fill_value=999999)

    """
    # Make sure that condition is a valid standard-type mask.
    cond = make_mask(condition)
    a = np.array(a, copy=copy, subok=True)

    (cshape, ashape) = (cond.shape, a.shape)
    if cshape and cshape != ashape:
        raise IndexError("Inconsistant shape between the condition and the input"
                         " (got %s and %s)" % (cshape, ashape))
    if hasattr(a, '_mask'):
        cond = mask_or(cond, a._mask)
        cls = type(a)
    else:
        cls = MaskedArray
    result = a.view(cls)
    result._mask = cond
    return result


def masked_greater(x, value, copy=True):
    """
    Mask an array where greater than a given value.

    This function is a shortcut to ``masked_where``, with
    `condition` = (x > value).

    See Also
    --------
    masked_where : Mask where a condition is met.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(4)
    >>> a
    array([0, 1, 2, 3])
    >>> ma.masked_greater(a, 2)
    masked_array(data = [0 1 2 --],
          mask = [False False False  True],
          fill_value=999999)

    """
    return masked_where(greater(x, value), x, copy=copy)


def masked_greater_equal(x, value, copy=True):
    """
    Mask an array where greater than or equal to a given value.

    This function is a shortcut to ``masked_where``, with
    `condition` = (x >= value).

    See Also
    --------
    masked_where : Mask where a condition is met.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(4)
    >>> a
    array([0, 1, 2, 3])
    >>> ma.masked_greater_equal(a, 2)
    masked_array(data = [0 1 -- --],
          mask = [False False  True  True],
          fill_value=999999)

    """
    return masked_where(greater_equal(x, value), x, copy=copy)


def masked_less(x, value, copy=True):
    """
    Mask an array where less than a given value.

    This function is a shortcut to ``masked_where``, with
    `condition` = (x < value).

    See Also
    --------
    masked_where : Mask where a condition is met.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(4)
    >>> a
    array([0, 1, 2, 3])
    >>> ma.masked_less(a, 2)
    masked_array(data = [-- -- 2 3],
          mask = [ True  True False False],
          fill_value=999999)

    """
    return masked_where(less(x, value), x, copy=copy)


def masked_less_equal(x, value, copy=True):
    """
    Mask an array where less than or equal to a given value.

    This function is a shortcut to ``masked_where``, with
    `condition` = (x <= value).

    See Also
    --------
    masked_where : Mask where a condition is met.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(4)
    >>> a
    array([0, 1, 2, 3])
    >>> ma.masked_less_equal(a, 2)
    masked_array(data = [-- -- -- 3],
          mask = [ True  True  True False],
          fill_value=999999)

    """
    return masked_where(less_equal(x, value), x, copy=copy)


def masked_not_equal(x, value, copy=True):
    """
    Mask an array where `not` equal to a given value.

    This function is a shortcut to ``masked_where``, with
    `condition` = (x != value).

    See Also
    --------
    masked_where : Mask where a condition is met.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(4)
    >>> a
    array([0, 1, 2, 3])
    >>> ma.masked_not_equal(a, 2)
    masked_array(data = [-- -- 2 --],
          mask = [ True  True False  True],
          fill_value=999999)

    """
    return masked_where(not_equal(x, value), x, copy=copy)


def masked_equal(x, value, copy=True):
    """
    Mask an array where equal to a given value.

    This function is a shortcut to ``masked_where``, with
    `condition` = (x == value).  For floating point arrays,
    consider using ``masked_values(x, value)``.

    See Also
    --------
    masked_where : Mask where a condition is met.
    masked_values : Mask using floating point equality.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(4)
    >>> a
    array([0, 1, 2, 3])
    >>> ma.masked_equal(a, 2)
    masked_array(data = [0 1 -- 3],
          mask = [False False  True False],
          fill_value=999999)

    """
    # An alternative implementation relies on filling first: probably not needed.
    # d = filled(x, 0)
    # c = umath.equal(d, value)
    # m = mask_or(c, getmask(x))
    # return array(d, mask=m, copy=copy)
    output = masked_where(equal(x, value), x, copy=copy)
    output.fill_value = value
    return output


def masked_inside(x, v1, v2, copy=True):
    """
    Mask an array inside a given interval.

    Shortcut to ``masked_where``, where `condition` is True for `x` inside
    the interval [v1,v2] (v1 <= x <= v2).  The boundaries `v1` and `v2`
    can be given in either order.

    See Also
    --------
    masked_where : Mask where a condition is met.

    Notes
    -----
    The array `x` is prefilled with its filling value.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> x = [0.31, 1.2, 0.01, 0.2, -0.4, -1.1]
    >>> ma.masked_inside(x, -0.3, 0.3)
    masked_array(data = [0.31 1.2 -- -- -0.4 -1.1],
          mask = [False False  True  True False False],
          fill_value=1e+20)

    The order of `v1` and `v2` doesn't matter.

    >>> ma.masked_inside(x, 0.3, -0.3)
    masked_array(data = [0.31 1.2 -- -- -0.4 -1.1],
          mask = [False False  True  True False False],
          fill_value=1e+20)

    """
    if v2 < v1:
        (v1, v2) = (v2, v1)
    xf = filled(x)
    condition = (xf >= v1) & (xf <= v2)
    return masked_where(condition, x, copy=copy)


def masked_outside(x, v1, v2, copy=True):
    """
    Mask an array outside a given interval.

    Shortcut to ``masked_where``, where `condition` is True for `x` outside
    the interval [v1,v2] (x < v1)|(x > v2).
    The boundaries `v1` and `v2` can be given in either order.

    See Also
    --------
    masked_where : Mask where a condition is met.

    Notes
    -----
    The array `x` is prefilled with its filling value.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> x = [0.31, 1.2, 0.01, 0.2, -0.4, -1.1]
    >>> ma.masked_outside(x, -0.3, 0.3)
    masked_array(data = [-- -- 0.01 0.2 -- --],
          mask = [ True  True False False  True  True],
          fill_value=1e+20)

    The order of `v1` and `v2` doesn't matter.

    >>> ma.masked_outside(x, 0.3, -0.3)
    masked_array(data = [-- -- 0.01 0.2 -- --],
          mask = [ True  True False False  True  True],
          fill_value=1e+20)

    """
    if v2 < v1:
        (v1, v2) = (v2, v1)
    xf = filled(x)
    condition = (xf < v1) | (xf > v2)
    return masked_where(condition, x, copy=copy)


def masked_object(x, value, copy=True, shrink=True):
    """
    Mask the array `x` where the data are exactly equal to value.

    This function is similar to `masked_values`, but only suitable
    for object arrays: for floating point, use `masked_values` instead.

    Parameters
    ----------
    x : array_like
        Array to mask
    value : object
        Comparison value
    copy : {True, False}, optional
        Whether to return a copy of `x`.
    shrink : {True, False}, optional
        Whether to collapse a mask full of False to nomask

    Returns
    -------
    result : MaskedArray
        The result of masking `x` where equal to `value`.

    See Also
    --------
    masked_where : Mask where a condition is met.
    masked_equal : Mask where equal to a given value (integers).
    masked_values : Mask using floating point equality.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> food = np.array(['green_eggs', 'ham'], dtype=object)
    >>> # don't eat spoiled food
    >>> eat = ma.masked_object(food, 'green_eggs')
    >>> print eat
    [-- ham]
    >>> # plain ol` ham is boring
    >>> fresh_food = np.array(['cheese', 'ham', 'pineapple'], dtype=object)
    >>> eat = ma.masked_object(fresh_food, 'green_eggs')
    >>> print eat
    [cheese ham pineapple]

    Note that `mask` is set to ``nomask`` if possible.

    >>> eat
    masked_array(data = [cheese ham pineapple],
          mask = False,
          fill_value=?)

    """
    if isMaskedArray(x):
        condition = umath.equal(x._data, value)
        mask = x._mask
    else:
        condition = umath.equal(np.asarray(x), value)
        mask = nomask
    mask = mask_or(mask, make_mask(condition, shrink=shrink))
    return masked_array(x, mask=mask, copy=copy, fill_value=value)


def masked_values(x, value, rtol=1e-5, atol=1e-8, copy=True, shrink=True):
    """
    Mask using floating point equality.

    Return a MaskedArray, masked where the data in array `x` are approximately
    equal to `value`, i.e. where the following condition is True

    (abs(x - value) <= atol+rtol*abs(value))

    The fill_value is set to `value` and the mask is set to ``nomask`` if
    possible.  For integers, consider using ``masked_equal``.

    Parameters
    ----------
    x : array_like
        Array to mask.
    value : float
        Masking value.
    rtol : float, optional
        Tolerance parameter.
    atol : float, optional
        Tolerance parameter (1e-8).
    copy : bool, optional
        Whether to return a copy of `x`.
    shrink : bool, optional
        Whether to collapse a mask full of False to ``nomask``.

    Returns
    -------
    result : MaskedArray
        The result of masking `x` where approximately equal to `value`.

    See Also
    --------
    masked_where : Mask where a condition is met.
    masked_equal : Mask where equal to a given value (integers).

    Examples
    --------
    >>> import numpy.ma as ma
    >>> x = np.array([1, 1.1, 2, 1.1, 3])
    >>> ma.masked_values(x, 1.1)
    masked_array(data = [1.0 -- 2.0 -- 3.0],
          mask = [False  True False  True False],
          fill_value=1.1)

    Note that `mask` is set to ``nomask`` if possible.

    >>> ma.masked_values(x, 1.5)
    masked_array(data = [ 1.   1.1  2.   1.1  3. ],
          mask = False,
          fill_value=1.5)

    For integers, the fill value will be different in general to the
    result of ``masked_equal``.

    >>> x = np.arange(5)
    >>> x
    array([0, 1, 2, 3, 4])
    >>> ma.masked_values(x, 2)
    masked_array(data = [0 1 -- 3 4],
          mask = [False False  True False False],
          fill_value=2)
    >>> ma.masked_equal(x, 2)
    masked_array(data = [0 1 -- 3 4],
          mask = [False False  True False False],
          fill_value=999999)

    """
    mabs = umath.absolute
    xnew = filled(x, value)
    if issubclass(xnew.dtype.type, np.floating):
        condition = umath.less_equal(mabs(xnew - value), atol + rtol * mabs(value))
        mask = getattr(x, '_mask', nomask)
    else:
        condition = umath.equal(xnew, value)
        mask = nomask
    mask = mask_or(mask, make_mask(condition, shrink=shrink))
    return masked_array(xnew, mask=mask, copy=copy, fill_value=value)


def masked_invalid(a, copy=True):
    """
    Mask an array where invalid values occur (NaNs or infs).

    This function is a shortcut to ``masked_where``, with
    `condition` = ~(np.isfinite(a)). Any pre-existing mask is conserved.
    Only applies to arrays with a dtype where NaNs or infs make sense
    (i.e. floating point types), but accepts any array_like object.

    See Also
    --------
    masked_where : Mask where a condition is met.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.arange(5, dtype=np.float)
    >>> a[2] = np.NaN
    >>> a[3] = np.PINF
    >>> a
    array([  0.,   1.,  NaN,  Inf,   4.])
    >>> ma.masked_invalid(a)
    masked_array(data = [0.0 1.0 -- -- 4.0],
          mask = [False False  True  True False],
          fill_value=1e+20)

    """
    a = np.array(a, copy=copy, subok=True)
    mask = getattr(a, '_mask', None)
    if mask is not None:
        condition = ~(np.isfinite(getdata(a)))
        if mask is not nomask:
            condition |= mask
        cls = type(a)
    else:
        condition = ~(np.isfinite(a))
        cls = MaskedArray
    result = a.view(cls)
    result._mask = condition
    return result


#####--------------------------------------------------------------------------
#---- --- Printing options ---
#####--------------------------------------------------------------------------

class _MaskedPrintOption:
    """
    Handle the string used to represent missing data in a masked array.

    """
    def __init__ (self, display):
        "Create the masked_print_option object."
        self._display = display
        self._enabled = True

    def display(self):
        "Display the string to print for masked values."
        return self._display

    def set_display (self, s):
        "Set the string to print for masked values."
        self._display = s

    def enabled(self):
        "Is the use of the display value enabled?"
        return self._enabled

    def enable(self, shrink=1):
        "Set the enabling shrink to `shrink`."
        self._enabled = shrink

    def __str__ (self):
        return str(self._display)

    __repr__ = __str__

#if you single index into a masked location you get this object.
masked_print_option = _MaskedPrintOption('--')


def _recursive_printoption(result, mask, printopt):
    """
    Puts printoptions in result where mask is True.
    Private function allowing for recursion
    """
    names = result.dtype.names
    for name in names:
        (curdata, curmask) = (result[name], mask[name])
        if curdata.dtype.names:
            _recursive_printoption(curdata, curmask, printopt)
        else:
            np.copyto(curdata, printopt, where=curmask)
    return

_print_templates = dict(long_std="""\
masked_%(name)s(data =
 %(data)s,
       %(nlen)s mask =
 %(mask)s,
 %(nlen)s fill_value = %(fill)s)
""",
                        short_std="""\
masked_%(name)s(data = %(data)s,
       %(nlen)s mask = %(mask)s,
%(nlen)s  fill_value = %(fill)s)
""",
                        long_flx="""\
masked_%(name)s(data =
 %(data)s,
       %(nlen)s mask =
 %(mask)s,
%(nlen)s  fill_value = %(fill)s,
      %(nlen)s dtype = %(dtype)s)
""",
                        short_flx="""\
masked_%(name)s(data = %(data)s,
%(nlen)s        mask = %(mask)s,
%(nlen)s  fill_value = %(fill)s,
%(nlen)s       dtype = %(dtype)s)
""")

#####--------------------------------------------------------------------------
#---- --- MaskedArray class ---
#####--------------------------------------------------------------------------

def _recursive_filled(a, mask, fill_value):
    """
    Recursively fill `a` with `fill_value`.
    Private function
    """
    names = a.dtype.names
    for name in names:
        current = a[name]
        if current.dtype.names:
            _recursive_filled(current, mask[name], fill_value[name])
        else:
            np.copyto(current, fill_value[name], where=mask[name])



def flatten_structured_array(a):
    """
    Flatten a structured array.

    The data type of the output is chosen such that it can represent all of the
    (nested) fields.

    Parameters
    ----------
    a : structured array

    Returns
    -------
    output : masked array or ndarray
        A flattened masked array if the input is a masked array, otherwise a
        standard ndarray.

    Examples
    --------
    >>> ndtype = [('a', int), ('b', float)]
    >>> a = np.array([(1, 1), (2, 2)], dtype=ndtype)
    >>> flatten_structured_array(a)
    array([[1., 1.],
           [2., 2.]])

    """
    #
    def flatten_sequence(iterable):
        """Flattens a compound of nested iterables."""
        for elm in iter(iterable):
            if hasattr(elm, '__iter__'):
                for f in flatten_sequence(elm):
                    yield f
            else:
                yield elm
    #
    a = np.asanyarray(a)
    inishape = a.shape
    a = a.ravel()
    if isinstance(a, MaskedArray):
        out = np.array([tuple(flatten_sequence(d.item())) for d in a._data])
        out = out.view(MaskedArray)
        out._mask = np.array([tuple(flatten_sequence(d.item()))
                              for d in getmaskarray(a)])
    else:
        out = np.array([tuple(flatten_sequence(d.item())) for d in a])
    if len(inishape) > 1:
        newshape = list(out.shape)
        newshape[0] = inishape
        out.shape = tuple(flatten_sequence(newshape))
    return out



class _arraymethod(object):
    """
    Define a wrapper for basic array methods.

    Upon call, returns a masked array, where the new ``_data`` array is
    the output of the corresponding method called on the original
    ``_data``.

    If `onmask` is True, the new mask is the output of the method called
    on the initial mask. Otherwise, the new mask is just a reference
    to the initial mask.

    Attributes
    ----------
    _onmask : bool
        Holds the `onmask` parameter.
    obj : object
        The object calling `_arraymethod`.

    Parameters
    ----------
    funcname : str
        Name of the function to apply on data.
    onmask : bool
        Whether the mask must be processed also (True) or left
        alone (False). Default is True. Make available as `_onmask`
        attribute.

    """
    def __init__(self, funcname, onmask=True):
        self.__name__ = funcname
        self._onmask = onmask
        self.obj = None
        self.__doc__ = self.getdoc()
    #
    def getdoc(self):
        "Return the doc of the function (from the doc of the method)."
        methdoc = getattr(ndarray, self.__name__, None) or \
                  getattr(np, self.__name__, None)
        if methdoc is not None:
            return methdoc.__doc__
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__(self, *args, **params):
        methodname = self.__name__
        instance = self.obj
        # Fallback : if the instance has not been initialized, use the first arg
        if instance is None:
            args = list(args)
            instance = args.pop(0)
        data = instance._data
        mask = instance._mask
        cls = type(instance)
        result = getattr(data, methodname)(*args, **params).view(cls)
        result._update_from(instance)
        if result.ndim:
            if not self._onmask:
                result.__setmask__(mask)
            elif mask is not nomask:
                result.__setmask__(getattr(mask, methodname)(*args, **params))
        else:
            if mask.ndim and (not mask.dtype.names and mask.all()):
                return masked
        return result


class MaskedIterator(object):
    """
    Flat iterator object to iterate over masked arrays.

    A `MaskedIterator` iterator is returned by ``x.flat`` for any masked array
    `x`. It allows iterating over the array as if it were a 1-D array,
    either in a for-loop or by calling its `next` method.

    Iteration is done in C-contiguous style, with the last index varying the
    fastest. The iterator can also be indexed using basic slicing or
    advanced indexing.

    See Also
    --------
    MaskedArray.flat : Return a flat iterator over an array.
    MaskedArray.flatten : Returns a flattened copy of an array.

    Notes
    -----
    `MaskedIterator` is not exported by the `ma` module. Instead of
    instantiating a `MaskedIterator` directly, use `MaskedArray.flat`.

    Examples
    --------
    >>> x = np.ma.array(arange(6).reshape(2, 3))
    >>> fl = x.flat
    >>> type(fl)
    <class 'numpy.ma.core.MaskedIterator'>
    >>> for item in fl:
    ...     print item
    ...
    0
    1
    2
    3
    4
    5

    Extracting more than a single element b indexing the `MaskedIterator`
    returns a masked array:

    >>> fl[2:4]
    masked_array(data = [2 3],
                 mask = False,
           fill_value = 999999)

    """
    def __init__(self, ma):
        self.ma = ma
        self.dataiter = ma._data.flat
        #
        if ma._mask is nomask:
            self.maskiter = None
        else:
            self.maskiter = ma._mask.flat

    def __iter__(self):
        return self

    def __getitem__(self, indx):
        result = self.dataiter.__getitem__(indx).view(type(self.ma))
        if self.maskiter is not None:
            _mask = self.maskiter.__getitem__(indx)
            if isinstance(_mask, ndarray):
                # set shape to match that of data; this is needed for matrices
                _mask.shape = result.shape
                result._mask = _mask
            elif isinstance(_mask, np.void):
                return mvoid(result, mask=_mask, hardmask=self.ma._hardmask)
            elif _mask:  # Just a scalar, masked
                return masked
        return result

    ### This won't work is ravel makes a copy
    def __setitem__(self, index, value):
        self.dataiter[index] = getdata(value)
        if self.maskiter is not None:
            self.maskiter[index] = getmaskarray(value)

    def __next__(self):
        """
        Return the next value, or raise StopIteration.

        Examples
        --------
        >>> x = np.ma.array([3, 2], mask=[0, 1])
        >>> fl = x.flat
        >>> fl.next()
        3
        >>> fl.next()
        masked_array(data = --,
                     mask = True,
               fill_value = 1e+20)
        >>> fl.next()
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "/home/ralf/python/numpy/numpy/ma/core.py", line 2243, in next
            d = self.dataiter.next()
        StopIteration

        """
        d = next(self.dataiter)
        if self.maskiter is not None:
            m = next(self.maskiter)
            if isinstance(m, np.void):
                return mvoid(d, mask=m, hardmask=self.ma._hardmask)
            elif m:  # Just a scalar, masked
                return masked
        return d

    next = __next__


class MaskedArray(ndarray):
    """
    An array class with possibly masked values.

    Masked values of True exclude the corresponding element from any
    computation.

    Construction::

      x = MaskedArray(data, mask=nomask, dtype=None,
                      copy=False, subok=True, ndmin=0, fill_value=None,
                      keep_mask=True, hard_mask=None, shrink=True)

    Parameters
    ----------
    data : array_like
        Input data.
    mask : sequence, optional
        Mask. Must be convertible to an array of booleans with the same
        shape as `data`. True indicates a masked (i.e. invalid) data.
    dtype : dtype, optional
        Data type of the output.
        If `dtype` is None, the type of the data argument (``data.dtype``)
        is used. If `dtype` is not None and different from ``data.dtype``,
        a copy is performed.
    copy : bool, optional
        Whether to copy the input data (True), or to use a reference instead.
        Default is False.
    subok : bool, optional
        Whether to return a subclass of `MaskedArray` if possible (True) or a
        plain `MaskedArray`. Default is True.
    ndmin : int, optional
        Minimum number of dimensions. Default is 0.
    fill_value : scalar, optional
        Value used to fill in the masked values when necessary.
        If None, a default based on the data-type is used.
    keep_mask : bool, optional
        Whether to combine `mask` with the mask of the input data, if any
        (True), or to use only `mask` for the output (False). Default is True.
    hard_mask : bool, optional
        Whether to use a hard mask or not. With a hard mask, masked values
        cannot be unmasked. Default is False.
    shrink : bool, optional
        Whether to force compression of an empty mask. Default is True.

    """

    __array_priority__ = 15
    _defaultmask = nomask
    _defaulthardmask = False
    _baseclass = ndarray

    def __new__(cls, data=None, mask=nomask, dtype=None, copy=False,
                subok=True, ndmin=0, fill_value=None,
                keep_mask=True, hard_mask=None, shrink=True,
                **options):
        """
    Create a new masked array from scratch.

    Notes
    -----
    A masked array can also be created by taking a .view(MaskedArray).

        """
        # Process data............
        _data = np.array(data, dtype=dtype, copy=copy, subok=True, ndmin=ndmin)
        _baseclass = getattr(data, '_baseclass', type(_data))
        # Check that we're not erasing the mask..........
        if isinstance(data, MaskedArray) and (data.shape != _data.shape):
            copy = True
        # Careful, cls might not always be MaskedArray...
        if not isinstance(data, cls) or not subok:
            _data = ndarray.view(_data, cls)
        else:
            _data = ndarray.view(_data, type(data))
        # Backwards compatibility w/ numpy.core.ma .......
        if hasattr(data, '_mask') and not isinstance(data, ndarray):
            _data._mask = data._mask
            _sharedmask = True
        # Process mask ...............................
        # Number of named fields (or zero if none)
        names_ = _data.dtype.names or ()
        # Type of the mask
        if names_:
            mdtype = make_mask_descr(_data.dtype)
        else:
            mdtype = MaskType
        # Case 1. : no mask in input ............
        if mask is nomask:
             # Erase the current mask ?
            if not keep_mask:
                # With a reduced version
                if shrink:
                    _data._mask = nomask
                # With full version
                else:
                    _data._mask = np.zeros(_data.shape, dtype=mdtype)
            # Check whether we missed something
            elif isinstance(data, (tuple, list)):
                try:
                    # If data is a sequence of masked array
                    mask = np.array([getmaskarray(m) for m in data],
                                    dtype=mdtype)
                except ValueError:
                    # If data is nested
                    mask = nomask
                # Force shrinking of the mask if needed (and possible)
                if (mdtype == MaskType) and mask.any():
                    _data._mask = mask
                    _data._sharedmask = False
            else:
                if copy:
                    _data._mask = _data._mask.copy()
                    _data._sharedmask = False
                    # Reset the shape of the original mask
                    if getmask(data) is not nomask:
                        data._mask.shape = data.shape
                else:
                    _data._sharedmask = True
        # Case 2. : With a mask in input ........
        else:
            # Read the mask with the current mdtype
            try:
                mask = np.array(mask, copy=copy, dtype=mdtype)
            # Or assume it's a sequence of bool/int
            except TypeError:
                mask = np.array([tuple([m] * len(mdtype)) for m in mask],
                                 dtype=mdtype)
            # Make sure the mask and the data have the same shape
            if mask.shape != _data.shape:
                (nd, nm) = (_data.size, mask.size)
                if nm == 1:
                    mask = np.resize(mask, _data.shape)
                elif nm == nd:
                    mask = np.reshape(mask, _data.shape)
                else:
                    msg = "Mask and data not compatible: data size is %i, " + \
                          "mask size is %i."
                    raise MaskError(msg % (nd, nm))
                copy = True
            # Set the mask to the new value
            if _data._mask is nomask:
                _data._mask = mask
                _data._sharedmask = not copy
            else:
                if not keep_mask:
                    _data._mask = mask
                    _data._sharedmask = not copy
                else:
                    if names_:
                        def _recursive_or(a, b):
                            "do a|=b on each field of a, recursively"
                            for name in a.dtype.names:
                                (af, bf) = (a[name], b[name])
                                if af.dtype.names:
                                    _recursive_or(af, bf)
                                else:
                                    af |= bf
                            return
                        _recursive_or(_data._mask, mask)
                    else:
                        _data._mask = np.logical_or(mask, _data._mask)
                    _data._sharedmask = False
        # Update fill_value.......
        if fill_value is None:
            fill_value = getattr(data, '_fill_value', None)
        # But don't run the check unless we have something to check....
        if fill_value is not None:
            _data._fill_value = _check_fill_value(fill_value, _data.dtype)
        # Process extra options ..
        if hard_mask is None:
            _data._hardmask = getattr(data, '_hardmask', False)
        else:
            _data._hardmask = hard_mask
        _data._baseclass = _baseclass
        return _data
    #
    def _update_from(self, obj):
        """Copies some attributes of obj to self.
        """
        if obj is not None and isinstance(obj, ndarray):
            _baseclass = type(obj)
        else:
            _baseclass = ndarray
        # We need to copy the _basedict to avoid backward propagation
        _optinfo = {}
        _optinfo.update(getattr(obj, '_optinfo', {}))
        _optinfo.update(getattr(obj, '_basedict', {}))
        if not isinstance(obj, MaskedArray):
            _optinfo.update(getattr(obj, '__dict__', {}))
        _dict = dict(_fill_value=getattr(obj, '_fill_value', None),
                     _hardmask=getattr(obj, '_hardmask', False),
                     _sharedmask=getattr(obj, '_sharedmask', False),
                     _isfield=getattr(obj, '_isfield', False),
                     _baseclass=getattr(obj, '_baseclass', _baseclass),
                     _optinfo=_optinfo,
                     _basedict=_optinfo)
        self.__dict__.update(_dict)
        self.__dict__.update(_optinfo)
        return


    def __array_finalize__(self, obj):
        """Finalizes the masked array.
        """
        # Get main attributes .........
        self._update_from(obj)
        if isinstance(obj, ndarray):
            odtype = obj.dtype
            if odtype.names:
                _mask = getattr(obj, '_mask', make_mask_none(obj.shape, odtype))
            else:
                _mask = getattr(obj, '_mask', nomask)
        else:
            _mask = nomask
        self._mask = _mask
        # Finalize the mask ...........
        if self._mask is not nomask:
            try:
                self._mask.shape = self.shape
            except ValueError:
                self._mask = nomask
            except (TypeError, AttributeError):
                # When _mask.shape is not writable (because it's a void)
                pass
        # Finalize the fill_value for structured arrays
        if self.dtype.names:
            if self._fill_value is None:
                self._fill_value = _check_fill_value(None, self.dtype)
        return


    def __array_wrap__(self, obj, context=None):
        """
        Special hook for ufuncs.
        Wraps the numpy array and sets the mask according to context.
        """
        result = obj.view(type(self))
        result._update_from(self)
        #..........
        if context is not None:
            result._mask = result._mask.copy()
            (func, args, _) = context
            m = reduce(mask_or, [getmaskarray(arg) for arg in args])
            # Get the domain mask................
            domain = ufunc_domain.get(func, None)
            if domain is not None:
                # Take the domain, and make sure it's a ndarray
                if len(args) > 2:
                    d = filled(reduce(domain, args), True)
                else:
                    d = filled(domain(*args), True)
                # Fill the result where the domain is wrong
                try:
                    # Binary domain: take the last value
                    fill_value = ufunc_fills[func][-1]
                except TypeError:
                    # Unary domain: just use this one
                    fill_value = ufunc_fills[func]
                except KeyError:
                    # Domain not recognized, use fill_value instead
                    fill_value = self.fill_value
                result = result.copy()
                np.copyto(result, fill_value, where=d)
                # Update the mask
                if m is nomask:
                    if d is not nomask:
                        m = d
                else:
                    # Don't modify inplace, we risk back-propagation
                    m = (m | d)
            # Make sure the mask has the proper size
            if result.shape == () and m:
                return masked
            else:
                result._mask = m
                result._sharedmask = False
        #....
        return result


    def view(self, dtype=None, type=None, fill_value=None):
        """
        Return a view of the MaskedArray data

        Parameters
        ----------
        dtype : data-type or ndarray sub-class, optional
            Data-type descriptor of the returned view, e.g., float32 or int16.
            The default, None, results in the view having the same data-type
            as `a`. As with ``ndarray.view``, dtype can also be specified as
            an ndarray sub-class, which then specifies the type of the
            returned object (this is equivalent to setting the ``type``
            parameter).
        type : Python type, optional
            Type of the returned view, e.g., ndarray or matrix.  Again, the
            default None results in type preservation.

        Notes
        -----

        ``a.view()`` is used two different ways:

        ``a.view(some_dtype)`` or ``a.view(dtype=some_dtype)`` constructs a view
        of the array's memory with a different data-type.  This can cause a
        reinterpretation of the bytes of memory.

        ``a.view(ndarray_subclass)`` or ``a.view(type=ndarray_subclass)`` just
        returns an instance of `ndarray_subclass` that looks at the same array
        (same shape, dtype, etc.)  This does not cause a reinterpretation of the
        memory.

        If `fill_value` is not specified, but `dtype` is specified (and is not
        an ndarray sub-class), the `fill_value` of the MaskedArray will be
        reset. If neither `fill_value` nor `dtype` are specified (or if
        `dtype` is an ndarray sub-class), then the fill value is preserved.
        Finally, if `fill_value` is specified, but `dtype` is not, the fill
        value is set to the specified value.

        For ``a.view(some_dtype)``, if ``some_dtype`` has a different number of
        bytes per entry than the previous dtype (for example, converting a
        regular array to a structured array), then the behavior of the view
        cannot be predicted just from the superficial appearance of ``a`` (shown
        by ``print(a)``). It also depends on exactly how ``a`` is stored in
        memory. Therefore if ``a`` is C-ordered versus fortran-ordered, versus
        defined as a slice or transpose, etc., the view may give different
        results.
        """

        if dtype is None:
            if type is None:
                output = ndarray.view(self)
            else:
                output = ndarray.view(self, type)
        elif type is None:
            try:
                if issubclass(dtype, ndarray):
                    output = ndarray.view(self, dtype)
                    dtype = None
                else:
                    output = ndarray.view(self, dtype)
            except TypeError:
                output = ndarray.view(self, dtype)
        else:
            output = ndarray.view(self, dtype, type)
        # Should we update the mask ?
        if (getattr(output, '_mask', nomask) is not nomask):
            if dtype is None:
                dtype = output.dtype
            mdtype = make_mask_descr(dtype)
            output._mask = self._mask.view(mdtype, ndarray)
            # Try to reset the shape of the mask (if we don't have a void)
            try:
                output._mask.shape = output.shape
            except (AttributeError, TypeError):
                pass
        # Make sure to reset the _fill_value if needed
        if getattr(output, '_fill_value', None) is not None:
            if fill_value is None:
                if dtype is None:
                    pass # leave _fill_value as is
                else:
                    output._fill_value = None
            else:
                output.fill_value = fill_value
        return output
    view.__doc__ = ndarray.view.__doc__


    def astype(self, newtype):
        """
        Returns a copy of the MaskedArray cast to given newtype.

        Returns
        -------
        output : MaskedArray
            A copy of self cast to input newtype.
            The returned record shape matches self.shape.

        Examples
        --------
        >>> x = np.ma.array([[1,2,3.1],[4,5,6],[7,8,9]], mask=[0] + [1,0]*4)
        >>> print x
        [[1.0 -- 3.1]
         [-- 5.0 --]
         [7.0 -- 9.0]]
        >>> print x.astype(int32)
        [[1 -- 3]
         [-- 5 --]
         [7 -- 9]]

        """
        newtype = np.dtype(newtype)
        output = self._data.astype(newtype).view(type(self))
        output._update_from(self)
        names = output.dtype.names
        if names is None:
            output._mask = self._mask.astype(bool)
        else:
            if self._mask is nomask:
                output._mask = nomask
            else:
                output._mask = self._mask.astype([(n, bool) for n in names])
        # Don't check _fill_value if it's None, that'll speed things up
        if self._fill_value is not None:
            output._fill_value = _check_fill_value(self._fill_value, newtype)
        return output


    def __getitem__(self, indx):
        """x.__getitem__(y) <==> x[y]

        Return the item described by i, as a masked array.

        """
        # This test is useful, but we should keep things light...
#        if getmask(indx) is not nomask:
#            msg = "Masked arrays must be filled before they can be used as indices!"
#            raise IndexError(msg)
        _data = ndarray.view(self, ndarray)
        dout = ndarray.__getitem__(_data, indx)
        # We could directly use ndarray.__getitem__ on self...
        # But then we would have to modify __array_finalize__ to prevent the
        # mask of being reshaped if it hasn't been set up properly yet...
        # So it's easier to stick to the current version
        _mask = self._mask
        if not getattr(dout, 'ndim', False):
            # A record ................
            if isinstance(dout, np.void):
                mask = _mask[indx]
                # We should always re-cast to mvoid, otherwise users can
                # change masks on rows that already have masked values, but not
                # on rows that have no masked values, which is inconsistent.
                dout = mvoid(dout, mask=mask, hardmask=self._hardmask)
            # Just a scalar............
            elif _mask is not nomask and _mask[indx]:
                return masked
        else:
            # Force dout to MA ........
            dout = dout.view(type(self))
            # Inherit attributes from self
            dout._update_from(self)
            # Check the fill_value ....
            if isinstance(indx, basestring):
                if self._fill_value is not None:
                    dout._fill_value = self._fill_value[indx]
                dout._isfield = True
            # Update the mask if needed
            if _mask is not nomask:
                dout._mask = _mask[indx]
                dout._sharedmask = True
#               Note: Don't try to check for m.any(), that'll take too long...
        return dout

    def __setitem__(self, indx, value):
        """x.__setitem__(i, y) <==> x[i]=y

        Set item described by index. If value is masked, masks those
        locations.

        """
        if self is masked:
            raise MaskError('Cannot alter the masked element.')
        # This test is useful, but we should keep things light...
#        if getmask(indx) is not nomask:
#            msg = "Masked arrays must be filled before they can be used as indices!"
#            raise IndexError(msg)
        _data = ndarray.view(self, ndarray.__getattribute__(self, '_baseclass'))
        _mask = ndarray.__getattribute__(self, '_mask')
        if isinstance(indx, basestring):
            ndarray.__setitem__(_data, indx, value)
            if _mask is nomask:
                self._mask = _mask = make_mask_none(self.shape, self.dtype)
            _mask[indx] = getmask(value)
            return
        #........................................
        _dtype = ndarray.__getattribute__(_data, 'dtype')
        nbfields = len(_dtype.names or ())
        #........................................
        if value is masked:
            # The mask wasn't set: create a full version...
            if _mask is nomask:
                _mask = self._mask = make_mask_none(self.shape, _dtype)
            # Now, set the mask to its value.
            if nbfields:
                _mask[indx] = tuple([True] * nbfields)
            else:
                _mask[indx] = True
            if not self._isfield:
                self._sharedmask = False
            return
        #........................................
        # Get the _data part of the new value
        dval = value
        # Get the _mask part of the new value
        mval = getattr(value, '_mask', nomask)
        if nbfields and mval is nomask:
            mval = tuple([False] * nbfields)
        if _mask is nomask:
            # Set the data, then the mask
            ndarray.__setitem__(_data, indx, dval)
            if mval is not nomask:
                _mask = self._mask = make_mask_none(self.shape, _dtype)
                ndarray.__setitem__(_mask, indx, mval)
        elif not self._hardmask:
            # Unshare the mask if necessary to avoid propagation
            if not self._isfield:
                self.unshare_mask()
                _mask = ndarray.__getattribute__(self, '_mask')
            # Set the data, then the mask
            ndarray.__setitem__(_data, indx, dval)
            ndarray.__setitem__(_mask, indx, mval)
        elif hasattr(indx, 'dtype') and (indx.dtype == MaskType):
            indx = indx * umath.logical_not(_mask)
            ndarray.__setitem__(_data, indx, dval)
        else:
            if nbfields:
                err_msg = "Flexible 'hard' masks are not yet supported..."
                raise NotImplementedError(err_msg)
            mindx = mask_or(_mask[indx], mval, copy=True)
            dindx = self._data[indx]
            if dindx.size > 1:
                np.copyto(dindx, dval, where=~mindx)
            elif mindx is nomask:
                dindx = dval
            ndarray.__setitem__(_data, indx, dindx)
            _mask[indx] = mindx
        return


    def __getslice__(self, i, j):
        """x.__getslice__(i, j) <==> x[i:j]

        Return the slice described by (i, j).  The use of negative
        indices is not supported.

        """
        return self.__getitem__(slice(i, j))

    def __setslice__(self, i, j, value):
        """x.__setslice__(i, j, value) <==> x[i:j]=value

    Set the slice (i,j) of a to value. If value is masked, mask
    those locations.

        """
        self.__setitem__(slice(i, j), value)


    def __setmask__(self, mask, copy=False):
        """Set the mask.

        """
        idtype = ndarray.__getattribute__(self, 'dtype')
        current_mask = ndarray.__getattribute__(self, '_mask')
        if mask is masked:
            mask = True
        # Make sure the mask is set
        if (current_mask is nomask):
            # Just don't do anything is there's nothing to do...
            if mask is nomask:
                return
            current_mask = self._mask = make_mask_none(self.shape, idtype)
        # No named fields.........
        if idtype.names is None:
            # Hardmask: don't unmask the data
            if self._hardmask:
                current_mask |= mask
            # Softmask: set everything to False
            # If it's obviously a compatible scalar, use a quick update
            # method...
            elif isinstance(mask, (int, float, np.bool_, np.number)):
                current_mask[...] = mask
            # ...otherwise fall back to the slower, general purpose way.
            else:
                current_mask.flat = mask
        # Named fields w/ ............
        else:
            mdtype = current_mask.dtype
            mask = np.array(mask, copy=False)
            # Mask is a singleton
            if not mask.ndim:
                # It's a boolean : make a record
                if mask.dtype.kind == 'b':
                    mask = np.array(tuple([mask.item()]*len(mdtype)),
                                    dtype=mdtype)
                # It's a record: make sure the dtype is correct
                else:
                    mask = mask.astype(mdtype)
            # Mask is a sequence
            else:
                # Make sure the new mask is a ndarray with the proper dtype
                try:
                    mask = np.array(mask, copy=copy, dtype=mdtype)
                # Or assume it's a sequence of bool/int
                except TypeError:
                    mask = np.array([tuple([m] * len(mdtype)) for m in mask],
                                    dtype=mdtype)
            # Hardmask: don't unmask the data
            if self._hardmask:
                for n in idtype.names:
                    current_mask[n] |= mask[n]
            # Softmask: set everything to False
            # If it's obviously a compatible scalar, use a quick update
            # method...
            elif isinstance(mask, (int, float, np.bool_, np.number)):
                current_mask[...] = mask
            # ...otherwise fall back to the slower, general purpose way.
            else:
                current_mask.flat = mask
        # Reshape if needed
        if current_mask.shape:
            current_mask.shape = self.shape
        return
    _set_mask = __setmask__
    #....
    def _get_mask(self):
        """Return the current mask.

        """
        # We could try to force a reshape, but that wouldn't work in some cases.
#        return self._mask.reshape(self.shape)
        return self._mask
    mask = property(fget=_get_mask, fset=__setmask__, doc="Mask")


    def _get_recordmask(self):
        """
    Return the mask of the records.
    A record is masked when all the fields are masked.

        """
        _mask = ndarray.__getattribute__(self, '_mask').view(ndarray)
        if _mask.dtype.names is None:
            return _mask
        return np.all(flatten_structured_array(_mask), axis= -1)


    def _set_recordmask(self):
        """Return the mask of the records.
    A record is masked when all the fields are masked.

        """
        raise NotImplementedError("Coming soon: setting the mask per records!")
    recordmask = property(fget=_get_recordmask)

    #............................................
    def harden_mask(self):
        """
        Force the mask to hard.

        Whether the mask of a masked array is hard or soft is determined by
        its `hardmask` property. `harden_mask` sets `hardmask` to True.

        See Also
        --------
        hardmask

        """
        self._hardmask = True
        return self

    def soften_mask(self):
        """
        Force the mask to soft.

        Whether the mask of a masked array is hard or soft is determined by
        its `hardmask` property. `soften_mask` sets `hardmask` to False.

        See Also
        --------
        hardmask

        """
        self._hardmask = False
        return self

    hardmask = property(fget=lambda self: self._hardmask,
                        doc="Hardness of the mask")


    def unshare_mask(self):
        """
        Copy the mask and set the sharedmask flag to False.

        Whether the mask is shared between masked arrays can be seen from
        the `sharedmask` property. `unshare_mask` ensures the mask is not shared.
        A copy of the mask is only made if it was shared.

        See Also
        --------
        sharedmask

        """
        if self._sharedmask:
            self._mask = self._mask.copy()
            self._sharedmask = False
        return self

    sharedmask = property(fget=lambda self: self._sharedmask,
                          doc="Share status of the mask (read-only).")

    def shrink_mask(self):
        """
        Reduce a mask to nomask when possible.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        >>> x = np.ma.array([[1,2 ], [3, 4]], mask=[0]*4)
        >>> x.mask
        array([[False, False],
               [False, False]], dtype=bool)
        >>> x.shrink_mask()
        >>> x.mask
        False

        """
        m = self._mask
        if m.ndim and not m.any():
            self._mask = nomask
        return self

    #............................................

    baseclass = property(fget=lambda self:self._baseclass,
                         doc="Class of the underlying data (read-only).")

    def _get_data(self):
        """Return the current data, as a view of the original
        underlying data.

        """
        return ndarray.view(self, self._baseclass)
    _data = property(fget=_get_data)
    data = property(fget=_get_data)

    def _get_flat(self):
        "Return a flat iterator."
        return MaskedIterator(self)
    #
    def _set_flat (self, value):
        "Set a flattened version of self to value."
        y = self.ravel()
        y[:] = value
    #
    flat = property(fget=_get_flat, fset=_set_flat,
                    doc="Flat version of the array.")


    def get_fill_value(self):
        """
        Return the filling value of the masked array.

        Returns
        -------
        fill_value : scalar
            The filling value.

        Examples
        --------
        >>> for dt in [np.int32, np.int64, np.float64, np.complex128]:
        ...     np.ma.array([0, 1], dtype=dt).get_fill_value()
        ...
        999999
        999999
        1e+20
        (1e+20+0j)

        >>> x = np.ma.array([0, 1.], fill_value=-np.inf)
        >>> x.get_fill_value()
        -inf

        """
        if self._fill_value is None:
            self._fill_value = _check_fill_value(None, self.dtype)
        return self._fill_value[()]

    def set_fill_value(self, value=None):
        """
        Set the filling value of the masked array.

        Parameters
        ----------
        value : scalar, optional
            The new filling value. Default is None, in which case a default
            based on the data type is used.

        See Also
        --------
        ma.set_fill_value : Equivalent function.

        Examples
        --------
        >>> x = np.ma.array([0, 1.], fill_value=-np.inf)
        >>> x.fill_value
        -inf
        >>> x.set_fill_value(np.pi)
        >>> x.fill_value
        3.1415926535897931

        Reset to default:

        >>> x.set_fill_value()
        >>> x.fill_value
        1e+20

        """
        target = _check_fill_value(value, self.dtype)
        _fill_value = self._fill_value
        if _fill_value is None:
            # Create the attribute if it was undefined
            self._fill_value = target
        else:
            # Don't overwrite the attribute, just fill it (for propagation)
            _fill_value[()] = target

    fill_value = property(fget=get_fill_value, fset=set_fill_value,
                          doc="Filling value.")


    def filled(self, fill_value=None):
        """
        Return a copy of self, with masked values filled with a given value.

        Parameters
        ----------
        fill_value : scalar, optional
            The value to use for invalid entries (None by default).
            If None, the `fill_value` attribute of the array is used instead.

        Returns
        -------
        filled_array : ndarray
            A copy of ``self`` with invalid entries replaced by *fill_value*
            (be it the function argument or the attribute of ``self``.

        Notes
        -----
        The result is **not** a MaskedArray!

        Examples
        --------
        >>> x = np.ma.array([1,2,3,4,5], mask=[0,0,1,0,1], fill_value=-999)
        >>> x.filled()
        array([1, 2, -999, 4, -999])
        >>> type(x.filled())
        <type 'numpy.ndarray'>

        Subclassing is preserved. This means that if the data part of the masked
        array is a matrix, `filled` returns a matrix:

        >>> x = np.ma.array(np.matrix([[1, 2], [3, 4]]), mask=[[0, 1], [1, 0]])
        >>> x.filled()
        matrix([[     1, 999999],
                [999999,      4]])

        """
        m = self._mask
        if m is nomask:
            return self._data
        #
        if fill_value is None:
            fill_value = self.fill_value
        else:
            fill_value = _check_fill_value(fill_value, self.dtype)
        #
        if self is masked_singleton:
            return np.asanyarray(fill_value)
        #
        if m.dtype.names:
            result = self._data.copy('K')
            _recursive_filled(result, self._mask, fill_value)
        elif not m.any():
            return self._data
        else:
            result = self._data.copy('K')
            try:
                np.copyto(result, fill_value, where=m)
            except (TypeError, AttributeError):
                fill_value = narray(fill_value, dtype=object)
                d = result.astype(object)
                result = np.choose(m, (d, fill_value))
            except IndexError:
                #ok, if scalar
                if self._data.shape:
                    raise
                elif m:
                    result = np.array(fill_value, dtype=self.dtype)
                else:
                    result = self._data
        return result

    def compressed(self):
        """
        Return all the non-masked data as a 1-D array.

        Returns
        -------
        data : ndarray
            A new `ndarray` holding the non-masked data is returned.

        Notes
        -----
        The result is **not** a MaskedArray!

        Examples
        --------
        >>> x = np.ma.array(np.arange(5), mask=[0]*2 + [1]*3)
        >>> x.compressed()
        array([0, 1])
        >>> type(x.compressed())
        <type 'numpy.ndarray'>

        """
        data = ndarray.ravel(self._data)
        if self._mask is not nomask:
            data = data.compress(np.logical_not(ndarray.ravel(self._mask)))
        return data


    def compress(self, condition, axis=None, out=None):
        """
        Return `a` where condition is ``True``.

        If condition is a `MaskedArray`, missing values are considered
        as ``False``.

        Parameters
        ----------
        condition : var
            Boolean 1-d array selecting which entries to return. If len(condition)
            is less than the size of a along the axis, then output is truncated
            to length of condition array.
        axis : {None, int}, optional
            Axis along which the operation must be performed.
        out : {None, ndarray}, optional
            Alternative output array in which to place the result. It must have
            the same shape as the expected output but the type will be cast if
            necessary.

        Returns
        -------
        result : MaskedArray
            A :class:`MaskedArray` object.

        Notes
        -----
        Please note the difference with :meth:`compressed` !
        The output of :meth:`compress` has a mask, the output of
        :meth:`compressed` does not.

        Examples
        --------
        >>> x = np.ma.array([[1,2,3],[4,5,6],[7,8,9]], mask=[0] + [1,0]*4)
        >>> print x
        [[1 -- 3]
         [-- 5 --]
         [7 -- 9]]
        >>> x.compress([1, 0, 1])
        masked_array(data = [1 3],
              mask = [False False],
              fill_value=999999)

        >>> x.compress([1, 0, 1], axis=1)
        masked_array(data =
         [[1 3]
         [-- --]
         [7 9]],
              mask =
         [[False False]
         [ True  True]
         [False False]],
              fill_value=999999)

        """
        # Get the basic components
        (_data, _mask) = (self._data, self._mask)
        # Force the condition to a regular ndarray (forget the missing values...)
        condition = np.array(condition, copy=False, subok=False)
        #
        _new = _data.compress(condition, axis=axis, out=out).view(type(self))
        _new._update_from(self)
        if _mask is not nomask:
            _new._mask = _mask.compress(condition, axis=axis)
        return _new

    #............................................
    def __str__(self):
        """String representation.

        """
        if masked_print_option.enabled():
            f = masked_print_option
            if self is masked:
                return str(f)
            m = self._mask
            if m is nomask:
                res = self._data
            else:
                if m.shape == ():
                    if m.dtype.names:
                        m = m.view((bool, len(m.dtype)))
                        if m.any():
                            return str(tuple((f if _m else _d) for _d, _m in
                                             zip(self._data.tolist(), m)))
                        else:
                            return str(self._data)
                    elif m:
                        return str(f)
                    else:
                        return str(self._data)
                # convert to object array to make filled work
                names = self.dtype.names
                if names is None:
                    res = self._data.astype("O")
                    res.view(ndarray)[m] = f
                else:
                    rdtype = _recursive_make_descr(self.dtype, "O")
                    res = self._data.astype(rdtype)
                    _recursive_printoption(res, m, f)
        else:
            res = self.filled(self.fill_value)
        return str(res)

    def __repr__(self):
        """Literal string representation.

        """
        n = len(self.shape)
        if self._baseclass is np.ndarray:
            name = 'array'
        else:
            name = self._baseclass.__name__

        parameters = dict(name=name, nlen=" " * len(name),
                          data=str(self), mask=str(self._mask),
                          fill=str(self.fill_value), dtype=str(self.dtype))
        if self.dtype.names:
            if n <= 1:
                return _print_templates['short_flx'] % parameters
            return _print_templates['long_flx'] % parameters
        elif n <= 1:
            return _print_templates['short_std'] % parameters
        return _print_templates['long_std'] % parameters

    def __eq__(self, other):
        "Check whether other equals self elementwise"
        if self is masked:
            return masked
        omask = getattr(other, '_mask', nomask)
        if omask is nomask:
            check = ndarray.__eq__(self.filled(0), other)
            try:
                check = check.view(type(self))
                check._mask = self._mask
            except AttributeError:
                # Dang, we have a bool instead of an array: return the bool
                return check
        else:
            odata = filled(other, 0)
            check = ndarray.__eq__(self.filled(0), odata).view(type(self))
            if self._mask is nomask:
                check._mask = omask
            else:
                mask = mask_or(self._mask, omask)
                if mask.dtype.names:
                    if mask.size > 1:
                        axis = 1
                    else:
                        axis = None
                    try:
                        mask = mask.view((bool_, len(self.dtype))).all(axis)
                    except ValueError:
                        mask = np.all([[f[n].all() for n in mask.dtype.names]
                                        for f in mask], axis=axis)
                check._mask = mask
        return check
    #
    def __ne__(self, other):
        "Check whether other doesn't equal self elementwise"
        if self is masked:
            return masked
        omask = getattr(other, '_mask', nomask)
        if omask is nomask:
            check = ndarray.__ne__(self.filled(0), other)
            try:
                check = check.view(type(self))
                check._mask = self._mask
            except AttributeError:
                # In case check is a boolean (or a numpy.bool)
                return check
        else:
            odata = filled(other, 0)
            check = ndarray.__ne__(self.filled(0), odata).view(type(self))
            if self._mask is nomask:
                check._mask = omask
            else:
                mask = mask_or(self._mask, omask)
                if mask.dtype.names:
                    if mask.size > 1:
                        axis = 1
                    else:
                        axis = None
                    try:
                        mask = mask.view((bool_, len(self.dtype))).all(axis)
                    except ValueError:
                        mask = np.all([[f[n].all() for n in mask.dtype.names]
                                        for f in mask], axis=axis)
                check._mask = mask
        return check
    #
    def __add__(self, other):
        "Add other to self, and return a new masked array."
        return add(self, other)
    #
    def __radd__(self, other):
        "Add other to self, and return a new masked array."
        return add(self, other)
    #
    def __sub__(self, other):
        "Subtract other to self, and return a new masked array."
        return subtract(self, other)
    #
    def __rsub__(self, other):
        "Subtract other to self, and return a new masked array."
        return subtract(other, self)
    #
    def __mul__(self, other):
        "Multiply other by self, and return a new masked array."
        return multiply(self, other)
    #
    def __rmul__(self, other):
        "Multiply other by self, and return a new masked array."
        return multiply(self, other)
    #
    def __div__(self, other):
        "Divide other into self, and return a new masked array."
        return divide(self, other)
    #
    def __truediv__(self, other):
        "Divide other into self, and return a new masked array."
        return true_divide(self, other)
    #
    def __rtruediv__(self, other):
        "Divide other into self, and return a new masked array."
        return true_divide(other, self)
    #
    def __floordiv__(self, other):
        "Divide other into self, and return a new masked array."
        return floor_divide(self, other)
    #
    def __rfloordiv__(self, other):
        "Divide other into self, and return a new masked array."
        return floor_divide(other, self)
    #
    def __pow__(self, other):
        "Raise self to the power other, masking the potential NaNs/Infs"
        return power(self, other)
    #
    def __rpow__(self, other):
        "Raise self to the power other, masking the potential NaNs/Infs"
        return power(other, self)
    #............................................
    def __iadd__(self, other):
        "Add other to self in-place."
        m = getmask(other)
        if self._mask is nomask:
            if m is not nomask and m.any():
                self._mask = make_mask_none(self.shape, self.dtype)
                self._mask += m
        else:
            if m is not nomask:
                self._mask += m
        ndarray.__iadd__(self._data, np.where(self._mask, 0, getdata(other)))
        return self
    #....
    def __isub__(self, other):
        "Subtract other from self in-place."
        m = getmask(other)
        if self._mask is nomask:
            if m is not nomask and m.any():
                self._mask = make_mask_none(self.shape, self.dtype)
                self._mask += m
        elif m is not nomask:
            self._mask += m
        ndarray.__isub__(self._data, np.where(self._mask, 0, getdata(other)))
        return self
    #....
    def __imul__(self, other):
        "Multiply self by other in-place."
        m = getmask(other)
        if self._mask is nomask:
            if m is not nomask and m.any():
                self._mask = make_mask_none(self.shape, self.dtype)
                self._mask += m
        elif m is not nomask:
            self._mask += m
        ndarray.__imul__(self._data, np.where(self._mask, 1, getdata(other)))
        return self
    #....
    def __idiv__(self, other):
        "Divide self by other in-place."
        other_data = getdata(other)
        dom_mask = _DomainSafeDivide().__call__(self._data, other_data)
        other_mask = getmask(other)
        new_mask = mask_or(other_mask, dom_mask)
        # The following 3 lines control the domain filling
        if dom_mask.any():
            (_, fval) = ufunc_fills[np.divide]
            other_data = np.where(dom_mask, fval, other_data)
#        self._mask = mask_or(self._mask, new_mask)
        self._mask |= new_mask
        ndarray.__idiv__(self._data, np.where(self._mask, 1, other_data))
        return self
    #....
    def __ifloordiv__(self, other):
        "Floor divide self by other in-place."
        other_data = getdata(other)
        dom_mask = _DomainSafeDivide().__call__(self._data, other_data)
        other_mask = getmask(other)
        new_mask = mask_or(other_mask, dom_mask)
        # The following 3 lines control the domain filling
        if dom_mask.any():
            (_, fval) = ufunc_fills[np.floor_divide]
            other_data = np.where(dom_mask, fval, other_data)
#        self._mask = mask_or(self._mask, new_mask)
        self._mask |= new_mask
        ndarray.__ifloordiv__(self._data, np.where(self._mask, 1, other_data))
        return self
    #....
    def __itruediv__(self, other):
        "True divide self by other in-place."
        other_data = getdata(other)
        dom_mask = _DomainSafeDivide().__call__(self._data, other_data)
        other_mask = getmask(other)
        new_mask = mask_or(other_mask, dom_mask)
        # The following 3 lines control the domain filling
        if dom_mask.any():
            (_, fval) = ufunc_fills[np.true_divide]
            other_data = np.where(dom_mask, fval, other_data)
#        self._mask = mask_or(self._mask, new_mask)
        self._mask |= new_mask
        ndarray.__itruediv__(self._data, np.where(self._mask, 1, other_data))
        return self
    #...
    def __ipow__(self, other):
        "Raise self to the power other, in place."
        other_data = getdata(other)
        other_mask = getmask(other)
        with np.errstate(divide='ignore', invalid='ignore'):
            ndarray.__ipow__(self._data, np.where(self._mask, 1, other_data))
        invalid = np.logical_not(np.isfinite(self._data))
        if invalid.any():
            if self._mask is not nomask:
                self._mask |= invalid
            else:
                self._mask = invalid
            np.copyto(self._data, self.fill_value, where=invalid)
        new_mask = mask_or(other_mask, invalid)
        self._mask = mask_or(self._mask, new_mask)
        return self
    #............................................
    def __float__(self):
        "Convert to float."
        if self.size > 1:
            raise TypeError("Only length-1 arrays can be converted "
                            "to Python scalars")
        elif self._mask:
            warnings.warn("Warning: converting a masked element to nan.")
            return np.nan
        return float(self.item())

    def __int__(self):
        "Convert to int."
        if self.size > 1:
            raise TypeError("Only length-1 arrays can be converted "
                            "to Python scalars")
        elif self._mask:
            raise MaskError('Cannot convert masked element to a Python int.')
        return int(self.item())


    def get_imag(self):
        """
        Return the imaginary part of the masked array.

        The returned array is a view on the imaginary part of the `MaskedArray`
        whose `get_imag` method is called.

        Parameters
        ----------
        None

        Returns
        -------
        result : MaskedArray
            The imaginary part of the masked array.

        See Also
        --------
        get_real, real, imag

        Examples
        --------
        >>> x = np.ma.array([1+1.j, -2j, 3.45+1.6j], mask=[False, True, False])
        >>> x.get_imag()
        masked_array(data = [1.0 -- 1.6],
                     mask = [False  True False],
               fill_value = 1e+20)

        """
        result = self._data.imag.view(type(self))
        result.__setmask__(self._mask)
        return result
    imag = property(fget=get_imag, doc="Imaginary part.")

    def get_real(self):
        """
        Return the real part of the masked array.

        The returned array is a view on the real part of the `MaskedArray`
        whose `get_real` method is called.

        Parameters
        ----------
        None

        Returns
        -------
        result : MaskedArray
            The real part of the masked array.

        See Also
        --------
        get_imag, real, imag

        Examples
        --------
        >>> x = np.ma.array([1+1.j, -2j, 3.45+1.6j], mask=[False, True, False])
        >>> x.get_real()
        masked_array(data = [1.0 -- 3.45],
                     mask = [False  True False],
               fill_value = 1e+20)

        """
        result = self._data.real.view(type(self))
        result.__setmask__(self._mask)
        return result
    real = property(fget=get_real, doc="Real part")


    #............................................
    def count(self, axis=None):
        """
        Count the non-masked elements of the array along the given axis.

        Parameters
        ----------
        axis : int, optional
            Axis along which to count the non-masked elements. If `axis` is
            `None`, all non-masked elements are counted.

        Returns
        -------
        result : int or ndarray
            If `axis` is `None`, an integer count is returned. When `axis` is
            not `None`, an array with shape determined by the lengths of the
            remaining axes, is returned.

        See Also
        --------
        count_masked : Count masked elements in array or along a given axis.

        Examples
        --------
        >>> import numpy.ma as ma
        >>> a = ma.arange(6).reshape((2, 3))
        >>> a[1, :] = ma.masked
        >>> a
        masked_array(data =
         [[0 1 2]
         [-- -- --]],
                     mask =
         [[False False False]
         [ True  True  True]],
               fill_value = 999999)
        >>> a.count()
        3

        When the `axis` keyword is specified an array of appropriate size is
        returned.

        >>> a.count(axis=0)
        array([1, 1, 1])
        >>> a.count(axis=1)
        array([3, 0])

        """
        m = self._mask
        s = self.shape
        if m is nomask:
            if axis is None:
                return self.size
            else:
                n = s[axis]
                t = list(s)
                del t[axis]
                return np.full(t, n, dtype=np.intp)
        n1 = np.size(m, axis)
        n2 = np.sum(m, axis=axis, dtype=np.intp)
        if axis is None:
            return (n1 - n2)
        else:
            return narray(n1 - n2)
    #............................................
    flatten = _arraymethod('flatten')
    #
    def ravel(self):
        """
        Returns a 1D version of self, as a view.

        Returns
        -------
        MaskedArray
            Output view is of shape ``(self.size,)`` (or
            ``(np.ma.product(self.shape),)``).

        Examples
        --------
        >>> x = np.ma.array([[1,2,3],[4,5,6],[7,8,9]], mask=[0] + [1,0]*4)
        >>> print x
        [[1 -- 3]
         [-- 5 --]
         [7 -- 9]]
        >>> print x.ravel()
        [1 -- 3 -- 5 -- 7 -- 9]

        """
        r = ndarray.ravel(self._data).view(type(self))
        r._update_from(self)
        if self._mask is not nomask:
            r._mask = ndarray.ravel(self._mask).reshape(r.shape)
        else:
            r._mask = nomask
        return r
    #
    repeat = _arraymethod('repeat')
    #
    def reshape (self, *s, **kwargs):
        """
        Give a new shape to the array without changing its data.

        Returns a masked array containing the same data, but with a new shape.
        The result is a view on the original array; if this is not possible, a
        ValueError is raised.

        Parameters
        ----------
        shape : int or tuple of ints
            The new shape should be compatible with the original shape. If an
            integer is supplied, then the result will be a 1-D array of that
            length.
        order : {'C', 'F'}, optional
            Determines whether the array data should be viewed as in C
            (row-major) or FORTRAN (column-major) order.

        Returns
        -------
        reshaped_array : array
            A new view on the array.

        See Also
        --------
        reshape : Equivalent function in the masked array module.
        numpy.ndarray.reshape : Equivalent method on ndarray object.
        numpy.reshape : Equivalent function in the NumPy module.

        Notes
        -----
        The reshaping operation cannot guarantee that a copy will not be made,
        to modify the shape in place, use ``a.shape = s``

        Examples
        --------
        >>> x = np.ma.array([[1,2],[3,4]], mask=[1,0,0,1])
        >>> print x
        [[-- 2]
         [3 --]]
        >>> x = x.reshape((4,1))
        >>> print x
        [[--]
         [2]
         [3]
         [--]]

        """
        kwargs.update(order=kwargs.get('order', 'C'))
        result = self._data.reshape(*s, **kwargs).view(type(self))
        result._update_from(self)
        mask = self._mask
        if mask is not nomask:
            result._mask = mask.reshape(*s, **kwargs)
        return result
    #
    def resize(self, newshape, refcheck=True, order=False):
        """
        .. warning::

            This method does nothing, except raise a ValueError exception. A
            masked array does not own its data and therefore cannot safely be
            resized in place. Use the `numpy.ma.resize` function instead.

        This method is difficult to implement safely and may be deprecated in
        future releases of NumPy.

        """
        # Note : the 'order' keyword looks broken, let's just drop it
#        try:
#            ndarray.resize(self, newshape, refcheck=refcheck)
#            if self.mask is not nomask:
#                self._mask.resize(newshape, refcheck=refcheck)
#        except ValueError:
#            raise ValueError("Cannot resize an array that has been referenced "
#                             "or is referencing another array in this way.\n"
#                             "Use the numpy.ma.resize function.")
#        return None
        errmsg = "A masked array does not own its data "\
                 "and therefore cannot be resized.\n" \
                 "Use the numpy.ma.resize function instead."
        raise ValueError(errmsg)
    #
    def put(self, indices, values, mode='raise'):
        """
        Set storage-indexed locations to corresponding values.

        Sets self._data.flat[n] = values[n] for each n in indices.
        If `values` is shorter than `indices` then it will repeat.
        If `values` has some masked values, the initial mask is updated
        in consequence, else the corresponding values are unmasked.

        Parameters
        ----------
        indices : 1-D array_like
            Target indices, interpreted as integers.
        values : array_like
            Values to place in self._data copy at target indices.
        mode : {'raise', 'wrap', 'clip'}, optional
            Specifies how out-of-bounds indices will behave.
            'raise' : raise an error.
            'wrap' : wrap around.
            'clip' : clip to the range.

        Notes
        -----
        `values` can be a scalar or length 1 array.

        Examples
        --------
        >>> x = np.ma.array([[1,2,3],[4,5,6],[7,8,9]], mask=[0] + [1,0]*4)
        >>> print x
        [[1 -- 3]
         [-- 5 --]
         [7 -- 9]]
        >>> x.put([0,4,8],[10,20,30])
        >>> print x
        [[10 -- 3]
         [-- 20 --]
         [7 -- 30]]

        >>> x.put(4,999)
        >>> print x
        [[10 -- 3]
         [-- 999 --]
         [7 -- 30]]

        """
        m = self._mask
        # Hard mask: Get rid of the values/indices that fall on masked data
        if self._hardmask and self._mask is not nomask:
            mask = self._mask[indices]
            indices = narray(indices, copy=False)
            values = narray(values, copy=False, subok=True)
            values.resize(indices.shape)
            indices = indices[~mask]
            values = values[~mask]
        #....
        self._data.put(indices, values, mode=mode)
        #....
        if m is nomask:
            m = getmask(values)
        else:
            m = m.copy()
            if getmask(values) is nomask:
                m.put(indices, False, mode=mode)
            else:
                m.put(indices, values._mask, mode=mode)
            m = make_mask(m, copy=False, shrink=True)
        self._mask = m
    #............................................
    def ids (self):
        """
        Return the addresses of the data and mask areas.

        Parameters
        ----------
        None

        Examples
        --------
        >>> x = np.ma.array([1, 2, 3], mask=[0, 1, 1])
        >>> x.ids()
        (166670640, 166659832)

        If the array has no mask, the address of `nomask` is returned. This address
        is typically not close to the data in memory:

        >>> x = np.ma.array([1, 2, 3])
        >>> x.ids()
        (166691080, 3083169284L)

        """
        if self._mask is nomask:
            return (self.ctypes.data, id(nomask))
        return (self.ctypes.data, self._mask.ctypes.data)

    def iscontiguous(self):
        """
        Return a boolean indicating whether the data is contiguous.

        Parameters
        ----------
        None

        Examples
        --------
        >>> x = np.ma.array([1, 2, 3])
        >>> x.iscontiguous()
        True

        `iscontiguous` returns one of the flags of the masked array:

        >>> x.flags
          C_CONTIGUOUS : True
          F_CONTIGUOUS : True
          OWNDATA : False
          WRITEABLE : True
          ALIGNED : True
          UPDATEIFCOPY : False

        """
        return self.flags['CONTIGUOUS']

    #............................................
    def all(self, axis=None, out=None):
        """
    Check if all of the elements of `a` are true.

    Performs a :func:`logical_and` over the given axis and returns the result.
    Masked values are considered as True during computation.
    For convenience, the output array is masked where ALL the values along the
    current axis are masked: if the output would have been a scalar and that
    all the values are masked, then the output is `masked`.

    Parameters
    ----------
    axis : {None, integer}
        Axis to perform the operation over.
        If None, perform over flattened array.
    out : {None, array}, optional
        Array into which the result can be placed. Its type is preserved
        and it must be of the right shape to hold the output.

    See Also
    --------
    all : equivalent function

    Examples
    --------
    >>> np.ma.array([1,2,3]).all()
    True
    >>> a = np.ma.array([1,2,3], mask=True)
    >>> (a.all() is np.ma.masked)
    True

        """
        mask = _check_mask_axis(self._mask, axis)
        if out is None:
            d = self.filled(True).all(axis=axis).view(type(self))
            if d.ndim:
                d.__setmask__(mask)
            elif mask:
                return masked
            return d
        self.filled(True).all(axis=axis, out=out)
        if isinstance(out, MaskedArray):
            if out.ndim or mask:
                out.__setmask__(mask)
        return out


    def any(self, axis=None, out=None):
        """
        Check if any of the elements of `a` are true.

        Performs a logical_or over the given axis and returns the result.
        Masked values are considered as False during computation.

        Parameters
        ----------
        axis : {None, integer}
            Axis to perform the operation over.
            If None, perform over flattened array and return a scalar.
        out : {None, array}, optional
            Array into which the result can be placed. Its type is preserved
            and it must be of the right shape to hold the output.

        See Also
        --------
        any : equivalent function

        """
        mask = _check_mask_axis(self._mask, axis)
        if out is None:
            d = self.filled(False).any(axis=axis).view(type(self))
            if d.ndim:
                d.__setmask__(mask)
            elif mask:
                d = masked
            return d
        self.filled(False).any(axis=axis, out=out)
        if isinstance(out, MaskedArray):
            if out.ndim or mask:
                out.__setmask__(mask)
        return out


    def nonzero(self):
        """
        Return the indices of unmasked elements that are not zero.

        Returns a tuple of arrays, one for each dimension, containing the
        indices of the non-zero elements in that dimension. The corresponding
        non-zero values can be obtained with::

            a[a.nonzero()]

        To group the indices by element, rather than dimension, use
        instead::

            np.transpose(a.nonzero())

        The result of this is always a 2d array, with a row for each non-zero
        element.

        Parameters
        ----------
        None

        Returns
        -------
        tuple_of_arrays : tuple
            Indices of elements that are non-zero.

        See Also
        --------
        numpy.nonzero :
            Function operating on ndarrays.
        flatnonzero :
            Return indices that are non-zero in the flattened version of the input
            array.
        ndarray.nonzero :
            Equivalent ndarray method.
        count_nonzero :
            Counts the number of non-zero elements in the input array.

        Examples
        --------
        >>> import numpy.ma as ma
        >>> x = ma.array(np.eye(3))
        >>> x
        masked_array(data =
         [[ 1.  0.  0.]
         [ 0.  1.  0.]
         [ 0.  0.  1.]],
              mask =
         False,
              fill_value=1e+20)
        >>> x.nonzero()
        (array([0, 1, 2]), array([0, 1, 2]))

        Masked elements are ignored.

        >>> x[1, 1] = ma.masked
        >>> x
        masked_array(data =
         [[1.0 0.0 0.0]
         [0.0 -- 0.0]
         [0.0 0.0 1.0]],
              mask =
         [[False False False]
         [False  True False]
         [False False False]],
              fill_value=1e+20)
        >>> x.nonzero()
        (array([0, 2]), array([0, 2]))

        Indices can also be grouped by element.

        >>> np.transpose(x.nonzero())
        array([[0, 0],
               [2, 2]])

        A common use for ``nonzero`` is to find the indices of an array, where
        a condition is True.  Given an array `a`, the condition `a` > 3 is a
        boolean array and since False is interpreted as 0, ma.nonzero(a > 3)
        yields the indices of the `a` where the condition is true.

        >>> a = ma.array([[1,2,3],[4,5,6],[7,8,9]])
        >>> a > 3
        masked_array(data =
         [[False False False]
         [ True  True  True]
         [ True  True  True]],
              mask =
         False,
              fill_value=999999)
        >>> ma.nonzero(a > 3)
        (array([1, 1, 1, 2, 2, 2]), array([0, 1, 2, 0, 1, 2]))

        The ``nonzero`` method of the condition array can also be called.

        >>> (a > 3).nonzero()
        (array([1, 1, 1, 2, 2, 2]), array([0, 1, 2, 0, 1, 2]))

        """
        return narray(self.filled(0), copy=False).nonzero()


    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None):
        """
        (this docstring should be overwritten)
        """
        #!!!: implement out + test!
        m = self._mask
        if m is nomask:
            result = super(MaskedArray, self).trace(offset=offset, axis1=axis1,
                                                    axis2=axis2, out=out)
            return result.astype(dtype)
        else:
            D = self.diagonal(offset=offset, axis1=axis1, axis2=axis2)
            return D.astype(dtype).filled(0).sum(axis=None, out=out)
    trace.__doc__ = ndarray.trace.__doc__

    def sum(self, axis=None, dtype=None, out=None):
        """
        Return the sum of the array elements over the given axis.
        Masked elements are set to 0 internally.

        Parameters
        ----------
        axis : {None, -1, int}, optional
            Axis along which the sum is computed. The default
            (`axis` = None) is to compute over the flattened array.
        dtype : {None, dtype}, optional
            Determines the type of the returned array and of the accumulator
            where the elements are summed. If dtype has the value None and
            the type of a is an integer type of precision less than the default
            platform integer, then the default platform integer precision is
            used.  Otherwise, the dtype is the same as that of a.
        out :  {None, ndarray}, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output
            but the type will be cast if necessary.

        Returns
        -------
        sum_along_axis : MaskedArray or scalar
            An array with the same shape as self, with the specified
            axis removed.   If self is a 0-d array, or if `axis` is None, a scalar
            is returned.  If an output array is specified, a reference to
            `out` is returned.

        Examples
        --------
        >>> x = np.ma.array([[1,2,3],[4,5,6],[7,8,9]], mask=[0] + [1,0]*4)
        >>> print x
        [[1 -- 3]
         [-- 5 --]
         [7 -- 9]]
        >>> print x.sum()
        25
        >>> print x.sum(axis=1)
        [4 5 16]
        >>> print x.sum(axis=0)
        [8 5 12]
        >>> print type(x.sum(axis=0, dtype=np.int64)[0])
        <type 'numpy.int64'>

        """
        _mask = ndarray.__getattribute__(self, '_mask')
        newmask = _check_mask_axis(_mask, axis)
        # No explicit output
        if out is None:
            result = self.filled(0).sum(axis, dtype=dtype)
            rndim = getattr(result, 'ndim', 0)
            if rndim:
                result = result.view(type(self))
                result.__setmask__(newmask)
            elif newmask:
                result = masked
            return result
        # Explicit output
        result = self.filled(0).sum(axis, dtype=dtype, out=out)
        if isinstance(out, MaskedArray):
            outmask = getattr(out, '_mask', nomask)
            if (outmask is nomask):
                outmask = out._mask = make_mask_none(out.shape)
            outmask.flat = newmask
        return out


    def cumsum(self, axis=None, dtype=None, out=None):
        """
        Return the cumulative sum of the elements along the given axis.
        The cumulative sum is calculated over the flattened array by
        default, otherwise over the specified axis.

        Masked values are set to 0 internally during the computation.
        However, their position is saved, and the result will be masked at
        the same locations.

        Parameters
        ----------
        axis : {None, -1, int}, optional
            Axis along which the sum is computed. The default (`axis` = None) is to
            compute over the flattened array. `axis` may be negative, in which case
            it counts from the   last to the first axis.
        dtype : {None, dtype}, optional
            Type of the returned array and of the accumulator in which the
            elements are summed.  If `dtype` is not specified, it defaults
            to the dtype of `a`, unless `a` has an integer dtype with a
            precision less than that of the default platform integer.  In
            that case, the default platform integer is used.
        out : ndarray, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output
            but the type will be cast if necessary.

        Returns
        -------
        cumsum : ndarray.
            A new array holding the result is returned unless ``out`` is
            specified, in which case a reference to ``out`` is returned.

        Notes
        -----
        The mask is lost if `out` is not a valid :class:`MaskedArray` !

        Arithmetic is modular when using integer types, and no error is
        raised on overflow.

        Examples
        --------
        >>> marr = np.ma.array(np.arange(10), mask=[0,0,0,1,1,1,0,0,0,0])
        >>> print marr.cumsum()
        [0 1 3 -- -- -- 9 16 24 33]

        """
        result = self.filled(0).cumsum(axis=axis, dtype=dtype, out=out)
        if out is not None:
            if isinstance(out, MaskedArray):
                out.__setmask__(self.mask)
            return out
        result = result.view(type(self))
        result.__setmask__(self._mask)
        return result


    def prod(self, axis=None, dtype=None, out=None):
        """
        Return the product of the array elements over the given axis.
        Masked elements are set to 1 internally for computation.

        Parameters
        ----------
        axis : {None, int}, optional
            Axis over which the product is taken. If None is used, then the
            product is over all the array elements.
        dtype : {None, dtype}, optional
            Determines the type of the returned array and of the accumulator
            where the elements are multiplied. If ``dtype`` has the value ``None``
            and the type of a is an integer type of precision less than the default
            platform integer, then the default platform integer precision is
            used.  Otherwise, the dtype is the same as that of a.
        out : {None, array}, optional
            Alternative output array in which to place the result. It must have
            the same shape as the expected output but the type will be cast if
            necessary.

        Returns
        -------
        product_along_axis : {array, scalar}, see dtype parameter above.
            Returns an array whose shape is the same as a with the specified
            axis removed. Returns a 0d array when a is 1d or axis=None.
            Returns a reference to the specified output array if specified.

        See Also
        --------
        prod : equivalent function

        Notes
        -----
        Arithmetic is modular when using integer types, and no error is raised
        on overflow.

        Examples
        --------
        >>> np.prod([1.,2.])
        2.0
        >>> np.prod([1.,2.], dtype=np.int32)
        2
        >>> np.prod([[1.,2.],[3.,4.]])
        24.0
        >>> np.prod([[1.,2.],[3.,4.]], axis=1)
        array([  2.,  12.])

        """
        _mask = ndarray.__getattribute__(self, '_mask')
        newmask = _check_mask_axis(_mask, axis)
        # No explicit output
        if out is None:
            result = self.filled(1).prod(axis, dtype=dtype)
            rndim = getattr(result, 'ndim', 0)
            if rndim:
                result = result.view(type(self))
                result.__setmask__(newmask)
            elif newmask:
                result = masked
            return result
        # Explicit output
        result = self.filled(1).prod(axis, dtype=dtype, out=out)
        if isinstance(out, MaskedArray):
            outmask = getattr(out, '_mask', nomask)
            if (outmask is nomask):
                outmask = out._mask = make_mask_none(out.shape)
            outmask.flat = newmask
        return out

    product = prod

    def cumprod(self, axis=None, dtype=None, out=None):
        """
        Return the cumulative product of the elements along the given axis.
        The cumulative product is taken over the flattened array by
        default, otherwise over the specified axis.

        Masked values are set to 1 internally during the computation.
        However, their position is saved, and the result will be masked at
        the same locations.

        Parameters
        ----------
        axis : {None, -1, int}, optional
            Axis along which the product is computed. The default
            (`axis` = None) is to compute over the flattened array.
        dtype : {None, dtype}, optional
            Determines the type of the returned array and of the accumulator
            where the elements are multiplied. If ``dtype`` has the value ``None``
            and the type of ``a`` is an integer type of precision less than the
            default platform integer, then the default platform integer precision
            is used.  Otherwise, the dtype is the same as that of ``a``.
        out : ndarray, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output
            but the type will be cast if necessary.

        Returns
        -------
        cumprod : ndarray
            A new array holding the result is returned unless out is specified,
            in which case a reference to out is returned.

        Notes
        -----
        The mask is lost if `out` is not a valid MaskedArray !

        Arithmetic is modular when using integer types, and no error is
        raised on overflow.

        """
        result = self.filled(1).cumprod(axis=axis, dtype=dtype, out=out)
        if out is not None:
            if isinstance(out, MaskedArray):
                out.__setmask__(self._mask)
            return out
        result = result.view(type(self))
        result.__setmask__(self._mask)
        return result


    def mean(self, axis=None, dtype=None, out=None):
        """
        Returns the average of the array elements.

        Masked entries are ignored.
        The average is taken over the flattened array by default, otherwise over
        the specified axis. Refer to `numpy.mean` for the full documentation.

        Parameters
        ----------
        a : array_like
            Array containing numbers whose mean is desired. If `a` is not an
            array, a conversion is attempted.
        axis : int, optional
            Axis along which the means are computed. The default is to compute
            the mean of the flattened array.
        dtype : dtype, optional
            Type to use in computing the mean. For integer inputs, the default
            is float64; for floating point, inputs it is the same as the input
            dtype.
        out : ndarray, optional
            Alternative output array in which to place the result. It must have
            the same shape as the expected output but the type will be cast if
            necessary.

        Returns
        -------
        mean : ndarray, see dtype parameter above
            If `out=None`, returns a new array containing the mean values,
            otherwise a reference to the output array is returned.

        See Also
        --------
        numpy.ma.mean : Equivalent function.
        numpy.mean : Equivalent function on non-masked arrays.
        numpy.ma.average: Weighted average.

        Examples
        --------
        >>> a = np.ma.array([1,2,3], mask=[False, False, True])
        >>> a
        masked_array(data = [1 2 --],
                     mask = [False False  True],
               fill_value = 999999)
        >>> a.mean()
        1.5

        """
        if self._mask is nomask:
            result = super(MaskedArray, self).mean(axis=axis, dtype=dtype)
        else:
            dsum = self.sum(axis=axis, dtype=dtype)
            cnt = self.count(axis=axis)
            if cnt.shape == () and (cnt == 0):
                result = masked
            else:
                result = dsum * 1. / cnt
        if out is not None:
            out.flat = result
            if isinstance(out, MaskedArray):
                outmask = getattr(out, '_mask', nomask)
                if (outmask is nomask):
                    outmask = out._mask = make_mask_none(out.shape)
                outmask.flat = getattr(result, '_mask', nomask)
            return out
        return result

    def anom(self, axis=None, dtype=None):
        """
        Compute the anomalies (deviations from the arithmetic mean)
        along the given axis.

        Returns an array of anomalies, with the same shape as the input and
        where the arithmetic mean is computed along the given axis.

        Parameters
        ----------
        axis : int, optional
            Axis over which the anomalies are taken.
            The default is to use the mean of the flattened array as reference.
        dtype : dtype, optional
            Type to use in computing the variance. For arrays of integer type
             the default is float32; for arrays of float types it is the same as
             the array type.

        See Also
        --------
        mean : Compute the mean of the array.

        Examples
        --------
        >>> a = np.ma.array([1,2,3])
        >>> a.anom()
        masked_array(data = [-1.  0.  1.],
                     mask = False,
               fill_value = 1e+20)

        """
        m = self.mean(axis, dtype)
        if not axis:
            return (self - m)
        else:
            return (self - expand_dims(m, axis))

    def var(self, axis=None, dtype=None, out=None, ddof=0):
        ""
        # Easy case: nomask, business as usual
        if self._mask is nomask:
            return self._data.var(axis=axis, dtype=dtype, out=out, ddof=ddof)
        # Some data are masked, yay!
        cnt = self.count(axis=axis) - ddof
        danom = self.anom(axis=axis, dtype=dtype)
        if iscomplexobj(self):
            danom = umath.absolute(danom) ** 2
        else:
            danom *= danom
        dvar = divide(danom.sum(axis), cnt).view(type(self))
        # Apply the mask if it's not a scalar
        if dvar.ndim:
            dvar._mask = mask_or(self._mask.all(axis), (cnt <= 0))
            dvar._update_from(self)
        elif getattr(dvar, '_mask', False):
        # Make sure that masked is returned when the scalar is masked.
            dvar = masked
            if out is not None:
                if isinstance(out, MaskedArray):
                    out.flat = 0
                    out.__setmask__(True)
                elif out.dtype.kind in 'biu':
                    errmsg = "Masked data information would be lost in one or "\
                             "more location."
                    raise MaskError(errmsg)
                else:
                    out.flat = np.nan
                return out
        # In case with have an explicit output
        if out is not None:
            # Set the data
            out.flat = dvar
            # Set the mask if needed
            if isinstance(out, MaskedArray):
                out.__setmask__(dvar.mask)
            return out
        return dvar
    var.__doc__ = np.var.__doc__


    def std(self, axis=None, dtype=None, out=None, ddof=0):
        ""
        dvar = self.var(axis=axis, dtype=dtype, out=out, ddof=ddof)
        if dvar is not masked:
            if out is not None:
                np.power(out, 0.5, out=out, casting='unsafe')
                return out
            dvar = sqrt(dvar)
        return dvar
    std.__doc__ = np.std.__doc__

    #............................................
    def round(self, decimals=0, out=None):
        """
        Return an array rounded a to the given number of decimals.

        Refer to `numpy.around` for full documentation.

        See Also
        --------
        numpy.around : equivalent function

        """
        result = self._data.round(decimals=decimals, out=out).view(type(self))
        result._mask = self._mask
        result._update_from(self)
        # No explicit output: we're done
        if out is None:
            return result
        if isinstance(out, MaskedArray):
            out.__setmask__(self._mask)
        return out
    round.__doc__ = ndarray.round.__doc__

    #............................................
    def argsort(self, axis=None, kind='quicksort', order=None, fill_value=None):
        """
        Return an ndarray of indices that sort the array along the
        specified axis.  Masked values are filled beforehand to
        `fill_value`.

        Parameters
        ----------
        axis : int, optional
            Axis along which to sort.  The default is -1 (last axis).
            If None, the flattened array is used.
        fill_value : var, optional
            Value used to fill the array before sorting.
            The default is the `fill_value` attribute of the input array.
        kind : {'quicksort', 'mergesort', 'heapsort'}, optional
            Sorting algorithm.
        order : list, optional
            When `a` is an array with fields defined, this argument specifies
            which fields to compare first, second, etc.  Not all fields need be
            specified.

        Returns
        -------
        index_array : ndarray, int
            Array of indices that sort `a` along the specified axis.
            In other words, ``a[index_array]`` yields a sorted `a`.

        See Also
        --------
        sort : Describes sorting algorithms used.
        lexsort : Indirect stable sort with multiple keys.
        ndarray.sort : Inplace sort.

        Notes
        -----
        See `sort` for notes on the different sorting algorithms.

        Examples
        --------
        >>> a = np.ma.array([3,2,1], mask=[False, False, True])
        >>> a
        masked_array(data = [3 2 --],
                     mask = [False False  True],
               fill_value = 999999)
        >>> a.argsort()
        array([1, 0, 2])

        """
        if fill_value is None:
            fill_value = default_fill_value(self)
        d = self.filled(fill_value).view(ndarray)
        return d.argsort(axis=axis, kind=kind, order=order)


    def argmin(self, axis=None, fill_value=None, out=None):
        """
        Return array of indices to the minimum values along the given axis.

        Parameters
        ----------
        axis : {None, integer}
            If None, the index is into the flattened array, otherwise along
            the specified axis
        fill_value : {var}, optional
            Value used to fill in the masked values.  If None, the output of
            minimum_fill_value(self._data) is used instead.
        out : {None, array}, optional
            Array into which the result can be placed. Its type is preserved
            and it must be of the right shape to hold the output.

        Returns
        -------
        {ndarray, scalar}
            If multi-dimension input, returns a new ndarray of indices to the
            minimum values along the given axis.  Otherwise, returns a scalar
            of index to the minimum values along the given axis.

        Examples
        --------
        >>> x = np.ma.array(arange(4), mask=[1,1,0,0])
        >>> x.shape = (2,2)
        >>> print x
        [[-- --]
         [2 3]]
        >>> print x.argmin(axis=0, fill_value=-1)
        [0 0]
        >>> print x.argmin(axis=0, fill_value=9)
        [1 1]

        """
        if fill_value is None:
            fill_value = minimum_fill_value(self)
        d = self.filled(fill_value).view(ndarray)
        return d.argmin(axis, out=out)


    def argmax(self, axis=None, fill_value=None, out=None):
        """
        Returns array of indices of the maximum values along the given axis.
        Masked values are treated as if they had the value fill_value.

        Parameters
        ----------
        axis : {None, integer}
            If None, the index is into the flattened array, otherwise along
            the specified axis
        fill_value : {var}, optional
            Value used to fill in the masked values.  If None, the output of
            maximum_fill_value(self._data) is used instead.
        out : {None, array}, optional
            Array into which the result can be placed. Its type is preserved
            and it must be of the right shape to hold the output.

        Returns
        -------
        index_array : {integer_array}

        Examples
        --------
        >>> a = np.arange(6).reshape(2,3)
        >>> a.argmax()
        5
        >>> a.argmax(0)
        array([1, 1, 1])
        >>> a.argmax(1)
        array([2, 2])

        """
        if fill_value is None:
            fill_value = maximum_fill_value(self._data)
        d = self.filled(fill_value).view(ndarray)
        return d.argmax(axis, out=out)


    def sort(self, axis= -1, kind='quicksort', order=None,
             endwith=True, fill_value=None):
        """
    Sort the array, in-place

    Parameters
    ----------
    a : array_like
        Array to be sorted.
    axis : int, optional
        Axis along which to sort. If None, the array is flattened before
        sorting. The default is -1, which sorts along the last axis.
    kind : {'quicksort', 'mergesort', 'heapsort'}, optional
        Sorting algorithm. Default is 'quicksort'.
    order : list, optional
        When `a` is a structured array, this argument specifies which fields
        to compare first, second, and so on.  This list does not need to
        include all of the fields.
    endwith : {True, False}, optional
        Whether missing values (if any) should be forced in the upper indices
        (at the end of the array) (True) or lower indices (at the beginning).
        When the array contains unmasked values of the largest (or smallest if
        False) representable value of the datatype the ordering of these values
        and the masked values is undefined.  To enforce the masked values are
        at the end (beginning) in this case one must sort the mask.
    fill_value : {var}, optional
        Value used internally for the masked values.
        If ``fill_value`` is not None, it supersedes ``endwith``.

    Returns
    -------
    sorted_array : ndarray
        Array of the same type and shape as `a`.

    See Also
    --------
    ndarray.sort : Method to sort an array in-place.
    argsort : Indirect sort.
    lexsort : Indirect stable sort on multiple keys.
    searchsorted : Find elements in a sorted array.

    Notes
    -----
    See ``sort`` for notes on the different sorting algorithms.

    Examples
    --------
    >>> a = ma.array([1, 2, 5, 4, 3],mask=[0, 1, 0, 1, 0])
    >>> # Default
    >>> a.sort()
    >>> print a
    [1 3 5 -- --]

    >>> a = ma.array([1, 2, 5, 4, 3],mask=[0, 1, 0, 1, 0])
    >>> # Put missing values in the front
    >>> a.sort(endwith=False)
    >>> print a
    [-- -- 1 3 5]

    >>> a = ma.array([1, 2, 5, 4, 3],mask=[0, 1, 0, 1, 0])
    >>> # fill_value takes over endwith
    >>> a.sort(endwith=False, fill_value=3)
    >>> print a
    [1 -- -- 3 5]

        """
        if self._mask is nomask:
            ndarray.sort(self, axis=axis, kind=kind, order=order)
        else:
            if self is masked:
                return self
            if fill_value is None:
                if endwith:
                    filler = minimum_fill_value(self)
                else:
                    filler = maximum_fill_value(self)
            else:
                filler = fill_value

            sidx = self.filled(filler).argsort(axis=axis, kind=kind,
                                               order=order)
            # save meshgrid memory for 1d arrays
            if self.ndim == 1:
                idx = sidx
            else:
                idx = np.meshgrid(*[np.arange(x) for x in self.shape], sparse=True,
                                  indexing='ij')
                idx[axis] = sidx
            tmp_mask = self._mask[idx].flat
            tmp_data = self._data[idx].flat
            self._data.flat = tmp_data
            self._mask.flat = tmp_mask
        return

    #............................................
    def min(self, axis=None, out=None, fill_value=None):
        """
    Return the minimum along a given axis.

    Parameters
    ----------
    axis : {None, int}, optional
        Axis along which to operate.  By default, ``axis`` is None and the
        flattened input is used.
    out : array_like, optional
        Alternative output array in which to place the result.  Must be of
        the same shape and buffer length as the expected output.
    fill_value : {var}, optional
        Value used to fill in the masked values.
        If None, use the output of `minimum_fill_value`.

    Returns
    -------
    amin : array_like
        New array holding the result.
        If ``out`` was specified, ``out`` is returned.

    See Also
    --------
    minimum_fill_value
        Returns the minimum filling value for a given datatype.

        """
        _mask = ndarray.__getattribute__(self, '_mask')
        newmask = _check_mask_axis(_mask, axis)
        if fill_value is None:
            fill_value = minimum_fill_value(self)
        # No explicit output
        if out is None:
            result = self.filled(fill_value).min(axis=axis, out=out).view(type(self))
            if result.ndim:
                # Set the mask
                result.__setmask__(newmask)
                # Get rid of Infs
                if newmask.ndim:
                    np.copyto(result, result.fill_value, where=newmask)
            elif newmask:
                result = masked
            return result
        # Explicit output
        result = self.filled(fill_value).min(axis=axis, out=out)
        if isinstance(out, MaskedArray):
            outmask = getattr(out, '_mask', nomask)
            if (outmask is nomask):
                outmask = out._mask = make_mask_none(out.shape)
            outmask.flat = newmask
        else:
            if out.dtype.kind in 'biu':
                errmsg = "Masked data information would be lost in one or more"\
                         " location."
                raise MaskError(errmsg)
            np.copyto(out, np.nan, where=newmask)
        return out

    def mini(self, axis=None):
        """
        Return the array minimum along the specified axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to find the minima. Default is None, in which case
            the minimum value in the whole array is returned.

        Returns
        -------
        min : scalar or MaskedArray
            If `axis` is None, the result is a scalar. Otherwise, if `axis` is
            given and the array is at least 2-D, the result is a masked array with
            dimension one smaller than the array on which `mini` is called.

        Examples
        --------
        >>> x = np.ma.array(np.arange(6), mask=[0 ,1, 0, 0, 0 ,1]).reshape(3, 2)
        >>> print x
        [[0 --]
         [2 3]
         [4 --]]
        >>> x.mini()
        0
        >>> x.mini(axis=0)
        masked_array(data = [0 3],
                     mask = [False False],
               fill_value = 999999)
        >>> print x.mini(axis=1)
        [0 2 4]

        """
        if axis is None:
            return minimum(self)
        else:
            return minimum.reduce(self, axis)

    #........................
    def max(self, axis=None, out=None, fill_value=None):
        """
        Return the maximum along a given axis.

        Parameters
        ----------
        axis : {None, int}, optional
            Axis along which to operate.  By default, ``axis`` is None and the
            flattened input is used.
        out : array_like, optional
            Alternative output array in which to place the result.  Must
            be of the same shape and buffer length as the expected output.
        fill_value : {var}, optional
            Value used to fill in the masked values.
            If None, use the output of maximum_fill_value().

        Returns
        -------
        amax : array_like
            New array holding the result.
            If ``out`` was specified, ``out`` is returned.

        See Also
        --------
        maximum_fill_value
            Returns the maximum filling value for a given datatype.

        """
        _mask = ndarray.__getattribute__(self, '_mask')
        newmask = _check_mask_axis(_mask, axis)
        if fill_value is None:
            fill_value = maximum_fill_value(self)
        # No explicit output
        if out is None:
            result = self.filled(fill_value).max(axis=axis, out=out).view(type(self))
            if result.ndim:
                # Set the mask
                result.__setmask__(newmask)
                # Get rid of Infs
                if newmask.ndim:
                    np.copyto(result, result.fill_value, where=newmask)
            elif newmask:
                result = masked
            return result
        # Explicit output
        result = self.filled(fill_value).max(axis=axis, out=out)
        if isinstance(out, MaskedArray):
            outmask = getattr(out, '_mask', nomask)
            if (outmask is nomask):
                outmask = out._mask = make_mask_none(out.shape)
            outmask.flat = newmask
        else:

            if out.dtype.kind in 'biu':
                errmsg = "Masked data information would be lost in one or more"\
                         " location."
                raise MaskError(errmsg)
            np.copyto(out, np.nan, where=newmask)
        return out

    def ptp(self, axis=None, out=None, fill_value=None):
        """
        Return (maximum - minimum) along the the given dimension
        (i.e. peak-to-peak value).

        Parameters
        ----------
        axis : {None, int}, optional
            Axis along which to find the peaks.  If None (default) the
            flattened array is used.
        out : {None, array_like}, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output
            but the type will be cast if necessary.
        fill_value : {var}, optional
            Value used to fill in the masked values.

        Returns
        -------
        ptp : ndarray.
            A new array holding the result, unless ``out`` was
            specified, in which case a reference to ``out`` is returned.

        """
        if out is None:
            result = self.max(axis=axis, fill_value=fill_value)
            result -= self.min(axis=axis, fill_value=fill_value)
            return result
        out.flat = self.max(axis=axis, out=out, fill_value=fill_value)
        min_value = self.min(axis=axis, fill_value=fill_value)
        np.subtract(out, min_value, out=out, casting='unsafe')
        return out

    def take(self, indices, axis=None, out=None, mode='raise'):
        """
        """
        (_data, _mask) = (self._data, self._mask)
        cls = type(self)
        # Make sure the indices are not masked
        maskindices = getattr(indices, '_mask', nomask)
        if maskindices is not nomask:
            indices = indices.filled(0)
        # Get the data
        if out is None:
            out = _data.take(indices, axis=axis, mode=mode).view(cls)
        else:
            np.take(_data, indices, axis=axis, mode=mode, out=out)
        # Get the mask
        if isinstance(out, MaskedArray):
            if _mask is nomask:
                outmask = maskindices
            else:
                outmask = _mask.take(indices, axis=axis, mode=mode)
                outmask |= maskindices
            out.__setmask__(outmask)
        return out


    # Array methods ---------------------------------------
    copy = _arraymethod('copy')
    diagonal = _arraymethod('diagonal')
    transpose = _arraymethod('transpose')
    T = property(fget=lambda self:self.transpose())
    swapaxes = _arraymethod('swapaxes')
    clip = _arraymethod('clip', onmask=False)
    copy = _arraymethod('copy')
    squeeze = _arraymethod('squeeze')
    #--------------------------------------------
    def tolist(self, fill_value=None):
        """
        Return the data portion of the masked array as a hierarchical Python list.

        Data items are converted to the nearest compatible Python type.
        Masked values are converted to `fill_value`. If `fill_value` is None,
        the corresponding entries in the output list will be ``None``.

        Parameters
        ----------
        fill_value : scalar, optional
            The value to use for invalid entries. Default is None.

        Returns
        -------
        result : list
            The Python list representation of the masked array.

        Examples
        --------
        >>> x = np.ma.array([[1,2,3], [4,5,6], [7,8,9]], mask=[0] + [1,0]*4)
        >>> x.tolist()
        [[1, None, 3], [None, 5, None], [7, None, 9]]
        >>> x.tolist(-999)
        [[1, -999, 3], [-999, 5, -999], [7, -999, 9]]

        """
        _mask = self._mask
        # No mask ? Just return .data.tolist ?
        if _mask is nomask:
            return self._data.tolist()
        # Explicit fill_value: fill the array and get the list
        if fill_value is not None:
            return self.filled(fill_value).tolist()
        # Structured array .............
        names = self.dtype.names
        if names:
            result = self._data.astype([(_, object) for _ in names])
            for n in names:
                result[n][_mask[n]] = None
            return result.tolist()
        # Standard arrays ...............
        if _mask is nomask:
            return [None]
        # Set temps to save time when dealing w/ marrays...
        inishape = self.shape
        result = np.array(self._data.ravel(), dtype=object)
        result[_mask.ravel()] = None
        result.shape = inishape
        return result.tolist()
#        if fill_value is not None:
#            return self.filled(fill_value).tolist()
#        result = self.filled().tolist()
#        # Set temps to save time when dealing w/ mrecarrays...
#        _mask = self._mask
#        if _mask is nomask:
#            return result
#        nbdims = self.ndim
#        dtypesize = len(self.dtype)
#        if nbdims == 0:
#            return tuple([None] * dtypesize)
#        elif nbdims == 1:
#            maskedidx = _mask.nonzero()[0].tolist()
#            if dtypesize:
#                nodata = tuple([None] * dtypesize)
#            else:
#                nodata = None
#            [operator.setitem(result, i, nodata) for i in maskedidx]
#        else:
#            for idx in zip(*[i.tolist() for i in _mask.nonzero()]):
#                tmp = result
#                for i in idx[:-1]:
#                    tmp = tmp[i]
#                tmp[idx[-1]] = None
#        return result
    #........................
    def tostring(self, fill_value=None, order='C'):
        """
        This function is a compatibility alias for tobytes. Despite its name it
        returns bytes not strings.
        """

        return self.tobytes(fill_value, order='C')
    #........................
    def tobytes(self, fill_value=None, order='C'):
        """
        Return the array data as a string containing the raw bytes in the array.

        The array is filled with a fill value before the string conversion.

        .. versionadded:: 1.9.0

        Parameters
        ----------
        fill_value : scalar, optional
            Value used to fill in the masked values. Deafult is None, in which
            case `MaskedArray.fill_value` is used.
        order : {'C','F','A'}, optional
            Order of the data item in the copy. Default is 'C'.

            - 'C'   -- C order (row major).
            - 'F'   -- Fortran order (column major).
            - 'A'   -- Any, current order of array.
            - None  -- Same as 'A'.

        See Also
        --------
        ndarray.tobytes
        tolist, tofile

        Notes
        -----
        As for `ndarray.tobytes`, information about the shape, dtype, etc.,
        but also about `fill_value`, will be lost.

        Examples
        --------
        >>> x = np.ma.array(np.array([[1, 2], [3, 4]]), mask=[[0, 1], [1, 0]])
        >>> x.tobytes()
        '\\x01\\x00\\x00\\x00?B\\x0f\\x00?B\\x0f\\x00\\x04\\x00\\x00\\x00'

        """
        return self.filled(fill_value).tobytes(order=order)
    #........................
    def tofile(self, fid, sep="", format="%s"):
        """
        Save a masked array to a file in binary format.

        .. warning::
          This function is not implemented yet.

        Raises
        ------
        NotImplementedError
            When `tofile` is called.

        """
        raise NotImplementedError("Not implemented yet, sorry...")

    def toflex(self):
        """
        Transforms a masked array into a flexible-type array.

        The flexible type array that is returned will have two fields:

        * the ``_data`` field stores the ``_data`` part of the array.
        * the ``_mask`` field stores the ``_mask`` part of the array.

        Parameters
        ----------
        None

        Returns
        -------
        record : ndarray
            A new flexible-type `ndarray` with two fields: the first element
            containing a value, the second element containing the corresponding
            mask boolean. The returned record shape matches self.shape.

        Notes
        -----
        A side-effect of transforming a masked array into a flexible `ndarray` is
        that meta information (``fill_value``, ...) will be lost.

        Examples
        --------
        >>> x = np.ma.array([[1,2,3],[4,5,6],[7,8,9]], mask=[0] + [1,0]*4)
        >>> print x
        [[1 -- 3]
         [-- 5 --]
         [7 -- 9]]
        >>> print x.toflex()
        [[(1, False) (2, True) (3, False)]
         [(4, True) (5, False) (6, True)]
         [(7, False) (8, True) (9, False)]]

        """
        # Get the basic dtype ....
        ddtype = self.dtype
        # Make sure we have a mask
        _mask = self._mask
        if _mask is None:
            _mask = make_mask_none(self.shape, ddtype)
        # And get its dtype
        mdtype = self._mask.dtype
        #
        record = np.ndarray(shape=self.shape,
                            dtype=[('_data', ddtype), ('_mask', mdtype)])
        record['_data'] = self._data
        record['_mask'] = self._mask
        return record
    torecords = toflex
    #--------------------------------------------
    # Pickling
    def __getstate__(self):
        """Return the internal state of the masked array, for pickling
        purposes.

        """
        cf = 'CF'[self.flags.fnc]
        state = (1,
                 self.shape,
                 self.dtype,
                 self.flags.fnc,
                 self._data.tobytes(cf),
                 #self._data.tolist(),
                 getmaskarray(self).tobytes(cf),
                 #getmaskarray(self).tolist(),
                 self._fill_value,
                 )
        return state
    #
    def __setstate__(self, state):
        """Restore the internal state of the masked array, for
        pickling purposes.  ``state`` is typically the output of the
        ``__getstate__`` output, and is a 5-tuple:

        - class name
        - a tuple giving the shape of the data
        - a typecode for the data
        - a binary string for the data
        - a binary string for the mask.

        """
        (_, shp, typ, isf, raw, msk, flv) = state
        ndarray.__setstate__(self, (shp, typ, isf, raw))
        self._mask.__setstate__((shp, make_mask_descr(typ), isf, msk))
        self.fill_value = flv
    #
    def __reduce__(self):
        """Return a 3-tuple for pickling a MaskedArray.

        """
        return (_mareconstruct,
                (self.__class__, self._baseclass, (0,), 'b',),
                self.__getstate__())
    #
    def __deepcopy__(self, memo=None):
        from copy import deepcopy
        copied = MaskedArray.__new__(type(self), self, copy=True)
        if memo is None:
            memo = {}
        memo[id(self)] = copied
        for (k, v) in self.__dict__.items():
            copied.__dict__[k] = deepcopy(v, memo)
        return copied


def _mareconstruct(subtype, baseclass, baseshape, basetype,):
    """Internal function that builds a new MaskedArray from the
    information stored in a pickle.

    """
    _data = ndarray.__new__(baseclass, baseshape, basetype)
    _mask = ndarray.__new__(ndarray, baseshape, make_mask_descr(basetype))
    return subtype.__new__(subtype, _data, mask=_mask, dtype=basetype,)






class mvoid(MaskedArray):
    """
    Fake a 'void' object to use for masked array with structured dtypes.
    """
    #
    def __new__(self, data, mask=nomask, dtype=None, fill_value=None,
                hardmask=False, copy=False, subok=True):
        _data = np.array(data, copy=copy, subok=subok, dtype=dtype)
        _data = _data.view(self)
        _data._hardmask = hardmask
        if mask is not nomask:
            if isinstance(mask, np.void):
                _data._mask = mask
            else:
                try:
                    # Mask is already a 0D array
                    _data._mask = np.void(mask)
                except TypeError:
                    # Transform the mask to a void
                    mdtype = make_mask_descr(dtype)
                    _data._mask = np.array(mask, dtype=mdtype)[()]
        if fill_value is not None:
            _data.fill_value = fill_value
        return _data

    def _get_data(self):
        # Make sure that the _data part is a np.void
        return self.view(ndarray)[()]
    _data = property(fget=_get_data)

    def __getitem__(self, indx):
        "Get the index..."
        m = self._mask
        if m is not nomask and m[indx]:
            return masked
        return self._data[indx]

    def __setitem__(self, indx, value):
        self._data[indx] = value
        if self._hardmask:
            self._mask[indx] |= getattr(value, "_mask", False)
        else:
            self._mask[indx] = getattr(value, "_mask", False)

    def __str__(self):
        m = self._mask
        if (m is nomask):
            return self._data.__str__()
        m = tuple(m)
        if (not any(m)):
            return self._data.__str__()
        r = self._data.tolist()
        p = masked_print_option
        if not p.enabled():
            p = 'N/A'
        else:
            p = str(p)
        r = [(str(_), p)[int(_m)] for (_, _m) in zip(r, m)]
        return "(%s)" % ", ".join(r)

    def __repr__(self):
        m = self._mask
        if (m is nomask):
            return self._data.__repr__()
        m = tuple(m)
        if not any(m):
            return self._data.__repr__()
        p = masked_print_option
        if not p.enabled():
            return self.filled(self.fill_value).__repr__()
        p = str(p)
        r = [(str(_), p)[int(_m)] for (_, _m) in zip(self._data.tolist(), m)]
        return "(%s)" % ", ".join(r)

    def __iter__(self):
        "Defines an iterator for mvoid"
        (_data, _mask) = (self._data, self._mask)
        if _mask is nomask:
            for d in _data:
                yield d
        else:
            for (d, m) in zip(_data, _mask):
                if m:
                    yield masked
                else:
                    yield d

    def __len__(self):
        return self._data.__len__()

    def filled(self, fill_value=None):
        """
        Return a copy with masked fields filled with a given value.

        Parameters
        ----------
        fill_value : scalar, optional
            The value to use for invalid entries (None by default).
            If None, the `fill_value` attribute is used instead.

        Returns
        -------
        filled_void
            A `np.void` object

        See Also
        --------
        MaskedArray.filled

    """
        return asarray(self).filled(fill_value)[()]

    def tolist(self):
        """
    Transforms the mvoid object into a tuple.

    Masked fields are replaced by None.

    Returns
    -------
    returned_tuple
        Tuple of fields
        """
        _mask = self._mask
        if _mask is nomask:
            return self._data.tolist()
        result = []
        for (d, m) in zip(self._data, self._mask):
            if m:
                result.append(None)
            else:
                # .item() makes sure we return a standard Python object
                result.append(d.item())
        return tuple(result)



#####--------------------------------------------------------------------------
#---- --- Shortcuts ---
#####---------------------------------------------------------------------------
def isMaskedArray(x):
    """
    Test whether input is an instance of MaskedArray.

    This function returns True if `x` is an instance of MaskedArray
    and returns False otherwise.  Any object is accepted as input.

    Parameters
    ----------
    x : object
        Object to test.

    Returns
    -------
    result : bool
        True if `x` is a MaskedArray.

    See Also
    --------
    isMA : Alias to isMaskedArray.
    isarray : Alias to isMaskedArray.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.eye(3, 3)
    >>> a
    array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])
    >>> m = ma.masked_values(a, 0)
    >>> m
    masked_array(data =
     [[1.0 -- --]
     [-- 1.0 --]
     [-- -- 1.0]],
          mask =
     [[False  True  True]
     [ True False  True]
     [ True  True False]],
          fill_value=0.0)
    >>> ma.isMaskedArray(a)
    False
    >>> ma.isMaskedArray(m)
    True
    >>> ma.isMaskedArray([0, 1, 2])
    False

    """
    return isinstance(x, MaskedArray)
isarray = isMaskedArray
isMA = isMaskedArray  #backward compatibility

# We define the masked singleton as a float for higher precedence...
# Note that it can be tricky sometimes w/ type comparison

class MaskedConstant(MaskedArray):
    #
    _data = data = np.array(0.)
    _mask = mask = np.array(True)
    _baseclass = ndarray
    #
    def __new__(self):
        return self._data.view(self)
    #
    def __array_finalize__(self, obj):
        return
    #
    def __array_wrap__(self, obj):
        return self
    #
    def __str__(self):
        return str(masked_print_option._display)
    #
    def __repr__(self):
        return 'masked'
    #
    def flatten(self):
        return masked_array([self._data], dtype=float, mask=[True])

    def __reduce__(self):
        """Override of MaskedArray's __reduce__.
        """
        return (self.__class__, ())


masked = masked_singleton = MaskedConstant()



masked_array = MaskedArray

def array(data, dtype=None, copy=False, order=False,
          mask=nomask, fill_value=None,
          keep_mask=True, hard_mask=False, shrink=True, subok=True, ndmin=0,
          ):
    """array(data, dtype=None, copy=False, order=False, mask=nomask,
             fill_value=None, keep_mask=True, hard_mask=False, shrink=True,
             subok=True, ndmin=0)

    Acts as shortcut to MaskedArray, with options in a different order
    for convenience.  And backwards compatibility...

    """
    #!!!: we should try to put 'order' somwehere
    return MaskedArray(data, mask=mask, dtype=dtype, copy=copy, subok=subok,
                       keep_mask=keep_mask, hard_mask=hard_mask,
                       fill_value=fill_value, ndmin=ndmin, shrink=shrink)
array.__doc__ = masked_array.__doc__

def is_masked(x):
    """
    Determine whether input has masked values.

    Accepts any object as input, but always returns False unless the
    input is a MaskedArray containing masked values.

    Parameters
    ----------
    x : array_like
        Array to check for masked values.

    Returns
    -------
    result : bool
        True if `x` is a MaskedArray with masked values, False otherwise.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> x = ma.masked_equal([0, 1, 0, 2, 3], 0)
    >>> x
    masked_array(data = [-- 1 -- 2 3],
          mask = [ True False  True False False],
          fill_value=999999)
    >>> ma.is_masked(x)
    True
    >>> x = ma.masked_equal([0, 1, 0, 2, 3], 42)
    >>> x
    masked_array(data = [0 1 0 2 3],
          mask = False,
          fill_value=999999)
    >>> ma.is_masked(x)
    False

    Always returns False if `x` isn't a MaskedArray.

    >>> x = [False, True, False]
    >>> ma.is_masked(x)
    False
    >>> x = 'a string'
    >>> ma.is_masked(x)
    False

    """
    m = getmask(x)
    if m is nomask:
        return False
    elif m.any():
        return True
    return False


#####---------------------------------------------------------------------------
#---- --- Extrema functions ---
#####---------------------------------------------------------------------------
class _extrema_operation(object):
    """
    Generic class for maximum/minimum functions.

    .. note::
      This is the base class for `_maximum_operation` and
      `_minimum_operation`.

    """
    def __call__(self, a, b=None):
        "Executes the call behavior."
        if b is None:
            return self.reduce(a)
        return where(self.compare(a, b), a, b)
    #.........
    def reduce(self, target, axis=None):
        "Reduce target along the given axis."
        target = narray(target, copy=False, subok=True)
        m = getmask(target)
        if axis is not None:
            kargs = { 'axis' : axis }
        else:
            kargs = {}
            target = target.ravel()
            if not (m is nomask):
                m = m.ravel()
        if m is nomask:
            t = self.ufunc.reduce(target, **kargs)
        else:
            target = target.filled(self.fill_value_func(target)).view(type(target))
            t = self.ufunc.reduce(target, **kargs)
            m = umath.logical_and.reduce(m, **kargs)
            if hasattr(t, '_mask'):
                t._mask = m
            elif m:
                t = masked
        return t
    #.........
    def outer (self, a, b):
        "Return the function applied to the outer product of a and b."
        ma = getmask(a)
        mb = getmask(b)
        if ma is nomask and mb is nomask:
            m = nomask
        else:
            ma = getmaskarray(a)
            mb = getmaskarray(b)
            m = logical_or.outer(ma, mb)
        result = self.ufunc.outer(filled(a), filled(b))
        if not isinstance(result, MaskedArray):
            result = result.view(MaskedArray)
        result._mask = m
        return result

#............................
class _minimum_operation(_extrema_operation):
    "Object to calculate minima"
    def __init__ (self):
        """minimum(a, b) or minimum(a)
In one argument case, returns the scalar minimum.
        """
        self.ufunc = umath.minimum
        self.afunc = amin
        self.compare = less
        self.fill_value_func = minimum_fill_value

#............................
class _maximum_operation(_extrema_operation):
    "Object to calculate maxima"
    def __init__ (self):
        """maximum(a, b) or maximum(a)
           In one argument case returns the scalar maximum.
        """
        self.ufunc = umath.maximum
        self.afunc = amax
        self.compare = greater
        self.fill_value_func = maximum_fill_value

#..........................................................
def min(obj, axis=None, out=None, fill_value=None):
    try:
        return obj.min(axis=axis, fill_value=fill_value, out=out)
    except (AttributeError, TypeError):
        # If obj doesn't have a min method,
        # ...or if the method doesn't accept a fill_value argument
        return asanyarray(obj).min(axis=axis, fill_value=fill_value, out=out)
min.__doc__ = MaskedArray.min.__doc__

def max(obj, axis=None, out=None, fill_value=None):
    try:
        return obj.max(axis=axis, fill_value=fill_value, out=out)
    except (AttributeError, TypeError):
        # If obj doesn't have a max method,
        # ...or if the method doesn't accept a fill_value argument
        return asanyarray(obj).max(axis=axis, fill_value=fill_value, out=out)
max.__doc__ = MaskedArray.max.__doc__

def ptp(obj, axis=None, out=None, fill_value=None):
    """a.ptp(axis=None) =  a.max(axis)-a.min(axis)"""
    try:
        return obj.ptp(axis, out=out, fill_value=fill_value)
    except (AttributeError, TypeError):
        # If obj doesn't have a ptp method,
        # ...or if the method doesn't accept a fill_value argument
        return asanyarray(obj).ptp(axis=axis, fill_value=fill_value, out=out)
ptp.__doc__ = MaskedArray.ptp.__doc__


#####---------------------------------------------------------------------------
#---- --- Definition of functions from the corresponding methods ---
#####---------------------------------------------------------------------------
class _frommethod:
    """
    Define functions from existing MaskedArray methods.

    Parameters
    ----------
    methodname : str
        Name of the method to transform.

    """
    def __init__(self, methodname, reversed=False):
        self.__name__ = methodname
        self.__doc__ = self.getdoc()
        self.reversed = reversed
    #
    def getdoc(self):
        "Return the doc of the function (from the doc of the method)."
        meth = getattr(MaskedArray, self.__name__, None) or\
               getattr(np, self.__name__, None)
        signature = self.__name__ + get_object_signature(meth)
        if meth is not None:
            doc = """    %s\n%s""" % (signature, getattr(meth, '__doc__', None))
            return doc
    #
    def __call__(self, a, *args, **params):
        if self.reversed:
            args = list(args)
            arr = args[0]
            args[0] = a
            a = arr
        # Get the method from the array (if possible)
        method_name = self.__name__
        method = getattr(a, method_name, None)
        if method is not None:
            return method(*args, **params)
        # Still here ? Then a is not a MaskedArray
        method = getattr(MaskedArray, method_name, None)
        if method is not None:
            return method(MaskedArray(a), *args, **params)
        # Still here ? OK, let's call the corresponding np function
        method = getattr(np, method_name)
        return method(a, *args, **params)

all = _frommethod('all')
anomalies = anom = _frommethod('anom')
any = _frommethod('any')
compress = _frommethod('compress', reversed=True)
cumprod = _frommethod('cumprod')
cumsum = _frommethod('cumsum')
copy = _frommethod('copy')
diagonal = _frommethod('diagonal')
harden_mask = _frommethod('harden_mask')
ids = _frommethod('ids')
maximum = _maximum_operation()
mean = _frommethod('mean')
minimum = _minimum_operation()
nonzero = _frommethod('nonzero')
prod = _frommethod('prod')
product = _frommethod('prod')
ravel = _frommethod('ravel')
repeat = _frommethod('repeat')
shrink_mask = _frommethod('shrink_mask')
soften_mask = _frommethod('soften_mask')
std = _frommethod('std')
sum = _frommethod('sum')
swapaxes = _frommethod('swapaxes')
#take = _frommethod('take')
trace = _frommethod('trace')
var = _frommethod('var')

def take(a, indices, axis=None, out=None, mode='raise'):
    """
    """
    a = masked_array(a)
    return a.take(indices, axis=axis, out=out, mode=mode)


#..............................................................................
def power(a, b, third=None):
    """
    Returns element-wise base array raised to power from second array.

    This is the masked array version of `numpy.power`. For details see
    `numpy.power`.

    See Also
    --------
    numpy.power

    Notes
    -----
    The *out* argument to `numpy.power` is not supported, `third` has to be
    None.

    """
    if third is not None:
        raise MaskError("3-argument power not supported.")
    # Get the masks
    ma = getmask(a)
    mb = getmask(b)
    m = mask_or(ma, mb)
    # Get the rawdata
    fa = getdata(a)
    fb = getdata(b)
    # Get the type of the result (so that we preserve subclasses)
    if isinstance(a, MaskedArray):
        basetype = type(a)
    else:
        basetype = MaskedArray
    # Get the result and view it as a (subclass of) MaskedArray
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(m, fa, umath.power(fa, fb)).view(basetype)
    result._update_from(a)
    # Find where we're in trouble w/ NaNs and Infs
    invalid = np.logical_not(np.isfinite(result.view(ndarray)))
    # Add the initial mask
    if m is not nomask:
        if not (result.ndim):
            return masked
        result._mask = np.logical_or(m, invalid)
    # Fix the invalid parts
    if invalid.any():
        if not result.ndim:
            return masked
        elif result._mask is nomask:
            result._mask = invalid
        result._data[invalid] = result.fill_value
    return result

#    if fb.dtype.char in typecodes["Integer"]:
#        return masked_array(umath.power(fa, fb), m)
#    m = mask_or(m, (fa < 0) & (fb != fb.astype(int)))
#    if m is nomask:
#        return masked_array(umath.power(fa, fb))
#    else:
#        fa = fa.copy()
#        if m.all():
#            fa.flat = 1
#        else:
#            np.copyto(fa, 1, where=m)
#        return masked_array(umath.power(fa, fb), m)

#..............................................................................
def argsort(a, axis=None, kind='quicksort', order=None, fill_value=None):
    "Function version of the eponymous method."
    if fill_value is None:
        fill_value = default_fill_value(a)
    d = filled(a, fill_value)
    if axis is None:
        return d.argsort(kind=kind, order=order)
    return d.argsort(axis, kind=kind, order=order)
argsort.__doc__ = MaskedArray.argsort.__doc__

def argmin(a, axis=None, fill_value=None):
    "Function version of the eponymous method."
    if fill_value is None:
        fill_value = default_fill_value(a)
    d = filled(a, fill_value)
    return d.argmin(axis=axis)
argmin.__doc__ = MaskedArray.argmin.__doc__

def argmax(a, axis=None, fill_value=None):
    "Function version of the eponymous method."
    if fill_value is None:
        fill_value = default_fill_value(a)
        try:
            fill_value = -fill_value
        except:
            pass
    d = filled(a, fill_value)
    return d.argmax(axis=axis)
argmax.__doc__ = MaskedArray.argmax.__doc__

def sort(a, axis= -1, kind='quicksort', order=None, endwith=True, fill_value=None):
    "Function version of the eponymous method."
    a = narray(a, copy=True, subok=True)
    if axis is None:
        a = a.flatten()
        axis = 0
    if fill_value is None:
        if endwith:
            filler = minimum_fill_value(a)
        else:
            filler = maximum_fill_value(a)
    else:
        filler = fill_value

    sindx = filled(a, filler).argsort(axis=axis, kind=kind, order=order)

    # save meshgrid memory for 1d arrays
    if a.ndim == 1:
        indx = sindx
    else:
        indx = np.meshgrid(*[np.arange(x) for x in a.shape], sparse=True,
                           indexing='ij')
        indx[axis] = sindx
    return a[indx]
sort.__doc__ = MaskedArray.sort.__doc__


def compressed(x):
    """
    Return all the non-masked data as a 1-D array.

    This function is equivalent to calling the "compressed" method of a
    `MaskedArray`, see `MaskedArray.compressed` for details.

    See Also
    --------
    MaskedArray.compressed
        Equivalent method.

    """
    if not isinstance(x, MaskedArray):
        x = asanyarray(x)
    return x.compressed()


def concatenate(arrays, axis=0):
    """
    Concatenate a sequence of arrays along the given axis.

    Parameters
    ----------
    arrays : sequence of array_like
        The arrays must have the same shape, except in the dimension
        corresponding to `axis` (the first, by default).
    axis : int, optional
        The axis along which the arrays will be joined. Default is 0.

    Returns
    -------
    result : MaskedArray
        The concatenated array with any masked entries preserved.

    See Also
    --------
    numpy.concatenate : Equivalent function in the top-level NumPy module.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = ma.arange(3)
    >>> a[1] = ma.masked
    >>> b = ma.arange(2, 5)
    >>> a
    masked_array(data = [0 -- 2],
                 mask = [False  True False],
           fill_value = 999999)
    >>> b
    masked_array(data = [2 3 4],
                 mask = False,
           fill_value = 999999)
    >>> ma.concatenate([a, b])
    masked_array(data = [0 -- 2 2 3 4],
                 mask = [False  True False False False False],
           fill_value = 999999)

    """
    d = np.concatenate([getdata(a) for a in arrays], axis)
    rcls = get_masked_subclass(*arrays)
    data = d.view(rcls)
    # Check whether one of the arrays has a non-empty mask...
    for x in arrays:
        if getmask(x) is not nomask:
            break
    else:
        return data
    # OK, so we have to concatenate the masks
    dm = np.concatenate([getmaskarray(a) for a in arrays], axis)
    # If we decide to keep a '_shrinkmask' option, we want to check that ...
    # ... all of them are True, and then check for dm.any()
#    shrink = numpy.logical_or.reduce([getattr(a,'_shrinkmask',True) for a in arrays])
#    if shrink and not dm.any():
    if not dm.dtype.fields and not dm.any():
        data._mask = nomask
    else:
        data._mask = dm.reshape(d.shape)
    return data

def count(a, axis=None):
    if isinstance(a, MaskedArray):
        return a.count(axis)
    return masked_array(a, copy=False).count(axis)
count.__doc__ = MaskedArray.count.__doc__


def diag(v, k=0):
    """
    Extract a diagonal or construct a diagonal array.

    This function is the equivalent of `numpy.diag` that takes masked
    values into account, see `numpy.diag` for details.

    See Also
    --------
    numpy.diag : Equivalent function for ndarrays.

    """
    output = np.diag(v, k).view(MaskedArray)
    if getmask(v) is not nomask:
        output._mask = np.diag(v._mask, k)
    return output


def expand_dims(x, axis):
    """
    Expand the shape of an array.

    Expands the shape of the array by including a new axis before the one
    specified by the `axis` parameter. This function behaves the same as
    `numpy.expand_dims` but preserves masked elements.

    See Also
    --------
    numpy.expand_dims : Equivalent function in top-level NumPy module.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> x = ma.array([1, 2, 4])
    >>> x[1] = ma.masked
    >>> x
    masked_array(data = [1 -- 4],
                 mask = [False  True False],
           fill_value = 999999)
    >>> np.expand_dims(x, axis=0)
    array([[1, 2, 4]])
    >>> ma.expand_dims(x, axis=0)
    masked_array(data =
     [[1 -- 4]],
                 mask =
     [[False  True False]],
           fill_value = 999999)

    The same result can be achieved using slicing syntax with `np.newaxis`.

    >>> x[np.newaxis, :]
    masked_array(data =
     [[1 -- 4]],
                 mask =
     [[False  True False]],
           fill_value = 999999)

    """
    result = n_expand_dims(x, axis)
    if isinstance(x, MaskedArray):
        new_shape = result.shape
        result = x.view()
        result.shape = new_shape
        if result._mask is not nomask:
            result._mask.shape = new_shape
    return result

#......................................
def left_shift (a, n):
    """
    Shift the bits of an integer to the left.

    This is the masked array version of `numpy.left_shift`, for details
    see that function.

    See Also
    --------
    numpy.left_shift

    """
    m = getmask(a)
    if m is nomask:
        d = umath.left_shift(filled(a), n)
        return masked_array(d)
    else:
        d = umath.left_shift(filled(a, 0), n)
        return masked_array(d, mask=m)

def right_shift (a, n):
    """
    Shift the bits of an integer to the right.

    This is the masked array version of `numpy.right_shift`, for details
    see that function.

    See Also
    --------
    numpy.right_shift

    """
    m = getmask(a)
    if m is nomask:
        d = umath.right_shift(filled(a), n)
        return masked_array(d)
    else:
        d = umath.right_shift(filled(a, 0), n)
        return masked_array(d, mask=m)

#......................................
def put(a, indices, values, mode='raise'):
    """
    Set storage-indexed locations to corresponding values.

    This function is equivalent to `MaskedArray.put`, see that method
    for details.

    See Also
    --------
    MaskedArray.put

    """
    # We can't use 'frommethod', the order of arguments is different
    try:
        return a.put(indices, values, mode=mode)
    except AttributeError:
        return narray(a, copy=False).put(indices, values, mode=mode)

def putmask(a, mask, values): #, mode='raise'):
    """
    Changes elements of an array based on conditional and input values.

    This is the masked array version of `numpy.putmask`, for details see
    `numpy.putmask`.

    See Also
    --------
    numpy.putmask

    Notes
    -----
    Using a masked array as `values` will **not** transform a `ndarray` into
    a `MaskedArray`.

    """
    # We can't use 'frommethod', the order of arguments is different
    if not isinstance(a, MaskedArray):
        a = a.view(MaskedArray)
    (valdata, valmask) = (getdata(values), getmask(values))
    if getmask(a) is nomask:
        if valmask is not nomask:
            a._sharedmask = True
            a._mask = make_mask_none(a.shape, a.dtype)
            np.copyto(a._mask, valmask, where=mask)
    elif a._hardmask:
        if valmask is not nomask:
            m = a._mask.copy()
            np.copyto(m, valmask, where=mask)
            a.mask |= m
    else:
        if valmask is nomask:
            valmask = getmaskarray(values)
        np.copyto(a._mask, valmask, where=mask)
    np.copyto(a._data, valdata, where=mask)
    return

def transpose(a, axes=None):
    """
    Permute the dimensions of an array.

    This function is exactly equivalent to `numpy.transpose`.

    See Also
    --------
    numpy.transpose : Equivalent function in top-level NumPy module.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> x = ma.arange(4).reshape((2,2))
    >>> x[1, 1] = ma.masked
    >>>> x
    masked_array(data =
     [[0 1]
     [2 --]],
                 mask =
     [[False False]
     [False  True]],
           fill_value = 999999)
    >>> ma.transpose(x)
    masked_array(data =
     [[0 2]
     [1 --]],
                 mask =
     [[False False]
     [False  True]],
           fill_value = 999999)

    """
    #We can't use 'frommethod', as 'transpose' doesn't take keywords
    try:
        return a.transpose(axes)
    except AttributeError:
        return narray(a, copy=False).transpose(axes).view(MaskedArray)

def reshape(a, new_shape, order='C'):
    """
    Returns an array containing the same data with a new shape.

    Refer to `MaskedArray.reshape` for full documentation.

    See Also
    --------
    MaskedArray.reshape : equivalent function

    """
    #We can't use 'frommethod', it whine about some parameters. Dmmit.
    try:
        return a.reshape(new_shape, order=order)
    except AttributeError:
        _tmp = narray(a, copy=False).reshape(new_shape, order=order)
        return _tmp.view(MaskedArray)

def resize(x, new_shape):
    """
    Return a new masked array with the specified size and shape.

    This is the masked equivalent of the `numpy.resize` function. The new
    array is filled with repeated copies of `x` (in the order that the
    data are stored in memory). If `x` is masked, the new array will be
    masked, and the new mask will be a repetition of the old one.

    See Also
    --------
    numpy.resize : Equivalent function in the top level NumPy module.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = ma.array([[1, 2] ,[3, 4]])
    >>> a[0, 1] = ma.masked
    >>> a
    masked_array(data =
     [[1 --]
     [3 4]],
                 mask =
     [[False  True]
     [False False]],
           fill_value = 999999)
    >>> np.resize(a, (3, 3))
    array([[1, 2, 3],
           [4, 1, 2],
           [3, 4, 1]])
    >>> ma.resize(a, (3, 3))
    masked_array(data =
     [[1 -- 3]
     [4 1 --]
     [3 4 1]],
                 mask =
     [[False  True False]
     [False False  True]
     [False False False]],
           fill_value = 999999)

    A MaskedArray is always returned, regardless of the input type.

    >>> a = np.array([[1, 2] ,[3, 4]])
    >>> ma.resize(a, (3, 3))
    masked_array(data =
     [[1 2 3]
     [4 1 2]
     [3 4 1]],
                 mask =
     False,
           fill_value = 999999)

    """
    # We can't use _frommethods here, as N.resize is notoriously whiny.
    m = getmask(x)
    if m is not nomask:
        m = np.resize(m, new_shape)
    result = np.resize(x, new_shape).view(get_masked_subclass(x))
    if result.ndim:
        result._mask = m
    return result


#................................................
def rank(obj):
    "maskedarray version of the numpy function."
    return np.rank(getdata(obj))
rank.__doc__ = np.rank.__doc__
#
def shape(obj):
    "maskedarray version of the numpy function."
    return np.shape(getdata(obj))
shape.__doc__ = np.shape.__doc__
#
def size(obj, axis=None):
    "maskedarray version of the numpy function."
    return np.size(getdata(obj), axis)
size.__doc__ = np.size.__doc__
#................................................

#####--------------------------------------------------------------------------
#---- --- Extra functions ---
#####--------------------------------------------------------------------------
def where (condition, x=None, y=None):
    """
    Return a masked array with elements from x or y, depending on condition.

    Returns a masked array, shaped like condition, where the elements
    are from `x` when `condition` is True, and from `y` otherwise.
    If neither `x` nor `y` are given, the function returns a tuple of
    indices where `condition` is True (the result of
    ``condition.nonzero()``).

    Parameters
    ----------
    condition : array_like, bool
        The condition to meet. For each True element, yield the corresponding
        element from `x`, otherwise from `y`.
    x, y : array_like, optional
        Values from which to choose. `x` and `y` need to have the same shape
        as condition, or be broadcast-able to that shape.

    Returns
    -------
    out : MaskedArray or tuple of ndarrays
        The resulting masked array if `x` and `y` were given, otherwise
        the result of ``condition.nonzero()``.

    See Also
    --------
    numpy.where : Equivalent function in the top-level NumPy module.

    Examples
    --------
    >>> x = np.ma.array(np.arange(9.).reshape(3, 3), mask=[[0, 1, 0],
    ...                                                    [1, 0, 1],
    ...                                                    [0, 1, 0]])
    >>> print x
    [[0.0 -- 2.0]
     [-- 4.0 --]
     [6.0 -- 8.0]]
    >>> np.ma.where(x > 5)    # return the indices where x > 5
    (array([2, 2]), array([0, 2]))

    >>> print np.ma.where(x > 5, x, -3.1416)
    [[-3.1416 -- -3.1416]
     [-- -3.1416 --]
     [6.0 -- 8.0]]

    """
    if x is None and y is None:
        return filled(condition, 0).nonzero()
    elif x is None or y is None:
        raise ValueError("Either both or neither x and y should be given.")
    # Get the condition ...............
    fc = filled(condition, 0).astype(MaskType)
    notfc = np.logical_not(fc)
    # Get the data ......................................
    xv = getdata(x)
    yv = getdata(y)
    if x is masked:
        ndtype = yv.dtype
    elif y is masked:
        ndtype = xv.dtype
    else:
        ndtype = np.find_common_type([xv.dtype, yv.dtype], [])
    # Construct an empty array and fill it
    d = np.empty(fc.shape, dtype=ndtype).view(MaskedArray)
    _data = d._data
    np.copyto(_data, xv.astype(ndtype), where=fc)
    np.copyto(_data, yv.astype(ndtype), where=notfc)
    # Create an empty mask and fill it
    _mask = d._mask = np.zeros(fc.shape, dtype=MaskType)
    np.copyto(_mask, getmask(x), where=fc)
    np.copyto(_mask, getmask(y), where=notfc)
    _mask |= getmaskarray(condition)
    if not _mask.any():
        d._mask = nomask
    return d

def choose (indices, choices, out=None, mode='raise'):
    """
    Use an index array to construct a new array from a set of choices.

    Given an array of integers and a set of n choice arrays, this method
    will create a new array that merges each of the choice arrays.  Where a
    value in `a` is i, the new array will have the value that choices[i]
    contains in the same place.

    Parameters
    ----------
    a : ndarray of ints
        This array must contain integers in ``[0, n-1]``, where n is the
        number of choices.
    choices : sequence of arrays
        Choice arrays. The index array and all of the choices should be
        broadcastable to the same shape.
    out : array, optional
        If provided, the result will be inserted into this array. It should
        be of the appropriate shape and `dtype`.
    mode : {'raise', 'wrap', 'clip'}, optional
        Specifies how out-of-bounds indices will behave.

        * 'raise' : raise an error
        * 'wrap' : wrap around
        * 'clip' : clip to the range

    Returns
    -------
    merged_array : array

    See Also
    --------
    choose : equivalent function

    Examples
    --------
    >>> choice = np.array([[1,1,1], [2,2,2], [3,3,3]])
    >>> a = np.array([2, 1, 0])
    >>> np.ma.choose(a, choice)
    masked_array(data = [3 2 1],
          mask = False,
          fill_value=999999)

    """
    def fmask (x):
        "Returns the filled array, or True if masked."
        if x is masked:
            return True
        return filled(x)
    def nmask (x):
        "Returns the mask, True if ``masked``, False if ``nomask``."
        if x is masked:
            return True
        return getmask(x)
    # Get the indices......
    c = filled(indices, 0)
    # Get the masks........
    masks = [nmask(x) for x in choices]
    data = [fmask(x) for x in choices]
    # Construct the mask
    outputmask = np.choose(c, masks, mode=mode)
    outputmask = make_mask(mask_or(outputmask, getmask(indices)),
                           copy=0, shrink=True)
    # Get the choices......
    d = np.choose(c, data, mode=mode, out=out).view(MaskedArray)
    if out is not None:
        if isinstance(out, MaskedArray):
            out.__setmask__(outputmask)
        return out
    d.__setmask__(outputmask)
    return d


def round_(a, decimals=0, out=None):
    """
    Return a copy of a, rounded to 'decimals' places.

    When 'decimals' is negative, it specifies the number of positions
    to the left of the decimal point.  The real and imaginary parts of
    complex numbers are rounded separately. Nothing is done if the
    array is not of float type and 'decimals' is greater than or equal
    to 0.

    Parameters
    ----------
    decimals : int
        Number of decimals to round to. May be negative.
    out : array_like
        Existing array to use for output.
        If not given, returns a default copy of a.

    Notes
    -----
    If out is given and does not have a mask attribute, the mask of a
    is lost!

    """
    if out is None:
        return np.round_(a, decimals, out)
    else:
        np.round_(getdata(a), decimals, out)
        if hasattr(out, '_mask'):
            out._mask = getmask(a)
        return out
round = round_

def inner(a, b):
    """
    Returns the inner product of a and b for arrays of floating point types.

    Like the generic NumPy equivalent the product sum is over the last dimension
    of a and b.

    Notes
    -----
    The first argument is not conjugated.

    """
    fa = filled(a, 0)
    fb = filled(b, 0)
    if len(fa.shape) == 0:
        fa.shape = (1,)
    if len(fb.shape) == 0:
        fb.shape = (1,)
    return np.inner(fa, fb).view(MaskedArray)
inner.__doc__ = doc_note(np.inner.__doc__,
                         "Masked values are replaced by 0.")
innerproduct = inner

def outer(a, b):
    "maskedarray version of the numpy function."
    fa = filled(a, 0).ravel()
    fb = filled(b, 0).ravel()
    d = np.outer(fa, fb)
    ma = getmask(a)
    mb = getmask(b)
    if ma is nomask and mb is nomask:
        return masked_array(d)
    ma = getmaskarray(a)
    mb = getmaskarray(b)
    m = make_mask(1 - np.outer(1 - ma, 1 - mb), copy=0)
    return masked_array(d, mask=m)
outer.__doc__ = doc_note(np.outer.__doc__,
                         "Masked values are replaced by 0.")
outerproduct = outer

def allequal (a, b, fill_value=True):
    """
    Return True if all entries of a and b are equal, using
    fill_value as a truth value where either or both are masked.

    Parameters
    ----------
    a, b : array_like
        Input arrays to compare.
    fill_value : bool, optional
        Whether masked values in a or b are considered equal (True) or not
        (False).

    Returns
    -------
    y : bool
        Returns True if the two arrays are equal within the given
        tolerance, False otherwise. If either array contains NaN,
        then False is returned.

    See Also
    --------
    all, any
    numpy.ma.allclose

    Examples
    --------
    >>> a = ma.array([1e10, 1e-7, 42.0], mask=[0, 0, 1])
    >>> a
    masked_array(data = [10000000000.0 1e-07 --],
          mask = [False False  True],
          fill_value=1e+20)

    >>> b = array([1e10, 1e-7, -42.0])
    >>> b
    array([  1.00000000e+10,   1.00000000e-07,  -4.20000000e+01])
    >>> ma.allequal(a, b, fill_value=False)
    False
    >>> ma.allequal(a, b)
    True

    """
    m = mask_or(getmask(a), getmask(b))
    if m is nomask:
        x = getdata(a)
        y = getdata(b)
        d = umath.equal(x, y)
        return d.all()
    elif fill_value:
        x = getdata(a)
        y = getdata(b)
        d = umath.equal(x, y)
        dm = array(d, mask=m, copy=False)
        return dm.filled(True).all(None)
    else:
        return False

def allclose (a, b, masked_equal=True, rtol=1e-5, atol=1e-8):
    """
    Returns True if two arrays are element-wise equal within a tolerance.

    This function is equivalent to `allclose` except that masked values
    are treated as equal (default) or unequal, depending on the `masked_equal`
    argument.

    Parameters
    ----------
    a, b : array_like
        Input arrays to compare.
    masked_equal : bool, optional
        Whether masked values in `a` and `b` are considered equal (True) or not
        (False). They are considered equal by default.
    rtol : float, optional
        Relative tolerance. The relative difference is equal to ``rtol * b``.
        Default is 1e-5.
    atol : float, optional
        Absolute tolerance. The absolute difference is equal to `atol`.
        Default is 1e-8.

    Returns
    -------
    y : bool
        Returns True if the two arrays are equal within the given
        tolerance, False otherwise. If either array contains NaN, then
        False is returned.

    See Also
    --------
    all, any
    numpy.allclose : the non-masked `allclose`.

    Notes
    -----
    If the following equation is element-wise True, then `allclose` returns
    True::

      absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    Return True if all elements of `a` and `b` are equal subject to
    given tolerances.

    Examples
    --------
    >>> a = ma.array([1e10, 1e-7, 42.0], mask=[0, 0, 1])
    >>> a
    masked_array(data = [10000000000.0 1e-07 --],
                 mask = [False False  True],
           fill_value = 1e+20)
    >>> b = ma.array([1e10, 1e-8, -42.0], mask=[0, 0, 1])
    >>> ma.allclose(a, b)
    False

    >>> a = ma.array([1e10, 1e-8, 42.0], mask=[0, 0, 1])
    >>> b = ma.array([1.00001e10, 1e-9, -42.0], mask=[0, 0, 1])
    >>> ma.allclose(a, b)
    True
    >>> ma.allclose(a, b, masked_equal=False)
    False

    Masked values are not compared directly.

    >>> a = ma.array([1e10, 1e-8, 42.0], mask=[0, 0, 1])
    >>> b = ma.array([1.00001e10, 1e-9, 42.0], mask=[0, 0, 1])
    >>> ma.allclose(a, b)
    True
    >>> ma.allclose(a, b, masked_equal=False)
    False

    """
    x = masked_array(a, copy=False)
    y = masked_array(b, copy=False)

    # make sure y is an inexact type to avoid abs(MIN_INT); will cause
    # casting of x later.
    dtype = np.result_type(y, 1.)
    if y.dtype != dtype:
        y = masked_array(y, dtype=dtype, copy=False)

    m = mask_or(getmask(x), getmask(y))
    xinf = np.isinf(masked_array(x, copy=False, mask=m)).filled(False)
    # If we have some infs, they should fall at the same place.
    if not np.all(xinf == filled(np.isinf(y), False)):
        return False
    # No infs at all
    if not np.any(xinf):
        d = filled(umath.less_equal(umath.absolute(x - y),
                                    atol + rtol * umath.absolute(y)),
                   masked_equal)
        return np.all(d)

    if not np.all(filled(x[xinf] == y[xinf], masked_equal)):
        return False
    x = x[~xinf]
    y = y[~xinf]

    d = filled(umath.less_equal(umath.absolute(x - y),
                                atol + rtol * umath.absolute(y)),
               masked_equal)

    return np.all(d)

#..............................................................................
def asarray(a, dtype=None, order=None):
    """
    Convert the input to a masked array of the given data-type.

    No copy is performed if the input is already an `ndarray`. If `a` is
    a subclass of `MaskedArray`, a base class `MaskedArray` is returned.

    Parameters
    ----------
    a : array_like
        Input data, in any form that can be converted to a masked array. This
        includes lists, lists of tuples, tuples, tuples of tuples, tuples
        of lists, ndarrays and masked arrays.
    dtype : dtype, optional
        By default, the data-type is inferred from the input data.
    order : {'C', 'F'}, optional
        Whether to use row-major ('C') or column-major ('FORTRAN') memory
        representation.  Default is 'C'.

    Returns
    -------
    out : MaskedArray
        Masked array interpretation of `a`.

    See Also
    --------
    asanyarray : Similar to `asarray`, but conserves subclasses.

    Examples
    --------
    >>> x = np.arange(10.).reshape(2, 5)
    >>> x
    array([[ 0.,  1.,  2.,  3.,  4.],
           [ 5.,  6.,  7.,  8.,  9.]])
    >>> np.ma.asarray(x)
    masked_array(data =
     [[ 0.  1.  2.  3.  4.]
     [ 5.  6.  7.  8.  9.]],
                 mask =
     False,
           fill_value = 1e+20)
    >>> type(np.ma.asarray(x))
    <class 'numpy.ma.core.MaskedArray'>

    """
    return masked_array(a, dtype=dtype, copy=False, keep_mask=True, subok=False)

def asanyarray(a, dtype=None):
    """
    Convert the input to a masked array, conserving subclasses.

    If `a` is a subclass of `MaskedArray`, its class is conserved.
    No copy is performed if the input is already an `ndarray`.

    Parameters
    ----------
    a : array_like
        Input data, in any form that can be converted to an array.
    dtype : dtype, optional
        By default, the data-type is inferred from the input data.
    order : {'C', 'F'}, optional
        Whether to use row-major ('C') or column-major ('FORTRAN') memory
        representation.  Default is 'C'.

    Returns
    -------
    out : MaskedArray
        MaskedArray interpretation of `a`.

    See Also
    --------
    asarray : Similar to `asanyarray`, but does not conserve subclass.

    Examples
    --------
    >>> x = np.arange(10.).reshape(2, 5)
    >>> x
    array([[ 0.,  1.,  2.,  3.,  4.],
           [ 5.,  6.,  7.,  8.,  9.]])
    >>> np.ma.asanyarray(x)
    masked_array(data =
     [[ 0.  1.  2.  3.  4.]
     [ 5.  6.  7.  8.  9.]],
                 mask =
     False,
           fill_value = 1e+20)
    >>> type(np.ma.asanyarray(x))
    <class 'numpy.ma.core.MaskedArray'>

    """
    return masked_array(a, dtype=dtype, copy=False, keep_mask=True, subok=True)


#####--------------------------------------------------------------------------
#---- --- Pickling ---
#####--------------------------------------------------------------------------
def dump(a, F):
    """
    Pickle a masked array to a file.

    This is a wrapper around ``cPickle.dump``.

    Parameters
    ----------
    a : MaskedArray
        The array to be pickled.
    F : str or file-like object
        The file to pickle `a` to. If a string, the full path to the file.

    """
    if not hasattr(F, 'readline'):
        F = open(F, 'w')
    return pickle.dump(a, F)

def dumps(a):
    """
    Return a string corresponding to the pickling of a masked array.

    This is a wrapper around ``cPickle.dumps``.

    Parameters
    ----------
    a : MaskedArray
        The array for which the string representation of the pickle is
        returned.

    """
    return pickle.dumps(a)

def load(F):
    """
    Wrapper around ``cPickle.load`` which accepts either a file-like object
    or a filename.

    Parameters
    ----------
    F : str or file
        The file or file name to load.

    See Also
    --------
    dump : Pickle an array

    Notes
    -----
    This is different from `numpy.load`, which does not use cPickle but loads
    the NumPy binary .npy format.

    """
    if not hasattr(F, 'readline'):
        F = open(F, 'r')
    return pickle.load(F)

def loads(strg):
    """
    Load a pickle from the current string.

    The result of ``cPickle.loads(strg)`` is returned.

    Parameters
    ----------
    strg : str
        The string to load.

    See Also
    --------
    dumps : Return a string corresponding to the pickling of a masked array.

    """
    return pickle.loads(strg)

################################################################################
def fromfile(file, dtype=float, count= -1, sep=''):
    raise NotImplementedError("Not yet implemented. Sorry")


def fromflex(fxarray):
    """
    Build a masked array from a suitable flexible-type array.

    The input array has to have a data-type with ``_data`` and ``_mask``
    fields. This type of array is output by `MaskedArray.toflex`.

    Parameters
    ----------
    fxarray : ndarray
        The structured input array, containing ``_data`` and ``_mask``
        fields. If present, other fields are discarded.

    Returns
    -------
    result : MaskedArray
        The constructed masked array.

    See Also
    --------
    MaskedArray.toflex : Build a flexible-type array from a masked array.

    Examples
    --------
    >>> x = np.ma.array(np.arange(9).reshape(3, 3), mask=[0] + [1, 0] * 4)
    >>> rec = x.toflex()
    >>> rec
    array([[(0, False), (1, True), (2, False)],
           [(3, True), (4, False), (5, True)],
           [(6, False), (7, True), (8, False)]],
          dtype=[('_data', '<i4'), ('_mask', '|b1')])
    >>> x2 = np.ma.fromflex(rec)
    >>> x2
    masked_array(data =
     [[0 -- 2]
     [-- 4 --]
     [6 -- 8]],
                 mask =
     [[False  True False]
     [ True False  True]
     [False  True False]],
           fill_value = 999999)

    Extra fields can be present in the structured array but are discarded:

    >>> dt = [('_data', '<i4'), ('_mask', '|b1'), ('field3', '<f4')]
    >>> rec2 = np.zeros((2, 2), dtype=dt)
    >>> rec2
    array([[(0, False, 0.0), (0, False, 0.0)],
           [(0, False, 0.0), (0, False, 0.0)]],
          dtype=[('_data', '<i4'), ('_mask', '|b1'), ('field3', '<f4')])
    >>> y = np.ma.fromflex(rec2)
    >>> y
    masked_array(data =
     [[0 0]
     [0 0]],
                 mask =
     [[False False]
     [False False]],
           fill_value = 999999)

    """
    return masked_array(fxarray['_data'], mask=fxarray['_mask'])



class _convert2ma:
    """
    Convert functions from numpy to numpy.ma.

    Parameters
    ----------
        _methodname : string
            Name of the method to transform.

    """
    __doc__ = None
    #
    def __init__(self, funcname, params=None):
        self._func = getattr(np, funcname)
        self.__doc__ = self.getdoc()
        self._extras = params or {}
    #
    def getdoc(self):
        "Return the doc of the function (from the doc of the method)."
        doc = getattr(self._func, '__doc__', None)
        sig = get_object_signature(self._func)
        if doc:
            # Add the signature of the function at the beginning of the doc
            if sig:
                sig = "%s%s\n" % (self._func.__name__, sig)
            doc = sig + doc
        return doc
    #
    def __call__(self, a, *args, **params):
        # Find the common parameters to the call and the definition
        _extras = self._extras
        common_params = set(params).intersection(_extras)
        # Drop the common parameters from the call
        for p in common_params:
            _extras[p] = params.pop(p)
        # Get the result
        result = self._func.__call__(a, *args, **params).view(MaskedArray)
        if "fill_value" in common_params:
            result.fill_value = _extras.get("fill_value", None)
        if "hardmask" in common_params:
            result._hardmask = bool(_extras.get("hard_mask", False))
        return result

arange = _convert2ma('arange', params=dict(fill_value=None, hardmask=False))
clip = np.clip
diff = np.diff
empty = _convert2ma('empty', params=dict(fill_value=None, hardmask=False))
empty_like = _convert2ma('empty_like')
frombuffer = _convert2ma('frombuffer')
fromfunction = _convert2ma('fromfunction')
identity = _convert2ma('identity', params=dict(fill_value=None, hardmask=False))
indices = np.indices
ones = _convert2ma('ones', params=dict(fill_value=None, hardmask=False))
ones_like = np.ones_like
squeeze = np.squeeze
zeros = _convert2ma('zeros', params=dict(fill_value=None, hardmask=False))
zeros_like = np.zeros_like

###############################################################################
def append(a, b, axis=None):
    """Append values to the end of an array.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    arr : array_like
        Values are appended to a copy of this array.
    values : array_like
        These values are appended to a copy of `arr`.  It must be of the
        correct shape (the same shape as `arr`, excluding `axis`).  If `axis`
        is not specified, `values` can be any shape and will be flattened
        before use.
    axis : int, optional
        The axis along which `values` are appended.  If `axis` is not given,
        both `arr` and `values` are flattened before use.

    Returns
    -------
    append : MaskedArray
        A copy of `arr` with `values` appended to `axis`.  Note that `append`
        does not occur in-place: a new array is allocated and filled.  If
        `axis` is None, the result is a flattened array.

    See Also
    --------
    numpy.append : Equivalent function in the top-level NumPy module.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = ma.masked_values([1, 2, 3], 2)
    >>> b = ma.masked_values([[4, 5, 6], [7, 8, 9]], 7)
    >>> print(ma.append(a, b))
    [1 -- 3 4 5 6 -- 8 9]
    """
    return concatenate([a, b], axis)
