from copy import deepcopy
import numpy as np

#------------------------------------------------------------------------
# Replace : slightly adapted from bottleneck

loop_template = 'for iINDEX%d in range(nINDEX%d):'
indent = '    '
#replace_op = ('%sif mask[INDEXALL]:\n'
#              '%s    a[INDEXALL] = new%s')

nonna_op = ('%sai = a[INDEXALL]\n'
            '%sif ai == old:\n'
            '%s    a[INDEXALL] = new%s')
na_op = ('%sai = a[INDEXALL]\n'
         '%sif ai != ai:\n'
         '%s    a[INDEXALL] = new%s')

generic_top = """
@cython.boundscheck(False)
@cython.wraparound(False)
def NAME_NDIMd_DTYPE_axisAXIS(np.ndarray[np.DTYPE_t, ndim=NDIM] a,
    double old, double new):
    "replace (inplace) specified elements of NDIMd array of dtype=DTYPE."
    cdef np.DTYPE_t ai
"""

int_check = """\
        oldint = <np.DTYPE_t>old
        newint = <np.DTYPE_t>new
        if oldint != old:
            raise ValueError('Cannot safely cast `old` to int.')
        if newint != new:
            raise ValueError('Cannot safely cast `new` to int.')
"""

def float_loop(ndims=3, type_suffix=''):
    loop = {}
    for n in range(1, ndims + 1):
        loop_str = indent + 'if old==old: \n'
        for i in range(n): # for i in range:
            loop_str += indent * (i + 2) + (loop_template % (i, i)) + '\n'

        dent = indent * (n + 2)
        loop_str += nonna_op % (dent, dent, dent, type_suffix)

        loop_str += '\n' + indent + 'else:\n'
        for i in range(n): # for i in range:
            loop_str += indent * (i + 2) + (loop_template % (i, i)) + '\n'

        dent = indent * (n + 2)
        loop_str += na_op % (dent, dent, dent, type_suffix)

        loop[n] = loop_str + '\n'
    return loop

def int_loop(ndims=3, type_suffix='int'):
    loop = {}
    for n in range(1, ndims + 1):
        loop_str = indent + 'if old==old: \n' + int_check
        for i in range(n): # for i in range:
            loop_str += indent * (i + 2) + (loop_template % (i, i)) + '\n'

        dent = indent * (n + 2)
        loop_str += nonna_op % (dent, dent, dent, type_suffix)
        loop[n] = loop_str + '\n'
    return loop


# float type functions
floats = {}
floats['dtypes'] = ['float32', 'float64']
floats['axisNone'] = True
floats['force_output_dtype'] = False
floats['reuse_non_nan_func'] = False
floats['top'] = generic_top
floats['loop'] = float_loop()

# int type functions
ints = deepcopy(floats)
ints['dtypes'] = ['int32', 'int64']
ints['top'] = generic_top + """
    cdef np.DTYPE_t oldint, newint
    newint = <np.DTYPE_t>new
    if newint != new:
        raise ValueError('Cannot safely cast `new` to int.')
"""
ints['loop'] = int_loop()

# Slow, unaccelerated ndim/dtype --------------------------------------------
def replace(arr, old, new):
    "Slow replace (inplace) used for unaccelerated ndim/dtype combinations."
    if type(arr) is not np.ndarray:
        raise TypeError("`arr` must be a numpy array.")
    if not issubclass(arr.dtype.type, np.inexact):
        if int(old) != old:
            raise ValueError("Cannot safely cast `old` to int.")
        if int(new) != new:
            raise ValueError("Cannot safely cast `new` to int.")
    if old != old:
        mask = np.isnan(arr)
    else:
        mask = arr == old
    np.putmask(arr, mask, new)

slow = {}
slow['name'] = "replace"
slow['signature'] = "arr, old, new"
slow['func'] = "slow_replace(arr, old, new)"

replace = {}
replace['name'] = 'replace'
replace['is_reducing_function'] = False
replace['cdef_output'] = False
replace['slow'] = slow
replace['templates'] = {}
replace['templates']['float_None'] = floats
replace['templates']['int_None'] = ints
replace['pyx_file'] = 'replace.pyx'

replace['main'] = '''"replace auto-generated from template"

def replace(arr, old, new):
    """
    Replace (inplace) given scalar values of an array with new values.

    similar to putmask but faster

    Parameters
    ----------
    arr : numpy.ndarray
        The input array, which is also the output array since this functions
        works inplace.
    old : scalar
    new : scalar
        All masked elements in `arr` will be replaced by `new`.

    Returns
    -------
    None, the operation is inplace.
    """
    func = replace_selector(arr)
    if np.isscalar(old):
        return func(arr, old, new)
    else:
        for o in old:
            func(arr, o, new)
        return arr

def replace_selector(arr):
    """
    Return replace function and array that matches `arr`.

    Under the hood Bottleneck uses a separate replace() Cython function for
    each combination of ndim and dtype. A lot of the overhead in bn.replace()
    is inselecting the low level function to use.

    You can get rid of the overhead by doing all this before you, for example,
    enter an inner loop, by using this function.

    Parameters
    ----------
    arr : numpy.ndarray
        Input array.

    Returns
    -------
    func : function
        The replace() function that matches the number of dimensions and dtype
        of the input array.
    """
    axis = None
    if type(arr) is not np.ndarray:
        raise TypeError("`arr` must be a numpy array.")
    cdef int ndim = PyArray_NDIM(arr)
    cdef int dtype = PyArray_TYPE(arr)
    cdef tuple key = (ndim, dtype, axis)
    try:
        func = replace_dict[key]
    except KeyError:
        try:
            func = replace_slow_dict[axis]
        except KeyError:
            tup = (str(ndim), str(arr.dtype), str(axis))
            raise TypeError("Unsupported ndim/dtype/axis (%s/%s/%s)." % tup)
    return func
'''
