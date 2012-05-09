"replace auto-generated from template"

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
            arr = func(arr, o, new)
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

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_1d_int32_axisNone(np.ndarray[np.int32_t, ndim=1] a,
    double old, double new):
    "replace (inplace) specified elements of 1d array of dtype=int32."
    cdef np.int32_t ai

    cdef np.int32_t oldint, newint
    newint = <np.int32_t>new
    if newint != new:
        raise ValueError('Cannot safely cast `new` to int.')
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    if old==old: 
        oldint = <np.int32_t>old
        newint = <np.int32_t>new
        if oldint != old:
            raise ValueError('Cannot safely cast `old` to int.')
        if newint != new:
            raise ValueError('Cannot safely cast `new` to int.')
        for i0 in range(n0):
            ai = a[i0]
            if ai == old:
                a[i0] = newint

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_1d_int64_axisNone(np.ndarray[np.int64_t, ndim=1] a,
    double old, double new):
    "replace (inplace) specified elements of 1d array of dtype=int64."
    cdef np.int64_t ai

    cdef np.int64_t oldint, newint
    newint = <np.int64_t>new
    if newint != new:
        raise ValueError('Cannot safely cast `new` to int.')
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    if old==old: 
        oldint = <np.int64_t>old
        newint = <np.int64_t>new
        if oldint != old:
            raise ValueError('Cannot safely cast `old` to int.')
        if newint != new:
            raise ValueError('Cannot safely cast `new` to int.')
        for i0 in range(n0):
            ai = a[i0]
            if ai == old:
                a[i0] = newint

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_2d_int32_axisNone(np.ndarray[np.int32_t, ndim=2] a,
    double old, double new):
    "replace (inplace) specified elements of 2d array of dtype=int32."
    cdef np.int32_t ai

    cdef np.int32_t oldint, newint
    newint = <np.int32_t>new
    if newint != new:
        raise ValueError('Cannot safely cast `new` to int.')
    cdef Py_ssize_t i0, i1
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef Py_ssize_t n1 = dim[1]
    if old==old: 
        oldint = <np.int32_t>old
        newint = <np.int32_t>new
        if oldint != old:
            raise ValueError('Cannot safely cast `old` to int.')
        if newint != new:
            raise ValueError('Cannot safely cast `new` to int.')
        for i0 in range(n0):
            for i1 in range(n1):
                ai = a[i0, i1]
                if ai == old:
                    a[i0, i1] = newint

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_2d_int64_axisNone(np.ndarray[np.int64_t, ndim=2] a,
    double old, double new):
    "replace (inplace) specified elements of 2d array of dtype=int64."
    cdef np.int64_t ai

    cdef np.int64_t oldint, newint
    newint = <np.int64_t>new
    if newint != new:
        raise ValueError('Cannot safely cast `new` to int.')
    cdef Py_ssize_t i0, i1
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef Py_ssize_t n1 = dim[1]
    if old==old: 
        oldint = <np.int64_t>old
        newint = <np.int64_t>new
        if oldint != old:
            raise ValueError('Cannot safely cast `old` to int.')
        if newint != new:
            raise ValueError('Cannot safely cast `new` to int.')
        for i0 in range(n0):
            for i1 in range(n1):
                ai = a[i0, i1]
                if ai == old:
                    a[i0, i1] = newint

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_3d_int32_axisNone(np.ndarray[np.int32_t, ndim=3] a,
    double old, double new):
    "replace (inplace) specified elements of 3d array of dtype=int32."
    cdef np.int32_t ai

    cdef np.int32_t oldint, newint
    newint = <np.int32_t>new
    if newint != new:
        raise ValueError('Cannot safely cast `new` to int.')
    cdef Py_ssize_t i0, i1, i2
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef Py_ssize_t n1 = dim[1]
    cdef Py_ssize_t n2 = dim[2]
    if old==old: 
        oldint = <np.int32_t>old
        newint = <np.int32_t>new
        if oldint != old:
            raise ValueError('Cannot safely cast `old` to int.')
        if newint != new:
            raise ValueError('Cannot safely cast `new` to int.')
        for i0 in range(n0):
            for i1 in range(n1):
                for i2 in range(n2):
                    ai = a[i0, i1, i2]
                    if ai == old:
                        a[i0, i1, i2] = newint

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_3d_int64_axisNone(np.ndarray[np.int64_t, ndim=3] a,
    double old, double new):
    "replace (inplace) specified elements of 3d array of dtype=int64."
    cdef np.int64_t ai

    cdef np.int64_t oldint, newint
    newint = <np.int64_t>new
    if newint != new:
        raise ValueError('Cannot safely cast `new` to int.')
    cdef Py_ssize_t i0, i1, i2
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef Py_ssize_t n1 = dim[1]
    cdef Py_ssize_t n2 = dim[2]
    if old==old: 
        oldint = <np.int64_t>old
        newint = <np.int64_t>new
        if oldint != old:
            raise ValueError('Cannot safely cast `old` to int.')
        if newint != new:
            raise ValueError('Cannot safely cast `new` to int.')
        for i0 in range(n0):
            for i1 in range(n1):
                for i2 in range(n2):
                    ai = a[i0, i1, i2]
                    if ai == old:
                        a[i0, i1, i2] = newint

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_1d_float32_axisNone(np.ndarray[np.float32_t, ndim=1] a,
    double old, double new):
    "replace (inplace) specified elements of 1d array of dtype=float32."
    cdef np.float32_t ai
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    if old==old: 
        for i0 in range(n0):
            ai = a[i0]
            if ai == old:
                a[i0] = new
    else:
        for i0 in range(n0):
            ai = a[i0]
            if ai != ai:
                a[i0] = new

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_1d_float64_axisNone(np.ndarray[np.float64_t, ndim=1] a,
    double old, double new):
    "replace (inplace) specified elements of 1d array of dtype=float64."
    cdef np.float64_t ai
    cdef Py_ssize_t i0
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    if old==old: 
        for i0 in range(n0):
            ai = a[i0]
            if ai == old:
                a[i0] = new
    else:
        for i0 in range(n0):
            ai = a[i0]
            if ai != ai:
                a[i0] = new

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_2d_float32_axisNone(np.ndarray[np.float32_t, ndim=2] a,
    double old, double new):
    "replace (inplace) specified elements of 2d array of dtype=float32."
    cdef np.float32_t ai
    cdef Py_ssize_t i0, i1
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef Py_ssize_t n1 = dim[1]
    if old==old: 
        for i0 in range(n0):
            for i1 in range(n1):
                ai = a[i0, i1]
                if ai == old:
                    a[i0, i1] = new
    else:
        for i0 in range(n0):
            for i1 in range(n1):
                ai = a[i0, i1]
                if ai != ai:
                    a[i0, i1] = new

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_2d_float64_axisNone(np.ndarray[np.float64_t, ndim=2] a,
    double old, double new):
    "replace (inplace) specified elements of 2d array of dtype=float64."
    cdef np.float64_t ai
    cdef Py_ssize_t i0, i1
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef Py_ssize_t n1 = dim[1]
    if old==old: 
        for i0 in range(n0):
            for i1 in range(n1):
                ai = a[i0, i1]
                if ai == old:
                    a[i0, i1] = new
    else:
        for i0 in range(n0):
            for i1 in range(n1):
                ai = a[i0, i1]
                if ai != ai:
                    a[i0, i1] = new

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_3d_float32_axisNone(np.ndarray[np.float32_t, ndim=3] a,
    double old, double new):
    "replace (inplace) specified elements of 3d array of dtype=float32."
    cdef np.float32_t ai
    cdef Py_ssize_t i0, i1, i2
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef Py_ssize_t n1 = dim[1]
    cdef Py_ssize_t n2 = dim[2]
    if old==old: 
        for i0 in range(n0):
            for i1 in range(n1):
                for i2 in range(n2):
                    ai = a[i0, i1, i2]
                    if ai == old:
                        a[i0, i1, i2] = new
    else:
        for i0 in range(n0):
            for i1 in range(n1):
                for i2 in range(n2):
                    ai = a[i0, i1, i2]
                    if ai != ai:
                        a[i0, i1, i2] = new

@cython.boundscheck(False)
@cython.wraparound(False)
def replace_3d_float64_axisNone(np.ndarray[np.float64_t, ndim=3] a,
    double old, double new):
    "replace (inplace) specified elements of 3d array of dtype=float64."
    cdef np.float64_t ai
    cdef Py_ssize_t i0, i1, i2
    cdef np.npy_intp *dim
    dim = PyArray_DIMS(a)
    cdef Py_ssize_t n0 = dim[0]
    cdef Py_ssize_t n1 = dim[1]
    cdef Py_ssize_t n2 = dim[2]
    if old==old: 
        for i0 in range(n0):
            for i1 in range(n1):
                for i2 in range(n2):
                    ai = a[i0, i1, i2]
                    if ai == old:
                        a[i0, i1, i2] = new
    else:
        for i0 in range(n0):
            for i1 in range(n1):
                for i2 in range(n2):
                    ai = a[i0, i1, i2]
                    if ai != ai:
                        a[i0, i1, i2] = new

cdef dict replace_dict = {}
replace_dict[(1, NPY_int32, 0)] = replace_1d_int32_axisNone
replace_dict[(1, NPY_int32, None)] = replace_1d_int32_axisNone
replace_dict[(1, NPY_int64, 0)] = replace_1d_int64_axisNone
replace_dict[(1, NPY_int64, None)] = replace_1d_int64_axisNone
replace_dict[(2, NPY_int32, None)] = replace_2d_int32_axisNone
replace_dict[(2, NPY_int64, None)] = replace_2d_int64_axisNone
replace_dict[(3, NPY_int32, None)] = replace_3d_int32_axisNone
replace_dict[(3, NPY_int64, None)] = replace_3d_int64_axisNone
replace_dict[(1, NPY_float32, 0)] = replace_1d_float32_axisNone
replace_dict[(1, NPY_float32, None)] = replace_1d_float32_axisNone
replace_dict[(1, NPY_float64, 0)] = replace_1d_float64_axisNone
replace_dict[(1, NPY_float64, None)] = replace_1d_float64_axisNone
replace_dict[(2, NPY_float32, None)] = replace_2d_float32_axisNone
replace_dict[(2, NPY_float64, None)] = replace_2d_float64_axisNone
replace_dict[(3, NPY_float32, None)] = replace_3d_float32_axisNone
replace_dict[(3, NPY_float64, None)] = replace_3d_float64_axisNone

def replace_slow_axis0(arr, old, new):
    "Unaccelerated (slow) replace along axis 0."
    return slow_replace(arr, old, new)

def replace_slow_axis1(arr, old, new):
    "Unaccelerated (slow) replace along axis 1."
    return slow_replace(arr, old, new)

def replace_slow_axis2(arr, old, new):
    "Unaccelerated (slow) replace along axis 2."
    return slow_replace(arr, old, new)

def replace_slow_axis3(arr, old, new):
    "Unaccelerated (slow) replace along axis 3."
    return slow_replace(arr, old, new)

def replace_slow_axis4(arr, old, new):
    "Unaccelerated (slow) replace along axis 4."
    return slow_replace(arr, old, new)

def replace_slow_axis5(arr, old, new):
    "Unaccelerated (slow) replace along axis 5."
    return slow_replace(arr, old, new)

def replace_slow_axis6(arr, old, new):
    "Unaccelerated (slow) replace along axis 6."
    return slow_replace(arr, old, new)

def replace_slow_axis7(arr, old, new):
    "Unaccelerated (slow) replace along axis 7."
    return slow_replace(arr, old, new)

def replace_slow_axis8(arr, old, new):
    "Unaccelerated (slow) replace along axis 8."
    return slow_replace(arr, old, new)

def replace_slow_axis9(arr, old, new):
    "Unaccelerated (slow) replace along axis 9."
    return slow_replace(arr, old, new)

def replace_slow_axis10(arr, old, new):
    "Unaccelerated (slow) replace along axis 10."
    return slow_replace(arr, old, new)

def replace_slow_axis11(arr, old, new):
    "Unaccelerated (slow) replace along axis 11."
    return slow_replace(arr, old, new)

def replace_slow_axis12(arr, old, new):
    "Unaccelerated (slow) replace along axis 12."
    return slow_replace(arr, old, new)

def replace_slow_axis13(arr, old, new):
    "Unaccelerated (slow) replace along axis 13."
    return slow_replace(arr, old, new)

def replace_slow_axis14(arr, old, new):
    "Unaccelerated (slow) replace along axis 14."
    return slow_replace(arr, old, new)

def replace_slow_axis15(arr, old, new):
    "Unaccelerated (slow) replace along axis 15."
    return slow_replace(arr, old, new)

def replace_slow_axis16(arr, old, new):
    "Unaccelerated (slow) replace along axis 16."
    return slow_replace(arr, old, new)

def replace_slow_axis17(arr, old, new):
    "Unaccelerated (slow) replace along axis 17."
    return slow_replace(arr, old, new)

def replace_slow_axis18(arr, old, new):
    "Unaccelerated (slow) replace along axis 18."
    return slow_replace(arr, old, new)

def replace_slow_axis19(arr, old, new):
    "Unaccelerated (slow) replace along axis 19."
    return slow_replace(arr, old, new)

def replace_slow_axis20(arr, old, new):
    "Unaccelerated (slow) replace along axis 20."
    return slow_replace(arr, old, new)

def replace_slow_axis21(arr, old, new):
    "Unaccelerated (slow) replace along axis 21."
    return slow_replace(arr, old, new)

def replace_slow_axis22(arr, old, new):
    "Unaccelerated (slow) replace along axis 22."
    return slow_replace(arr, old, new)

def replace_slow_axis23(arr, old, new):
    "Unaccelerated (slow) replace along axis 23."
    return slow_replace(arr, old, new)

def replace_slow_axis24(arr, old, new):
    "Unaccelerated (slow) replace along axis 24."
    return slow_replace(arr, old, new)

def replace_slow_axis25(arr, old, new):
    "Unaccelerated (slow) replace along axis 25."
    return slow_replace(arr, old, new)

def replace_slow_axis26(arr, old, new):
    "Unaccelerated (slow) replace along axis 26."
    return slow_replace(arr, old, new)

def replace_slow_axis27(arr, old, new):
    "Unaccelerated (slow) replace along axis 27."
    return slow_replace(arr, old, new)

def replace_slow_axis28(arr, old, new):
    "Unaccelerated (slow) replace along axis 28."
    return slow_replace(arr, old, new)

def replace_slow_axis29(arr, old, new):
    "Unaccelerated (slow) replace along axis 29."
    return slow_replace(arr, old, new)

def replace_slow_axis30(arr, old, new):
    "Unaccelerated (slow) replace along axis 30."
    return slow_replace(arr, old, new)

def replace_slow_axis31(arr, old, new):
    "Unaccelerated (slow) replace along axis 31."
    return slow_replace(arr, old, new)

def replace_slow_axis32(arr, old, new):
    "Unaccelerated (slow) replace along axis 32."
    return slow_replace(arr, old, new)

def replace_slow_axisNone(arr, old, new):
    "Unaccelerated (slow) replace along axis None."
    return slow_replace(arr, old, new)


cdef dict replace_slow_dict = {}
replace_slow_dict[0] = replace_slow_axis0
replace_slow_dict[1] = replace_slow_axis1
replace_slow_dict[2] = replace_slow_axis2
replace_slow_dict[3] = replace_slow_axis3
replace_slow_dict[4] = replace_slow_axis4
replace_slow_dict[5] = replace_slow_axis5
replace_slow_dict[6] = replace_slow_axis6
replace_slow_dict[7] = replace_slow_axis7
replace_slow_dict[8] = replace_slow_axis8
replace_slow_dict[9] = replace_slow_axis9
replace_slow_dict[10] = replace_slow_axis10
replace_slow_dict[11] = replace_slow_axis11
replace_slow_dict[12] = replace_slow_axis12
replace_slow_dict[13] = replace_slow_axis13
replace_slow_dict[14] = replace_slow_axis14
replace_slow_dict[15] = replace_slow_axis15
replace_slow_dict[16] = replace_slow_axis16
replace_slow_dict[17] = replace_slow_axis17
replace_slow_dict[18] = replace_slow_axis18
replace_slow_dict[19] = replace_slow_axis19
replace_slow_dict[20] = replace_slow_axis20
replace_slow_dict[21] = replace_slow_axis21
replace_slow_dict[22] = replace_slow_axis22
replace_slow_dict[23] = replace_slow_axis23
replace_slow_dict[24] = replace_slow_axis24
replace_slow_dict[25] = replace_slow_axis25
replace_slow_dict[26] = replace_slow_axis26
replace_slow_dict[27] = replace_slow_axis27
replace_slow_dict[28] = replace_slow_axis28
replace_slow_dict[29] = replace_slow_axis29
replace_slow_dict[30] = replace_slow_axis30
replace_slow_dict[31] = replace_slow_axis31
replace_slow_dict[32] = replace_slow_axis32
replace_slow_dict[None] = replace_slow_axisNone