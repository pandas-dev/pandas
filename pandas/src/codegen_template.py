"Copied from bottleneck: Turn templates into Cython pyx files."
import os.path

def template(func):
    "'Convert template dictionary `func` to a pyx file.'\n"
    codes = []
    codes.append(func['main'])
    select = Selector(func['name'])
    for key in func['templates']:
        f = func['templates'][key]
        code = subtemplate(name=func['name'],
                           top=f['top'],
                           loop=f['loop'],
                           axisNone=f['axisNone'],
                           dtypes=f['dtypes'],
                           force_output_dtype=f['force_output_dtype'],
                           reuse_non_nan_func=f['reuse_non_nan_func'],
                           is_reducing_function=func['is_reducing_function'],
                           cdef_output=func['cdef_output'],
                           select=select)
        codes.append(code)
    codes.append('\n' + str(select))
    if 'slow' in func:
        if func['slow'] is not None:
            slow = func['slow']
            code1 = slow_selector(slow['name'])
            code2 = slow_functions(slow['name'],
                                   slow['signature'],
                                   slow['func'])
            codes.append(code2)
            codes.append(code1)
    modpath = os.path.dirname(__file__)
    fid = open(os.path.join(modpath, func['pyx_file']), 'w')
    fid.write(''.join(codes))
    fid.close()

def subtemplate(name, top, loop, axisNone, dtypes, force_output_dtype,
                reuse_non_nan_func, is_reducing_function, cdef_output, select):
    "Assemble template"
    ndims = sorted(loop.keys())
    funcs = []
    for ndim in ndims:
        if axisNone:
            axes = [None]
        else:
            axes = list(range(ndim))
        for dtype in dtypes:
            for axis in axes:

                if reuse_non_nan_func:

                    select.append(ndim, dtype, axis, True)

                else:

                    # Code template
                    func = top

                    # loop
                    if force_output_dtype is not False:
                        ydtype = force_output_dtype
                    else:
                        ydtype = dtype
                    func += loop_cdef(ndim, ydtype, axis, is_reducing_function,
                                      cdef_output)
                    func += looper(loop[ndim], ndim, axis)

                    # name, ndim, dtype, axis
                    func = func.replace('NAME', name)
                    func = func.replace('NDIM', str(ndim))
                    func = func.replace('DTYPE', dtype)
                    func = func.replace('AXIS', str(axis))

                    funcs.append(func)
                    select.append(ndim, dtype, axis)

    return ''.join(funcs)

def looper(loop, ndim, axis):
    """
    Given loop template, expand index markers for given `ndim` and `axis`.

    Parameters
    ----------
    loop : str
        Code of loop where the following template markers will be expanded
        (example given is for 3d input, similarly for other nd):

        ================= =================================================
        INDEXALL          Replace with i0, i1, i2
        INDEXPOP          If axis=1, e.g., replace with i0, i2
        INDEXN            If N=1, e.g., replace with 1
        INDEXREPLACE|exp| If exp = 'k - window' and axis=1, e.g., replace
                          with i0, k - window, i2
        NREPLACE|exp|     If exp = 'n - window' and axis=1, e.g., replace
                          with n0, n - window, n2
        ================= =================================================
    ndim : int
        Number of dimensions in the loop.
    axis : {int, None}
        Axis over which the loop is evaluated.

    Returns
    -------
    code : str
        Code for the loop with templated index markers expanded.

    Examples
    --------
    Make a 3d loop template:

    >>> loop = '''
    .... for iINDEX0 in range(nINDEX0):
    ....    for iINDEX1 in range(nINDEX1):
    ....        amin = MAXDTYPE
    ....        for iINDEX2 in range(nINDEX2):
    ....            ai = a[INDEXALL]
    ....            if ai <= amin:
    ....                amin = ai
    ....         y[INDEXPOP] = amin
    .... '''

    Import the looper function:

    >>> from bottleneck.src.template.template import looper

    Make a loop over axis=0:

    >>> print(looper(loop, ndim=3, axis=0))
    for i1 in range(n1):
        for i2 in range(n2):
            amin = MAXDTYPE
            for i0 in range(n0):
                ai = a[i0, i1, i2]
                if ai <= amin:
                    amin = ai
            y[i1, i2] = amin

    Make a loop over axis=1:

    >>> print(looper(loop, ndim=3, axis=1))
    for i0 in range(n0):
        for i2 in range(n2):
            amin = MAXDTYPE
            for i1 in range(n1):
                ai = a[i0, i1, i2]
                if ai <= amin:
                    amin = ai
            y[i0, i2] = amin

    Make a loop over axis=2:

    >>> print(looper(loop, ndim=3, axis=2))
    for i0 in range(n0):
        for i1 in range(n1):
            amin = MAXDTYPE
            for i2 in range(n2):
                ai = a[i0, i1, i2]
                if ai <= amin:
                    amin = ai
            y[i0, i1] = amin

    """

    if ndim < 1:
        raise ValueError("ndim(=%d) must be and integer greater than 0" % ndim)
    if axis is not None:
        if axis < 0:
            raise ValueError("`axis` must be a non-negative integer or None")
        elif axis >= ndim:
            raise ValueError("`axis` must be less then `ndim`")

    # INDEXALL
    INDEXALL = ', '.join('i' + str(i) for i in range(ndim))
    code = loop.replace('INDEXALL', INDEXALL)

    # INDEXPOP
    idx = list(range(ndim))
    if axis is not None:
        idx.pop(axis)
    INDEXPOP = ', '.join(['i' + str(i) for i in idx])
    code = code.replace('INDEXPOP', INDEXPOP)

    # INDEXN
    idx = list(range(ndim))
    if axis is not None:
        idxpop = idx.pop(axis)
        idx.append(idxpop)
    for i, j in enumerate(idx):
        code = code.replace('INDEX%d' % i, '%d' % j)

    # INDEXREPLACE|x|
    mark = 'INDEXREPLACE|'
    nreplace = code.count(mark)
    if (nreplace > 0) and (axis is None):
        raise ValueError("`INDEXREPLACE` cannot be used when axis is None.")
    while mark in code:
        idx0 = code.index(mark)
        idx1 = idx0 + len(mark)
        idx2 = idx1 + code[idx1:].index('|')
        if (idx0 >= idx1) or (idx1 >= idx2):
            raise RuntimeError("Parsing error or poorly formatted input.")
        replacement = code[idx1:idx2]
        idx = ['i' + str(i) for i in range(ndim)]
        idx[axis] = replacement
        idx = ', '.join(idx)
        code = code[:idx0] + idx + code[idx2+1:]

    # NREPLACE|x|
    mark = 'NREPLACE|'
    nreplace = code.count(mark)
    # TODO: reuse while loop above, only difference is 'i' --> 'n'
    while mark in code:
        idx0 = code.index(mark)
        idx1 = idx0 + len(mark)
        idx2 = idx1 + code[idx1:].index('|')
        if (idx0 >= idx1) or (idx1 >= idx2):
            raise RuntimeError("Parsing error or poorly formatted input.")
        replacement = code[idx1:idx2]
        idx = ['n' + str(i) for i in range(ndim)]
        idx[axis] = replacement
        idx = ', '.join(idx)
        code = code[:idx0] + idx + code[idx2+1:]

    return code

def loop_cdef(ndim, dtype, axis, is_reducing_function, cdef_output=True):
    """
    String of code that initializes variables needed in a for loop.

    The output string contains code for: index array counters, one for each
    dimension (cdef Py_size_t i0, i1, i2, ....); the length along each
    dimension of the input array, `a` (cdef Py_ssize_t n0 = a.shape[0],...);
    the initialized, empty output array, `y`.

    Parameters
    ----------
    ndim = int
        Number of dimensions.
    dtype : str
        The data type of the output. Used for initilizing the empty output
        array, `y`.
    is_reducing_function : bool
        If True then remove the dimension given by `axis` when initializing
        the output array, `y`.
    cdef_output : bool, optional
        If False then only initialize indices (i) and shapes (n). If True
        (default) then also intialized output array `y`.

    Returns
    -------
    cdefs : str
        String of code to use to initialize variables needed for loop.

    Examples
    --------
    Define parameters:

    >>> ndim = 3
    >>> dtype = 'float64'
    >>> axis = 1
    >>> is_reducing_function = True

    Import loop_cdef:

    >>> from bottleneck.src.template.template import loop_cdef

    Make loop initialization code:

    >>> print(loop_cdef(ndim, dtype, axis, is_reducing_function))
        cdef Py_ssize_t i0, i1, i2
        cdef np.npy_intp *dim
        dim = PyArray_DIMS(a)
        Py_ssize_t n0 = dim[0]
        Py_ssize_t n1 = dim[1]
        Py_ssize_t n2 = dim[2]
        cdef np.npy_intp *dims = [n0, n2]
        cdef np.ndarray[np.float64_t, ndim=2] y = PyArray_EMPTY(2, dims,
                                                  NPY_float64, 0)

    Repeat, but this time make the output non-reducing:

    >>> is_reducing_function = False
    >>> print(loop_cdef(ndim, dtype, axis, is_reducing_function))
        cdef Py_ssize_t i0, i1, i2
        cdef np.npy_intp *dim
        dim = PyArray_DIMS(a)
        Py_ssize_t n0 = dim[0]
        Py_ssize_t n1 = dim[1]
        Py_ssize_t n2 = dim[2]
        cdef np.npy_intp *dims = [n0, n1, n2]
        cdef np.ndarray[np.float64_t, ndim=3] y = PyArray_EMPTY(3, dims,
                                                  NPY_float64, 0)

    """

    if ndim < 1:
        raise ValueError("ndim(=%d) must be and integer greater than 0" % ndim)
    if axis is not None:
        if axis < 0:
            raise ValueError("`axis` must be a non-negative integer or None")
        elif axis >= ndim:
            raise ValueError("`axis` must be less then `ndim`")

    tab = '    '
    cdefs = []

    # cdef loop indices
    idx = ', '.join('i'+str(i) for i in range(ndim))
    cdefs.append(tab + 'cdef Py_ssize_t ' + idx)

    # Length along each dimension
    cdefs.append(tab + "cdef np.npy_intp *dim")
    cdefs.append(tab + "dim = PyArray_DIMS(a)")
    for dim in range(ndim):
        cdefs.append(tab + "cdef Py_ssize_t n%d = dim[%d]" % (dim, dim))

    if not cdef_output:
        return '\n'.join(cdefs) + '\n'

    # cdef initialize output
    if is_reducing_function:
        if (ndim > 1) and (axis is not None):
            idx = list(range(ndim))
            del idx[axis]
            ns = ', '.join(['n'+str(i) for i in idx])
            cdefs.append("%scdef np.npy_intp *dims = [%s]" % (tab, ns))
            y = "%scdef np.ndarray[np.%s_t, ndim=%d] "
            y += "y = PyArray_EMPTY(%d, dims,"
            y += "\n                                              NPY_%s, 0)"
            cdefs.append(y % (tab, dtype, ndim-1, ndim-1, dtype))
    else:
        idx = list(range(ndim))
        ns = ', '.join('n'+str(i) for i in idx)
        cdefs.append("%scdef np.npy_intp *dims = [%s]" % (tab, ns))
        y = "%scdef np.ndarray[np.%s_t, ndim=%d] "
        y += "y = PyArray_EMPTY(%d, dims,"
        y += "\n                                              NPY_%s, 0)"
        cdefs.append(y % (tab, dtype, ndim, ndim, dtype))

    return '\n'.join(cdefs) + '\n'

class Selector(object):
    "String of code for dictionary that maps dtype to cython function."

    def __init__(self, name):
        self.name = name
        self.data = []

    def append(self, ndim, dtype, axis, reuse=False):
        self.data.append((ndim, dtype, axis, reuse))

    def __str__(self):
        fmt = "%s_dict[(%s, NPY_%s, %s)] = %s_%sd_%s_axis%s"
        src = []
        src.append("cdef dict %s_dict = {}" % self.name)
        for ndim, dtype, axis, reuse in self.data:
            name = self.name
            if reuse:
                name = name.replace('nan', '')
            if (ndim == 1) and (axis is None):
                tup = (self.name, str(ndim), str(dtype), str(0),
                       name, str(ndim), str(dtype), str(axis))
                src.append(fmt % tup)
            tup = (self.name, str(ndim), str(dtype), str(axis),
                   name, str(ndim), str(dtype), str(axis))
            src.append(fmt % tup)
        return '\n'.join(src)

def slow_selector(name, maxaxis=32):
    "String of code for slow function mapping dictionary."
    axes = list(range(maxaxis+1)) + [None]
    src = ['\n']
    src.append("cdef dict %s_slow_dict = {}" % name)
    fmt = "%s_slow_dict[%s] = %s_slow_axis%s"
    for axis in axes:
        tup = 2 * (name, str(axis))
        src.append(fmt % tup)
    return '\n'.join(src)

def slow_functions(name, signature, func, maxaxis=32):
    "String of code for slow functions."
    axes = list(range(maxaxis+1)) + [None]
    tab = '    '
    sig = "def %s_slow_axis%s(%s):"
    doc = '%s"Unaccelerated (slow) %s along axis %s."'
    function = "%sreturn %s\n"
    src = ['\n']
    for axis in axes:

        axis = str(axis)

        # signature
        code = sig % (name, axis, signature)
        code = code.replace('AXIS', axis)
        src.append(code)

        # docstring
        code = doc % (tab, name, axis)
        code = code.replace('AXIS', axis)
        src.append(code)

        # function
        code = function % (tab, func)
        code = code.replace('AXIS', axis)
        src.append(code)

    return '\n'.join(src)
