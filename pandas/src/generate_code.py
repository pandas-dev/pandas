from cStringIO import StringIO

take_1d_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_%(name)s(ndarray[%(c_type)s] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[%(c_type)s] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            %(na_action)s
        else:
            outbuf[i] = values[idx]

"""

take_2d_axis0_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_%(name)s(ndarray[%(c_type)s, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[%(c_type)s, ndim=2] outbuf

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]

        if idx == -1:
            for j from 0 <= j < k:
                %(na_action)s
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

"""

take_2d_axis1_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_%(name)s(ndarray[%(c_type)s, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[%(c_type)s, ndim=2] outbuf

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    for j from 0 <= j < k:
        idx = indexer[j]

        if idx == -1:
            for i from 0 <= i < n:
                %(na_action)s
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

"""

set_na = "outbuf[i] = NaN"
set_na_2d = "outbuf[i, j] = NaN"
raise_on_na = "raise ValueError('No NA values allowed')"

merge_indexer_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def merge_indexer_%(name)s(ndarray[%(c_type)s] values, dict oldMap):
    cdef int i, j, length, newLength
    cdef %(c_type)s idx
    cdef ndarray[int32_t] fill_vec

    newLength = len(values)
    fill_vec = np.empty(newLength, dtype=np.int32)
    for i from 0 <= i < newLength:
        idx = values[i]
        if idx in oldMap:
            fill_vec[i] = oldMap[idx]
        else:
            fill_vec[i] = -1

    return fill_vec

"""

# name, ctype, capable of holding NA
function_list = [
    ('float64', 'float64_t', True),
    ('object', 'object', True),
    ('int32', 'int32_t', False),
    ('int64', 'int64_t', False),
    ('bool', 'uint8_t', False)
]

def generate_from_template(template, ndim=1):
    output = StringIO()
    for name, c_type, can_hold_na in function_list:
        if ndim == 1:
            na_action = set_na if can_hold_na else raise_on_na
        elif ndim == 2:
            na_action = set_na_2d if can_hold_na else raise_on_na
        func = template % {'name' : name, 'c_type' : c_type,
                           'na_action' : na_action}
        output.write(func)
    return output.getvalue()

def generate_take_cython_file(path='generated.pyx'):
    with open(path, 'w') as f:
        print >> f, generate_from_template(merge_indexer_template)
        print >> f, generate_from_template(take_1d_template)
        print >> f, generate_from_template(take_2d_axis0_template, ndim=2)
        print >> f, generate_from_template(take_2d_axis1_template, ndim=2)

if __name__ == '__main__':
    generate_take_cython_file()
