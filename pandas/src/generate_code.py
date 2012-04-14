from cStringIO import StringIO

take_1d_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_%(name)s(ndarray[%(c_type)s] values,
                     ndarray[int32_t] indexer,
                     out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[%(c_type)s] outbuf
        %(c_type)s fv

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    if %(raise_on_na)s and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                raise ValueError('No NA values allowed')
            else:
                outbuf[i] = values[idx]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                outbuf[i] = fv
            else:
                outbuf[i] = values[idx]

"""

take_2d_axis0_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_%(name)s(ndarray[%(c_type)s, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[%(c_type)s, ndim=2] outbuf
        %(c_type)s fv

    n = len(indexer)
    k = values.shape[1]

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if %(raise_on_na)s and _checknan(fill_value):
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j from 0 <= j < k:
                    raise ValueError('No NA values allowed')
            else:
                for j from 0 <= j < k:
                    outbuf[i, j] = values[idx, j]
    else:
        fv = fill_value
        for i in range(n):
            idx = indexer[i]
            if idx == -1:
                for j from 0 <= j < k:
                    outbuf[i, j] = fv
            else:
                for j from 0 <= j < k:
                    outbuf[i, j] = values[idx, j]

"""

take_2d_axis1_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_%(name)s(ndarray[%(c_type)s, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None, fill_value=np.nan):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[%(c_type)s, ndim=2] outbuf
        %(c_type)s fv

    n = len(values)
    k = len(indexer)

    if out is None:
        outbuf = np.empty((n, k), dtype=values.dtype)
    else:
        outbuf = out

    if %(raise_on_na)s and _checknan(fill_value):
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    raise ValueError('No NA values allowed')
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]
    else:
        fv = fill_value
        for j in range(k):
            idx = indexer[j]

            if idx == -1:
                for i in range(n):
                    outbuf[i, j] = fv
            else:
                for i in range(n):
                    outbuf[i, j] = values[i, idx]

"""

def set_na(na ="NaN"):
    return "outbuf[i] = %s" % na

def set_na_2d(na = "NaN"):
    return "outbuf[i, j] = %s" % na

raise_on_na = "raise ValueError('No NA values allowed')"

'''
Backfilling logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .        0               1
         .        0               1
         .        0               1
A        A        0               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
         .        1               1
B        B        1               1
         .        2               1
         .        2               1
         .        2               1
C        C        2               1
         .                        0
         .                        0
D
'''

backfill_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_%(name)s(ndarray[%(c_type)s] old, ndarray[%(c_type)s] new,
                      limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int32_t, ndim=1] indexer
    cdef %(c_type)s cur, prev
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int32)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[0] > old[nleft - 1]:
        return indexer

    i = nleft - 1
    j = nright - 1

    cur = old[nleft - 1]

    while j >= 0 and new[j] > cur:
        j -= 1

    while True:
        if j < 0:
            break

        if i == 0:
            while j >= 0:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] < cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j -= 1
            break

        prev = old[i - 1]

        while j >= 0 and prev < new[j] <= cur:
            if new[j] == cur:
                indexer[j] = i
            elif new[j] < cur and fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j -= 1

        fill_count = 0
        i -= 1
        cur = prev

    return indexer

"""


pad_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def pad_%(name)s(ndarray[%(c_type)s] old, ndarray[%(c_type)s] new,
                   limit=None):
    cdef Py_ssize_t i, j, nleft, nright
    cdef ndarray[int32_t, ndim=1] indexer
    cdef %(c_type)s cur, next
    cdef int lim, fill_count = 0

    nleft = len(old)
    nright = len(new)
    indexer = np.empty(nright, dtype=np.int32)
    indexer.fill(-1)

    if limit is None:
        lim = nright
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    if nleft == 0 or nright == 0 or new[nright - 1] < old[0]:
        return indexer

    i = j = 0

    cur = old[0]

    while j <= nright - 1 and new[j] < cur:
        j += 1

    while True:
        if j == nright:
            break

        if i == nleft - 1:
            while j < nright:
                if new[j] == cur:
                    indexer[j] = i
                elif new[j] > cur and fill_count < lim:
                    indexer[j] = i
                    fill_count += 1
                j += 1
            break

        next = old[i + 1]

        while j < nright and cur <= new[j] < next:
            if new[j] == cur:
                indexer[j] = i
            elif fill_count < lim:
                indexer[j] = i
                fill_count += 1
            j += 1

        fill_count = 0
        i += 1
        cur = next

    return indexer

"""

pad_1d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def pad_inplace_%(name)s(ndarray[%(c_type)s] values,
                         ndarray[uint8_t, cast=True] mask,
                         limit=None):
    cdef Py_ssize_t i, N
    cdef %(c_type)s val
    cdef int lim, fill_count = 0

    N = len(values)

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[0]
    for i in range(N):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]

"""

pad_2d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def pad_2d_inplace_%(name)s(ndarray[%(c_type)s, ndim=2] values,
                            ndarray[uint8_t, ndim=2] mask,
                            limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef %(c_type)s val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, 0]
        for i in range(N):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]
"""

backfill_2d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_2d_inplace_%(name)s(ndarray[%(c_type)s, ndim=2] values,
                                 ndarray[uint8_t, ndim=2] mask,
                                 limit=None):
    cdef Py_ssize_t i, j, N, K
    cdef %(c_type)s val
    cdef int lim, fill_count = 0

    K, N = (<object> values).shape

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    for j in range(K):
        fill_count = 0
        val = values[j, N - 1]
        for i in range(N - 1, -1 , -1):
            if mask[j, i]:
                if fill_count >= lim:
                    continue
                fill_count += 1
                values[j, i] = val
            else:
                fill_count = 0
                val = values[j, i]
"""

backfill_1d_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_inplace_%(name)s(ndarray[%(c_type)s] values,
                              ndarray[uint8_t, cast=True] mask,
                              limit=None):
    cdef Py_ssize_t i, N
    cdef %(c_type)s val
    cdef int lim, fill_count = 0

    N = len(values)

    if limit is None:
        lim = N
    else:
        if limit < 0:
            raise ValueError('Limit must be non-negative')
        lim = limit

    val = values[N - 1]
    for i in range(N - 1, -1 , -1):
        if mask[i]:
            if fill_count >= lim:
                continue
            fill_count += 1
            values[i] = val
        else:
            fill_count = 0
            val = values[i]
"""

is_monotonic_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_%(name)s(ndarray[%(c_type)s] arr):
    '''
    Returns
    -------
    is_monotonic, is_unique
    '''
    cdef:
        Py_ssize_t i, n
        %(c_type)s prev, cur
        bint is_unique = 1

    n = len(arr)

    if n < 2:
        return True, True

    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if cur < prev:
            return False, None
        elif cur == prev:
            is_unique = 0
        prev = cur
    return True, is_unique
"""

map_indices_template = """@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_%(name)s(ndarray[%(c_type)s] index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i in range(length):
        result[index[i]] = i

    return result

"""

groupby_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_%(name)s(ndarray[%(c_type)s] index, ndarray labels):
    cdef dict result = {}
    cdef Py_ssize_t i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    for i in range(length):
        key = util.get_value_1d(labels, i)

        if _checknull(key):
            continue

        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

"""

arrmap_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_%(name)s(ndarray[%(c_type)s] index, object func):
    cdef Py_ssize_t length = index.shape[0]
    cdef Py_ssize_t i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    for i in range(length):
        result[i] = func(index[i])

    return maybe_convert_objects(result)

"""

#----------------------------------------------------------------------
# Joins on ordered, unique indices

# right might contain non-unique values

left_join_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_%(name)s(ndarray[%(c_type)s] left,
                             ndarray[%(c_type)s] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int32_t] indexer
        %(c_type)s lval, rval

    i = 0
    j = 0
    nleft = len(left)
    nright = len(right)

    indexer = np.empty(nleft, dtype=np.int32)
    while True:
        if i == nleft:
            break

        if j == nright:
            indexer[i] = -1
            i += 1
            continue

        rval = right[j]

        while i < nleft - 1 and left[i] == rval:
            indexer[i] = j
            i += 1

        if left[i] == right[j]:
            indexer[i] = j
            i += 1
            while i < nleft - 1 and left[i] == rval:
                indexer[i] = j
                i += 1
            j += 1
        elif left[i] > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer

"""

inner_join_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_%(name)s(ndarray[%(c_type)s] left,
                              ndarray[%(c_type)s] right):
    '''
    Two-pass algorithm?
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        %(c_type)s lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[%(c_type)s] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    while True:
        if i == nleft or j == nright:
             break
        else:
            lval = left[i]
            rval = right[j]
            if lval == rval:
                i += 1
                j += 1
                count += 1
            elif lval < rval:
                i += 1
            else:
                j += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int32)
    rindexer = np.empty(count, dtype=np.int32)
    result = np.empty(count, dtype=%(dtype)s)

    i = 0
    j = 0
    count = 0
    while True:
        if i == nleft or j == nright:
             break
        else:
            lval = left[i]
            rval = right[j]
            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                i += 1
                j += 1
                count += 1
            elif lval < rval:
                i += 1
            else:
                j += 1

    return result, lindexer, rindexer

"""

outer_join_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_%(name)s(ndarray[%(c_type)s] left,
                                ndarray[%(c_type)s] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        %(c_type)s lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[%(c_type)s] result

    nleft = len(left)
    nright = len(right)

    i = 0
    j = 0
    count = 0
    while True:
        if i == nleft:
            if j == nright:
                # we are done
                break
            else:
                while j < nright:
                    j += 1
                    count += 1
                break
        elif j == nright:
            while i < nleft:
                i += 1
                count += 1
            break
        else:
            if left[i] == right[j]:
                i += 1
                j += 1
            elif left[i] < right[j]:
                i += 1
            else:
                j += 1

            count += 1

    lindexer = np.empty(count, dtype=np.int32)
    rindexer = np.empty(count, dtype=np.int32)
    result = np.empty(count, dtype=%(dtype)s)

    # do it again, but populate the indexers / result

    i = 0
    j = 0
    count = 0
    while True:
        if i == nleft:
            if j == nright:
                # we are done
                break
            else:
                while j < nright:
                    lindexer[count] = -1
                    rindexer[count] = j
                    result[count] = right[j]
                    j += 1
                    count += 1
                break
        elif j == nright:
            while i < nleft:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = left[i]
                i += 1
                count += 1
            break
        else:
            lval = left[i]
            rval = right[j]
            if lval == rval:
                lindexer[count] = i
                rindexer[count] = j
                result[count] = lval
                i += 1
                j += 1
            elif lval < rval:
                lindexer[count] = i
                rindexer[count] = -1
                result[count] = lval
                i += 1
            else:
                lindexer[count] = -1
                rindexer[count] = j
                result[count] = rval
                j += 1

            count += 1

    return result, lindexer, rindexer

"""

#----------------------------------------------------------------------
# Fast "put" logic for speeding up interleaving logic

put2d_template = """
def put2d_%(name)s_%(dest_type)s(ndarray[%(c_type)s, ndim=2, cast=True] values,
                              ndarray[int32_t] indexer, Py_ssize_t loc,
                              ndarray[%(dest_type2)s] out):
    cdef:
        Py_ssize_t i, j, k

    k = len(values)
    for j from 0 <= j < k:
        i = indexer[j]
        out[i] = values[j, loc]
"""

def generate_put_functions():
    function_list = [
        ('float64', 'float64_t', 'object'),
        ('float64', 'float64_t', 'float64_t'),
        ('object', 'object', 'object'),
        ('int32', 'int32_t', 'int64_t'),
        ('int32', 'int32_t', 'float64_t'),
        ('int32', 'int32_t', 'object'),
        ('int64', 'int64_t', 'int64_t'),
        ('int64', 'int64_t', 'float64_t'),
        ('int64', 'int64_t', 'object'),
        ('bool', 'uint8_t', 'uint8_t'),
        ('bool', 'uint8_t', 'object')
    ]

    output = StringIO()
    for name, c_type, dest_type in function_list:
        func = put2d_template % {'name' : name, 'c_type' : c_type,
                                 'dest_type' : dest_type.replace('_t', ''),
                                 'dest_type2' : dest_type}
        output.write(func)
    return output.getvalue()


# name, ctype, capable of holding NA
function_list = [
    ('float64', 'float64_t', 'np.float64', True),
    ('object', 'object', 'object', True),
    ('int32', 'int32_t', 'np.int32', False),
    ('int64', 'int64_t', 'np.int64', False),
    ('bool', 'uint8_t', 'np.bool', False)
]

def generate_from_template(template, ndim=1, exclude=None):
    output = StringIO()
    for name, c_type, dtype, can_hold_na in function_list:
        if exclude is not None and name in exclude:
            continue

        func = template % {'name': name, 'c_type': c_type,
                           'dtype': dtype,
                           'raise_on_na': 'False' if can_hold_na else 'True'}
        output.write(func)
    return output.getvalue()

templates_1d = [map_indices_template,
                pad_template,
                backfill_template,
                pad_1d_template,
                backfill_1d_template,
                pad_2d_template,
                backfill_2d_template,
                take_1d_template,
                is_monotonic_template,
                groupby_template,
                arrmap_template]

nobool_1d_templates = [left_join_template,
                       outer_join_template,
                       inner_join_template]

templates_2d = [take_2d_axis0_template,
                take_2d_axis1_template]


# templates_1d_datetime = [take_1d_template]
# templates_2d_datetime = [take_2d_axis0_template,
#                          take_2d_axis1_template]


def generate_take_cython_file(path='generated.pyx'):
    with open(path, 'w') as f:
        for template in templates_1d:
            print >> f, generate_from_template(template)

        for template in templates_2d:
            print >> f, generate_from_template(template, ndim=2)

        # for template in templates_1d_datetime:
        #     print >> f, generate_from_template_datetime(template)

        # for template in templates_2d_datetime:
        #     print >> f, generate_from_template_datetime(template, ndim=2)

        for template in nobool_1d_templates:
            print >> f, generate_from_template(template, exclude=['bool'])

        # print >> f, generate_put_functions()

if __name__ == '__main__':
    generate_take_cython_file()
