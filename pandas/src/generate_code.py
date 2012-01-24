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

    for i in range(n):
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

    for i in range(n):
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

    for j in range(k):
        idx = indexer[j]

        if idx == -1:
            for i in range(n):
                %(na_action)s
        else:
            for i in range(n):
                outbuf[i, j] = values[i, idx]

"""

set_na = "outbuf[i] = NaN"
set_na_2d = "outbuf[i, j] = NaN"
raise_on_na = "raise ValueError('No NA values allowed')"

merge_indexer_template = """@cython.wraparound(False)
@cython.boundscheck(False)
def merge_indexer_%(name)s(ndarray[%(c_type)s] values, dict oldMap):
    cdef Py_ssize_t i, j, length, newLength
    cdef %(c_type)s idx
    cdef ndarray[int32_t] fill_vec

    newLength = len(values)
    fill_vec = np.empty(newLength, dtype=np.int32)
    for i in range(newLength):
        idx = values[i]
        if idx in oldMap:
            fill_vec[i] = oldMap[idx]
        else:
            fill_vec[i] = -1

    return fill_vec

"""

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
def backfill_%(name)s(ndarray[%(c_type)s] oldIndex,
                      ndarray[%(c_type)s] newIndex,
                      dict oldMap, dict newMap):
    cdef Py_ssize_t i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef Py_ssize_t newPos, oldPos
    cdef %(c_type)s prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

    if oldLength == 0 or newLength == 0:
        return fill_vec

    oldPos = oldLength - 1
    newPos = newLength - 1

    if newIndex[0] > oldIndex[oldLength - 1]:
        return fill_vec

    while newPos >= 0:
        curOld = oldIndex[oldPos]

        while newIndex[newPos] > curOld:
            newPos -= 1
            if newPos < 0:
                break

        curLoc = oldMap[curOld]

        if oldPos == 0:
            if newIndex[newPos] <= curOld:
                fill_vec[:newPos + 1] = curLoc
            break
        else:
            prevOld = oldIndex[oldPos - 1]

            while newIndex[newPos] > prevOld:
                fill_vec[newPos] = curLoc

                newPos -= 1
                if newPos < 0:
                    break
        oldPos -= 1

    return fill_vec

"""

'''
Padding logic for generating fill vector

Diagram of what's going on

Old      New    Fill vector    Mask
         .                        0
         .                        0
         .                        0
A        A        0               1
         .        0               1
         .        0               1
         .        0               1
         .        0               1
         .        0               1
B        B        1               1
         .        1               1
         .        1               1
         .        1               1
C        C        2               1
'''


pad_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def pad_%(name)s(ndarray[%(c_type)s] oldIndex,
                 ndarray[%(c_type)s] newIndex,
                 dict oldMap, dict newMap):
    cdef Py_ssize_t i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef Py_ssize_t newPos, oldPos
    cdef %(c_type)s prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

    if oldLength == 0 or newLength == 0:
        return fill_vec

    oldPos = 0
    newPos = 0

    if newIndex[newLength - 1] < oldIndex[0]:
        return fill_vec

    while newPos < newLength:
        curOld = oldIndex[oldPos]

        while newIndex[newPos] < curOld:
            newPos += 1
            if newPos > newLength - 1:
                break

        curLoc = oldMap[curOld]

        if oldPos == oldLength - 1:
            if newIndex[newPos] >= curOld:
                fill_vec[newPos:] = curLoc
            break
        else:
            nextOld = oldIndex[oldPos + 1]
            done = 0

            while newIndex[newPos] < nextOld:
                fill_vec[newPos] = curLoc
                newPos += 1

                if newPos > newLength - 1:
                    done = 1
                    break

            if done:
                break

        oldPos += 1

    return fill_vec

"""

is_monotonic_template = """@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_%(name)s(ndarray[%(c_type)s] arr):
    cdef:
        Py_ssize_t i, n
        %(c_type)s prev, cur

    n = len(arr)

    if n < 2:
        return True

    prev = arr[0]
    for i in range(1, n):
        cur = arr[i]
        if cur < prev:
            return False
        prev = cur
    return True

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

        if ndim == 1:
            na_action = set_na if can_hold_na else raise_on_na
        elif ndim == 2:
            na_action = set_na_2d if can_hold_na else raise_on_na
        func = template % {'name' : name, 'c_type' : c_type,
                           'dtype' : dtype, 'na_action' : na_action}
        output.write(func)
    return output.getvalue()

templates_1d = [map_indices_template,
                merge_indexer_template,
                pad_template,
                backfill_template,
                take_1d_template,
                is_monotonic_template,
                groupby_template,
                arrmap_template]

nobool_1d_templates = [left_join_template,
                       outer_join_template,
                       inner_join_template]

templates_2d = [take_2d_axis0_template,
                take_2d_axis1_template]

def generate_take_cython_file(path='generated.pyx'):
    with open(path, 'w') as f:
        for template in templates_1d:
            print >> f, generate_from_template(template)

        for template in templates_2d:
            print >> f, generate_from_template(template, ndim=2)

        for template in nobool_1d_templates:
            print >> f, generate_from_template(template, exclude=['bool'])

        # print >> f, generate_put_functions()

if __name__ == '__main__':
    generate_take_cython_file()
