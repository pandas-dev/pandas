
#-------------------------------------------------------------------------------
# Groupby-related functions

cdef inline _isnan(object o):
    return o != o

@cython.boundscheck(False)
def arrmap(ndarray[object, ndim=1] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object, ndim=1] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = trycall(func, index[i])

    return result

@cython.boundscheck(False)
def groupby(object index, object mapper, output=None):
    cdef dict result
    cdef ndarray[object, ndim=1] mapped_index
    cdef ndarray[object, ndim=1] index_buf
    cdef ndarray[int8_t, ndim=1] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)

    if output is None:
        result = {}
    else:
        result = <dict> output

    index_buf = np.asarray(index)
    mapped_index = arrmap(index_buf, mapper)

    mask = _isnullobj(mapped_index)
    nullkeys = index_buf[mask.astype(bool)]

    if nullkeys.any():
        result[np.NaN] = nullkeys

    for i from 0 <= i < length:
        if mask[i]:
            continue
        key = mapped_index[i]
        idx = index_buf[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.boundscheck(False)
def groupby_indices(object index, object mapper):
    cdef dict result
    cdef ndarray[object, ndim=1] mapped_index
    cdef ndarray[int8_t, ndim=1] mask
    cdef int i, length
    cdef list members, null_list
    cdef object key

    length = len(index)

    result = {}
    index = np.asarray(index)
    mapped_index = arrmap(index, mapper)

    mask = _isnullobj(mapped_index)

    if mask.astype(bool).any():
        null_list = []
        result[np.NaN] = null_list

    for i from 0 <= i < length:
        if mask[i]:
            null_list.append(i)
        key = mapped_index[i]
        if key in result:
            (<list> result[key]).append(i)
        else:
            result[key] = [i]

    return result
