
#-------------------------------------------------------------------------------
# Groupby-related functions

@cython.boundscheck(False)
def arrmap(ndarray[object, ndim=1] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object, ndim=1] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result

@cython.boundscheck(False)
def groupby_withnull_old(ndarray[object, ndim = 1] index, object keyfunc):
    cdef dict groups
    cdef int length = index.shape[0]
    cdef object idx
    cdef object curKey, key
    cdef list members

    groups = PyDict_New()

    if length != index.shape[0]:
        raise Exception('Dates and values were not the same length!')

    cdef ndarray[object, ndim=1] mapped_index = arrmap(index, keyfunc)

    cdef ndarray[npy_int8, ndim=1] null_mask = _isnullobj(mapped_index)

    bool_mask = null_mask.astype(bool)

    null_values = np.asarray(index)[bool_mask]

    if null_values.any():
        PyDict_SetItem(groups, np.NaN, null_values)

    cdef int i = 0
    idx = index[0]
    key = mapped_index[0]

    # Algorithm notes
    #   - Tries to reduce the number of calls to PyDict_GetItem,
    #   'lazily' evaluates

    while i < length:
        if not PyDict_Contains(groups, key):
            members = [idx]
            PyDict_SetItem(groups, key, members)
            i += 1
            curKey = key
            while i < length:
                if null_mask[i]:
                    i += 1
                    continue

                idx = index[i]
                key = mapped_index[i]
                if key == curKey:
                    members.append(idx)
                    i += 1
                else:
                    break
        else:
            members = <list> PyDict_GetItem(groups, key)
            members.append(idx)
            i += 1
            curKey = key
            while null_mask[i] and i < length:
                i += 1

            while i < length:
                if null_mask[i]:
                    i += 1
                    continue

                idx = index[i]
                key = mapped_index[i]
                if key == curKey:
                    members.append(idx)
                    i += 1
                else:
                    break

    return groups

@cython.boundscheck(False)
def groupby_withnull(ndarray[object, ndim = 1] index, object keyfunc):
    cdef dict groups
    cdef int length = index.shape[0]
    cdef object idx
    cdef object curKey, key
    cdef list members

    groups = PyDict_New()

    if length != index.shape[0]:
        raise Exception('Dates and values were not the same length!')

    cdef ndarray[object, ndim=1] mapped_index = arrmap(index, keyfunc)

    cdef ndarray[npy_int8, ndim=1] null_mask = _isnullobj(mapped_index)

    bool_mask = null_mask.astype(bool)

    null_values = np.asarray(index)[bool_mask]

    if null_values.any():
        PyDict_SetItem(groups, np.NaN, null_values)

    cdef int i = 0
    idx = index[0]
    key = mapped_index[0]

    # Algorithm notes
    #   - Tries to reduce the number of calls to PyDict_GetItem,
    #   'lazily' evaluates

    while i < length:
        if key not in groups:
            members = [idx]
            groups[key] = members
            i += 1
            curKey = key
            while i < length:
                if null_mask[i]:
                    i += 1
                    continue

                idx = index[i]
                key = mapped_index[i]
                if key == curKey:
                    members.append(idx)
                    i += 1
                else:
                    break
        else:
            members = <list> groups[key]
            members.append(idx)
            i += 1
            curKey = key
            while null_mask[i] and i < length:
                i += 1

            while i < length:
                if null_mask[i]:
                    i += 1
                    continue

                idx = index[i]
                key = mapped_index[i]
                if key == curKey:
                    members.append(idx)
                    i += 1
                else:
                    break

    return groups

@cython.boundscheck(False)
def groupby(ndarray[object, ndim = 1] index, object keyfunc):
    cdef dict groups
    cdef int length = index.shape[0]
    cdef object idx
    cdef object curKey, key
    cdef list members

    groups = PyDict_New()

    if length != index.shape[0]:
        raise Exception('Dates and values were not the same length!')

    cdef int i = 0
    idx = index[i]
    key = keyfunc(idx)

    # Algorithm notes
    #   - Tries to reduce the number of calls to PyDict_GetItem, 'lazily' evaluates

    while i < length:
        if not PyDict_Contains(groups, key):
            members = [idx]
            PyDict_SetItem(groups, key, members)
            i += 1
            curKey = key
            while i < length:
                idx = index[i]
                key = trycall(keyfunc, idx)
                if key == curKey:
                    members.append(idx)
                    i += 1
                else:
                    break
        else:
            members = <list> PyDict_GetItem(groups, key)
            members.append(idx)
            i += 1
            curKey = key
            while i < length:
                idx = index[i]
                key = trycall(keyfunc, idx)
                if key == curKey:
                    members.append(idx)
                    i += 1
                else:
                    break

    return groups

@cython.boundscheck(False)
def groupbyfunc(ndarray[object, ndim = 1] index,
                ndarray[npy_float64, ndim = 1] values,
                object keyfunc, object applyfunc):
    '''
    Doing this proper in Cython
    Not sure how much it will really speed things up
    '''
    cdef dict groups
    cdef int length = values.shape[0]
    cdef object idx
    cdef object curKey, key
    cdef list members, grouplist

    groups = PyDict_New()

    if length != index.shape[0]:
        raise Exception('Dates and values were not the same length!')

    cdef int i = 0
    idx = index[i]
    key = trycall(keyfunc, idx)

    # Algorithm notes
    #   - Tries to reduce the number of calls to PyDict_GetItem,
    #   'lazily' evaluates

    while i < length:
        if not PyDict_Contains(groups, key):
            members = [values[i]]
            PyDict_SetItem(groups, key, members)
            i += 1
            curKey = key
            while i < length:
                idx = index[i]
                key = trycall(keyfunc, idx)
                if key == curKey:
                    members.append(values[i])
                    i += 1
                else:
                    break
        else:
            members = <list> PyDict_GetItem(groups, key)
            members.append(values[i])
            i += 1
            curKey = key
            while i < length:
                idx = index[i]
                key = trycall(keyfunc, idx)
                if key == curKey:
                    members.append(values[i])
                    i += 1
                else:
                    break

    grouplist = PyDict_Keys(groups)

    i = 0
    length = len(grouplist)
    for i from 0 <= i < length:
        key = grouplist[i]
        members = <list> PyDict_GetItem(groups, key)
        PyDict_SetItem(groups, key, applyfunc(np.asarray(members)))

    return groups
