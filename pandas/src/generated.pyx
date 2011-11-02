@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_float64(ndarray[float64_t] index):
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

    for i from 0 <= i < length:
        result[index[i]] = i

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_object(ndarray[object] index):
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

    for i from 0 <= i < length:
        result[index[i]] = i

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_int32(ndarray[int32_t] index):
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

    for i from 0 <= i < length:
        result[index[i]] = i

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_int64(ndarray[int64_t] index):
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

    for i from 0 <= i < length:
        result[index[i]] = i

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef map_indices_bool(ndarray[uint8_t] index):
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

    for i from 0 <= i < length:
        result[index[i]] = i

    return result


@cython.wraparound(False)
@cython.boundscheck(False)
def merge_indexer_float64(ndarray[float64_t] values, dict oldMap):
    cdef int i, j, length, newLength
    cdef float64_t idx
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

@cython.wraparound(False)
@cython.boundscheck(False)
def merge_indexer_object(ndarray[object] values, dict oldMap):
    cdef int i, j, length, newLength
    cdef object idx
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

@cython.wraparound(False)
@cython.boundscheck(False)
def merge_indexer_int32(ndarray[int32_t] values, dict oldMap):
    cdef int i, j, length, newLength
    cdef int32_t idx
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

@cython.wraparound(False)
@cython.boundscheck(False)
def merge_indexer_int64(ndarray[int64_t] values, dict oldMap):
    cdef int i, j, length, newLength
    cdef int64_t idx
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

@cython.wraparound(False)
@cython.boundscheck(False)
def merge_indexer_bool(ndarray[uint8_t] values, dict oldMap):
    cdef int i, j, length, newLength
    cdef uint8_t idx
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


@cython.boundscheck(False)
@cython.wraparound(False)
def pad_float64(ndarray[float64_t] oldIndex,
                 ndarray[float64_t] newIndex,
                 dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef float64_t prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def pad_object(ndarray[object] oldIndex,
                 ndarray[object] newIndex,
                 dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef object prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def pad_int32(ndarray[int32_t] oldIndex,
                 ndarray[int32_t] newIndex,
                 dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef int32_t prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def pad_int64(ndarray[int64_t] oldIndex,
                 ndarray[int64_t] newIndex,
                 dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef int64_t prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def pad_bool(ndarray[uint8_t] oldIndex,
                 ndarray[uint8_t] newIndex,
                 dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef uint8_t prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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


@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_float64(ndarray[float64_t] oldIndex,
                      ndarray[float64_t] newIndex,
                      dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef float64_t prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_object(ndarray[object] oldIndex,
                      ndarray[object] newIndex,
                      dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef object prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_int32(ndarray[int32_t] oldIndex,
                      ndarray[int32_t] newIndex,
                      dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef int32_t prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_int64(ndarray[int64_t] oldIndex,
                      ndarray[int64_t] newIndex,
                      dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef int64_t prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def backfill_bool(ndarray[uint8_t] oldIndex,
                      ndarray[uint8_t] newIndex,
                      dict oldMap, dict newMap):
    cdef int i, j, oldLength, newLength, curLoc
    cdef ndarray[int32_t, ndim=1] fill_vec
    cdef int newPos, oldPos
    cdef uint8_t prevOld, curOld

    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fill_vec = np.empty(len(newIndex), dtype = np.int32)
    fill_vec.fill(-1)

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


@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_float64(ndarray[float64_t] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[float64_t] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            outbuf[i] = NaN
        else:
            outbuf[i] = values[idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_object(ndarray[object] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[object] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            outbuf[i] = NaN
        else:
            outbuf[i] = values[idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_int32(ndarray[int32_t] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[int32_t] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            raise ValueError('No NA values allowed')
        else:
            outbuf[i] = values[idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_int64(ndarray[int64_t] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[int64_t] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            raise ValueError('No NA values allowed')
        else:
            outbuf[i] = values[idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_1d_bool(ndarray[uint8_t] values, ndarray[int32_t] indexer,
                     out=None):
    cdef:
        Py_ssize_t i, n, idx
        ndarray[uint8_t] outbuf

    n = len(indexer)

    if out is None:
        outbuf = np.empty(n, dtype=values.dtype)
    else:
        outbuf = out

    for i from 0 <= i < n:
        idx = indexer[i]
        if idx == -1:
            raise ValueError('No NA values allowed')
        else:
            outbuf[i] = values[idx]


@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_float64(ndarray[float64_t] arr):
    cdef:
        Py_ssize_t i, n
        float64_t prev, cur

    n = len(arr)

    if n < 2:
        return True

    prev = arr[0]
    for i from 1 <= i < n:
        cur = arr[i]
        if cur < prev:
            return False
        prev = cur
    return True

@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_object(ndarray[object] arr):
    cdef:
        Py_ssize_t i, n
        object prev, cur

    n = len(arr)

    if n < 2:
        return True

    prev = arr[0]
    for i from 1 <= i < n:
        cur = arr[i]
        if cur < prev:
            return False
        prev = cur
    return True

@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_int32(ndarray[int32_t] arr):
    cdef:
        Py_ssize_t i, n
        int32_t prev, cur

    n = len(arr)

    if n < 2:
        return True

    prev = arr[0]
    for i from 1 <= i < n:
        cur = arr[i]
        if cur < prev:
            return False
        prev = cur
    return True

@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_int64(ndarray[int64_t] arr):
    cdef:
        Py_ssize_t i, n
        int64_t prev, cur

    n = len(arr)

    if n < 2:
        return True

    prev = arr[0]
    for i from 1 <= i < n:
        cur = arr[i]
        if cur < prev:
            return False
        prev = cur
    return True

@cython.boundscheck(False)
@cython.wraparound(False)
def is_monotonic_bool(ndarray[uint8_t] arr):
    cdef:
        Py_ssize_t i, n
        uint8_t prev, cur

    n = len(arr)

    if n < 2:
        return True

    prev = arr[0]
    for i from 1 <= i < n:
        cur = arr[i]
        if cur < prev:
            return False
        prev = cur
    return True


@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_float64(ndarray[float64_t] index, ndarray[object] labels):
    cdef dict result = {}
    cdef ndarray[uint8_t] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)
    mask = isnullobj(labels).view(np.uint8)

    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = labels[i]
        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_object(ndarray[object] index, ndarray[object] labels):
    cdef dict result = {}
    cdef ndarray[uint8_t] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)
    mask = isnullobj(labels).view(np.uint8)

    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = labels[i]
        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_int32(ndarray[int32_t] index, ndarray[object] labels):
    cdef dict result = {}
    cdef ndarray[uint8_t] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)
    mask = isnullobj(labels).view(np.uint8)

    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = labels[i]
        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_int64(ndarray[int64_t] index, ndarray[object] labels):
    cdef dict result = {}
    cdef ndarray[uint8_t] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)
    mask = isnullobj(labels).view(np.uint8)

    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = labels[i]
        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def groupby_bool(ndarray[uint8_t] index, ndarray[object] labels):
    cdef dict result = {}
    cdef ndarray[uint8_t] mask
    cdef int i, length
    cdef list members
    cdef object idx, key

    length = len(index)
    mask = isnullobj(labels).view(np.uint8)

    for i from 0 <= i < length:
        if mask[i]:
            continue

        key = labels[i]
        idx = index[i]
        if key in result:
            members = result[key]
            members.append(idx)
        else:
            result[key] = [idx]

    return result


@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_float64(ndarray[float64_t] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_object(ndarray[object] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_int32(ndarray[int32_t] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_int64(ndarray[int64_t] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def arrmap_bool(ndarray[uint8_t] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_float64(ndarray[float64_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[float64_t, ndim=2] outbuf

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
                outbuf[i, j] = NaN
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_object(ndarray[object, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[object, ndim=2] outbuf

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
                outbuf[i, j] = NaN
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int32(ndarray[int32_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int32_t, ndim=2] outbuf

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
                raise ValueError('No NA values allowed')
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_int64(ndarray[int64_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t, ndim=2] outbuf

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
                raise ValueError('No NA values allowed')
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis0_bool(ndarray[uint8_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[uint8_t, ndim=2] outbuf

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
                raise ValueError('No NA values allowed')
        else:
            for j from 0 <= j < k:
                outbuf[i, j] = values[idx, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_float64(ndarray[float64_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[float64_t, ndim=2] outbuf

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
                outbuf[i, j] = NaN
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_object(ndarray[object, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[object, ndim=2] outbuf

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
                outbuf[i, j] = NaN
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int32(ndarray[int32_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int32_t, ndim=2] outbuf

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
                raise ValueError('No NA values allowed')
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_int64(ndarray[int64_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[int64_t, ndim=2] outbuf

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
                raise ValueError('No NA values allowed')
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]

@cython.wraparound(False)
@cython.boundscheck(False)
def take_2d_axis1_bool(ndarray[uint8_t, ndim=2] values,
                           ndarray[int32_t] indexer,
                           out=None):
    cdef:
        Py_ssize_t i, j, k, n, idx
        ndarray[uint8_t, ndim=2] outbuf

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
                raise ValueError('No NA values allowed')
        else:
            for i from 0 <= i < n:
                outbuf[i, j] = values[i, idx]


@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_float64(ndarray[float64_t] left,
                             ndarray[float64_t] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int32_t] indexer
        float64_t lval, rval

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

        lval = left[i]
        rval = right[j]

        if lval == right[j]:
            indexer[i] = j
            i += 1
            j += 1
        elif lval > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_object(ndarray[object] left,
                             ndarray[object] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int32_t] indexer
        object lval, rval

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

        lval = left[i]
        rval = right[j]

        if lval == right[j]:
            indexer[i] = j
            i += 1
            j += 1
        elif lval > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_int32(ndarray[int32_t] left,
                             ndarray[int32_t] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int32_t] indexer
        int32_t lval, rval

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

        lval = left[i]
        rval = right[j]

        if lval == right[j]:
            indexer[i] = j
            i += 1
            j += 1
        elif lval > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_indexer_int64(ndarray[int64_t] left,
                             ndarray[int64_t] right):
    cdef:
        Py_ssize_t i, j, nleft, nright
        ndarray[int32_t] indexer
        int64_t lval, rval

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

        lval = left[i]
        rval = right[j]

        if lval == right[j]:
            indexer[i] = j
            i += 1
            j += 1
        elif lval > rval:
            indexer[i] = -1
            j += 1
        else:
            indexer[i] = -1
            i += 1
    return indexer


@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_float64(ndarray[float64_t] left,
                                ndarray[float64_t] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        float64_t lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[float64_t] result

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
    result = np.empty(count, dtype=np.float64)

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

@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_object(ndarray[object] left,
                                ndarray[object] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        object lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[object] result

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
    result = np.empty(count, dtype=object)

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

@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_int32(ndarray[int32_t] left,
                                ndarray[int32_t] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        int32_t lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[int32_t] result

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
    result = np.empty(count, dtype=np.int32)

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

@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer_int64(ndarray[int64_t] left,
                                ndarray[int64_t] right):
    cdef:
        Py_ssize_t i, j, nright, nleft, count
        int64_t lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[int64_t] result

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
    result = np.empty(count, dtype=np.int64)

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


@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_float64(ndarray[float64_t] left,
                              ndarray[float64_t] right):
    '''
    Two-pass algorithm?
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        float64_t lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[float64_t] result

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
    result = np.empty(count, dtype=np.float64)

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

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_object(ndarray[object] left,
                              ndarray[object] right):
    '''
    Two-pass algorithm?
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        object lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[object] result

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
    result = np.empty(count, dtype=object)

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

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_int32(ndarray[int32_t] left,
                              ndarray[int32_t] right):
    '''
    Two-pass algorithm?
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        int32_t lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[int32_t] result

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
    result = np.empty(count, dtype=np.int32)

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

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer_int64(ndarray[int64_t] left,
                              ndarray[int64_t] right):
    '''
    Two-pass algorithm?
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        int64_t lval, rval
        ndarray[int32_t] lindexer, rindexer
        ndarray[int64_t] result

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
    result = np.empty(count, dtype=np.int64)

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


