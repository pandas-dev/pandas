def getFillVec(ndarray oldIndex, ndarray newIndex, dict oldMap, dict newMap,
               kind=None):

    if kind is None:
        fillVec, maskVec = getMergeVec(newIndex, oldMap)
    elif kind == 'PAD':
        fillVec, maskVec = _pad(oldIndex, newIndex, oldMap, newMap)
    elif kind == 'BACKFILL':
        fillVec, maskVec = _backfill(oldIndex, newIndex, oldMap, newMap)
    else:
        raise Exception("Don't recognize method: %s" % kind)

    return fillVec, maskVec.astype(np.bool)

@cython.wraparound(False)
def _backfill(ndarray[object] oldIndex, ndarray[object] newIndex,
              dict oldMap, dict newMap):
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
    cdef int i, j, oldLength, newLength, curLoc
    # Make empty vectors
    cdef ndarray[int32_t, ndim=1] fillVec
    cdef ndarray[int8_t, ndim=1] mask
    cdef int newPos, oldPos
    cdef object prevOld, curOld

    # Get the size
    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fillVec = np.empty(len(newIndex), dtype = np.int32)
    fillVec.fill(-1)

    mask = np.zeros(len(newIndex), dtype = np.int8)

    # Current positions
    oldPos = oldLength - 1
    newPos = newLength - 1

    # corner case, no filling possible
    if newIndex[0] > oldIndex[oldLength - 1]:
        return fillVec, mask

    while newPos >= 0:
        curOld = oldIndex[oldPos]

        # Until we reach a point where we are before the curOld point
        while newIndex[newPos] > curOld:
            newPos -= 1
            if newPos < 0:
                break

        # Get the location in the old index
        curLoc = oldMap[curOld]

        # At the beginning of the old index
        if oldPos == 0:
            # Make sure we are before the curOld index
            if newIndex[newPos] <= curOld:
                fillVec[:newPos + 1] = curLoc
                mask[:newPos + 1] = 1
            # Exit the main loop
            break
        else:
            # Get the index there
            prevOld = oldIndex[oldPos - 1]

            # Until we reach the previous index
            while newIndex[newPos] > prevOld:
                # Set the current fill location
                fillVec[newPos] = curLoc
                mask[newPos] = 1

                newPos -= 1
                if newPos < 0:
                    break

        # Move one period back
        oldPos -= 1

    return (fillVec, mask)

@cython.wraparound(False)
def _pad(ndarray[object] oldIndex, ndarray[object] newIndex,
         dict oldMap, dict newMap):
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
    cdef int i, j, oldLength, newLength, curLoc
    # Make empty vectors
    cdef ndarray[int32_t, ndim=1] fillVec
    cdef ndarray[int8_t, ndim=1] mask
    cdef int newPos, oldPos
    cdef object prevOld, curOld

    # Get the size
    oldLength = len(oldIndex)
    newLength = len(newIndex)

    fillVec = np.empty(len(newIndex), dtype = np.int32)
    fillVec.fill(-1)

    mask = np.zeros(len(newIndex), dtype = np.int8)

    oldPos = 0
    newPos = 0

    # corner case, no filling possible
    if newIndex[newLength - 1] < oldIndex[0]:
        return fillVec, mask

    while newPos < newLength:
        curOld = oldIndex[oldPos]

        # At beginning, keep going until we go exceed the
        # first OLD index in the NEW index
        while newIndex[newPos] < curOld:
            newPos += 1
            if newPos > newLength - 1:
                break

        # We got there, get the current location in the old index
        curLoc = oldMap[curOld]

        # We're at the end of the road, need to propagate this value to the end
        if oldPos == oldLength - 1:
            if newIndex[newPos] >= curOld:
                fillVec[newPos:] = curLoc
                mask[newPos:] = 1
            break
        else:
            # Not at the end, need to go about filling

            # Get the next index so we know when to stop propagating this value
            nextOld = oldIndex[oldPos + 1]

            done = 0

            # Until we reach the next OLD value in the NEW index
            while newIndex[newPos] < nextOld:
                # Use this location to fill
                fillVec[newPos] = curLoc

                # Set mask to be 1 so will not be NaN'd
                mask[newPos] = 1
                newPos += 1

                # We got to the end of the new index
                if newPos > newLength - 1:
                    done = 1
                    break

            # We got to the end of the new index
            if done:
                break

        # We already advanced the iterold pointer to the next value,
        # inc the count
        oldPos += 1

    return fillVec, mask

def pad_inplace_float64(ndarray[float64_t] values,
                        ndarray[np.uint8_t, cast=True] mask):
    '''
    mask: True if needs to be padded otherwise False

    e.g.
    pad_inplace_float64(values, isnull(values))
    '''
    cdef:
        Py_ssize_t i, n
        float64_t val

    n = len(values)
    val = NaN
    for i from 0 <= i < n:
        if mask[i]:
            values[i] = val
        else:
            val = values[i]

def get_pad_indexer(ndarray[np.uint8_t, cast=True] mask):
    '''
    mask: True if needs to be padded otherwise False

    e.g.
    pad_inplace_float64(values, isnull(values))
    '''
    cdef:
        Py_ssize_t i, n
        int32_t idx
        ndarray[int32_t] indexer

    n = len(mask)
    indexer = np.empty(n, dtype=np.int32)

    idx = 0
    for i from 0 <= i < n:
        if not mask[i]:
            idx = i
        indexer[i] = idx

    return indexer

def get_backfill_indexer(ndarray[np.uint8_t, cast=True] mask):
    '''
    mask: True if needs to be padded otherwise False

    e.g.
    pad_inplace_float64(values, isnull(values))
    '''
    cdef:
        Py_ssize_t i, n
        int32_t idx
        ndarray[int32_t] indexer

    n = len(mask)
    indexer = np.empty(n, dtype=np.int32)

    idx = n - 1
    i = n - 1
    while i >= 0:
        if not mask[i]:
            idx = i
        indexer[i] = idx
        i -= 1

    return indexer

def backfill_inplace_float64(ndarray[float64_t] values,
                             ndarray[np.uint8_t, cast=True] mask):
    '''
    mask: True if needs to be backfilled otherwise False
    '''
    cdef:
        Py_ssize_t i, n
        float64_t val

    n = len(values)
    val = NaN
    i = n - 1
    while i >= 0:
        if mask[i]:
            values[i] = val
        else:
            val = values[i]
        i -= 1

@cython.wraparound(False)
@cython.boundscheck(False)
def getMergeVec(ndarray[object] values, dict oldMap):
    cdef int i, j, length, newLength
    cdef object idx
    cdef ndarray[int32_t] fillVec
    cdef ndarray[int8_t] mask

    newLength = len(values)
    fillVec = np.empty(newLength, dtype=np.int32)
    mask = np.zeros(newLength, dtype=np.int8)
    for i from 0 <= i < newLength:
        idx = values[i]
        if idx in oldMap:
            fillVec[i] = oldMap[idx]
            mask[i] = 1

    for i from 0 <= i < newLength:
        if mask[i] == 0:
            fillVec[i] = -1

    return fillVec, mask.astype(bool)

def ordered_left_join(ndarray[object] left, ndarray[object] right):
    # cdef dict right_map = map_indices_buf(right)
    # return getMergeVec(left, right_map)
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[int32_t] indexer
        ndarray[uint8_t, cast=True] mask
        object val

    i = 0
    j = 0
    n = len(left)
    k = len(right)

    indexer = np.zeros(n, dtype=np.int32)
    mask = np.ones(n, dtype=np.bool)

    for i from 0 <= i < n:
        val = left[i]

        while j < k and right[j] < val:
            j += 1

        if j == k:
            break

        if val == right[j]:
            indexer[i] = j
            mask[i] = 0

    return indexer, mask

@cython.wraparound(False)
@cython.boundscheck(False)
def ordered_left_join_int64(ndarray[int64_t] left, ndarray[int64_t] right):
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[int32_t] indexer
        ndarray[uint8_t, cast=True] mask
        int64_t val

    i = 0
    j = 0
    n = len(left)
    k = len(right)

    indexer = np.zeros(n, dtype=np.int32)
    mask = np.ones(n, dtype=np.bool)

    for i from 0 <= i < n:
        val = left[i]

        while j < k and right[j] < val:
            j += 1

        if j == k:
            break

        if val == right[j]:
            indexer[i] = j
            mask[i] = 0

    return indexer, mask

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_2d(ndarray[int64_t] left, ndarray[int64_t] right,
                 ndarray[float64_t, ndim=2] lvalues,
                 ndarray[float64_t, ndim=2] rvalues,
                 ndarray[float64_t, ndim=2] out):
    cdef:
        Py_ssize_t i, j, k, nright, nleft, kright, kleft
        int64_t val

    nleft, kleft = (<object> lvalues).shape
    nright, kright = (<object> rvalues).shape

    j = 0
    for i from 0 <= i < nleft:
        for k from 0 <= k < kleft:
            out[i, k] = lvalues[i, k]

        val = left[i]

        while j < nright and right[j] < val:
            j += 1

        if j == nright:
            for k from kleft <= k < kleft + kright:
                out[i, k] = NaN
            continue

        if val == right[j]:
            for k from kleft <= k < kleft + kright:
                out[i, k] = rvalues[j, k - kleft]
        else:
            for k from kleft <= k < kleft + kright:
                out[i, k] = NaN

@cython.wraparound(False)
@cython.boundscheck(False)
def left_join_1d(ndarray[int64_t] left, ndarray[int64_t] right,
                 ndarray[float64_t] lvalues,
                 ndarray[float64_t] rvalues,
                 ndarray[float64_t, ndim=2] out):
    cdef:
        Py_ssize_t i, j, nright, nleft
        int64_t val

    nleft = len(lvalues)
    nright = len(rvalues)

    j = 0
    for i from 0 <= i < nleft:
        out[i, 0] = lvalues[i]

        val = left[i]

        while j < nright and right[j] < val:
            j += 1

        if j == nright:
            out[i, 1] = NaN
            continue

        if val == right[j]:
            out[i, 1] = rvalues[j]
        else:
            out[i, 1] = NaN

@cython.wraparound(False)
@cython.boundscheck(False)
def inner_join_indexer(ndarray[int64_t] left, ndarray[int64_t] right):
    '''
    Two-pass algorithm?
    '''
    cdef:
        Py_ssize_t i, j, k, nright, nleft, count
        int64_t val
        ndarray[int32_t] lindexer, rindexer
        ndarray[int64_t] result

    nleft = len(left)
    nright = len(right)

    assert(left.flags.contiguous)
    assert(right.flags.contiguous)

    cdef int64_t *lptr = <int64_t*> left.data
    cdef int64_t *rptr = <int64_t*> right.data

    i = 0
    j = 0
    count = 0
    while i < nleft:
        while j < nright and rptr[j] < lptr[i]:
            j += 1

        if j == nright:
            break

        if lptr[i] == rptr[j]:
            count += 1
            i += 1
            j += 1
        else:
            while lptr[i] < rptr[j]:
                i += 1

    # do it again now that result size is known

    lindexer = np.empty(count, dtype=np.int32)
    rindexer = np.empty(count, dtype=np.int32)
    result = np.empty(count, dtype=np.int64)

    cdef int32_t *liptr = <int32_t*> lindexer.data
    cdef int32_t *riptr = <int32_t*> rindexer.data
    cdef int64_t *resptr = <int64_t*> result.data

    i = 0
    j = 0
    count = 0
    while i < nleft:
        val = lptr[i]
        while j < nright and rptr[j] < val:
            j += 1

        if j == nright:
            break

        if val == rptr[j]:
            liptr[count] = i
            riptr[count] = j
            resptr[count] = val
            count += 1
            i += 1
            j += 1
        else:
            while lptr[i] < rptr[j]:
                i += 1

    return result, lindexer, rindexer

def _inner_join_count(ndarray[int64_t] left, ndarray[int64_t] right):
    pass

@cython.wraparound(False)
@cython.boundscheck(False)
def outer_join_indexer(ndarray[int64_t] left, ndarray[int64_t] right):
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
                result[count] = left[j]
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
def take_axis0(ndarray[float64_t, ndim=2] values,
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
def take_join_contiguous(ndarray[float64_t, ndim=2] lvalues,
                         ndarray[float64_t, ndim=2] rvalues,
                         ndarray[int32_t] lindexer,
                         ndarray[int32_t] rindexer,
                         ndarray out):
    cdef:
        Py_ssize_t i, j, rk, lk, n, lidx, ridx
        float64_t *outbuf

    assert(out.flags.contiguous)

    outbuf = <float64_t*> out.data

    n = len(lindexer)
    lk = lvalues.shape[1]
    rk = rvalues.shape[1]

    for i from 0 <= i < n:
        lidx = lindexer[i]
        ridx = rindexer[i]

        if lidx == -1:
            for j from 0 <= j < lk:
                outbuf[0] = NaN
                outbuf = outbuf + 1
        else:
            for j from 0 <= j < lk:
                outbuf[0] = lvalues[lidx, j]
                outbuf = outbuf + 1

        if lidx == -1:
            for j from 0 <= j < rk:
                outbuf[0] = NaN
                outbuf = outbuf + 1
        else:
            for j from 0 <= j < rk:
                outbuf[0] = rvalues[ridx, j]
                outbuf = outbuf + 1

@cython.wraparound(False)
@cython.boundscheck(False)
def take_axis1(ndarray[float64_t, ndim=2] values,
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
def take_1d(ndarray[float64_t] values, ndarray[int32_t] indexer,
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

def ordered_put_indexer(ndarray[int64_t] left, ndarray[int64_t] right,
                        ndarray[float64_t, ndim=2] lvalues,
                        ndarray[float64_t, ndim=2] rvalues,
                        ndarray[float64_t, ndim=2] out):
    pass

def ordered_outer_join(ndarray[int64_t] left, ndarray[int64_t] right):
    cdef:
        Py_ssize_t i, j, k, nright, nleft, kright, kleft
        int64_t val
    pass


def ordered_inner_join(ndarray[object] left, ndarray[object] right):
    pass

