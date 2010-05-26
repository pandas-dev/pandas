def getFillVec(ndarray oldIndex, ndarray newIndex, dict oldMap, dict newMap,
               object kind):

    if kind is None:
        fillVec, maskVec = getMergeVec(newIndex, oldMap)
    elif kind == 'PAD':
        fillVec, maskVec = _pad(oldIndex, newIndex, oldMap, newMap)
    elif kind == 'BACKFILL':
        fillVec, maskVec = _backfill(oldIndex, newIndex, oldMap, newMap)
    else:
        raise Exception("Don't recognize fillMethod: %s" % kind)

    return fillVec, maskVec.astype(np.bool)

@cython.wraparound(False)
def _backfill(ndarray[object, ndim=1] oldIndex,
              ndarray[object, ndim=1] newIndex,
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
def _pad(ndarray[object, ndim=1] oldIndex,
         ndarray[object, ndim=1] newIndex,
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

@cython.boundscheck(False)
def getMergeVec(ndarray values, dict oldMap):
    cdef int i, j, length, newLength

    cdef flatiter iternew
    cdef object idx
    cdef ndarray[int32_t, ndim=1] fillVec
    cdef ndarray[int8_t, ndim=1] mask

    newLength = len(values)
    fillVec = np.empty(newLength, dtype=np.int32)
    mask = np.zeros(newLength, dtype=np.int8)

    iternew = <flatiter> PyArray_IterNew(values)

    for i from 0 <= i < newLength:
        idx = PyArray_GETITEM(values, PyArray_ITER_DATA(iternew))

        if idx in oldMap:
            fillVec[i] = oldMap[idx]
            mask[i] = 1

        PyArray_ITER_NEXT(iternew)

    for i from 0 <= i < newLength:
        if mask[i] == 0:
            fillVec[i] = -1

    return fillVec, mask.astype(bool)

