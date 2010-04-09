
cdef tuple _backfill(ndarray oldIndex, ndarray newIndex, dict oldMap, dict newMap):
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
    cdef ndarray fillVec
    cdef ndarray maskVec
    fillVec = <ndarray> np.empty(len(newIndex), dtype = np.int32)
    maskVec = <ndarray> np.zeros(len(newIndex), dtype = np.int8)

    # Get references
    cdef int *fillLocs
    cdef char *mask
    fillLocs = <int *> fillVec.data
    mask = <char *> maskVec.data

    # Create the iterators
    cdef flatiter iterold, iternew
    iterold = <flatiter> PyArray_IterNew(oldIndex)
    iternew = <flatiter> PyArray_IterNew(newIndex)

    # Get the size
    oldLength = PyArray_SIZE(oldIndex)
    newLength = PyArray_SIZE(newIndex)

    # Current positions
    cdef int newPos, oldPos
    oldPos = oldLength - 1
    newPos = newLength - 1

    # References holding indices
    cdef object prevOld, curOld

    while newPos >= 0:
        # Move to the current position
        PyArray_ITER_GOTO1D(iternew, newPos)
        PyArray_ITER_GOTO1D(iterold, oldPos)

        # Get the current index
        curOld = PyArray_GETITEM(oldIndex, PyArray_ITER_DATA(iterold))

        # Until we reach a point where we are before the curOld point
        while PyArray_GETITEM(newIndex, PyArray_ITER_DATA(iternew)) > curOld:
            newPos -= 1
            if newPos < 0:
                break
            PyArray_ITER_GOTO1D(iternew, newPos)

        # Get the location in the old index
        curLoc = oldMap[curOld]

        # At the beginning of the old index
        if oldPos == 0:

            # Make sure we are before the curOld index
            if PyArray_GETITEM(newIndex, PyArray_ITER_DATA(iternew)) <= curOld:
                fillVec[:newPos + 1] = curLoc
                maskVec[:newPos + 1] = 1

            # Exit the main loop
            break

        else:
            # Move one position back
            PyArray_ITER_GOTO1D(iterold, oldPos - 1)

            # Get the index there
            prevOld = PyArray_GETITEM(oldIndex, PyArray_ITER_DATA(iterold))

            # Until we reach the previous index
            while PyArray_GETITEM(newIndex, PyArray_ITER_DATA(iternew)) > prevOld:

                # Set the current fill location
                fillLocs[newPos] = curLoc
                mask[newPos] = 1

                newPos -= 1
                if newPos < 0:
                    break

                # Move the iterator back
                PyArray_ITER_GOTO1D(iternew, newPos)

        # Move one period back
        oldPos -= 1

    for i from 0 <= i < newLength:
        if mask[i] == 0:
            # Fill from some generic location
            fillLocs[i] = -1

    return (fillVec, maskVec)

cdef tuple _pad(ndarray oldIndex, ndarray newIndex, dict oldMap, dict newMap):
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

    # Declare variables
    cdef ndarray fillVec
    cdef ndarray maskVec
    cdef int *fillLocs
    cdef char *mask
    cdef int i, j, oldLength, newLength, curLoc, newPos, oldPos
    cdef flatiter iterold, iternew
    cdef object nextOld, curOld
    cdef char done

    # Make empty fill vector and mask vector, cast to ndarray
    fillVec = <ndarray> np.empty(len(newIndex), dtype = np.int32)
    maskVec = <ndarray> np.zeros(len(newIndex), dtype = np.int8)

    # Get reference to the arrays inside
    fillLocs = <int *> fillVec.data
    mask = <char *> maskVec.data

    # Create simple ndarray iterators using C API
    iterold = <flatiter> PyArray_IterNew(oldIndex)
    iternew = <flatiter> PyArray_IterNew(newIndex)

    # Length of each index
    oldLength = PyArray_SIZE(oldIndex)
    newLength = PyArray_SIZE(newIndex)

    oldPos = 0
    newPos = 0
    while newPos < newLength:
        curOld = PyArray_GETITEM(oldIndex, PyArray_ITER_DATA(iterold))

        # At beginning, keep going until we go exceed the
        # first OLD index in the NEW index
        while PyArray_GETITEM(newIndex, PyArray_ITER_DATA(iternew)) < curOld:
            newPos += 1
            if newPos > newLength - 1:
                break
            PyArray_ITER_NEXT(iternew)

        # We got there, get the current location in the old index
        curLoc = oldMap[curOld]

        # We're at the end of the road, need to propagate this value to the end
        if oldPos == oldLength - 1:
            if PyArray_GETITEM(newIndex, PyArray_ITER_DATA(iternew)) >= curOld:
                fillVec[newPos:] = curLoc
                maskVec[newPos:] = 1
            break
        else:
            # Not at the end, need to go about filling

            # Get the next index so we know when to stop propagating this value
            PyArray_ITER_NEXT(iterold)
            nextOld = PyArray_GETITEM(oldIndex, PyArray_ITER_DATA(iterold))

            done = 0

            # Until we reach the next OLD value in the NEW index
            while PyArray_GETITEM(newIndex, PyArray_ITER_DATA(iternew)) < nextOld:

                # Use this location to fill
                fillLocs[newPos] = curLoc

                # Set mask to be 1 so will not be NaN'd
                mask[newPos] = 1
                newPos += 1

                # We got to the end of the new index
                if newPos > newLength - 1:
                    done = 1
                    break

                # Advance the pointer
                PyArray_ITER_NEXT(iternew)

            # We got to the end of the new index
            if done:
                break

        # We already advanced the iterold pointer to the next value,
        # inc the count
        oldPos += 1

    # Places where the mask is 0, fill with an arbitrary value
    # (will be NA'd out)
    for i from 0 <= i < newLength:
        if mask[i] == 0:
            fillLocs[i] = -1

    return fillVec, maskVec

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

@cython.boundscheck(False)
def getMergeVec(ndarray values, dict oldMap):
    cdef int i, j, length, newLength

    cdef flatiter iternew
    cdef object idx
    cdef ndarray[int32_t, ndim=1] fillLocs
    cdef ndarray[int8_t, ndim=1] mask

    newLength = len(values)
    fillLocs = np.empty(newLength, dtype=np.int32)
    mask = np.zeros(newLength, dtype=np.int8)

    iternew = <flatiter> PyArray_IterNew(values)

    for i from 0 <= i < newLength:
        idx = PyArray_GETITEM(values, PyArray_ITER_DATA(iternew))

        if idx in oldMap:
            fillLocs[i] = oldMap[idx]
            mask[i] = 1

        PyArray_ITER_NEXT(iternew)

    for i from 0 <= i < newLength:
        if mask[i] == 0:
            fillLocs[i] = -1

    return fillLocs, mask.astype(bool)

