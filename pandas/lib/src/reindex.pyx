
def reindex(ndarray index, ndarray arr, dict idxMap):
    '''
    Using the provided new index, a given array, and a mapping of index-value
    correpondences in the value array, return a new ndarray conforming to
    the new index.

    This is significantly faster than doing it in pure Python.
    '''
    cdef ndarray result
    cdef double *result_data
    cdef int i, length
    cdef flatiter itera, iteridx
    cdef double nan
    cdef object idx

    nan = <double> np.NaN

    length = PyArray_SIZE(index)

    result = <ndarray> np.empty(length, np.float64)

    result_data = <double *> result.data

    itera = <flatiter> PyArray_IterNew(arr)
    iteridx = <flatiter> PyArray_IterNew(index)

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, PyArray_ITER_DATA(iteridx))
        PyArray_ITER_NEXT(iteridx)
        if idx not in idxMap:
            result_data[i] = nan
            continue
        PyArray_ITER_GOTO1D(itera, idxMap[idx])
        result_data[i] = (<double *>PyArray_ITER_DATA(itera))[0]

    return result

def reindexObj(ndarray index, ndarray arr, dict idxMap):
    '''
    Using the provided new index, a given array, and a mapping of index-value
    correpondences in the value array, return a new ndarray conforming to
    the new index.

    This is significantly faster than doing it in pure Python.
    '''
    cdef ndarray result
    cdef int i, length
    cdef flatiter itera, iteridx, iterresult
    cdef object idx, nan, obj

    nan = np.NaN
    length = PyArray_SIZE(index)

    result = <ndarray> np.empty(length, dtype=np.object_)

    itera = <flatiter> PyArray_IterNew(arr)
    iteridx = <flatiter> PyArray_IterNew(index)
    iterresult = <flatiter> PyArray_IterNew(result)

    cdef int res

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, PyArray_ITER_DATA(iteridx))
        PyArray_ITER_NEXT(iteridx)

        if idx not in idxMap:
            PyArray_SETITEM(result, PyArray_ITER_DATA(iterresult), nan)
            PyArray_ITER_NEXT(iterresult)
            continue

        PyArray_ITER_GOTO1D(itera, idxMap[idx])
        obj = PyArray_GETITEM(arr, PyArray_ITER_DATA(itera))

        res = PyArray_SETITEM(result, PyArray_ITER_DATA(iterresult), obj)
        PyArray_ITER_NEXT(iterresult)

    return result

@cython.boundscheck(False)
def reindexObject(ndarray[object, ndim=1] index,
                  ndarray[object, ndim=1] arr,
                  dict idxMap):
    '''
    Using the provided new index, a given array, and a mapping of index-value
    correpondences in the value array, return a new ndarray conforming to
    the new index.

    Returns
    -------
    ndarray
    '''
    cdef ndarray[object, ndim = 1] result
    cdef int i, loc, length
    cdef object idx, value
    cdef object nan = np.NaN

    length = index.shape[0]
    result = np.empty(length, dtype=object)

    for i from 0 <= i < length:
        idx = index[i]
        if idx not in idxMap:
            result[i] = nan
            continue

        result[i] = arr[idxMap[idx]]

    return result

cdef tuple _nofill(ndarray oldIndex, ndarray newIndex, dict oldMap, dict newMap):
    cdef int *fillLocs
    cdef char *mask







    cdef int i, j, length, newLength

    cdef flatiter iterold
    cdef object idx
    cdef ndarray fillVec
    cdef ndarray maskVec

    fillVec = <ndarray> np.empty(len(newIndex), dtype = np.int32)
    maskVec = <ndarray> np.zeros(len(newIndex), dtype = np.int8)

    fillLocs = <int *> fillVec.data
    mask = <char *> maskVec.data

    newLength = PyArray_SIZE(fillVec)

    length = PyArray_SIZE(oldIndex)
    iterold = <flatiter> PyArray_IterNew(oldIndex)

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(oldIndex, PyArray_ITER_DATA(iterold))
        if i < length - 1:
           PyArray_ITER_NEXT(iterold)
        if idx in newMap:
            j = newMap[idx]
            fillLocs[j] = i
            mask[j] = 1

    for i from 0 <= i < newLength:
        if mask[i] == 0:
            fillLocs[i] = -1

    return fillVec, maskVec

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

    if kind == '':
        fillVec, maskVec = _nofill(oldIndex, newIndex, oldMap, newMap)
    elif kind == 'PAD':
        fillVec, maskVec = _pad(oldIndex, newIndex, oldMap, newMap)
    elif kind == 'BACKFILL':
        fillVec, maskVec = _backfill(oldIndex, newIndex, oldMap, newMap)
    else:
        raise Exception("Don't recognize fillMethod: %s" % fillMethod)

    return fillVec, maskVec.astype(np.bool)

def getMergeVec(ndarray values, dict indexMap):
    cdef int *fillLocs
    cdef char *mask
    cdef int i, j, length

    cdef flatiter itervals
    cdef object val
    cdef ndarray fillVec
    cdef ndarray maskVec

    cdef int newLength = len(values)

    fillVec = <ndarray> np.empty(newLength, dtype = np.int32)
    maskVec = <ndarray> np.zeros(newLength, dtype = np.int8)

    fillLocs = <int *> fillVec.data
    mask = <char *> maskVec.data

    length = PyArray_SIZE(values)
    itervals = <flatiter> PyArray_IterNew(values)

    for i from 0 <= i < length:
        val = PyArray_GETITEM(values, PyArray_ITER_DATA(itervals))
        if val in indexMap:
            j = indexMap[val]
            fillLocs[i] = j
            mask[i] = 1

        PyArray_ITER_NEXT(itervals)

    for i from 0 <= i < newLength:
        if mask[i] == 0:
            fillLocs[i] = -1

    return fillVec, maskVec.astype(np.bool)

