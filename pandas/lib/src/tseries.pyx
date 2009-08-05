include "numpy.pxi"
include "datetime.pxi"
include "Python.pxi"

# initialize numpy
import_array()

import numpy as np
cimport numpy as np

isnan = np.isnan
cdef double NaN = <double> np.NaN

from datetime import datetime as pydatetime

from python_dict cimport *
from numpy cimport ndarray, npy_float64, npy_int32, npy_int8, npy_float128

cimport cython

cdef inline object trycall(object func, object arg):
    cdef object result
    try:
        result = func(arg)
    except:
        raise Exception('Error calling func on index %s' % arg)
    return result

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a >= b else b

def map_indices(ndarray index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.
    
    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}
        
    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef int i, length
    cdef flatiter iter
    cdef dict result
    cdef object idx
        
    result = {}
    
    iter = PyArray_IterNew(index)

    length = PyArray_SIZE(index)
        
    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, <void *> iter.dataptr)
        result[idx] = i
        PyArray_ITER_NEXT(iter)
        
    return result

def match(ndarray A, ndarray B):
    '''
    --> match(a, b)
    
    Close equivalent of R's match function.
    
    For given input index A, find matching locations for values of A in B.
    
    Example:
    >>> b
    array([[ 0.        ,  0.26929312],
           [ 1.        ,  0.49540359],
           [ 2.        ,  0.66389941],
           [ 3.        ,  0.66235806],
           [ 4.        ,  0.97993956],
           [ 5.        ,  0.83804732],
           [ 6.        ,  0.75033074],
           [ 7.        ,  0.10250388],
           [ 8.        ,  0.66591799],
           [ 9.        ,  0.18337242]])
    >>> a
        array([1, 3, 6, 8, 4, 5, 7, 0, 2, 9])
    
    # Now with match we can realign b based on a
    
    >>> b[match(a, b[:,0]),:]
    array([[ 1.        ,  0.49540359],
           [ 3.        ,  0.66235806],
           [ 6.        ,  0.75033074],
           [ 8.        ,  0.66591799],
           [ 4.        ,  0.97993956],
           [ 5.        ,  0.83804732],
           [ 7.        ,  0.10250388],
           [ 0.        ,  0.26929312],
           [ 2.        ,  0.66389941],
           [ 9.        ,  0.18337242]])
   
    '''
    
    cdef int i, length
    cdef flatiter itera
    cdef dict bmap
    cdef double *result_data
    cdef double nan
    cdef object idx
    cdef ndarray result

    nan = <double> np.NaN
    
    bmap = map_indices(B)
        
    itera = PyArray_IterNew(A)
    length = PyArray_SIZE(A)
    
    result = <ndarray> np.empty(length, np.float64)

    result_data = <double *> result.data    
    
    for i from 0 <= i < length:
        idx = PyArray_GETITEM(A, <void *> itera.dataptr)
        if idx in bmap:
            result_data[i] = <double> bmap[idx]
        else:
            result_data[i] = nan
            
        PyArray_ITER_NEXT(itera)
    
    return result.astype(int)
    
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

    itera = PyArray_IterNew(arr)
    iteridx = PyArray_IterNew(index)

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, <void *> iteridx.dataptr)
        PyArray_ITER_NEXT(iteridx)
        if idx not in idxMap:
            result_data[i] = nan
            continue
        PyArray_ITER_GOTO1D(itera, idxMap[idx])
        result_data[i] = (<double *>(itera.dataptr))[0]

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

    itera = PyArray_IterNew(arr)
    iteridx = PyArray_IterNew(index)
    iterresult = PyArray_IterNew(result)

    cdef int res

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, <void *> iteridx.dataptr)
        PyArray_ITER_NEXT(iteridx)
        
        if idx not in idxMap:
            PyArray_SETITEM(result, <void *> iterresult.dataptr, nan)
            PyArray_ITER_NEXT(iterresult)
            continue
            
        PyArray_ITER_GOTO1D(itera, idxMap[idx])
        obj = PyArray_GETITEM(arr, <void *> itera.dataptr)        
        
        res = PyArray_SETITEM(result, <void *> iterresult.dataptr, obj)
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
    '''
    cdef int j, loc, length
    cdef object idx, value
    cdef object nan = np.NaN

    length = index.shape[0]
    cdef ndarray[object, ndim = 1] result = np.empty(length, dtype=object)
    
    loc = 0
    cdef int i = 0
    for i from 0 <= i < length:
        idx = index[i]
        if not PyDict_Contains(idxMap, idx):
            result[i] = nan
            continue
        value = arr[idxMap[idx]]
        result[i] = value
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
    iterold = PyArray_IterNew(oldIndex)

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(oldIndex, <void *> iterold.dataptr)
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
    iterold = PyArray_IterNew(oldIndex)
    iternew = PyArray_IterNew(newIndex)
    
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
        curOld = PyArray_GETITEM(oldIndex, <void *> iterold.dataptr)
        
        # Until we reach a point where we are before the curOld point
        while PyArray_GETITEM(newIndex, <void *> iternew.dataptr) > curOld:
            newPos -= 1
            if newPos < 0:
                break
            PyArray_ITER_GOTO1D(iternew, newPos)
        
        # Get the location in the old index
        curLoc = oldMap[curOld]
        
        # At the beginning of the old index
        if oldPos == 0:

            # Make sure we are before the curOld index
            if PyArray_GETITEM(newIndex, <void *> iternew.dataptr) <= curOld:
                fillVec[:newPos + 1] = curLoc
                maskVec[:newPos + 1] = 1
            
            # Exit the main loop
            break

        else:
            # Move one position back
            PyArray_ITER_GOTO1D(iterold, oldPos - 1)
            
            # Get the index there
            prevOld = PyArray_GETITEM(oldIndex, <void *> iterold.dataptr)
            
            # Until we reach the previous index
            while PyArray_GETITEM(newIndex, <void *> iternew.dataptr) > prevOld:

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
    iterold = PyArray_IterNew(oldIndex)
    iternew = PyArray_IterNew(newIndex)
    
    # Length of each index
    oldLength = PyArray_SIZE(oldIndex)
    newLength = PyArray_SIZE(newIndex)

    oldPos = 0
    newPos = 0
    while newPos < newLength:
        curOld = PyArray_GETITEM(oldIndex, <void *> iterold.dataptr)

        # At beginning, keep going until we go exceed the 
        # first OLD index in the NEW index
        while PyArray_GETITEM(newIndex, <void *> iternew.dataptr) < curOld:
            newPos += 1
            if newPos > newLength - 1:
                break
            PyArray_ITER_NEXT(iternew)

        # We got there, get the current location in the old index
        curLoc = oldMap[curOld]

        # We're at the end of the road, need to propagate this value to the end
        if oldPos == oldLength - 1:
            if PyArray_GETITEM(newIndex, <void *> iternew.dataptr) >= curOld:
                fillVec[newPos:] = curLoc
                maskVec[newPos:] = 1
            break
        else:
            # Not at the end, need to go about filling

            # Get the next index so we know when to stop propagating this value
            PyArray_ITER_NEXT(iterold)
            nextOld = PyArray_GETITEM(oldIndex, <void *> iterold.dataptr)

            done = 0
            
            # Until we reach the next OLD value in the NEW index
            while PyArray_GETITEM(newIndex, <void *> iternew.dataptr) < nextOld:
                
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
    itervals = PyArray_IterNew(values)

    for i from 0 <= i < length:
        val = PyArray_GETITEM(values, <void *> itervals.dataptr)
        if val in indexMap:
            j = indexMap[val]
            fillLocs[i] = j
            mask[i] = 1

        PyArray_ITER_NEXT(itervals)
            
    for i from 0 <= i < newLength:
        if mask[i] == 0:
            fillLocs[i] = -1

    return fillVec, maskVec.astype(np.bool)

cdef double INF = <double> np.inf
cdef double NEGINF = -INF

cdef inline _checknull(object val):
    return val is None or val != val or val == INF or val == NEGINF    

cdef ndarray _isnullobj(input):
    cdef int i, length
    cdef object val
    cdef ndarray[npy_int8, ndim=1] result    
    cdef flatiter iter 

    length = PyArray_SIZE(input)
    
    result = <ndarray> np.zeros(length, dtype=np.int8)
    
    iter= PyArray_IterNew(input)
            
    for i from 0 <= i < length:
        val = PyArray_GETITEM(input, <void *> iter.dataptr)
        
        if _checknull(val):
            result[i] = 1

        PyArray_ITER_NEXT(iter)
            
    return result
    
def isnull(input):    
    '''
    Replacement for numpy.isnan / -numpy.isfinite which is suitable
    for use on object arrays.

    Parameters
    ----------
    arr: ndarray or object value
    
    Returns
    -------
    boolean ndarray or boolean
    '''
    cdef ndarray[npy_int8, ndim=1] result
    
    if isinstance(input, np.ndarray):
        if input.dtype.kind in ('O', 'S'):
            result = _isnullobj(input)
            
            return result.astype(np.bool)
        else:
            return -np.isfinite(input)
    else:
        return _checknull(input)
    
def notnull(input):    
    '''
    Replacement for numpy.isfinite / -numpy.isnan which is suitable
    for use on object arrays.
    
    Parameters
    ----------
    arr: ndarray or object value
    
    Returns
    -------
    boolean ndarray or boolean
    '''
    if isinstance(input, np.ndarray):
        return -isnull(input)
    else:
        return not bool(_checknull(input))
    
def reindexNew(ndarray index, ndarray arr, dict idxMap):
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

    itera = PyArray_IterNew(arr)
    iteridx = PyArray_IterNew(index)

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, <void *> iteridx.dataptr)
        PyArray_ITER_NEXT(iteridx)
        if idx not in idxMap:
            result_data[i] = nan
            continue
        PyArray_ITER_GOTO1D(itera, idxMap[idx])
        result_data[i] = (<double *>(itera.dataptr))[0]

    return result
    
cdef double __add(double a, double b):
    return a + b
cdef double __sub(double a, double b):
    return a - b
cdef double __div(double a, double b):
    return a / b
cdef double __mul(double a, double b):
    return a * b
cdef double __eq(double a, double b):
    return a == b
cdef double __ne(double a, double b):
    return a != b
cdef double __lt(double a, double b):
    return a < b
cdef double __gt(double a, double b):
    return a > b
cdef double __pow(double a, double b):
    return a ** b

ctypedef double (* double_func)(double a, double b)

cdef ndarray _applyFunc(double_func func, ndarray index, object ao, 
                        object bo, dict aMap, dict bMap):
    '''
    C function taking a function pointer for quickly adding two Series objects.
    '''
    cdef ndarray A, B, result
    cdef double *result_data
    cdef int i, length
    cdef flatiter itera, iterb, iteridx
    cdef double nan
    cdef object idx
    
    # This is EXTREMELY important, otherwise you will get very 
    # undesired results
    A = PyArray_ContiguousFromAny(ao, NPY_DOUBLE, 1, 1)
    B = PyArray_ContiguousFromAny(bo, NPY_DOUBLE, 1, 1)

    nan = <double> np.NaN
    length = PyArray_SIZE(index)
    
    result = <ndarray> np.empty(length, np.float64)
    result_data = <double *>result.data
    
    itera = <flatiter> PyArray_IterNew(A)
    iterb = <flatiter> PyArray_IterNew(B)
    iteridx = PyArray_IterNew(index)
    
    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, <void *> iteridx.dataptr)
        PyArray_ITER_NEXT(iteridx)
        
        if idx not in aMap or idx not in bMap:
            result_data[i] = nan
            continue

        result_data[i] = func((<double *>A.data)[aMap[idx]], 
                            (<double *>B.data)[bMap[idx]])
                                         
    return result
    
def combineFunc(object name, ndarray index, object ao, 
                object bo, dict aMap, dict bMap):
    '''
    Combine two series (values and index maps for each passed in) using the 
    indicated function.
    '''
    if name == "__add__":
        return _applyFunc(__add, index, ao, bo, aMap, bMap)
    elif name == "__sub__":
        return _applyFunc(__sub, index, ao, bo, aMap, bMap)
    elif name == "__div__":
        return _applyFunc(__div, index, ao, bo, aMap, bMap)
    elif name == "__mul__":
        return _applyFunc(__mul, index, ao, bo, aMap, bMap)
    elif name == "__eq__":
        return _applyFunc(__eq, index, ao, bo, aMap, bMap)
    elif name == "__ne__":
        return _applyFunc(__ne, index, ao, bo, aMap, bMap)
    elif name == "__lt__":
        return _applyFunc(__lt, index, ao, bo, aMap, bMap)
    elif name == "__gt__":
        return _applyFunc(__gt, index, ao, bo, aMap, bMap)
    elif name == "__pow__":
        return _applyFunc(__pow, index, ao, bo, aMap, bMap)
    else:
        raise Exception('bad funcname requested of Cython code')
        
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

