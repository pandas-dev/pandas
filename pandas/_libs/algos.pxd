from pandas._libs.util cimport numeric


cdef inline Py_ssize_t swap(numeric *a, numeric *b) nogil:
    cdef:
        numeric t

    # cython doesn't allow pointer dereference so use array syntax
    t = a[0]
    a[0] = b[0]
    b[0] = t
    return 0


cdef inline numeric kth_smallest_c(numeric* arr, Py_ssize_t k, Py_ssize_t n) nogil:
    """
    Compute the kth smallest value in an array

    Parameters
    ----------
    arr: numeric* arr
        Pointer to the start of the array
    k: Py_ssize_t
    n: Number of values in arr to consider (no more than len(arr))

    Returns
    -------
    numeric
        The kth smallest value in arr
    """
    cdef:
        Py_ssize_t i, j, l, m
        numeric x

    l = 0
    m = n - 1

    while l < m:
        x = arr[k]
        i = l
        j = m

        while 1:
            while arr[i] < x: i += 1
            while x < arr[j]: j -= 1
            if i <= j:
                swap(&arr[i], &arr[j])
                i += 1; j -= 1

            if i > j: break

        if j < k: l = i
        if k < i: m = j
    return arr[k]


cdef enum TiebreakEnumType:
    TIEBREAK_AVERAGE
    TIEBREAK_MIN,
    TIEBREAK_MAX
    TIEBREAK_FIRST
    TIEBREAK_FIRST_DESCENDING
    TIEBREAK_DENSE
