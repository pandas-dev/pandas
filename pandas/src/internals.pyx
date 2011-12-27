def get_reverse_indexer(ndarray[int32_t] indexer, Py_ssize_t length):
    cdef:
        Py_ssize_t i
        ndarray[int32_t] rev_indexer
        int32_t idx

    rev_indexer = np.empty(length, dtype='i4')
    for i in range(len(indexer)):
        idx = indexer[i]
        if idx != -1:
            rev_indexer[idx] = i

    return rev_indexer
