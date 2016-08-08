"""
Template for each `dtype` helper function for hashtable

WARNING: DO NOT edit .pxi FILE directly, .pxi is generated from .pxi.in
"""

#----------------------------------------------------------------------
# asof_join_by
#----------------------------------------------------------------------


from hashtable cimport *


def asof_join_int64_t_by_object(ndarray[int64_t] left_values,
        ndarray[int64_t] right_values,
        ndarray[object] left_by_values,
        ndarray[object] right_by_values,
        bint allow_exact_matches=1,
        tolerance=None):

    cdef:
        Py_ssize_t left_pos, right_pos, left_size, right_size, found_right_pos
        ndarray[int64_t] left_indexer, right_indexer
        bint has_tolerance = 0
        int64_t tolerance_
        PyObjectHashTable hash_table
        object by_value

    # if we are using tolerance, set our objects
    if tolerance is not None:
        has_tolerance = 1
        tolerance_ = tolerance

    left_size = len(left_values)
    right_size = len(right_values)

    left_indexer = np.empty(left_size, dtype=np.int64)
    right_indexer = np.empty(left_size, dtype=np.int64)

    hash_table = PyObjectHashTable(right_size)

    right_pos = 0
    for left_pos in range(left_size):
        # restart right_pos if it went negative in a previous iteration
        if right_pos < 0:
            right_pos = 0

        # find last position in right whose value is less than left's value
        if allow_exact_matches:
            while right_pos < right_size and\
                right_values[right_pos] <= left_values[left_pos]:
                hash_table.set_item(right_by_values[right_pos], right_pos)
                right_pos += 1
        else:
            while right_pos < right_size and\
                right_values[right_pos] < left_values[left_pos]:
                hash_table.set_item(right_by_values[right_pos], right_pos)
                right_pos += 1
        right_pos -= 1

        # save positions as the desired index
        by_value = left_by_values[left_pos]
        found_right_pos = hash_table.get_item(by_value)\
                          if by_value in hash_table else -1
        left_indexer[left_pos] = left_pos
        right_indexer[left_pos] = found_right_pos

        # if needed, verify that tolerance is met
        if has_tolerance and found_right_pos != -1:
            diff = left_values[left_pos] - right_values[found_right_pos]
            if diff > tolerance_:
                right_indexer[left_pos] = -1

    return left_indexer, right_indexer


def asof_join_double_by_object(ndarray[double] left_values,
        ndarray[double] right_values,
        ndarray[object] left_by_values,
        ndarray[object] right_by_values,
        bint allow_exact_matches=1,
        tolerance=None):

    cdef:
        Py_ssize_t left_pos, right_pos, left_size, right_size, found_right_pos
        ndarray[int64_t] left_indexer, right_indexer
        bint has_tolerance = 0
        double tolerance_
        PyObjectHashTable hash_table
        object by_value

    # if we are using tolerance, set our objects
    if tolerance is not None:
        has_tolerance = 1
        tolerance_ = tolerance

    left_size = len(left_values)
    right_size = len(right_values)

    left_indexer = np.empty(left_size, dtype=np.int64)
    right_indexer = np.empty(left_size, dtype=np.int64)

    hash_table = PyObjectHashTable(right_size)

    right_pos = 0
    for left_pos in range(left_size):
        # restart right_pos if it went negative in a previous iteration
        if right_pos < 0:
            right_pos = 0

        # find last position in right whose value is less than left's value
        if allow_exact_matches:
            while right_pos < right_size and\
                right_values[right_pos] <= left_values[left_pos]:
                hash_table.set_item(right_by_values[right_pos], right_pos)
                right_pos += 1
        else:
            while right_pos < right_size and\
                right_values[right_pos] < left_values[left_pos]:
                hash_table.set_item(right_by_values[right_pos], right_pos)
                right_pos += 1
        right_pos -= 1

        # save positions as the desired index
        by_value = left_by_values[left_pos]
        found_right_pos = hash_table.get_item(by_value)\
                          if by_value in hash_table else -1
        left_indexer[left_pos] = left_pos
        right_indexer[left_pos] = found_right_pos

        # if needed, verify that tolerance is met
        if has_tolerance and found_right_pos != -1:
            diff = left_values[left_pos] - right_values[found_right_pos]
            if diff > tolerance_:
                right_indexer[left_pos] = -1

    return left_indexer, right_indexer


def asof_join_int64_t_by_int64_t(ndarray[int64_t] left_values,
        ndarray[int64_t] right_values,
        ndarray[int64_t] left_by_values,
        ndarray[int64_t] right_by_values,
        bint allow_exact_matches=1,
        tolerance=None):

    cdef:
        Py_ssize_t left_pos, right_pos, left_size, right_size, found_right_pos
        ndarray[int64_t] left_indexer, right_indexer
        bint has_tolerance = 0
        int64_t tolerance_
        Int64HashTable hash_table
        int64_t by_value

    # if we are using tolerance, set our objects
    if tolerance is not None:
        has_tolerance = 1
        tolerance_ = tolerance

    left_size = len(left_values)
    right_size = len(right_values)

    left_indexer = np.empty(left_size, dtype=np.int64)
    right_indexer = np.empty(left_size, dtype=np.int64)

    hash_table = Int64HashTable(right_size)

    right_pos = 0
    for left_pos in range(left_size):
        # restart right_pos if it went negative in a previous iteration
        if right_pos < 0:
            right_pos = 0

        # find last position in right whose value is less than left's value
        if allow_exact_matches:
            while right_pos < right_size and\
                right_values[right_pos] <= left_values[left_pos]:
                hash_table.set_item(right_by_values[right_pos], right_pos)
                right_pos += 1
        else:
            while right_pos < right_size and\
                right_values[right_pos] < left_values[left_pos]:
                hash_table.set_item(right_by_values[right_pos], right_pos)
                right_pos += 1
        right_pos -= 1

        # save positions as the desired index
        by_value = left_by_values[left_pos]
        found_right_pos = hash_table.get_item(by_value)\
                          if by_value in hash_table else -1
        left_indexer[left_pos] = left_pos
        right_indexer[left_pos] = found_right_pos

        # if needed, verify that tolerance is met
        if has_tolerance and found_right_pos != -1:
            diff = left_values[left_pos] - right_values[found_right_pos]
            if diff > tolerance_:
                right_indexer[left_pos] = -1

    return left_indexer, right_indexer


def asof_join_double_by_int64_t(ndarray[double] left_values,
        ndarray[double] right_values,
        ndarray[int64_t] left_by_values,
        ndarray[int64_t] right_by_values,
        bint allow_exact_matches=1,
        tolerance=None):

    cdef:
        Py_ssize_t left_pos, right_pos, left_size, right_size, found_right_pos
        ndarray[int64_t] left_indexer, right_indexer
        bint has_tolerance = 0
        double tolerance_
        Int64HashTable hash_table
        int64_t by_value

    # if we are using tolerance, set our objects
    if tolerance is not None:
        has_tolerance = 1
        tolerance_ = tolerance

    left_size = len(left_values)
    right_size = len(right_values)

    left_indexer = np.empty(left_size, dtype=np.int64)
    right_indexer = np.empty(left_size, dtype=np.int64)

    hash_table = Int64HashTable(right_size)

    right_pos = 0
    for left_pos in range(left_size):
        # restart right_pos if it went negative in a previous iteration
        if right_pos < 0:
            right_pos = 0

        # find last position in right whose value is less than left's value
        if allow_exact_matches:
            while right_pos < right_size and\
                right_values[right_pos] <= left_values[left_pos]:
                hash_table.set_item(right_by_values[right_pos], right_pos)
                right_pos += 1
        else:
            while right_pos < right_size and\
                right_values[right_pos] < left_values[left_pos]:
                hash_table.set_item(right_by_values[right_pos], right_pos)
                right_pos += 1
        right_pos -= 1

        # save positions as the desired index
        by_value = left_by_values[left_pos]
        found_right_pos = hash_table.get_item(by_value)\
                          if by_value in hash_table else -1
        left_indexer[left_pos] = left_pos
        right_indexer[left_pos] = found_right_pos

        # if needed, verify that tolerance is met
        if has_tolerance and found_right_pos != -1:
            diff = left_values[left_pos] - right_values[found_right_pos]
            if diff > tolerance_:
                right_indexer[left_pos] = -1

    return left_indexer, right_indexer


#----------------------------------------------------------------------
# asof_join
#----------------------------------------------------------------------


def asof_join_int64_t(ndarray[int64_t] left_values,
        ndarray[int64_t] right_values,
        bint allow_exact_matches=1,
        tolerance=None):

    cdef:
        Py_ssize_t left_pos, right_pos, left_size, right_size
        ndarray[int64_t] left_indexer, right_indexer
        bint has_tolerance = 0
        int64_t tolerance_

    # if we are using tolerance, set our objects
    if tolerance is not None:
        has_tolerance = 1
        tolerance_ = tolerance

    left_size = len(left_values)
    right_size = len(right_values)

    left_indexer = np.empty(left_size, dtype=np.int64)
    right_indexer = np.empty(left_size, dtype=np.int64)

    right_pos = 0
    for left_pos in range(left_size):
        # restart right_pos if it went negative in a previous iteration
        if right_pos < 0:
            right_pos = 0

        # find last position in right whose value is less than left's value
        if allow_exact_matches:
            while right_pos < right_size and\
                right_values[right_pos] <= left_values[left_pos]:
                right_pos += 1
        else:
            while right_pos < right_size and\
                right_values[right_pos] < left_values[left_pos]:
                right_pos += 1
        right_pos -= 1

        # save positions as the desired index
        left_indexer[left_pos] = left_pos
        right_indexer[left_pos] = right_pos

        # if needed, verify that tolerance is met
        if has_tolerance and right_pos != -1:
            diff = left_values[left_pos] - right_values[right_pos]
            if diff > tolerance_:
                right_indexer[left_pos] = -1

    return left_indexer, right_indexer


def asof_join_double(ndarray[double] left_values,
        ndarray[double] right_values,
        bint allow_exact_matches=1,
        tolerance=None):

    cdef:
        Py_ssize_t left_pos, right_pos, left_size, right_size
        ndarray[int64_t] left_indexer, right_indexer
        bint has_tolerance = 0
        double tolerance_

    # if we are using tolerance, set our objects
    if tolerance is not None:
        has_tolerance = 1
        tolerance_ = tolerance

    left_size = len(left_values)
    right_size = len(right_values)

    left_indexer = np.empty(left_size, dtype=np.int64)
    right_indexer = np.empty(left_size, dtype=np.int64)

    right_pos = 0
    for left_pos in range(left_size):
        # restart right_pos if it went negative in a previous iteration
        if right_pos < 0:
            right_pos = 0

        # find last position in right whose value is less than left's value
        if allow_exact_matches:
            while right_pos < right_size and\
                right_values[right_pos] <= left_values[left_pos]:
                right_pos += 1
        else:
            while right_pos < right_size and\
                right_values[right_pos] < left_values[left_pos]:
                right_pos += 1
        right_pos -= 1

        # save positions as the desired index
        left_indexer[left_pos] = left_pos
        right_indexer[left_pos] = right_pos

        # if needed, verify that tolerance is met
        if has_tolerance and right_pos != -1:
            diff = left_values[left_pos] - right_values[right_pos]
            if diff > tolerance_:
                right_indexer[left_pos] = -1

    return left_indexer, right_indexer
