"""
Template for each `dtype` helper function using groupby

WARNING: DO NOT edit .pxi FILE directly, .pxi is generated from .pxi.in
"""

cdef extern from "numpy/npy_math.h":
    double NAN "NPY_NAN"
_int64_max = np.iinfo(np.int64).max

#----------------------------------------------------------------------
# group_add, group_prod, group_var, group_mean, group_ohlc
#----------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
def group_add_float64(ndarray[float64_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float64_t, ndim=2] values,
                       ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    with nogil:

        if K > 1:

            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val:
                        nobs[lab, j] += 1
                        sumx[lab, j] += val

        else:

            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val:
                    nobs[lab, 0] += 1
                    sumx[lab, 0] += val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = sumx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_prod_float64(ndarray[float64_t, ndim=2] out,
                        ndarray[int64_t] counts,
                        ndarray[float64_t, ndim=2] values,
                        ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float64_t val, count
        ndarray[float64_t, ndim=2] prodx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val:
                        nobs[lab, j] += 1
                        prodx[lab, j] *= val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val:
                    nobs[lab, 0] += 1
                    prodx[lab, 0] *= val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = prodx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def group_var_float64(ndarray[float64_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float64_t, ndim=2] values,
                       ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float64_t val, ct, oldmean
        ndarray[float64_t, ndim=2] nobs, mean

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    mean = np.zeros_like(out)

    N, K = (<object> values).shape

    out[:, :] = 0.0

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1

            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    oldmean = mean[lab, j]
                    mean[lab, j] += (val - oldmean) / nobs[lab, j]
                    out[lab, j] += (val - mean[lab, j]) * (val - oldmean)

        for i in range(ncounts):
            for j in range(K):
                ct = nobs[i, j]
                if ct < 2:
                    out[i, j] = NAN
                else:
                    out[i, j] /= (ct - 1)
# add passing bin edges, instead of labels


@cython.wraparound(False)
@cython.boundscheck(False)
def group_mean_float64(ndarray[float64_t, ndim=2] out,
                        ndarray[int64_t] counts,
                        ndarray[float64_t, ndim=2] values,
                        ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float64_t val, count
        ndarray[float64_t, ndim=2] sumx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]
                    # not nan
                    if val == val:
                        nobs[lab, j] += 1
                        sumx[lab, j] += val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]
                # not nan
                if val == val:
                    nobs[lab, 0] += 1
                    sumx[lab, 0] += val

        for i in range(ncounts):
            for j in range(K):
                count = nobs[i, j]
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = sumx[i, j] / count


@cython.wraparound(False)
@cython.boundscheck(False)
def group_ohlc_float64(ndarray[float64_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float64_t, ndim=2] values,
                  ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab
        float64_t val, count
        Py_ssize_t ngroups = len(counts)

    if len(labels) == 0:
        return

    N, K = (<object> values).shape

    if out.shape[1] != 4:
        raise ValueError('Output array must have 4 columns')

    if K > 1:
        raise NotImplementedError("Argument 'values' must have only "
                                  "one dimension")
    out.fill(np.nan)

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab == -1:
                continue

            counts[lab] += 1
            val = values[i, 0]
            if val != val:
                continue

            if out[lab, 0] != out[lab, 0]:
                out[lab, 0] = out[lab, 1] = out[lab, 2] = out[lab, 3] = val
            else:
                out[lab, 1] = max(out[lab, 1], val)
                out[lab, 2] = min(out[lab, 2], val)
                out[lab, 3] = val


@cython.wraparound(False)
@cython.boundscheck(False)
def group_add_float32(ndarray[float32_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float32_t, ndim=2] values,
                       ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float32_t val, count
        ndarray[float32_t, ndim=2] sumx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    with nogil:

        if K > 1:

            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val:
                        nobs[lab, j] += 1
                        sumx[lab, j] += val

        else:

            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val:
                    nobs[lab, 0] += 1
                    sumx[lab, 0] += val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = sumx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_prod_float32(ndarray[float32_t, ndim=2] out,
                        ndarray[int64_t] counts,
                        ndarray[float32_t, ndim=2] values,
                        ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float32_t val, count
        ndarray[float32_t, ndim=2] prodx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    prodx = np.ones_like(out)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val:
                        nobs[lab, j] += 1
                        prodx[lab, j] *= val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val:
                    nobs[lab, 0] += 1
                    prodx[lab, 0] *= val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = prodx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def group_var_float32(ndarray[float32_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float32_t, ndim=2] values,
                       ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float32_t val, ct, oldmean
        ndarray[float32_t, ndim=2] nobs, mean

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    mean = np.zeros_like(out)

    N, K = (<object> values).shape

    out[:, :] = 0.0

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1

            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val:
                    nobs[lab, j] += 1
                    oldmean = mean[lab, j]
                    mean[lab, j] += (val - oldmean) / nobs[lab, j]
                    out[lab, j] += (val - mean[lab, j]) * (val - oldmean)

        for i in range(ncounts):
            for j in range(K):
                ct = nobs[i, j]
                if ct < 2:
                    out[i, j] = NAN
                else:
                    out[i, j] /= (ct - 1)
# add passing bin edges, instead of labels


@cython.wraparound(False)
@cython.boundscheck(False)
def group_mean_float32(ndarray[float32_t, ndim=2] out,
                        ndarray[int64_t] counts,
                        ndarray[float32_t, ndim=2] values,
                        ndarray[int64_t] labels):
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float32_t val, count
        ndarray[float32_t, ndim=2] sumx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)
    sumx = np.zeros_like(out)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]
                    # not nan
                    if val == val:
                        nobs[lab, j] += 1
                        sumx[lab, j] += val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]
                # not nan
                if val == val:
                    nobs[lab, 0] += 1
                    sumx[lab, 0] += val

        for i in range(ncounts):
            for j in range(K):
                count = nobs[i, j]
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = sumx[i, j] / count


@cython.wraparound(False)
@cython.boundscheck(False)
def group_ohlc_float32(ndarray[float32_t, ndim=2] out,
                  ndarray[int64_t] counts,
                  ndarray[float32_t, ndim=2] values,
                  ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab
        float32_t val, count
        Py_ssize_t ngroups = len(counts)

    if len(labels) == 0:
        return

    N, K = (<object> values).shape

    if out.shape[1] != 4:
        raise ValueError('Output array must have 4 columns')

    if K > 1:
        raise NotImplementedError("Argument 'values' must have only "
                                  "one dimension")
    out.fill(np.nan)

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab == -1:
                continue

            counts[lab] += 1
            val = values[i, 0]
            if val != val:
                continue

            if out[lab, 0] != out[lab, 0]:
                out[lab, 0] = out[lab, 1] = out[lab, 2] = out[lab, 3] = val
            else:
                out[lab, 1] = max(out[lab, 1], val)
                out[lab, 2] = min(out[lab, 2], val)
                out[lab, 3] = val

#----------------------------------------------------------------------
# group_nth, group_last
#----------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
def group_last_float64(ndarray[float64_t, ndim=2] out,
                        ndarray[int64_t] counts,
                        ndarray[float64_t, ndim=2] values,
                        ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float64_t val, count
        ndarray[float64_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val and val != NAN:
                    nobs[lab, j] += 1
                    resx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = resx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_nth_float64(ndarray[float64_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float64_t, ndim=2] values,
                       ndarray[int64_t] labels, int64_t rank):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float64_t val, count
        ndarray[float64_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val and val != NAN:
                    nobs[lab, j] += 1
                    if nobs[lab, j] == rank:
                        resx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = resx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_last_float32(ndarray[float32_t, ndim=2] out,
                        ndarray[int64_t] counts,
                        ndarray[float32_t, ndim=2] values,
                        ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float32_t val, count
        ndarray[float32_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val and val != NAN:
                    nobs[lab, j] += 1
                    resx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = resx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_nth_float32(ndarray[float32_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float32_t, ndim=2] values,
                       ndarray[int64_t] labels, int64_t rank):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float32_t val, count
        ndarray[float32_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val and val != NAN:
                    nobs[lab, j] += 1
                    if nobs[lab, j] == rank:
                        resx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = resx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_last_int64(ndarray[int64_t, ndim=2] out,
                        ndarray[int64_t] counts,
                        ndarray[int64_t, ndim=2] values,
                        ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        int64_t val, count
        ndarray[int64_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val and val != iNaT:
                    nobs[lab, j] += 1
                    resx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = iNaT
                else:
                    out[i, j] = resx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_nth_int64(ndarray[int64_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[int64_t, ndim=2] values,
                       ndarray[int64_t] labels, int64_t rank):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        int64_t val, count
        ndarray[int64_t, ndim=2] resx
        ndarray[int64_t, ndim=2] nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros((<object> out).shape, dtype=np.int64)
    resx = np.empty_like(out)

    N, K = (<object> values).shape

    with nogil:
        for i in range(N):
            lab = labels[i]
            if lab < 0:
                continue

            counts[lab] += 1
            for j in range(K):
                val = values[i, j]

                # not nan
                if val == val and val != iNaT:
                    nobs[lab, j] += 1
                    if nobs[lab, j] == rank:
                        resx[lab, j] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = iNaT
                else:
                    out[i, j] = resx[i, j]

#----------------------------------------------------------------------
# group_min, group_max
#----------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_float64(ndarray[float64_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float64_t, ndim=2] values,
                       ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float64_t val, count
        ndarray[float64_t, ndim=2] maxx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    maxx = np.empty_like(out)
    maxx.fill(-np.inf)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val and val != NAN:
                        nobs[lab, j] += 1
                        if val > maxx[lab, j]:
                            maxx[lab, j] = val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val and val != NAN:
                    nobs[lab, 0] += 1
                    if val > maxx[lab, 0]:
                        maxx[lab, 0] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = maxx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_float64(ndarray[float64_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float64_t, ndim=2] values,
                       ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float64_t val, count
        ndarray[float64_t, ndim=2] minx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(np.inf)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val and val != NAN:

                        nobs[lab, j] += 1
                        if val < minx[lab, j]:
                            minx[lab, j] = val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val and val != NAN:
                    nobs[lab, 0] += 1
                    if val < minx[lab, 0]:
                        minx[lab, 0] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = minx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_float32(ndarray[float32_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float32_t, ndim=2] values,
                       ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float32_t val, count
        ndarray[float32_t, ndim=2] maxx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    maxx = np.empty_like(out)
    maxx.fill(-np.inf)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val and val != NAN:
                        nobs[lab, j] += 1
                        if val > maxx[lab, j]:
                            maxx[lab, j] = val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val and val != NAN:
                    nobs[lab, 0] += 1
                    if val > maxx[lab, 0]:
                        maxx[lab, 0] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = maxx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_float32(ndarray[float32_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[float32_t, ndim=2] values,
                       ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        float32_t val, count
        ndarray[float32_t, ndim=2] minx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(np.inf)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val and val != NAN:

                        nobs[lab, j] += 1
                        if val < minx[lab, j]:
                            minx[lab, j] = val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val and val != NAN:
                    nobs[lab, 0] += 1
                    if val < minx[lab, 0]:
                        minx[lab, 0] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = NAN
                else:
                    out[i, j] = minx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_max_int64(ndarray[int64_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[int64_t, ndim=2] values,
                       ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        int64_t val, count
        ndarray[int64_t, ndim=2] maxx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    maxx = np.empty_like(out)
    maxx.fill(-_int64_max)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val and val != iNaT:
                        nobs[lab, j] += 1
                        if val > maxx[lab, j]:
                            maxx[lab, j] = val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val and val != iNaT:
                    nobs[lab, 0] += 1
                    if val > maxx[lab, 0]:
                        maxx[lab, 0] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = iNaT
                else:
                    out[i, j] = maxx[i, j]


@cython.wraparound(False)
@cython.boundscheck(False)
def group_min_int64(ndarray[int64_t, ndim=2] out,
                       ndarray[int64_t] counts,
                       ndarray[int64_t, ndim=2] values,
                       ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, lab, ncounts = len(counts)
        int64_t val, count
        ndarray[int64_t, ndim=2] minx, nobs

    if not len(values) == len(labels):
        raise AssertionError("len(index) != len(labels)")

    nobs = np.zeros_like(out)

    minx = np.empty_like(out)
    minx.fill(_int64_max)

    N, K = (<object> values).shape

    with nogil:
        if K > 1:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                for j in range(K):
                    val = values[i, j]

                    # not nan
                    if val == val and val != iNaT:

                        nobs[lab, j] += 1
                        if val < minx[lab, j]:
                            minx[lab, j] = val
        else:
            for i in range(N):
                lab = labels[i]
                if lab < 0:
                    continue

                counts[lab] += 1
                val = values[i, 0]

                # not nan
                if val == val and val != iNaT:
                    nobs[lab, 0] += 1
                    if val < minx[lab, 0]:
                        minx[lab, 0] = val

        for i in range(ncounts):
            for j in range(K):
                if nobs[i, j] == 0:
                    out[i, j] = iNaT
                else:
                    out[i, j] = minx[i, j]

#----------------------------------------------------------------------
# other grouping functions not needing a template
#----------------------------------------------------------------------


def group_median_float64(ndarray[float64_t, ndim=2] out,
                         ndarray[int64_t] counts,
                         ndarray[float64_t, ndim=2] values,
                         ndarray[int64_t] labels):
    """
    Only aggregates on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, ngroups, size
        ndarray[int64_t] _counts
        ndarray data
        float64_t* ptr
    ngroups = len(counts)
    N, K = (<object> values).shape

    indexer, _counts = groupsort_indexer(labels, ngroups)
    counts[:] = _counts[1:]

    data = np.empty((K, N), dtype=np.float64)
    ptr = <float64_t*> data.data

    take_2d_axis1_float64_float64(values.T, indexer, out=data)

    for i in range(K):
        # exclude NA group
        ptr += _counts[0]
        for j in range(ngroups):
            size = _counts[j + 1]
            out[j, i] = _median_linear(ptr, size)
            ptr += size


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cumprod_float64(float64_t[:, :] out,
                          float64_t[:, :] values,
                          int64_t[:] labels,
                          float64_t[:, :] accum):
    """
    Only transforms on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, size
        float64_t val
        int64_t lab

    N, K = (<object> values).shape
    accum = np.ones_like(accum)

    with nogil:
        for i in range(N):
            lab = labels[i]

            if lab < 0:
                continue
            for j in range(K):
                val = values[i, j]
                if val == val:
                    accum[lab, j] *= val
                    out[i, j] = accum[lab, j]


@cython.boundscheck(False)
@cython.wraparound(False)
def group_cumsum(numeric[:, :] out,
                 numeric[:, :] values,
                 int64_t[:] labels,
                 numeric[:, :] accum):
    """
    Only transforms on axis=0
    """
    cdef:
        Py_ssize_t i, j, N, K, size
        numeric val
        int64_t lab

    N, K = (<object> values).shape
    accum = np.zeros_like(accum)

    with nogil:
        for i in range(N):
            lab = labels[i]

            if lab < 0:
                continue
            for j in range(K):
                val = values[i, j]
                if val == val:
                    accum[lab, j] += val
                    out[i, j] = accum[lab, j]


@cython.boundscheck(False)
@cython.wraparound(False)
def group_shift_indexer(int64_t[:] out, int64_t[:] labels,
                        int ngroups, int periods):
    cdef:
        Py_ssize_t N, i, j, ii
        int offset, sign
        int64_t lab, idxer, idxer_slot
        int64_t[:] label_seen = np.zeros(ngroups, dtype=np.int64)
        int64_t[:, :] label_indexer

    N, = (<object> labels).shape

    if periods < 0:
        periods = -periods
        offset = N - 1
        sign = -1
    elif periods > 0:
        offset = 0
        sign = 1

    if periods == 0:
        with nogil:
            for i in range(N):
                out[i] = i
    else:
        # array of each previous indexer seen
        label_indexer = np.zeros((ngroups, periods), dtype=np.int64)
        with nogil:
            for i in range(N):
                ## reverse iterator if shifting backwards
                ii = offset + sign * i
                lab = labels[ii]

                # Skip null keys
                if lab == -1:
                    out[ii] = -1
                    continue

                label_seen[lab] += 1

                idxer_slot = label_seen[lab] % periods
                idxer = label_indexer[lab, idxer_slot]

                if label_seen[lab] > periods:
                    out[ii] = idxer
                else:
                    out[ii] = -1

                label_indexer[lab, idxer_slot] = ii
