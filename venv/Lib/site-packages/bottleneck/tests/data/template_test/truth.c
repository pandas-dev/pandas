#line 1 "{DIRPATH}/test_template.c"
// Copyright 2010-2019 Keith Goodman
// Copyright 2019 Bottleneck Developers
#include "bottleneck.h"
#include "iterators.h"

/* nanmin, nanmax -------------------------------------------------------- */

#line 14
FOO(nanmin, float64) {
    npy_float64 bar = BN_INFINITY;
    if (bar <= 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}

#line 14
FOO(nanmin, float32) {
    npy_float32 bar = BN_INFINITY;
    if (bar <= 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}

#line 24
FOO(nanmin, int64) {
    npy_int64 bar = BN_INFINITY;
    if (bar <= 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}

#line 24
FOO(nanmin, int32) {
    npy_int32 bar = BN_INFINITY;
    if (bar <= 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}

REDUCE_MAIN(nanmin, 0)

#line 14
FOO(nanmax, float64) {
    npy_float64 bar = -BN_INFINITY;
    if (bar >= 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}

#line 14
FOO(nanmax, float32) {
    npy_float32 bar = -BN_INFINITY;
    if (bar >= 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}

#line 24
FOO(nanmax, int64) {
    npy_int64 bar = -BN_INFINITY;
    if (bar >= 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}

#line 24
FOO(nanmax, int32) {
    npy_int32 bar = -BN_INFINITY;
    if (bar >= 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}

REDUCE_MAIN(nanmax, 0)
