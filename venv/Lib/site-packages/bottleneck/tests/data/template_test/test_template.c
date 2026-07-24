// Copyright 2010-2019 Keith Goodman
// Copyright 2019 Bottleneck Developers
#include "bottleneck.h"
#include "iterators.h"


/* nanmin, nanmax -------------------------------------------------------- */

/* repeat = {'NAME':      ['nanmin',         'nanmax'],
             'COMPARE':   ['<=',             '>='],
             'BIG_FLOAT': ['BN_INFINITY',    '-BN_INFINITY'],
             'BIG_INT':   ['NPY_MAX_DTYPE0', 'NPY_MIN_DTYPE0']} */
/* dtype = [['float64'], ['float32']] */
FOO(NAME, DTYPE0) {
    npy_DTYPE0 bar = BIG_FLOAT;
    if (bar COMPARE 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}
/* dtype end */

/* dtype = [['int64'], ['int32']] */
FOO(NAME, DTYPE0) {
    npy_DTYPE0 bar = BIG_FLOAT;
    if (bar COMPARE 0) {
        bar = 0;
    }
    return PyFloat_FromDouble(bar);
}
/* dtype end */

REDUCE_MAIN(NAME, 0)
/* repeat end */
