#ifndef NUMBA_ARYSTRUCT_H_
#define NUMBA_ARYSTRUCT_H_
/*
 * Fill in the *arystruct* with information from the Numpy array *obj*.
 * *arystruct*'s layout is defined in numba.targets.arrayobj (look
 * for the ArrayTemplate class).
 */

typedef struct {
    void     *meminfo;  /* see _nrt_python.c and nrt.h in numba/core/runtime */
    PyObject *parent;
    npy_intp nitems;
    npy_intp itemsize;
    void *data;

    npy_intp shape_and_strides[];
} arystruct_t;


#endif  /* NUMBA_ARYSTRUCT_H_ */

