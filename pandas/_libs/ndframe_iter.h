// TODO: could use better naming than "PdBlockIter"
// showing that this is really "deconstructed" block data
// in native form

#ifndef PANDAS__NDFRAME_ITER
#define PANDAS__NDFRAME_ITER

#define PY_SSIZE_T_CLEAN
#include <Python.h>


// Struct containing a pointer to len # of NDArrays
// The order of each item in data should match the
// order specified by the BlockManager
typedef struct {
  Py_ssize_t len;
  PyObject **ndarrays;  // TODO: need some kind of destructor
} PdBlocksIter;

// Provided a DataFrame and axis deconstructs the block
// data to match the order represented by the BlockManager
// Returns NULL on error
PdBlocksIter *PdFrameIter_New(PyObject *df, int axis);

#endif
