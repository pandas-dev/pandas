#include "ndframe_iter.h"

// returns -1 on error
// performs no bounds checking, so will likely segfault if you pass
// and inappropriate dim for the object
static Py_ssize_t getDimLength(PyObject *df, Py_ssize_t dim) {
  PyObject *shape;
  Py_ssize_t ncols;
  
  shape = PyObject_GetAttrString(df, "shape");
  if (shape == NULL) {
    return -1;
  }
  
  if (!PyTuple_Check(shape)) {
    Py_DECREF(shape);
    return -1;
  }

  ncols = PyLong_AsLongLong(PyTuple_GET_ITEM(shape, dim));
  Py_DECREF(shape);

  return ncols;
}

// Return the blocks object associated with a dataframe
// Checks that the return value is a Tuple
static PyObject *getBlocksTuple(PyObject *df) {
  PyObject *blockManager, *blocks;

  blockManager = PyObject_GetAttrString(df, "_data");
  if (blockManager == NULL) {
    return NULL;
  }

  blocks = PyObject_GetAttrString(blockManager, "blocks");
  Py_DECREF(blockManager);
  if (blocks == NULL) {
    return NULL;
  }

  if (!PyTuple_Check(blocks)) {
    // TODO: Set error message here
    Py_DECREF(blocks);
    return NULL;
  }

  return blocks;
}

// Return the integer placements of the blocks columns in its owning frame
// Checks that the return value is a List
static PyObject *getManagerLocationsAsList(PyObject *block) {
  PyObject *managerLocs, *ndarray, *list;

  managerLocs = PyObject_GetAttrString(block, "mgr_locs");
  if (managerLocs == NULL) {
    return NULL;
  }

  // TODO: we could probably just supply managerLocs to the list()
  // built-in instead of going from mgr_locs->as_array->tolist
  ndarray = PyObject_GetAttrString(managerLocs, "as_array");
  Py_DECREF(managerLocs);
  if (ndarray == NULL) {
    return NULL;
  }

  list = PyObject_CallMethod(ndarray, "tolist", NULL);
  Py_DECREF(ndarray);
  if (list == NULL) {
    return NULL;
  } else if (!PyList_Check(list)) {
    Py_DECREF(list);
    return NULL;
  }
  
  return list;
}


PdBlocksIter *PdFrameIter_New(PyObject *df, int axis) {
  PyObject *blocks, *block, *managerLocs, *blockValues;
  char ***data;  // individual ndarrays; C-order should match column order
  Py_ssize_t i, j, loc, ncols;
  NpyIter *iter;
  NpyIter_IterNextFunc *iternext;

  printf("we are in!");
  ncols = getDimLength(df, 1);
  blocks = getBlocksTuple(df);
  if (blocks == NULL) {
    return NULL;
  }

  printf("down here\n"); 
  data = PyObject_Malloc(sizeof(char **) * ncols);
  if (data == NULL) {
    Py_DECREF(blocks);
    return NULL;
  }

  for (i = 0; i < PyTuple_GET_SIZE(blocks); i++) {
    block = PyTuple_GET_ITEM(blocks, i);
    
    blockValues = PyObject_CallMethod(block, "get_block_values", NULL);
    if (blockValues == NULL) {
      Py_DECREF(blocks);
      PyObject_Free(data);
      return NULL;
    }

    managerLocs = getManagerLocationsAsList(block);
    if (managerLocs == NULL) {
      Py_DECREF(blockValues);
      Py_DECREF(blocks);
      PyObject_Free(data);
      return NULL;
    }

    // TODO: we have the array data and we know the indices
    // of where they are located in the dataframe, so figure out how we
    // should iterate and store that knowledge
    Py_DECREF(managerLocs);
    Py_DECREF(blockValues);
  }

  Py_DECREF(blocks);

  return NULL;
}
