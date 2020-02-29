#include "assert.h"
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


// Given a DataFrame, this will create a struct containing both
// the length of ndarrays the frame uses along with an array of
// PyArrayObjects, appearing in C-order along the requested axis.
//
// For example, if we had a DataFrame that looks as follows
//    a     b   c
// 0  1   2.0   3
// 1  4   5.0   6
//
// This would deconstruct the blocks and provide a struct
// with len 3; accessing
// ndarrays[0] will yield ndarray([1, 4])
// ndarrays[1] will yield ndarray([2., 5.])
// ndarrays[2] will yield ndarray([3, 6])
//
// The goal of this is to provide a performant way in C to
// maintain axis order and dtype information of the
// supplied DataFrame.
PdOrderedArrays *PdOrderedArrays_New(PyObject *df, int axis) {
  PyObject *blocks, *block, *managerLocs, *blockValues, *ndarr, *key;
  PyObject **ndarrays;
  Py_ssize_t i, j, loc, ncols;
  PdOrderedArrays *result;

  // 
  assert(axis == 0);

  ncols = getDimLength(df, 1);
  blocks = getBlocksTuple(df);
  if (blocks == NULL) {
    return NULL;
  }

  ndarrays = PyObject_Malloc(sizeof(PyObject *) * ncols);
  if (ndarrays == NULL) {
    Py_DECREF(blocks);
    return NULL;
  }


  // Iterate over all of the blocks, getting a slice from
  // each block sequentially and storing it in the appropriate order
  // in `arrays`
  for (i = 0; i < PyTuple_GET_SIZE(blocks); i++) {
    block = PyTuple_GET_ITEM(blocks, i);
    
    blockValues = PyObject_CallMethod(block, "get_block_values", NULL);
    if (blockValues == NULL) {
      Py_DECREF(blocks);
      PyObject_Free(ndarrays);
      return NULL;
    }

    managerLocs = getManagerLocationsAsList(block);
    if (managerLocs == NULL) {
      Py_DECREF(blockValues);
      Py_DECREF(blocks);
      PyObject_Free(ndarrays);
      return NULL;
    }

    // managerLocs tells use where each column of the block
    // exists in the maintaining dataframe, so iterate
    // managerLocs sequentially and assign, get the column of
    // data and assign that reference to the appropriate position
    // in `ndarrays`
    for (j = 0; j < PyList_GET_SIZE(managerLocs); j++) {
      loc = PyLong_AsLongLong(PyList_GET_ITEM(managerLocs, j));
      key = PyLong_FromLongLong(j);
      if (key == NULL) {
        goto LOOP_ERROR;
      }

      // TODO: We should actually be grabbing an ndarray reference here
      // and storing that in ndarrays; I haven't quite figured out
      // with the Numpy C-API how to do that though
      //
      // For illustration, let's assume we have a DataFrame that looks like:
      //     colA  colB  colC
      //  0     0     X     1
      //  1     2     Y     3
      //
      // If we are currently iterating over the int block, it would look
      // something like this:
      //
      //  [[0, 1],
      //   [2, 3]]
      //
      // and managerLocs would be [0, 2] (note position of int dtypes in dataframe)
      //
      // So in this loop we would ideally create a ndarray that references [0, 2]
      // (i.e. the first column of the above block) and store a pointer to that
      // at arrays[0]. We would then get an array that references [1, 3] and store
      // that at arrays[2].
      //
      // Iteration over the string block should populate an ndarray at arrays[1]
      // that holds a reference to ["X", "Y"]
      ndarr = PyObject_GetItem(blockValues, key);  // TODO: no GetItem; construct ndarray somehow
      Py_DECREF(key);
      if (ndarr == NULL) {
        goto LOOP_ERROR;
      }

      ndarrays[loc] = ndarr;
      continue;
    LOOP_ERROR:
      Py_DECREF(managerLocs);
      Py_DECREF(blockValues);
      Py_DECREF(blocks);
      PyObject_Free(ndarrays);
      return NULL;
    }
    
    Py_DECREF(managerLocs);
    Py_DECREF(blockValues);
  }

  Py_DECREF(blocks);

  result = PyObject_Malloc(sizeof(PdOrderedArrays));
  if (result == NULL) {
    return NULL;
  }

  result->len = ncols;
  result->ndarrays = ndarrays;

  return result;
}

void PdOrderedArrays_Destroy(PdOrderedArrays *orderedArrays) {
  Py_ssize_t i;

  for (i = 0; i < orderedArrays->len; i++) {
    Py_DECREF(orderedArrays->ndarrays++);
  }
  
  PyObject_Free(orderedArrays->ndarrays);
  PyObject_Free(orderedArrays);
}
