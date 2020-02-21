// TODO: could use better naming than "PdBlockIter"
// showing that this is really "deconstructed" block data
// in native form

// Struct containing a pointer to len # of NDArrays
// The order of each item in data should match the
// order specified by the BlockManager
typedef struct {
  Py_ssize_t len;
  char ***data;  // each of these points to an array containing numpy data
} PdBlocksIter;

// Provided a DataFrame and axis deconstructs the block
// data to match the order represented by the BlockManager
// Returns NULL on error
PdBlocksIter *PdFrameIter_New(PyObject *df, int axis);

