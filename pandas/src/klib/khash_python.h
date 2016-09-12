#include <Python.h>

#include "khash.h"

// Previously we were using the built in cpython hash function for doubles
// python 2.7 https://github.com/python/cpython/blob/2.7/Objects/object.c#L1021
// python 3.5 https://github.com/python/cpython/blob/3.5/Python/pyhash.c#L85

// The python 3 hash function has the invariant hash(x) == hash(int(x)) == hash(decimal(x))
// and the size of hash may be different by platform / version (long in py2, Py_ssize_t in py3).
// We don't need those invariants because types will be cast before hashing, and if Py_ssize_t
// is 64 bits the truncation causes collission issues.  Given all that, we use our own
// simple hash, viewing the double bytes as an int64 and using khash's default
// hash for 64 bit integers.
// GH 13436
khint64_t PANDAS_INLINE asint64(double key) {
  return *(khint64_t *)(&key);
}
#define kh_float64_hash_func(key) (khint32_t)((asint64(key))>>33^(asint64(key))^(asint64(key))<<11)
#define kh_float64_hash_equal(a, b) ((a) == (b) || ((b) != (b) && (a) != (a)))

#define KHASH_MAP_INIT_FLOAT64(name, khval_t)								\
	KHASH_INIT(name, khfloat64_t, khval_t, 1, kh_float64_hash_func, kh_float64_hash_equal)

KHASH_MAP_INIT_FLOAT64(float64, size_t)


int PANDAS_INLINE pyobject_cmp(PyObject* a, PyObject* b) {
	int result = PyObject_RichCompareBool(a, b, Py_EQ);
	if (result < 0) {
		PyErr_Clear();
		return 0;
	}
	return result;
}


#define kh_python_hash_func(key) (PyObject_Hash(key))
#define kh_python_hash_equal(a, b) (pyobject_cmp(a, b))


// Python object

typedef PyObject* kh_pyobject_t;

#define KHASH_MAP_INIT_PYOBJECT(name, khval_t)							\
	KHASH_INIT(name, kh_pyobject_t, khval_t, 1,						\
			   kh_python_hash_func, kh_python_hash_equal)

KHASH_MAP_INIT_PYOBJECT(pymap, Py_ssize_t)

#define KHASH_SET_INIT_PYOBJECT(name)                                  \
	KHASH_INIT(name, kh_pyobject_t, char, 0,     \
			   kh_python_hash_func, kh_python_hash_equal)

KHASH_SET_INIT_PYOBJECT(pyset)

#define kh_exist_pymap(h, k) (kh_exist(h, k))
#define kh_exist_pyset(h, k) (kh_exist(h, k))

KHASH_MAP_INIT_STR(strbox, kh_pyobject_t)
