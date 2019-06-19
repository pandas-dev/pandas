#include <string.h>
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
  khint64_t val;
  memcpy(&val, &key, sizeof(double));
  return val;
}

// correct for all inputs but not -0.0 and NaNs
#define kh_float64_hash_func_0_NAN(key) (khint32_t)((asint64(key))>>33^(asint64(key))^(asint64(key))<<11)

// correct for all inputs but not NaNs
#define kh_float64_hash_func_NAN(key) ((key) == 0.0 ?                       \
                                        kh_float64_hash_func_0_NAN(0.0) : \
                                        kh_float64_hash_func_0_NAN(key))

// correct for all
#define kh_float64_hash_func(key) ((key) != (key) ?                       \
                                   kh_float64_hash_func_NAN(Py_NAN) :     \
                                   kh_float64_hash_func_NAN(key))

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
    if (result == 0) {  // still could be two NaNs
        return PyFloat_CheckExact(a) &&
               PyFloat_CheckExact(b) &&
               Py_IS_NAN(PyFloat_AS_DOUBLE(a)) &&
               Py_IS_NAN(PyFloat_AS_DOUBLE(b));
    }
	return result;
}

// For PyObject_Hash holds:
//    hash(0.0) == 0 == hash(-0.0)
//    hash(X) == 0 if X is a NaN-value
// so it is OK to use it directly
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

typedef struct {
	kh_str_t *table;
	int starts[256];
} kh_str_starts_t;

typedef kh_str_starts_t* p_kh_str_starts_t;

p_kh_str_starts_t PANDAS_INLINE kh_init_str_starts(void) {
	kh_str_starts_t *result = (kh_str_starts_t*)calloc(1, sizeof(kh_str_starts_t));
	result->table = kh_init_str();
	return result;
}

khint_t PANDAS_INLINE kh_put_str_starts_item(kh_str_starts_t* table, char* key, int* ret) {
    khint_t result = kh_put_str(table->table, key, ret);
	if (*ret != 0) {
		table->starts[(unsigned char)key[0]] = 1;
	}
    return result;
}

khint_t PANDAS_INLINE kh_get_str_starts_item(kh_str_starts_t* table, char* key) {
    unsigned char ch = *key;
	if (table->starts[ch]) {
		if (ch == '\0' || kh_get_str(table->table, key) != table->table->n_buckets) return 1;
	}
    return 0;
}

void PANDAS_INLINE kh_destroy_str_starts(kh_str_starts_t* table) {
	kh_destroy_str(table->table);
	free(table);
}

void PANDAS_INLINE kh_resize_str_starts(kh_str_starts_t* table, khint_t val) {
	kh_resize_str(table->table, val);
}