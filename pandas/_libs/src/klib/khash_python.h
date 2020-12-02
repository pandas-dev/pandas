#include <string.h>
#include <Python.h>

// khash should report usage to tracemalloc
#if PY_VERSION_HEX >= 0x03060000
#include <pymem.h>
#if PY_VERSION_HEX < 0x03070000
#define PyTraceMalloc_Track _PyTraceMalloc_Track
#define PyTraceMalloc_Untrack _PyTraceMalloc_Untrack
#endif
#else
#define PyTraceMalloc_Track(...)
#define PyTraceMalloc_Untrack(...)
#endif


static const int KHASH_TRACE_DOMAIN = 424242;
void *traced_malloc(size_t size){
    void * ptr = malloc(size);
    if(ptr!=NULL){
        PyTraceMalloc_Track(KHASH_TRACE_DOMAIN, (uintptr_t)ptr, size);
    }
    return ptr;
}

void *traced_calloc(size_t num, size_t size){
    void * ptr = calloc(num, size);
    if(ptr!=NULL){
        PyTraceMalloc_Track(KHASH_TRACE_DOMAIN, (uintptr_t)ptr, num*size);
    }
    return ptr;
}

void *traced_realloc(void* old_ptr, size_t size){
    void * ptr = realloc(old_ptr, size);
    if(ptr!=NULL){
        if(old_ptr != ptr){
            PyTraceMalloc_Untrack(KHASH_TRACE_DOMAIN, (uintptr_t)old_ptr);
        }
        PyTraceMalloc_Track(KHASH_TRACE_DOMAIN, (uintptr_t)ptr, size);
    }
    return ptr;
}

void traced_free(void* ptr){
    if(ptr!=NULL){
        PyTraceMalloc_Untrack(KHASH_TRACE_DOMAIN, (uintptr_t)ptr);
    }
    free(ptr);
}


#define KHASH_MALLOC traced_malloc
#define KHASH_REALLOC traced_realloc
#define KHASH_CALLOC traced_calloc
#define KHASH_FREE traced_free
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
// GH 13436 showed that _Py_HashDouble doesn't work well with khash
// GH 28303 showed, that the simple xoring-version isn't good enough
// See GH 36729 for evaluation of the currently used murmur2-hash version
// An interesting alternative to expensive murmur2-hash would be to change
// the probing strategy and use e.g. the probing strategy from CPython's
// implementation of dicts, which shines for smaller sizes but is more
// predisposed to superlinear running times (see GH 36729 for comparison)


khint64_t PANDAS_INLINE asint64(double key) {
    khint64_t val;
    memcpy(&val, &key, sizeof(double));
    return val;
}

khint32_t PANDAS_INLINE asint32(float key) {
    khint32_t val;
    memcpy(&val, &key, sizeof(float));
    return val;
}

#define ZERO_HASH 0
#define NAN_HASH  0

khint32_t PANDAS_INLINE kh_float64_hash_func(double val){
    // 0.0 and -0.0 should have the same hash:
    if (val == 0.0){
        return ZERO_HASH;
    }
    // all nans should have the same hash:
    if ( val!=val ){
        return NAN_HASH;
    }
    khint64_t as_int = asint64(val);
    return murmur2_64to32(as_int);
}

khint32_t PANDAS_INLINE kh_float32_hash_func(float val){
    // 0.0 and -0.0 should have the same hash:
    if (val == 0.0f){
        return ZERO_HASH;
    }
    // all nans should have the same hash:
    if ( val!=val ){
        return NAN_HASH;
    }
    khint32_t as_int = asint32(val);
    return murmur2_32to32(as_int);
}

#define kh_floats_hash_equal(a, b) ((a) == (b) || ((b) != (b) && (a) != (a)))

#define KHASH_MAP_INIT_FLOAT64(name, khval_t)								\
	KHASH_INIT(name, khfloat64_t, khval_t, 1, kh_float64_hash_func, kh_floats_hash_equal)

KHASH_MAP_INIT_FLOAT64(float64, size_t)

#define KHASH_MAP_INIT_FLOAT32(name, khval_t)								\
	KHASH_INIT(name, khfloat32_t, khval_t, 1, kh_float32_hash_func, kh_floats_hash_equal)

KHASH_MAP_INIT_FLOAT32(float32, size_t)


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
	kh_str_starts_t *result = (kh_str_starts_t*)KHASH_CALLOC(1, sizeof(kh_str_starts_t));
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

khint_t PANDAS_INLINE kh_get_str_starts_item(const kh_str_starts_t* table, const char* key) {
    unsigned char ch = *key;
	if (table->starts[ch]) {
		if (ch == '\0' || kh_get_str(table->table, key) != table->table->n_buckets) return 1;
	}
    return 0;
}

void PANDAS_INLINE kh_destroy_str_starts(kh_str_starts_t* table) {
	kh_destroy_str(table->table);
	KHASH_FREE(table);
}

void PANDAS_INLINE kh_resize_str_starts(kh_str_starts_t* table, khint_t val) {
	kh_resize_str(table->table, val);
}
