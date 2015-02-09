#ifndef _KLIB_KHASH_PYTHON_H_
#define _KLIB_KHASH_PYTHON_H_

#include <Python.h>

#ifndef PANDAS_INLINE
  #if defined(__GNUC__)
    #define PANDAS_INLINE __inline__
  #elif defined(_MSC_VER)
    #define PANDAS_INLINE __inline
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define PANDAS_INLINE inline
  #else
    #define PANDAS_INLINE
  #endif
#endif

#define kh_inline PANDAS_INLINE
#include "khash.h"

#define kh_exist_str(h, k) (kh_exist(h, k))
#define kh_exist_float64(h, k) (kh_exist(h, k))
#define kh_exist_int64(h, k) (kh_exist(h, k))
#define kh_exist_int32(h, k) (kh_exist(h, k))

#include "xxhash/xxhash.h"

/*
 * By default khash uses crappy x31 hash function which puts strings that
 * differ only in the last character into neighbouring buckets which is not
 * good given that quadratic probing tries small steps first.
 *
 * xxhash gives better bucket distribution and performance-wise is great for
 * long-ish strings, but it is a bit slower than x31 on the shortest ones
 * (turns out at length == 2 the difference is already negligible).
 *
 * Inlining will hinder merging in upstream releases, but 1-character strings
 * are a valid use case for pandas, so let's pre-calculate a vector of 256
 * values to avoid calling two functions (strlen and XXH32) if there's only one
 * character to hash.
 *
 * This table was generated with the following code.  Feel free to re-run it if
 * an update comes in:

#include <stdio.h>
#include "xxhash.h"

int main(int argc, char *argv[])
{
  printf("static khint_t XXH32_EMPTY_HASH = 0x%08x;\n",
         XXH32("", 0, 0xdeadbeef));
  printf("static khint_t XXH32_ONECHAR_HASH[256] = {");
  unsigned char s[2] = {0};
  for (int i = 0; i < 256; ++i) {
    if (i % 8 == 0) {
      printf("\n    ");
    }
    s[0] = i;
    printf("0x%08x", XXH32(s, 1, 0xdeadbeef));
    if (i < 255) {
      printf(", ");
    }
  }
  printf("\n};\n");
  return 0;
}
*/

static khint_t XXH32_EMPTY_HASH = 0xc372c6cb;
static khint_t XXH32_ONECHAR_HASH[256] = {
    0x39110451, 0xd3efa134, 0xea8d6dc4, 0xe59a066b, 0x89f3a4f5, 0xdcce5bc9, 0x44be0c3e, 0x96469248, 
    0x7885ddeb, 0x24417b24, 0xb77b30b2, 0xa83d21eb, 0x6f6ba52b, 0x7315bbe5, 0xce858701, 0x52299f26, 
    0x440ec810, 0xd02a934f, 0xf873d394, 0xd168a8e1, 0x31c30198, 0x37c3967b, 0xc1bdbdf8, 0x3ddaf3cc, 
    0xb7222f4a, 0x96625cdf, 0xabf92a2f, 0x69e97975, 0x55f24523, 0x6b1abaa0, 0xe5b033ab, 0x9e21842c, 
    0x3ac2a339, 0x827b0af2, 0xd7ea0f97, 0x72317ee6, 0xe6bd4439, 0xb0b183f1, 0xca90e5e0, 0x57960753, 
    0x6eefe374, 0xb9c9c5b5, 0x57396d1f, 0x6db79351, 0xab55c12d, 0x32229df4, 0xbfa3a164, 0x58f9f4ba, 
    0x5987c643, 0xffbfa961, 0x1080d4eb, 0xc5c3d846, 0x16a7fd8e, 0xed29fd3a, 0x8d78613d, 0xd088b720, 
    0x8d597f4c, 0x2df1ce8f, 0x79bc5215, 0x749d67c1, 0xa9ad300c, 0x60c6237d, 0xeeb080e7, 0xb74eef62, 
    0x6ddba2f2, 0x3d9f18cf, 0x0b6ad1bd, 0xc7a33d19, 0x3cb6352f, 0x872839f9, 0x259ced1e, 0x0f9d713b, 
    0x6816620f, 0x8d2c96a7, 0x377fb2f9, 0x2616b5b5, 0x9bae3a05, 0x8368a004, 0x3a67fd94, 0x312529c4, 
    0xc9238f87, 0x3e85e142, 0x973dedc6, 0xcbc3d4ba, 0xd2629b58, 0x2aae9a6d, 0x82ffc598, 0x4a8512b3, 
    0x51146ceb, 0x85ddc3f4, 0xa83b942f, 0x55769a32, 0xf7fa3fdf, 0xfbe35842, 0x342ff574, 0x848400a6, 
    0x92707153, 0x48cd58fd, 0xbdae4a11, 0x701bbadb, 0x4a5b37c4, 0x98770eeb, 0xfc1b98fc, 0x05dd6894, 
    0xd3ba005c, 0x453bc774, 0xfe186d14, 0xa25acde2, 0xcc738313, 0x1dbdefa7, 0x83ed6f1e, 0xf9d8e195, 
    0x5f10c546, 0xf22c5a0f, 0x31da5f5e, 0x5341c163, 0xabd3f750, 0x882e33d8, 0x4d8105cd, 0xc1f6f3d9, 
    0x347e1d5c, 0xdb06193c, 0x64841a53, 0x3991a6e6, 0x0abdd625, 0xedcf00f7, 0xa8e64229, 0x2fc9029b, 
    0x4fc5ca41, 0x1f5aaae5, 0x29bdda91, 0x55446dae, 0x1566ec40, 0x9ac8391e, 0xcd4d6ab1, 0x0f3807f6, 
    0xf3be6887, 0x9f4b88bd, 0x33c401df, 0xaa9df64f, 0xce5c70ac, 0x9ee55a87, 0x4cb91c84, 0x8c322b3d, 
    0x8e40fb24, 0x3af430fb, 0xeea567c2, 0xe80c7dc2, 0x6f619449, 0xe0ca8048, 0x984c626e, 0x50bf1281, 
    0x4895cbee, 0x5d016a96, 0xe58b8980, 0x3457ef7c, 0x2a24f819, 0x0641cc30, 0xbddc5f84, 0x03ce4656, 
    0xbcb73c9c, 0xcd29be82, 0x0930d945, 0xf3fc8e3c, 0xbed775cd, 0xd6668fae, 0x6876f949, 0xcf34fbd7, 
    0x0537d916, 0x7efd5f26, 0xb2d32520, 0x10d58995, 0x19d64e1c, 0xacae767c, 0xf23a4e7d, 0xdcb654fe, 
    0xe1ec9a9f, 0x3061302b, 0x453a0b7c, 0xe845436e, 0xb2b690df, 0x245c17b5, 0x756a9374, 0x470998f5, 
    0xe31a5f5b, 0x60dbad02, 0xf738299d, 0x0db8b11a, 0xd34cb801, 0xb2f3597d, 0xa627e466, 0xda4f9935, 
    0x5c58e1df, 0x4b5319d6, 0x48acc08f, 0xce18d68e, 0xeb995e7f, 0x11a07cba, 0x025127b2, 0xd1325331, 
    0x55d76240, 0x281bba14, 0xb9ac069d, 0x25e60bcc, 0xf077fbd3, 0xe460ece9, 0x725a9971, 0xa6b5c6b4, 
    0xe5f216a3, 0xbee80d71, 0x1a049114, 0x851012d4, 0xa6e175cc, 0x6ec98c95, 0x56a77202, 0x7e2ab05f, 
    0x4850279c, 0x1b009afe, 0xf71e36b6, 0x9cadc37a, 0x43a167da, 0x5d75b5f3, 0xc432215c, 0x93ff1905, 
    0x8764d057, 0xf44cd35d, 0x03d3a324, 0xd65a5047, 0xe872b4d8, 0x8dcb9a23, 0xfebf9113, 0x59701be9, 
    0xdf9f6090, 0xce9b2907, 0x664c6a5a, 0x81bfefc4, 0x13829979, 0xda98b6ab, 0x7b7e9ff0, 0x13c24005, 
    0xcee61b6b, 0x15737a85, 0xe2f95e48, 0xf2136570, 0xd1ccfdab, 0xa9adfb16, 0x1f7339a9, 0x83247f43, 
    0x68c6c8bf, 0x5046f6fc, 0x2d3dea84, 0x79a0be74, 0x39dd7eb3, 0x4d5cc636, 0xe4e1352d, 0xd1317a99
};

/* Seed value is chosen arbitrarily. */
static khint_t XXH32_SEED = 0xdeadbeef;

static khint_t PANDAS_INLINE str_xxhash_hash_func(kh_cstr_t key) {
    if (!key[0]) {
        return XXH32_EMPTY_HASH;
    }
    if (!key[1]) {
        return XXH32_ONECHAR_HASH[(uint8_t)key[0]];
    }
    return XXH32(key, strlen(key), XXH32_SEED);
}

KHASH_INIT(str, kh_cstr_t, size_t, 1,
           str_xxhash_hash_func, kh_str_hash_equal)

KHASH_MAP_INIT_INT(int32, size_t)
KHASH_MAP_INIT_INT64(int64, size_t)

// kludge

#define kh_float64_hash_func _Py_HashDouble
#define kh_float64_hash_equal kh_int64_hash_equal

#define KHASH_MAP_INIT_FLOAT64(name, khval_t)								\
	KHASH_INIT(name, khfloat64_t, khval_t, 1, kh_float64_hash_func, kh_float64_hash_equal)

KHASH_MAP_INIT_FLOAT64(float64, size_t)


static int PANDAS_INLINE pyobject_cmp(PyObject* a, PyObject* b) {
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

KHASH_INIT(strbox, kh_cstr_t, kh_pyobject_t, 1,
           str_xxhash_hash_func, kh_str_hash_equal)

/* Plain old C buffer structure */
typedef struct {
    kh_cstr_t buf;
    Py_ssize_t len;
} cbuf_t;

static khint_t PANDAS_INLINE cbuf_xxhash(cbuf_t val) {
    switch (val.len) {
    case 0:
        return XXH32_EMPTY_HASH;
    case 1:
        return XXH32_ONECHAR_HASH[(uint8_t)val.buf[0]];
    default:
        return XXH32(val.buf, val.len, XXH32_SEED);
    }
}

static int PANDAS_INLINE cbuf_equal(cbuf_t a, cbuf_t b) {
    int i;
    if (a.len != b.len) {
        return 0;
    }
    if (a.buf == b.buf) {
        return 1;
    }
    for (i = 0; i < a.len; ++i) {
        if (a.buf[i] != b.buf[i]) {
            return 0;
        }
    }
    return 1;
}

/* [cbuf_t -> size_t] hash map */
KHASH_INIT(cbuf_map, cbuf_t, size_t, 1, cbuf_xxhash, cbuf_equal)
#define kh_exist_cbuf_map(h, k) (kh_exist(h, k))

#endif /* _KLIB_KHASH_PYTHON_H_ */
