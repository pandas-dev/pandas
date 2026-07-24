#include "pythoncapi_compat.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#include "mypyc_util.h"
#include "CPy.h"
#include "librt_random.h"

//
// ChaCha8 PRNG with forward secrecy
//

#define CHACHA8_RESEED_INTERVAL 16

typedef struct {
    uint32_t seed[8];       // 256-bit key
    uint32_t buf[16];       // output buffer: one ChaCha8 block
    uint32_t counter;       // block counter
    uint8_t  used;          // index into buf
    uint8_t  n;             // usable values in buf (8 or 16)
    uint8_t  blocks_left;   // blocks until next reseed
} chacha8_rng;

static inline uint32_t
rotl32(uint32_t x, int n) {
    return (x << n) | (x >> (32 - n));
}

#define QUARTERROUND(a, b, c, d) \
    do { \
        a += b; d ^= a; d = rotl32(d, 16); \
        c += d; b ^= c; b = rotl32(b, 12); \
        a += b; d ^= a; d = rotl32(d, 8);  \
        c += d; b ^= c; b = rotl32(b, 7);  \
    } while (0)

static void
chacha8_block(const uint32_t seed[8], uint32_t counter, uint32_t out[16])
{
    // "expand 32-byte k"
    uint32_t s[16] = {
        0x61707865, 0x3320646e, 0x79622d32, 0x6b206574,
        seed[0], seed[1], seed[2], seed[3],
        seed[4], seed[5], seed[6], seed[7],
        counter, 0, 0, 0   // counter (low 32), counter (high 32), nonce
    };

    memcpy(out, s, sizeof(uint32_t) * 16);

    // 4 double-rounds = 8 rounds
    for (int i = 0; i < 4; i++) {
        // Column rounds
        QUARTERROUND(out[0], out[4], out[ 8], out[12]);
        QUARTERROUND(out[1], out[5], out[ 9], out[13]);
        QUARTERROUND(out[2], out[6], out[10], out[14]);
        QUARTERROUND(out[3], out[7], out[11], out[15]);
        // Diagonal rounds
        QUARTERROUND(out[0], out[5], out[10], out[15]);
        QUARTERROUND(out[1], out[6], out[11], out[12]);
        QUARTERROUND(out[2], out[7], out[ 8], out[13]);
        QUARTERROUND(out[3], out[4], out[ 9], out[14]);
    }

    // Add original state back (standard ChaCha finalization)
    for (int i = 0; i < 16; i++)
        out[i] += s[i];
}

// Fill entropy from OS via os.urandom(), which handles short reads,
// EINTR, and platform differences internally.
// Returns 0 on success, -1 on failure (with Python exception set).
static int
fill_os_entropy(void *buf, size_t len)
{
    PyObject *os_mod = PyImport_ImportModule("os");
    if (os_mod == NULL)
        return -1;
    PyObject *bytes = PyObject_CallMethod(os_mod, "urandom", "n", (Py_ssize_t)len);
    Py_DECREF(os_mod);
    if (bytes == NULL)
        return -1;
    memcpy(buf, PyBytes_AS_STRING(bytes), len);
    Py_DECREF(bytes);
    return 0;
}

static void
chacha8_refill(chacha8_rng *rng)
{
    chacha8_block(rng->seed, rng->counter, rng->buf);
    rng->counter++;
    rng->used = 0;
    rng->blocks_left--;

    if (unlikely(rng->blocks_left == 0)) {
        // Forward secrecy reseed: steal last 8 words as new key
        memcpy(rng->seed, rng->buf + 8, sizeof(uint32_t) * 8);
        rng->n = 8;  // only 8 words usable this block
        rng->counter = 0;
        rng->blocks_left = CHACHA8_RESEED_INTERVAL;
    } else {
        rng->n = 16;
    }
}

static inline uint32_t
chacha8_next(chacha8_rng *rng)
{
    if (unlikely(rng->used >= rng->n))
        chacha8_refill(rng);
    return rng->buf[rng->used++];
}

// Return 64 bits of randomness (two consecutive 32-bit words, single bounds check).
static inline uint64_t
chacha8_next64(chacha8_rng *rng)
{
    // Need 2 words available; if fewer than 2, refill first.
    if (unlikely(rng->used + 1 >= rng->n))
        // Use two separate calls to handle block boundary correctly.
        return ((uint64_t)chacha8_next(rng) << 32) | chacha8_next(rng);
    uint32_t hi = rng->buf[rng->used++];
    uint32_t lo = rng->buf[rng->used++];
    return ((uint64_t)hi << 32) | lo;
}

// Return a uniformly distributed random value in [0, range).
// Use Lemire's nearly divisionless method for small ranges, and a portable
// rejection sampler for larger ranges to avoid non-standard 128-bit arithmetic.
static inline uint64_t
chacha8_next_ranged(chacha8_rng *rng, uint64_t range)
{
    assert(range != 0);
    if (likely(range <= UINT32_MAX)) {
        // 32-bit Lemire: multiply r * range to get 64-bit product,
        // upper 32 bits are the result in [0, range).
        uint64_t m = (uint64_t)chacha8_next(rng) * range;
        uint32_t lo = (uint32_t)m;
        if (unlikely(lo < range)) {
            uint32_t thresh = (uint32_t)(-(uint32_t)range) % (uint32_t)range;
            while (lo < thresh) {
                m = (uint64_t)chacha8_next(rng) * range;
                lo = (uint32_t)m;
            }
        }
        return m >> 32;
    }
    // If range is a power of two, masking produces an unbiased result.
    if ((range & (range - 1)) == 0) {
        return chacha8_next64(rng) & (range - 1);
    }
    uint64_t r;
    // In unsigned arithmetic, -range is 2**64 - range, so this computes
    // 2**64 % range. Rejecting values below this threshold leaves exactly
    // floor(2**64 / range) full buckets of size range, avoiding modulo bias.
    uint64_t thresh = -range % range;
    do {
        r = chacha8_next64(rng);
    } while (unlikely(r < thresh));
    return r % range;
}

// Return a random i64 starting at 'start', with 'range' possible values.
// A zero range represents the full 2**64 i64 domain.
static inline int64_t
random_i64_from_range(chacha8_rng *rng, int64_t start, uint64_t range)
{
    uint64_t offset = range == 0 ? chacha8_next64(rng) : chacha8_next_ranged(rng, range);
    return (int64_t)((uint64_t)start + offset);
}

static void
chacha8_reset(chacha8_rng *rng)
{
    rng->counter = 0;
    rng->used = 16;  // force immediate refill on first call
    rng->n = 16;
    rng->blocks_left = CHACHA8_RESEED_INTERVAL;
}

static int
chacha8_init(chacha8_rng *rng)
{
    if (fill_os_entropy(rng->seed, sizeof(rng->seed)) < 0)
        return -1;
    chacha8_reset(rng);
    return 0;
}

// Seed from an integer by hashing it through ChaCha8 to fill the 256-bit key.
static void
chacha8_seed_int(chacha8_rng *rng, int64_t seed_val)
{
    // Use the integer to construct a simple initial key, then run one
    // ChaCha8 block to diffuse it across all 256 bits.
    memset(rng->seed, 0, sizeof(rng->seed));
    rng->seed[0] = (uint32_t)(seed_val & 0xFFFFFFFF);
    rng->seed[1] = (uint32_t)((uint64_t)seed_val >> 32);

    uint32_t out[16];
    chacha8_block(rng->seed, 0, out);
    memcpy(rng->seed, out, sizeof(rng->seed));
    chacha8_reset(rng);
}

//
// Thread-local global RNG for module-level random()/randint()
//
// thread_local pointer for fast access (direct %fs/%gs-relative load),
// platform TLS key with destructor for cleanup on thread exit.
//

#ifdef _WIN32
static __declspec(thread) chacha8_rng *tls_rng = NULL;
#else
static __thread chacha8_rng *tls_rng = NULL;
#endif

#ifdef _WIN32
static DWORD tls_key = FLS_OUT_OF_INDEXES;

static void NTAPI
tls_rng_destructor(void *ptr)
{
    if (ptr != NULL) {
        memset(ptr, 0, sizeof(chacha8_rng));
        PyMem_RawFree(ptr);
    }
}
#else
static pthread_key_t tls_key;

static void
tls_rng_destructor(void *ptr)
{
    if (ptr != NULL) {
        memset(ptr, 0, sizeof(chacha8_rng));
        PyMem_RawFree(ptr);
    }
}
#endif

static int tls_key_created = 0;

static int
ensure_tls_key(void)
{
    if (likely(tls_key_created))
        return 0;
#ifdef _WIN32
    tls_key = FlsAlloc(tls_rng_destructor);
    if (tls_key == FLS_OUT_OF_INDEXES) {
        PyErr_SetString(PyExc_OSError, "FlsAlloc failed");
        return -1;
    }
#else
    if (pthread_key_create(&tls_key, tls_rng_destructor) != 0) {
        PyErr_SetString(PyExc_OSError, "pthread_key_create failed");
        return -1;
    }
#endif
    tls_key_created = 1;
    return 0;
}

// Get the thread-local RNG, initializing on first use.
// Returns NULL with Python exception set on failure.
static inline chacha8_rng *
get_thread_rng(void)
{
    chacha8_rng *rng = tls_rng;
    if (likely(rng != NULL))
        return rng;

    // First use on this thread — allocate and seed
    rng = PyMem_RawMalloc(sizeof(chacha8_rng));
    if (rng == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    if (chacha8_init(rng) < 0) {
        PyMem_RawFree(rng);
        return NULL;
    }

    // Register with platform TLS for destructor
#ifdef _WIN32
    FlsSetValue(tls_key, rng);
#else
    pthread_setspecific(tls_key, rng);
#endif

    tls_rng = rng;
    return rng;
}

// Return a random double in [0.0, 1.0) with 53 bits of mantissa precision.
static inline double
random_double_impl(chacha8_rng *rng)
{
    uint64_t r = chacha8_next64(rng);
    return (double)(r >> 11) * (1.0 / 9007199254740992.0);  // 1/2^53
}

//
// Module-level random() and randint()
//

static PyObject*
module_random(PyObject *module, PyObject *Py_UNUSED(ignored))
{
    chacha8_rng *rng = get_thread_rng();
    if (rng == NULL)
        return NULL;
    return PyFloat_FromDouble(random_double_impl(rng));
}

// Generate random integer in [a, b] using the given RNG.
static inline PyObject*
randint_impl(chacha8_rng *rng, int64_t a, int64_t b)
{
    uint64_t range = (uint64_t)b - (uint64_t)a + 1;
    return PyLong_FromLongLong(random_i64_from_range(rng, a, range));
}

static PyObject*
module_randint(PyObject *module, PyObject *const *args, Py_ssize_t nargs)
{
    if (nargs != 2) {
        PyErr_Format(PyExc_TypeError,
                     "randint() takes exactly 2 arguments (%zd given)", nargs);
        return NULL;
    }

    int64_t a = CPyLong_AsInt64(args[0]);
    if (unlikely(a == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;

    int64_t b = CPyLong_AsInt64(args[1]);
    if (unlikely(b == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;

    if (a > b) {
        PyErr_SetString(PyExc_ValueError,
                        "empty range for randint()");
        return NULL;
    }

    chacha8_rng *rng = get_thread_rng();
    if (rng == NULL)
        return NULL;

    return randint_impl(rng, a, b);
}

// Parse 1 or 2 int args for randrange([start,] stop).
// Sets *a to start (default 0), *b to stop-1.
// Returns 0 on success, -1 on error (with exception set).
static int
parse_randrange_args(PyObject *const *args, Py_ssize_t nargs,
                     int64_t *a, int64_t *b)
{
    if (nargs == 1) {
        *a = 0;
        int64_t stop = CPyLong_AsInt64(args[0]);
        if (unlikely(stop == CPY_LL_INT_ERROR && PyErr_Occurred()))
            return -1;
        if (stop <= 0) {
            PyErr_SetString(PyExc_ValueError, "empty range for randrange()");
            return -1;
        }
        *b = stop - 1;
    } else if (nargs == 2) {
        *a = CPyLong_AsInt64(args[0]);
        if (unlikely(*a == CPY_LL_INT_ERROR && PyErr_Occurred()))
            return -1;
        int64_t stop = CPyLong_AsInt64(args[1]);
        if (unlikely(stop == CPY_LL_INT_ERROR && PyErr_Occurred()))
            return -1;
        if (*a >= stop) {
            PyErr_SetString(PyExc_ValueError, "empty range for randrange()");
            return -1;
        }
        *b = stop - 1;
    } else {
        PyErr_Format(PyExc_TypeError,
                     "randrange() takes 1 or 2 arguments (%zd given)", nargs);
        return -1;
    }
    return 0;
}

static PyObject*
module_randrange(PyObject *module, PyObject *const *args, Py_ssize_t nargs)
{
    int64_t a, b;
    if (parse_randrange_args(args, nargs, &a, &b) < 0)
        return NULL;

    chacha8_rng *rng = get_thread_rng();
    if (rng == NULL)
        return NULL;

    return randint_impl(rng, a, b);
}

static PyObject*
module_seed(PyObject *module, PyObject *const *args, Py_ssize_t nargs)
{
    if (nargs != 1) {
        PyErr_Format(PyExc_TypeError,
                     "seed() takes exactly 1 argument (%zd given)", nargs);
        return NULL;
    }
    int64_t seed_val = CPyLong_AsInt64(args[0]);
    if (unlikely(seed_val == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;

    chacha8_rng *rng = get_thread_rng();
    if (rng == NULL)
        return NULL;

    chacha8_seed_int(rng, seed_val);
    Py_RETURN_NONE;
}

//
// Random Python type
//

typedef struct {
    PyObject_HEAD
    chacha8_rng rng;
} RandomObject;

static PyTypeObject RandomType;

static PyObject*
Random_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    if (type != &RandomType) {
        PyErr_SetString(PyExc_TypeError, "Random cannot be subclassed");
        return NULL;
    }

    RandomObject *self = (RandomObject *)type->tp_alloc(type, 0);
    // Seeding is done in tp_init
    return (PyObject *)self;
}

static int
Random_init(RandomObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *seed_obj = NULL;

    if (!PyArg_ParseTuple(args, "|O", &seed_obj)) {
        return -1;
    }

    if (kwds != NULL && PyDict_Size(kwds) > 0) {
        PyErr_SetString(PyExc_TypeError,
                        "Random() takes no keyword arguments");
        return -1;
    }

    if (seed_obj == NULL || seed_obj == Py_None) {
        if (chacha8_init(&self->rng) < 0)
            return -1;
    } else {
        int64_t seed_val = CPyLong_AsInt64(seed_obj);
        if (unlikely(seed_val == CPY_LL_INT_ERROR && PyErr_Occurred()))
            return -1;
        chacha8_seed_int(&self->rng, seed_val);
    }

    return 0;
}

// Internal constructors for capsule API (bypass tp_new/tp_init)

static PyObject *
Random_internal(void) {
    RandomObject *self = (RandomObject *)RandomType.tp_alloc(&RandomType, 0);
    if (self == NULL)
        return NULL;
    if (chacha8_init(&self->rng) < 0) {
        Py_DECREF(self);
        return NULL;
    }
    return (PyObject *)self;
}

static PyObject *
Random_from_seed_internal(int64_t seed_val) {
    RandomObject *self = (RandomObject *)RandomType.tp_alloc(&RandomType, 0);
    if (self == NULL)
        return NULL;
    chacha8_seed_int(&self->rng, seed_val);
    return (PyObject *)self;
}

static PyTypeObject *
Random_type_internal(void) {
    return &RandomType;
}

static int64_t
Random_randrange1_internal(PyObject *self, int64_t stop) {
    if (unlikely(stop <= 0)) {
        PyErr_SetString(PyExc_ValueError, "empty range for randrange()");
        return CPY_LL_INT_ERROR;
    }
    return (int64_t)chacha8_next_ranged(&((RandomObject *)self)->rng, (uint64_t)stop);
}

static int64_t
Random_randrange2_internal(PyObject *self, int64_t start, int64_t stop) {
    if (unlikely(start >= stop)) {
        PyErr_SetString(PyExc_ValueError, "empty range for randrange()");
        return CPY_LL_INT_ERROR;
    }
    uint64_t range = (uint64_t)stop - (uint64_t)start;
    return random_i64_from_range(&((RandomObject *)self)->rng, start, range);
}

static int64_t
Random_randint_internal(PyObject *self, int64_t a, int64_t b) {
    if (unlikely(a > b)) {
        PyErr_SetString(PyExc_ValueError, "empty range for randint()");
        return CPY_LL_INT_ERROR;
    }
    uint64_t range = (uint64_t)b - (uint64_t)a + 1;
    return random_i64_from_range(&((RandomObject *)self)->rng, a, range);
}

static double
Random_random_internal(PyObject *self) {
    return random_double_impl(&((RandomObject *)self)->rng);
}

static PyObject*
Random_randint(RandomObject *self, PyObject *const *args, Py_ssize_t nargs) {
    if (nargs != 2) {
        PyErr_Format(PyExc_TypeError,
                     "randint() takes exactly 2 arguments (%zd given)", nargs);
        return NULL;
    }

    int64_t a = CPyLong_AsInt64(args[0]);
    if (unlikely(a == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;

    int64_t b = CPyLong_AsInt64(args[1]);
    if (unlikely(b == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;

    if (a > b) {
        PyErr_SetString(PyExc_ValueError,
                        "empty range for randint()");
        return NULL;
    }

    return randint_impl(&self->rng, a, b);
}

static PyObject*
Random_randrange(RandomObject *self, PyObject *const *args, Py_ssize_t nargs) {
    int64_t a, b;
    if (parse_randrange_args(args, nargs, &a, &b) < 0)
        return NULL;
    return randint_impl(&self->rng, a, b);
}

static PyObject*
Random_random(RandomObject *self, PyObject *Py_UNUSED(ignored)) {
    return PyFloat_FromDouble(random_double_impl(&self->rng));
}

static PyObject*
Random_seed(RandomObject *self, PyObject *const *args, Py_ssize_t nargs) {
    if (nargs != 1) {
        PyErr_Format(PyExc_TypeError,
                     "seed() takes exactly 1 argument (%zd given)", nargs);
        return NULL;
    }
    int64_t seed_val = CPyLong_AsInt64(args[0]);
    if (unlikely(seed_val == CPY_LL_INT_ERROR && PyErr_Occurred()))
        return NULL;
    chacha8_seed_int(&self->rng, seed_val);
    Py_RETURN_NONE;
}

static PyMethodDef Random_methods[] = {
    {"randint", (PyCFunction) Random_randint, METH_FASTCALL,
     PyDoc_STR("Return random integer in range [a, b], including both end points.")
    },
    {"randrange", (PyCFunction) Random_randrange, METH_FASTCALL,
     PyDoc_STR("Return random integer in range [start, stop).")
    },
    {"random", (PyCFunction) Random_random, METH_NOARGS,
     PyDoc_STR("Return random float in [0.0, 1.0).")
    },
    {"seed", (PyCFunction) Random_seed, METH_FASTCALL,
     PyDoc_STR("Seed the random number generator with an integer.")
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject RandomType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Random",
    .tp_doc = PyDoc_STR("Fast random number generator using ChaCha8"),
    .tp_basicsize = sizeof(RandomObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = Random_new,
    .tp_init = (initproc) Random_init,
    .tp_methods = Random_methods,
};

// Module definition

static PyMethodDef librt_random_module_methods[] = {
    {"random", (PyCFunction) module_random, METH_NOARGS,
     PyDoc_STR("Return random float in [0.0, 1.0) using thread-local RNG.")
    },
    {"randint", (PyCFunction) module_randint, METH_FASTCALL,
     PyDoc_STR("Return random integer in range [a, b] using thread-local RNG.")
    },
    {"randrange", (PyCFunction) module_randrange, METH_FASTCALL,
     PyDoc_STR("Return random integer in range [start, stop) using thread-local RNG.")
    },
    {"seed", (PyCFunction) module_seed, METH_FASTCALL,
     PyDoc_STR("Seed the thread-local RNG with an integer.")
    },
    {NULL, NULL, 0, NULL}
};

// Module-level internal functions for mypyc primitives (use thread-local RNG)

static double
module_random_internal(void) {
    chacha8_rng *rng = get_thread_rng();
    if (rng == NULL)
        return CPY_FLOAT_ERROR;
    return random_double_impl(rng);
}

static int64_t
module_randint_internal(int64_t a, int64_t b) {
    if (unlikely(a > b)) {
        PyErr_SetString(PyExc_ValueError, "empty range for randint()");
        return CPY_LL_INT_ERROR;
    }
    chacha8_rng *rng = get_thread_rng();
    if (rng == NULL)
        return CPY_LL_INT_ERROR;
    uint64_t range = (uint64_t)b - (uint64_t)a + 1;
    return random_i64_from_range(rng, a, range);
}

static int64_t
module_randrange1_internal(int64_t stop) {
    if (unlikely(stop <= 0)) {
        PyErr_SetString(PyExc_ValueError, "empty range for randrange()");
        return CPY_LL_INT_ERROR;
    }
    chacha8_rng *rng = get_thread_rng();
    if (rng == NULL)
        return CPY_LL_INT_ERROR;
    return (int64_t)chacha8_next_ranged(rng, (uint64_t)stop);
}

static int64_t
module_randrange2_internal(int64_t start, int64_t stop) {
    if (unlikely(start >= stop)) {
        PyErr_SetString(PyExc_ValueError, "empty range for randrange()");
        return CPY_LL_INT_ERROR;
    }
    chacha8_rng *rng = get_thread_rng();
    if (rng == NULL)
        return CPY_LL_INT_ERROR;
    uint64_t range = (uint64_t)stop - (uint64_t)start;
    return random_i64_from_range(rng, start, range);
}

static int
random_abi_version(void) {
    return LIBRT_RANDOM_ABI_VERSION;
}

static int
random_api_version(void) {
    return LIBRT_RANDOM_API_VERSION;
}

static int
librt_random_module_exec(PyObject *m)
{
    if (ensure_tls_key() < 0) {
        return -1;
    }
    if (PyType_Ready(&RandomType) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "Random", (PyObject *) &RandomType) < 0) {
        return -1;
    }
    // Export mypyc internal C API via capsule
    static void *librt_random_api[LIBRT_RANDOM_API_LEN] = {
        (void *)random_abi_version,
        (void *)random_api_version,
        (void *)Random_internal,
        (void *)Random_from_seed_internal,
        (void *)Random_type_internal,
        (void *)Random_random_internal,
        (void *)Random_randint_internal,
        (void *)Random_randrange1_internal,
        (void *)Random_randrange2_internal,
        (void *)module_random_internal,
        (void *)module_randint_internal,
        (void *)module_randrange1_internal,
        (void *)module_randrange2_internal,
    };
    PyObject *c_api_object = PyCapsule_New((void *)librt_random_api, "librt.random._C_API", NULL);
    if (PyModule_Add(m, "_C_API", c_api_object) < 0) {
        return -1;
    }
    return 0;
}

static PyModuleDef_Slot librt_random_module_slots[] = {
    {Py_mod_exec, librt_random_module_exec},
#ifdef Py_MOD_GIL_NOT_USED
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static PyModuleDef librt_random_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "random",
    .m_doc = "Fast random number generation using ChaCha8",
    .m_size = 0,
    .m_methods = librt_random_module_methods,
    .m_slots = librt_random_module_slots,
};

PyMODINIT_FUNC
PyInit_random(void)
{
    return PyModuleDef_Init(&librt_random_module);
}
