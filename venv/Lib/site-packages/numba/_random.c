/*
 * PRNG support.
 */

#ifdef _MSC_VER
#define HAVE_PTHREAD_ATFORK 0
#else
#define HAVE_PTHREAD_ATFORK 1
#include <pthread.h>
#endif


/* Magic Mersenne Twister constants */
#define MT_N 624
#define MT_M 397
#define MT_MATRIX_A 0x9908b0dfU
#define MT_UPPER_MASK 0x80000000U
#define MT_LOWER_MASK 0x7fffffffU

/*
 * Note this structure is accessed in numba.targets.randomimpl,
 * any changes here should be reflected there too.
 */
typedef struct {
    int index;
    /* unsigned int is sufficient on modern machines as we only need 32 bits */
    unsigned int mt[MT_N];
    int has_gauss;
    double gauss;
    int is_initialized;
} rnd_state_t;

/* Some code portions below from CPython's _randommodule.c, some others
   from Numpy's and Jean-Sebastien Roy's randomkit.c. */

NUMBA_EXPORT_FUNC(void)
numba_rnd_shuffle(rnd_state_t *state)
{
    int i;
    unsigned int y;

    for (i = 0; i < MT_N - MT_M; i++) {
        y = (state->mt[i] & MT_UPPER_MASK) | (state->mt[i+1] & MT_LOWER_MASK);
        state->mt[i] = state->mt[i+MT_M] ^ (y >> 1) ^
                       (-(int) (y & 1) & MT_MATRIX_A);
    }
    for (; i < MT_N - 1; i++) {
        y = (state->mt[i] & MT_UPPER_MASK) | (state->mt[i+1] & MT_LOWER_MASK);
        state->mt[i] = state->mt[i+(MT_M-MT_N)] ^ (y >> 1) ^
                       (-(int) (y & 1) & MT_MATRIX_A);
    }
    y = (state->mt[MT_N - 1] & MT_UPPER_MASK) | (state->mt[0] & MT_LOWER_MASK);
    state->mt[MT_N - 1] = state->mt[MT_M - 1] ^ (y >> 1) ^
                          (-(int) (y & 1) & MT_MATRIX_A);
}

/* Initialize mt[] with an integer seed */
NUMBA_EXPORT_FUNC(void)
numba_rnd_init(rnd_state_t *state, unsigned int seed)
{
    unsigned int pos;
    seed &= 0xffffffffU;

    /* Knuth's PRNG as used in the Mersenne Twister reference implementation */
    for (pos = 0; pos < MT_N; pos++) {
        state->mt[pos] = seed;
        seed = (1812433253U * (seed ^ (seed >> 30)) + pos + 1) & 0xffffffffU;
    }
    state->index = MT_N;
    state->has_gauss = 0;
    state->gauss = 0.0;
    state->is_initialized = 1;
}

/* Perturb mt[] with a key array */
static void
rnd_init_by_array(rnd_state_t *state, unsigned int init_key[], size_t key_length)
{
    size_t i, j, k;
    unsigned int *mt = state->mt;

    numba_rnd_init(state, 19650218U);
    i = 1; j = 0;
    k = (MT_N > key_length ? MT_N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525U))
                 + init_key[j] + (unsigned int) j; /* non linear */
        mt[i] &= 0xffffffffU;
        i++; j++;
        if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
        if (j >= key_length) j = 0;
    }
    for (k = MT_N - 1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941U))
                 - (unsigned int) i; /* non linear */
        mt[i] &= 0xffffffffU;
        i++;
        if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i=1; }
    }

    mt[0] = 0x80000000U; /* MSB is 1; ensuring non-zero initial array */
    state->index = MT_N;
    state->has_gauss = 0;
    state->gauss = 0.0;
    state->is_initialized = 1;
}

/*
 * Management of thread-local random state.
 */

static int rnd_globally_initialized;

#ifdef _MSC_VER
#define THREAD_LOCAL(ty) __declspec(thread) ty
#else
/* Non-standard C99 extension that's understood by gcc and clang */
#define THREAD_LOCAL(ty) __thread ty
#endif

static THREAD_LOCAL(rnd_state_t) numba_py_random_state;
static THREAD_LOCAL(rnd_state_t) numba_np_random_state;
static THREAD_LOCAL(rnd_state_t) numba_internal_random_state;

/* Seed the state with random bytes */
static int
rnd_seed_with_bytes(rnd_state_t *state, Py_buffer *buf)
{
    unsigned int *keys;
    unsigned char *bytes;
    size_t i, nkeys;

    nkeys = buf->len / sizeof(unsigned int);
    keys = (unsigned int *) PyMem_Malloc(nkeys * sizeof(unsigned int));
    if (keys == NULL) {
        PyBuffer_Release(buf);
        return -1;
    }
    bytes = (unsigned char *) buf->buf;
    /* Convert input bytes to int32 keys, without violating alignment
     * constraints.
     */
    for (i = 0; i < nkeys; i++, bytes += 4) {
        keys[i] =
            ((unsigned int)bytes[3] << 24) +
            ((unsigned int)bytes[2] << 16) +
            ((unsigned int)bytes[1] << 8) +
            ((unsigned int)bytes[0] << 0);
    }
    PyBuffer_Release(buf);
    rnd_init_by_array(state, keys, nkeys);
    PyMem_Free(keys);
    return 0;
}

#if HAVE_PTHREAD_ATFORK
/* After a fork(), the child should reseed its random states.
 * Since only the main thread survives in the child, it's enough to mark
 * the current thread-local states as uninitialized.
 */
static void
rnd_atfork_child(void)
{
    numba_py_random_state.is_initialized = 0;
    numba_np_random_state.is_initialized = 0;
    numba_internal_random_state.is_initialized = 0;
}
#endif

/* Global initialization routine.  It must be called as early as possible.
 */
NUMBA_EXPORT_FUNC(void)
numba_rnd_ensure_global_init(void)
{
    if (!rnd_globally_initialized) {
#if HAVE_PTHREAD_ATFORK
        pthread_atfork(NULL, NULL, rnd_atfork_child);
#endif
        numba_py_random_state.is_initialized = 0;
        numba_np_random_state.is_initialized = 0;
        numba_internal_random_state.is_initialized = 0;
        rnd_globally_initialized = 1;
    }
}

/* First-time init a random state */
static void
rnd_implicit_init(rnd_state_t *state)
{
    /* Initialize with random bytes.  The easiest way to get good-quality
     * cross-platform random bytes is still to call os.urandom()
     * using the Python interpreter...
     */
    PyObject *module, *bufobj;
    Py_buffer buf;
    PyGILState_STATE gilstate = PyGILState_Ensure();

    module = PyImport_ImportModule("os");
    if (module == NULL)
        goto error;
    /* Read as many bytes as necessary to get the full entropy
     * exploitable by the MT generator.
     */
    bufobj = PyObject_CallMethod(module, "urandom", "i",
                                 (int) (MT_N * sizeof(unsigned int)));
    Py_DECREF(module);
    if (bufobj == NULL)
        goto error;
    if (PyObject_GetBuffer(bufobj, &buf, PyBUF_SIMPLE))
        goto error;
    Py_DECREF(bufobj);
    if (rnd_seed_with_bytes(state, &buf))
        goto error;
    /* state->is_initialized is set now */

    PyGILState_Release(gilstate);
    return;

error:
    /* In normal conditions, os.urandom() and PyMem_Malloc() shouldn't fail,
     * and we don't want the caller to deal with errors, so just bail out.
     */
    if (PyErr_Occurred())
        PyErr_Print();
    Py_FatalError(NULL);
}

/* Functions returning the thread-local random state pointer.
 * The LLVM JIT doesn't support thread-local variables so we rely
 * on the C compiler instead.
 */

NUMBA_EXPORT_FUNC(rnd_state_t *)
numba_get_py_random_state(void)
{
    rnd_state_t *state = &numba_py_random_state;
    if (!state->is_initialized)
        rnd_implicit_init(state);
    return state;
}

NUMBA_EXPORT_FUNC(rnd_state_t *)
numba_get_np_random_state(void)
{
    rnd_state_t *state = &numba_np_random_state;
    if (!state->is_initialized)
        rnd_implicit_init(state);
    return state;
}

NUMBA_EXPORT_FUNC(rnd_state_t *)
numba_get_internal_random_state(void)
{
    rnd_state_t *state = &numba_internal_random_state;
    if (!state->is_initialized)
        rnd_implicit_init(state);
    return state;
}

/*
 * Python-exposed helpers for state management and testing.
 */
static int
rnd_state_converter(PyObject *obj, rnd_state_t **state)
{
    *state = (rnd_state_t *) PyLong_AsVoidPtr(obj);
    return (*state != NULL || !PyErr_Occurred());
}

NUMBA_EXPORT_FUNC(PyObject *)
_numba_rnd_get_py_state_ptr(PyObject *self)
{
    return PyLong_FromVoidPtr(numba_get_py_random_state());
}

NUMBA_EXPORT_FUNC(PyObject *)
_numba_rnd_get_np_state_ptr(PyObject *self)
{
    return PyLong_FromVoidPtr(numba_get_np_random_state());
}

NUMBA_EXPORT_FUNC(PyObject *)
_numba_rnd_shuffle(PyObject *self, PyObject *arg)
{
    rnd_state_t *state;
    if (!rnd_state_converter(arg, &state))
        return NULL;
    numba_rnd_shuffle(state);
    Py_RETURN_NONE;
}

NUMBA_EXPORT_FUNC(PyObject *)
_numba_rnd_set_state(PyObject *self, PyObject *args)
{
    int i, index;
    rnd_state_t *state;
    PyObject *tuplearg, *intlist;

    if (!PyArg_ParseTuple(args, "O&O!:rnd_set_state",
                          rnd_state_converter, &state,
                          &PyTuple_Type, &tuplearg))
        return NULL;
    if (!PyArg_ParseTuple(tuplearg, "iO!", &index, &PyList_Type, &intlist))
        return NULL;
    if (PyList_GET_SIZE(intlist) != MT_N) {
        PyErr_SetString(PyExc_ValueError, "list object has wrong size");
        return NULL;
    }
    state->index = index;
    for (i = 0; i < MT_N; i++) {
        PyObject *v = PyList_GET_ITEM(intlist, i);
        unsigned long x = PyLong_AsUnsignedLong(v);
        if (x == (unsigned long) -1 && PyErr_Occurred())
            return NULL;
        state->mt[i] = (unsigned int) x;
    }
    state->has_gauss = 0;
    state->gauss = 0.0;
    state->is_initialized = 1;
    Py_RETURN_NONE;
}

NUMBA_EXPORT_FUNC(PyObject *)
_numba_rnd_get_state(PyObject *self, PyObject *arg)
{
    PyObject *intlist;
    int i;
    rnd_state_t *state;
    if (!rnd_state_converter(arg, &state))
        return NULL;

    intlist = PyList_New(MT_N);
    if (intlist == NULL)
        return NULL;
    for (i = 0; i < MT_N; i++) {
        PyObject *v = PyLong_FromUnsignedLong(state->mt[i]);
        if (v == NULL) {
            Py_DECREF(intlist);
            return NULL;
        }
        PyList_SET_ITEM(intlist, i, v);
    }
    return Py_BuildValue("iN", state->index, intlist);
}

NUMBA_EXPORT_FUNC(PyObject *)
_numba_rnd_seed(PyObject *self, PyObject *args)
{
    unsigned int seed;
    rnd_state_t *state;

    if (!PyArg_ParseTuple(args, "O&I:rnd_seed",
                          rnd_state_converter, &state, &seed)) {
        /* rnd_seed_*(bytes-like object) */
        Py_buffer buf;

        PyErr_Clear();
        if (!PyArg_ParseTuple(args, "O&s*:rnd_seed",
                              rnd_state_converter, &state, &buf))
            return NULL;

        if (rnd_seed_with_bytes(state, &buf))
            return NULL;
        else
            Py_RETURN_NONE;
    }
    else {
        /* rnd_seed_*(int32) */
        numba_rnd_init(state, seed);
        Py_RETURN_NONE;
    }
}

/*
 * Random distribution helpers.
 * Most code straight from Numpy's distributions.c.
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

NUMBA_EXPORT_FUNC(unsigned int)
get_next_int32(rnd_state_t *state)
{
    unsigned int y;

    if (state->index == MT_N) {
        numba_rnd_shuffle(state);
        state->index = 0;
    }
    y = state->mt[state->index++];
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680U;
    y ^= (y << 15) & 0xefc60000U;
    y ^= (y >> 18);
    return y;
}

NUMBA_EXPORT_FUNC(double)
get_next_double(rnd_state_t *state)
{
    double a = get_next_int32(state) >> 5;
    double b = get_next_int32(state) >> 6;
    return (a * 67108864.0 + b) / 9007199254740992.0;
}

NUMBA_EXPORT_FUNC(double)
loggam(double x)
{
    double x0, x2, xp, gl, gl0;
    long k, n;

    static double a[10] = {8.333333333333333e-02,-2.777777777777778e-03,
         7.936507936507937e-04,-5.952380952380952e-04,
         8.417508417508418e-04,-1.917526917526918e-03,
         6.410256410256410e-03,-2.955065359477124e-02,
         1.796443723688307e-01,-1.39243221690590e+00};
    x0 = x;
    n = 0;
    if ((x == 1.0) || (x == 2.0))
    {
        return 0.0;
    }
    else if (x <= 7.0)
    {
        n = (long)(7 - x);
        x0 = x + n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2*M_PI;
    gl0 = a[9];
    for (k=8; k>=0; k--)
    {
        gl0 *= x2;
        gl0 += a[k];
    }
    gl = gl0/x0 + 0.5*log(xp) + (x0-0.5)*log(x0) - x0;
    if (x <= 7.0)
    {
        for (k=1; k<=n; k++)
        {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}


NUMBA_EXPORT_FUNC(int64_t)
numba_poisson_ptrs(rnd_state_t *state, double lam)
{
    /* This method is invoked only if the parameter lambda of this
     * distribution is big enough ( >= 10 ). The algorithm used is
     * described in "Hörmann, W. 1992. 'The Transformed Rejection
     * Method for Generating Poisson Random Variables'.
     * The implementation comes straight from Numpy.
     */
    int64_t k;
    double U, V, slam, loglam, a, b, invalpha, vr, us;

    slam = sqrt(lam);
    loglam = log(lam);
    b = 0.931 + 2.53*slam;
    a = -0.059 + 0.02483*b;
    invalpha = 1.1239 + 1.1328/(b-3.4);
    vr = 0.9277 - 3.6224/(b-2);

    while (1)
    {
        U = get_next_double(state) - 0.5;
        V = get_next_double(state);
        us = 0.5 - fabs(U);
        k = (int64_t) floor((2*a/us + b)*U + lam + 0.43);
        if ((us >= 0.07) && (V <= vr))
        {
            return k;
        }
        if ((k < 0) ||
            ((us < 0.013) && (V > us)))
        {
            continue;
        }
        if ((log(V) + log(invalpha) - log(a/(us*us)+b)) <=
            (-lam + (double) k*loglam - loggam((double) k+1)))
        {
            return k;
        }
    }
}
