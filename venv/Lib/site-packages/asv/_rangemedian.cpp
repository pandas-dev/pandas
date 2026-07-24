// Licensed under a 3-clause BSD style license - see LICENSE.rst

// Fast range median distance computations for dataset `y`:
//
//    mu(l, r) = median(y[l:r+1])
//    dist(l, r) = sum(abs(x - mu(l, r)) for x in y[l:r+1])
//
// and an implementation of the find-best-partition dynamic program.
//
// We don't implement a rolling median computation, on the assumption that
// accesses are concentrated on small windows in the data.

#include <vector>
#include <queue>
#include <limits>
#include <map>
#include <mutex>
#include <algorithm>
#include <utility>

#include <Python.h>



//
// Median computation.
//

template <class const_iterator>
void compute_weighted_median(const_iterator start, const_iterator end,
                             double *mu, double *dist)
{
    std::vector<std::pair<double,double> > tmp;
    std::vector<std::pair<double,double> >::iterator it;
    double midpoint, wsum;

    if (start == end) {
        *mu = 0;
        *dist = 0;
        return;
    }

    tmp.insert(tmp.end(), start, end);
    std::sort(tmp.begin(), tmp.end());

    midpoint = 0;
    for (it = tmp.begin(); it != tmp.end(); ++it) {
        midpoint += it->second;
    }
    midpoint /= 2;

    wsum = 0;
    for (it = tmp.begin(); it != tmp.end(); ++it) {
        wsum += it->second;
        if (wsum >= midpoint) {
            break;
        }
    }

    if (it != tmp.end()) {
        *mu = it->first;
        if (wsum == midpoint) {
            ++it;
            if (it != tmp.end()) {
                *mu = (it->first + *mu) / 2;
            }
        }
    }
    else {
        // Error condition, maybe some floating point summation issue
        --it;
        *mu = it->first;
    }

    *dist = 0;
    for (const_iterator it = start; it < end; ++it) {
        *dist += it->second * fabs(it->first - *mu);
    }
}


//
// Cache for cache[left,right] == (mu, dist)
//

class Cache
{
private:
    struct Item
    {
        double mu, dist;
    };

    std::map<std::pair<size_t, size_t>, Item> items_;

public:
    Cache() {};
    bool get(size_t left, size_t right, double *mu, double *dist) const {
        auto it = items_.find(std::make_pair(left, right));
        if (it != items_.end()) {
            *mu = it->second.mu;
            *dist = it->second.dist;
            return true;
        }
        return false;
    }

    void set(size_t left, size_t right, double mu, double dist) {
        items_[std::make_pair(left, right)] = { mu, dist };
    }
};

// Module state
typedef struct {
    PyObject *RangeMedian_Type;    // Xxo class
} rangemedian_state;


//
// RangeMedian object.
//

typedef struct {
    PyObject_HEAD
    std::vector<std::pair<double,double> > *y;
    Cache *cache;
#ifdef Py_GIL_DISABLED
    std::mutex *use_mtx;
#endif
} RangeMedianObject;

#define RangeMedianObject_CAST(op)  ((RangeMedianObject *)(op))

PyObject *RangeMedian_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    RangeMedianObject *self;
    allocfunc rmalloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    self = (RangeMedianObject*)rmalloc(type, 0);
    if (self == NULL) {
        return NULL;
    }
    self->y = NULL;
    self->cache = NULL;
#ifdef Py_GIL_DISABLED
    self->use_mtx = NULL;
    try {
        self->use_mtx = new std::mutex();
    }
    catch (const std::bad_alloc&) {
        PyErr_SetString(PyExc_MemoryError, "Allocating memory failed");
        freefunc free = (freefunc)PyType_GetSlot(type, Py_tp_free);
        free((PyObject*)self);
        return NULL;
    }
#endif
    return (PyObject*)self;
}


#ifdef Py_GIL_DISABLED
class RangeMedianUseGuard
{
private:
    std::unique_lock<std::mutex> lock_;

public:
    explicit RangeMedianUseGuard(std::mutex *mtx) : lock_(*mtx, std::try_to_lock) {}

    bool acquired() const {
        return lock_.owns_lock();
    }
};

static bool
RangeMedian_try_acquire_use_lock(RangeMedianObject *self, RangeMedianUseGuard *guard)
{
    if (!guard->acquired()) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "RangeMedian objects cannot be used concurrently across threads"
        );
        return false;
    }
    return true;
}
#endif


int RangeMedian_init(PyObject *op, PyObject *args, PyObject *kwds)
{
    RangeMedianObject *self = RangeMedianObject_CAST(op);
    static const char *kwlist[] = {"y", "w", NULL};
    PyObject *y_obj, *w_obj;
    Py_ssize_t size, wsize, k;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!", (char**)kwlist,
                                     &PyList_Type, &y_obj,
                                     &PyList_Type, &w_obj)) {
        return -1;
    }

#ifdef Py_GIL_DISABLED
    RangeMedianUseGuard guard(self->use_mtx);
    if (!RangeMedian_try_acquire_use_lock(self, &guard)) {
        return -1;
    }
#endif

    size = PyList_Size(y_obj);
    wsize = PyList_Size(w_obj);

    if (wsize != size) {
        PyErr_SetString(PyExc_ValueError, "y and w must have same length");
        return -1;
    }

    // Free any previous data from a prior __init__ call
    delete self->y;
    delete self->cache;
    self->y = NULL;
    self->cache = NULL;

    try {
        self->y = new std::vector<std::pair<double,double> >(size);
        self->cache = new Cache();
    }
    catch (const std::bad_alloc&) {
        PyErr_SetString(PyExc_MemoryError, "Allocating memory failed");
        return -1;
    }

    for (k = 0; k < size; ++k) {
        PyObject *x, *wx;

        x = PyNumber_Float(PyList_GetItem(y_obj, k));
        if (x == NULL || !PyFloat_Check(x)) {
            Py_XDECREF(x);
            return -1;
        }

        wx = PyNumber_Float(PyList_GetItem(w_obj, k));
        if (wx == NULL || !PyFloat_Check(wx)) {
            Py_XDECREF(x);
            Py_XDECREF(wx);
            return -1;
        }

        (*self->y)[k] = std::make_pair(PyFloat_AsDouble(x),
                                       PyFloat_AsDouble(wx));
        Py_DECREF(x);
        Py_DECREF(wx);
    }

    return 0;
}
/* finalization.
 *
 * Types that store references to other PyObjects generally need to implement
 * the GC slots: traverse, clear, dealloc, and (optionally) finalize.
 */

// traverse: Visit all references from an object, including its type
static int
RangeMedian_traverse(PyObject *op, visitproc visit, void *arg)
{
    // Visit the type
    Py_VISIT(Py_TYPE(op));
    return 0;
}

// clear: drop references in order to break all reference cycles
static int
RangeMedian_clear(PyObject *op)
{
    return 0;
}


static void RangeMedian_finalize(void *op)
{
    RangeMedianObject *self = RangeMedianObject_CAST(op);
    delete self->y;
    delete self->cache;
#ifdef Py_GIL_DISABLED
    delete self->use_mtx;
#endif
}


static void RangeMedian_dealloc(void *self)
{
    PyObject_GC_UnTrack(self);
    RangeMedian_finalize(self);
    PyTypeObject *tp = Py_TYPE(self);
    freefunc free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    free(self);
    Py_DECREF(tp);
}


static int RangeMedian_mu_dist(PyObject *op, Py_ssize_t left, Py_ssize_t right,
                               double *mu, double *dist)
{
    RangeMedianObject *self = RangeMedianObject_CAST(op);
    Py_ssize_t size = (Py_ssize_t)self->y->size();

    if (left < 0 || right < 0 || left >= size || right >= size) {
        PyErr_SetString(PyExc_ValueError, "argument out of range");
        return -1;
    }

    if (!self->cache->get(left, right, mu, dist)) {
        compute_weighted_median(self->y->begin() + left, self->y->begin() + right + 1, mu, dist);
        self->cache->set(left, right, *mu, *dist);
    }

    return 0;
}


static PyObject *RangeMedian_mu(PyObject *op, PyObject *args)
{
    Py_ssize_t left, right;
    double mu = 0, dist;

    if (!PyArg_ParseTuple(args, "nn", &left, &right)) {
        return NULL;
    }

#ifdef Py_GIL_DISABLED
    RangeMedianObject *self = RangeMedianObject_CAST(op);
    RangeMedianUseGuard guard(self->use_mtx);
    if (!RangeMedian_try_acquire_use_lock(self, &guard)) {
        return NULL;
    }
#endif

    if (RangeMedian_mu_dist(op, left, right, &mu, &dist) == -1) {
        return NULL;
    }

    return PyFloat_FromDouble(mu);
}


static PyObject *RangeMedian_dist(PyObject *op, PyObject *args)
{
    Py_ssize_t left, right;
    double mu, dist = 0;

    if (!PyArg_ParseTuple(args, "nn", &left, &right)) {
        return NULL;
    }

#ifdef Py_GIL_DISABLED
    RangeMedianObject *self = RangeMedianObject_CAST(op);
    RangeMedianUseGuard guard(self->use_mtx);
    if (!RangeMedian_try_acquire_use_lock(self, &guard)) {
        return NULL;
    }
#endif

    if (RangeMedian_mu_dist(op, left, right, &mu, &dist) == -1) {
        return NULL;
    }

    return PyFloat_FromDouble(dist);
}


static PyObject *RangeMedian_find_best_partition(PyObject *op, PyObject *args)
{
    RangeMedianObject *self = RangeMedianObject_CAST(op);
    Py_ssize_t min_size, max_size, min_pos, max_pos;
    double gamma;
    Py_ssize_t size;

    if (!PyArg_ParseTuple(args, "dnnnn", &gamma, &min_size, &max_size, &min_pos, &max_pos)) {
        return NULL;
    }

#ifdef Py_GIL_DISABLED
    RangeMedianUseGuard guard(self->use_mtx);
    if (!RangeMedian_try_acquire_use_lock(self, &guard)) {
        return NULL;
    }
#endif

    size = self->y->size();

    if (!(0 < min_size && min_size <= max_size &&
          0 <= min_pos && min_pos <= max_pos && max_pos <= size)) {
        PyErr_SetString(PyExc_ValueError, "invalid input indices");
        return NULL;
    }

    double inf = std::numeric_limits<double>::infinity();

    std::vector<double> B(max_pos - min_pos + 1);
    std::vector<Py_ssize_t> p(max_pos - min_pos);

    B[0] = -gamma;

    for (Py_ssize_t right = min_pos; right < max_pos; ++right) {
        B[right + 1 - min_pos] = inf;

        Py_ssize_t aa = std::max(right + 1 - max_size, min_pos);
        Py_ssize_t bb = std::max(right + 1 - min_size + 1, min_pos);
        for (Py_ssize_t left = aa; left < bb; ++left) {
            double mu, dist;
            if (RangeMedian_mu_dist(op, left, right, &mu, &dist) == -1) {
                return NULL;
            }

            double b = B[left - min_pos] + gamma + dist;
            if (b <= B[right + 1 - min_pos]) {
                B[right + 1 - min_pos] = b;
                p[right - min_pos] = left - 1;
            }
        }
    }

    PyObject *p_list;

    p_list = PyList_New(p.size());
    if (p_list == NULL) {
        return NULL;
    }

    for (Py_ssize_t k = 0; k < (Py_ssize_t)p.size(); ++k) {
        PyObject *num = PyLong_FromSsize_t(p[k]);
        if (num == NULL) {
            Py_DECREF(p_list);
            return NULL;
        }
        PyList_SetItem(p_list, k, num);
    }

    return p_list;
}


//
// RangeMedian type.
//

static PyMethodDef RangeMedian_methods[] = {
    {"mu", (PyCFunction)RangeMedian_mu, METH_VARARGS, NULL},
    {"dist", (PyCFunction)RangeMedian_dist, METH_VARARGS, NULL},
    {"find_best_partition", (PyCFunction)RangeMedian_find_best_partition, METH_VARARGS, NULL},
    {NULL, NULL}
};

static PyType_Slot RangeMedian_Type_slots[] = {
    {Py_tp_dealloc, (void*)RangeMedian_dealloc},
    {Py_tp_traverse, (void*)RangeMedian_traverse},
    {Py_tp_clear, (void*)RangeMedian_clear},
    {Py_tp_finalize, (void*)RangeMedian_finalize},
    {Py_tp_methods, RangeMedian_methods},
    {Py_tp_init, (void*)RangeMedian_init},
    {Py_tp_new, (void*)RangeMedian_new},
    {0, 0},  /* sentinel */
};

static PyType_Spec RangeMedian_Type_spec = {
    "_rangemedian.RangeMedian",  // name
    sizeof(RangeMedianObject),   // basicsize
    0,                           // itemsize
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, // flags
    RangeMedian_Type_slots     // slots
};


PyDoc_STRVAR(module_doc,
"expose the RangeMedian type");


static int
rangemedian_modexec(PyObject *m)
{
    rangemedian_state *state = (rangemedian_state *)PyModule_GetState(m);

    state->RangeMedian_Type = PyType_FromSpec(&RangeMedian_Type_spec);
    if (state->RangeMedian_Type == NULL) {
        return -1;
    }
    if (PyType_Ready((PyTypeObject *)state->RangeMedian_Type) < 0) {
        return -1;
    }
    Py_INCREF(state->RangeMedian_Type);
    if (PyModule_AddObject(m, "RangeMedian", (PyObject *)state->RangeMedian_Type) < 0) {
        return -1;
    }

    return 0;
}

static PyModuleDef_Slot rangemedian_slots[] = {

    /* exec function to initialize the module (called as part of import
     * after the object was added to sys.modules)
     */
    {Py_mod_exec, (void *)rangemedian_modexec},

#ifdef Py_mod_gil
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif

    {0, NULL}
};


// Module finalization: modules that hold references in their module state
// need to implement the fullowing GC hooks. They're similar to the ones for
// types (see "finalization").

static int
rangemedian_traverse(PyObject *module, visitproc visit, void *arg)
{
    rangemedian_state *state = (rangemedian_state *)PyModule_GetState(module);
    Py_VISIT(state->RangeMedian_Type);
    return 0;
}

static int
rangemedian_clear(PyObject *module)
{
    rangemedian_state *state = (rangemedian_state *)PyModule_GetState(module);
    Py_CLEAR(state->RangeMedian_Type);
    return 0;
}

static void
rangemedian_free(void *module)
{
    // allow modexec to omit calling clear on error
    (void)rangemedian_clear((PyObject *)module);
}

//
// Module initialization.
//

extern "C" {

static struct PyModuleDef moduledef = {
   PyModuleDef_HEAD_INIT,
   "_rangemedian",      // m_name
   module_doc,          // m_doc
   sizeof(rangemedian_state), // m_size
   NULL,                 // m_methods
   rangemedian_slots,    // m_slots
   rangemedian_traverse, // m_traverse
   rangemedian_clear,    // m_clear
   rangemedian_free      // m_free
};

PyMODINIT_FUNC
PyInit__rangemedian(void)
{
    return PyModuleDef_Init(&moduledef);
}

} // extern "C"
