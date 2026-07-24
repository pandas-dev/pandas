#ifndef STATIC_DATA_H
#define STATIC_DATA_H

#include <Python.h>

#ifdef __cplusplus
extern "C" {
#endif

// Adopted from numpy 2.4.0: numpy/_core/src/multiarry/npy_static_data.h

int intern_strings(void);

typedef struct mypyc_interned_str_struct {
    PyObject *__init_subclass__;
    PyObject *__module__;
    PyObject *__mro_entries__;
    PyObject *__mypyc_attrs__;
    PyObject *__orig_bases__;
    PyObject *__qualname__;
    PyObject *__slots__;
    PyObject *__name__;
    PyObject *__radd__;
    PyObject *__rsub__;
    PyObject *__rmul__;
    PyObject *__rtruediv__;
    PyObject *__rmod__;
    PyObject *__rdivmod__;
    PyObject *__rfloordiv__;
    PyObject *__rpow__;
    PyObject *__rmatmul__;
    PyObject *__rand__;
    PyObject *__ror__;
    PyObject *__rxor__;
    PyObject *__rlshift__;
    PyObject *__rrshift__;
    PyObject *__eq__;
    PyObject *__ne__;
    PyObject *__gt__;
    PyObject *__le__;
    PyObject *__lt__;
    PyObject *__ge__;
    PyObject *clear;
    PyObject *close_;
    PyObject *copy;
    PyObject *dispatch_cache;
    PyObject *endswith;
    PyObject *get_type_hints;
    PyObject *keys;
    PyObject *lower;
    PyObject *items;
    PyObject *join;
    PyObject *register_;
    PyObject *registry;
    PyObject *send;
    PyObject *setdefault;
    PyObject *startswith;
    PyObject *super;
    PyObject *throw_;
    PyObject *translate;
    PyObject *update;
    PyObject *upper;
    PyObject *values;
} mypyc_interned_str_struct;

extern mypyc_interned_str_struct mypyc_interned_str;

#ifdef __cplusplus
}
#endif

#endif
