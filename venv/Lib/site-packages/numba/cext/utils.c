#include "cext.h"

/* Align size *sz* to pointer width */
Py_ssize_t
aligned_size(Py_ssize_t sz) {
    Py_ssize_t alignment = sizeof(void*);
    return sz + (alignment - sz % alignment) % alignment;
}
