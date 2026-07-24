#ifndef LIBRT_RANDOM_API_H
#define LIBRT_RANDOM_API_H

#include <stdbool.h>
#include <stdint.h>
#include <Python.h>
#include "librt_random.h"

int
import_librt_random(void);

extern void *LibRTRandom_API[LIBRT_RANDOM_API_LEN];

#define LibRTRandom_ABIVersion (*(int (*)(void)) LibRTRandom_API[0])
#define LibRTRandom_APIVersion (*(int (*)(void)) LibRTRandom_API[1])
#define LibRTRandom_Random_internal (*(PyObject* (*)(void)) LibRTRandom_API[2])
#define LibRTRandom_Random_from_seed_internal (*(PyObject* (*)(int64_t)) LibRTRandom_API[3])
#define LibRTRandom_Random_type_internal (*(PyTypeObject* (*)(void)) LibRTRandom_API[4])
#define LibRTRandom_Random_random_internal (*(double (*)(PyObject*)) LibRTRandom_API[5])
#define LibRTRandom_Random_randint_internal (*(int64_t (*)(PyObject*, int64_t, int64_t)) LibRTRandom_API[6])
#define LibRTRandom_Random_randrange1_internal (*(int64_t (*)(PyObject*, int64_t)) LibRTRandom_API[7])
#define LibRTRandom_Random_randrange2_internal (*(int64_t (*)(PyObject*, int64_t, int64_t)) LibRTRandom_API[8])
#define LibRTRandom_module_random_internal (*(double (*)(void)) LibRTRandom_API[9])
#define LibRTRandom_module_randint_internal (*(int64_t (*)(int64_t, int64_t)) LibRTRandom_API[10])
#define LibRTRandom_module_randrange1_internal (*(int64_t (*)(int64_t)) LibRTRandom_API[11])
#define LibRTRandom_module_randrange2_internal (*(int64_t (*)(int64_t, int64_t)) LibRTRandom_API[12])

static inline bool CPyRandom_Check(PyObject *obj) {
    return Py_TYPE(obj) == LibRTRandom_Random_type_internal();
}

#endif  // LIBRT_RANDOM_API_H
