#ifndef LIBRT_VECS_API_H
#define LIBRT_VECS_API_H

#include "librt_vecs.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>

int
import_librt_vecs(void);

// Global API pointers initialized by import_librt_vecs()
extern VecCapsule *VecApi;
extern VecI64API VecI64Api;
extern VecI32API VecI32Api;
extern VecI16API VecI16Api;
extern VecU8API VecU8Api;
extern VecFloatAPI VecFloatApi;
extern VecBoolAPI VecBoolApi;
extern VecTAPI VecTApi;
extern VecNestedAPI VecNestedApi;

#endif  // LIBRT_VECS_API_H
