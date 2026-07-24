#include "librt_vecs_api.h"

VecCapsule *VecApi = NULL;
VecI64API VecI64Api = {0};
VecI32API VecI32Api = {0};
VecI16API VecI16Api = {0};
VecU8API VecU8Api = {0};
VecFloatAPI VecFloatApi = {0};
VecBoolAPI VecBoolApi = {0};
VecTAPI VecTApi = {0};
VecNestedAPI VecNestedApi = {0};

int
import_librt_vecs(void)
{
    PyObject *mod = PyImport_ImportModule("librt.vecs");
    if (mod == NULL)
        return -1;
    Py_DECREF(mod);  // we import just for the side effect of making the below work.
    VecCapsule *capsule = PyCapsule_Import("librt.vecs._C_API", 0);
    if (!capsule)
        return -1;
    if (capsule->abi_version() != LIBRT_VECS_ABI_VERSION) {
        char err[128];
        snprintf(err, sizeof(err),
                 "ABI version conflict for librt.vecs, expected %d, found %d",
                 LIBRT_VECS_ABI_VERSION,
                 capsule->abi_version());
        PyErr_SetString(PyExc_ValueError, err);
        return -1;
    }
    if (capsule->api_version() < LIBRT_VECS_API_VERSION) {
        char err[128];
        snprintf(err, sizeof(err),
                 "API version conflict for librt.vecs, expected %d or newer, found %d (hint: upgrade librt)",
                 LIBRT_VECS_API_VERSION,
                 capsule->api_version());
        PyErr_SetString(PyExc_ValueError, err);
        return -1;
    }
    VecApi = capsule;
    VecI64Api = *VecApi->i64;
    VecI32Api = *VecApi->i32;
    VecI16Api = *VecApi->i16;
    VecU8Api = *VecApi->u8;
    VecFloatApi = *VecApi->float_;
    VecBoolApi = *VecApi->bool_;
    VecTApi = *VecApi->t;
    VecNestedApi = *VecApi->nested;
    return 0;
}
