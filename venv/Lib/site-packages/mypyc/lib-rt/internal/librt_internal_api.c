#include "librt_internal_api.h"

void *NativeInternal_API[LIBRT_INTERNAL_API_LEN] = {0};

int
import_librt_internal(void)
{
    PyObject *mod = PyImport_ImportModule("librt.internal");
    if (mod == NULL)
        return -1;
    Py_DECREF(mod);  // we import just for the side effect of making the below work.
    void **capsule = (void **)PyCapsule_Import("librt.internal._C_API", 0);
    if (capsule == NULL)
        return -1;

    // Each librt capsule gives at least 20 entries. Validate version before copying the table.
    int (*abi_version)(void) = (int (*)(void))capsule[13];
    if (abi_version() != LIBRT_INTERNAL_ABI_VERSION) {
        char err[128];
        snprintf(err, sizeof(err), "ABI version conflict for librt.internal, expected %d, found %d",
            LIBRT_INTERNAL_ABI_VERSION,
            abi_version()
        );
        PyErr_SetString(PyExc_ValueError, err);
        return -1;
    }
    int (*api_version)(void) = (int (*)(void))capsule[19];
    if (api_version() < LIBRT_INTERNAL_API_VERSION) {
        char err[128];
        snprintf(err, sizeof(err),
                 "API version conflict for librt.internal, expected %d or newer, found %d (hint: upgrade librt)",
            LIBRT_INTERNAL_API_VERSION,
            api_version()
        );
        PyErr_SetString(PyExc_ValueError, err);
        return -1;
    }
    // Provider API version is >= our expected version, which (by the API
    // compatibility contract) means it has at least LIBRT_INTERNAL_API_LEN
    // entries, so this copy is safe.
    memcpy(NativeInternal_API, capsule, sizeof(NativeInternal_API));
    return 0;
}
