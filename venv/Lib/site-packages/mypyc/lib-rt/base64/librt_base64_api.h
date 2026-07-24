#ifndef LIBRT_BASE64_API_H
#define LIBRT_BASE64_API_H

#include "librt_base64.h"

extern void *LibRTBase64_API[LIBRT_BASE64_API_LEN];

#define LibRTBase64_ABIVersion (*(int (*)(void)) LibRTBase64_API[0])
#define LibRTBase64_APIVersion (*(int (*)(void)) LibRTBase64_API[1])
#define LibRTBase64_b64encode_internal (*(PyObject* (*)(PyObject *source, bool urlsafe)) LibRTBase64_API[2])
#define LibRTBase64_b64decode_internal (*(PyObject* (*)(PyObject *source, bool urlsafe)) LibRTBase64_API[3])

int import_librt_base64(void);

#endif  // LIBRT_BASE64_API_H
