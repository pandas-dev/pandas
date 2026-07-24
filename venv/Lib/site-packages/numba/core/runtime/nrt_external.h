#ifndef NUMBA_NRT_EXTERNAL_H_
#define NUMBA_NRT_EXTERNAL_H_

#include <stdlib.h>

typedef struct MemInfo NRT_MemInfo;

typedef void NRT_managed_dtor(void *data);

typedef void *(*NRT_external_malloc_func)(size_t size, void *opaque_data);
typedef void *(*NRT_external_realloc_func)(void *ptr, size_t new_size, void *opaque_data);
typedef void (*NRT_external_free_func)(void *ptr, void *opaque_data);

struct ExternalMemAllocator {
    NRT_external_malloc_func malloc;
    NRT_external_realloc_func realloc;
    NRT_external_free_func free;
    void *opaque_data;
};

typedef struct ExternalMemAllocator NRT_ExternalAllocator;

typedef struct {
    /* Methods to create MemInfos.

    MemInfos are like smart pointers for objects that are managed by the Numba.
    */

    /* Allocate memory

    *nbytes* is the number of bytes to be allocated

    Returning a new reference.
    */
    NRT_MemInfo* (*allocate)(size_t nbytes);
    /* Allocates memory using an external allocator but still using Numba's MemInfo.
     *
     * NOTE: An externally provided allocator must behave the same way as C99
     *       stdlib.h's "malloc" function with respect to return value
     *       (including the behaviour that occurs when requesting an allocation
     *        of size 0 bytes).
     */
    NRT_MemInfo* (*allocate_external)(size_t nbytes, NRT_ExternalAllocator *allocator);

    /* Convert externally allocated memory into a MemInfo.

    *data* is the memory pointer
    *dtor* is the deallocator of the memory
    */
    NRT_MemInfo* (*manage_memory)(void *data, NRT_managed_dtor dtor);

    /* Acquire a reference */
    void (*acquire)(NRT_MemInfo* mi);

    /* Release a reference */
    void (*release)(NRT_MemInfo* mi);

    /* Get MemInfo data pointer */
    void* (*get_data)(NRT_MemInfo* mi);

} NRT_api_functions;



#endif /* NUMBA_NRT_EXTERNAL_H_ */
