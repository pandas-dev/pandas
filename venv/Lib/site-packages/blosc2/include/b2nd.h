/*********************************************************************
  Blosc - Blocked Shuffling and Compression Library

  Copyright (c) 2021  Blosc Development Team <blosc@blosc.org>
  https://blosc.org
  License: BSD 3-Clause (see LICENSE.txt)

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

/** @file b2nd.h
 * @brief Blosc2 NDim header file.
 *
 * This file contains Blosc2 NDim public API and the structures needed to use it.
 * @author Blosc Development Team <blosc@blosc.org>
 */

#ifndef BLOSC_B2ND_H
#define BLOSC_B2ND_H

#ifdef __cplusplus
extern "C" {
#endif
#include "blosc2/blosc2-export.h"
#ifdef __cplusplus
}
#endif

#include "blosc2.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_MSC_VER)
#define B2ND_DEPRECATED(msg) __declspec(deprecated(msg))
#elif defined(__GNUC__) || defined(__clang__)
#define B2ND_DEPRECATED(msg) __attribute__((deprecated(msg)))
#else
#define B2ND_DEPRECATED(msg)
#endif

/* The version for metalayer format; starts from 0 and it must not exceed 127 */
#define B2ND_METALAYER_VERSION 0

/* The maximum number of dimensions for b2nd arrays */
#define B2ND_MAX_DIM 16

/* The maximum number of metalayers for b2nd arrays */
#define B2ND_MAX_METALAYERS (BLOSC2_MAX_METALAYERS - 1)

/* NumPy dtype format
 * https://numpy.org/doc/stable/reference/arrays.dtypes.html#arrays-dtypes-constructing
 */
#define DTYPE_NUMPY_FORMAT 0

/* The default data type */
#define B2ND_DEFAULT_DTYPE "|u1"
/* The default data format */
#define B2ND_DEFAULT_DTYPE_FORMAT DTYPE_NUMPY_FORMAT

/**
 * @brief An *optional* cache for a single block.
 *
 * When a chunk is needed, it is copied into this cache. In this way, if the same chunk is needed
 * again afterwards, it is not necessary to recover it because it is already in the cache.
 */
struct chunk_cache_s {
  uint8_t *data;
  //!< The chunk data.
  int64_t nchunk;
  //!< The chunk number in cache. If @p nchunk equals to -1, it means that the cache is empty.
};

/**
 * @brief General parameters needed for the creation of a b2nd array.
 */
typedef struct b2nd_context_s b2nd_context_t;   /* opaque type */

/**
 * @brief A multidimensional array of data that can be compressed.
 */
typedef struct {
  blosc2_schunk *sc;
  //!< Pointer to a Blosc super-chunk
  int64_t shape[B2ND_MAX_DIM];
  //!< Shape of original data.
  int32_t chunkshape[B2ND_MAX_DIM];
  //!< Shape of each chunk.
  int64_t extshape[B2ND_MAX_DIM];
  //!< Shape of padded data.
  int32_t blockshape[B2ND_MAX_DIM];
  //!< Shape of each block.
  int64_t extchunkshape[B2ND_MAX_DIM];
  //!< Shape of padded chunk.
  int64_t nitems;
  //!< Number of items in original data.
  int32_t chunknitems;
  //!< Number of items in each chunk.
  int64_t extnitems;
  //!< Number of items in padded data.
  int32_t blocknitems;
  //!< Number of items in each block.
  int64_t extchunknitems;
  //!< Number of items in a padded chunk.
  int8_t ndim;
  //!< Data dimensions.
  struct chunk_cache_s chunk_cache;
  //!< A partition cache.
  int64_t item_array_strides[B2ND_MAX_DIM];
  //!< Item - shape strides.
  int64_t item_chunk_strides[B2ND_MAX_DIM];
  //!< Item - shape strides.
  int64_t item_extchunk_strides[B2ND_MAX_DIM];
  //!< Item - shape strides.
  int64_t item_block_strides[B2ND_MAX_DIM];
  //!< Item - shape strides.
  int64_t block_chunk_strides[B2ND_MAX_DIM];
  //!< Item - shape strides.
  int64_t chunk_array_strides[B2ND_MAX_DIM];
  //!< Item - shape strides.
  char *dtype;
  //!< Data type. Different formats can be supported (see dtype_format).
  int8_t dtype_format;
  //!< The format of the data type.  Default is DTYPE_NUMPY_FORMAT.
  int64_t last_tick;
  //!< The schunk change_tick this struct was last synced with (see blosc2_schunk).
} b2nd_array_t;


/**
 * @brief Create b2nd params.
 *
 * @param b2_storage The Blosc2 storage params.
 * @param ndim The dimensions.
 * @param shape The shape.
 * @param chunkshape The chunk shape.
 * @param blockshape The block shape.
 * @param dtype The data type expressed as a string version.
 * @param dtype_format The data type format; DTYPE_NUMPY_FORMAT should be chosen for NumPy compatibility.
 * @param metalayers The memory pointer to the list of the metalayers desired.
 * @param nmetalayers The number of metalayers.
 *
 * @return A pointer to the new b2nd params. NULL is returned if this fails.
 *
 * @note The pointer returned must be freed when not used anymore with #b2nd_free_ctx.
 *
 */
BLOSC_EXPORT b2nd_context_t *
b2nd_create_ctx(const blosc2_storage *b2_storage, int8_t ndim, const int64_t *shape, const int32_t *chunkshape,
                const int32_t *blockshape, const char *dtype, int8_t dtype_format, const blosc2_metalayer *metalayers,
                int32_t nmetalayers);


/**
 * @brief Free the resources associated with b2nd_context_t.
 *
 * @param ctx The b2nd context to free.
 *
 * @return An error code.
 *
 * @note This is safe in the sense that it will not free the schunk pointer in internal cparams.
 *
 */
BLOSC_EXPORT int b2nd_free_ctx(b2nd_context_t *ctx);


/**
 * @brief Create an uninitialized array.
 *
 * @param ctx The b2nd context for the new array.
 * @param array The memory pointer where the array will be created.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_uninit(b2nd_context_t *ctx, b2nd_array_t **array);


/**
 * @brief Create an empty array.
 *
 * @param ctx The b2nd context for the new array.
 * @param array The memory pointer where the array will be created.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_empty(b2nd_context_t *ctx, b2nd_array_t **array);


/**
 * Create an array, with zero being used as the default value for
 * uninitialized portions of the array.
 *
 * @param ctx The b2nd context for the new array.
 * @param array The memory pointer where the array will be created.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_zeros(b2nd_context_t *ctx, b2nd_array_t **array);


/**
 * Create an array, with NaN being used as the default value for
 * uninitialized portions of the array. Should only be used with type sizes
 * of either 4 or 8. Other sizes generate an error.
 *
 * @param ctx The b2nd context for the new array.
 * @param array The memory pointer where the array will be created.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_nans(b2nd_context_t *ctx, b2nd_array_t **array);


/**
 * Create an array, with @p fill_value being used as the default value for
 * uninitialized portions of the array.
 *
 * @param ctx The b2nd context for the new array.
 * @param array The memory pointer where the array will be created.
 * @param fill_value Default value for uninitialized portions of the array.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_full(b2nd_context_t *ctx, b2nd_array_t **array, const void *fill_value);

/**
 * @brief Free an array.
 *
 * @param array The array.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_free(b2nd_array_t *array);

/**
 * @brief Create a b2nd array from a super-chunk. It can only be used if the array
 * is backed by a blosc super-chunk.
 *
 * @param schunk The blosc super-chunk where the b2nd array is stored.
 * @param array The memory pointer where the array will be created.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_from_schunk(blosc2_schunk *schunk, b2nd_array_t **array);

/**
 * Create a serialized super-chunk from a b2nd array.
 *
 * @param array The b2nd array to be serialized.
 * @param cframe The pointer of the buffer where the in-memory array will be copied.
 * @param cframe_len The length of the in-memory array buffer.
 * @param needs_free Whether the buffer should be freed or not.
 *
 * @return An error code
 */
BLOSC_EXPORT int b2nd_to_cframe(const b2nd_array_t *array, uint8_t **cframe,
                                int64_t *cframe_len, bool *needs_free);

/**
 * @brief Create a b2nd array from a serialized super-chunk.
 *
 * @param cframe The buffer of the in-memory array.
 * @param cframe_len The size (in bytes) of the in-memory array.
 * @param copy Whether b2nd should make a copy of the cframe data or not. The copy will be made to an internal sparse frame.
 * @param array The memory pointer where the array will be created.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_from_cframe(uint8_t *cframe, int64_t cframe_len, bool copy, b2nd_array_t **array);

/**
 * @brief Open a b2nd array from a file.
 *
 * @param urlpath The path of the b2nd array on disk.
 * @param array The memory pointer where the array info will be stored.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_open(const char *urlpath, b2nd_array_t **array);

/**
 * @brief Open a b2nd array from a file using an offset.
 *
 * @param urlpath The path of the b2nd array on disk.
 * @param array The memory pointer where the array info will be stored.
 * @param offset The offset in the file where the b2nd array frame starts.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_open_offset(const char *urlpath, b2nd_array_t **array, int64_t offset);

/**
 * @brief Save b2nd array into a specific urlpath.
 *
 * @param array The array to be saved.
 * @param urlpath The urlpath where the array will be stored.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_save(const b2nd_array_t *array, char *urlpath);

/**
 * @brief Append a b2nd array into a file.
 *
 * @param array The array to write.
 * @param urlpath The path for persistent storage.
 *
 * @return If successful, return the offset where @p array has been appended in @p urlpath.
 * Else, a negative value.
 */
BLOSC_EXPORT int64_t b2nd_save_append(const b2nd_array_t *array, const char *urlpath);

/**
 * @brief Create a b2nd array from a C buffer.
 *
 * @param ctx The b2nd context for the new array.
 * @param array The memory pointer where the array will be created.
 * @param buffer The buffer where source data is stored.
 * @param buffersize The size (in bytes) of the buffer.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_from_cbuffer(b2nd_context_t *ctx, b2nd_array_t **array, const void *buffer, int64_t buffersize);

/**
 * @brief Extract the data from a b2nd array into a C buffer.
 *
 * @param array The b2nd array.
 * @param buffer The buffer where the data will be stored.
 * @param buffersize Size (in bytes) of the buffer.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_to_cbuffer(const b2nd_array_t *array, void *buffer, int64_t buffersize);

/**
 * @brief Gather items from a b2nd array of any dimensionality by flat logical coordinates.
 *
 * Fully supports arrays of any number of dimensions (1-D, 2-D, ... N-D).
 * Regardless of the array's dimensionality, every element is addressed by a
 * single flat, C-order logical index in the range ``[0, array->nitems)``
 * -- the same linearisation used by ``b2nd_from_cbuffer`` / ``b2nd_to_cbuffer``.
 * For example, element ``(i, j)`` of a 2-D array with shape ``(R, C)`` has
 * flat index ``i * C + j``; element ``(i, j, k)`` of a 3-D array with shape
 * ``(A, B, C)`` has flat index ``i * B * C + j * C + k``; and so on for
 * higher dimensions.
 *
 * The function translates each flat logical coordinate to the padded storage
 * coordinate used internally by b2nd (chunk → block → item within block) and
 * then gathers the requested items into ``buffer`` in the same order as
 * ``coords``.  Repeated and unsorted coordinates are allowed.
 *
 * This is a low-level sparse gather primitive.  It is most useful for fancy
 * indexing and as a building block for higher-level ``take`` operations.
 *
 * @param array The source b2nd array (any ndim).
 * @param ncoords Number of coordinates to gather.  May be zero.
 * @param coords Flat C-order logical coordinates in ``[0, array->nitems)``.
 *               Must not be NULL when ``ncoords > 0``.
 * @param buffer Destination buffer.  Must have room for at least
 *               ``ncoords * array->sc->typesize`` bytes and must not be NULL
 *               when ``ncoords > 0``.
 * @param buffersize Destination buffer size in bytes.
 *
 * @return A non-negative value on success, a negative error code otherwise.
 */
BLOSC_EXPORT int b2nd_get_sparse_cbuffer(const b2nd_array_t *array, int64_t ncoords,
                                         const int64_t *coords, void *buffer, int64_t buffersize);


/**
 * @brief Get a slice from an array and store it into a new array.
 *
 * @param ctx The b2nd context for the new array.
 * @param array The memory pointer where the array will be created.
 * @param src The array from which the slice will be extracted
 * @param start The coordinates where the slice will begin.
 * @param stop The coordinates where the slice will end.
 *
 * @return An error code.
 *
 * @note The ndim and shape from ctx will be overwritten by the src and stop-start respectively.
 *
 */
BLOSC_EXPORT int b2nd_get_slice(b2nd_context_t *ctx, b2nd_array_t **array, const b2nd_array_t *src,
                                const int64_t *start, const int64_t *stop);

/**
 * @brief Squeeze a b2nd array
 *
 * This function remove selected single-dimensional entries from the shape of a
 b2nd array.
 *
 * @param array The b2nd array.
 * @param view The memory pointer where the new view will be created.
 * @param index Indexes of the single-dimensional entries to remove.
 *
 * @return An error code
 */
BLOSC_EXPORT int b2nd_squeeze_index(b2nd_array_t *array, b2nd_array_t **view, const bool *index);

/**
 * @brief Squeeze a b2nd array
 *
 * This function remove single-dimensional entries from the shape of a b2nd array.
 *
 * @param array The b2nd array.
 *  @param view The memory pointer where the new view will be created.
 *
 * @return An error code
 */
BLOSC_EXPORT int b2nd_squeeze(b2nd_array_t *array, b2nd_array_t **view);

/**
 * @brief Add a newaxis to a b2nd array at location @p axis.
 *
 * @param array The b2nd array to be expanded.
 * @param axis The axes where the new dimensions will be added.
 * @param view The memory pointer where the new view will be created.
 * @param final_dims The final number of dimensions. Should be same as the number of elements in @p axis.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_expand_dims(const b2nd_array_t *array, b2nd_array_t **view, const bool *axis,
    const uint8_t final_dims);

/**
 * @brief Get a slice from an array and store it into a C buffer.
 *
 * @param array The array from which the slice will be extracted.
 * @param start The coordinates where the slice will begin.
 * @param stop The coordinates where the slice will end.
 * @param buffershape The shape of the buffer.
 * @param buffer The buffer where the data will be stored.
 * @param buffersize The size (in bytes) of the buffer.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_get_slice_cbuffer(const b2nd_array_t *array, const int64_t *start, const int64_t *stop,
                                        void *buffer, const int64_t *buffershape, int64_t buffersize);

/**
 * @brief Set a slice in a b2nd array using a C buffer.
 *
 * @param buffer The buffer where the slice data is.
 * @param buffershape The shape of the buffer.
 * @param buffersize The size (in bytes) of the buffer.
 * @param start The coordinates where the slice will begin.
 * @param stop The coordinates where the slice will end.
 * @param array The b2nd array where the slice will be set
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_set_slice_cbuffer(const void *buffer, const int64_t *buffershape, int64_t buffersize,
                                        const int64_t *start, const int64_t *stop, b2nd_array_t *array);

/**
 * @brief Make a copy of the array data. The copy is done into a new b2nd array.
 *
 * @param ctx The b2nd context for the new array.
 * @param src The array from which data is copied.
 * @param array The memory pointer where the array will be created.
 *
 * @return An error code
 *
 * @note The ndim and shape in ctx will be overwritten by the src ctx.
 *
 */
BLOSC_EXPORT int b2nd_copy(b2nd_context_t *ctx, const b2nd_array_t *src, b2nd_array_t **array);

/**
 * @brief Concatenate arrays. The result is stored in a new b2nd array, or an enlarged one.
 *
 * @param ctx The b2nd context for the new array.
 * @param src1 The first array from which data is copied.
 * @param src2 The second array from which data is copied.
 * @param axis The axis along which the arrays will be concatenated.
 * @param copy Whether the data should be copied or not. If false, the @p src1 array
 *   will be expanded as needed to keep the result.
 * @param array The memory pointer where the array will be created.  It will have the same
 *   metalayers of @p src1, except for the b2nd metalayer, which will be updated with the
 *   new shape.
 *
 * @ note The two arrays must have the same shape in all dimensions except the concatenation axis.
 * Also, the typesize of the two arrays must be the same.
 *
 * @return An error code
 *
 * @note The ndim and shape in ctx will be overwritten by the src1 ctx.
 *
 */
BLOSC_EXPORT int b2nd_concatenate(b2nd_context_t *ctx, const b2nd_array_t *src1, const b2nd_array_t *src2,
                                  int8_t axis, bool copy, b2nd_array_t **array);

/**
 * @brief Print metalayer parameters.
 *
 * @param array The array where the metalayer is stored.
 *
 * @return An error code
 */
BLOSC_EXPORT int b2nd_print_meta(const b2nd_array_t *array);

/**
 * @brief Resize the shape of an array
 *
 * @param array The array to be resized.
 * @param new_shape The new shape from the array.
 * @param start The position in which the array will be extended or shrunk.
 *
 * @return An error code
 */
BLOSC_EXPORT int b2nd_resize(b2nd_array_t *array, const int64_t *new_shape, const int64_t *start);


/**
 * @brief Re-sync the cached geometry (shape, strides, metalayers) of a
 * disk-based array when another handle — possibly in another process — has
 * changed it (e.g. via @ref b2nd_resize).
 *
 * This gives a deterministic sync point for readers in a single-writer,
 * multiple-readers (SWMR) workflow: after calling it, the shape in @p array
 * reflects the current on-disk state.  Data-access functions (get/set slice
 * and friends) already do this implicitly, so calling it is only needed to
 * observe shape changes without reading data (e.g. when polling).
 *
 * A no-op for arrays not backed by an on-disk frame.  Only shape changes are
 * followed: if ndim, chunkshape or blockshape changed behind this handle,
 * BLOSC2_ERROR_DATA is returned.
 *
 * @param array The array to refresh.
 *
 * @return 1 if the geometry was re-synced, 0 if it was already current, or a
 * negative error code.
 */
BLOSC_EXPORT int b2nd_refresh(b2nd_array_t *array);


/**
 * @brief Insert given buffer in an array extending the given axis.
 *
 * @param array The array to insert the data in.
 * @param buffer The buffer data to be inserted.
 * @param buffersize The size (in bytes) of the buffer.
 * @param axis The axis that will be extended.
 * @param insert_start The position inside the axis to start inserting the data.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_insert(b2nd_array_t *array, const void *buffer, int64_t buffersize,
                             int8_t axis, int64_t insert_start);

/**
 * Append a buffer at the end of a b2nd array.
 *
 * @param array The array to append the data in.
 * @param buffer The buffer data to be appended.
 * @param buffersize Size (in bytes) of the buffer.
 * @param axis The axis that will be extended to append the data.
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_append(b2nd_array_t *array, const void *buffer, int64_t buffersize,
                             int8_t axis);

/**
 * @brief Delete shrinking the given axis delete_len items.
 *
 * @param array The array to shrink.
 * @param axis The axis to shrink.
 * @param delete_start The start position from the axis to start deleting chunks.
 * @param delete_len The number of items to delete to the array->shape[axis].
 *   The newshape[axis] will be the old array->shape[axis] - delete_len
 *
 * @return An error code.
 *
 * @note See also b2nd_resize
 */
BLOSC_EXPORT int b2nd_delete(b2nd_array_t *array, int8_t axis,
                             int64_t delete_start, int64_t delete_len);


// Indexing section

/**
 * @brief Get an element selection along each dimension of an array independently.
 *
 * @param array The array to get the data from.
 * @param selection The elements along each dimension.
 * @param selection_size The size of the selection along each dimension.
 * @param buffer The buffer for getting the data.
 * @param buffershape The shape of the buffer.
 * @param buffersize The buffer size (in bytes).
 *
 * @return An error code.
 *
 * @note See also b2nd_set_orthogonal_selection.
 */
BLOSC_EXPORT int b2nd_get_orthogonal_selection(const b2nd_array_t *array, int64_t **selection,
                                               int64_t *selection_size, void *buffer,
                                               int64_t *buffershape, int64_t buffersize);

/**
 * @brief Set an element selection along each dimension of an array independently.
 *
 * @param array The array to set the data to.
 * @param selection The elements along each dimension.
 * @param selection_size The size of the selection along each dimension.
 * @param buffer The buffer with the data for setting.
 * @param buffershape The shape of the buffer.
 * @param buffersize The buffer size (in bytes).
 *
 * @return An error code.
 *
 * @note See also b2nd_get_orthogonal_selection.
 */
BLOSC_EXPORT int b2nd_set_orthogonal_selection(b2nd_array_t *array, int64_t **selection,
                                               int64_t *selection_size, const void *buffer,
                                               int64_t *buffershape, int64_t buffersize);


/**
 * @brief Create the metainfo for the b2nd metalayer.
 *
 * @param ndim The number of dimensions in the array.
 * @param shape The shape of the array.
 * @param chunkshape The shape of the chunks in the array.
 * @param blockshape The shape of the blocks in the array.
 * @param dtype A string representation of the data type of the array.
 * @param dtype_format The format of the dtype representation. 0 means NumPy.
 * @param smeta The msgpack buffer (output).
 *
 * @return An error code.
 */
BLOSC_EXPORT int b2nd_serialize_meta(int8_t ndim, const int64_t *shape, const int32_t *chunkshape,
                                     const int32_t *blockshape, const char *dtype,
                                     int8_t dtype_format, uint8_t **smeta);

/**
 * @brief Read the metainfo in the b2nd metalayer (inline version).
 *
 * @param smeta The msgpack buffer (input).
 * @param smeta_len The length of the smeta buffer (input).
 * @param ndim The number of dimensions in the array (output).
 * @param shape The shape of the array (output).
 * @param chunkshape The shape of the chunks in the array (output).
 * @param blockshape The shape of the blocks in the array (output).
 * @param dtype A string representation of the data type of the array (output).
 * @param dtype_format The format of the dtype representation (output). 0 means NumPy (the default).
 *
 * @note This function is inlined so that external codec/filter plugins (like
 * blosc2_grok) can use it without linking against libblosc2.  This avoids
 * pulling all of libblosc2's symbols (e.g. internal ZFP, Zstd) into the
 * global namespace at load time, which would otherwise shadow symbols from
 * other libraries that need differently-configured builds of the same
 * dependencies.
 *
 * @return An error code.
 */
static inline int b2nd_deserialize_meta_inline(const uint8_t *smeta, int32_t smeta_len, int8_t *ndim, int64_t *shape,
                                                int32_t *chunkshape, int32_t *blockshape, char **dtype, int8_t *dtype_format) {
  BLOSC_ERROR_NULL(smeta, BLOSC2_ERROR_NULL_POINTER);
  BLOSC_ERROR_NULL(ndim, BLOSC2_ERROR_NULL_POINTER);
  BLOSC_ERROR_NULL(shape, BLOSC2_ERROR_NULL_POINTER);
  BLOSC_ERROR_NULL(chunkshape, BLOSC2_ERROR_NULL_POINTER);
  BLOSC_ERROR_NULL(blockshape, BLOSC2_ERROR_NULL_POINTER);
  if (dtype != NULL) {
    *dtype = NULL;
  }
  if (dtype_format != NULL) {
    *dtype_format = 0;
  }
  if (smeta_len <= 0) {
    BLOSC_TRACE_ERROR("Malformed b2nd metalayer: empty metadata");
    return BLOSC2_ERROR_FAILURE;
  }

  const uint8_t *pmeta = smeta;

#define B2ND_REQUIRE_META_NBYTES(nbytes)                                             \
  do {                                                                                \
    size_t consumed = (size_t)(pmeta - smeta);                                        \
    size_t total = (size_t)smeta_len;                                                 \
    if (consumed > total || (total - consumed) < (size_t)(nbytes)) {                 \
      BLOSC_TRACE_ERROR("Malformed b2nd metalayer: truncated metadata");            \
      return BLOSC2_ERROR_FAILURE;                                                    \
    }                                                                                 \
  } while (0)

  // Check that we have an array with 7 entries (version, ndim, shape, chunkshape, blockshape, dtype_format, dtype)
  B2ND_REQUIRE_META_NBYTES(1);
  pmeta += 1;

  // version entry
  // int8_t version = (int8_t)pmeta[0];  // positive fixnum (7-bit positive integer) commented to avoid warning
  B2ND_REQUIRE_META_NBYTES(1);
  pmeta += 1;

  // ndim entry
  B2ND_REQUIRE_META_NBYTES(1);
  *ndim = (int8_t) pmeta[0];
  int8_t ndim_aux = *ndim;  // positive fixnum (7-bit positive integer)
  if (ndim_aux < 0 || ndim_aux > B2ND_MAX_DIM) {
    BLOSC_TRACE_ERROR("ndim %d is out of range", ndim_aux);
    return BLOSC2_ERROR_FAILURE;
  }
  pmeta += 1;

  // shape entry
  // Initialize to ones, as required by b2nd
  for (int i = 0; i < ndim_aux; i++) shape[i] = 1;
  B2ND_REQUIRE_META_NBYTES(1);
  pmeta += 1;
  for (int8_t i = 0; i < ndim_aux; i++) {
    B2ND_REQUIRE_META_NBYTES(1 + sizeof(int64_t));
    pmeta += 1;
    swap_store(shape + i, pmeta, sizeof(int64_t));
    pmeta += sizeof(int64_t);
  }

  // chunkshape entry
  // Initialize to ones, as required by b2nd
  for (int i = 0; i < ndim_aux; i++) chunkshape[i] = 1;
  B2ND_REQUIRE_META_NBYTES(1);
  pmeta += 1;
  for (int8_t i = 0; i < ndim_aux; i++) {
    B2ND_REQUIRE_META_NBYTES(1 + sizeof(int32_t));
    pmeta += 1;
    swap_store(chunkshape + i, pmeta, sizeof(int32_t));
    pmeta += sizeof(int32_t);
  }

  // blockshape entry
  // Initialize to ones, as required by b2nd
  for (int i = 0; i < ndim_aux; i++) blockshape[i] = 1;
  B2ND_REQUIRE_META_NBYTES(1);
  pmeta += 1;
  for (int8_t i = 0; i < ndim_aux; i++) {
    B2ND_REQUIRE_META_NBYTES(1 + sizeof(int32_t));
    pmeta += 1;
    swap_store(blockshape + i, pmeta, sizeof(int32_t));
    pmeta += sizeof(int32_t);
  }

  // dtype entry
  if (dtype_format == NULL || dtype == NULL) {
    return (int32_t)(pmeta - smeta);
  }
  if (pmeta - smeta < smeta_len) {
    // dtype info is here
    B2ND_REQUIRE_META_NBYTES(1 + 1 + sizeof(int32_t));
    *dtype_format = (int8_t) *(pmeta++);
    if (*pmeta != 0xdb) {
      BLOSC_TRACE_ERROR("Malformed b2nd metalayer: invalid dtype MsgPack marker");
      return BLOSC2_ERROR_FAILURE;
    }
    pmeta += 1;
    int32_t dtype_len;
    swap_store(&dtype_len, pmeta, sizeof(int32_t));
    pmeta += sizeof(int32_t);
    if (dtype_len < 0) {
      BLOSC_TRACE_ERROR("Malformed b2nd metalayer: negative dtype length");
      return BLOSC2_ERROR_FAILURE;
    }
    B2ND_REQUIRE_META_NBYTES(dtype_len);
    size_t dtype_len_ = (size_t)dtype_len;
    *dtype = (char*)malloc(dtype_len_ + 1);
    BLOSC_ERROR_NULL(*dtype, BLOSC2_ERROR_MEMORY_ALLOC);
    char* dtype_ = *dtype;
    memcpy(dtype_, (char*)pmeta, dtype_len_);
    dtype_[dtype_len_] = '\0';
    pmeta += dtype_len_;
  }
  else {
    // dtype is mandatory in b2nd metalayer, but this is mainly meant as
    // a fall-back for deprecated caterva headers
    *dtype = NULL;
    *dtype_format = 0;
  }

#undef B2ND_REQUIRE_META_NBYTES

  int32_t slen = (int32_t) (pmeta - smeta);
  return (int)slen;
}

/**
 * @brief Read the metainfo in the b2nd metalayer (ABI entry point).
 *
 * This is a wrapper around @ref b2nd_deserialize_meta_inline for ABI compatibility.
 * New code that only includes headers should use b2nd_deserialize_meta_inline()
 * directly to avoid the need to link against libblosc2.
 *
 * @see b2nd_deserialize_meta_inline
 */
BLOSC_EXPORT int b2nd_deserialize_meta(const uint8_t *smeta, int32_t smeta_len, int8_t *ndim, int64_t *shape,
                                       int32_t *chunkshape, int32_t *blockshape, char **dtype, int8_t *dtype_format);

// Utilities for C buffers representing multidimensional arrays

/**
 * @brief Copy a slice of a source array into another array. The arrays have
 * the same number of dimensions (though their shapes may differ), the same
 * item size, and they are stored as C buffers with contiguous data (any
 * padding is considered part of the array).
 *
 * @param ndim The number of dimensions in both arrays.
 * @param itemsize The size of the individual data item in both arrays.
 * @param src The buffer for getting the data from the source array.
 * @param src_pad_shape The shape of the source array, including padding.
 * @param src_start The source coordinates where the slice will begin.
 * @param src_stop The source coordinates where the slice will end.
 * @param dst The buffer for setting the data into the destination array.
 * @param dst_pad_shape The shape of the destination array, including padding.
 * @param dst_start The destination coordinates where the slice will be placed.
 *
 * @return An error code.
 *
 * @note This is kept for backward compatibility with existing code out there.  New code should use
 * b2nd_copy_buffer2 instead.
 *
 * @note Please make sure that slice boundaries fit within the source and
 * destination arrays before using this function, as it does not perform these
 * checks itself.
 */
B2ND_DEPRECATED("Use b2nd_copy_buffer2 instead.")
BLOSC_EXPORT int b2nd_copy_buffer(int8_t ndim,
                                  uint8_t itemsize,
                                  const void *src, const int64_t *src_pad_shape,
                                  const int64_t *src_start, const int64_t *src_stop,
                                  void *dst, const int64_t *dst_pad_shape,
                                  const int64_t *dst_start);

/**
 * @brief Copy a slice of a source array into another array. The arrays have
 * the same number of dimensions (though their shapes may differ), the same
 * item size, and they are stored as C buffers with contiguous data (any
 * padding is considered part of the array).
 *
 * @param ndim The number of dimensions in both arrays.
 * @param itemsize The size of the individual data item in both arrays.
 * @param src The buffer for getting the data from the source array.
 * @param src_pad_shape The shape of the source array, including padding.
 * @param src_start The source coordinates where the slice will begin.
 * @param src_stop The source coordinates where the slice will end.
 * @param dst The buffer for setting the data into the destination array.
 * @param dst_pad_shape The shape of the destination array, including padding.
 * @param dst_start The destination coordinates where the slice will be placed.
 *
 * @return An error code.
 *
 * @note This is a version of (now deprecated) b2nd_copy_buffer() that uses
 * signed 32-bit integers for copying data. This is useful when data is stored
 * in a buffer that uses itemsizes that are larger than 255 bytes.
 *
 * @note Please make sure that slice boundaries fit within the source and
 * destination arrays before using this function, as it does not perform these
 * checks itself.
 */
BLOSC_EXPORT int b2nd_copy_buffer2(int8_t ndim,
                                   int32_t itemsize,
                                   const void *src, const int64_t *src_pad_shape,
                                   const int64_t *src_start, const int64_t *src_stop,
                                   void *dst, const int64_t *dst_pad_shape,
                                   const int64_t *dst_start);


#ifdef __cplusplus
}
#endif

#endif  /* BLOSC_B2ND_H */
