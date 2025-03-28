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


/* The version for metalayer format; starts from 0 and it must not exceed 127 */
#define B2ND_METALAYER_VERSION 0

/* The maximum number of dimensions for b2nd arrays */
#define B2ND_MAX_DIM 8

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
 * @param index Indexes of the single-dimensional entries to remove.
 *
 * @return An error code
 */
BLOSC_EXPORT int b2nd_squeeze_index(b2nd_array_t *array, const bool *index);

/**
 * @brief Squeeze a b2nd array
 *
 * This function remove single-dimensional entries from the shape of a b2nd array.
 *
 * @param array The b2nd array.
 *
 * @return An error code
 */
BLOSC_EXPORT int b2nd_squeeze(b2nd_array_t *array);

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
 * @brief Read the metainfo in the b2nd metalayer.
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
 * @note This function is inlined and available even when not linking with libblosc2.
 *
 * @return An error code.
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
 * @note Please make sure that slice boundaries fit within the source and
 * destination arrays before using this function, as it does not perform these
 * checks itself.
 */
BLOSC_EXPORT int b2nd_copy_buffer(int8_t ndim,
                                  uint8_t itemsize,
                                  const void *src, const int64_t *src_pad_shape,
                                  const int64_t *src_start, const int64_t *src_stop,
                                  void *dst, const int64_t *dst_pad_shape,
                                  const int64_t *dst_start);


#ifdef __cplusplus
}
#endif

#endif  /* BLOSC_B2ND_H */
