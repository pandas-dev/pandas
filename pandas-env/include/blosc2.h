/*********************************************************************
  Blosc - Blocked Shuffling and Compression Library

  Copyright (c) 2021  Blosc Development Team <blosc@blosc.org>
  https://blosc.org
  License: BSD 3-Clause (see LICENSE.txt)

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

/*********************************************************************
  @file blosc2.h
  @brief Blosc2 header file.

  This file contains Blosc2 public API and the structures needed to use it.
  @author Blosc Development Team <blosc@blosc.org>
**********************************************************************/

#ifndef BLOSC_BLOSC2_H
#define BLOSC_BLOSC2_H

#include "blosc2/blosc2-export.h"
#include "blosc2/blosc2-common.h"
#include "blosc2/blosc2-stdio.h"

#if defined(_WIN32) && !defined(__MINGW32__)
#include <windows.h>
#include <malloc.h>
#include <process.h>
#define getpid _getpid
#endif

#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// For compatibility with the Blosc 1.x series
#ifdef BLOSC1_COMPAT
  // Blosc2 symbols that should be accessible from Blosc 1.x API
  #define BLOSC_VERSION_MAJOR BLOSC2_VERSION_MAJOR
  #define BLOSC_VERSION_MINOR BLOSC2_VERSION_MINOR
  #define BLOSC_VERSION_RELEASE BLOSC2_VERSION_RELEASE
  #define BLOSC_VERSION_STRING BLOSC2_VERSION_STRING
  #define BLOSC_VERSION_DATE BLOSC2_VERSION_DATE
  #define BLOSC_MAX_OVERHEAD BLOSC2_MAX_OVERHEAD
  #define BLOSC_MAX_BUFFERSIZE BLOSC2_MAX_BUFFERSIZE

  // API that changed to blosc1_ prefix
  #define blosc_compress blosc1_compress
  #define blosc_decompress blosc1_decompress
  #define blosc_getitem blosc1_getitem
  #define blosc_get_compressor blosc1_get_compressor
  #define blosc_set_compressor blosc1_set_compressor
  #define blosc_cbuffer_sizes blosc1_cbuffer_sizes
  #define blosc_cbuffer_validate blosc1_cbuffer_validate
  #define blosc_cbuffer_metainfo blosc1_cbuffer_metainfo
  #define blosc_get_blocksize blosc1_get_blocksize
  #define blosc_set_blocksize blosc1_set_blocksize
  #define blosc_set_splitmode blosc1_set_splitmode

  // API that changed to blosc2_ prefix
  #define blosc_init blosc2_init
  #define blosc_destroy blosc2_destroy
  #define blosc_free_resources blosc2_free_resources
  #define blosc_get_nthreads blosc2_get_nthreads
  #define blosc_set_nthreads blosc2_set_nthreads
  #define blosc_compcode_to_compname blosc2_compcode_to_compname
  #define blosc_compname_to_compcode blosc2_compname_to_compcode
  #define blosc_list_compressors blosc2_list_compressors
  #define blosc_get_version_string blosc2_get_version_string
  #define blosc_get_complib_info blosc2_get_complib_info
  #define blosc_cbuffer_versions blosc2_cbuffer_versions
  #define blosc_cbuffer_complib blosc2_cbuffer_complib
#endif


/* Version numbers */
#define BLOSC2_VERSION_MAJOR    2    /* for major interface/format changes  */
#define BLOSC2_VERSION_MINOR    15   /* for minor interface/format changes  */
#define BLOSC2_VERSION_RELEASE  2    /* for tweaks, bug-fixes, or development */

#define BLOSC2_VERSION_STRING   "2.15.2.dev"  /* string version.  Sync with above! */
#define BLOSC2_VERSION_DATE     "$Date:: 2024-07-30 #$"    /* date version year-month-day */


/* The maximum number of dimensions for Blosc2 NDim arrays */
#define BLOSC2_MAX_DIM 8


/* Tracing macros */
#define BLOSC_TRACE_ERROR(msg, ...) BLOSC_TRACE(error, msg, ##__VA_ARGS__)
#define BLOSC_TRACE_WARNING(msg, ...) BLOSC_TRACE(warning, msg, ##__VA_ARGS__)
#define BLOSC_TRACE_INFO(msg, ...) BLOSC_TRACE(info, msg, ##__VA_ARGS__)
#define BLOSC_TRACE(cat, msg, ...)                  \
    do {                                            \
        const char *__e = getenv("BLOSC_TRACE");    \
        if (!__e) { break; }                        \
        fprintf(stderr, "[%s] - " msg " (%s:%d)\n", #cat, ##__VA_ARGS__, __FILE__, __LINE__); \
    } while(0)

#define BLOSC_ERROR_NULL(pointer, rc)               \
    do {                                            \
        if ((pointer) == NULL) {                    \
            BLOSC_TRACE_ERROR("Pointer is null");   \
            return (rc);                            \
        }                                           \
    } while (0)
#define BLOSC_ERROR(rc)                             \
    do {                                            \
        int rc_ = (rc);                             \
        if (rc_ < BLOSC2_ERROR_SUCCESS) {           \
            char *error_msg = print_error(rc_);     \
            BLOSC_TRACE_ERROR("%s", error_msg);     \
            return rc_;                             \
        }                                           \
    } while (0)

#define BLOSC_INFO(msg, ...)                        \
    do {                                            \
        const char *__e = getenv("BLOSC_INFO");     \
        if (!__e) { break; }                        \
        fprintf(stderr, "[INFO] - " msg "\n", ##__VA_ARGS__); \
    } while(0)


/* The VERSION_FORMAT symbols below should be just 1-byte long */
enum {
  /* Blosc format version, starting at 1
     1 -> Blosc pre-1.0
     2 -> Blosc 1.x stable series
     3 -> Blosc 2-alpha.x series
     4 -> Blosc 2.x beta.1 series
     5 -> Blosc 2.x stable series
     */
  BLOSC1_VERSION_FORMAT_PRE1 = 1,
  BLOSC1_VERSION_FORMAT = 2,
  BLOSC2_VERSION_FORMAT_ALPHA = 3,
  BLOSC2_VERSION_FORMAT_BETA1 = 4,
  BLOSC2_VERSION_FORMAT_STABLE = 5,
  BLOSC2_VERSION_FORMAT = BLOSC2_VERSION_FORMAT_STABLE,
};


/* The FRAME_FORMAT_VERSION symbols below should be just 4-bit long */
enum {
  /* Blosc format version
   *  1 -> First version (introduced in beta.2)
   *  2 -> Second version (introduced in rc.1)
   *
   */
  BLOSC2_VERSION_FRAME_FORMAT_BETA2 = 1,  // for 2.0.0-beta2 and after
  BLOSC2_VERSION_FRAME_FORMAT_RC1 = 2,    // for 2.0.0-rc1 and after
  BLOSC2_VERSION_FRAME_FORMAT = BLOSC2_VERSION_FRAME_FORMAT_RC1,
};


//!< Struct for storing data from instrumentation of codecs
// This can be flexible because it is typically used mainly for development
typedef struct {
  float cratio;
  float cspeed;
  float filter_speed;
  //float memory;
  //float power;
  uint8_t flags[4];
} blosc2_instr;


enum {
#ifndef BLOSC_H
  BLOSC_MIN_HEADER_LENGTH = 16,
  //!< Minimum header length (Blosc1)
#endif // BLOSC_H
  BLOSC_EXTENDED_HEADER_LENGTH = 32,
  //!< Extended header length (Blosc2, see README_HEADER)
  BLOSC2_MAX_OVERHEAD = BLOSC_EXTENDED_HEADER_LENGTH,
  //!< The maximum overhead during compression in bytes. This equals
  //!< to @ref BLOSC_EXTENDED_HEADER_LENGTH now, but can be higher in future
  //!< implementations.
  BLOSC2_MAX_BUFFERSIZE = (INT_MAX - BLOSC2_MAX_OVERHEAD),
  //!< Maximum source buffer size to be compressed
#ifndef BLOSC_H
  BLOSC_MAX_TYPESIZE = UINT8_MAX,
  //!< Maximum typesize before considering source buffer as a stream of bytes.
  //!< Cannot be larger than 255.
#endif // BLOSC_H
  BLOSC_MIN_BUFFERSIZE = 32,
  //!< Minimum buffer size to be compressed.
};

enum {
    BLOSC2_DEFINED_TUNER_START = 0,
    BLOSC2_DEFINED_TUNER_STOP = 31,
    //!< Blosc-defined tuners must be between 0 - 31.
    BLOSC2_GLOBAL_REGISTERED_TUNER_START = 32,
    BLOSC2_GLOBAL_REGISTERED_TUNER_STOP = 159,
    //!< Blosc-registered tuners must be between 31 - 159.
    BLOSC2_GLOBAL_REGISTERED_TUNERS = 0,
    //!< Number of Blosc-registered tuners at the moment.
    BLOSC2_USER_REGISTERED_TUNER_START = 160,
    BLOSC2_USER_REGISTERED_TUNER_STOP = 255,
    //!< User-defined tuners must be between 160 - 255.
};

/**
 * @brief Codes for the different tuners shipped with Blosc
 */
enum {
    BLOSC_STUNE = 0,
    BLOSC_LAST_TUNER = 1,
    //!< Determine the last tuner defined by Blosc.
    BLOSC_LAST_REGISTERED_TUNE = BLOSC2_GLOBAL_REGISTERED_TUNER_START + BLOSC2_GLOBAL_REGISTERED_TUNERS - 1,
    //!< Determine the last registered tuner. It is used to check if a tuner is registered or not.
};

enum {
  BLOSC2_DEFINED_FILTERS_START = 0,
  BLOSC2_DEFINED_FILTERS_STOP = 31,
  //!< Blosc-defined filters must be between 0 - 31.
  BLOSC2_GLOBAL_REGISTERED_FILTERS_START = 32,
  BLOSC2_GLOBAL_REGISTERED_FILTERS_STOP = 159,
  //!< Blosc-registered filters must be between 32 - 159.
  BLOSC2_GLOBAL_REGISTERED_FILTERS = 5,
  //!< Number of Blosc-registered filters at the moment.
  BLOSC2_USER_REGISTERED_FILTERS_START = 160,
  BLOSC2_USER_REGISTERED_FILTERS_STOP = 255,
  //!< User-defined filters must be between 128 - 255.
  BLOSC2_MAX_FILTERS = 6,
  //!< Maximum number of filters in the filter pipeline
  BLOSC2_MAX_UDFILTERS = 16,
  //!< Maximum number of filters that a user can register.
};


/**
 * @brief Codes for filters.
 *
 * @sa #blosc1_compress
 */
enum {
#ifndef BLOSC_H
  BLOSC_NOSHUFFLE = 0,
  //!< No shuffle (for compatibility with Blosc1).
  BLOSC_NOFILTER = 0,
  //!< No filter.
  BLOSC_SHUFFLE = 1,
  //!< Byte-wise shuffle. `filters_meta` does not have any effect here.
  BLOSC_BITSHUFFLE = 2,
  //!< Bit-wise shuffle. `filters_meta` does not have any effect here.
#endif // BLOSC_H
  BLOSC_DELTA = 3,
  //!< Delta filter. `filters_meta` does not have any effect here.
  BLOSC_TRUNC_PREC = 4,
  //!< Truncate mantissa precision.
  //!< Positive values in `filters_meta` will keep bits; negative values will zero bits.
  BLOSC_LAST_FILTER = 5,
  //!< sentinel
  BLOSC_LAST_REGISTERED_FILTER = BLOSC2_GLOBAL_REGISTERED_FILTERS_START + BLOSC2_GLOBAL_REGISTERED_FILTERS - 1,
  //!< Determine the last registered filter. It is used to check if a filter is registered or not.
};

/**
 * @brief Codes for internal flags (see blosc1_cbuffer_metainfo)
 */
enum {
#ifndef BLOSC_H
  BLOSC_DOSHUFFLE = 0x1,     //!< byte-wise shuffle
  BLOSC_MEMCPYED = 0x2,      //!< plain copy
  BLOSC_DOBITSHUFFLE = 0x4,  //!< bit-wise shuffle
#endif // BLOSC_H
  BLOSC_DODELTA = 0x8,       //!< delta coding
};

/**
 * @brief Codes for new internal flags in Blosc2
 */
enum {
  BLOSC2_USEDICT = 0x1,          //!< use dictionaries with codec
  BLOSC2_BIGENDIAN = 0x2,        //!< data is in big-endian ordering
  BLOSC2_INSTR_CODEC = 0x80,     //!< codec is instrumented (mainly for development)
};

/**
 * @brief Values for different Blosc2 capabilities
 */
enum {
  BLOSC2_MAXDICTSIZE = 128 * 1024, //!< maximum size for compression dicts
  BLOSC2_MAXBLOCKSIZE = 536866816  //!< maximum size for blocks
};


enum {
  BLOSC2_DEFINED_CODECS_START = 0,
  BLOSC2_DEFINED_CODECS_STOP = 31,
  //!< Blosc-defined codecs must be between 0 - 31.
  BLOSC2_GLOBAL_REGISTERED_CODECS_START = 32,
  BLOSC2_GLOBAL_REGISTERED_CODECS_STOP = 159,
  //!< Blosc-registered codecs must be between 31 - 159.
  BLOSC2_GLOBAL_REGISTERED_CODECS = 5,
    //!< Number of Blosc-registered codecs at the moment.
  BLOSC2_USER_REGISTERED_CODECS_START = 160,
  BLOSC2_USER_REGISTERED_CODECS_STOP = 255,
  //!< User-defined codecs must be between 160 - 255.
};

/**
 * @brief Codes for the different compressors shipped with Blosc
 */
enum {
#ifndef BLOSC_H
  BLOSC_BLOSCLZ = 0,
  BLOSC_LZ4 = 1,
  BLOSC_LZ4HC = 2,
  BLOSC_ZLIB = 4,
  BLOSC_ZSTD = 5,
#endif // BLOSC_H
  BLOSC_LAST_CODEC = 6,
  //!< Determine the last codec defined by Blosc.
  BLOSC_LAST_REGISTERED_CODEC = BLOSC2_GLOBAL_REGISTERED_CODECS_START + BLOSC2_GLOBAL_REGISTERED_CODECS - 1,
  //!< Determine the last registered codec. It is used to check if a codec is registered or not.
};


// Names for the different compressors shipped with Blosc

#ifndef BLOSC_H
#define BLOSC_BLOSCLZ_COMPNAME   "blosclz"
#define BLOSC_LZ4_COMPNAME       "lz4"
#define BLOSC_LZ4HC_COMPNAME     "lz4hc"
#define BLOSC_ZLIB_COMPNAME      "zlib"
#define BLOSC_ZSTD_COMPNAME      "zstd"
#endif // BLOSC_H

/**
 * @brief Codes for compression libraries shipped with Blosc (code must be < 8)
 */
enum {
#ifndef BLOSC_H
  BLOSC_BLOSCLZ_LIB = 0,
  BLOSC_LZ4_LIB = 1,
  BLOSC_ZLIB_LIB = 3,
  BLOSC_ZSTD_LIB = 4,
#endif // BLOSC_H
  BLOSC_UDCODEC_LIB = 6,
  BLOSC_SCHUNK_LIB = 7,   //!< compressor library in super-chunk header
};

/**
 * @brief Names for the different compression libraries shipped with Blosc
 */
#ifndef BLOSC_H
#define BLOSC_BLOSCLZ_LIBNAME   "BloscLZ"
#define BLOSC_LZ4_LIBNAME       "LZ4"
#define BLOSC_ZLIB_LIBNAME      "Zlib"
#define BLOSC_ZSTD_LIBNAME      "Zstd"
#endif // BLOSC_H

/**
 * @brief The codes for compressor formats shipped with Blosc
 */
enum {
#ifndef BLOSC_H
  BLOSC_BLOSCLZ_FORMAT = BLOSC_BLOSCLZ_LIB,
  BLOSC_LZ4_FORMAT = BLOSC_LZ4_LIB,
  //!< LZ4HC and LZ4 share the same format
  BLOSC_LZ4HC_FORMAT = BLOSC_LZ4_LIB,
  BLOSC_ZLIB_FORMAT = BLOSC_ZLIB_LIB,
  BLOSC_ZSTD_FORMAT = BLOSC_ZSTD_LIB,
#endif // BLOSC_H
  BLOSC_UDCODEC_FORMAT = BLOSC_UDCODEC_LIB,
};

/**
 * @brief The version formats for compressors shipped with Blosc.
 * All versions here starts at 1
 */
enum {
#ifndef BLOSC_H
  BLOSC_BLOSCLZ_VERSION_FORMAT = 1,
  BLOSC_LZ4_VERSION_FORMAT = 1,
  BLOSC_LZ4HC_VERSION_FORMAT = 1,  /* LZ4HC and LZ4 share the same format */
  BLOSC_ZLIB_VERSION_FORMAT = 1,
  BLOSC_ZSTD_VERSION_FORMAT = 1,
#endif // BLOSC_H
  BLOSC_UDCODEC_VERSION_FORMAT = 1,
};

/**
 * @brief Split mode for blocks.
 * NEVER and ALWAYS are for experimenting with compression ratio.
 * AUTO for nearly optimal behaviour (based on heuristics).
 * FORWARD_COMPAT provides best forward compatibility (default).
 */
#ifndef BLOSC_H
enum {
  BLOSC_ALWAYS_SPLIT = 1,
  BLOSC_NEVER_SPLIT = 2,
  BLOSC_AUTO_SPLIT = 3,
  BLOSC_FORWARD_COMPAT_SPLIT = 4,
};
#endif // BLOSC_H

/**
 * @brief Offsets for fields in Blosc2 chunk header.
 */
enum {
  BLOSC2_CHUNK_VERSION = 0x0,       //!< the version for the chunk format
  BLOSC2_CHUNK_VERSIONLZ = 0x1,     //!< the version for the format of internal codec
  BLOSC2_CHUNK_FLAGS = 0x2,         //!< flags and codec info
  BLOSC2_CHUNK_TYPESIZE = 0x3,      //!< (uint8) the number of bytes of the atomic type
  BLOSC2_CHUNK_NBYTES = 0x4,        //!< (int32) uncompressed size of the buffer (this header is not included)
  BLOSC2_CHUNK_BLOCKSIZE = 0x8,     //!< (int32) size of internal blocks
  BLOSC2_CHUNK_CBYTES = 0xc,        //!< (int32) compressed size of the buffer (including this header)
  BLOSC2_CHUNK_FILTER_CODES = 0x10, //!< the codecs for the filter pipeline (1 byte per code)
  BLOSC2_CHUNK_FILTER_META = 0x18,  //!< meta info for the filter pipeline (1 byte per code)
  BLOSC2_CHUNK_BLOSC2_FLAGS = 0x1F, //!< flags specific for Blosc2 functionality
};

/**
 * @brief Run lengths for special values for chunks/frames
 */
enum {
  BLOSC2_NO_SPECIAL = 0x0,       //!< no special value
  BLOSC2_SPECIAL_ZERO = 0x1,     //!< zero special value
  BLOSC2_SPECIAL_NAN = 0x2,      //!< NaN special value
  BLOSC2_SPECIAL_VALUE = 0x3,    //!< generic special value
  BLOSC2_SPECIAL_UNINIT = 0x4,   //!< non initialized values
  BLOSC2_SPECIAL_LASTID = 0x4,   //!< last valid ID for special value (update this adequately)
  BLOSC2_SPECIAL_MASK = 0x7      //!< special value mask (prev IDs cannot be larger than this)
};

/**
 * @brief Error codes
 * Each time an error code is added here, its corresponding message error should be added in
 * print_error()
 */
enum {
  BLOSC2_ERROR_SUCCESS = 0,           //<! Success
  BLOSC2_ERROR_FAILURE = -1,          //<! Generic failure
  BLOSC2_ERROR_STREAM = -2,           //<! Bad stream
  BLOSC2_ERROR_DATA = -3,             //<! Invalid data
  BLOSC2_ERROR_MEMORY_ALLOC = -4,     //<! Memory alloc/realloc failure
  BLOSC2_ERROR_READ_BUFFER = -5,      //!< Not enough space to read
  BLOSC2_ERROR_WRITE_BUFFER = -6,     //!< Not enough space to write
  BLOSC2_ERROR_CODEC_SUPPORT = -7,    //!< Codec not supported
  BLOSC2_ERROR_CODEC_PARAM = -8,      //!< Invalid parameter supplied to codec
  BLOSC2_ERROR_CODEC_DICT = -9,       //!< Codec dictionary error
  BLOSC2_ERROR_VERSION_SUPPORT = -10, //!< Version not supported
  BLOSC2_ERROR_INVALID_HEADER = -11,  //!< Invalid value in header
  BLOSC2_ERROR_INVALID_PARAM = -12,   //!< Invalid parameter supplied to function
  BLOSC2_ERROR_FILE_READ = -13,       //!< File read failure
  BLOSC2_ERROR_FILE_WRITE = -14,      //!< File write failure
  BLOSC2_ERROR_FILE_OPEN = -15,       //!< File open failure
  BLOSC2_ERROR_NOT_FOUND = -16,       //!< Not found
  BLOSC2_ERROR_RUN_LENGTH = -17,      //!< Bad run length encoding
  BLOSC2_ERROR_FILTER_PIPELINE = -18, //!< Filter pipeline error
  BLOSC2_ERROR_CHUNK_INSERT = -19,    //!< Chunk insert failure
  BLOSC2_ERROR_CHUNK_APPEND = -20,    //!< Chunk append failure
  BLOSC2_ERROR_CHUNK_UPDATE = -21,    //!< Chunk update failure
  BLOSC2_ERROR_2GB_LIMIT = -22,       //!< Sizes larger than 2gb not supported
  BLOSC2_ERROR_SCHUNK_COPY = -23,     //!< Super-chunk copy failure
  BLOSC2_ERROR_FRAME_TYPE = -24,      //!< Wrong type for frame
  BLOSC2_ERROR_FILE_TRUNCATE = -25,   //!< File truncate failure
  BLOSC2_ERROR_THREAD_CREATE = -26,   //!< Thread or thread context creation failure
  BLOSC2_ERROR_POSTFILTER = -27,      //!< Postfilter failure
  BLOSC2_ERROR_FRAME_SPECIAL = -28,   //!< Special frame failure
  BLOSC2_ERROR_SCHUNK_SPECIAL = -29,  //!< Special super-chunk failure
  BLOSC2_ERROR_PLUGIN_IO = -30,       //!< IO plugin error
  BLOSC2_ERROR_FILE_REMOVE = -31,     //!< Remove file failure
  BLOSC2_ERROR_NULL_POINTER = -32,    //!< Pointer is null
  BLOSC2_ERROR_INVALID_INDEX = -33,   //!< Invalid index
  BLOSC2_ERROR_METALAYER_NOT_FOUND = -34,   //!< Metalayer has not been found
  BLOSC2_ERROR_MAX_BUFSIZE_EXCEEDED = -35,  //!< Max buffer size exceeded
  BLOSC2_ERROR_TUNER = -36,           //!< Tuner failure
};


#ifdef __GNUC__
#define BLOSC_ATTRIBUTE_UNUSED __attribute__((unused))
#else
#define BLOSC_ATTRIBUTE_UNUSED
#endif

static char *print_error(int rc) BLOSC_ATTRIBUTE_UNUSED;
static char *print_error(int rc) {
  switch (rc) {
    case BLOSC2_ERROR_FAILURE:
      return (char *) "Generic failure";
    case BLOSC2_ERROR_STREAM:
      return (char *) "Bad stream";
    case BLOSC2_ERROR_DATA:
      return (char *) "Invalid data";
    case BLOSC2_ERROR_MEMORY_ALLOC:
      return (char *) "Memory alloc/realloc failure";
    case BLOSC2_ERROR_READ_BUFFER:
      return (char *) "Not enough space to read";
    case BLOSC2_ERROR_WRITE_BUFFER:
      return (char *) "Not enough space to write";
    case BLOSC2_ERROR_CODEC_SUPPORT:
      return (char *) "Codec not supported";
    case BLOSC2_ERROR_CODEC_PARAM:
      return (char *) "Invalid parameter supplied to codec";
    case BLOSC2_ERROR_CODEC_DICT:
      return (char *) "Codec dictionary error";
    case BLOSC2_ERROR_VERSION_SUPPORT:
      return (char *) "Version not supported";
    case BLOSC2_ERROR_INVALID_HEADER:
      return (char *) "Invalid value in header";
    case BLOSC2_ERROR_INVALID_PARAM:
      return (char *) "Invalid parameter supplied to function";
    case BLOSC2_ERROR_FILE_READ:
      return (char *) "File read failure";
    case BLOSC2_ERROR_FILE_WRITE:
      return (char *) "File write failure";
    case BLOSC2_ERROR_FILE_OPEN:
      return (char *) "File open failure";
    case BLOSC2_ERROR_NOT_FOUND:
      return (char *) "Not found";
    case BLOSC2_ERROR_RUN_LENGTH:
      return (char *) "Bad run length encoding";
    case BLOSC2_ERROR_FILTER_PIPELINE:
      return (char *) "Filter pipeline error";
    case BLOSC2_ERROR_CHUNK_INSERT:
      return (char *) "Chunk insert failure";
    case BLOSC2_ERROR_CHUNK_APPEND:
      return (char *) "Chunk append failure";
    case BLOSC2_ERROR_CHUNK_UPDATE:
      return (char *) "Chunk update failure";
    case BLOSC2_ERROR_2GB_LIMIT:
      return (char *) "Sizes larger than 2gb not supported";
    case BLOSC2_ERROR_SCHUNK_COPY:
      return (char *) "Super-chunk copy failure";
    case BLOSC2_ERROR_FRAME_TYPE:
      return (char *) "Wrong type for frame";
    case BLOSC2_ERROR_FILE_TRUNCATE:
      return (char *) "File truncate failure";
    case BLOSC2_ERROR_THREAD_CREATE:
      return (char *) "Thread or thread context creation failure";
    case BLOSC2_ERROR_POSTFILTER:
      return (char *) "Postfilter failure";
    case BLOSC2_ERROR_FRAME_SPECIAL:
      return (char *) "Special frame failure";
    case BLOSC2_ERROR_SCHUNK_SPECIAL:
      return (char *) "Special super-chunk failure";
    case BLOSC2_ERROR_PLUGIN_IO:
      return (char *) "IO plugin error";
    case BLOSC2_ERROR_FILE_REMOVE:
      return (char *) "Remove file failure";
    case BLOSC2_ERROR_NULL_POINTER:
      return (char *) "Pointer is null";
    case BLOSC2_ERROR_INVALID_INDEX:
      return (char *) "Invalid index";
    case BLOSC2_ERROR_METALAYER_NOT_FOUND:
      return (char *) "Metalayer has not been found";
    case BLOSC2_ERROR_MAX_BUFSIZE_EXCEEDED:
      return (char *) "Maximum buffersize exceeded";
    default:
      return (char *) "Unknown error";
  }
}


/**
 * @brief Initialize the Blosc library environment.
 *
 * You must call this previous to any other Blosc call, unless you want
 * Blosc to be used simultaneously in a multi-threaded environment, in
 * which case you can use the #blosc2_compress_ctx #blosc2_decompress_ctx pair.
 *
 * @sa #blosc2_destroy
 */
BLOSC_EXPORT void blosc2_init(void);


/**
 * @brief Destroy the Blosc library environment.
 *
 * You must call this after to you are done with all the Blosc calls,
 * unless you have not used blosc2_init() before.
 *
 * @sa #blosc2_init
 */
BLOSC_EXPORT void blosc2_destroy(void);


/**
 * @brief Compress a block of data in the @p src buffer and returns the size of
 * compressed block.
 *
 * @remark Compression is memory safe and guaranteed not to write @p dest
 * more than what is specified in @p destsize.
 * There is not a minimum for @p src buffer size @p nbytes.
 * Equivalent to #blosc2_compress.
 *
 * @warning The @p src buffer and the @p dest buffer can not overlap.
 *
 * @param clevel The desired compression level and must be a number
 * between 0 (no compression) and 9 (maximum compression).
 * @param doshuffle Specifies whether the shuffle compression preconditioner
 * should be applied or not. #BLOSC_NOFILTER means not applying filters,
 * #BLOSC_SHUFFLE means applying shuffle at a byte level and
 * #BLOSC_BITSHUFFLE at a bit level (slower but *may* achieve better
 * compression).
 * @param typesize Is the number of bytes for the atomic type in binary
 * @p src buffer.  This is mainly useful for the shuffle preconditioner.
 * For implementation reasons, only a 1 < typesize < 256 will allow the
 * shuffle filter to work.  When typesize is not in this range, shuffle
 * will be silently disabled.
 * @param nbytes The number of bytes to compress in the @p src buffer.
 * @param src The buffer containing the data to compress.
 * @param dest The buffer where the compressed data will be put,
 * must have at least the size of @p destsize.
 * @param destsize The size of the dest buffer. Blosc
 * guarantees that if you set @p destsize to, at least,
 * (@p nbytes + #BLOSC2_MAX_OVERHEAD), the compression will always succeed.
 *
 * @return The number of bytes compressed.
 * If @p src buffer cannot be compressed into @p destsize, the return
 * value is zero and you should discard the contents of the @p dest
 * buffer. A negative return value means that an internal error happened. This
 * should never happen. If you see this, please report it back
 * together with the buffer data causing this and compression settings.
 *
 *
 * @par Environment variables
 * @parblock
 *
 * This function honors different environment variables to control
 * internal parameters without the need of doing that programmatically.
 * Here are the ones supported:
 *
 * * **BLOSC_CLEVEL=(INTEGER)**: This will overwrite the @p clevel parameter
 * before the compression process starts.
 *
 * * **BLOSC_SHUFFLE=[NOSHUFFLE | SHUFFLE | BITSHUFFLE]**: This will
 * overwrite the @p doshuffle parameter before the compression process
 * starts.
 *
 * * **BLOSC_DELTA=(1|0)**: This will call #blosc2_set_delta() before the
 * compression process starts.
 *
 * * **BLOSC_TYPESIZE=(INTEGER)**: This will overwrite the @p typesize
 * parameter before the compression process starts.
 *
 * * **BLOSC_COMPRESSOR=[BLOSCLZ | LZ4 | LZ4HC | ZLIB | ZSTD]**:
 * This will call #blosc1_set_compressor before the compression process starts.
 *
 * * **BLOSC_NTHREADS=(INTEGER)**: This will call
 * #blosc2_set_nthreads before the compression process starts.
 *
 * * **BLOSC_SPLITMODE=(ALWAYS | NEVER | AUTO | FORWARD_COMPAT)**:
 * This will call #blosc1_set_splitmode() before the compression process starts.
 *
 * * **BLOSC_BLOCKSIZE=(INTEGER)**: This will call
 * #blosc1_set_blocksize before the compression process starts.
 * *NOTE:* The *blocksize* is a critical parameter with
 * important restrictions in the allowed values, so use this with care.
 *
 * * **BLOSC_NOLOCK=(ANY VALUE)**: This will call #blosc2_compress_ctx under
 * the hood, with the *compressor*, *blocksize* and
 * *numinternalthreads* parameters set to the same as the last calls to
 * #blosc1_set_compressor, #blosc1_set_blocksize and
 * #blosc2_set_nthreads. *BLOSC_CLEVEL*, *BLOSC_SHUFFLE*, *BLOSC_DELTA* and
 * *BLOSC_TYPESIZE* environment vars will also be honored.
 *
 * @endparblock
 *
 * @sa #blosc1_decompress
 */
BLOSC_EXPORT int blosc1_compress(int clevel, int doshuffle, size_t typesize,
                                 size_t nbytes, const void* src, void* dest,
                                 size_t destsize);


/**
 * @brief Decompress a block of compressed data in @p src, put the result in
 * @p dest and returns the size of the decompressed block.
 *
 * @warning The @p src buffer and the @p dest buffer can not overlap.
 *
 * @remark Decompression is memory safe and guaranteed not to write the @p dest
 * buffer more than what is specified in @p destsize.
 * Similar to #blosc2_decompress.
 *
 * @remark In case you want to keep under control the number of bytes read from
 * source, you can call #blosc1_cbuffer_sizes first to check whether the
 * @p nbytes (i.e. the number of bytes to be read from @p src buffer by this
 * function) in the compressed buffer is ok with you.
 *
 * @param src The buffer to be decompressed.
 * @param dest The buffer where the decompressed data will be put.
 * @param destsize The size of the @p dest buffer.
 *
 * @return The number of bytes decompressed.
 * If an error occurs, e.g. the compressed data is corrupted or the
 * output buffer is not large enough, then a negative value
 * will be returned instead.
 *
 * @par Environment variables
 * @parblock
 * This function honors different environment variables to control
 * internal parameters without the need of doing that programmatically.
 * Here are the ones supported:
 *
 * * **BLOSC_NTHREADS=(INTEGER)**: This will call
 * #blosc2_set_nthreads before the proper decompression
 * process starts.
 *
 * * **BLOSC_NOLOCK=(ANY VALUE)**: This will call #blosc2_decompress_ctx
 * under the hood, with the *numinternalthreads* parameter set to the
 * same value as the last call to #blosc2_set_nthreads.
 *
 * @endparblock
 *
 * @sa #blosc1_compress
 */
BLOSC_EXPORT int blosc1_decompress(const void* src, void* dest, size_t destsize);


/**
 * @brief Get @p nitems (of @p typesize size) in @p src buffer starting in @p start.
 * The items are returned in @p dest buffer, which has to have enough
 * space for storing all items.
 *
 * @remark The function #blosc2_getitem is a more complete and secure version.
 *
 * @param src The compressed buffer from data will be decompressed.
 * @param start The position of the first item (of @p typesize size) from where data
 * will be retrieved.
 * @param nitems The number of items (of @p typesize size) that will be retrieved.
 * @param dest The buffer where the decompressed data retrieved will be put.
 *
 * @return The number of bytes copied to @p dest or a negative value if
 * some error happens.
 */
BLOSC_EXPORT int blosc1_getitem(const void* src, int start, int nitems, void* dest);

/**
 * @brief Get @p nitems (of @p typesize size) in @p src buffer starting in @p start.
 * The items are returned in @p dest buffer. The dest buffer should have enough space
 * for storing all items. This function is a more secure version of #blosc1_getitem.
 *
 * @param src The compressed buffer holding the data to be retrieved.
 * @param srcsize Size of the compressed buffer.
 * @param start The position of the first item (of @p typesize size) from where data
 * will be retrieved.
 * @param nitems The number of items (of @p typesize size) that will be retrieved.
 * @param dest The buffer where the retrieved data will be stored decompressed.
 * @param destsize Size of the buffer where retrieved data will be stored.
 *
 * @return The number of bytes copied to @p dest or a negative value if
 * some error happens.
 */
BLOSC_EXPORT int blosc2_getitem(const void* src, int32_t srcsize, int start, int nitems,
                                void* dest, int32_t destsize);

/**
  Pointer to a callback function that executes `dojob(jobdata + i*jobdata_elsize)` for `i = 0 to numjobs-1`,
  possibly in parallel threads (but not returning until all `dojob` calls have returned).   This allows the
  caller to provide a custom threading backend as an alternative to the default Blosc-managed threads.
  `callback_data` is passed through from `blosc2_set_threads_callback`.
 */
typedef void (*blosc_threads_callback)(void *callback_data, void (*dojob)(void *), int numjobs, size_t jobdata_elsize, void *jobdata);

/**
  * @brief Set the threading backend for parallel compression/decompression to use @p callback to execute work
  * instead of using the Blosc-managed threads.  The @p callback_data argument is passed through to the callback.

  * @param callback: the callback to use. Passing `NULL` uses the default Blosc threading backend.
  * @param callback_data: the callback data.
  *
  * @warning This function is *not* thread-safe and should be called
  * before any other Blosc function: it affects all Blosc contexts.
  * @sa https://github.com/Blosc/c-blosc2/pull/81
  */
BLOSC_EXPORT void blosc2_set_threads_callback(blosc_threads_callback callback, void *callback_data);


/**
 * @brief Returns the current number of threads that are used for
 * compression/decompression.
 */
BLOSC_EXPORT int16_t blosc2_get_nthreads(void);


/**
 * @brief Initialize a pool of threads for compression/decompression. If
 * @p nthreads is 1, then the serial version is chosen and a possible
 * previous existing pool is ended. If this is not called, @p nthreads
 * is set to 1 internally.
 *
 * @param nthreads The number of threads to use.
 *
 * @return The previous number of threads.
 */
BLOSC_EXPORT int16_t blosc2_set_nthreads(int16_t nthreads);


/**
 * @brief Get the current compressor that is used for compression.
 *
 * @return The string identifying the compressor being used.
 */
BLOSC_EXPORT const char* blosc1_get_compressor(void);


/**
 * @brief Select the compressor to be used. The supported ones are "blosclz",
 * "lz4", "lz4hc", "zlib" and "ztsd". If this function is not
 * called, then "blosclz" will be used.
 *
 * @param compname The name identifier of the compressor to be set.
 *
 * @return The code for the compressor (>=0). In case the compressor
 * is not recognized, or there is not support for it in this build,
 * it returns a -1.
 */
BLOSC_EXPORT int blosc1_set_compressor(const char* compname);


/**
 * @brief Select the delta coding filter to be used.
 *
 * @param dodelta A value >0 will activate the delta filter.
 * If 0, it will be de-activated
 *
 * This call should always succeed.
 */
BLOSC_EXPORT void blosc2_set_delta(int dodelta);


/**
 * @brief Get the compressor name associated with the compressor code.
 *
 * @param compcode The code identifying the compressor
 * @param compname The pointer to a string where the compressor name will be put.
 *
 * @return The compressor code. If the compressor code is not recognized,
 * or there is not support for it in this build, -1 is returned.
 */
BLOSC_EXPORT int blosc2_compcode_to_compname(int compcode, const char** compname);


/**
 * @brief Get the compressor code associated with the compressor name.
 *
 * @param compname The string containing the compressor name.
 *
 * @return The compressor code. If the compressor name is not recognized,
 * or there is not support for it in this build, -1 is returned instead.
 */
BLOSC_EXPORT int blosc2_compname_to_compcode(const char* compname);


/**
 * @brief Get a list of compressors supported in the current build.
 *
 * @return The comma separated string with the list of compressor names
 * supported.
 *
 * This function does not leak, so you should not free() the returned
 * list.
 *
 * This function should always succeed.
 */
BLOSC_EXPORT const char* blosc2_list_compressors(void);


/**
 * @brief Get the version of Blosc in string format.
 *
 * @return The string with the current Blosc version.
 * Useful for dynamic libraries.
 */
BLOSC_EXPORT const char* blosc2_get_version_string(void);


/**
 * @brief Get info from compression libraries included in the current build.
 *
 * @param compname The compressor name that you want info from.
 * @param complib The pointer to a string where the
 * compression library name, if available, will be put.
 * @param version The pointer to a string where the
 * compression library version, if available, will be put.
 *
 * @warning You are in charge of the @p complib and @p version strings,
 * you should free() them so as to avoid leaks.
 *
 * @return The code for the compression library (>=0). If it is not supported,
 * this function returns -1.
 */
BLOSC_EXPORT int blosc2_get_complib_info(const char* compname, char** complib,
                                         char** version);


/**
 * @brief Free possible memory temporaries and thread resources. Use this
 * when you are not going to use Blosc for a long while.
 *
 * @return A 0 if succeeds, in case of problems releasing the resources,
 * it returns a negative number.
 */
BLOSC_EXPORT int blosc2_free_resources(void);


/**
 * @brief Get information about a compressed buffer, namely the number of
 * uncompressed bytes (@p nbytes) and compressed (@p cbytes). It also
 * returns the @p blocksize (which is used internally for doing the
 * compression by blocks).
 *
 * @remark Equivalent to function #blosc2_cbuffer_sizes.
 *
 * @param cbuffer The buffer of compressed data.
 * @param nbytes The pointer where the number of uncompressed bytes will be put.
 * @param cbytes The pointer where the number of compressed bytes will be put.
 * @param blocksize The pointer where the block size will be put.
 *
 * You only need to pass the first BLOSC_MIN_HEADER_LENGTH bytes of a
 * compressed buffer for this call to work.
 *
 * This function should always succeed.
 */
BLOSC_EXPORT void blosc1_cbuffer_sizes(const void* cbuffer, size_t* nbytes,
                                       size_t* cbytes, size_t* blocksize);
/**
 * @brief Get information about a compressed buffer, namely the number of
 * uncompressed bytes (@p nbytes) and compressed (@p cbytes). It also
 * returns the @p blocksize (which is used internally for doing the
 * compression by blocks).
 *
 * @param cbuffer The buffer of compressed data.
 * @param nbytes The pointer where the number of uncompressed bytes will be put.
 * @param cbytes The pointer where the number of compressed bytes will be put.
 * @param blocksize The pointer where the block size will be put.
 *
 * @note: if any of the nbytes, cbytes or blocksize is NULL, it will not be returned.
 *
 * You only need to pass the first BLOSC_MIN_HEADER_LENGTH bytes of a
 * compressed buffer for this call to work.
 *
 * @return On failure, returns negative value.
 */
BLOSC_EXPORT int blosc2_cbuffer_sizes(const void* cbuffer, int32_t* nbytes,
                                      int32_t* cbytes, int32_t* blocksize);

/**
 * @brief Checks that the compressed buffer starting at @p cbuffer of length @p cbytes
 * may contain valid blosc compressed data, and that it is safe to call
 * blosc1_decompress/blosc1_getitem.
 * On success, returns 0 and sets @p nbytes to the size of the uncompressed data.
 * This does not guarantee that the decompression function won't return an error,
 * but does guarantee that it is safe to attempt decompression.
 *
 * @param cbuffer The buffer of compressed data.
 * @param cbytes The number of compressed bytes.
 * @param nbytes The pointer where the number of uncompressed bytes will be put.
 *
 * @return On failure, returns negative value.
 */
BLOSC_EXPORT int blosc1_cbuffer_validate(const void* cbuffer, size_t cbytes,
                                         size_t* nbytes);

/**
 * @brief Get information about a compressed buffer, namely the type size
 * (@p typesize), as well as some internal @p flags.
 *
 * @param cbuffer The buffer of compressed data.
 * @param typesize The pointer where the type size will be put.
 * @param flags The pointer of the integer where the additional info is encoded.
 * The @p flags is a set of bits, where the currently used ones are:
 *   * bit 0: whether the shuffle filter has been applied or not
 *   * bit 1: whether the internal buffer is a pure memcpy or not
 *   * bit 2: whether the bitshuffle filter has been applied or not
 *   * bit 3: whether the delta coding filter has been applied or not
 *
 * You can use the @p BLOSC_DOSHUFFLE, @p BLOSC_DOBITSHUFFLE, @p BLOSC_DODELTA
 * and @p BLOSC_MEMCPYED symbols for extracting the interesting bits
 * (e.g. @p flags & @p BLOSC_DOSHUFFLE says whether the buffer is byte-shuffled
 * or not).
 *
 * This function should always succeed.
 */
BLOSC_EXPORT void blosc1_cbuffer_metainfo(const void* cbuffer, size_t* typesize,
                                          int* flags);


/**
 * @brief Get information about a compressed buffer, namely the internal
 * Blosc format version (@p version) and the format for the internal
 * Lempel-Ziv compressor used (@p versionlz).
 *
 * @param cbuffer The buffer of compressed data.
 * @param version The pointer where the Blosc format version will be put.
 * @param versionlz The pointer where the Lempel-Ziv version will be put.
 *
 * This function should always succeed.
 */
BLOSC_EXPORT void blosc2_cbuffer_versions(const void* cbuffer, int* version,
                                          int* versionlz);


/**
 * @brief Get the compressor library/format used in a compressed buffer.
 *
 * @param cbuffer The buffer of compressed data.
 *
 * @return The string identifying the compressor library/format used.
 *
 * This function should always succeed.
 */
BLOSC_EXPORT const char* blosc2_cbuffer_complib(const void* cbuffer);

/*********************************************************************
  Structures and functions related with user-defined input/output.
*********************************************************************/

enum {
  BLOSC2_IO_FILESYSTEM = 0,
  BLOSC2_IO_FILESYSTEM_MMAP = 1,
  BLOSC_IO_LAST_BLOSC_DEFINED = 2,  // sentinel
  BLOSC_IO_LAST_REGISTERED = 32,  // sentinel
};

enum {
  BLOSC2_IO_BLOSC_DEFINED = 32,
  BLOSC2_IO_REGISTERED = 160,
  BLOSC2_IO_USER_DEFINED = 256
};

typedef void*   (*blosc2_open_cb)(const char *urlpath, const char *mode, void *params);
typedef int     (*blosc2_close_cb)(void *stream);
typedef int64_t (*blosc2_size_cb)(void *stream);
typedef int64_t (*blosc2_write_cb)(const void *ptr, int64_t size, int64_t nitems, int64_t position, void *stream);
typedef int64_t (*blosc2_read_cb)(void **ptr, int64_t size, int64_t nitems, int64_t position, void *stream);
typedef int     (*blosc2_truncate_cb)(void *stream, int64_t size);
typedef int     (*blosc2_destroy_cb)(void *params);


/*
 * Input/Output callbacks.
 */
typedef struct {
  uint8_t id;
  //!< The IO identifier.
  char* name;
  //!< The IO name.
  bool is_allocation_necessary;
  //!< If true, the caller needs to allocate data for the read function (ptr argument). If false, the read function
  //!< takes care of memory allocation and stores the address in the allocated_ptr argument.
  blosc2_open_cb open;
  //!< The IO open callback.
  blosc2_close_cb close;
  //!< The IO close callback.
  blosc2_size_cb size;
  //!< The IO size callback.
  blosc2_write_cb write;
  //!< The IO write callback.
  blosc2_read_cb read;
  //!< The IO read callback.
  blosc2_truncate_cb truncate;
  //!< The IO truncate callback.
  blosc2_destroy_cb destroy;
  //!< The IO destroy callback (called in the end when finished with the schunk).
} blosc2_io_cb;


/*
 * Input/Output parameters.
 */
typedef struct {
  uint8_t id;
  const char *name;
  //!< The IO identifier.
  void *params;
  //!< The IO parameters.
} blosc2_io;

static const blosc2_io BLOSC2_IO_DEFAULTS = {
  /* .id = */ BLOSC2_IO_FILESYSTEM,
  /* .name = */ "filesystem",
  /* .params = */ NULL,
};


/**
 * @brief Register a user-defined input/output callbacks in Blosc.
 *
 * @param io The callbacks API to register.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_register_io_cb(const blosc2_io_cb *io);

/**
 * @brief Get a user-defined input/output callback in Blosc.
 *
 * @param id The id of the callback to get.
 *
 * @return A pointer containing the desired callback if success. Else a NULL pointer is returned.
 */
BLOSC_EXPORT blosc2_io_cb *blosc2_get_io_cb(uint8_t id);

/*********************************************************************
  Structures and functions related with contexts.
*********************************************************************/

typedef struct blosc2_context_s blosc2_context;   /* opaque type */

typedef struct {
  int (*init)(void * config, blosc2_context* cctx, blosc2_context* dctx);
  //!< Initialize tuner. Keep in mind dctx may be NULL. This should memcpy the cctx->tuner_params.
  int (*next_blocksize)(blosc2_context * context);
  //!< Only compute the next blocksize. Only it is executed if tuner is not initialized.
  int (*next_cparams)(blosc2_context * context);
  //!< Compute the next cparams. Only is executed if tuner is initialized.
  int (*update)(blosc2_context * context, double ctime);
  //!< Update the tuner parameters.
  int (*free)(blosc2_context * context);
  //!< Free the tuner.
  int id;
  //!< The tuner id
  char *name;
  //!< The tuner name
} blosc2_tuner;


/**
 * @brief Register locally a user-defined tuner in Blosc.
 *
 * @param tuner The tuner to register.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_register_tuner(blosc2_tuner *tuner);


/**
 * @brief The parameters for a prefilter function.
 *
 */
typedef struct {
  void *user_data;  // user-provided info (optional)
  const uint8_t *input;  // the input buffer
  uint8_t *output;  // the output buffer
  int32_t output_size;  // the output size (in bytes)
  int32_t output_typesize;  // the output typesize
  int32_t output_offset; // offset to reach the start of the output buffer
  int64_t nchunk;  // the current nchunk in associated schunk (if exists; if not -1)
  int32_t nblock;  // the current nblock in associated chunk
  int32_t tid;  // thread id
  uint8_t *ttmp;  // a temporary that is able to hold several blocks for the output and is private for each thread
  size_t ttmp_nbytes;  // the size of the temporary in bytes
  blosc2_context *ctx;  // the compression context
} blosc2_prefilter_params;

/**
 * @brief The parameters for a postfilter function.
 *
 */
typedef struct {
  void *user_data;  // user-provided info (optional)
  const uint8_t *input;  // the input buffer
  uint8_t *output;  // the output buffer
  int32_t size;  // the input size (in bytes)
  int32_t typesize;  // the input typesize
  int32_t offset;  // offset to reach the start of the input buffer
  int64_t nchunk;  // the current nchunk in associated schunk (if exists; if not -1)
  int32_t nblock;  // the current nblock in associated chunk
  int32_t tid;  // thread id
  uint8_t *ttmp;  // a temporary that is able to hold several blocks for the output and is private for each thread
  size_t ttmp_nbytes;  // the size of the temporary in bytes
  blosc2_context *ctx;  // the decompression context
} blosc2_postfilter_params;

/**
 * @brief The type of the prefilter function.
 *
 * If the function call is successful, the return value should be 0; else, a negative value.
 */
typedef int (*blosc2_prefilter_fn)(blosc2_prefilter_params* params);

/**
 * @brief The type of the postfilter function.
 *
 * If the function call is successful, the return value should be 0; else, a negative value.
 */
typedef int (*blosc2_postfilter_fn)(blosc2_postfilter_params* params);

/**
 * @brief The parameters for creating a context for compression purposes.
 *
 * In parenthesis it is shown the default value used internally when a 0
 * (zero) in the fields of the struct is passed to a function.
 */
typedef struct {
  uint8_t compcode;
  //!< The compressor codec.
  uint8_t compcode_meta;
  //!< The metadata for the compressor codec.
  uint8_t clevel;
  //!< The compression level (5).
  int use_dict;
  //!< Use dicts or not when compressing (only for ZSTD).
  int32_t typesize;
  //!< The type size (8).
  int16_t nthreads;
  //!< The number of threads to use internally (1).
  int32_t blocksize;
  //!< The requested size of the compressed blocks (0 means automatic).
  int32_t splitmode;
  //!< Whether the blocks should be split or not.
  void* schunk;
  //!< The associated schunk, if any (NULL).
  uint8_t filters[BLOSC2_MAX_FILTERS];
  //!< The (sequence of) filters.
  uint8_t filters_meta[BLOSC2_MAX_FILTERS];
  //!< The metadata for filters.
  blosc2_prefilter_fn prefilter;
  //!< The prefilter function.
  blosc2_prefilter_params *preparams;
  //!< The prefilter parameters.
  void *tuner_params;
  //!< Tune configuration.
  int tuner_id;
  //!< The tuner id.
  bool instr_codec;
  //!< Whether the codec is instrumented or not
  void *codec_params;
  //!< User defined parameters for the codec
  void *filter_params[BLOSC2_MAX_FILTERS];
  //!< User defined parameters for the filters
} blosc2_cparams;

/**
 * @brief Default struct for compression params meant for user initialization.
 */
static const blosc2_cparams BLOSC2_CPARAMS_DEFAULTS = {
        BLOSC_BLOSCLZ, 0, 5, 0, 8, 1, 0,
        BLOSC_FORWARD_COMPAT_SPLIT, NULL,
        {0, 0, 0, 0, 0, BLOSC_SHUFFLE},
        {0, 0, 0, 0, 0, 0},
        NULL, NULL, NULL, 0, 0,
        NULL, {NULL, NULL, NULL, NULL, NULL, NULL}
        };


/**
  @brief The parameters for creating a context for decompression purposes.

  In parenthesis it is shown the default value used internally when a 0
  (zero) in the fields of the struct is passed to a function.
 */
typedef struct {
  int16_t nthreads;
  //!< The number of threads to use internally (1).
  void* schunk;
  //!< The associated schunk, if any (NULL).
  blosc2_postfilter_fn postfilter;
  //!< The postfilter function.
  blosc2_postfilter_params *postparams;
  //!< The postfilter parameters.
} blosc2_dparams;

/**
 * @brief Default struct for decompression params meant for user initialization.
 */
static const blosc2_dparams BLOSC2_DPARAMS_DEFAULTS = {1, NULL, NULL, NULL};

/**
 * @brief Create a context for @a *_ctx() compression functions.
 *
 * @param cparams The blosc2_cparams struct with the compression parameters.
 *
 * @return A pointer to the new context. NULL is returned if this fails.
 *
 * @note This supports the same environment variables than #blosc2_compress
 * for overriding the programmatic compression values.
 *
 * @sa #blosc2_compress
 */
BLOSC_EXPORT blosc2_context* blosc2_create_cctx(blosc2_cparams cparams);

/**
 * @brief Create a context for *_ctx() decompression functions.
 *
 * @param dparams The blosc2_dparams struct with the decompression parameters.
 *
 * @return A pointer to the new context. NULL is returned if this fails.
 *
 * @note This supports the same environment variables than #blosc2_decompress
 * for overriding the programmatic decompression values.
 *
 * @sa #blosc2_decompress
 *
 */
BLOSC_EXPORT blosc2_context* blosc2_create_dctx(blosc2_dparams dparams);

/**
 * @brief Free the resources associated with a context.
 *
 * @param context The context to free.
 *
 * This function should always succeed and is valid for contexts meant for
 * both compression and decompression.
 */
BLOSC_EXPORT void blosc2_free_ctx(blosc2_context* context);

/**
 * @brief Create a @p cparams associated to a context.
 *
 * @param ctx The context from where to extract the compression parameters.
 * @param cparams The pointer where the compression params will be stored.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_ctx_get_cparams(blosc2_context *ctx, blosc2_cparams *cparams);

/**
 * @brief Create a @p dparams associated to a context.
 *
 * @param ctx The context from where to extract the decompression parameters.
 * @param dparams The pointer where the decompression params will be stored.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_ctx_get_dparams(blosc2_context *ctx, blosc2_dparams *dparams);

/**
 * @brief Set a maskout so as to avoid decompressing specified blocks.
 *
 * @param ctx The decompression context to update.
 *
 * @param maskout The boolean mask for the blocks where decompression
 * is to be avoided.
 *
 * @remark The maskout is valid for contexts *only* meant for decompressing
 * a chunk via #blosc2_decompress_ctx.  Once a call to #blosc2_decompress_ctx
 * is done, this mask is reset so that next call to #blosc2_decompress_ctx
 * will decompress the whole chunk.
 *
 * @param nblocks The number of blocks in maskout above.
 *
 * @return If success, a 0 is returned.  An error is signaled with a negative int.
 *
 */
BLOSC_EXPORT int blosc2_set_maskout(blosc2_context *ctx, bool *maskout, int nblocks);

/**
 * @brief Compress a block of data in the @p src buffer and returns the size of
 * compressed block.
 *
 * @remark Compression is memory safe and guaranteed not to write @p dest
 * more than what is specified in @p destsize.
 * There is not a minimum for @p src buffer size @p nbytes.
 *
 * @warning The @p src buffer and the @p dest buffer can not overlap.
 *
 * @param clevel The desired compression level and must be a number
 * between 0 (no compression) and 9 (maximum compression).
 * @param doshuffle Specifies whether the shuffle compression preconditioner
 * should be applied or not. #BLOSC_NOFILTER means not applying filters,
 * #BLOSC_SHUFFLE means applying shuffle at a byte level and
 * #BLOSC_BITSHUFFLE at a bit level (slower but *may* achieve better
 * compression).
 * @param typesize Is the number of bytes for the atomic type in binary
 * @p src buffer.  This is mainly useful for the shuffle preconditioner.
 * For implementation reasons, only a 1 < typesize < 256 will allow the
 * shuffle filter to work.  When typesize is not in this range, shuffle
 * will be silently disabled.
 * @param src The buffer containing the data to compress.
 * @param srcsize The number of bytes to compress in the @p src buffer.
 * @param dest The buffer where the compressed data will be put,
 * must have at least the size of @p destsize.
 * @param destsize The size of the dest buffer. Blosc
 * guarantees that if you set @p destsize to, at least,
 * (@p nbytes + #BLOSC2_MAX_OVERHEAD), the compression will always succeed.
 *
 * @return The number of bytes compressed.
 * If @p src buffer cannot be compressed into @p destsize, the return
 * value is zero and you should discard the contents of the @p dest
 * buffer. A negative return value means that either a parameter is not correct
 * or that an internal error happened. Set the BLOSC_TRACE environment variable
 * for getting more info on what is happening. If the error is not related with
 * wrong params, please report it back together with the buffer data causing this,
 * as well as the compression params used.
*/
/*
 * Environment variables
 * _____________________
 *
 * *blosc2_compress()* honors different environment variables to control
 * internal parameters without the need of doing that programmatically.
 * Here are the ones supported:
 *
 * **BLOSC_CLEVEL=(INTEGER)**: This will overwrite the @p clevel parameter
 * before the compression process starts.
 *
 * **BLOSC_SHUFFLE=[NOSHUFFLE | SHUFFLE | BITSHUFFLE]**: This will
 * overwrite the *doshuffle* parameter before the compression process
 * starts.
 *
 * **BLOSC_DELTA=(1|0)**: This will call *blosc2_set_delta()^* before the
 * compression process starts.
 *
 * **BLOSC_TYPESIZE=(INTEGER)**: This will overwrite the *typesize*
 * parameter before the compression process starts.
 *
 * **BLOSC_COMPRESSOR=[BLOSCLZ | LZ4 | LZ4HC | ZLIB | ZSTD]**:
 * This will call #blosc_set_compressor before the compression process starts.
 *
 * **BLOSC_NTHREADS=(INTEGER)**: This will call
 * #blosc_set_nthreads before the compression process
 * starts.
 *
 * **BLOSC_SPLITMODE=(ALWAYS | NEVER | AUTO | FORWARD_COMPAT)**:
 * This will call #blosc1_set_splitmode() before the compression process starts.
 *
 * **BLOSC_BLOCKSIZE=(INTEGER)**: This will call
 * #blosc_set_blocksize before the compression process starts.
 * *NOTE:* The blocksize is a critical parameter with
 * important restrictions in the allowed values, so use this with care.
 *
 * **BLOSC_NOLOCK=(ANY VALUE)**: This will call #blosc2_compress_ctx under
 * the hood, with the *compressor*, *blocksize* and
 * *numinternalthreads* parameters set to the same as the last calls to
 * #blosc1_set_compressor, #blosc1_set_blocksize and
 * #blosc2_set_nthreads. *BLOSC_CLEVEL*, *BLOSC_SHUFFLE*, *BLOSC_DELTA* and
 * *BLOSC_TYPESIZE* environment vars will also be honored.
 *
 */
BLOSC_EXPORT int blosc2_compress(int clevel, int doshuffle, int32_t typesize,
                                 const void* src, int32_t srcsize, void* dest,
                                 int32_t destsize);


/**
 * @brief Decompress a block of compressed data in @p src, put the result in
 * @p dest and returns the size of the decompressed block.
 *
 * @warning The @p src buffer and the @p dest buffer can not overlap.
 *
 * @remark Decompression is memory safe and guaranteed not to write the @p dest
 * buffer more than what is specified in @p destsize.
 *
 * @remark In case you want to keep under control the number of bytes read from
 * source, you can call #blosc1_cbuffer_sizes first to check whether the
 * @p nbytes (i.e. the number of bytes to be read from @p src buffer by this
 * function) in the compressed buffer is ok with you.
 *
 * @param src The buffer to be decompressed.
 * @param srcsize The size of the buffer to be decompressed.
 * @param dest The buffer where the decompressed data will be put.
 * @param destsize The size of the @p dest buffer.
 *
 * @return The number of bytes decompressed.
 * If an error occurs, e.g. the compressed data is corrupted or the
 * output buffer is not large enough, then a negative value
 * will be returned instead.
*/
/*
 * Environment variables
 * _____________________
 *
 * *blosc1_decompress* honors different environment variables to control
 * internal parameters without the need of doing that programmatically.
 * Here are the ones supported:
 *
 * **BLOSC_NTHREADS=(INTEGER)**: This will call
 * *blosc_set_nthreads(BLOSC_NTHREADS)* before the proper decompression
 * process starts.
 *
 * **BLOSC_NOLOCK=(ANY VALUE)**: This will call *blosc2_decompress_ctx*
 * under the hood, with the *numinternalthreads* parameter set to the
 * same value as the last call to *blosc2_set_nthreads*.
 *
 */
BLOSC_EXPORT int blosc2_decompress(const void* src, int32_t srcsize,
                                   void* dest, int32_t destsize);

/**
 * @brief Context interface to Blosc compression. This does not require a call
 * to #blosc2_init and can be called from multithreaded applications
 * without the global lock being used, so allowing Blosc be executed
 * simultaneously in those scenarios.
 *
 * @param context A blosc2_context struct with the different compression params.
 * @param src The buffer containing the data to be compressed.
 * @param srcsize The number of bytes to be compressed from the @p src buffer.
 * @param dest The buffer where the compressed data will be put.
 * @param destsize The size in bytes of the @p dest buffer.
 *
 * @return The number of bytes compressed.
 * If @p src buffer cannot be compressed into @p destsize, the return
 * value is zero and you should discard the contents of the @p dest
 * buffer.  A negative return value means that an internal error happened.
 * It could happen that context is not meant for compression (which is stated in stderr).
 * Otherwise, please report it back together with the buffer data causing this
 * and compression settings.
 */
BLOSC_EXPORT int blosc2_compress_ctx(
        blosc2_context* context, const void* src, int32_t srcsize, void* dest,
        int32_t destsize);


/**
 * @brief Context interface to Blosc decompression. This does not require a
 * call to #blosc2_init and can be called from multithreaded
 * applications without the global lock being used, so allowing Blosc
 * be executed simultaneously in those scenarios.
 *
 * @param context The blosc2_context struct with the different compression params.
 * @param src The buffer of compressed data.
 * @param srcsize The length of buffer of compressed data.
 * @param dest The buffer where the decompressed data will be put.
 * @param destsize The size in bytes of the @p dest buffer.
 *
 * @warning The @p src buffer and the @p dest buffer can not overlap.
 *
 * @remark Decompression is memory safe and guaranteed not to write the @p dest
 * buffer more than what is specified in @p destsize.
 *
 * @remark In case you want to keep under control the number of bytes read from
 * source, you can call #blosc1_cbuffer_sizes first to check the @p nbytes
 * (i.e. the number of bytes to be read from @p src buffer by this function)
 * in the compressed buffer.
 *
 * @remark If #blosc2_set_maskout is called prior to this function, its
 * @p block_maskout parameter will be honored for just *one single* shot;
 * i.e. the maskout in context will be automatically reset to NULL, so
 * mask won't be used next time (unless #blosc2_set_maskout is called again).
 *
 * @return The number of bytes decompressed (i.e. the maskout blocks are not
 * counted). If an error occurs, e.g. the compressed data is corrupted,
 * @p destsize is not large enough or context is not meant for decompression,
 * then a negative value will be returned instead.
 */
BLOSC_EXPORT int blosc2_decompress_ctx(blosc2_context* context, const void* src,
                                       int32_t srcsize, void* dest, int32_t destsize);

/**
 * @brief Create a chunk made of zeros.
 *
 * @param cparams The compression parameters.
 * @param nbytes The size (in bytes) of the chunk.
 * @param dest The buffer where the data chunk will be put.
 * @param destsize The size (in bytes) of the @p dest buffer;
 * must be BLOSC_EXTENDED_HEADER_LENGTH at least.
 *
 * @return The number of bytes compressed (BLOSC_EXTENDED_HEADER_LENGTH).
 * If negative, there has been an error and @p dest is unusable.
 * */
BLOSC_EXPORT int blosc2_chunk_zeros(blosc2_cparams cparams, int32_t nbytes,
                                    void* dest, int32_t destsize);


/**
 * @brief Create a chunk made of nans.
 *
 * @param cparams The compression parameters;
 * only 4 bytes (float) and 8 bytes (double) are supported.
 * @param nbytes The size (in bytes) of the chunk.
 * @param dest The buffer where the data chunk will be put.
 * @param destsize The size (in bytes) of the @p dest buffer;
 * must be BLOSC_EXTENDED_HEADER_LENGTH at least.
 *
 * @note Whether the NaNs are floats or doubles will be given by the typesize.
 *
 * @return The number of bytes compressed (BLOSC_EXTENDED_HEADER_LENGTH).
 * If negative, there has been an error and @p dest is unusable.
 * */
BLOSC_EXPORT int blosc2_chunk_nans(blosc2_cparams cparams, int32_t nbytes,
                                   void* dest, int32_t destsize);


/**
 * @brief Create a chunk made of repeated values.
 *
 * @param cparams The compression parameters.
 * @param nbytes The size (in bytes) of the chunk.
 * @param dest The buffer where the data chunk will be put.
 * @param destsize The size (in bytes) of the @p dest buffer.
 * @param repeatval A pointer to the repeated value (little endian).
 * The size of the value is given by @p cparams.typesize param.
 *
 * @return The number of bytes compressed (BLOSC_EXTENDED_HEADER_LENGTH + typesize).
 * If negative, there has been an error and @p dest is unusable.
 * */
BLOSC_EXPORT int blosc2_chunk_repeatval(blosc2_cparams cparams, int32_t nbytes,
                                        void* dest, int32_t destsize, const void* repeatval);


/**
 * @brief Create a chunk made of uninitialized values.
 *
 * @param cparams The compression parameters.
 * @param nbytes The size (in bytes) of the chunk.
 * @param dest The buffer where the data chunk will be put.
 * @param destsize The size (in bytes) of the @p dest buffer;
 * must be BLOSC_EXTENDED_HEADER_LENGTH at least.
 *
 * @return The number of bytes compressed (BLOSC_EXTENDED_HEADER_LENGTH).
 * If negative, there has been an error and @p dest is unusable.
 * */
BLOSC_EXPORT int blosc2_chunk_uninit(blosc2_cparams cparams, int32_t nbytes,
                                     void* dest, int32_t destsize);


/**
 * @brief Context interface counterpart for #blosc1_getitem.
 *
 * @param context Context pointer.
 * @param src The compressed buffer from data will be decompressed.
 * @param srcsize Compressed buffer length.
 * @param start The position of the first item (of @p typesize size) from where data
 * will be retrieved.
 * @param nitems The number of items (of @p typesize size) that will be retrieved.
 * @param dest The buffer where the decompressed data retrieved will be put.
 * @param destsize Output buffer length.
 *
 * @return The number of bytes copied to @p dest or a negative value if
 * some error happens.
 */
BLOSC_EXPORT int blosc2_getitem_ctx(blosc2_context* context, const void* src,
                                    int32_t srcsize, int start, int nitems, void* dest,
                                    int32_t destsize);


/*********************************************************************
  Super-chunk related structures and functions.
*********************************************************************/

#define BLOSC2_MAX_METALAYERS 16
#define BLOSC2_METALAYER_NAME_MAXLEN 31

// Allow for a reasonable number of vl metalayers
// max is 64 * 1024 due to msgpack map 16 in frame
// mem usage 8 * 1024 entries for blosc2_schunk.vlmetalayers[] is 64 KB
#define BLOSC2_MAX_VLMETALAYERS (8 * 1024)
#define BLOSC2_VLMETALAYERS_NAME_MAXLEN BLOSC2_METALAYER_NAME_MAXLEN

/**
 * @brief This struct is meant for holding storage parameters for a
 * for a blosc2 container, allowing to specify, for example, how to interpret
 * the contents included in the schunk.
 */
typedef struct {
    bool contiguous;
    //!< Whether the chunks are contiguous or sparse.
    char* urlpath;
    //!< The path for persistent storage. If NULL, that means in-memory.
    blosc2_cparams* cparams;
    //!< The compression params when creating a schunk.
    //!< If NULL, sensible defaults are used depending on the context.
    blosc2_dparams* dparams;
    //!< The decompression params when creating a schunk.
    //!< If NULL, sensible defaults are used depending on the context.
    blosc2_io *io;
    //!< Input/output backend.
} blosc2_storage;

/**
 * @brief Default struct for #blosc2_storage meant for user initialization.
 */
static const blosc2_storage BLOSC2_STORAGE_DEFAULTS = {false, NULL, NULL, NULL, NULL};

/**
 * @brief Get default struct for compression params meant for user initialization.
 */
BLOSC_EXPORT blosc2_cparams blosc2_get_blosc2_cparams_defaults(void);

/**
 * @brief Get default struct for decompression params meant for user initialization.
 */
BLOSC_EXPORT blosc2_dparams blosc2_get_blosc2_dparams_defaults(void);

/**
 * @brief Get default struct for #blosc2_storage meant for user initialization.
 */
BLOSC_EXPORT blosc2_storage blosc2_get_blosc2_storage_defaults(void);

/**
 * @brief Get default struct for #blosc2_io meant for user initialization.
 */
BLOSC_EXPORT blosc2_io blosc2_get_blosc2_io_defaults(void);

/**
 * @brief Get default struct for #blosc2_stdio_mmap meant for user initialization.
 */
BLOSC_EXPORT blosc2_stdio_mmap blosc2_get_blosc2_stdio_mmap_defaults(void);

typedef struct blosc2_frame_s blosc2_frame;   /* opaque type */

/**
 * @brief This struct is meant to store metadata information inside
 * a #blosc2_schunk, allowing to specify, for example, how to interpret
 * the contents included in the schunk.
 */
typedef struct blosc2_metalayer {
  char* name;          //!< The metalayer identifier for Blosc client (e.g. Blosc2 NDim).
  uint8_t* content;    //!< The serialized (msgpack preferably) content of the metalayer.
  int32_t content_len; //!< The length in bytes of the content.
} blosc2_metalayer;

/**
 * @brief This struct is the standard container for Blosc 2 compressed data.
 *
 * This is essentially a container for Blosc 1 chunks of compressed data,
 * and it allows to overcome the 32-bit limitation in Blosc 1. Optionally,
 * a #blosc2_frame can be attached so as to store the compressed chunks contiguously.
 */
typedef struct blosc2_schunk {
  uint8_t version;
  uint8_t compcode;
  //!< The default compressor. Each chunk can override this.
  uint8_t compcode_meta;
  //!< The default compressor metadata. Each chunk can override this.
  uint8_t clevel;
  //!< The compression level and other compress params.
  uint8_t splitmode;
  //!< The split mode.
  int32_t typesize;
  //!< The type size.
  int32_t blocksize;
  //!< The requested size of the compressed blocks (0; meaning automatic).
  int32_t chunksize;
  //!< Size of each chunk. 0 if not a fixed chunksize.
  uint8_t filters[BLOSC2_MAX_FILTERS];
  //!< The (sequence of) filters.  8-bit per filter.
  uint8_t filters_meta[BLOSC2_MAX_FILTERS];
  //!< Metadata for filters. 8-bit per meta-slot.
  int64_t nchunks;
  //!< Number of chunks in super-chunk.
  int64_t current_nchunk;
  //!< The current chunk that is being accessed
  int64_t nbytes;
  //!< The data size (uncompressed).
  int64_t cbytes;
  //!< The data size + chunks header size (compressed).
  uint8_t** data;
  //!< Pointer to chunk data pointers buffer.
  size_t data_len;
  //!< Length of the chunk data pointers buffer.
  blosc2_storage* storage;
  //!< Pointer to storage info.
  blosc2_frame* frame;
  //!< Pointer to frame used as store for chunks.
  //!<uint8_t* ctx;
  //!< Context for the thread holder. NULL if not acquired.
  blosc2_context* cctx;
  //!< Context for compression
  blosc2_context* dctx;
  //!< Context for decompression.
  struct blosc2_metalayer *metalayers[BLOSC2_MAX_METALAYERS];
  //!< The array of metalayers.
  uint16_t nmetalayers;
  //!< The number of metalayers in the super-chunk
  struct blosc2_metalayer *vlmetalayers[BLOSC2_MAX_VLMETALAYERS];
  //<! The array of variable-length metalayers.
  int16_t nvlmetalayers;
  //!< The number of variable-length metalayers.
  void *tuner_params;
  //!< Tune configuration.
  int tuner_id;
  //<! Id for tuner
  int8_t ndim;
  //<! The ndim (mainly for ZFP usage)
  int64_t *blockshape;
  //<! The blockshape (mainly for ZFP usage)
} blosc2_schunk;


/**
 * @brief Create a new super-chunk.
 *
 * @param storage The storage properties.
 *
 * @remark In case that storage.urlpath is not NULL, the data is stored
 * on-disk.  If the data file(s) exist, they are *overwritten*.
 *
 * @return The new super-chunk.
 */
BLOSC_EXPORT blosc2_schunk* blosc2_schunk_new(blosc2_storage *storage);

/**
 * Create a copy of a super-chunk.
 *
 * @param schunk The super-chunk to be copied.
 * @param storage The storage properties.
 *
 * @return The new super-chunk.
 */
BLOSC_EXPORT blosc2_schunk* blosc2_schunk_copy(blosc2_schunk *schunk, blosc2_storage *storage);

/**
 * @brief Create a super-chunk out of a contiguous frame buffer.
 *
 * @param cframe The buffer of the in-memory frame.
 * @param copy Whether the super-chunk should make a copy of
 * the @p cframe data or not.  The copy will be made to an internal
 * sparse frame.
 *
 * @remark If copy is false, the @p cframe buffer passed will be owned
 * by the super-chunk and will be automatically freed when
 * blosc2_schunk_free() is called.  If the user frees it after the
 * opening, bad things will happen.  Don't do that (or set @p copy).
 *
 * @param len The length of the buffer (in bytes).
 *
 * @return The new super-chunk.
 */
BLOSC_EXPORT blosc2_schunk* blosc2_schunk_from_buffer(uint8_t *cframe, int64_t len, bool copy);

/**
 * @brief Set the private `avoid_cframe_free` field in a frame.
 *
 * @param schunk The super-chunk referencing the frame.
 * @param avoid_cframe_free The value to set in the blosc2_frame_s structure.
 *
 * @warning If you set it to `true` you will be responsible of freeing it.
 */
BLOSC_EXPORT void blosc2_schunk_avoid_cframe_free(blosc2_schunk *schunk, bool avoid_cframe_free);

/**
 * @brief Open an existing super-chunk that is on-disk (frame). No in-memory copy is made.
 *
 * @param urlpath The file name.
 *
 * @return The new super-chunk.  NULL if not found or not in frame format.
 */
BLOSC_EXPORT blosc2_schunk* blosc2_schunk_open(const char* urlpath);

/**
 * @brief Open an existing super-chunk that is on-disk (frame). No in-memory copy is made.
 *
 * @param urlpath The file name.
 *
 * @param offset The frame offset.
 *
 * @return The new super-chunk.  NULL if not found or not in frame format.
 */
BLOSC_EXPORT blosc2_schunk* blosc2_schunk_open_offset(const char* urlpath, int64_t offset);

/**
 * @brief Open an existing super-chunk (no copy is made) using a user-defined I/O interface.
 *
 * @param urlpath The file name.
 *
 * @param udio The user-defined I/O interface.
 *
 * @return The new super-chunk.
 */
BLOSC_EXPORT blosc2_schunk* blosc2_schunk_open_udio(const char* urlpath, const blosc2_io *udio);

/**
 * @brief Open an existing super-chunk (no copy is made) using a user-defined I/O interface.
 *
 * @param urlpath The file name.
 *
 * @param offset The frame offset.
 *
 * @param udio The user-defined I/O interface.
 *
 * @return The new super-chunk.
 */
BLOSC_EXPORT blosc2_schunk* blosc2_schunk_open_offset_udio(const char* urlpath, int64_t offset, const blosc2_io *udio);

/* @brief Convert a super-chunk into a contiguous frame buffer.
 *
 * @param schunk The super-chunk to convert.
 * @param cframe The address of the destination buffer (output).
 * @param needs_free The pointer to a boolean indicating if it is the user's
 * responsibility to free the resulting @p cframe buffer or not.
 *
 * @return If successful, return the size of the (frame) @p cframe buffer.
 * Else, a negative value.
 */
BLOSC_EXPORT int64_t blosc2_schunk_to_buffer(blosc2_schunk* schunk, uint8_t** cframe, bool* needs_free);

/* @brief Store a super-chunk into a file.
 *
 * @param schunk The super-chunk to write.
 * @param urlpath The path for persistent storage.
 *
 * @return If successful, return the size of the (fileframe) in @p urlpath.
 * Else, a negative value.
 */
BLOSC_EXPORT int64_t blosc2_schunk_to_file(blosc2_schunk* schunk, const char* urlpath);

/* @brief Append a super-chunk into a file.
 *
 * @param schunk The super-chunk to write.
 * @param urlpath The path for persistent storage.
 *
 * @return If successful, return the offset where @p schunk has been appended in @p urlpath.
 * Else, a negative value.
 */
BLOSC_EXPORT int64_t blosc2_schunk_append_file(blosc2_schunk* schunk, const char* urlpath);

/**
 * @brief Release resources from a super-chunk.
 *
 * @param schunk The super-chunk to be freed.
 *
 * @remark All the memory resources attached to the super-chunk are freed.
 * If the super-chunk is on-disk, the data continues there for a later
 * re-opening.
 *
 * @return 0 if success.
 */
BLOSC_EXPORT int blosc2_schunk_free(blosc2_schunk *schunk);

/**
 * @brief Append an existing @p chunk to a super-chunk.
 *
 * @param schunk The super-chunk where the chunk will be appended.
 * @param chunk The @p chunk to append.  An internal copy is made, so @p chunk can be reused or
 * freed if desired.
 * @param copy Whether the chunk should be copied internally or can be used as-is.
 *
 * @return The number of chunks in super-chunk. If some problem is
 * detected, this number will be negative.
 */
BLOSC_EXPORT int64_t blosc2_schunk_append_chunk(blosc2_schunk *schunk, uint8_t *chunk, bool copy);

/**
  * @brief Update a chunk at a specific position in a super-chunk.
  *
  * @param schunk The super-chunk where the chunk will be updated.
  * @param nchunk The position where the chunk will be updated.
  * @param chunk The new @p chunk. If an internal copy is made, the @p chunk can be reused or
  * freed if desired.
  * @param copy Whether the chunk should be copied internally or can be used as-is.
  *
  * @return The number of chunks in super-chunk. If some problem is
  * detected, this number will be negative.
  */
BLOSC_EXPORT int64_t blosc2_schunk_update_chunk(blosc2_schunk *schunk, int64_t nchunk, uint8_t *chunk, bool copy);

/**
 * @brief Insert a chunk at a specific position in a super-chunk.
 *
 * @param schunk The super-chunk where the chunk will be appended.
 * @param nchunk The position where the chunk will be inserted.
 * @param chunk The @p chunk to insert. If an internal copy is made, the @p chunk can be reused or
 * freed if desired.
 * @param copy Whether the chunk should be copied internally or can be used as-is.
 *
 * @return The number of chunks in super-chunk. If some problem is
 * detected, this number will be negative.
 */
BLOSC_EXPORT int64_t blosc2_schunk_insert_chunk(blosc2_schunk *schunk, int64_t nchunk, uint8_t *chunk, bool copy);

/**
 * @brief Delete a chunk at a specific position in a super-chunk.
 *
 * @param schunk The super-chunk where the chunk will be deleted.
 * @param nchunk The position where the chunk will be deleted.
 *
 * @return The number of chunks in super-chunk. If some problem is
 * detected, this number will be negative.
 */
BLOSC_EXPORT int64_t blosc2_schunk_delete_chunk(blosc2_schunk *schunk, int64_t nchunk);

/**
 * @brief Append a @p src data buffer to a super-chunk.
 *
 * @param schunk The super-chunk where data will be appended.
 * @param src The buffer of data to compress.
 * @param nbytes The size of the @p src buffer.
 *
 * @return The number of chunks in super-chunk. If some problem is
 * detected, this number will be negative.
 */
BLOSC_EXPORT int64_t blosc2_schunk_append_buffer(blosc2_schunk *schunk, const void *src, int32_t nbytes);

/**
 * @brief Decompress and return the @p nchunk chunk of a super-chunk.
 *
 * If the chunk is uncompressed successfully, it is put in the @p *dest
 * pointer.
 *
 * @param schunk The super-chunk from where the chunk will be decompressed.
 * @param nchunk The chunk to be decompressed (0 indexed).
 * @param dest The buffer where the decompressed data will be put.
 * @param nbytes The size of the area pointed by @p *dest.
 *
 * @warning You must make sure that you have enough space to store the
 * uncompressed data.
 *
 * @return The size of the decompressed chunk or 0 if it is non-initialized. If some problem is
 * detected, a negative code is returned instead.
 */
BLOSC_EXPORT int blosc2_schunk_decompress_chunk(blosc2_schunk *schunk, int64_t nchunk, void *dest, int32_t nbytes);

/**
 * @brief Return a compressed chunk that is part of a super-chunk in the @p chunk parameter.
 *
 * @param schunk The super-chunk from where to extract a chunk.
 * @param nchunk The chunk to be extracted (0 indexed).
 * @param chunk The pointer to the chunk of compressed data.
 * @param needs_free The pointer to a boolean indicating if it is the user's
 * responsibility to free the chunk returned or not.
 *
 * @warning If the super-chunk is backed by a frame that is disk-based, a buffer is allocated for the
 * (compressed) chunk, and hence a free is needed.
 * You can check whether the chunk requires a free with the @p needs_free parameter.
 * If the chunk does not need a free, it means that a pointer to the location in the super-chunk
 * (or the backing in-memory frame) is returned in the @p chunk parameter.
 *
 * @return The size of the (compressed) chunk or 0 if it is non-initialized. If some problem is
 * detected, a negative code is returned instead.
 */
BLOSC_EXPORT int blosc2_schunk_get_chunk(blosc2_schunk *schunk, int64_t nchunk, uint8_t **chunk,
                                         bool *needs_free);

/**
 * @brief Return a (lazy) compressed chunk that is part of a super-chunk in the @p chunk parameter.
 *
 * @param schunk The super-chunk from where to extract a chunk.
 * @param nchunk The chunk to be extracted (0 indexed).
 * @param chunk The pointer to the (lazy) chunk of compressed data.
 * @param needs_free The pointer to a boolean indicating if it is the user's
 * responsibility to free the chunk returned or not.
 *
 * @note For disk-based frames, a lazy chunk is always returned.
 *
 * @warning Currently, a lazy chunk can only be used by #blosc2_decompress_ctx and #blosc2_getitem_ctx.
 *
 * @warning If the super-chunk is backed by a frame that is disk-based, a buffer is allocated for the
 * (compressed) chunk, and hence a free is needed.
 * You can check whether requires a free with the @p needs_free parameter.
 * If the chunk does not need a free, it means that a pointer to the location in the super-chunk
 * (or the backing in-memory frame) is returned in the @p chunk parameter.  In this case the returned
 * chunk is not lazy.
 *
 * @return The size of the (compressed) chunk or 0 if it is non-initialized. If some problem is
 * detected, a negative code is returned instead.  Note that a lazy chunk is somewhat larger than
 * a regular chunk because of the trailer section (for details see `README_CHUNK_FORMAT.rst`).
 */
BLOSC_EXPORT int blosc2_schunk_get_lazychunk(blosc2_schunk *schunk, int64_t nchunk, uint8_t **chunk,
                                             bool *needs_free);

/**
 * @brief Fill buffer with a schunk slice.
 *
 * @param schunk The super-chunk from where to extract a slice.
 * @param start Index (0-based) where the slice begins.
 * @param stop The first index (0-based) that is not in the selected slice.
 * @param buffer The buffer where the data will be stored.
 *
 * @warning You must make sure that you have enough space in buffer to store the
 * uncompressed data.
 *
 * @return An error code.
 */
BLOSC_EXPORT int blosc2_schunk_get_slice_buffer(blosc2_schunk *schunk, int64_t start, int64_t stop, void *buffer);

/**
 * @brief Update a schunk slice from buffer.
 *
 * @param schunk The super-chunk where to set the slice.
 * @param start Index (0-based) where the slice begins.
 * @param stop The first index (0-based) that is not in the selected slice.
 * @param buffer The buffer containing the data to set.
 *
 *
 * @return An error code.
 */
BLOSC_EXPORT int blosc2_schunk_set_slice_buffer(blosc2_schunk *schunk, int64_t start, int64_t stop, void *buffer);

/**
 * @brief Return the @p cparams associated to a super-chunk.
 *
 * @param schunk The super-chunk from where to extract the compression parameters.
 * @param cparams The pointer where the compression params will be returned.
 *
 * @warning A new struct is allocated, and the user should free it after use.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_schunk_get_cparams(blosc2_schunk *schunk, blosc2_cparams **cparams);

/**
 * @brief Return the @p dparams struct associated to a super-chunk.
 *
 * @param schunk The super-chunk from where to extract the decompression parameters.
 * @param dparams The pointer where the decompression params will be returned.
 *
 * @warning A new struct is allocated, and the user should free it after use.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_schunk_get_dparams(blosc2_schunk *schunk, blosc2_dparams **dparams);

/**
 * @brief Reorder the chunk offsets of an existing super-chunk.
 *
 * @param schunk The super-chunk whose chunk offsets are to be reordered.
 * @param offsets_order The new order of the chunk offsets.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_schunk_reorder_offsets(blosc2_schunk *schunk, int64_t *offsets_order);

/**
 * @brief Get the length (in bytes) of the internal frame of the super-chunk.
 *
 * @param schunk The super-chunk.
 *
 * @return The length (in bytes) of the internal frame.
 * If there is not an internal frame, an estimate of the length is provided.
 */
BLOSC_EXPORT int64_t blosc2_schunk_frame_len(blosc2_schunk* schunk);

/**
 * @brief Quickly fill an empty frame with special values (zeros, NaNs, uninit).
 *
 * @param schunk The super-chunk to be filled.  This must be empty initially.
 * @param nitems The number of items to fill.
 * @param special_value The special value to use for filling.  The only values
 * supported for now are BLOSC2_SPECIAL_ZERO, BLOSC2_SPECIAL_NAN and BLOSC2_SPECIAL_UNINIT.
 * @param chunksize The chunksize for the chunks that are to be added to the super-chunk.
 *
 * @return The total number of chunks that have been added to the super-chunk.
 * If there is an error, a negative value is returned.
 */
BLOSC_EXPORT int64_t blosc2_schunk_fill_special(blosc2_schunk* schunk, int64_t nitems,
                                                int special_value, int32_t chunksize);


/*********************************************************************
  Functions related with fixed-length metalayers.
*********************************************************************/

/**
 * @brief Find whether the schunk has a metalayer or not.
 *
 * @param schunk The super-chunk from which the metalayer will be checked.
 * @param name The name of the metalayer to be checked.
 *
 * @return If successful, return the index of the metalayer. Else, return a negative value.
 */
static inline int blosc2_meta_exists(blosc2_schunk *schunk, const char *name) {
  if (strlen(name) > BLOSC2_METALAYER_NAME_MAXLEN) {
    BLOSC_TRACE_ERROR("Metalayers cannot be larger than %d chars.", BLOSC2_METALAYER_NAME_MAXLEN);
    return BLOSC2_ERROR_INVALID_PARAM;
  }

  if (schunk == NULL) {
    BLOSC_TRACE_ERROR("Schunk must not be NUll.");
    return BLOSC2_ERROR_INVALID_PARAM;
  }

  for (int nmetalayer = 0; nmetalayer < schunk->nmetalayers; nmetalayer++) {
    if (strcmp(name, schunk->metalayers[nmetalayer]->name) == 0) {
      return nmetalayer;
    }
  }
  return BLOSC2_ERROR_NOT_FOUND;
}

/**
 * @brief Add content into a new metalayer.
 *
 * @param schunk The super-chunk to which the metalayer should be added.
 * @param name The name of the metalayer.
 * @param content The content of the metalayer.
 * @param content_len The length of the content.
 *
 * @return If successful, the index of the new metalayer. Else, return a negative value.
 */
BLOSC_EXPORT int blosc2_meta_add(blosc2_schunk *schunk, const char *name, uint8_t *content,
                                 int32_t content_len);

/**
 * @brief Update the content of an existing metalayer.
 *
 * @param schunk The frame containing the metalayer.
 * @param name The name of the metalayer to be updated.
 * @param content The new content of the metalayer.
 * @param content_len The length of the content.
 *
 * @note Contrarily to #blosc2_meta_add the updates to metalayers
 * are automatically serialized into a possible attached frame.
 *
 * @return If successful, the index of the metalayer. Else, return a negative value.
 */
BLOSC_EXPORT int blosc2_meta_update(blosc2_schunk *schunk, const char *name, uint8_t *content,
                                    int32_t content_len);

/**
 * @brief Get the content out of a metalayer.
 *
 * @param schunk The frame containing the metalayer.
 * @param name The name of the metalayer.
 * @param content The pointer where the content will be put.
 * @param content_len The length of the content.
 *
 * @warning The @p **content receives a malloc'ed copy of the content.
 * The user is responsible of freeing it.
 *
 * @note This function is inlined and available even when not linking with libblosc2.
 *
 * @return If successful, the index of the new metalayer. Else, return a negative value.
 */
static inline int blosc2_meta_get(blosc2_schunk *schunk, const char *name, uint8_t **content,
                                  int32_t *content_len) {
  int nmetalayer = blosc2_meta_exists(schunk, name);
  if (nmetalayer < 0) {
    BLOSC_TRACE_WARNING("Metalayer \"%s\" not found.", name);
    return nmetalayer;
  }
  *content_len = schunk->metalayers[nmetalayer]->content_len;
  *content = (uint8_t*)malloc((size_t)*content_len);
  memcpy(*content, schunk->metalayers[nmetalayer]->content, (size_t)*content_len);
  return nmetalayer;
}


/*********************************************************************
  Variable-length metalayers functions.
*********************************************************************/

/**
 * @brief Find whether the schunk has a variable-length metalayer or not.
 *
 * @param schunk The super-chunk from which the variable-length metalayer will be checked.
 * @param name The name of the variable-length metalayer to be checked.
 *
 * @return If successful, return the index of the variable-length metalayer. Else, return a negative value.
 */
BLOSC_EXPORT int blosc2_vlmeta_exists(blosc2_schunk *schunk, const char *name);

/**
 * @brief Add content into a new variable-length metalayer.
 *
 * @param schunk The super-chunk to which the variable-length metalayer should be added.
 * @param name The name of the variable-length metalayer.
 * @param content The content to be added.
 * @param content_len The length of the content.
 * @param cparams The parameters for compressing the variable-length metalayer content. If NULL,
 * the `BLOSC2_CPARAMS_DEFAULTS` will be used.
 *
 * @return If successful, the index of the new variable-length metalayer. Else, return a negative value.
 */
BLOSC_EXPORT int blosc2_vlmeta_add(blosc2_schunk *schunk, const char *name,
                                   uint8_t *content, int32_t content_len,
                                   blosc2_cparams *cparams);

/**
 * @brief Update the content of an existing variable-length metalayer.
 *
 * @param schunk The super-chunk containing the variable-length metalayer.
 * @param name The name of the variable-length metalayer to be updated.
 * @param content The new content of the variable-length metalayer.
 * @param content_len The length of the content.
 * @param cparams The parameters for compressing the variable-length metalayer content. If NULL,
 * the `BLOSC2_CPARAMS_DEFAULTS` will be used.
 *
 * @return If successful, the index of the variable-length metalayer. Else, return a negative value.
 */
BLOSC_EXPORT int blosc2_vlmeta_update(blosc2_schunk *schunk, const char *name,
                                      uint8_t *content, int32_t content_len,
                                      blosc2_cparams *cparams);

/**
 * @brief Get the content out of a variable-length metalayer.
 *
 * @param schunk The super-chunk containing the variable-length metalayer.
 * @param name The name of the variable-length metalayer.
 * @param content The pointer where the content will be put.
 * @param content_len The pointer where the length of the content will be put.
 *
 * @warning The @p **content receives a malloc'ed copy of the content.
 * The user is responsible of freeing it.
 *
 * @return If successful, the index of the new variable-length metalayer. Else, return a negative value.
 */
BLOSC_EXPORT int blosc2_vlmeta_get(blosc2_schunk *schunk, const char *name,
                                   uint8_t **content, int32_t *content_len);

/**
 * @brief Delete the variable-length metalayer from the super-chunk.
 *
 * @param schunk The super-chunk containing the variable-length metalayer.
 * @param name The name of the variable-length metalayer.
 *
 * @return If successful, the number of the variable-length metalayers in the super-chunk. Else, return a negative value.
 */
BLOSC_EXPORT int blosc2_vlmeta_delete(blosc2_schunk *schunk, const char *name);

/**
 * @brief Get a list of all the variable-length metalayer names.
 *
 * @param schunk The super-chunk containing the variable-length metalayers.
 * @param names The pointer to a char** to store the name pointers. This should
 * be of size *schunk->nvlmetalayers * sizeof(char*).
 *
 * @return The number of the variable-length metalayers in the super-chunk.
 * This cannot fail unless the user does not pass a @p names which is large enough to
 * keep pointers to all names, in which case funny things (seg faults and such) will happen.
 */
BLOSC_EXPORT int blosc2_vlmeta_get_names(blosc2_schunk *schunk, char **names);


/*********************************************************************
  Time measurement utilities.
*********************************************************************/

#if defined(_WIN32)
/* For QueryPerformanceCounter(), etc. */
  #include <windows.h>
#elif defined(__MACH__)
#include <mach/clock.h>
#include <mach/mach.h>
#include <time.h>
#elif defined(__unix__)
#if defined(__linux__)
    #include <time.h>
  #else
    #include <sys/time.h>
  #endif
#else
  #error Unable to detect platform.
#endif

/* The type of timestamp used on this system. */
#if defined(_WIN32)
typedef LARGE_INTEGER blosc_timestamp_t;
#else
typedef struct timespec blosc_timestamp_t;
#endif

/*
 * @brief Set a timestamp.
 *
 * @param timestamp
 *
 */
BLOSC_EXPORT void blosc_set_timestamp(blosc_timestamp_t* timestamp);

/*
 * @brief Get the nanoseconds between 2 timestamps.
 *
 * @param start_time
 * @param end_time
 *
 * @return The nanoseconds between start_time and end_time.
 */
BLOSC_EXPORT double blosc_elapsed_nsecs(blosc_timestamp_t start_time,
                                        blosc_timestamp_t end_time);

/*
 * @brief Get the seconds between 2 timestamps.
 *
 * @param start_time
 * @param end_time
 *
 * @return The seconds between start_time and end_time.
 */
BLOSC_EXPORT double blosc_elapsed_secs(blosc_timestamp_t start_time,
                                       blosc_timestamp_t end_time);


/*********************************************************************
  Low-level functions follows.  Use them only if you are an expert!
*********************************************************************/

/**
 * @brief Get the internal blocksize to be used during compression. 0 means
 * that an automatic blocksize is computed internally.
 *
 * @return The size in bytes of the internal block size.
 */
BLOSC_EXPORT int blosc1_get_blocksize(void);

/**
 * @brief Force the use of a specific blocksize. If 0, an automatic
 * blocksize will be used (the default).
 *
 * @warning The blocksize is a critical parameter with important
 * restrictions in the allowed values, so use this with care.
 */
BLOSC_EXPORT void blosc1_set_blocksize(size_t blocksize);


/**
  * @brief Set the split mode.

  * @param splitmode It can take the next values:
  *  BLOSC_FORWARD_COMPAT_SPLIT
  *  BLOSC_AUTO_SPLIT
  *  BLOSC_NEVER_SPLIT
  *  BLOSC_ALWAYS_SPLIT
  *
  * BLOSC_FORWARD_COMPAT offers reasonably forward compatibility,
  * BLOSC_AUTO_SPLIT is for nearly optimal results (based on heuristics),
  * BLOSC_NEVER_SPLIT and BLOSC_ALWAYS_SPLIT are for the user experimenting
  *  when trying to get best compression ratios and/or speed.
  *
  * If not called, the default mode is BLOSC_FORWARD_COMPAT_SPLIT.
  *
  * This function should always succeed.
 */
BLOSC_EXPORT void blosc1_set_splitmode(int splitmode);


/**
 * @brief Get the offsets of a frame in a super-chunk.
 *
 * @param schunk The super-chunk containing the frame.
 *
 * @return If successful, return a pointer to a buffer of the decompressed offsets.
 * The number of offsets is equal to schunk->nchunks; the user is
 * responsible to free this buffer. Else, return a NULL value.
 */
BLOSC_EXPORT int64_t* blosc2_frame_get_offsets(blosc2_schunk *schunk);


/*********************************************************************
  Structures and functions related with compression codecs.
*********************************************************************/

typedef int (* blosc2_codec_encoder_cb) (const uint8_t *input, int32_t input_len, uint8_t *output, int32_t output_len,
            uint8_t meta, blosc2_cparams *cparams, const void* chunk);
typedef int (* blosc2_codec_decoder_cb) (const uint8_t *input, int32_t input_len, uint8_t *output, int32_t output_len,
            uint8_t meta, blosc2_dparams *dparams, const void* chunk);

typedef struct {
  uint8_t compcode;
  //!< The codec identifier.
  char *compname;
  //!< The codec name.
  uint8_t complib;
  //!< The codec library format.
  uint8_t version;
  //!< The codec version.
  blosc2_codec_encoder_cb encoder;
  //!< The codec encoder that is used during compression.
  blosc2_codec_decoder_cb decoder;
  //!< The codec decoder that is used during decompression.
} blosc2_codec;

/**
 * @brief Register locally a user-defined codec in Blosc.
 *
 * @param codec The codec to register.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_register_codec(blosc2_codec *codec);


/*********************************************************************
  Structures and functions related with filters plugins.
*********************************************************************/

typedef int (* blosc2_filter_forward_cb)  (const uint8_t *, uint8_t *, int32_t, uint8_t, blosc2_cparams *,
                                           uint8_t);
typedef int (* blosc2_filter_backward_cb) (const uint8_t *, uint8_t *, int32_t, uint8_t, blosc2_dparams *,
                                           uint8_t);

/**
 * @brief The parameters for a user-defined filter.
 */
typedef struct {
  uint8_t id;
  //!< The filter identifier.
  char * name;
  //!< The filter name.
  uint8_t version;
  //!< The filter version.
  blosc2_filter_forward_cb forward;
  //!< The filter function that is used during compression.
  blosc2_filter_backward_cb backward;
  //!< The filter function that is used during decompression.
} blosc2_filter;

/**
 * @brief Register locally a user-defined filter in Blosc.
 *
 * @param filter The filter to register.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_register_filter(blosc2_filter *filter);

/*********************************************************************
  Directory utilities.
*********************************************************************/

/*
 * @brief Remove a directory and its files.
 *
 * @param path The directory to remove.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_remove_dir(const char *path);

/*
 * @brief Remove a file or a directory given by path.
 *
 * @param path The file or directory to remove.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_remove_urlpath(const char *path);

/*
 * @brief Rename a file or a directory given by old_urlpath to new_path.
 *
 * @param old_urlpath The original path to the directory or file.
 * @param new_path The new path to the directory or file.
 *
 * @return 0 if succeeds. Else a negative code is returned.
 */
BLOSC_EXPORT int blosc2_rename_urlpath(char* old_urlpath, char* new_path);


/*********************************************************************
  Index utilities.
*********************************************************************/

/*
 * @brief Convert a sequential index into a multidimensional index
 */
BLOSC_EXPORT void blosc2_unidim_to_multidim(uint8_t ndim, int64_t *shape, int64_t i, int64_t *index);

/*
 * @brief Convert a multidimensional index into a sequential index
 */
BLOSC_EXPORT void blosc2_multidim_to_unidim(const int64_t *index, int8_t ndim, const int64_t *strides, int64_t *i);

/*
 * @brief Get the unidimensional chunk indexes needed to get a slice of a schunk or a b2nd array
 *
 * @param schunk The super-chunk (of b2nd array or not).
 * @param start Index (0-based if it is a schunk) where the slice begins.
 * @param stop The first index (0-based if it is a schunk) that is not in the selected slice.
 * @param chunks_idx The pointer to the buffer where the indexes will be written. It is the user responsibility
 * to free the buffer.
 *
 * @return The number of chunks needed to get the slice. If some problem is
 * detected, a negative code is returned instead.
 */
BLOSC_EXPORT int blosc2_get_slice_nchunks(blosc2_schunk* schunk, int64_t *start, int64_t *stop, int64_t **chunks_idx);


/*********************************************************************
  Private functions, these are here for convenience,
  and are not meant to be included in public docs
*********************************************************************/

// Private function needed in b2nd.h for deserializing meta
static inline void swap_store(void *dest, const void *pa, int size) {
  uint8_t *pa_ = (uint8_t *) pa;
  uint8_t *pa2_ = (uint8_t*)malloc((size_t) size);
  int i = 1; /* for big/little endian detection */
  char *p = (char *) &i;

  if (p[0] == 1) {
    /* little endian */
    switch (size) {
      case 8:
        pa2_[0] = pa_[7];
        pa2_[1] = pa_[6];
        pa2_[2] = pa_[5];
        pa2_[3] = pa_[4];
        pa2_[4] = pa_[3];
        pa2_[5] = pa_[2];
        pa2_[6] = pa_[1];
        pa2_[7] = pa_[0];
        break;
      case 4:
        pa2_[0] = pa_[3];
        pa2_[1] = pa_[2];
        pa2_[2] = pa_[1];
        pa2_[3] = pa_[0];
        break;
      case 2:
        pa2_[0] = pa_[1];
        pa2_[1] = pa_[0];
        break;
      case 1:
        pa2_[0] = pa_[0];
        break;
      default:
        fprintf(stderr, "Unhandled nitems: %d\n", size);
    }
  }
  memcpy(dest, pa2_, size);
  free(pa2_);
}

#ifdef __cplusplus
}
#endif

#endif /* BLOSC_BLOSC2_H */
