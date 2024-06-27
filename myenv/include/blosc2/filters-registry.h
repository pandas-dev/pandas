/*********************************************************************
  Blosc - Blocked Shuffling and Compression Library

  Copyright (c) 2021  Blosc Development Team <blosc@blosc.org>
  https://blosc.org
  License: BSD 3-Clause (see LICENSE.txt)

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

#ifndef BLOSC_BLOSC2_FILTERS_REGISTRY_H
#define BLOSC_BLOSC2_FILTERS_REGISTRY_H

#ifdef __cplusplus
extern "C" {
#endif

enum {
    BLOSC_FILTER_NDCELL = 32,
    //!< Simple filter for grouping NDim cell data together.
    //!< See https://github.com/Blosc/c-blosc2/blob/main/plugins/filters/ndcell/README.md
    BLOSC_FILTER_NDMEAN = 33,
    //!< Simple filter for replacing content of a NDim cell with its mean value.
    //!< See https://github.com/Blosc/c-blosc2/blob/main/plugins/filters/ndmean/README.md
    BLOSC_FILTER_BYTEDELTA_BUGGY = 34,
    // buggy version. See #524
    BLOSC_FILTER_BYTEDELTA = 35,
    //!< Byteshuffle + delta.  The typesize should be specified in the `filters_meta` slot.
    //!<  Sometimes this can represent an advantage over
    //!< @ref BLOSC_SHUFFLE or @ref BLOSC_BITSHUFFLE.
    //!< See https://www.blosc.org/posts/bytedelta-enhance-compression-toolset/
    BLOSC_FILTER_INT_TRUNC = 36,
    //!< Truncate int precision; positive values in `filters_meta` slot will keep bits;
    //!< negative values will remove (set to zero) bits.
    //!< This is similar to @ref BLOSC_TRUNC_PREC, but for integers instead of floating point data.
};

void register_filters(void);

// For dynamically loaded filters
typedef struct {
    char *forward;
    char *backward;
} filter_info;

#ifdef __cplusplus
}
#endif

#endif /* BLOSC_BLOSC2_FILTERS_REGISTRY_H */
