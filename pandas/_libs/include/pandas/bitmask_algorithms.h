#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "pandas/vendored/nanoarrow.h"

/*
  Concatenates the data from an array of bitmaps with size nbitmaps
  into a buffer "out". Order is preserved and out is assumed to have
  enough bytes to hold all elements.
*/
void ConcatenateBitmapData(const struct ArrowBitmap **bitmaps, size_t nbitmaps,
                           uint8_t *out);

bool BitmapAny(const struct ArrowBitmap *bitmap);
bool BitmapAll(const struct ArrowBitmap *bitmap);

/* Returns -1 on failure. On success returns 0 and writes to out */
int BitmapOr(const struct ArrowBitmap *bitmap1,
             const struct ArrowBitmap *bitmap2, struct ArrowBitmap *out);

/* Returns -1 on failure. On success returns 0 and writes to out */
int BitmapXor(const struct ArrowBitmap *bitmap1,
              const struct ArrowBitmap *bitmap2, struct ArrowBitmap *out);

/* Returns -1 on failure. On success returns 0 and writes to out */
int BitmapAnd(const struct ArrowBitmap *bitmap1,
              const struct ArrowBitmap *bitmap2, struct ArrowBitmap *out);

/* Returns -1 on failure. On success returns 0 and writes to out */
int BitmapInvert(const struct ArrowBitmap *bitmap, struct ArrowBitmap *out);

/* Returns -1 on failure. On success returns 0 and writes to out */
int BitmapTake(const struct ArrowBitmap *bitmap, const int64_t *indices,
               size_t nindices, struct ArrowBitmap *out);

/* Returns -1 on failure. On success returns 0 and writes to out */
int BitmapPutFromBufferMask(struct ArrowBitmap *bitmap, const uint8_t *buf,
                            size_t n, uint8_t value);
