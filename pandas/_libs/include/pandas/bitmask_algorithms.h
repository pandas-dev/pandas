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

bool BitmapAny(const struct ArrowBitmap* bitmap);
