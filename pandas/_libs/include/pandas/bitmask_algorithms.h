#include <stddef.h>
#include <stdint.h>

#include "pandas/vendored/nanoarrow.h"

/*
  Ordered concatenation of bitmasks. Masks is the data itself,
  nmasks is the number of masks to concatenate, mask_nbits is the
  number of bits within each mask to concatenate.

  Concatenation preserves order.

  out is assumed to have enough bytes to hold all elements.
*/
void ConcatenateBitmapData(struct ArrowBitmap *bitmaps, size_t nbitmaps,
                           uint8_t *out);
