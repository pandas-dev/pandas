#include <string.h>

#include "pandas/bitmask_algorithms.h"

static const uint8_t clear_mask[8] = {0x0, 0x1,  0x3,  0x7,
                                      0xf, 0x1f, 0x3f, 0x7f};

void ConcatenateBitmapData(struct ArrowBitmap *bitmaps, size_t nbitmaps,
                           uint8_t *out) {
  if (nbitmaps == 0) {
    return;
  }

  uint8_t *out_cursor = out;
  // As we loop through each array, any time we end up starting
  // on a word boundary we can simply use memcpy. If we are not
  // so lucky we fall back to bit shifting each element
  size_t start_bit_pos = 0;
  for (size_t i = 0; i < nbitmaps; i++) {
    struct ArrowBitmap bitmap = bitmaps[i];
    int64_t nbytes = bitmap.buffer.size_bytes;
    size_t trailing_nbits = bitmap.size_bits % 8;

    if (start_bit_pos == 0) {
      memcpy(out_cursor, bitmap.buffer.data, nbytes);
    } else {
      for (size_t j = 0; j < nbytes; j++) {
        uint8_t lshifted = bitmap.buffer.data[j] << start_bit_pos;
        out_cursor[j] = (out_cursor[j] & clear_mask[start_bit_pos]) | lshifted;

        uint8_t rshifted = bitmap.buffer.data[j] >> (8 - start_bit_pos);
        out_cursor[j + 1] = rshifted;
      }
    }

    start_bit_pos = (start_bit_pos + trailing_nbits) % 8;
    if (start_bit_pos == 0) {
      out_cursor += nbytes;
    } else {
      out_cursor += nbytes - 1;
    }
  }
}
