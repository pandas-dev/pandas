#include <string.h>

#include "pandas/bitmask_algorithms.h"

static const uint8_t clear_mask[8] = {0x0, 0x1,  0x3,  0x7,
                                      0xf, 0x1f, 0x3f, 0x7f};

void ConcatenateBitmapData(const struct ArrowBitmap **bitmaps, size_t nbitmaps,
                           struct ArrowBitmap *out) {
  if (nbitmaps == 0) {
    return;
  }

  int64_t bits_processed = 0;
  uint8_t *out_cursor = out->buffer.data;
  size_t start_bit_pos = 0;
  for (size_t i = 0; i < nbitmaps; i++) {
    const struct ArrowBitmap *bitmap = bitmaps[i];
    const int64_t nbytes = bitmap->buffer.size_bytes;
    const size_t trailing_nbits = bitmap->size_bits % 8;

    // As we loop through each array, any time we end up starting
    // on a word boundary we can simply use memcpy. If we are not
    // so lucky we fall back to bit shifting each element
    if (start_bit_pos == 0) {
      memcpy(out_cursor, bitmap->buffer.data, nbytes);
    } else {
      for (size_t j = 0; j < nbytes; j++) {
        const uint8_t lshifted = bitmap->buffer.data[j] << start_bit_pos;
        out_cursor[j] = (out_cursor[j] & clear_mask[start_bit_pos]) | lshifted;

        const uint8_t rshifted = bitmap->buffer.data[j] >> (8 - start_bit_pos);
        out_cursor[j + 1] = rshifted;
      }
    }

    start_bit_pos = (start_bit_pos + trailing_nbits) % 8;
    if (start_bit_pos == 0) {
      out_cursor += nbytes;
    } else {
      out_cursor += nbytes - 1;
    }

    bits_processed += bitmap->size_bits;
  }

  out->size_bits = bits_processed;
  out->buffer.size_bytes = bits_processed / 8;
  if ((bits_processed % 8) > 0) {
    out->buffer.size_bytes += 1;
  }
}

bool BitmapAny(const struct ArrowBitmap *bitmap) {
  const size_t nbits = bitmap->size_bits;
  const size_t size_bytes = bitmap->buffer.size_bytes;
  const uint8_t *buf = bitmap->buffer.data;

  if (nbits < 1) {
    return false;
  }

  for (size_t i = 0; i < size_bytes - 1; i++) {
    if (buf[i] > 0) {
      return true;
    }
  }

  const size_t bits_remaining = nbits % 8;
  for (size_t i = 0; i < bits_remaining; i++) {
    if (ArrowBitGet(buf, nbits - i - 1)) {
      return true;
    }
  }

  return false;
}

bool BitmapAll(const struct ArrowBitmap *bitmap) {
  const size_t nbits = bitmap->size_bits;
  const size_t size_bytes = bitmap->buffer.size_bytes;
  const uint8_t *buf = bitmap->buffer.data;

  if (nbits < 1) {
    return true;
  }

  for (size_t i = 0; i < size_bytes - 1; i++) {
    if (buf[i] != 0xff) {
      return false;
    }
  }

  const size_t bits_remaining = nbits % 8;
  for (size_t i = 0; i < bits_remaining; i++) {
    if (ArrowBitGet(buf, nbits - i - 1) == 0) {
      return false;
    }
  }

  return true;
}

int BitmapOr(const struct ArrowBitmap *bitmap1,
             const struct ArrowBitmap *bitmap2, struct ArrowBitmap *out) {
  if (bitmap1->size_bits != bitmap2->size_bits) {
    return -1;
  } else if (!(out->buffer.capacity_bytes >= bitmap1->buffer.size_bytes)) {
    return -1;
  }

  for (size_t i = 0; i < bitmap1->buffer.size_bytes; i++) {
    out->buffer.data[i] = bitmap1->buffer.data[i] | bitmap2->buffer.data[i];
  }

  out->size_bits = bitmap1->size_bits;
  out->buffer.size_bytes = bitmap1->buffer.size_bytes;

  return 0;
}

int BitmapAnd(const struct ArrowBitmap *bitmap1,
              const struct ArrowBitmap *bitmap2, struct ArrowBitmap *out) {
  if (bitmap1->size_bits != bitmap2->size_bits) {
    return -1;
  } else if (!(out->buffer.capacity_bytes >= bitmap1->buffer.size_bytes)) {
    return -1;
  }

  for (size_t i = 0; i < bitmap1->buffer.size_bytes; i++) {
    out->buffer.data[i] = bitmap1->buffer.data[i] & bitmap2->buffer.data[i];
  }

  out->size_bits = bitmap1->size_bits;
  out->buffer.size_bytes = bitmap1->buffer.size_bytes;

  return 0;
}

int BitmapXor(const struct ArrowBitmap *bitmap1,
              const struct ArrowBitmap *bitmap2, struct ArrowBitmap *out) {
  if (bitmap1->size_bits != bitmap2->size_bits) {
    return -1;
  } else if (!(out->buffer.capacity_bytes >= bitmap1->buffer.size_bytes)) {
    return -1;
  }

  for (size_t i = 0; i < bitmap1->buffer.size_bytes; i++) {
    out->buffer.data[i] = bitmap1->buffer.data[i] ^ bitmap2->buffer.data[i];
  }

  out->size_bits = bitmap1->size_bits;
  out->buffer.size_bytes = bitmap1->buffer.size_bytes;

  return 0;
}

int BitmapInvert(const struct ArrowBitmap *bitmap, struct ArrowBitmap *out) {
  if (!(out->buffer.capacity_bytes >= bitmap->buffer.size_bytes)) {
    return -1;
  }

  for (size_t i = 0; i < bitmap->buffer.size_bytes; i++) {
    out->buffer.data[i] = ~bitmap->buffer.data[i];
  }

  out->size_bits = bitmap->size_bits;
  out->buffer.size_bytes = bitmap->buffer.size_bytes;

  return 0;
}

int BitmapTake(const struct ArrowBitmap *bitmap, const int64_t *indices,
               size_t nindices, struct ArrowBitmap *out) {
  int64_t bytes_needed = nindices / 8;
  if ((nindices % 8) > 0) {
    bytes_needed += 1;
  }

  if (!(out->buffer.capacity_bytes >= bytes_needed)) {
    return -1;
  }

  for (size_t i = 0; i < nindices; i++) {
    int64_t index = indices[i];
    if (index < 0) {
      return -1;
    }

    int8_t value = ArrowBitGet(bitmap->buffer.data, index);
    ArrowBitmapAppendUnsafe(out, value, 1);
  }

  return 0;
}

int BitmapPutFromBufferMask(struct ArrowBitmap *bitmap, const uint8_t *buf,
                            size_t n, uint8_t value) {
  int64_t bytes_needed = n / 8;
  if ((n % 8) > 0) {
    bytes_needed += 1;
  }

  if (bytes_needed > bitmap->buffer.capacity_bytes) {
    return -1;
  }

  for (size_t i = 0; i < n; i++) {
    if (buf[i]) {
      ArrowBitSetTo(bitmap->buffer.data, i, value);
    }
  }

  return 0;
}
