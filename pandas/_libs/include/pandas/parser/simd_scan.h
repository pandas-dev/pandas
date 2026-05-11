/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Minimum bytes the scanner can process in one call. Callers should
// fall through to the scalar path when fewer bytes remain.
#define PD_SCAN_MIN_BYTES 16

typedef struct pd_scanner pd_scanner;

// Build a scanner that halts on any of `n` special bytes. Supported
// values for `n` are 2 (quoted-field scan) and 6 (unquoted-field scan).
// Returns NULL on allocation failure or unsupported `n`.
pd_scanner *pd_scanner_create(const char *chars, int n);

// Free a scanner. Accepts NULL.
void pd_scanner_destroy(pd_scanner *scanner);

// Returns the byte offset of the first special char in data[0..len),
// or `len` if no special char was found within full SIMD chunks. The
// trailing <PD_SCAN_MIN_BYTES bytes are not scanned; the caller's
// scalar fallback handles them.
size_t pd_scanner_scan(const pd_scanner *scanner, const char *data, size_t len);

#ifdef __cplusplus
}
#endif
