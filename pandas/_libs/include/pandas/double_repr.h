/*
Copyright (c) 2026, PyData Development Team
All rights reserved.
Distributed under the terms of the BSD Simplified License.
*/
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/*
Shortest round-trip float formatting matching numpy scalar str().

Notation choice replicates numpy (value-based, not digit-exponent-based):
fixed iff 1e-4 <= |val| < 1e16 for float64. For float32 the upper cutoff is
supplied by the caller (fixed_hi) because numpy changed it in 2.3.0 (1e16 ->
1e6). Compared in double precision. Scientific style matches Python/numpy:
"d[.digits]e[+-]XX" with a sign and at least two exponent digits.
*/

/* Enough for sign + 17 digits + dot + zeros/exponent, with headroom. */
#define PD_DOUBLE_REPR_BUFSIZE 40

/*
Write the repr of val to out (not NUL-terminated), returning the length
written, or -1 on error. out must have room for PD_DOUBLE_REPR_BUFSIZE bytes.
When floating-point std::to_chars is unavailable (PD_HAVE_FP_TO_CHARS not
defined) this falls back to PyOS_double_to_string and requires the GIL.
*/
int pd_double_repr(double val, char *out);

/*
float32 analogue of pd_double_repr; only usable when
pd_float32_repr_available() returns nonzero (requires std::to_chars).
fixed_hi is the value-based fixed/scientific cutoff (1e16 for numpy < 2.3,
1e6 for numpy >= 2.3).
*/
int pd_float32_repr(float val, double fixed_hi, char *out);
int pd_float32_repr_available(void);

#ifdef __cplusplus
}
#endif
