/*
Copyright (c) 2026, PyData Development Team
All rights reserved.
Distributed under the terms of the BSD Simplified License.

Shortest round-trip float formatting matching numpy scalar str(): shortest
digits (via std::to_chars, or PyOS_double_to_string as fallback) laid out
with numpy's value-based fixed/scientific cutoffs.
*/

#include "pandas/double_repr.h"

#include <cmath>
#include <cstring>

#ifdef PD_HAVE_FP_TO_CHARS
#  include <charconv>
#  include <system_error>
#else
#  define PY_SSIZE_T_CLEAN
#  include <Python.h>
#endif

namespace {

// Handle nan/inf/zero; returns length written or 0 if val needs digits.
int emit_special(double val, char *out) {
  if (std::isnan(val)) {
    std::memcpy(out, "nan", 3);
    return 3;
  }
  char *p = out;
  if (std::signbit(val)) {
    *p++ = '-';
  }
  if (std::isinf(val)) {
    std::memcpy(p, "inf", 3);
    return static_cast<int>(p - out) + 3;
  }
  if (val == 0.0) {
    std::memcpy(p, "0.0", 3);
    return static_cast<int>(p - out) + 3;
  }
  return 0;
}

// Lay out shortest-round-trip digits (value = d.igits x 10^exp10) in numpy
// scalar-str notation.
int emit_digits(bool neg, const char *digits, int ndigits, int exp10,
                bool use_fixed, char *out) {
  char *p = out;
  if (neg) {
    *p++ = '-';
  }
  if (use_fixed) {
    if (exp10 >= 0) {
      int int_digits = exp10 + 1;
      if (int_digits >= ndigits) {
        std::memcpy(p, digits, ndigits);
        p += ndigits;
        for (int i = ndigits; i < int_digits; i++) {
          *p++ = '0';
        }
        *p++ = '.';
        *p++ = '0';
      } else {
        std::memcpy(p, digits, int_digits);
        p += int_digits;
        *p++ = '.';
        std::memcpy(p, digits + int_digits, ndigits - int_digits);
        p += ndigits - int_digits;
      }
    } else {
      *p++ = '0';
      *p++ = '.';
      for (int i = 0; i < -exp10 - 1; i++) {
        *p++ = '0';
      }
      std::memcpy(p, digits, ndigits);
      p += ndigits;
    }
  } else {
    *p++ = digits[0];
    if (ndigits > 1) {
      *p++ = '.';
      std::memcpy(p, digits + 1, ndigits - 1);
      p += ndigits - 1;
    }
    *p++ = 'e';
    int aexp = exp10;
    if (exp10 < 0) {
      *p++ = '-';
      aexp = -exp10;
    } else {
      *p++ = '+';
    }
    char ebuf[8];
    int elen = 0;
    do {
      ebuf[elen++] = static_cast<char>('0' + aexp % 10);
      aexp /= 10;
    } while (aexp != 0);
    while (elen < 2) {
      ebuf[elen++] = '0';
    }
    while (elen > 0) {
      *p++ = ebuf[--elen];
    }
  }
  return static_cast<int>(p - out);
}

#ifdef PD_HAVE_FP_TO_CHARS

// Format via shortest-scientific to_chars, then re-layout.
template <typename T>
int repr_to_chars(T val, double fixed_lo, double fixed_hi, char *out) {
  int res = emit_special(static_cast<double>(val), out);
  if (res != 0) {
    return res;
  }
  char sci[PD_DOUBLE_REPR_BUFSIZE];
  const auto conv =
      std::to_chars(sci, sci + sizeof(sci), val, std::chars_format::scientific);
  if (conv.ec != std::errc()) {
    return -1;
  }
  const char *q = sci;
  const bool neg = (*q == '-');
  if (neg) {
    q++;
  }
  char digits[24];
  int ndigits = 0;
  digits[ndigits++] = *q++;
  if (*q == '.') {
    q++;
    while (q < conv.ptr && *q != 'e') {
      digits[ndigits++] = *q++;
    }
  }
  // q now points at 'e'
  q++;
  bool exp_neg = false;
  if (*q == '+' || *q == '-') {
    exp_neg = (*q == '-');
    q++;
  }
  int exp10 = 0;
  while (q < conv.ptr) {
    exp10 = exp10 * 10 + (*q++ - '0');
  }
  if (exp_neg) {
    exp10 = -exp10;
  }
  const double aval = std::fabs(static_cast<double>(val));
  const bool use_fixed = (aval >= fixed_lo) && (aval < fixed_hi);
  return emit_digits(neg, digits, ndigits, exp10, use_fixed, out);
}

#endif

} // namespace

#ifdef PD_HAVE_FP_TO_CHARS

int pd_double_repr(double val, char *out) {
  return repr_to_chars(val, 1e-4, 1e16, out);
}

int pd_float32_repr(float val, char *out) {
  return repr_to_chars(val, 1e-4, 1e6, out);
}

int pd_float32_repr_available(void) { return 1; }

#else

// Fallback: PyOS_double_to_string repr mode gives the same shortest digits;
// re-layout them with numpy's value-based notation choice. Requires the GIL.
int pd_double_repr(double val, char *out) {
  int res = emit_special(val, out);
  if (res != 0) {
    return res;
  }
  char *pyrepr = PyOS_double_to_string(val, 'r', 0, 0, NULL);
  if (pyrepr == NULL) {
    return -1;
  }
  const char *q = pyrepr;
  const bool neg = (*q == '-');
  if (neg) {
    q++;
  }
  char digits[32];
  int ndigits = 0;
  int int_len = -1;
  int exp_part = 0;
  bool exp_neg = false;
  for (; *q != '\0'; q++) {
    if (*q == '.') {
      int_len = ndigits;
    } else if (*q == 'e') {
      q++;
      if (*q == '+' || *q == '-') {
        exp_neg = (*q == '-');
        q++;
      }
      for (; *q != '\0'; q++) {
        exp_part = exp_part * 10 + (*q - '0');
      }
      break;
    } else {
      digits[ndigits++] = *q;
    }
  }
  PyMem_Free(pyrepr);
  if (int_len == -1) {
    int_len = ndigits;
  }
  if (exp_neg) {
    exp_part = -exp_part;
  }
  int exp10 = exp_part + int_len - 1;
  // normalize: strip leading zeros ("0.0001" -> digits "1", exp10 -4) and
  // trailing zeros ("1.0" -> digits "1")
  int start = 0;
  while (start < ndigits - 1 && digits[start] == '0') {
    start++;
    exp10--;
  }
  while (ndigits - start > 1 && digits[ndigits - 1] == '0') {
    ndigits--;
  }
  const double aval = std::fabs(val);
  const bool use_fixed = (aval >= 1e-4) && (aval < 1e16);
  return emit_digits(neg, digits + start, ndigits - start, exp10, use_fixed,
                     out);
}

int pd_float32_repr(float val, char *out) {
  (void)val;
  (void)out;
  return -1;
}

int pd_float32_repr_available(void) { return 0; }

#endif
