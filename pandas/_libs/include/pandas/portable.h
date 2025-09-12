/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#include <string.h>

#if defined(_MSC_VER)
#define strcasecmp(s1, s2) _stricmp(s1, s2)
#endif

// GH-23516 - works around locale perf issues
// from MUSL libc, licence at LICENSES/MUSL_LICENSE
#define isdigit_ascii(c) (((unsigned)(c) - '0') < 10u)
#define getdigit_ascii(c, default)                                             \
  (isdigit_ascii(c) ? ((int)((c) - '0')) : default)
#define isspace_ascii(c) (((c) == ' ') || (((unsigned)(c) - '\t') < 5))
#define toupper_ascii(c) ((((unsigned)(c) - 'a') < 26) ? ((c) & 0x5f) : (c))
#define tolower_ascii(c) ((((unsigned)(c) - 'A') < 26) ? ((c) | 0x20) : (c))

#if defined(_WIN32)
#define PD_FALLTHROUGH                                                         \
  do {                                                                         \
  } while (0) /* fallthrough */
#elif __has_attribute(__fallthrough__)
#define PD_FALLTHROUGH __attribute__((__fallthrough__))
#else
#define PD_FALLTHROUGH                                                         \
  do {                                                                         \
  } while (0) /* fallthrough */
#endif
