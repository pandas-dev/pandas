#pragma once

#include <string.h>

#if defined(_MSC_VER)
#define strcasecmp( s1, s2 ) _stricmp( s1, s2 )
#endif

// GH-23516 - works around locale perf issues
// from MUSL libc, MIT Licensed - see LICENSES
#define isdigit_ascii(c) (((unsigned)(c) - '0') < 10u)
#define getdigit_ascii(c, default) (isdigit_ascii(c) ? ((int)((c) - '0')) : default)
