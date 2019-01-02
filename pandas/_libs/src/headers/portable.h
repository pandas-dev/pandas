#ifndef _PANDAS_PORTABLE_H_
#define _PANDAS_PORTABLE_H_

#if defined(_MSC_VER)
#define strcasecmp( s1, s2 ) _stricmp( s1, s2 )
#endif

// GH-23516 - works around locale perf issues
// from MUSL libc, MIT Licensed - see LICENSES
#define isdigit_ascii(c) ((unsigned)c - '0' < 10)
#define isspace_ascii(c) (c == ' ' || (unsigned)c-'\t' < 5)
#define toupper_ascii(c) (((unsigned)c-'a' < 26) ? (c & 0x5f) : c)

#endif
