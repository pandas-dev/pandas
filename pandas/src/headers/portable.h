#ifndef _PANDAS_PORTABLE_H_
#define _PANDAS_PORTABLE_H_

#if defined(_MSC_VER)
#define strcasecmp( s1, s2 ) _stricmp( s1, s2 )
#endif

#endif
