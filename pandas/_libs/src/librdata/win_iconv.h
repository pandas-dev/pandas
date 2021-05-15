/*

win-iconv - iconv implementation using Win32 API to convert.
Written in 2009-2016 by Yukihiro Nakadaira <https://github.com/ynkdir>
and contributors to win-iconv <https://github.com/win-iconv/win-iconv>

To the extent possible under law, the author(s) have dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

You should have received a copy of the CC0 Public Domain Dedication along with
this software. If not, see http://creativecommons.org/publicdomain/zero/1.0/.

 */

#ifndef PANDAS__LIBS_SRC_LIBRDATA_WIN_ICONV_H_
#define PANDAS__LIBS_SRC_LIBRDATA_WIN_ICONV_H_

//  #ifndef _LIBICONV_H
    #define _LIBICONV_H
    #include <stddef.h>
    #ifndef WINICONV_CONST
        # ifdef ICONV_CONST
            #  define WINICONV_CONST ICONV_CONST
        # else
            #  define WINICONV_CONST const
            # endif
    #endif
    #ifdef __cplusplus
        extern "C" {
    #endif

    typedef void* iconv_t;
    iconv_t iconv_open(const char *tocode, const char *fromcode);
    int iconv_close(iconv_t cd);
    size_t iconv(
        iconv_t cd,
        WINICONV_CONST char **inbuf,
        size_t *inbytesleft,
        char **outbuf,
        size_t *outbytesleft);

    #ifdef __cplusplus
        }
    #endif
//  #endif

#endif  // PANDAS__LIBS_SRC_LIBRDATA_WIN_ICONV_H_
