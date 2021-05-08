/*
 * No Copyright.
 *
 * iconv implementation using Win32 API to convert.
 * This file is placed in the public domain.
 */

#ifndef PANDAS_IO_RDATA_LIBRDATA_WIN_ICONV_H_
#define PANDAS_IO_RDATA_LIBRDATA_WIN_ICONV_H_

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

#endif  // PANDAS_IO_RDATA_LIBRDATA_WIN_ICONV_H_
