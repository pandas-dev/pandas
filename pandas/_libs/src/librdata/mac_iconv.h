/* Copyright (C) 1999-2003, 2005-2006 Free Software Foundation, Inc.
   This file is part of the GNU LIBICONV Library.

   The GNU LIBICONV Library is free software; you can redistribute it
   and/or modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   The GNU LIBICONV Library is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with the GNU LIBICONV Library; see the file COPYING.LIB.
   If not, write to the Free Software Foundation, Inc., 51 Franklin Street,
   Fifth Floor, Boston, MA 02110-1301, USA.  */

/* When installed, this file is called "iconv.h". */

#ifndef PANDAS__LIBS_SRC_LIBRDATA_MAC_ICONV_H_
#define PANDAS__LIBS_SRC_LIBRDATA_MAC_ICONV_H_

#ifndef _LIBICONV_H
#define _LIBICONV_H

#include <sys/cdefs.h>
#include <_types.h>
#include <sys/_types/_size_t.h>

#define _LIBICONV_VERSION 0x010B    /* version number: (major<<8) + minor */

#if BUILDING_LIBICONV
#define __LIBICONV_DLL_EXPORTED __attribute__((__visibility__("default")))
#else
#define __LIBICONV_DLL_EXPORTED
#endif
extern __LIBICONV_DLL_EXPORTED  int _libiconv_version; /* Likewise */

/* We would like to #include any system header file which could define
   iconv_t, 1. in order to eliminate the risk that the user gets compilation
   errors because some other system header file includes /usr/include/iconv.h
   which defines iconv_t or declares iconv after this file, 2. when compiling
   for LIBICONV_PLUG, we need the proper iconv_t type in order to produce
   binary compatible code.
   But gcc's #include_next is not portable. Thus, once libiconv's iconv.h
   has been installed in /usr/local/include, there is no way any more to
   include the original /usr/include/iconv.h. We simply have to get away
   without it.
   Ad 1. The risk that a system header file does
   #include "iconv.h"  or  #include_next "iconv.h"
   is small. They all do #include <iconv.h>.
   Ad 2. The iconv_t type is a pointer type in all cases I have seen. (It
   has to be a scalar type because (iconv_t)(-1) is a possible return value
   from iconv_open().) */

/* Define iconv_t ourselves. */
#ifndef _ICONV_T
#define _ICONV_T
typedef void* iconv_t;
#endif


#ifdef __cplusplus
extern "C" {
#endif


/* Allocates descriptor for code conversion from encoding `fromcode' to
   encoding `tocode'. */
extern __LIBICONV_DLL_EXPORTED iconv_t iconv_open(
    const char* __tocode,
    const char* __fromcode);

/* Converts, using conversion descriptor `cd', at most `*inbytesleft' bytes
   starting at `*inbuf', writing at most `*outbytesleft' bytes starting at
   `*outbuf'.
   Decrements `*inbytesleft' and increments `*inbuf' by the same amount.
   Decrements `*outbytesleft' and increments `*outbuf' by the same amount. */
extern __LIBICONV_DLL_EXPORTED size_t iconv(
    iconv_t __cd,
    char* * __restrict __inbuf,
    size_t * __restrict __inbytesleft,
    char* * __restrict __outbuf,
    size_t * __restrict __outbytesleft);

/* Frees resources allocated for conversion descriptor `cd'. */
extern __LIBICONV_DLL_EXPORTED int iconv_close(iconv_t _cd);

#if !defined(_POSIX_C_SOURCE) || defined(_DARWIN_C_SOURCE)

/* Nonstandard extensions. */

#include <sys/_types/_wchar_t.h>

/* Control of attributes. */
extern __LIBICONV_DLL_EXPORTED int iconvctl(
    iconv_t cd,
    int request,
    void* argument);

/* Hook performed after every successful conversion of a Unicode character. */
typedef void (*iconv_unicode_char_hook)(unsigned int uc, void* data);
/* Hook performed after every successful conversion of a wide character. */
typedef void (*iconv_wide_char_hook)(wchar_t wc, void* data);
/* Set of hooks. */
struct iconv_hooks {
  iconv_unicode_char_hook uc_hook;
  iconv_wide_char_hook wc_hook;
  void* data;
};

/* Fallback function.  Invoked when a small number of bytes could not be
   converted to a Unicode character.  This function should process all
   bytes from inbuf and may produce replacement Unicode characters by calling
   the write_replacement callback repeatedly.  */
typedef void (*iconv_unicode_mb_to_uc_fallback)(
    const char* inbuf, size_t inbufsize,
    void (*write_replacement)(
        const unsigned int *buf,
        size_t buflen,
        void* callback_arg),
    void* callback_arg,
    void* data);
/* Fallback function.  Invoked when a Unicode character could not be converted
   to the target encoding.  This function should process the character and
   may produce replacement bytes (in the target encoding) by calling the
   write_replacement callback repeatedly.  */
typedef void (*iconv_unicode_uc_to_mb_fallback)(
    unsigned int code,
    void (*write_replacement)(
        const char *buf,
        size_t buflen,
        void* callback_arg),
    void* callback_arg,
    void* data);
#if 1
/* Fallback function.  Invoked when a number of bytes could not be converted to
   a wide character.  This function should process all bytes from inbuf and may
   produce replacement wide characters by calling the write_replacement
   callback repeatedly.  */
typedef void (*iconv_wchar_mb_to_wc_fallback)(
    const char* inbuf, size_t inbufsize,
    void (*write_replacement)(
        const wchar_t *buf,
        size_t buflen,
        void* callback_arg),
    void* callback_arg,
    void* data);
/* Fallback function.  Invoked when a wide character could not be converted to
   the target encoding.  This function should process the character and may
   produce replacement bytes (in the target encoding) by calling the
   write_replacement callback repeatedly.  */
typedef void (*iconv_wchar_wc_to_mb_fallback)(
    wchar_t code,
    void (*write_replacement)(
        const char *buf,
        size_t buflen,
        void* callback_arg),
    void* callback_arg,
    void* data);
#else
/* If the wchar_t type does not exist, these two fallback functions are never
   invoked.  Their argument list therefore does not matter.  */
typedef void (*iconv_wchar_mb_to_wc_fallback) ();
typedef void (*iconv_wchar_wc_to_mb_fallback) ();
#endif
/* Set of fallbacks. */
struct iconv_fallbacks {
  iconv_unicode_mb_to_uc_fallback mb_to_uc_fallback;
  iconv_unicode_uc_to_mb_fallback uc_to_mb_fallback;
  iconv_wchar_mb_to_wc_fallback mb_to_wc_fallback;
  iconv_wchar_wc_to_mb_fallback wc_to_mb_fallback;
  void* data;
};

/* Requests for iconvctl. */
#define ICONV_TRIVIALP            0  /* int *argument */
#define ICONV_GET_TRANSLITERATE   1  /* int *argument */
#define ICONV_SET_TRANSLITERATE   2  /* const int *argument */
#define ICONV_GET_DISCARD_ILSEQ   3  /* int *argument */
#define ICONV_SET_DISCARD_ILSEQ   4  /* const int *argument */
/* const struct iconv_hooks *argument */
#define ICONV_SET_HOOKS           5
/* const struct iconv_fallbacks *argument */
#define ICONV_SET_FALLBACKS       6

/* Listing of locale independent encodings. */
extern __LIBICONV_DLL_EXPORTED void iconvlist(
    int (*do_one)(
        unsigned int namescount,
        const char * const * names,
        void* data),
    void* data);

/* Canonicalize an encoding name.
   The result is either a canonical encoding name, or name itself. */
extern __LIBICONV_DLL_EXPORTED const char * iconv_canonicalize(
    const char * name);

/* Support for relocatable packages.  */

/* Sets the original and the current installation prefix of the package.
   Relocation simply replaces a pathname starting with the original prefix
   by the corresponding pathname with the current prefix instead.  Both
   prefixes should be directory names without trailing slash (i.e. use ""
   instead of "/").  */
extern __LIBICONV_DLL_EXPORTED void libiconv_set_relocation_prefix(
    const char *orig_prefix,
    const char *curr_prefix);

#endif /* (!_POSIX_C_SOURCE || _DARWIN_C_SOURCE) */


#ifdef __cplusplus
}
#endif


#endif /* _LIBICONV_H */

#endif  // PANDAS__LIBS_SRC_LIBRDATA_MAC_ICONV_H_
