/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_INLINE_HELPER_H_
#define PANDAS__LIBS_SRC_INLINE_HELPER_H_

#ifndef PANDAS_INLINE
  #if defined(__GNUC__)
    #define PANDAS_INLINE static __inline__
  #elif defined(_MSC_VER)
    #define PANDAS_INLINE static __inline
  #elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define PANDAS_INLINE static inline
  #else
    #define PANDAS_INLINE
  #endif
#endif

#endif  // PANDAS__LIBS_SRC_INLINE_HELPER_H_
