#ifndef C_HELPER_H
#define C_HELPER_H

#ifndef PANDAS_INLINE
  #if defined(__GNUC__)
    #define PANDAS_INLINE static __inline__
  #elif defined(_MSC_VER)
    #define PANDAS_INLINE static __inline
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define PANDAS_INLINE static inline
  #else
    #define PANDAS_INLINE
  #endif
#endif

#endif
