#ifndef BASE64_CONFIG_H
#define BASE64_CONFIG_H

#if !defined(__APPLE__) && ((defined(__x86_64__) && defined(__LP64__)) || defined(_M_X64))
  #define HAVE_SSSE3 1
  #define HAVE_SSE41 1
  #define HAVE_SSE42 1
  #define HAVE_AVX 1
  #define HAVE_AVX2 1
  #define HAVE_AVX512 0
#elif (defined(__APPLE__) && defined(__aarch64__))
  #define HAVE_NEON64 1
#elif (defined(__wasm__) && defined(__wasm_simd128__))
  #include "emscripten/version.h"
  #if __EMSCRIPTEN_major__ == 3
    #define HAVE_NEON32 1
  #elif __EMSCRIPTEN_major__ > 3
    #define HAVE_NEON64 1
  #endif
#endif

#endif // BASE64_CONFIG_H
