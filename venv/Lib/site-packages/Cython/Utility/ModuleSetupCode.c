/////////////// InitLimitedAPI ///////////////

#if defined(Py_LIMITED_API)
  #if !defined(CYTHON_LIMITED_API)
  // Use Py_LIMITED_API as the main control for Cython's limited API mode.
  // However it's still possible to define CYTHON_LIMITED_API alone to
  // force Cython to use Limited-API code without enforcing it in Python.
  #define CYTHON_LIMITED_API 1
  #endif
#elif defined(CYTHON_LIMITED_API)
  #ifdef _MSC_VER
  #pragma message ("Limited API usage is enabled with 'CYTHON_LIMITED_API' but 'Py_LIMITED_API' does not define a Python target version. Consider setting 'Py_LIMITED_API' instead.")
  #else
  #warning Limited API usage is enabled with 'CYTHON_LIMITED_API' but 'Py_LIMITED_API' does not define a Python target version. Consider setting 'Py_LIMITED_API' instead.
  #endif
#endif

/////////////// CModulePreamble ///////////////

#include <stddef.h> /* For offsetof */
#ifndef offsetof
  #define offsetof(type, member) ( (size_t) & ((type*)0) -> member )
#endif

#if !defined(_WIN32) && !defined(WIN32) && !defined(MS_WINDOWS)
  #ifndef __stdcall
    #define __stdcall
  #endif
  #ifndef __cdecl
    #define __cdecl
  #endif
  #ifndef __fastcall
    #define __fastcall
  #endif
#endif

#ifndef DL_IMPORT
  #define DL_IMPORT(t) t
#endif
#ifndef DL_EXPORT
  #define DL_EXPORT(t) t
#endif

// For use in DL_IMPORT/DL_EXPORT macros.
#define __PYX_COMMA ,

#ifndef PY_LONG_LONG
  #define PY_LONG_LONG LONG_LONG
#endif

#ifndef Py_HUGE_VAL
  #define Py_HUGE_VAL HUGE_VAL
#endif

// For the limited API it often makes sense to use Py_LIMITED_API rather than PY_VERSION_HEX
// when doing version checks.
#define __PYX_LIMITED_VERSION_HEX PY_VERSION_HEX

#if defined(GRAALVM_PYTHON)
  /* For very preliminary testing purposes. Most variables are set the same as PyPy.
     The existence of this section does not imply that anything works or is even tested */
  // GRAALVM_PYTHON test comes before PyPy test because GraalPython unhelpfully defines PYPY_VERSION
  #define CYTHON_COMPILING_IN_PYPY 0
  #define CYTHON_COMPILING_IN_CPYTHON 0
  #define CYTHON_COMPILING_IN_LIMITED_API 0
  #define CYTHON_COMPILING_IN_GRAAL 1

  #define CYTHON_COMPILING_IN_CPYTHON_FREETHREADING 0

  #undef CYTHON_USE_TYPE_SLOTS
  #define CYTHON_USE_TYPE_SLOTS 0
  #undef CYTHON_USE_TYPE_SPECS
  #define CYTHON_USE_TYPE_SPECS 0
  #undef CYTHON_USE_PYTYPE_LOOKUP
  #define CYTHON_USE_PYTYPE_LOOKUP 0
  #undef CYTHON_USE_PYLIST_INTERNALS
  #define CYTHON_USE_PYLIST_INTERNALS 0
  #undef CYTHON_USE_UNICODE_INTERNALS
  #define CYTHON_USE_UNICODE_INTERNALS 0
  #undef CYTHON_USE_UNICODE_WRITER
  #define CYTHON_USE_UNICODE_WRITER 0
  #undef CYTHON_USE_PYLONG_INTERNALS
  #define CYTHON_USE_PYLONG_INTERNALS 0
  #undef CYTHON_AVOID_BORROWED_REFS
  #define CYTHON_AVOID_BORROWED_REFS 1
  #undef CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
  #define CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS 0
  #undef CYTHON_ASSUME_SAFE_MACROS
  #define CYTHON_ASSUME_SAFE_MACROS 0
  #undef CYTHON_ASSUME_SAFE_SIZE
  #define CYTHON_ASSUME_SAFE_SIZE 0
  #undef CYTHON_UNPACK_METHODS
  #define CYTHON_UNPACK_METHODS 0
  #undef CYTHON_FAST_THREAD_STATE
  #define CYTHON_FAST_THREAD_STATE 0
  #undef CYTHON_FAST_GIL
  #define CYTHON_FAST_GIL 0
  #undef CYTHON_METH_FASTCALL
  #define CYTHON_METH_FASTCALL 0
  #undef CYTHON_FAST_PYCALL
  #define CYTHON_FAST_PYCALL 0
  #ifndef CYTHON_PEP487_INIT_SUBCLASS
    #define CYTHON_PEP487_INIT_SUBCLASS 1
  #endif
  #undef CYTHON_PEP489_MULTI_PHASE_INIT
  #define CYTHON_PEP489_MULTI_PHASE_INIT 1
  #undef CYTHON_USE_MODULE_STATE
  #define CYTHON_USE_MODULE_STATE 0
  #undef CYTHON_USE_SYS_MONITORING
  #define CYTHON_USE_SYS_MONITORING 0
  #undef CYTHON_USE_TP_FINALIZE
  #define CYTHON_USE_TP_FINALIZE 0
  #undef CYTHON_USE_AM_SEND
  #define CYTHON_USE_AM_SEND 0
  #undef CYTHON_USE_DICT_VERSIONS
  #define CYTHON_USE_DICT_VERSIONS 0
  #undef CYTHON_USE_EXC_INFO_STACK
  #define CYTHON_USE_EXC_INFO_STACK 1
  #ifndef CYTHON_UPDATE_DESCRIPTOR_DOC
    #define CYTHON_UPDATE_DESCRIPTOR_DOC 0
  #endif
  #undef CYTHON_USE_FREELISTS
  #define CYTHON_USE_FREELISTS 0
  #undef CYTHON_IMMORTAL_CONSTANTS
  #define CYTHON_IMMORTAL_CONSTANTS 0

#elif defined(PYPY_VERSION)
  #define CYTHON_COMPILING_IN_PYPY 1
  #define CYTHON_COMPILING_IN_CPYTHON 0
  #define CYTHON_COMPILING_IN_LIMITED_API 0
  #define CYTHON_COMPILING_IN_GRAAL 0

  #define CYTHON_COMPILING_IN_CPYTHON_FREETHREADING 0

  #undef CYTHON_USE_TYPE_SLOTS
  #define CYTHON_USE_TYPE_SLOTS 1
  #ifndef CYTHON_USE_TYPE_SPECS
    #define CYTHON_USE_TYPE_SPECS 0
  #endif
  #undef CYTHON_USE_PYTYPE_LOOKUP
  #define CYTHON_USE_PYTYPE_LOOKUP 0
  #undef CYTHON_USE_PYLIST_INTERNALS
  #define CYTHON_USE_PYLIST_INTERNALS 0
  #undef CYTHON_USE_UNICODE_INTERNALS
  #define CYTHON_USE_UNICODE_INTERNALS 0
  #undef CYTHON_USE_UNICODE_WRITER
  #define CYTHON_USE_UNICODE_WRITER 0
  #undef CYTHON_USE_PYLONG_INTERNALS
  #define CYTHON_USE_PYLONG_INTERNALS 0
  #undef CYTHON_AVOID_BORROWED_REFS
  #define CYTHON_AVOID_BORROWED_REFS 1
  #undef CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
  #define CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS 1
  #undef CYTHON_ASSUME_SAFE_MACROS
  #define CYTHON_ASSUME_SAFE_MACROS 0
  #ifndef CYTHON_ASSUME_SAFE_SIZE
    #define CYTHON_ASSUME_SAFE_SIZE 1
  #endif
  #undef CYTHON_UNPACK_METHODS
  #define CYTHON_UNPACK_METHODS 0
  #undef CYTHON_FAST_THREAD_STATE
  #define CYTHON_FAST_THREAD_STATE 0
  #undef CYTHON_FAST_GIL
  #define CYTHON_FAST_GIL 0
  #undef CYTHON_METH_FASTCALL
  #define CYTHON_METH_FASTCALL 0
  #undef CYTHON_FAST_PYCALL
  #define CYTHON_FAST_PYCALL 0
  #ifndef CYTHON_PEP487_INIT_SUBCLASS
    #define CYTHON_PEP487_INIT_SUBCLASS 1
  #endif
  #if PY_VERSION_HEX < 0x03090000
    #undef CYTHON_PEP489_MULTI_PHASE_INIT
    #define CYTHON_PEP489_MULTI_PHASE_INIT 0
  #elif !defined(CYTHON_PEP489_MULTI_PHASE_INIT)
    #define CYTHON_PEP489_MULTI_PHASE_INIT 1
  #endif
  #undef CYTHON_USE_MODULE_STATE
  #define CYTHON_USE_MODULE_STATE 0
  #undef CYTHON_USE_SYS_MONITORING
  #define CYTHON_USE_SYS_MONITORING 0
  #ifndef CYTHON_USE_TP_FINALIZE
    #define CYTHON_USE_TP_FINALIZE (PYPY_VERSION_NUM >= 0x07030C00)
  #endif
  #undef CYTHON_USE_AM_SEND
  #define CYTHON_USE_AM_SEND 0
  #undef CYTHON_USE_DICT_VERSIONS
  #define CYTHON_USE_DICT_VERSIONS 0
  #undef CYTHON_USE_EXC_INFO_STACK
  #define CYTHON_USE_EXC_INFO_STACK 0
  #ifndef CYTHON_UPDATE_DESCRIPTOR_DOC
    #define CYTHON_UPDATE_DESCRIPTOR_DOC (PYPY_VERSION_NUM >= 0x07031100)
  #endif
  #undef CYTHON_USE_FREELISTS
  #define CYTHON_USE_FREELISTS 0
  #undef CYTHON_IMMORTAL_CONSTANTS
  #define CYTHON_IMMORTAL_CONSTANTS 0

#elif defined(CYTHON_LIMITED_API)
  // EXPERIMENTAL !!
  #ifdef Py_LIMITED_API
    #undef __PYX_LIMITED_VERSION_HEX
    #define __PYX_LIMITED_VERSION_HEX Py_LIMITED_API
  #endif
  #define CYTHON_COMPILING_IN_PYPY 0
  #define CYTHON_COMPILING_IN_CPYTHON 0
  #define CYTHON_COMPILING_IN_LIMITED_API 1
  #define CYTHON_COMPILING_IN_GRAAL 0

  #define CYTHON_COMPILING_IN_CPYTHON_FREETHREADING 0

  #undef CYTHON_USE_TYPE_SLOTS
  #define CYTHON_USE_TYPE_SLOTS 0
  #undef CYTHON_USE_TYPE_SPECS
  #define CYTHON_USE_TYPE_SPECS 1
  #undef CYTHON_USE_PYTYPE_LOOKUP
  #define CYTHON_USE_PYTYPE_LOOKUP 0
  #undef CYTHON_USE_PYLIST_INTERNALS
  #define CYTHON_USE_PYLIST_INTERNALS 0
  #undef CYTHON_USE_UNICODE_INTERNALS
  #define CYTHON_USE_UNICODE_INTERNALS 0
  #ifndef CYTHON_USE_UNICODE_WRITER
    #define CYTHON_USE_UNICODE_WRITER 0
  #endif
  #undef CYTHON_USE_PYLONG_INTERNALS
  #define CYTHON_USE_PYLONG_INTERNALS 0
  #ifndef CYTHON_AVOID_BORROWED_REFS
    #define CYTHON_AVOID_BORROWED_REFS 0
  #endif
  #ifndef CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
    #define CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS 0
  #endif
  #undef CYTHON_ASSUME_SAFE_MACROS
  #define CYTHON_ASSUME_SAFE_MACROS 0
  #undef CYTHON_ASSUME_SAFE_SIZE
  #define CYTHON_ASSUME_SAFE_SIZE 0
  #undef CYTHON_UNPACK_METHODS
  #define CYTHON_UNPACK_METHODS 0
  #undef CYTHON_FAST_THREAD_STATE
  #define CYTHON_FAST_THREAD_STATE 0
  #undef CYTHON_FAST_GIL
  #define CYTHON_FAST_GIL 0
  #undef CYTHON_METH_FASTCALL
  #define CYTHON_METH_FASTCALL (__PYX_LIMITED_VERSION_HEX >= 0x030C0000)
  #undef CYTHON_FAST_PYCALL
  #define CYTHON_FAST_PYCALL 0
  #ifndef CYTHON_PEP487_INIT_SUBCLASS
    #define CYTHON_PEP487_INIT_SUBCLASS 1
  #endif
  #ifndef CYTHON_PEP489_MULTI_PHASE_INIT
    #define CYTHON_PEP489_MULTI_PHASE_INIT 1
  #endif
  #ifndef CYTHON_USE_MODULE_STATE
    #define CYTHON_USE_MODULE_STATE 0
  #endif
  #undef CYTHON_USE_SYS_MONITORING
  #define CYTHON_USE_SYS_MONITORING 0
  #ifndef CYTHON_USE_TP_FINALIZE
    // PyObject_CallFinalizerFromDealloc is missing and not easily replaced
    #define CYTHON_USE_TP_FINALIZE 0
  #endif
  #ifndef CYTHON_USE_AM_SEND
    #define CYTHON_USE_AM_SEND (__PYX_LIMITED_VERSION_HEX >= 0x030A0000)
  #endif
  #undef CYTHON_USE_DICT_VERSIONS
  #define CYTHON_USE_DICT_VERSIONS 0
  #undef CYTHON_USE_EXC_INFO_STACK
  #define CYTHON_USE_EXC_INFO_STACK 0
  #ifndef CYTHON_UPDATE_DESCRIPTOR_DOC
    #define CYTHON_UPDATE_DESCRIPTOR_DOC 0
  #endif
  #ifndef CYTHON_USE_FREELISTS
  #define CYTHON_USE_FREELISTS 1
  #endif
  #undef CYTHON_IMMORTAL_CONSTANTS
  #define CYTHON_IMMORTAL_CONSTANTS 0

#else
  #define CYTHON_COMPILING_IN_PYPY 0
  #define CYTHON_COMPILING_IN_CPYTHON 1
  #define CYTHON_COMPILING_IN_LIMITED_API 0
  #define CYTHON_COMPILING_IN_GRAAL 0

  #ifdef Py_GIL_DISABLED
    #define CYTHON_COMPILING_IN_CPYTHON_FREETHREADING 1
  #else
    #define CYTHON_COMPILING_IN_CPYTHON_FREETHREADING 0
  #endif

  #if PY_VERSION_HEX < 0x030A0000
    // Before Py3.10, PyObject_GetSlot() rejects static (non-heap) types.
    #undef CYTHON_USE_TYPE_SLOTS
    #define CYTHON_USE_TYPE_SLOTS 1
  #elif !defined(CYTHON_USE_TYPE_SLOTS)
    #define CYTHON_USE_TYPE_SLOTS 1
  #endif
  #ifndef CYTHON_USE_TYPE_SPECS
    #define CYTHON_USE_TYPE_SPECS 0
  #endif
  #ifndef CYTHON_USE_PYTYPE_LOOKUP
    #define CYTHON_USE_PYTYPE_LOOKUP 1
  #endif
  #ifndef CYTHON_USE_PYLONG_INTERNALS
    #define CYTHON_USE_PYLONG_INTERNALS 1
  #endif
  #if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    #undef CYTHON_USE_PYLIST_INTERNALS
    // Use thread-safe CPython C API calls to manipulate list contents
    #define CYTHON_USE_PYLIST_INTERNALS 0
  #elif !defined(CYTHON_USE_PYLIST_INTERNALS)
    #define CYTHON_USE_PYLIST_INTERNALS 1
  #endif
  #ifndef CYTHON_USE_UNICODE_INTERNALS
    #define CYTHON_USE_UNICODE_INTERNALS 1
  #endif
  #if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING || PY_VERSION_HEX >= 0x030B00A2
    // Python 3.11a2 hid _PyLong_FormatAdvancedWriter and _PyFloat_FormatAdvancedWriter
    // therefore disable unicode writer until a better alternative appears
    #undef CYTHON_USE_UNICODE_WRITER
    #define CYTHON_USE_UNICODE_WRITER 0
  #elif !defined(CYTHON_USE_UNICODE_WRITER)
    #define CYTHON_USE_UNICODE_WRITER 1
  #endif
  // CYTHON_AVOID_BORROWED_REFS - Avoid borrowed references and always request owned references directly instead.
  #ifndef CYTHON_AVOID_BORROWED_REFS
    #define CYTHON_AVOID_BORROWED_REFS 0
  #endif
  // CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS - Avoid borrowed references that are not thread-safe in the free-threaded build of CPython.
  #if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    #undef CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
    #define CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS 1
  #elif !defined(CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS)
    #define CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS 0
  #endif
  // CYTHON_ASSUME_SAFE_MACROS - Assume that macro calls do not fail and do not raise exceptions.
  #ifndef CYTHON_ASSUME_SAFE_MACROS
    #define CYTHON_ASSUME_SAFE_MACROS 1
  #endif
  // CYTHON_ASSUME_SAFE_SIZE - Assume that Py*_GET_SIZE() calls do not fail and do not raise exceptions.
  #ifndef CYTHON_ASSUME_SAFE_SIZE
    #define CYTHON_ASSUME_SAFE_SIZE 1
  #endif
  #ifndef CYTHON_UNPACK_METHODS
    #define CYTHON_UNPACK_METHODS 1
  #endif
  #ifndef CYTHON_FAST_THREAD_STATE
    #define CYTHON_FAST_THREAD_STATE 1
  #endif
  #if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    #undef CYTHON_FAST_GIL
    #define CYTHON_FAST_GIL 0
  #elif !defined(CYTHON_FAST_GIL)
    // FIXME: FastGIL can probably be supported also in CPython 3.12 but needs to be adapted.
    // The gain is unclear, however, since the GIL handling itself became faster in recent CPython versions.
    #define CYTHON_FAST_GIL (PY_VERSION_HEX < 0x030C00A6)
  #endif
  #ifndef CYTHON_METH_FASTCALL
    // CPython 3.6 introduced METH_FASTCALL but with slightly different
    // semantics. It became stable starting from CPython 3.7.
    #define CYTHON_METH_FASTCALL 1
  #endif
  #ifndef CYTHON_FAST_PYCALL
    #define CYTHON_FAST_PYCALL 1
  #endif
  #ifndef CYTHON_PEP487_INIT_SUBCLASS
    #define CYTHON_PEP487_INIT_SUBCLASS 1
  #endif
  #ifndef CYTHON_PEP489_MULTI_PHASE_INIT
    #define CYTHON_PEP489_MULTI_PHASE_INIT 1
  #endif
  // CYTHON_USE_MODULE_STATE - Use a module state/globals struct tied to the module object.
  #ifndef CYTHON_USE_MODULE_STATE
    // EXPERIMENTAL !!
    #define CYTHON_USE_MODULE_STATE 0
  #endif
  #ifndef CYTHON_USE_SYS_MONITORING
    #define CYTHON_USE_SYS_MONITORING (PY_VERSION_HEX >= 0x030d00B1)
  #endif
  #ifndef CYTHON_USE_TP_FINALIZE
    #define CYTHON_USE_TP_FINALIZE 1
  #endif
  #ifndef CYTHON_USE_AM_SEND
    #define CYTHON_USE_AM_SEND 1
  #endif
  #if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    #undef CYTHON_USE_DICT_VERSIONS
    #define CYTHON_USE_DICT_VERSIONS 0
  #elif !defined(CYTHON_USE_DICT_VERSIONS)
    // Python 3.12a5 deprecated "ma_version_tag"
    // and we use static variables with dict versions so it's incompatible with module state
    #define CYTHON_USE_DICT_VERSIONS  (PY_VERSION_HEX < 0x030C00A5 && !CYTHON_USE_MODULE_STATE)
  #endif
  #ifndef CYTHON_USE_EXC_INFO_STACK
    #define CYTHON_USE_EXC_INFO_STACK 1
  #endif
  #ifndef CYTHON_UPDATE_DESCRIPTOR_DOC
    #define CYTHON_UPDATE_DESCRIPTOR_DOC 1
  #endif
  #ifndef CYTHON_USE_FREELISTS
    #define CYTHON_USE_FREELISTS (!CYTHON_COMPILING_IN_CPYTHON_FREETHREADING)
  #endif
  #if defined(CYTHON_IMMORTAL_CONSTANTS) && PY_VERSION_HEX < 0x030C0000
    #undef CYTHON_IMMORTAL_CONSTANTS
    #define CYTHON_IMMORTAL_CONSTANTS 0  // definitely won't work
  #elif !defined(CYTHON_IMMORTAL_CONSTANTS)
    // The assumption is that immortal constants mainly help reduce contention when accessed from many threads
    // (hence freethreading only), and that with module state the module may be unloaded and reloaded (so they
    // could be a memory leak). However the user can always override these assumptions.
    #define CYTHON_IMMORTAL_CONSTANTS (PY_VERSION_HEX >= 0x030C0000 && !CYTHON_USE_MODULE_STATE && CYTHON_COMPILING_IN_CPYTHON_FREETHREADING)
  #endif
#endif

#ifndef CYTHON_COMPRESS_STRINGS
  #define CYTHON_COMPRESS_STRINGS 1
#endif

#ifndef CYTHON_FAST_PYCCALL
#define CYTHON_FAST_PYCCALL  CYTHON_FAST_PYCALL
#endif

#ifndef CYTHON_VECTORCALL
#if CYTHON_COMPILING_IN_LIMITED_API
// Possibly needs a bit of clearing up, however:
//  the limited API doesn't define CYTHON_FAST_PYCCALL (because that involves
//  a lot of access to internals) but does define CYTHON_VECTORCALL because
//  that's available cleanly from Python 3.12. Note that only VectorcallDict isn't
//  available though.
#define CYTHON_VECTORCALL  (__PYX_LIMITED_VERSION_HEX >= 0x030C0000)
#else
#define CYTHON_VECTORCALL  (CYTHON_FAST_PYCCALL)
#endif
#endif

#if CYTHON_USE_PYLONG_INTERNALS
  /* These short defines from the PyLong header can easily conflict with other code */
  #undef SHIFT
  #undef BASE
  #undef MASK
  /* Compile-time sanity check that these are indeed equal.  Github issue #2670. */
  #ifdef SIZEOF_VOID_P
    enum { __pyx_check_sizeof_voidp = 1 / (int)(SIZEOF_VOID_P == sizeof(void*)) };
  #endif
#endif

#ifndef __has_attribute
  #define __has_attribute(x) 0
#endif

#ifndef __has_cpp_attribute
  #define __has_cpp_attribute(x) 0
#endif

// restrict
#ifndef CYTHON_RESTRICT
  #if defined(__GNUC__)
    #define CYTHON_RESTRICT __restrict__
  #elif defined(_MSC_VER) && _MSC_VER >= 1400
    #define CYTHON_RESTRICT __restrict
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define CYTHON_RESTRICT restrict
  #else
    #define CYTHON_RESTRICT
  #endif
#endif

// unused attribute
#ifndef CYTHON_UNUSED
  #if defined(__cplusplus)
    /* for clang __has_cpp_attribute(maybe_unused) is true even before C++17
     * but leads to warnings with -pedantic, since it is a C++17 feature */
    #if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus >= 201703L)
      #if __has_cpp_attribute(maybe_unused)
        #define CYTHON_UNUSED [[maybe_unused]]
      #endif
    #endif
  #endif
#endif
#ifndef CYTHON_UNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define CYTHON_UNUSED __attribute__ ((__unused__))
#   else
#     define CYTHON_UNUSED
#   endif
# elif defined(__ICC) || (defined(__INTEL_COMPILER) && !defined(_MSC_VER))
#   define CYTHON_UNUSED __attribute__ ((__unused__))
# else
#   define CYTHON_UNUSED
# endif
#endif

#ifndef CYTHON_UNUSED_VAR
#  if defined(__cplusplus)
     template<class T> void CYTHON_UNUSED_VAR( const T& ) { }
#  else
#    define CYTHON_UNUSED_VAR(x) (void)(x)
#  endif
#endif

#ifndef CYTHON_MAYBE_UNUSED_VAR
  #define CYTHON_MAYBE_UNUSED_VAR(x) CYTHON_UNUSED_VAR(x)
#endif

#ifndef CYTHON_NCP_UNUSED
# if CYTHON_COMPILING_IN_CPYTHON && !CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
#  define CYTHON_NCP_UNUSED
# else
#  define CYTHON_NCP_UNUSED CYTHON_UNUSED
# endif
#endif

#ifndef CYTHON_USE_CPP_STD_MOVE
  // msvc doesn't set __cplusplus to a useful value
  #if defined(__cplusplus) && ( \
    __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1600))
    #define CYTHON_USE_CPP_STD_MOVE 1
  #else
    #define CYTHON_USE_CPP_STD_MOVE 0
  #endif
#endif

#define __Pyx_void_to_None(void_result) ((void)(void_result), Py_INCREF(Py_None), Py_None)

#include <stdint.h>
typedef uintptr_t  __pyx_uintptr_t;

#ifndef CYTHON_FALLTHROUGH
  #if defined(__cplusplus)
    /* for clang __has_cpp_attribute(fallthrough) is true even before C++17
     * but leads to warnings with -pedantic, since it is a C++17 feature */
    #if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus >= 201703L)
      #if __has_cpp_attribute(fallthrough)
        #define CYTHON_FALLTHROUGH [[fallthrough]]
      #endif
    #endif

    #ifndef CYTHON_FALLTHROUGH
      #if __has_cpp_attribute(clang::fallthrough)
        #define CYTHON_FALLTHROUGH [[clang::fallthrough]]
      #elif __has_cpp_attribute(gnu::fallthrough)
        #define CYTHON_FALLTHROUGH [[gnu::fallthrough]]
      #endif
    #endif
  #endif

  #ifndef CYTHON_FALLTHROUGH
    #if __has_attribute(fallthrough)
      #define CYTHON_FALLTHROUGH __attribute__((fallthrough))
    #else
      #define CYTHON_FALLTHROUGH
    #endif
  #endif

  #if defined(__clang__) && defined(__apple_build_version__)
    #if __apple_build_version__ < 7000000 /* Xcode < 7.0 */
      #undef  CYTHON_FALLTHROUGH
      #define CYTHON_FALLTHROUGH
    #endif
  #endif
#endif

#ifndef Py_UNREACHABLE
  #define Py_UNREACHABLE()  assert(0); abort()
#endif

#ifdef __cplusplus
  template <typename T>
  struct __PYX_IS_UNSIGNED_IMPL {static const bool value = T(0) < T(-1);};
  #define __PYX_IS_UNSIGNED(type) (__PYX_IS_UNSIGNED_IMPL<type>::value)
#else
  #define __PYX_IS_UNSIGNED(type) (((type)-1) > 0)
#endif

#if CYTHON_COMPILING_IN_PYPY == 1
  #define __PYX_NEED_TP_PRINT_SLOT  (PY_VERSION_HEX < 0x030A0000)
#else
  #define __PYX_NEED_TP_PRINT_SLOT  (PY_VERSION_HEX < 0x03090000)
#endif
// reinterpret

// TODO: refactor existing code to use those macros
#define __PYX_REINTERPRET_FUNCION(func_pointer, other_pointer) ((func_pointer)(void(*)(void))(other_pointer))
// #define __PYX_REINTERPRET_POINTER(pointer_type, pointer) ((pointer_type)(void *)(pointer))
// #define __PYX_RUNTIME_REINTERPRET(type, var) (*(type *)(&var))


/////////////// CInitCode ///////////////

// inline attribute
#ifndef CYTHON_INLINE
  #if defined(__clang__)
    #define CYTHON_INLINE __inline__ __attribute__ ((__unused__))
  #elif defined(__GNUC__)
    #define CYTHON_INLINE __inline__
  #elif defined(_MSC_VER)
    #define CYTHON_INLINE __inline
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define CYTHON_INLINE inline
  #else
    #define CYTHON_INLINE
  #endif
#endif


/////////////// CppInitCode ///////////////

#ifndef __cplusplus
  #error "Cython files generated with the C++ option must be compiled with a C++ compiler."
#endif

// inline attribute
#ifndef CYTHON_INLINE
  #if defined(__clang__)
    #define CYTHON_INLINE __inline__ __attribute__ ((__unused__))
  #else
    #define CYTHON_INLINE inline
  #endif
#endif

// Work around clang bug https://stackoverflow.com/questions/21847816/c-invoke-nested-template-class-destructor
// (even without the clang bug, the need not to know the typename is generally a benefit)
template<typename T>
void __Pyx_call_destructor(T& x) {
    x.~T();
}

// Used for temporary variables of "reference" type.
template<typename T>
class __Pyx_FakeReference {
  public:
    __Pyx_FakeReference() : ptr(NULL) { }
    // __Pyx_FakeReference(T& ref) : ptr(&ref) { }
    // Const version needed as Cython doesn't know about const overloads (e.g. for stl containers).
    __Pyx_FakeReference(const T& ref) : ptr(const_cast<T*>(&ref)) { }
    T *operator->() { return ptr; }
    T *operator&() { return ptr; }
    operator T&() { return *ptr; }
    // TODO(robertwb): Delegate all operators (or auto-generate unwrapping code where needed).
    template<typename U> bool operator ==(const U& other) const { return *ptr == other; }
    template<typename U> bool operator !=(const U& other) const { return *ptr != other; }
    template<typename U> bool operator==(const __Pyx_FakeReference<U>& other) const { return *ptr == *other.ptr; }
    template<typename U> bool operator!=(const __Pyx_FakeReference<U>& other) const { return *ptr != *other.ptr; }
  private:
    T *ptr;
};


/////////////// PythonCompatibility ///////////////
//@substitute: naming

#define __PYX_BUILD_PY_SSIZE_T "n"
#define CYTHON_FORMAT_SSIZE_T "z"

// TODO: remove this block
#define __Pyx_BUILTIN_MODULE_NAME "builtins"
#define __Pyx_DefaultClassType PyType_Type

#if CYTHON_COMPILING_IN_LIMITED_API
    // Cython uses these constants but they are not available in the limited API.
    // Therefore define them as static variables and look them up at module init.
    #ifndef CO_OPTIMIZED
    static int CO_OPTIMIZED;
    #endif
    #ifndef CO_NEWLOCALS
    static int CO_NEWLOCALS;
    #endif
    #ifndef CO_VARARGS
    static int CO_VARARGS;
    #endif
    #ifndef CO_VARKEYWORDS
    static int CO_VARKEYWORDS;
    #endif
    #ifndef CO_ASYNC_GENERATOR
    static int CO_ASYNC_GENERATOR;
    #endif
    #ifndef CO_GENERATOR
    static int CO_GENERATOR;
    #endif
    #ifndef CO_COROUTINE
    static int CO_COROUTINE;
    #endif
#else
    #ifndef CO_COROUTINE
      #define CO_COROUTINE 0x80
    #endif
    #ifndef CO_ASYNC_GENERATOR
      #define CO_ASYNC_GENERATOR 0x200
    #endif
#endif
static int __Pyx_init_co_variables(void); /* proto */

#if PY_VERSION_HEX >= 0x030900A4 || defined(Py_IS_TYPE)
  #define __Pyx_IS_TYPE(ob, type) Py_IS_TYPE(ob, type)
#else
  #define __Pyx_IS_TYPE(ob, type) (((const PyObject*)ob)->ob_type == (type))
#endif

#if PY_VERSION_HEX >= 0x030A00B1 || defined(Py_Is)
  #define __Pyx_Py_Is(x, y)  Py_Is(x, y)
#else
  #define __Pyx_Py_Is(x, y) ((x) == (y))
#endif
#if PY_VERSION_HEX >= 0x030A00B1 || defined(Py_IsNone)
  #define __Pyx_Py_IsNone(ob) Py_IsNone(ob)
#else
  #define __Pyx_Py_IsNone(ob) __Pyx_Py_Is((ob), Py_None)
#endif
#if PY_VERSION_HEX >= 0x030A00B1 || defined(Py_IsTrue)
  #define __Pyx_Py_IsTrue(ob) Py_IsTrue(ob)
#else
  #define __Pyx_Py_IsTrue(ob) __Pyx_Py_Is((ob), Py_True)
#endif
#if PY_VERSION_HEX >= 0x030A00B1 || defined(Py_IsFalse)
  #define __Pyx_Py_IsFalse(ob) Py_IsFalse(ob)
#else
  #define __Pyx_Py_IsFalse(ob) __Pyx_Py_Is((ob), Py_False)
#endif
#define __Pyx_NoneAsNull(obj)  (__Pyx_Py_IsNone(obj) ? NULL : (obj))

#if PY_VERSION_HEX >= 0x030900F0 && !CYTHON_COMPILING_IN_PYPY
  #define __Pyx_PyObject_GC_IsFinalized(o) PyObject_GC_IsFinalized(o)
#else
  #define __Pyx_PyObject_GC_IsFinalized(o) _PyGC_FINALIZED(o)
#endif

#ifndef Py_TPFLAGS_CHECKTYPES
  #define Py_TPFLAGS_CHECKTYPES 0
#endif
#ifndef Py_TPFLAGS_HAVE_INDEX
  #define Py_TPFLAGS_HAVE_INDEX 0
#endif
#ifndef Py_TPFLAGS_HAVE_NEWBUFFER
  #define Py_TPFLAGS_HAVE_NEWBUFFER 0
#endif
#ifndef Py_TPFLAGS_HAVE_FINALIZE
  #define Py_TPFLAGS_HAVE_FINALIZE 0
#endif
#ifndef Py_TPFLAGS_SEQUENCE
  #define Py_TPFLAGS_SEQUENCE 0
#endif
#ifndef Py_TPFLAGS_MAPPING
  #define Py_TPFLAGS_MAPPING 0
#endif
#ifndef Py_TPFLAGS_IMMUTABLETYPE
  #define Py_TPFLAGS_IMMUTABLETYPE (1UL << 8)
#endif
#ifndef Py_TPFLAGS_DISALLOW_INSTANTIATION
  #define Py_TPFLAGS_DISALLOW_INSTANTIATION (1UL << 7)
#endif

#ifndef METH_STACKLESS
  // already defined for Stackless Python (all versions) and C-Python >= 3.7
  // value if defined: Stackless Python < 3.6: 0x80 else 0x100
  #define METH_STACKLESS 0
#endif
#ifndef METH_FASTCALL
  // new in CPython 3.6, but changed in 3.7 - see
  // positional-only parameters:
  //   https://bugs.python.org/issue29464
  // const args:
  //   https://bugs.python.org/issue32240
  #ifndef METH_FASTCALL
     #define METH_FASTCALL 0x80
  #endif
  typedef PyObject *(*__Pyx_PyCFunctionFast) (PyObject *self, PyObject *const *args, Py_ssize_t nargs);
  // new in CPython 3.7, used to be old signature of _PyCFunctionFast() in 3.6
  typedef PyObject *(*__Pyx_PyCFunctionFastWithKeywords) (PyObject *self, PyObject *const *args,
                                                          Py_ssize_t nargs, PyObject *kwnames);
#else
  #if PY_VERSION_HEX >= 0x030d00A4
  #  define __Pyx_PyCFunctionFast PyCFunctionFast
  #  define __Pyx_PyCFunctionFastWithKeywords PyCFunctionFastWithKeywords
  #else
  #  define __Pyx_PyCFunctionFast _PyCFunctionFast
  #  define __Pyx_PyCFunctionFastWithKeywords _PyCFunctionFastWithKeywords
  #endif
#endif

#if CYTHON_METH_FASTCALL
  #define __Pyx_METH_FASTCALL METH_FASTCALL
  #define __Pyx_PyCFunction_FastCall __Pyx_PyCFunctionFast
  #define __Pyx_PyCFunction_FastCallWithKeywords __Pyx_PyCFunctionFastWithKeywords
#else
  #define __Pyx_METH_FASTCALL METH_VARARGS
  #define __Pyx_PyCFunction_FastCall PyCFunction
  #define __Pyx_PyCFunction_FastCallWithKeywords PyCFunctionWithKeywords
#endif

#if CYTHON_VECTORCALL
  #define __pyx_vectorcallfunc vectorcallfunc
  #define __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET  PY_VECTORCALL_ARGUMENTS_OFFSET
  #define __Pyx_PyVectorcall_NARGS(n)  PyVectorcall_NARGS((size_t)(n))
#else
  #define __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET  0
  #define __Pyx_PyVectorcall_NARGS(n)  ((Py_ssize_t)(n))
#endif

// These PyCFunction related macros get redefined in CythonFunction.c.
// We need our own copies because the inline functions in CPython have a type-check assert
// that breaks with a CyFunction in debug mode.
#if PY_VERSION_HEX >= 0x030900B1
#define __Pyx_PyCFunction_CheckExact(func)  PyCFunction_CheckExact(func)
#else
#define __Pyx_PyCFunction_CheckExact(func)  PyCFunction_Check(func)
#endif
#define __Pyx_CyOrPyCFunction_Check(func)  PyCFunction_Check(func)

#if CYTHON_COMPILING_IN_CPYTHON
#define __Pyx_CyOrPyCFunction_GET_FUNCTION(func)  (((PyCFunctionObject*)(func))->m_ml->ml_meth)
#elif !CYTHON_COMPILING_IN_LIMITED_API
// It's probably easier for non-CPythons to support PyCFunction_GET_FUNCTION() than the object struct layout.
#define __Pyx_CyOrPyCFunction_GET_FUNCTION(func)  PyCFunction_GET_FUNCTION(func)
// Unused in CYTHON_COMPILING_IN_LIMITED_API.
#endif
#if CYTHON_COMPILING_IN_CPYTHON
#define __Pyx_CyOrPyCFunction_GET_FLAGS(func)  (((PyCFunctionObject*)(func))->m_ml->ml_flags)
static CYTHON_INLINE PyObject* __Pyx_CyOrPyCFunction_GET_SELF(PyObject *func) {
    return (__Pyx_CyOrPyCFunction_GET_FLAGS(func) & METH_STATIC) ? NULL : ((PyCFunctionObject*)func)->m_self;
}
// Only used if CYTHON_COMPILING_IN_CPYTHON.
#endif
static CYTHON_INLINE int __Pyx__IsSameCFunction(PyObject *func, void (*cfunc)(void)) {
#if CYTHON_COMPILING_IN_LIMITED_API
    return PyCFunction_Check(func) && PyCFunction_GetFunction(func) == (PyCFunction) cfunc;
#else
    return PyCFunction_Check(func) && PyCFunction_GET_FUNCTION(func) == (PyCFunction) cfunc;
#endif
}
#define __Pyx_IsSameCFunction(func, cfunc)   __Pyx__IsSameCFunction(func, cfunc)

// PEP-573: PyCFunction holds reference to defining class (PyCMethodObject)
#if PY_VERSION_HEX < 0x03090000 || (CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000)
  #define __Pyx_PyType_FromModuleAndSpec(m, s, b)  ((void)m, PyType_FromSpecWithBases(s, b))
  typedef PyObject *(*__Pyx_PyCMethod)(PyObject *, PyTypeObject *, PyObject *const *, size_t, PyObject *);
#else
  #define __Pyx_PyType_FromModuleAndSpec(m, s, b)  PyType_FromModuleAndSpec(m, s, b)
  #define __Pyx_PyCMethod  PyCMethod
#endif
#ifndef METH_METHOD
  #define METH_METHOD 0x200
#endif

#if CYTHON_COMPILING_IN_PYPY && !defined(PyObject_Malloc)
  #define PyObject_Malloc(s)   PyMem_Malloc(s)
  #define PyObject_Free(p)     PyMem_Free(p)
  #define PyObject_Realloc(p)  PyMem_Realloc(p)
#endif

#if CYTHON_COMPILING_IN_LIMITED_API
  // __Pyx_PyCode_HasFreeVars isn't easily emulated in the limited API (but isn't really necessary)
  #define __Pyx_PyFrame_SetLineNumber(frame, lineno)
#elif CYTHON_COMPILING_IN_GRAAL && defined(GRAALPY_VERSION_NUM) && GRAALPY_VERSION_NUM > 0x19000000
  #define __Pyx_PyCode_HasFreeVars(co)  (PyCode_GetNumFree(co) > 0)
  #define __Pyx_PyFrame_SetLineNumber(frame, lineno) GraalPyFrame_SetLineNumber((frame), (lineno))
#elif CYTHON_COMPILING_IN_GRAAL
  // Remove when GraalPy 24 goes EOL
  #define __Pyx_PyCode_HasFreeVars(co)  (PyCode_GetNumFree(co) > 0)
  #define __Pyx_PyFrame_SetLineNumber(frame, lineno) _PyFrame_SetLineNumber((frame), (lineno))
#else
  #define __Pyx_PyCode_HasFreeVars(co)  (PyCode_GetNumFree(co) > 0)
  #define __Pyx_PyFrame_SetLineNumber(frame, lineno)  (frame)->f_lineno = (lineno)
#endif

#if CYTHON_COMPILING_IN_LIMITED_API
  #define __Pyx_PyThreadState_Current PyThreadState_Get()
#elif !CYTHON_FAST_THREAD_STATE
  #define __Pyx_PyThreadState_Current PyThreadState_GET()
#elif PY_VERSION_HEX >= 0x030d00A1
  #define __Pyx_PyThreadState_Current PyThreadState_GetUnchecked()
#else
  #define __Pyx_PyThreadState_Current _PyThreadState_UncheckedGet()
#endif

#if CYTHON_USE_MODULE_STATE
static CYTHON_INLINE void *__Pyx__PyModule_GetState(PyObject *op)
{
    void *result;

    result = PyModule_GetState(op);
    if (!result)
        Py_FatalError("Couldn't find the module state");
    return result;
}
// Define a macro with a cast because the modulestate type isn't known yet and
// is a typedef struct so impossible to forward declare
#define __Pyx_PyModule_GetState(o) ($modulestatetype_cname *)__Pyx__PyModule_GetState(o)
#else
#define __Pyx_PyModule_GetState(op) ((void)op,$modulestateglobal_cname)
#endif

// The "Try" variants may return NULL on static types with the Limited API on earlier versions
// so should be used for optimization rather than where a result is required.
#define __Pyx_PyObject_GetSlot(obj, name, func_ctype)  __Pyx_PyType_GetSlot(Py_TYPE((PyObject *) obj), name, func_ctype)
#define __Pyx_PyObject_TryGetSlot(obj, name, func_ctype) __Pyx_PyType_TryGetSlot(Py_TYPE(obj), name, func_ctype)
#define __Pyx_PyObject_GetSubSlot(obj, sub, name, func_ctype) __Pyx_PyType_GetSubSlot(Py_TYPE(obj), sub, name, func_ctype)
#define __Pyx_PyObject_TryGetSubSlot(obj, sub, name, func_ctype) __Pyx_PyType_TryGetSubSlot(Py_TYPE(obj), sub, name, func_ctype)
#if CYTHON_USE_TYPE_SLOTS
  #define __Pyx_PyType_GetSlot(type, name, func_ctype)  ((type)->name)
  #define __Pyx_PyType_TryGetSlot(type, name, func_ctype) __Pyx_PyType_GetSlot(type, name, func_ctype)
  #define __Pyx_PyType_GetSubSlot(type, sub, name, func_ctype) (((type)->sub) ? ((type)->sub->name) : NULL)
  #define __Pyx_PyType_TryGetSubSlot(type, sub, name, func_ctype) __Pyx_PyType_GetSubSlot(type, sub, name, func_ctype)
#else
  #define __Pyx_PyType_GetSlot(type, name, func_ctype)  ((func_ctype) PyType_GetSlot((type), Py_##name))
  #define __Pyx_PyType_TryGetSlot(type, name, func_ctype) \
    ((__PYX_LIMITED_VERSION_HEX >= 0x030A0000 || \
     (PyType_GetFlags(type) & Py_TPFLAGS_HEAPTYPE) || __Pyx_get_runtime_version() >= 0x030A0000) ? \
     __Pyx_PyType_GetSlot(type, name, func_ctype) : NULL)
  #define __Pyx_PyType_GetSubSlot(obj, sub, name, func_ctype) __Pyx_PyType_GetSlot(obj, name, func_ctype)
  #define __Pyx_PyType_TryGetSubSlot(obj, sub, name, func_ctype) __Pyx_PyType_TryGetSlot(obj, name, func_ctype)
#endif

#if CYTHON_COMPILING_IN_CPYTHON || defined(_PyDict_NewPresized)
#define __Pyx_PyDict_NewPresized(n)  ((n <= 8) ? PyDict_New() : _PyDict_NewPresized(n))
#else
#define __Pyx_PyDict_NewPresized(n)  PyDict_New()
#endif

#define __Pyx_PyNumber_Divide(x,y)         PyNumber_TrueDivide(x,y)
#define __Pyx_PyNumber_InPlaceDivide(x,y)  PyNumber_InPlaceTrueDivide(x,y)

#if CYTHON_COMPILING_IN_CPYTHON && CYTHON_USE_UNICODE_INTERNALS
// _PyDict_GetItem_KnownHash() existed from CPython 3.5 to 3.12, but it was
// dropping exceptions in 3.5. Since 3.6, exceptions are kept.
#define __Pyx_PyDict_GetItemStrWithError(dict, name)  _PyDict_GetItem_KnownHash(dict, name, ((PyASCIIObject *) name)->hash)
static CYTHON_INLINE PyObject * __Pyx_PyDict_GetItemStr(PyObject *dict, PyObject *name) {
    PyObject *res = __Pyx_PyDict_GetItemStrWithError(dict, name);
    if (res == NULL) PyErr_Clear();
    return res;
}
#elif !CYTHON_COMPILING_IN_PYPY || PYPY_VERSION_NUM >= 0x07020000
#define __Pyx_PyDict_GetItemStrWithError  PyDict_GetItemWithError
#define __Pyx_PyDict_GetItemStr           PyDict_GetItem
#else
static CYTHON_INLINE PyObject * __Pyx_PyDict_GetItemStrWithError(PyObject *dict, PyObject *name) {
    // This is tricky - we should return a borrowed reference but not swallow non-KeyError exceptions. 8-|
    // But: this function is only used in Py2 and older PyPys,
    // and currently only for argument parsing and other non-correctness-critical lookups
    // and we know that 'name' is an interned 'str' with pre-calculated hash value (only comparisons can fail),
    // thus, performance matters more than correctness here, especially in the "not found" case.
#if CYTHON_COMPILING_IN_PYPY
    // So we ignore any exceptions in old PyPys ...
    return PyDict_GetItem(dict, name);
#else
    // and hack together a stripped-down and modified PyDict_GetItem() in CPython 2.
    PyDictEntry *ep;
    PyDictObject *mp = (PyDictObject*) dict;
    long hash = ((PyStringObject *) name)->ob_shash;
    assert(hash != -1); /* hash values of interned strings are always initialised */
    ep = (mp->ma_lookup)(mp, name, hash);
    if (ep == NULL) {
        // error occurred
        return NULL;
    }
    // found or not found
    return ep->me_value;
#endif
}
#define __Pyx_PyDict_GetItemStr           PyDict_GetItem
#endif

/* Type slots */

#if CYTHON_USE_TYPE_SLOTS
  #define __Pyx_PyType_GetFlags(tp)   (((PyTypeObject *)tp)->tp_flags)
  #define __Pyx_PyType_HasFeature(type, feature)  ((__Pyx_PyType_GetFlags(type) & (feature)) != 0)
#else
  #define __Pyx_PyType_GetFlags(tp)   (PyType_GetFlags((PyTypeObject *)tp))
  #define __Pyx_PyType_HasFeature(type, feature)  PyType_HasFeature(type, feature)
#endif

// There is no replacement for "Py_TYPE(obj)->iternext" in the C-API.
// PyIter_Next() discards the StopIteration, unlike Python's "next()".
#define __Pyx_PyObject_GetIterNextFunc(iterator)  __Pyx_PyObject_GetSlot(iterator, tp_iternext, iternextfunc)

#if CYTHON_USE_TYPE_SPECS
// In Py3.8+, instances of heap types need to decref their type on deallocation.
// https://bugs.python.org/issue35810
#define __Pyx_PyHeapTypeObject_GC_Del(obj)  { \
    PyTypeObject *type = Py_TYPE((PyObject*)obj); \
    assert(__Pyx_PyType_HasFeature(type, Py_TPFLAGS_HEAPTYPE)); \
    PyObject_GC_Del(obj); \
    Py_DECREF(type); \
}
#else
#define __Pyx_PyHeapTypeObject_GC_Del(obj)  PyObject_GC_Del(obj)
#endif

#if CYTHON_COMPILING_IN_LIMITED_API
  #define __Pyx_PyUnicode_READY(op)       (0)
  #define __Pyx_PyUnicode_READ_CHAR(u, i) PyUnicode_ReadChar(u, i)
  #define __Pyx_PyUnicode_MAX_CHAR_VALUE(u)   ((void)u, 1114111U)
  #define __Pyx_PyUnicode_KIND(u)         ((void)u, (0))
  // __Pyx_PyUnicode_DATA() and __Pyx_PyUnicode_READ() must go together, e.g. for iteration.
  #define __Pyx_PyUnicode_DATA(u)         ((void*)u)
  #define __Pyx_PyUnicode_READ(k, d, i)   ((void)k, PyUnicode_ReadChar((PyObject*)(d), i))
  //#define __Pyx_PyUnicode_WRITE(k, d, i, ch)  /* not available */
  #define __Pyx_PyUnicode_IS_TRUE(u)      (0 != PyUnicode_GetLength(u))
#else
  #if PY_VERSION_HEX >= 0x030C0000
    // Py3.12 / PEP-623 removed wstr type unicode strings and all of the PyUnicode_READY() machinery.
    #define __Pyx_PyUnicode_READY(op)       (0)
  #else
    #define __Pyx_PyUnicode_READY(op)       (likely(PyUnicode_IS_READY(op)) ? \
                                                0 : _PyUnicode_Ready((PyObject *)(op)))
  #endif

  #define __Pyx_PyUnicode_READ_CHAR(u, i) PyUnicode_READ_CHAR(u, i)
  #define __Pyx_PyUnicode_MAX_CHAR_VALUE(u)   PyUnicode_MAX_CHAR_VALUE(u)
  #define __Pyx_PyUnicode_KIND(u)         ((int)PyUnicode_KIND(u))
  #define __Pyx_PyUnicode_DATA(u)         PyUnicode_DATA(u)
  #define __Pyx_PyUnicode_READ(k, d, i)   PyUnicode_READ(k, d, i)
  #define __Pyx_PyUnicode_WRITE(k, d, i, ch)  PyUnicode_WRITE(k, d, i, (Py_UCS4) ch)
  #if PY_VERSION_HEX >= 0x030C0000
    #define __Pyx_PyUnicode_IS_TRUE(u)      (0 != PyUnicode_GET_LENGTH(u))
  #else
    #if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX >= 0x03090000
    // Avoid calling deprecated C-API functions in Py3.9+ that PEP-623 schedules for removal in Py3.12.
    // https://www.python.org/dev/peps/pep-0623/
    #define __Pyx_PyUnicode_IS_TRUE(u)      (0 != (likely(PyUnicode_IS_READY(u)) ? PyUnicode_GET_LENGTH(u) : ((PyCompactUnicodeObject *)(u))->wstr_length))
    #else
    #define __Pyx_PyUnicode_IS_TRUE(u)      (0 != (likely(PyUnicode_IS_READY(u)) ? PyUnicode_GET_LENGTH(u) : PyUnicode_GET_SIZE(u)))
    #endif
  #endif
#endif

#if CYTHON_COMPILING_IN_PYPY
  #define __Pyx_PyUnicode_Concat(a, b)      PyNumber_Add(a, b)
  #define __Pyx_PyUnicode_ConcatSafe(a, b)  PyNumber_Add(a, b)
#else
  #define __Pyx_PyUnicode_Concat(a, b)      PyUnicode_Concat(a, b)
  #define __Pyx_PyUnicode_ConcatSafe(a, b)  ((unlikely((a) == Py_None) || unlikely((b) == Py_None)) ? \
      PyNumber_Add(a, b) : __Pyx_PyUnicode_Concat(a, b))
#endif

#if CYTHON_COMPILING_IN_PYPY
  #if !defined(PyUnicode_DecodeUnicodeEscape)
    #define PyUnicode_DecodeUnicodeEscape(s, size, errors)  PyUnicode_Decode(s, size, "unicode_escape", errors)
  #endif
  #if !defined(PyUnicode_Contains)
    #define PyUnicode_Contains(u, s)  PySequence_Contains(u, s)
  #endif
  #if !defined(PyByteArray_Check)
    #define PyByteArray_Check(obj)  PyObject_TypeCheck(obj, &PyByteArray_Type)
  #endif
  #if !defined(PyObject_Format)
    #define PyObject_Format(obj, fmt)  PyObject_CallMethod(obj, "__format__", "O", fmt)
  #endif
#endif

// ("..." % x)  must call PyNumber_Remainder() if x is a string subclass that implements "__rmod__()".
#define __Pyx_PyUnicode_FormatSafe(a, b)  ((unlikely((a) == Py_None || (PyUnicode_Check(b) && !PyUnicode_CheckExact(b)))) ? PyNumber_Remainder(a, b) : PyUnicode_Format(a, b))

#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX >= 0x030E0000
  #define __Pyx_PySequence_ListKeepNew(obj) \
    (likely(PyList_CheckExact(obj) && PyUnstable_Object_IsUniquelyReferenced(obj)) ? __Pyx_NewRef(obj) : PySequence_List(obj))
#elif CYTHON_COMPILING_IN_CPYTHON
  #define __Pyx_PySequence_ListKeepNew(obj) \
    (likely(PyList_CheckExact(obj) && Py_REFCNT(obj) == 1) ? __Pyx_NewRef(obj) : PySequence_List(obj))
#else
  #define __Pyx_PySequence_ListKeepNew(obj)  PySequence_List(obj)
#endif

#ifndef PySet_CheckExact
  #define PySet_CheckExact(obj)        __Pyx_IS_TYPE(obj, &PySet_Type)
#endif

#if PY_VERSION_HEX >= 0x030900A4
  #define __Pyx_SET_REFCNT(obj, refcnt) Py_SET_REFCNT(obj, refcnt)
  #define __Pyx_SET_SIZE(obj, size) Py_SET_SIZE(obj, size)
#else
  #define __Pyx_SET_REFCNT(obj, refcnt) Py_REFCNT(obj) = (refcnt)
  #define __Pyx_SET_SIZE(obj, size) Py_SIZE(obj) = (size)
#endif

enum __Pyx_ReferenceSharing {
  __Pyx_ReferenceSharing_DefinitelyUnique, // We created it so we know it's unshared - no need to check
  __Pyx_ReferenceSharing_OwnStrongReference,
  __Pyx_ReferenceSharing_FunctionArgument,
  __Pyx_ReferenceSharing_SharedReference, // Never trust it to be unshared because it's a global or similar
};

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && PY_VERSION_HEX >= 0x030E0000
#define __Pyx_IS_UNIQUELY_REFERENCED(o, sharing) \
    (sharing == __Pyx_ReferenceSharing_DefinitelyUnique ? 1 : \
      (sharing == __Pyx_ReferenceSharing_FunctionArgument ? PyUnstable_Object_IsUniqueReferencedTemporary(o) : \
      (sharing == __Pyx_ReferenceSharing_OwnStrongReference ? PyUnstable_Object_IsUniquelyReferenced(o) : 0)))
#elif (CYTHON_COMPILING_IN_CPYTHON && !CYTHON_COMPILING_IN_CPYTHON_FREETHREADING) || CYTHON_COMPILING_IN_LIMITED_API
#define __Pyx_IS_UNIQUELY_REFERENCED(o, sharing) (((void)sharing), Py_REFCNT(o) == 1)
#else
// On other platforms we don't trust the refcount so don't optimize based on it
#define __Pyx_IS_UNIQUELY_REFERENCED(o, sharing) (((void)o), ((void)sharing), 0)
#endif


#if CYTHON_AVOID_BORROWED_REFS || CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
  #if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
    #define __Pyx_PyList_GetItemRef(o, i) PyList_GetItemRef(o, i)
  #elif CYTHON_COMPILING_IN_LIMITED_API || !CYTHON_ASSUME_SAFE_MACROS
    #define __Pyx_PyList_GetItemRef(o, i) (likely((i) >= 0) ? PySequence_GetItem(o, i) : (PyErr_SetString(PyExc_IndexError, "list index out of range"), (PyObject*)NULL))
  #else
    #define __Pyx_PyList_GetItemRef(o, i) PySequence_ITEM(o, i)
  #endif
#elif CYTHON_COMPILING_IN_LIMITED_API || !CYTHON_ASSUME_SAFE_MACROS
  #if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
    #define __Pyx_PyList_GetItemRef(o, i) PyList_GetItemRef(o, i)
  #else
    #define __Pyx_PyList_GetItemRef(o, i) __Pyx_XNewRef(PyList_GetItem(o, i))
  #endif
#else
  #define __Pyx_PyList_GetItemRef(o, i) __Pyx_NewRef(PyList_GET_ITEM(o, i))
#endif


#if CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS && !CYTHON_COMPILING_IN_LIMITED_API && CYTHON_ASSUME_SAFE_MACROS
  #define __Pyx_PyList_GetItemRefFast(o, i, unsafe_shared) (__Pyx_IS_UNIQUELY_REFERENCED(o, unsafe_shared) ? \
    __Pyx_NewRef(PyList_GET_ITEM(o, i)) : __Pyx_PyList_GetItemRef(o, i))
#else
  #define __Pyx_PyList_GetItemRefFast(o, i, unsafe_shared) __Pyx_PyList_GetItemRef(o, i)
#endif

#if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
#define __Pyx_PyDict_GetItemRef(dict, key, result) PyDict_GetItemRef(dict, key, result)
#elif CYTHON_AVOID_BORROWED_REFS || CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
static CYTHON_INLINE int __Pyx_PyDict_GetItemRef(PyObject *dict, PyObject *key, PyObject **result) {
  *result = PyObject_GetItem(dict, key);
  if (*result == NULL) {
    if (PyErr_ExceptionMatches(PyExc_KeyError)) {
      PyErr_Clear();
      return 0;
    }
    return -1;
  }
  return 1;
}
#else
static CYTHON_INLINE int __Pyx_PyDict_GetItemRef(PyObject *dict, PyObject *key, PyObject **result) {
  *result = PyDict_GetItemWithError(dict, key);
  if (*result == NULL) {
    return PyErr_Occurred() ? -1 : 0;
  }
  Py_INCREF(*result);
  return 1;
}
#endif

// No-op macro for calling Py_VISIT() on known constants that can never participate in reference cycles.
// Users can define "CYTHON_DEBUG_VISIT_CONST=1" to help in debugging reference issues.
#if defined(CYTHON_DEBUG_VISIT_CONST) && CYTHON_DEBUG_VISIT_CONST
  #define __Pyx_VISIT_CONST(obj)  Py_VISIT(obj)
#else
  #define __Pyx_VISIT_CONST(obj)
#endif

#if CYTHON_ASSUME_SAFE_MACROS
  #define __Pyx_PySequence_ITEM(o, i) PySequence_ITEM(o, i)
  #define __Pyx_PySequence_SIZE(seq)  Py_SIZE(seq)
  #define __Pyx_PyTuple_SET_ITEM(o, i, v) (PyTuple_SET_ITEM(o, i, v), (0))
  #define __Pyx_PyTuple_GET_ITEM(o, i) PyTuple_GET_ITEM(o, i)
  #define __Pyx_PyList_SET_ITEM(o, i, v) (PyList_SET_ITEM(o, i, v), (0))
  #define __Pyx_PyList_GET_ITEM(o, i) PyList_GET_ITEM(o, i)
#else
  #define __Pyx_PySequence_ITEM(o, i) PySequence_GetItem(o, i)
  // NOTE: might fail with exception => check for -1
  #define __Pyx_PySequence_SIZE(seq)  PySequence_Size(seq)
  // NOTE: this doesn't leak a reference to whatever is at o[i]
  #define __Pyx_PyTuple_SET_ITEM(o, i, v) PyTuple_SetItem(o, i, v)
  #define __Pyx_PyTuple_GET_ITEM(o, i) PyTuple_GetItem(o, i)
  #define __Pyx_PyList_SET_ITEM(o, i, v) PyList_SetItem(o, i, v)
  #define __Pyx_PyList_GET_ITEM(o, i) PyList_GetItem(o, i)
#endif

#if CYTHON_ASSUME_SAFE_SIZE
  #define __Pyx_PyTuple_GET_SIZE(o) PyTuple_GET_SIZE(o)
  #define __Pyx_PyList_GET_SIZE(o) PyList_GET_SIZE(o)
  #define __Pyx_PySet_GET_SIZE(o) PySet_GET_SIZE(o)
  #define __Pyx_PyBytes_GET_SIZE(o) PyBytes_GET_SIZE(o)
  #define __Pyx_PyByteArray_GET_SIZE(o) PyByteArray_GET_SIZE(o)
  #define __Pyx_PyUnicode_GET_LENGTH(o) PyUnicode_GET_LENGTH(o)
#else
  // These all need exception checks for -1.
  #define __Pyx_PyTuple_GET_SIZE(o) PyTuple_Size(o)
  #define __Pyx_PyList_GET_SIZE(o) PyList_Size(o)
  #define __Pyx_PySet_GET_SIZE(o) PySet_Size(o)
  #define __Pyx_PyBytes_GET_SIZE(o) PyBytes_Size(o)
  #define __Pyx_PyByteArray_GET_SIZE(o) PyByteArray_Size(o)
  #define __Pyx_PyUnicode_GET_LENGTH(o) PyUnicode_GetLength(o)
#endif

#if CYTHON_COMPILING_IN_PYPY && !defined(PyUnicode_InternFromString)
  #define PyUnicode_InternFromString(s) PyUnicode_FromString(s)
#endif

#define __Pyx_PyLong_FromHash_t PyLong_FromSsize_t
#define __Pyx_PyLong_AsHash_t   __Pyx_PyIndex_AsSsize_t


// backport of PyAsyncMethods from Py3.10 to older Py3.x versions
#if __PYX_LIMITED_VERSION_HEX >= 0x030A0000
    #define __Pyx_PySendResult PySendResult
#else
    typedef enum {
        PYGEN_RETURN = 0,
        PYGEN_ERROR = -1,
        PYGEN_NEXT = 1,
    } __Pyx_PySendResult;
#endif

#if CYTHON_COMPILING_IN_LIMITED_API || PY_VERSION_HEX < 0x030A00A3
  typedef __Pyx_PySendResult (*__Pyx_pyiter_sendfunc)(PyObject *iter, PyObject *value, PyObject **result);
#else
  #define __Pyx_pyiter_sendfunc sendfunc
#endif

// "Py_am_send" requires Py3.10 when using type specs (which utility code types do).
#if !CYTHON_USE_AM_SEND
#define __PYX_HAS_PY_AM_SEND 0
#elif __PYX_LIMITED_VERSION_HEX >= 0x030A0000
#define __PYX_HAS_PY_AM_SEND 1
#else
#define __PYX_HAS_PY_AM_SEND 2  // our own backported implementation
#endif

#if __PYX_HAS_PY_AM_SEND < 2
    #define __Pyx_PyAsyncMethodsStruct PyAsyncMethods
#else
    // PyAsyncMethods in Py<3.10 lacks "am_send"
    typedef struct {
        unaryfunc am_await;
        unaryfunc am_aiter;
        unaryfunc am_anext;
        __Pyx_pyiter_sendfunc am_send;
    } __Pyx_PyAsyncMethodsStruct;

    #define __Pyx_SlotTpAsAsync(s) ((PyAsyncMethods*)(s))
#endif

// Use a flag in Py < 3.10 to mark coroutines that have the "am_send" field.
#if CYTHON_USE_AM_SEND && PY_VERSION_HEX < 0x030A00F0
    #define __Pyx_TPFLAGS_HAVE_AM_SEND (1UL << 21)
#else
    #define __Pyx_TPFLAGS_HAVE_AM_SEND (0)
#endif

#if PY_VERSION_HEX >= 0x03090000
#define __Pyx_PyInterpreterState_Get() PyInterpreterState_Get()
#else
#define __Pyx_PyInterpreterState_Get() PyThreadState_Get()->interp
#endif

#if CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX < 0x030A0000
// PyMem_Calloc *is* in the Stable ABI in all the Limited API versions we care about.
// However, it is omitted from the Python headers which means that C incorrectly
// assumes it returns an int (and generates dubious code on based on that assumption).
// Therefore, copy the prototype.
#ifdef __cplusplus
extern "C"
#endif
PyAPI_FUNC(void *) PyMem_Calloc(size_t nelem, size_t elsize); /* proto */
#endif

#if CYTHON_COMPILING_IN_LIMITED_API
// returns 1 for success and 0 for failure to enable it to be chained in an &&
static int __Pyx_init_co_variable(PyObject *inspect, const char* name, int *write_to) {
    int value;
    PyObject *py_value = PyObject_GetAttrString(inspect, name);
    if (!py_value) return 0;
    // There's a small chance of overflow here, but it'd only happen if inspect was set up wrongly.
    value = (int) PyLong_AsLong(py_value);
    Py_DECREF(py_value);
    *write_to = value;
    return value != -1 || !PyErr_Occurred();
}

// Returns 0 on success and -1 on failure for normal error handling
static int __Pyx_init_co_variables(void) {
    PyObject *inspect;
    int result;
    inspect = PyImport_ImportModule("inspect");

    result =
#if !defined(CO_OPTIMIZED)
        __Pyx_init_co_variable(inspect, "CO_OPTIMIZED", &CO_OPTIMIZED) &&
#endif
#if !defined(CO_NEWLOCALS)
        __Pyx_init_co_variable(inspect, "CO_NEWLOCALS", &CO_NEWLOCALS) &&
#endif
#if !defined(CO_VARARGS)
        __Pyx_init_co_variable(inspect, "CO_VARARGS", &CO_VARARGS) &&
#endif
#if !defined(CO_VARKEYWORDS)
        __Pyx_init_co_variable(inspect, "CO_VARKEYWORDS", &CO_VARKEYWORDS) &&
#endif
#if !defined(CO_ASYNC_GENERATOR)
        __Pyx_init_co_variable(inspect, "CO_ASYNC_GENERATOR", &CO_ASYNC_GENERATOR) &&
#endif
#if !defined(CO_GENERATOR)
        __Pyx_init_co_variable(inspect, "CO_GENERATOR", &CO_GENERATOR) &&
#endif
#if !defined(CO_COROUTINE)
        __Pyx_init_co_variable(inspect, "CO_COROUTINE", &CO_COROUTINE) &&
#endif
        1;

    Py_DECREF(inspect);
    return result ? 0 : -1;
}
#else
static int __Pyx_init_co_variables(void) {
    return 0;  // It's a limited API-only feature
}
#endif

/////////////// CythonABIVersion.proto ///////////////
//@proto_block: module_declarations
// This needs to go after the utility code 'proto' section but before user code and utility impl.

#if CYTHON_COMPILING_IN_LIMITED_API
    // The limited API makes some significant changes to data structures, so we don't
    // want to share the implementations compiled with and without the limited API.
    #if CYTHON_METH_FASTCALL
        #define __PYX_FASTCALL_ABI_SUFFIX  "_fastcall"
    #else
        #define __PYX_FASTCALL_ABI_SUFFIX
    #endif

    #define __PYX_LIMITED_ABI_SUFFIX "limited" __PYX_FASTCALL_ABI_SUFFIX __PYX_AM_SEND_ABI_SUFFIX
#else
    #define __PYX_LIMITED_ABI_SUFFIX
#endif

#if __PYX_HAS_PY_AM_SEND == 1
    #define __PYX_AM_SEND_ABI_SUFFIX
#elif __PYX_HAS_PY_AM_SEND == 2
    #define __PYX_AM_SEND_ABI_SUFFIX "amsendbackport"
#else
    #define __PYX_AM_SEND_ABI_SUFFIX "noamsend"
#endif

#ifndef __PYX_MONITORING_ABI_SUFFIX
    #define __PYX_MONITORING_ABI_SUFFIX
#endif

#if CYTHON_USE_TP_FINALIZE
    #define __PYX_TP_FINALIZE_ABI_SUFFIX
#else
    // affects destruction of async generator/coroutines
    #define __PYX_TP_FINALIZE_ABI_SUFFIX "nofinalize"
#endif

#if CYTHON_USE_FREELISTS || !defined(__Pyx_AsyncGen_USED)
    #define __PYX_FREELISTS_ABI_SUFFIX
#else
    // affects allocation/deallocation of async generator objects.
    #define __PYX_FREELISTS_ABI_SUFFIX "nofreelists"
#endif

#define CYTHON_ABI  __PYX_ABI_VERSION __PYX_LIMITED_ABI_SUFFIX __PYX_MONITORING_ABI_SUFFIX __PYX_TP_FINALIZE_ABI_SUFFIX __PYX_FREELISTS_ABI_SUFFIX __PYX_AM_SEND_ABI_SUFFIX

#define __PYX_ABI_MODULE_NAME "_cython_" CYTHON_ABI
#define __PYX_TYPE_MODULE_PREFIX __PYX_ABI_MODULE_NAME "."

/////////////// PythonCompatibility.init ///////////////

if (likely(__Pyx_init_co_variables() == 0)); else


/////////////// IncludeStructmemberH.proto ///////////////
//@proto_block: utility_code_proto_before_types

#include <structmember.h>


/////////////// SmallCodeConfig ///////////////

#ifndef CYTHON_SMALL_CODE
#if defined(__clang__)
    #define CYTHON_SMALL_CODE
#elif defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3))
    #define CYTHON_SMALL_CODE __attribute__((cold))
#else
    #define CYTHON_SMALL_CODE
#endif
#endif


/////////////// PyModInitFuncType ///////////////

#ifndef CYTHON_NO_PYINIT_EXPORT
  #define __Pyx_PyMODINIT_FUNC PyMODINIT_FUNC
#else
  // define this to PyObject * manually because PyMODINIT_FUNC adds __declspec(dllexport) to it's definition.
  #ifdef __cplusplus
  #define __Pyx_PyMODINIT_FUNC extern "C" PyObject *
  #else
  #define __Pyx_PyMODINIT_FUNC PyObject *
  #endif
#endif


/////////////// FastTypeChecks.proto ///////////////

#if CYTHON_COMPILING_IN_CPYTHON
#define __Pyx_TypeCheck(obj, type) __Pyx_IsSubtype(Py_TYPE(obj), (PyTypeObject *)type)
#define __Pyx_TypeCheck2(obj, type1, type2) __Pyx_IsAnySubtype2(Py_TYPE(obj), (PyTypeObject *)type1, (PyTypeObject *)type2)
static CYTHON_INLINE int __Pyx_IsSubtype(PyTypeObject *a, PyTypeObject *b);/*proto*/
static CYTHON_INLINE int __Pyx_IsAnySubtype2(PyTypeObject *cls, PyTypeObject *a, PyTypeObject *b);/*proto*/
static CYTHON_INLINE int __Pyx_PyErr_GivenExceptionMatches(PyObject *err, PyObject *type);/*proto*/
static CYTHON_INLINE int __Pyx_PyErr_GivenExceptionMatches2(PyObject *err, PyObject *type1, PyObject *type2);/*proto*/
#else
#define __Pyx_TypeCheck(obj, type) PyObject_TypeCheck(obj, (PyTypeObject *)type)
#define __Pyx_TypeCheck2(obj, type1, type2) (PyObject_TypeCheck(obj, (PyTypeObject *)type1) || PyObject_TypeCheck(obj, (PyTypeObject *)type2))
#define __Pyx_PyErr_GivenExceptionMatches(err, type) PyErr_GivenExceptionMatches(err, type)
static CYTHON_INLINE int __Pyx_PyErr_GivenExceptionMatches2(PyObject *err, PyObject *type1, PyObject *type2) {
    return PyErr_GivenExceptionMatches(err, type1) || PyErr_GivenExceptionMatches(err, type2);
}
#endif
#define __Pyx_PyErr_ExceptionMatches2(err1, err2)  __Pyx_PyErr_GivenExceptionMatches2(__Pyx_PyErr_CurrentExceptionType(), err1, err2)

#define __Pyx_PyException_Check(obj) __Pyx_TypeCheck(obj, PyExc_Exception)
#ifdef PyExceptionInstance_Check
  #define __Pyx_PyBaseException_Check(obj) PyExceptionInstance_Check(obj)
#else
  #define __Pyx_PyBaseException_Check(obj) __Pyx_TypeCheck(obj, PyExc_BaseException)
#endif


/////////////// FastTypeChecks ///////////////
//@requires: Exceptions.c::PyThreadStateGet
//@requires: Exceptions.c::PyErrFetchRestore

#if CYTHON_COMPILING_IN_CPYTHON
static int __Pyx_InBases(PyTypeObject *a, PyTypeObject *b) {
    while (a) {
        a = __Pyx_PyType_GetSlot(a, tp_base, PyTypeObject*);
        if (a == b)
            return 1;
    }
    return b == &PyBaseObject_Type;
}

static CYTHON_INLINE int __Pyx_IsSubtype(PyTypeObject *a, PyTypeObject *b) {
    PyObject *mro;
    if (a == b) return 1;
    mro = a->tp_mro;
    if (likely(mro)) {
        Py_ssize_t i, n;
        n = PyTuple_GET_SIZE(mro);
        for (i = 0; i < n; i++) {
            if (PyTuple_GET_ITEM(mro, i) == (PyObject *)b)
                return 1;
        }
        return 0;
    }
    // should only get here for incompletely initialised types, i.e. never under normal usage patterns
    return __Pyx_InBases(a, b);
}

static CYTHON_INLINE int __Pyx_IsAnySubtype2(PyTypeObject *cls, PyTypeObject *a, PyTypeObject *b) {
    PyObject *mro;
    if (cls == a || cls == b) return 1;
    mro = cls->tp_mro;
    if (likely(mro)) {
        Py_ssize_t i, n;
        n = PyTuple_GET_SIZE(mro);
        for (i = 0; i < n; i++) {
            PyObject *base = PyTuple_GET_ITEM(mro, i);
            if (base == (PyObject *)a || base == (PyObject *)b)
                return 1;
        }
        return 0;
    }
    // should only get here for incompletely initialised types, i.e. never under normal usage patterns
    return __Pyx_InBases(cls, a) || __Pyx_InBases(cls, b);
}


static CYTHON_INLINE int __Pyx_inner_PyErr_GivenExceptionMatches2(PyObject *err, PyObject* exc_type1, PyObject *exc_type2) {
    if (exc_type1) {
        return __Pyx_IsAnySubtype2((PyTypeObject*)err, (PyTypeObject*)exc_type1, (PyTypeObject*)exc_type2);
    } else {
        return __Pyx_IsSubtype((PyTypeObject*)err, (PyTypeObject*)exc_type2);
    }
}

// so far, we only call PyErr_GivenExceptionMatches() with an exception type (not instance) as first argument
// => optimise for that case

static int __Pyx_PyErr_GivenExceptionMatchesTuple(PyObject *exc_type, PyObject *tuple) {
    Py_ssize_t i, n;
    assert(PyExceptionClass_Check(exc_type));
    n = PyTuple_GET_SIZE(tuple);
    // the tight subtype checking in Py3 allows faster out-of-order comparison
    for (i=0; i<n; i++) {
        if (exc_type == PyTuple_GET_ITEM(tuple, i)) return 1;
    }
    for (i=0; i<n; i++) {
        PyObject *t = PyTuple_GET_ITEM(tuple, i);
        if (likely(PyExceptionClass_Check(t))) {
            if (__Pyx_inner_PyErr_GivenExceptionMatches2(exc_type, NULL, t)) return 1;
        } else {
            // FIXME: Py3: PyErr_SetString(PyExc_TypeError, "catching classes that do not inherit from BaseException is not allowed");
        }
    }
    return 0;
}

static CYTHON_INLINE int __Pyx_PyErr_GivenExceptionMatches(PyObject *err, PyObject* exc_type) {
    if (likely(err == exc_type)) return 1;
    if (likely(PyExceptionClass_Check(err))) {
        if (likely(PyExceptionClass_Check(exc_type))) {
            return __Pyx_inner_PyErr_GivenExceptionMatches2(err, NULL, exc_type);
        } else if (likely(PyTuple_Check(exc_type))) {
            return __Pyx_PyErr_GivenExceptionMatchesTuple(err, exc_type);
        } else {
            // FIXME: Py3: PyErr_SetString(PyExc_TypeError, "catching classes that do not inherit from BaseException is not allowed");
        }
    }
    return PyErr_GivenExceptionMatches(err, exc_type);
}

static CYTHON_INLINE int __Pyx_PyErr_GivenExceptionMatches2(PyObject *err, PyObject *exc_type1, PyObject *exc_type2) {
    // Only used internally with known exception types => pure safety check assertions.
    assert(PyExceptionClass_Check(exc_type1));
    assert(PyExceptionClass_Check(exc_type2));
    if (likely(err == exc_type1 || err == exc_type2)) return 1;
    if (likely(PyExceptionClass_Check(err))) {
        return __Pyx_inner_PyErr_GivenExceptionMatches2(err, exc_type1, exc_type2);
    }
    return (PyErr_GivenExceptionMatches(err, exc_type1) || PyErr_GivenExceptionMatches(err, exc_type2));
}

#endif


/////////////// MathInitCode ///////////////

#if defined(_WIN32) || defined(WIN32) || defined(MS_WINDOWS)
  #ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
  #endif
#endif
#include <math.h>

#if defined(__CYGWIN__) && defined(_LDBL_EQ_DBL)
#define __Pyx_truncl trunc
#else
#define __Pyx_truncl truncl
#endif

/////////////// ForceInitThreads.proto ///////////////
//@proto_block: utility_code_proto_before_types

#ifndef __PYX_FORCE_INIT_THREADS
  #define __PYX_FORCE_INIT_THREADS 0
#endif


/////////////// ModuleCreationPEP489 ///////////////
//@substitute: naming

#if CYTHON_COMPILING_IN_LIMITED_API && (__PYX_LIMITED_VERSION_HEX < 0x03090000 \
      || ((defined(_WIN32) || defined(WIN32) || defined(MS_WINDOWS)) && __PYX_LIMITED_VERSION_HEX < 0x030A0000))
// Probably won't work before 3.8, but we don't use restricted API to find that out.
static PY_INT64_T __Pyx_GetCurrentInterpreterId(void) {
    {
        PyObject *module = PyImport_ImportModule("_interpreters"); // 3.13+ I think
        if (!module) {
            PyErr_Clear(); // just try the 3.8-3.12 version
            module = PyImport_ImportModule("_xxsubinterpreters");
            if (!module) goto bad;
        }
        PyObject *current = PyObject_CallMethod(module, "get_current", NULL);
        Py_DECREF(module);
        if (!current) goto bad;
        if (PyTuple_Check(current)) {
            // I think 3.13+ returns a tuple of (ID, whence),
            // but it's obviously a private module so the API changes a bit.
            PyObject *new_current = PySequence_GetItem(current, 0);
            Py_DECREF(current);
            current = new_current;
            if (!new_current) goto bad;
        }
        long long as_c_int = PyLong_AsLongLong(current);
        Py_DECREF(current);
        return as_c_int;
    }
  bad:
    PySys_WriteStderr("__Pyx_GetCurrentInterpreterId failed. Try setting the C define CYTHON_PEP489_MULTI_PHASE_INIT=0\n");
    return -1;
}
#endif

//#if CYTHON_PEP489_MULTI_PHASE_INIT
#if !CYTHON_USE_MODULE_STATE
static CYTHON_SMALL_CODE int __Pyx_check_single_interpreter(void) {
    static PY_INT64_T main_interpreter_id = -1;
#if CYTHON_COMPILING_IN_GRAAL && defined(GRAALPY_VERSION_NUM) && GRAALPY_VERSION_NUM > 0x19000000
    PY_INT64_T current_id = GraalPyInterpreterState_GetIDFromThreadState(PyThreadState_Get());
#elif CYTHON_COMPILING_IN_GRAAL
    PY_INT64_T current_id = PyInterpreterState_GetIDFromThreadState(PyThreadState_Get());
#elif CYTHON_COMPILING_IN_LIMITED_API && (__PYX_LIMITED_VERSION_HEX < 0x03090000 \
      || ((defined(_WIN32) || defined(WIN32) || defined(MS_WINDOWS)) && __PYX_LIMITED_VERSION_HEX < 0x030A0000))
    // Although PyInterpreterState_Get is part of the Stable ABI from Python 3.9,
    // it was somehow omitted from the Windows dll in that version.
    PY_INT64_T current_id = __Pyx_GetCurrentInterpreterId();
#elif CYTHON_COMPILING_IN_LIMITED_API
    PY_INT64_T current_id = PyInterpreterState_GetID(PyInterpreterState_Get());
#else
    PY_INT64_T current_id = PyInterpreterState_GetID(PyThreadState_Get()->interp);
#endif
    if (unlikely(current_id == -1)) {
        return -1;
    }
    if (main_interpreter_id == -1) {
        main_interpreter_id = current_id;
        return 0;
    } else if (unlikely(main_interpreter_id != current_id)) {
        PyErr_SetString(
            PyExc_ImportError,
            "Interpreter change detected - this module can only be loaded into one interpreter per process.");
        return -1;
    }
    return 0;
}
#endif

static CYTHON_SMALL_CODE int __Pyx_copy_spec_to_module(PyObject *spec, PyObject *moddict, const char* from_name, const char* to_name, int allow_none)
{
    PyObject *value = PyObject_GetAttrString(spec, from_name);
    int result = 0;
    if (likely(value)) {
        if (allow_none || value != Py_None) {
            result = PyDict_SetItemString(moddict, to_name, value);
        }
        Py_DECREF(value);
    } else if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
        PyErr_Clear();
    } else {
        result = -1;
    }
    return result;
}

static CYTHON_SMALL_CODE PyObject* ${pymodule_create_func_cname}(PyObject *spec, PyModuleDef *def) {
    PyObject *module = NULL, *moddict, *modname;
    CYTHON_UNUSED_VAR(def);

    #if !CYTHON_USE_MODULE_STATE
    // For now, we only have exactly one module instance.
    if (__Pyx_check_single_interpreter())
        return NULL;
    #endif
    if (${module_cname})
        return __Pyx_NewRef(${module_cname});

    modname = PyObject_GetAttrString(spec, "name");
    if (unlikely(!modname)) goto bad;

    module = PyModule_NewObject(modname);
    Py_DECREF(modname);
    if (unlikely(!module)) goto bad;

    moddict = PyModule_GetDict(module);
    if (unlikely(!moddict)) goto bad;
    // moddict is a borrowed reference

    if (unlikely(__Pyx_copy_spec_to_module(spec, moddict, "loader", "__loader__", 1) < 0)) goto bad;
    if (unlikely(__Pyx_copy_spec_to_module(spec, moddict, "origin", "__file__", 1) < 0)) goto bad;
    if (unlikely(__Pyx_copy_spec_to_module(spec, moddict, "parent", "__package__", 1) < 0)) goto bad;
    if (unlikely(__Pyx_copy_spec_to_module(spec, moddict, "submodule_search_locations", "__path__", 0) < 0)) goto bad;

    return module;
bad:
    Py_XDECREF(module);
    return NULL;
}
//#endif


/////////////// CodeObjectCache.proto ///////////////
//@requires: Synchronization.c::Atomics

#if CYTHON_COMPILING_IN_LIMITED_API
typedef PyObject __Pyx_CachedCodeObjectType;
#else
typedef PyCodeObject __Pyx_CachedCodeObjectType;
#endif

typedef struct {
    __Pyx_CachedCodeObjectType* code_object;
    int code_line;
} __Pyx_CodeObjectCacheEntry;

struct __Pyx_CodeObjectCache {
    int count;
    int max_count;
    __Pyx_CodeObjectCacheEntry* entries;
  #if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    // 0 for none, +ve for readers, -ve for writers.
    //
    __pyx_atomic_int_type accessor_count;
  #endif
};

static int __pyx_bisect_code_objects(__Pyx_CodeObjectCacheEntry* entries, int count, int code_line);

static __Pyx_CachedCodeObjectType *__pyx_find_code_object(int code_line);
static void __pyx_insert_code_object(int code_line, __Pyx_CachedCodeObjectType* code_object);

/////////////// CodeObjectCache.module_state_decls ////////////////

struct __Pyx_CodeObjectCache __pyx_code_cache;

/////////////// CodeObjectCache ///////////////
// Note that errors are simply ignored in the code below.
// This is just a cache, if a lookup or insertion fails - so what?

static int __pyx_bisect_code_objects(__Pyx_CodeObjectCacheEntry* entries, int count, int code_line) {
    int start = 0, mid = 0, end = count - 1;
    if (end >= 0 && code_line > entries[end].code_line) {
        return count;
    }
    while (start < end) {
        mid = start + (end - start) / 2;
        if (code_line < entries[mid].code_line) {
            end = mid;
        } else if (code_line > entries[mid].code_line) {
             start = mid + 1;
        } else {
            return mid;
        }
    }
    if (code_line <= entries[mid].code_line) {
        return mid;
    } else {
        return mid + 1;
    }
}

static __Pyx_CachedCodeObjectType *__pyx__find_code_object(struct __Pyx_CodeObjectCache *code_cache, int code_line) {
    __Pyx_CachedCodeObjectType* code_object;
    int pos;
    if (unlikely(!code_line) || unlikely(!code_cache->entries)) {
        return NULL;
    }
    pos = __pyx_bisect_code_objects(code_cache->entries, code_cache->count, code_line);
    if (unlikely(pos >= code_cache->count) || unlikely(code_cache->entries[pos].code_line != code_line)) {
        return NULL;
    }
    code_object = code_cache->entries[pos].code_object;
    Py_INCREF(code_object);
    return code_object;
}

static __Pyx_CachedCodeObjectType *__pyx_find_code_object(int code_line) {
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && !CYTHON_ATOMICS
    (void)__pyx__find_code_object;
    return NULL; // Most implementation should have atomics. But otherwise, don't make it thread-safe, just miss.
#else
    struct __Pyx_CodeObjectCache *code_cache = &CGLOBAL(__pyx_code_cache);
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    __pyx_nonatomic_int_type old_count = __pyx_atomic_incr_acq_rel(&code_cache->accessor_count);
    if (old_count < 0) {
        // It's being written so currently unreadable.
        __pyx_atomic_decr_acq_rel(&code_cache->accessor_count);
        return NULL;
    }
#endif
    __Pyx_CachedCodeObjectType *result = __pyx__find_code_object(code_cache, code_line);
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    __pyx_atomic_decr_acq_rel(&code_cache->accessor_count);
#endif
    return result;
#endif
}


static void __pyx__insert_code_object(struct __Pyx_CodeObjectCache *code_cache, int code_line, __Pyx_CachedCodeObjectType* code_object)
{
    int pos, i;
    __Pyx_CodeObjectCacheEntry* entries = code_cache->entries;
    if (unlikely(!code_line)) {
        return;
    }
    if (unlikely(!entries)) {
        entries = (__Pyx_CodeObjectCacheEntry*)PyMem_Malloc(64*sizeof(__Pyx_CodeObjectCacheEntry));
        if (likely(entries)) {
            code_cache->entries = entries;
            code_cache->max_count = 64;
            code_cache->count = 1;
            entries[0].code_line = code_line;
            entries[0].code_object = code_object;
            Py_INCREF(code_object);
        }
        return;
    }
    pos = __pyx_bisect_code_objects(code_cache->entries, code_cache->count, code_line);
    if ((pos < code_cache->count) && unlikely(code_cache->entries[pos].code_line == code_line)) {
        __Pyx_CachedCodeObjectType* tmp = entries[pos].code_object;
        entries[pos].code_object = code_object;
        Py_INCREF(code_object);
        Py_DECREF(tmp);
        return;
    }
    if (code_cache->count == code_cache->max_count) {
        int new_max = code_cache->max_count + 64;
        entries = (__Pyx_CodeObjectCacheEntry*)PyMem_Realloc(
            code_cache->entries, ((size_t)new_max) * sizeof(__Pyx_CodeObjectCacheEntry));
        if (unlikely(!entries)) {
            return;
        }
        code_cache->entries = entries;
        code_cache->max_count = new_max;
    }
    for (i=code_cache->count; i>pos; i--) {
        entries[i] = entries[i-1];
    }
    entries[pos].code_line = code_line;
    entries[pos].code_object = code_object;
    code_cache->count++;
    Py_INCREF(code_object);
}

static void __pyx_insert_code_object(int code_line, __Pyx_CachedCodeObjectType* code_object) {
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && !CYTHON_ATOMICS
    (void)__pyx__insert_code_object;
    return; // Most implementation should have atomics. But otherwise, don't make it thread-safe, just fail.
#else
    struct __Pyx_CodeObjectCache *code_cache = &CGLOBAL(__pyx_code_cache);
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    __pyx_nonatomic_int_type expected = 0;
    if (!__pyx_atomic_int_cmp_exchange(&code_cache->accessor_count, &expected, INT_MIN)) {
        // it's being written or read, Either way we can't do anything
        return;
    }
#endif
    __pyx__insert_code_object(code_cache, code_line, code_object);
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    __pyx_atomic_sub(&code_cache->accessor_count, INT_MIN);
#endif
#endif
}

/////////////// CodeObjectCache.cleanup ///////////////

{
  struct __Pyx_CodeObjectCache *code_cache = &CGLOBAL(__pyx_code_cache);
  if (code_cache->entries) {
      __Pyx_CodeObjectCacheEntry* entries = code_cache->entries;
      int i, count = code_cache->count;
      code_cache->count = 0;
      code_cache->max_count = 0;
      code_cache->entries = NULL;
      for (i=0; i<count; i++) {
          Py_DECREF(entries[i].code_object);
      }
      PyMem_Free(entries);
  }
}

/////////////// GetRuntimeVersion.proto ///////////////

#if __PYX_LIMITED_VERSION_HEX < 0x030b0000
static unsigned long __Pyx_cached_runtime_version = 0;

static void __Pyx_init_runtime_version(void);
#else
#define __Pyx_init_runtime_version()
#endif

// Does not require the GIL to be held
static unsigned long __Pyx_get_runtime_version(void);

/////////////// GetRuntimeVersion ///////////////

#if __PYX_LIMITED_VERSION_HEX < 0x030b0000
void __Pyx_init_runtime_version(void) {
    if (__Pyx_cached_runtime_version == 0) {
        const char* rt_version = Py_GetVersion();
        unsigned long version = 0;
        unsigned long factor = 0x01000000UL;
        unsigned int digit = 0;
        int i = 0;
        while (factor) {
            while ('0' <= rt_version[i] && rt_version[i] <= '9') {
                digit = digit * 10 + (unsigned int) (rt_version[i] - '0');
                ++i;
            }
            version += factor * digit;
            if (rt_version[i] != '.')
                break;
            digit = 0;
            factor >>= 8;
            ++i;
        }
        __Pyx_cached_runtime_version = version;
    }
}
#endif

static unsigned long __Pyx_get_runtime_version(void) {
    // We will probably never need the alpha/beta status, so avoid the complexity to parse it.
#if __PYX_LIMITED_VERSION_HEX >= 0x030b0000
    return Py_Version & ~0xFFUL;
#else
    return __Pyx_cached_runtime_version;
#endif
}

/////////////// CheckBinaryVersion.proto ///////////////

static int __Pyx_check_binary_version(unsigned long ct_version, unsigned long rt_version, int allow_newer);

/////////////// CheckBinaryVersion ///////////////

static int __Pyx_check_binary_version(unsigned long ct_version, unsigned long rt_version, int allow_newer) {
    // runtime version is: -1 => older, 0 => equal, 1 => newer
    const unsigned long MAJOR_MINOR = 0xFFFF0000UL;
    if ((rt_version & MAJOR_MINOR) == (ct_version & MAJOR_MINOR))
        return 0;
    if (likely(allow_newer && (rt_version & MAJOR_MINOR) > (ct_version & MAJOR_MINOR)))
        return 1;

    {
        char message[200];
        PyOS_snprintf(message, sizeof(message),
                      "compile time Python version %d.%d "
                      "of module '%.100s' "
                      "%s "
                      "runtime version %d.%d",
                       (int) (ct_version >> 24), (int) ((ct_version >> 16) & 0xFF),
                       __Pyx_MODULE_NAME,
                       (allow_newer) ? "was newer than" : "does not match",
                       (int) (rt_version >> 24), (int) ((rt_version >> 16) & 0xFF)
       );
        // returns 0 or -1
        return PyErr_WarnEx(NULL, message, 1);
    }
}

/////////////// IsLittleEndian.proto ///////////////

static CYTHON_INLINE int __Pyx_Is_Little_Endian(void);

/////////////// IsLittleEndian ///////////////

static CYTHON_INLINE int __Pyx_Is_Little_Endian(void)
{
  union {
    uint32_t u32;
    uint8_t u8[4];
  } S;
  S.u32 = 0x01020304;
  return S.u8[0] == 4;
}

/////////////// Refnanny.proto ///////////////

#ifndef CYTHON_REFNANNY
  #define CYTHON_REFNANNY 0
#endif

#if CYTHON_REFNANNY
  typedef struct {
    void (*INCREF)(void*, PyObject*, Py_ssize_t);
    void (*DECREF)(void*, PyObject*, Py_ssize_t);
    void (*GOTREF)(void*, PyObject*, Py_ssize_t);
    void (*GIVEREF)(void*, PyObject*, Py_ssize_t);
    void* (*SetupContext)(const char*, Py_ssize_t, const char*);
    void (*FinishContext)(void**);
  } __Pyx_RefNannyAPIStruct;
  static __Pyx_RefNannyAPIStruct *__Pyx_RefNanny = NULL;
  static __Pyx_RefNannyAPIStruct *__Pyx_RefNannyImportAPI(const char *modname); /*proto*/
  #define __Pyx_RefNannyDeclarations void *__pyx_refnanny = NULL;
  #define __Pyx_RefNannySetupContext(name, acquire_gil) \
          if (acquire_gil) { \
              PyGILState_STATE __pyx_gilstate_save = PyGILState_Ensure(); \
              __pyx_refnanny = __Pyx_RefNanny->SetupContext((name), (__LINE__), (__FILE__)); \
              PyGILState_Release(__pyx_gilstate_save); \
          } else { \
              __pyx_refnanny = __Pyx_RefNanny->SetupContext((name), (__LINE__), (__FILE__)); \
          }
  #define __Pyx_RefNannyFinishContextNogil() { \
              PyGILState_STATE __pyx_gilstate_save = PyGILState_Ensure(); \
              __Pyx_RefNannyFinishContext(); \
              PyGILState_Release(__pyx_gilstate_save); \
          }
  #define __Pyx_RefNannyFinishContextNogil() { \
              PyGILState_STATE __pyx_gilstate_save = PyGILState_Ensure(); \
              __Pyx_RefNannyFinishContext(); \
              PyGILState_Release(__pyx_gilstate_save); \
          }
  #define __Pyx_RefNannyFinishContext() \
          __Pyx_RefNanny->FinishContext(&__pyx_refnanny)
  #define __Pyx_INCREF(r)  __Pyx_RefNanny->INCREF(__pyx_refnanny, (PyObject *)(r), (__LINE__))
  #define __Pyx_DECREF(r)  __Pyx_RefNanny->DECREF(__pyx_refnanny, (PyObject *)(r), (__LINE__))
  #define __Pyx_GOTREF(r)  __Pyx_RefNanny->GOTREF(__pyx_refnanny, (PyObject *)(r), (__LINE__))
  #define __Pyx_GIVEREF(r) __Pyx_RefNanny->GIVEREF(__pyx_refnanny, (PyObject *)(r), (__LINE__))
  #define __Pyx_XINCREF(r)  do { if((r) == NULL); else {__Pyx_INCREF(r); }} while(0)
  #define __Pyx_XDECREF(r)  do { if((r) == NULL); else {__Pyx_DECREF(r); }} while(0)
  #define __Pyx_XGOTREF(r)  do { if((r) == NULL); else {__Pyx_GOTREF(r); }} while(0)
  #define __Pyx_XGIVEREF(r) do { if((r) == NULL); else {__Pyx_GIVEREF(r);}} while(0)
#else
  #define __Pyx_RefNannyDeclarations
  #define __Pyx_RefNannySetupContext(name, acquire_gil)
  #define __Pyx_RefNannyFinishContextNogil()
  #define __Pyx_RefNannyFinishContext()
  #define __Pyx_INCREF(r) Py_INCREF(r)
  #define __Pyx_DECREF(r) Py_DECREF(r)
  #define __Pyx_GOTREF(r)
  #define __Pyx_GIVEREF(r)
  #define __Pyx_XINCREF(r) Py_XINCREF(r)
  #define __Pyx_XDECREF(r) Py_XDECREF(r)
  #define __Pyx_XGOTREF(r)
  #define __Pyx_XGIVEREF(r)
#endif /* CYTHON_REFNANNY */

#define __Pyx_Py_XDECREF_SET(r, v) do {                         \
        PyObject *tmp = (PyObject *) r;                         \
        r = v; Py_XDECREF(tmp);                                 \
    } while (0)
#define __Pyx_XDECREF_SET(r, v) do {                            \
        PyObject *tmp = (PyObject *) r;                         \
        r = v; __Pyx_XDECREF(tmp);                              \
    } while (0)
#define __Pyx_DECREF_SET(r, v) do {                             \
        PyObject *tmp = (PyObject *) r;                         \
        r = v; __Pyx_DECREF(tmp);                               \
    } while (0)

#define __Pyx_CLEAR(r)    do { PyObject* tmp = ((PyObject*)(r)); r = NULL; __Pyx_DECREF(tmp);} while(0)
#define __Pyx_XCLEAR(r)   do { if((r) != NULL) {PyObject* tmp = ((PyObject*)(r)); r = NULL; __Pyx_DECREF(tmp);}} while(0)

/////////////// Refnanny ///////////////

#if CYTHON_REFNANNY
static __Pyx_RefNannyAPIStruct *__Pyx_RefNannyImportAPI(const char *modname) {
    PyObject *m = NULL, *p = NULL;
    void *r = NULL;
    m = PyImport_ImportModule(modname);
    if (!m) goto end;
    p = PyObject_GetAttrString(m, "RefNannyAPI");
    if (!p) goto end;
    r = PyLong_AsVoidPtr(p);
end:
    Py_XDECREF(p);
    Py_XDECREF(m);
    return (__Pyx_RefNannyAPIStruct *)r;
}
#endif /* CYTHON_REFNANNY */


/////////////// ImportRefnannyAPI ///////////////

#if CYTHON_REFNANNY
__Pyx_RefNanny = __Pyx_RefNannyImportAPI("refnanny");
if (!__Pyx_RefNanny) {
  PyErr_Clear();
  __Pyx_RefNanny = __Pyx_RefNannyImportAPI("Cython.Runtime.refnanny");
  if (!__Pyx_RefNanny)
      Py_FatalError("failed to import 'refnanny' module");
}
#endif


/////////////// RegisterModuleCleanup.proto ///////////////
//@substitute: naming

static void ${cleanup_cname}(PyObject *self); /*proto*/

#if CYTHON_COMPILING_IN_PYPY
static int __Pyx_RegisterCleanup(void); /*proto*/
#else
#define __Pyx_RegisterCleanup() (0)
#endif

/////////////// RegisterModuleCleanup ///////////////
//@substitute: naming

#if CYTHON_COMPILING_IN_PYPY
static PyObject* ${cleanup_cname}_atexit(PyObject *module, PyObject *unused) {
    CYTHON_UNUSED_VAR(unused);
    ${cleanup_cname}(module);
    Py_INCREF(Py_None); return Py_None;
}

static int __Pyx_RegisterCleanup(void) {
    // Don't use Py_AtExit because that has a 32-call limit and is called
    // after python finalization.
    // Also, we try to prepend the cleanup function to "atexit._exithandlers"
    // in Py2 because CPython runs them last-to-first. Being run last allows
    // user exit code to run before us that may depend on the globals
    // and cached objects that we are about to clean up.

    static PyMethodDef cleanup_def = {
        "__cleanup", (PyCFunction)${cleanup_cname}_atexit, METH_NOARGS, 0};

    PyObject *cleanup_func = 0;
    PyObject *atexit = 0;
    PyObject *reg = 0;
    PyObject *args = 0;
    PyObject *res = 0;
    int ret = -1;

    cleanup_func = PyCFunction_New(&cleanup_def, 0);
    if (!cleanup_func)
        goto bad;

    atexit = PyImport_ImportModule("atexit");
    if (!atexit)
        goto bad;
    reg = PyObject_GetAttrString(atexit, "_exithandlers");
    if (reg && PyList_Check(reg)) {
        PyObject *a, *kw;
        a = PyTuple_New(0);
        kw = PyDict_New();
        if (!a || !kw) {
            Py_XDECREF(a);
            Py_XDECREF(kw);
            goto bad;
        }
        args = PyTuple_Pack(3, cleanup_func, a, kw);
        Py_DECREF(a);
        Py_DECREF(kw);
        if (!args)
            goto bad;
        ret = PyList_Insert(reg, 0, args);
    } else {
        if (!reg)
            PyErr_Clear();
        Py_XDECREF(reg);
        reg = PyObject_GetAttrString(atexit, "register");
        if (!reg)
            goto bad;
        args = PyTuple_Pack(1, cleanup_func);
        if (!args)
            goto bad;
        res = PyObject_CallObject(reg, args);
        if (!res)
            goto bad;
        ret = 0;
    }
bad:
    Py_XDECREF(cleanup_func);
    Py_XDECREF(atexit);
    Py_XDECREF(reg);
    Py_XDECREF(args);
    Py_XDECREF(res);
    return ret;
}
#endif

/////////////// FastGil.init ///////////////
__Pyx_FastGilFuncInit();

/////////////// NoFastGil.proto ///////////////
//@proto_block: utility_code_proto_before_types

#define __Pyx_PyGILState_Ensure PyGILState_Ensure
#define __Pyx_PyGILState_Release PyGILState_Release
#define __Pyx_FastGIL_Remember()
#define __Pyx_FastGIL_Forget()
#define __Pyx_FastGilFuncInit()

/////////////// FastGil.proto ///////////////
//@proto_block: utility_code_proto_before_types

#if CYTHON_FAST_GIL

struct __Pyx_FastGilVtab {
  PyGILState_STATE (*Fast_PyGILState_Ensure)(void);
  void (*Fast_PyGILState_Release)(PyGILState_STATE oldstate);
  void (*FastGIL_Remember)(void);
  void (*FastGIL_Forget)(void);
};

static void __Pyx_FastGIL_Noop(void) {}
static struct __Pyx_FastGilVtab __Pyx_FastGilFuncs = {
  PyGILState_Ensure,
  PyGILState_Release,
  __Pyx_FastGIL_Noop,
  __Pyx_FastGIL_Noop
};

static void __Pyx_FastGilFuncInit(void);

#define __Pyx_PyGILState_Ensure __Pyx_FastGilFuncs.Fast_PyGILState_Ensure
#define __Pyx_PyGILState_Release __Pyx_FastGilFuncs.Fast_PyGILState_Release
#define __Pyx_FastGIL_Remember __Pyx_FastGilFuncs.FastGIL_Remember
#define __Pyx_FastGIL_Forget __Pyx_FastGilFuncs.FastGIL_Forget

#ifndef CYTHON_THREAD_LOCAL
  #if defined(__cplusplus) && __cplusplus >= 201103L
    #define CYTHON_THREAD_LOCAL thread_local
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 201112
    #define CYTHON_THREAD_LOCAL _Thread_local
  #elif defined(__GNUC__)
    #define CYTHON_THREAD_LOCAL __thread
  #elif defined(_MSC_VER)
    #define CYTHON_THREAD_LOCAL __declspec(thread)
  #endif
#endif

#else
#define __Pyx_PyGILState_Ensure PyGILState_Ensure
#define __Pyx_PyGILState_Release PyGILState_Release
#define __Pyx_FastGIL_Remember()
#define __Pyx_FastGIL_Forget()
#define __Pyx_FastGilFuncInit()
#endif

/////////////// FastGil ///////////////
//@requires: AddModuleRef
// The implementations of PyGILState_Ensure/Release calls PyThread_get_key_value
// several times which is turns out to be quite slow (slower in fact than
// acquiring the GIL itself).  Simply storing it in a thread local for the
// common case is much faster.
// To make optimal use of this thread local, we attempt to share it between
// modules.

#if CYTHON_FAST_GIL

#define __Pyx_FastGIL_ABI_module __PYX_ABI_MODULE_NAME
#define __Pyx_FastGIL_PyCapsuleName "FastGilFuncs"
#define __Pyx_FastGIL_PyCapsule \
    __Pyx_FastGIL_ABI_module "." __Pyx_FastGIL_PyCapsuleName

#ifdef CYTHON_THREAD_LOCAL

#include "pythread.h"
#include "pystate.h"

static CYTHON_THREAD_LOCAL PyThreadState *__Pyx_FastGil_tcur = NULL;
static CYTHON_THREAD_LOCAL int __Pyx_FastGil_tcur_depth = 0;
static int __Pyx_FastGil_autoTLSkey = -1;

static CYTHON_INLINE void __Pyx_FastGIL_Remember0(void) {
  ++__Pyx_FastGil_tcur_depth;
}

static CYTHON_INLINE void __Pyx_FastGIL_Forget0(void) {
  if (--__Pyx_FastGil_tcur_depth == 0) {
    __Pyx_FastGil_tcur = NULL;
  }
}

static CYTHON_INLINE PyThreadState *__Pyx_FastGil_get_tcur(void) {
  PyThreadState *tcur = __Pyx_FastGil_tcur;
  if (tcur == NULL) {
    tcur = __Pyx_FastGil_tcur = (PyThreadState*)PyThread_get_key_value(__Pyx_FastGil_autoTLSkey);
  }
  return tcur;
}

static PyGILState_STATE __Pyx_FastGil_PyGILState_Ensure(void) {
  int current;
  PyThreadState *tcur;
  __Pyx_FastGIL_Remember0();
  tcur = __Pyx_FastGil_get_tcur();
  if (tcur == NULL) {
    // Uninitialized, need to initialize now.
    return PyGILState_Ensure();
  }
  current = tcur == __Pyx_PyThreadState_Current;
  if (current == 0) {
    PyEval_RestoreThread(tcur);
  }
  ++tcur->gilstate_counter;
  return current ? PyGILState_LOCKED : PyGILState_UNLOCKED;
}

static void __Pyx_FastGil_PyGILState_Release(PyGILState_STATE oldstate) {
  PyThreadState *tcur = __Pyx_FastGil_get_tcur();
  __Pyx_FastGIL_Forget0();
  if (tcur->gilstate_counter == 1) {
    // This is the last lock, do all the cleanup as well.
    PyGILState_Release(oldstate);
  } else {
    --tcur->gilstate_counter;
    if (oldstate == PyGILState_UNLOCKED) {
      PyEval_SaveThread();
    }
  }
}

static void __Pyx_FastGilFuncInit0(void) {
  /* Try to detect autoTLSkey. */
  int key;
  void* this_thread_state = (void*) PyGILState_GetThisThreadState();
  for (key = 0; key < 100; key++) {
    if (PyThread_get_key_value(key) == this_thread_state) {
      __Pyx_FastGil_autoTLSkey = key;
      break;
    }
  }
  if (__Pyx_FastGil_autoTLSkey != -1) {
    PyObject* capsule = NULL;
    PyObject* abi_module = NULL;
    __Pyx_PyGILState_Ensure = __Pyx_FastGil_PyGILState_Ensure;
    __Pyx_PyGILState_Release = __Pyx_FastGil_PyGILState_Release;
    __Pyx_FastGIL_Remember = __Pyx_FastGIL_Remember0;
    __Pyx_FastGIL_Forget = __Pyx_FastGIL_Forget0;
    capsule = PyCapsule_New(&__Pyx_FastGilFuncs, __Pyx_FastGIL_PyCapsule, NULL);
    if (capsule) {
        abi_module = __Pyx_PyImport_AddModuleRef(__Pyx_FastGIL_ABI_module);
        if (abi_module) {
          PyObject_SetAttrString(abi_module, __Pyx_FastGIL_PyCapsuleName, capsule);
          Py_DECREF(abi_module);
        }
    }
    Py_XDECREF(capsule);
  }
}

#else

static void __Pyx_FastGilFuncInit0(void) {
}

#endif

static void __Pyx_FastGilFuncInit(void) {
  struct __Pyx_FastGilVtab* shared = (struct __Pyx_FastGilVtab*)PyCapsule_Import(__Pyx_FastGIL_PyCapsule, 1);
  if (shared) {
    __Pyx_FastGilFuncs = *shared;
  } else {
   PyErr_Clear();
    __Pyx_FastGilFuncInit0();
  }
}

#endif

///////////////////// PretendToInitialize ////////////////////////

#ifdef __cplusplus
// In C++ a variable must actually be initialized to make returning
// it defined behaviour, and there doesn't seem to be a viable compiler trick to
// avoid that.
#if __cplusplus > 201103L
#include <type_traits>
#endif

template <typename T>
static void __Pyx_pretend_to_initialize(T* ptr) {
    // In C++11 we have enough introspection to work out which types it's actually
    // necessary to apply this to (non-trivial types will have been initialized by
    // the definition). Below C++11 just initialize everything.
#if __cplusplus > 201103L
    if ((std::is_trivially_default_constructible<T>::value))
#endif
        *ptr = T();
    (void)ptr;
}
#else
// For C,  taking an address of a variable is enough to make returning it
// defined behaviour.
static CYTHON_INLINE void __Pyx_pretend_to_initialize(void* ptr) { (void)ptr; }
#endif

///////////////////// UtilityCodePragmas /////////////////////////

#ifdef _MSC_VER
#pragma warning( push )
/* Warning 4127: conditional expression is constant
 * Cython uses constant conditional expressions to allow in inline functions to be optimized at
 * compile-time, so this warning is not useful
 */
#pragma warning( disable : 4127 )
#endif

///////////////////// UtilityCodePragmasEnd //////////////////////

#ifdef _MSC_VER
#pragma warning( pop )  /* undo whatever Cython has done to warnings */
#endif


//////////////////// NewCodeObj.proto ////////////////////////
//@proto_block: init_codeobjects

static PyObject* __Pyx_PyCode_New(
        //int argcount,
        //int num_posonly_args,
        //int num_kwonly_args,
        //int nlocals,
        // int s,
        //int flags,
        //int first_line,
        const __Pyx_PyCode_New_function_description descr,
        // PyObject *code,
        // PyObject *consts,
        // PyObject* n,
        // PyObject *varnames_tuple,
        PyObject * const *varnames,
        // PyObject *freevars,
        // PyObject *cellvars,
        PyObject *filename,
        PyObject *funcname,
        PyObject *line_table,
        PyObject *tuple_dedup_map
);/*proto*/

//////////////////// NewCodeObj ////////////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
    // Note that the limited API doesn't know about PyCodeObject, so the type of this
    // is PyObject (unlike for the main API)
    static PyObject* __Pyx__PyCode_New(int a, int p, int k, int l, int s, int f,
                                       PyObject *code, PyObject *c, PyObject* n, PyObject *v,
                                       PyObject *fv, PyObject *cell, PyObject* fn,
                                       PyObject *name, int fline, PyObject *lnos) {
        // Backup option for generating a code object.
        // PyCode_NewEmpty isn't in the limited API. Therefore the two options are
        //  1. Python call of the code type with a long list of positional args.
        //  2. Generate a code object by compiling some trivial code, and customize.
        // We use the first option here.
        PyObject *exception_table = NULL;
        PyObject *types_module=NULL, *code_type=NULL, *result=NULL;
        #if __PYX_LIMITED_VERSION_HEX < 0x030b0000
        PyObject *version_info;  /* borrowed */
        PyObject *py_minor_version = NULL;
        #endif
        long minor_version = 0;
        PyObject *type, *value, *traceback;

        // we must be able to call this while an exception is happening - thus clear then restore the state
        PyErr_Fetch(&type, &value, &traceback);

        #if __PYX_LIMITED_VERSION_HEX >= 0x030b0000
        minor_version = 11;
        // we don't yet need to distinguish between versions > 11
        // Note that from 3.13, when we do, we can use Py_Version
        #else
        if (!(version_info = PySys_GetObject("version_info"))) goto end;
        if (!(py_minor_version = PySequence_GetItem(version_info, 1))) goto end;
        minor_version = PyLong_AsLong(py_minor_version);
        Py_DECREF(py_minor_version);
        if (minor_version == -1 && PyErr_Occurred()) goto end;
        #endif

        if (!(types_module = PyImport_ImportModule("types"))) goto end;
        if (!(code_type = PyObject_GetAttrString(types_module, "CodeType"))) goto end;

        if (minor_version <= 7) {
            // 3.7:
            // code(argcount, kwonlyargcount, nlocals, stacksize, flags, codestring,
            //        constants, names, varnames, filename, name, firstlineno,
            //        lnotab[, freevars[, cellvars]])
            (void)p;
            result = PyObject_CallFunction(code_type, "iiiiiOOOOOOiOOO", a, k, l, s, f, code,
                          c, n, v, fn, name, fline, lnos, fv, cell);
        } else if (minor_version <= 10) {
            // 3.8, 3.9, 3.10
            // code(argcount, posonlyargcount, kwonlyargcount, nlocals, stacksize,
            //    flags, codestring, constants, names, varnames, filename, name,
            //    firstlineno, lnotab[, freevars[, cellvars]])
            // 3.10 switches lnotab for linetable, but is otherwise the same
            result = PyObject_CallFunction(code_type, "iiiiiiOOOOOOiOOO", a,p, k, l, s, f, code,
                          c, n, v, fn, name, fline, lnos, fv, cell);
        } else {
            // 3.11, 3.12
            // code(argcount, posonlyargcount, kwonlyargcount, nlocals, stacksize,
            //    flags, codestring, constants, names, varnames, filename, name,
            //    qualname, firstlineno, linetable, exceptiontable, freevars=(), cellvars=(), /)
            // We use name and qualname for simplicity
            if (!(exception_table = PyBytes_FromStringAndSize(NULL, 0))) goto end;
            result = PyObject_CallFunction(code_type, "iiiiiiOOOOOOOiOOOO", a,p, k, l, s, f, code,
                          c, n, v, fn, name, name, fline, lnos, exception_table, fv, cell);
        }

    end:
        Py_XDECREF(code_type);
        Py_XDECREF(exception_table);
        Py_XDECREF(types_module);
        if (type) {
            PyErr_Restore(type, value, traceback);
        }
        return result;
    }

#elif PY_VERSION_HEX >= 0x030B0000
  static PyCodeObject* __Pyx__PyCode_New(int a, int p, int k, int l, int s, int f,
                                         PyObject *code, PyObject *c, PyObject* n, PyObject *v,
                                         PyObject *fv, PyObject *cell, PyObject* fn,
                                         PyObject *name, int fline, PyObject *lnos) {
    // As earlier versions, but
    //  1. pass an empty bytes string as exception_table
    //  2. pass name as qualname (TODO this might implementing properly in future)
    PyCodeObject *result;
    result =
      #if PY_VERSION_HEX >= 0x030C0000
        PyUnstable_Code_NewWithPosOnlyArgs
      #else
        PyCode_NewWithPosOnlyArgs
      #endif
        (a, p, k, l, s, f, code, c, n, v, fv, cell, fn, name, name, fline, lnos, EMPTY(bytes));

    #if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX >= 0x030c00A1
    if (likely(result))
        result->_co_firsttraceable = 0;
    #endif
    return result;
  }
#elif !CYTHON_COMPILING_IN_PYPY
  #define __Pyx__PyCode_New(a, p, k, l, s, f, code, c, n, v, fv, cell, fn, name, fline, lnos) \
          PyCode_NewWithPosOnlyArgs(a, p, k, l, s, f, code, c, n, v, fv, cell, fn, name, fline, lnos)
#else
  #define __Pyx__PyCode_New(a, p, k, l, s, f, code, c, n, v, fv, cell, fn, name, fline, lnos) \
          PyCode_New(a, k, l, s, f, code, c, n, v, fv, cell, fn, name, fline, lnos)
#endif

// This is a specialised helper function for creating Cython's function code objects.
// It only receives the arguments that differ between the Cython functions of the module.
// This minimises the calling code in the module init function.
static PyObject* __Pyx_PyCode_New(
        //int argcount,
        //int num_posonly_args,
        //int num_kwonly_args,
        //int nlocals,
        // int s,
        //int flags,
        //int first_line,
        const __Pyx_PyCode_New_function_description descr,
        // PyObject *code,
        // PyObject *consts,
        // PyObject* n,
        // PyObject *varnames_tuple,
        PyObject * const *varnames,
        // PyObject *freevars,
        // PyObject *cellvars,
        PyObject *filename,
        PyObject *funcname,
        // line table replaced lnotab in Py3.11 (PEP-626)
        PyObject *line_table,
        PyObject *tuple_dedup_map
) {

    PyObject *code_obj = NULL, *varnames_tuple_dedup = NULL, *code_bytes = NULL;
    Py_ssize_t var_count = (Py_ssize_t) descr.nlocals;

    PyObject *varnames_tuple = PyTuple_New(var_count);
    if (unlikely(!varnames_tuple)) return NULL;
    for (Py_ssize_t i=0; i < var_count; i++) {
        Py_INCREF(varnames[i]);
        if (__Pyx_PyTuple_SET_ITEM(varnames_tuple, i, varnames[i]) != (0)) goto done;
    }

    #if CYTHON_COMPILING_IN_LIMITED_API
    varnames_tuple_dedup = PyDict_GetItem(tuple_dedup_map, varnames_tuple);
    if (!varnames_tuple_dedup) {
        if (unlikely(PyDict_SetItem(tuple_dedup_map, varnames_tuple, varnames_tuple) < 0)) goto done;
        varnames_tuple_dedup = varnames_tuple;
    }
    #else
    varnames_tuple_dedup = PyDict_SetDefault(tuple_dedup_map, varnames_tuple, varnames_tuple);
    if (unlikely(!varnames_tuple_dedup)) goto done;
    #endif

    #if CYTHON_AVOID_BORROWED_REFS
    Py_INCREF(varnames_tuple_dedup);
    #endif

    if (__PYX_LIMITED_VERSION_HEX >= (0x030b0000) && line_table != NULL && !CYTHON_COMPILING_IN_GRAAL) {
        // Allocate a "byte code" array (oversized) to match the addresses in the line table.
        // Length and alignment must be a multiple of sizeof(_Py_CODEUNIT), which is CPython specific but currently 2.
        // CPython makes a copy of the code array internally, so make sure it's somewhat short (but not too short).
        Py_ssize_t line_table_length = __Pyx_PyBytes_GET_SIZE(line_table);
        #if !CYTHON_ASSUME_SAFE_SIZE
        if (unlikely(line_table_length == -1)) goto done;
        #endif
        Py_ssize_t code_len = (line_table_length * 2 + 4) & ~3LL;
        code_bytes = PyBytes_FromStringAndSize(NULL, code_len);
        if (unlikely(!code_bytes)) goto done;
        char* c_code_bytes = PyBytes_AsString(code_bytes);
        if (unlikely(!c_code_bytes)) goto done;
        // We initialise the code array to '\0' even though a NOP would be more accurate,
        // but NOP changes its byte code ID across Python versions/implementations.
        memset(c_code_bytes, 0, (size_t) code_len);
    }

    code_obj = (PyObject*) __Pyx__PyCode_New(
        (int) descr.argcount,
        (int) descr.num_posonly_args,
        (int) descr.num_kwonly_args,
        (int) descr.nlocals,
        0,
        (int) descr.flags,
        code_bytes ? code_bytes : EMPTY(bytes),
        EMPTY(tuple),
        EMPTY(tuple),
        varnames_tuple_dedup,
        EMPTY(tuple),
        EMPTY(tuple),
        filename,
        funcname,
        (int) descr.first_line,
        (__PYX_LIMITED_VERSION_HEX >= (0x030b0000) && line_table) ? line_table : EMPTY(bytes)
    );

done:
    Py_XDECREF(code_bytes);
    #if CYTHON_AVOID_BORROWED_REFS
    Py_XDECREF(varnames_tuple_dedup);
    #endif
    Py_DECREF(varnames_tuple);
    return code_obj;
}


////////////////////////// MultiPhaseInitModuleState.proto /////////////

#if CYTHON_PEP489_MULTI_PHASE_INIT && CYTHON_USE_MODULE_STATE
// This defines an ad-hoc, single module version of PyState_FindModule that
// works for multi-phase init modules. It's intended to be the last option
// when all the other official ways of getting the module are unavailable.
static PyObject *__Pyx_State_FindModule(void*); /* proto */
static int __Pyx_State_AddModule(PyObject* module, void*); /* proto */
static int __Pyx_State_RemoveModule(void*); /* proto */

#elif CYTHON_USE_MODULE_STATE
#define __Pyx_State_FindModule PyState_FindModule
#define __Pyx_State_AddModule PyState_AddModule
#define __Pyx_State_RemoveModule PyState_RemoveModule
#endif

////////////////////////// MultiPhaseInitModuleState /////////////
//@requires: Synchronization.c::Atomics


// Code to maintain a mapping between (sub)interpreters and the module instance that they imported.
// This is used to find the correct module state for the current interpreter.

#if CYTHON_PEP489_MULTI_PHASE_INIT && CYTHON_USE_MODULE_STATE

#ifndef CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
// If you're using multiple interpreters but a single GIL then
// this can be undefined for a bit of speed.
// Isolated subinterpreters were added in 3.12, and nogil in 3.13, so before that
// we can safely assume that we're protected by the GIL.
// TODO - turn this off as appropriate when the user is able to set
// Py_MOD_PER_INTERPRETER_GIL_SUPPORTED explicitly.
#if (CYTHON_COMPILING_IN_LIMITED_API || PY_VERSION_HEX >= 0x030C0000)
  #define CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE 1
#else
  #define CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE 0
#endif
#endif

#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE && !CYTHON_ATOMICS
#error "Module state with PEP489 requires atomics. Currently that's one of\
 C11, C++11, gcc atomic intrinsics or MSVC atomic intrinsics"
#endif

// We also need a lock. In order of preference:
// - PyMutex
// - a language standard library
// - pthreads
// - msvc
// - PyThread_lock isn't acceptable since we can't initialize it in a thread safe way
#if !CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE

#define __Pyx_ModuleStateLookup_Lock()
#define __Pyx_ModuleStateLookup_Unlock()

#elif !CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x030d0000

static PyMutex __Pyx_ModuleStateLookup_mutex = {0};
#define __Pyx_ModuleStateLookup_Lock() PyMutex_Lock(&__Pyx_ModuleStateLookup_mutex)
#define __Pyx_ModuleStateLookup_Unlock() PyMutex_Unlock(&__Pyx_ModuleStateLookup_mutex)

#elif defined(__cplusplus) && __cplusplus >= 201103L

#include <mutex>
static std::mutex __Pyx_ModuleStateLookup_mutex;
#define __Pyx_ModuleStateLookup_Lock() __Pyx_ModuleStateLookup_mutex.lock()
#define __Pyx_ModuleStateLookup_Unlock() __Pyx_ModuleStateLookup_mutex.unlock()

#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ > 201112L) && !defined(__STDC_NO_THREADS__)
#include <threads.h>
static mtx_t __Pyx_ModuleStateLookup_mutex;
static once_flag __Pyx_ModuleStateLookup_mutex_once_flag = ONCE_FLAG_INIT;
static void __Pyx_ModuleStateLookup_initialize_mutex(void) {
    mtx_init(&__Pyx_ModuleStateLookup_mutex, mtx_plain);
}
#define __Pyx_ModuleStateLookup_Lock() \
  call_once(&__Pyx_ModuleStateLookup_mutex_once_flag, __Pyx_ModuleStateLookup_initialize_mutex); \
  mtx_lock(&__Pyx_ModuleStateLookup_mutex)
#define __Pyx_ModuleStateLookup_Unlock() mtx_unlock(&__Pyx_ModuleStateLookup_mutex)

// HAVE_PTHREAD_H comes from pyconfig.h
#elif defined(HAVE_PTHREAD_H)

#include <pthread.h>
static pthread_mutex_t __Pyx_ModuleStateLookup_mutex = PTHREAD_MUTEX_INITIALIZER;
#define __Pyx_ModuleStateLookup_Lock() pthread_mutex_lock(&__Pyx_ModuleStateLookup_mutex)
#define __Pyx_ModuleStateLookup_Unlock() pthread_mutex_unlock(&__Pyx_ModuleStateLookup_mutex)

#elif defined(_WIN32)

#include <Windows.h>  // synchapi.h on its own doesn't work

// Using a slim-read-write lock (instead of a mutex/critical section)
// because it can be statically initialized.
static SRWLOCK __Pyx_ModuleStateLookup_mutex = SRWLOCK_INIT;
#define __Pyx_ModuleStateLookup_Lock() AcquireSRWLockExclusive(&__Pyx_ModuleStateLookup_mutex)
#define __Pyx_ModuleStateLookup_Unlock() ReleaseSRWLockExclusive(&__Pyx_ModuleStateLookup_mutex)

#else
#error "No suitable lock available for CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE.\
 Requires C standard >= C11, or C++ standard >= C++11,\
 or pthreads, or the Windows 32 API, or Python >= 3.13."
#endif


typedef struct {
    int64_t id;
    PyObject *module;
} __Pyx_InterpreterIdAndModule;

typedef struct {
    char interpreter_id_as_index;
    Py_ssize_t count;
    Py_ssize_t allocated;
    __Pyx_InterpreterIdAndModule table[1];
} __Pyx_ModuleStateLookupData;

#define __PYX_MODULE_STATE_LOOKUP_SMALL_SIZE 32
// "interpreter_id_as_index" above means "the maximum interpreter ID ever seen is smaller than
// __PYX_MODULE_STATE_LOOKUP_SMALL_SIZE and thus they're stored in an array
// where the index corresponds to interpreter ID, and __Pyx_ModuleStateLookup_count
// is the size of the array.

#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
static __pyx_atomic_int_type __Pyx_ModuleStateLookup_read_counter = 0;
#endif

// A sorted list of (sub)interpreter IDs and the module that was imported into that interpreter.
// For now look this up via binary search.
#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
static __pyx_atomic_ptr_type __Pyx_ModuleStateLookup_data = 0;
#else
static __Pyx_ModuleStateLookupData* __Pyx_ModuleStateLookup_data = NULL;
#endif


static __Pyx_InterpreterIdAndModule* __Pyx_State_FindModuleStateLookupTableLowerBound(
        __Pyx_InterpreterIdAndModule* table,
        Py_ssize_t count,
        int64_t interpreterId) {
    __Pyx_InterpreterIdAndModule* begin = table;
    __Pyx_InterpreterIdAndModule* end = begin + count;

    // fairly likely - e.g. single interpreter
    if (begin->id == interpreterId) {
        return begin;
    }

    while ((end - begin) > __PYX_MODULE_STATE_LOOKUP_SMALL_SIZE) {
        __Pyx_InterpreterIdAndModule* halfway = begin + (end - begin)/2;
        if (halfway->id == interpreterId) {
            return halfway;
        }
        if (halfway->id < interpreterId) {
            begin = halfway;
        } else {
            end = halfway;
        }
    }

    // Assume that for small ranges, it's quicker to do a linear search
    for (; begin < end; ++begin) {
        if (begin->id >= interpreterId) return begin;
    }
    return begin;
}

static PyObject *__Pyx_State_FindModule(CYTHON_UNUSED void* dummy) {
    int64_t interpreter_id = PyInterpreterState_GetID(__Pyx_PyInterpreterState_Get());
    if (interpreter_id == -1) return NULL;

#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
    __Pyx_ModuleStateLookupData* data = (__Pyx_ModuleStateLookupData*)__pyx_atomic_pointer_load_relaxed(&__Pyx_ModuleStateLookup_data);
    {
        // Thread sanitizer says that this is OK relaxed, but I think it needs to be acquire-release
        __pyx_atomic_incr_acq_rel(&__Pyx_ModuleStateLookup_read_counter);
        // data == NULL can either mean we're writing, or it's uninitialized.
        // Uninitialized only happens infrequently on the first few calls, so it's fine
        // to be on the slow path.
        if (likely(data)) {
            __Pyx_ModuleStateLookupData* new_data = (__Pyx_ModuleStateLookupData*)__pyx_atomic_pointer_load_acquire(&__Pyx_ModuleStateLookup_data);
            if (likely(data == new_data)) {
                // Nothing has written the data between incrementing the read counter and loading the pointer.
                goto read_finished;
            }
        }
        // In principle DW believes this could be "relaxed", but it's on the unlikely slow path anyway
        // so let's not add more macros.
        // Undo our addition to the read counter.
        __pyx_atomic_decr_acq_rel(&__Pyx_ModuleStateLookup_read_counter);
        // Wait for the write to finish and try again
        __Pyx_ModuleStateLookup_Lock();
        __pyx_atomic_incr_relaxed(&__Pyx_ModuleStateLookup_read_counter);
        data = (__Pyx_ModuleStateLookupData*)__pyx_atomic_pointer_load_relaxed(&__Pyx_ModuleStateLookup_data);
        __Pyx_ModuleStateLookup_Unlock();
    }
  read_finished:;

#else
    __Pyx_ModuleStateLookupData* data = __Pyx_ModuleStateLookup_data;
#endif

    __Pyx_InterpreterIdAndModule* found = NULL;

    // There's one "already imported" check that'll hit this
    if (unlikely(!data)) goto end;

    if (data->interpreter_id_as_index) {
        if (interpreter_id < data->count) {
            found = data->table+interpreter_id;
        }
    } else {
        found = __Pyx_State_FindModuleStateLookupTableLowerBound(
            data->table, data->count, interpreter_id);
    }

  end:
    {
        PyObject *result=NULL;

        if (found && found->id == interpreter_id) {
            result = found->module;
        }
#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
        __pyx_atomic_decr_acq_rel(&__Pyx_ModuleStateLookup_read_counter);
#endif
        return result;
    }
}

#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
static void __Pyx_ModuleStateLookup_wait_until_no_readers(void) {
    // Wait for any readers still working on the old data. Spin-lock is
    // fine because readers should be much faster than memory allocation.
    while (__pyx_atomic_load(&__Pyx_ModuleStateLookup_read_counter) != 0);
}
#else
#define __Pyx_ModuleStateLookup_wait_until_no_readers()
#endif

static int __Pyx_State_AddModuleInterpIdAsIndex(__Pyx_ModuleStateLookupData **old_data, PyObject* module, int64_t interpreter_id) {
    Py_ssize_t to_allocate = (*old_data)->allocated;
    while (to_allocate <= interpreter_id) {
        if (to_allocate == 0) to_allocate = 1;
        else to_allocate *= 2;
    }
    __Pyx_ModuleStateLookupData *new_data = *old_data;
    if (to_allocate != (*old_data)->allocated) {
         new_data = (__Pyx_ModuleStateLookupData *)realloc(
            *old_data,
            sizeof(__Pyx_ModuleStateLookupData)+(to_allocate-1)*sizeof(__Pyx_InterpreterIdAndModule));
        if (!new_data) {
            PyErr_NoMemory();
            return -1;
        }
        for (Py_ssize_t i = new_data->allocated; i < to_allocate; ++i) {
            new_data->table[i].id = i;
            new_data->table[i].module = NULL;
        }
        new_data->allocated = to_allocate;
    }
    new_data->table[interpreter_id].module = module;
    if (new_data->count < interpreter_id+1) {
        new_data->count = interpreter_id+1;
    }
    *old_data = new_data;
    return 0;
}

static void __Pyx_State_ConvertFromInterpIdAsIndex(__Pyx_ModuleStateLookupData *data) {
    __Pyx_InterpreterIdAndModule *read = data->table;
    __Pyx_InterpreterIdAndModule *write = data->table;
    __Pyx_InterpreterIdAndModule *end = read + data->count;

    for (; read<end; ++read) {
        if (read->module) {
            write->id = read->id;
            write->module = read->module;
            ++write;
        }
        // Otherwise empty; don't copy
    }
    data->count = write - data->table;
    for (; write<end; ++write) {
        // clear rest of array
        write->id = 0;
        write->module = NULL;
    }
    data->interpreter_id_as_index = 0;
}

static int __Pyx_State_AddModule(PyObject* module, CYTHON_UNUSED void* dummy) {
    int64_t interpreter_id = PyInterpreterState_GetID(__Pyx_PyInterpreterState_Get());
    if (interpreter_id == -1) return -1;

    int result = 0;

    __Pyx_ModuleStateLookup_Lock();

    // Adding modules is the slow path so I've not thought about memory ordering much and
    // just made it strict.
#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
    // we're working and maybe modifying it, swap for 0
    __Pyx_ModuleStateLookupData *old_data = (__Pyx_ModuleStateLookupData *)
            __pyx_atomic_pointer_exchange(&__Pyx_ModuleStateLookup_data, 0);
#else
    __Pyx_ModuleStateLookupData *old_data = __Pyx_ModuleStateLookup_data;
#endif
    __Pyx_ModuleStateLookupData *new_data = old_data;

    if (!new_data) {
        // If we don't yet have anything, initialize
        new_data = (__Pyx_ModuleStateLookupData *)calloc(1, sizeof(__Pyx_ModuleStateLookupData));
        if (!new_data) {
            result = -1;
            PyErr_NoMemory();
            goto end;
        }
        new_data->allocated = 1;
        new_data->interpreter_id_as_index = 1;
    }

    // Pretty much everything from here modifies the data, and so requires us to wait
    // until all existing readers have finished in order to be thread-safe.
    __Pyx_ModuleStateLookup_wait_until_no_readers();
    if (new_data->interpreter_id_as_index) {
        if (interpreter_id < __PYX_MODULE_STATE_LOOKUP_SMALL_SIZE) {
            result = __Pyx_State_AddModuleInterpIdAsIndex(&new_data, module, interpreter_id);
            goto end;
        }
        // otherwise we have to convert then proceed with a normal insertion
        __Pyx_State_ConvertFromInterpIdAsIndex(new_data);
    }
    {
        Py_ssize_t insert_at = 0;
        {
            __Pyx_InterpreterIdAndModule* lower_bound = __Pyx_State_FindModuleStateLookupTableLowerBound(
                new_data->table, new_data->count, interpreter_id);

            assert(lower_bound);

            insert_at = lower_bound - new_data->table;

            if (unlikely(insert_at < new_data->count && lower_bound->id == interpreter_id)) {
                lower_bound->module = module;
                goto end;  // already in table, nothing more to do
            }

        }

        if (new_data->count+1 >= new_data->allocated) {
            // Use C realloc. PyMem_RawMalloc is added to the limited API fairly late (3.13)
            // and we want allocation independent of the interpreter which I think excludes PyMem_Malloc.
            Py_ssize_t to_allocate = (new_data->count+1)*2;
            new_data =
                (__Pyx_ModuleStateLookupData*)realloc(
                    new_data,
                    sizeof(__Pyx_ModuleStateLookupData) +
                    (to_allocate-1)*sizeof(__Pyx_InterpreterIdAndModule));
            if (!new_data) {
                result = -1;
                new_data = old_data;
                PyErr_NoMemory();
                goto end;
            }
            new_data->allocated = to_allocate;
        }

        ++new_data->count;

        int64_t last_id = interpreter_id;
        PyObject *last_module = module;
        for (Py_ssize_t i=insert_at; i<new_data->count; ++i) {
            int64_t current_id = new_data->table[i].id;
            new_data->table[i].id = last_id;
            last_id = current_id;
            PyObject *current_module = new_data->table[i].module;
            new_data->table[i].module = last_module;
            last_module = current_module;
        }
    }

  end:
#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
    __pyx_atomic_pointer_exchange(&__Pyx_ModuleStateLookup_data, new_data);
#else
    __Pyx_ModuleStateLookup_data = new_data;
#endif

    __Pyx_ModuleStateLookup_Unlock();
    return result;
}

static int __Pyx_State_RemoveModule(CYTHON_UNUSED void* dummy) {
    int64_t interpreter_id = PyInterpreterState_GetID(__Pyx_PyInterpreterState_Get());
    if (interpreter_id == -1) return -1;

    __Pyx_ModuleStateLookup_Lock();

#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
    __Pyx_ModuleStateLookupData *data = (__Pyx_ModuleStateLookupData *)
            __pyx_atomic_pointer_exchange(&__Pyx_ModuleStateLookup_data, 0);
#else
    __Pyx_ModuleStateLookupData *data = __Pyx_ModuleStateLookup_data;
#endif

    if (data->interpreter_id_as_index) {
        if (interpreter_id < data->count) {
            data->table[interpreter_id].module = NULL;
        }
        goto done;
    }
    {
        __Pyx_ModuleStateLookup_wait_until_no_readers();

        __Pyx_InterpreterIdAndModule* lower_bound = __Pyx_State_FindModuleStateLookupTableLowerBound(
            data->table, data->count, interpreter_id);

        // TODO Errors here?
        if (!lower_bound) goto done;
        if (lower_bound->id != interpreter_id) goto done;

        __Pyx_InterpreterIdAndModule *end = data->table+data->count;
        for (;lower_bound<end-1; ++lower_bound) {
            lower_bound->id = (lower_bound+1)->id;
            lower_bound->module = (lower_bound+1)->module;
        }
    }
    --data->count;
    if (data->count == 0) {
        free(data);
        data = NULL;
    }
    // For now, never shrink the allocated table.
  done:
#if CYTHON_MODULE_STATE_LOOKUP_THREAD_SAFE
    __pyx_atomic_pointer_exchange(&__Pyx_ModuleStateLookup_data, data);
#else
    __Pyx_ModuleStateLookup_data = data;
#endif
    __Pyx_ModuleStateLookup_Unlock();
    return 0;
}

#endif

////////////////////// IncludeStdlibH.proto //////////////////////

#include <stdlib.h>

////////////////////// ReleaseUnknownGil.proto ///////////////////

#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
typedef struct {
  PyThreadState* ts;
  PyGILState_STATE gil_state;
} __Pyx_UnknownThreadState;
#else
#define __Pyx_UnknownThreadState PyThreadState*
#endif

static __Pyx_UnknownThreadState __Pyx_SaveUnknownThread(void); /* proto */
static void __Pyx_RestoreUnknownThread(__Pyx_UnknownThreadState state); /* proto */

// Note - may return false negatives for some versions of the Limited API
static CYTHON_INLINE int __Pyx_UnknownThreadStateDefinitelyHadGil(__Pyx_UnknownThreadState state); /* proto */
// Note - may return false positives for some versions of the Limited API
static CYTHON_INLINE int __Pyx_UnknownThreadStateMayHaveHadGil(__Pyx_UnknownThreadState state); /* proto */

////////////////////// ReleaseUnknownGil ///////////////////

static __Pyx_UnknownThreadState __Pyx_SaveUnknownThread(void) {
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
    if (__Pyx_get_runtime_version() >= 0x030d0000) {
        // On Python 3.13+ we can just swap the thread state for NULL, and that's fine,
        // whether we have the GIL or not.
        PyThreadState *ts = PyThreadState_Swap(NULL);
        __Pyx_UnknownThreadState out = { ts, PyGILState_UNLOCKED /* arbitrary value */ };
        return out;
    } else {
        // On Python <3.13 there's no way of telling if we have the GIL in the Limited API.
        // Therefore, the only safe thing to do is to acquire it then release it.
        PyGILState_STATE gil_state = PyGILState_Ensure();
        PyThreadState *ts = PyEval_SaveThread();
        __Pyx_UnknownThreadState out = { ts, gil_state };
        return out;
    }
#elif __PYX_LIMITED_VERSION_HEX >= 0x030d0000
    return PyThreadState_Swap(NULL);
#else
  #if CYTHON_COMPILING_IN_PYPY || PY_VERSION_HEX < 0x030C0000
    if (PyGILState_Check())
  #else
    if (_PyThreadState_UncheckedGet())  // UncheckedGet is a reliable check for the GIL from 3.12 upwards
  #endif
    {
      return PyEval_SaveThread();
    }
    return NULL; // Nothing to release - we don't have the GIL
#endif
}

static void __Pyx_RestoreUnknownThread(__Pyx_UnknownThreadState state) {
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
    if (!state.ts) return;
    PyEval_RestoreThread(state.ts);
    if (__Pyx_get_runtime_version() < 0x030d0000) {
        PyGILState_Release(state.gil_state);
    }
#else
    if (state) {
        PyEval_RestoreThread(state);
    }
#endif
}

static CYTHON_INLINE int __Pyx_UnknownThreadStateDefinitelyHadGil(__Pyx_UnknownThreadState state) {
  #if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
  // For older Limited API versions we can't tell if we had the GIL
  return ((state.ts != NULL) && (__Pyx_get_runtime_version() >= 0x030d0000));
  #else
  return state != NULL;
  #endif
}

static CYTHON_INLINE int __Pyx_UnknownThreadStateMayHaveHadGil(__Pyx_UnknownThreadState state) {
  #if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
  return state.ts != NULL;
  #else
  return state != NULL;
  #endif
}

///////////////////// AddModuleRef.proto ///////////////////////

#if ((CYTHON_COMPILING_IN_CPYTHON_FREETHREADING /* && __PYX_LIMITED_VERSION_HEX < some future value */) || \
     __PYX_LIMITED_VERSION_HEX < 0x030d0000) 
  // https://github.com/python/cpython/issues/137422 - PyImport_AddModule(Ref) isn't thread safe!
  static PyObject *__Pyx_PyImport_AddModuleRef(const char *name); /* proto */
#else
  #define __Pyx_PyImport_AddModuleRef(name) PyImport_AddModuleRef(name)
#endif

///////////////////// AddModuleRef ////////////////////////////

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING /* && __PYX_LIMITED_VERSION_HEX < some future value */
  // https://github.com/python/cpython/issues/137422 - PyImport_AddModule(Ref) isn't thread safe!
  static PyObject *__Pyx_PyImport_AddModuleObjectRef(PyObject *name) {
      // We're going to assume nobody is swapping the module dict out from under us
      // (even though they're allowed to) because we really can't write code that's
      // safe against that.
      PyObject *module_dict = PyImport_GetModuleDict();
      PyObject *m;
      if (PyMapping_GetOptionalItem(module_dict, name, &m) < 0) { 
          return NULL; 
      } 
      if (m != NULL && PyModule_Check(m)) { 
          return m; 
      }
      Py_XDECREF(m); 
      m = PyModule_NewObject(name);
      if (m == NULL) 
          return NULL; 
      if (PyDict_CheckExact(module_dict)) {
          PyObject *new_m;
          (void)PyDict_SetDefaultRef(module_dict, name, m, &new_m);
          Py_DECREF(m);
          return new_m;
      } else {
            // For non-dict sys-modules I don't think it's possible to reliably make thread-safe.
           if (PyObject_SetItem(module_dict, name, m) != 0) { 
                Py_DECREF(m); 
                return NULL; 
            }         
            return m; 
      }
  }
  static PyObject *__Pyx_PyImport_AddModuleRef(const char *name) {
      PyObject *py_name = PyUnicode_FromString(name);
      if (!py_name) return NULL;
      PyObject *module = __Pyx_PyImport_AddModuleObjectRef(py_name); 
      Py_DECREF(py_name);
      return module;
  }
#elif __PYX_LIMITED_VERSION_HEX >= 0x030d0000
  #define __Pyx_PyImport_AddModuleRef(name) PyImport_AddModuleRef(name)
#else
  static PyObject *__Pyx_PyImport_AddModuleRef(const char *name) {
      PyObject *module = PyImport_AddModule(name);
      Py_XINCREF(module);
      return module;
  }
#endif
