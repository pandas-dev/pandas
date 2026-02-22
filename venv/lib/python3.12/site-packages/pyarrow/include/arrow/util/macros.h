// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#pragma once

#include <cstdint>

#define ARROW_EXPAND(x) x
#define ARROW_STRINGIFY(x) #x
#define ARROW_CONCAT(x, y) x##y

// From Google gutil
#ifndef ARROW_DISALLOW_COPY_AND_ASSIGN
#  define ARROW_DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&) = delete;            \
    void operator=(const TypeName&) = delete
#endif

#ifndef ARROW_DEFAULT_MOVE_AND_ASSIGN
#  define ARROW_DEFAULT_MOVE_AND_ASSIGN(TypeName) \
    TypeName(TypeName&&) = default;               \
    TypeName& operator=(TypeName&&) = default
#endif

// With ARROW_PREDICT_FALSE, GCC and clang can be told that a certain branch is
// not likely to be taken (for instance, a CHECK failure), and use that information in
// static analysis. Giving the compiler this information can affect the generated code
// layout in the absence of better information (i.e. -fprofile-arcs). [1] explains how
// this feature can be used to improve code generation. It was written as a positive
// comment to a negative article about the use of these annotations.
//
// ARROW_COMPILER_ASSUME allows the compiler to assume that a given expression is
// true, without evaluating it, and to optimise based on this assumption [2]. If this
// condition is violated at runtime, the behavior is undefined. This can be useful to
// generate both faster and smaller code in compute kernels.
//
// IMPORTANT: Different optimisers are likely to react differently to this annotation!
// It should be used with care when we can prove by some means that the assumption
// is (1) guaranteed to always hold and (2) is useful for optimization [3]. If the
// assumption is pessimistic, it might even block the compiler from decisions that
// could lead to better code [4]. If you have a good intuition for what the compiler
// can do with assumptions [5], you can use this macro to guide it and end up with
// results you would only get with more complex code transformations.
// `clang -S -emit-llvm` can be used to check how the generated code changes with
// your specific use of this macro.
//
// [1] https://lobste.rs/s/uwgtkt/don_t_use_likely_unlikely_attributes#c_xi3wmc
// [2] "Portable assumptions"
//     https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p1774r4.pdf
// [3] "Assertions Are Pessimistic, Assumptions Are Optimistic"
//     https://blog.regehr.org/archives/1096
// [4] https://discourse.llvm.org/t/llvm-assume-blocks-optimization/71609
// [5] J. Doerfert et al. 2019. "Performance Exploration Through Optimistic Static
//     Program Annotations". https://github.com/jdoerfert/PETOSPA/blob/master/ISC19.pdf
#define ARROW_UNUSED(x) (void)(x)
#ifdef ARROW_WARN_DOCUMENTATION
#  define ARROW_ARG_UNUSED(x) x
#else
#  define ARROW_ARG_UNUSED(x)
#endif
#if defined(__GNUC__)  // GCC and compatible compilers (clang, Intel ICC)
#  define ARROW_NORETURN __attribute__((noreturn))
#  define ARROW_NOINLINE __attribute__((noinline))
#  define ARROW_FORCE_INLINE __attribute__((always_inline))
#  define ARROW_PREDICT_FALSE(x) (__builtin_expect(!!(x), 0))
#  define ARROW_PREDICT_TRUE(x) (__builtin_expect(!!(x), 1))
#  define ARROW_RESTRICT __restrict
#  if defined(__clang__)  // clang-specific
#    define ARROW_COMPILER_ASSUME(expr) __builtin_assume(expr)
#  else  // GCC-specific
#    if __GNUC__ >= 13
#      define ARROW_COMPILER_ASSUME(expr) __attribute__((assume(expr)))
#    else
// GCC does not have a built-in assume intrinsic before GCC 13, so we use an
// if statement and __builtin_unreachable() to achieve the same effect [2].
// Unlike clang's __builtin_assume and C++23's [[assume(expr)]], using this
// on GCC won't warn about side-effects in the expression, so make sure expr
// is side-effect free when working with GCC versions before 13 (Jan-2024),
// otherwise clang/MSVC builds will fail in CI.
#      define ARROW_COMPILER_ASSUME(expr) \
        if (expr) {                       \
        } else {                          \
          __builtin_unreachable();        \
        }
#    endif  // __GNUC__ >= 13
#  endif
#elif defined(_MSC_VER)  // MSVC
#  define ARROW_NORETURN __declspec(noreturn)
#  define ARROW_NOINLINE __declspec(noinline)
#  define ARROW_FORCE_INLINE __forceinline
#  define ARROW_PREDICT_FALSE(x) (x)
#  define ARROW_PREDICT_TRUE(x) (x)
#  define ARROW_RESTRICT __restrict
#  define ARROW_COMPILER_ASSUME(expr) __assume(expr)
#else
#  define ARROW_NORETURN
#  define ARROW_NOINLINE
#  define ARROW_FORCE_INLINE
#  define ARROW_PREDICT_FALSE(x) (x)
#  define ARROW_PREDICT_TRUE(x) (x)
#  define ARROW_RESTRICT
#  define ARROW_COMPILER_ASSUME(expr)
#endif

// ----------------------------------------------------------------------
// C++/CLI support macros (see ARROW-1134)

#ifndef NULLPTR

#  ifdef __cplusplus_cli
#    define NULLPTR __nullptr
#  else
#    define NULLPTR nullptr
#  endif

#endif  // ifndef NULLPTR

// ----------------------------------------------------------------------

// clang-format off
// [[deprecated]] is only available in C++14, use this for the time being
// This macro takes an optional deprecation message
#ifdef __COVERITY__
#  define ARROW_DEPRECATED(...)
#else
#  define ARROW_DEPRECATED(...) [[deprecated(__VA_ARGS__)]]
#endif

#ifdef __COVERITY__
#  define ARROW_DEPRECATED_ENUM_VALUE(...)
#else
#  define ARROW_DEPRECATED_ENUM_VALUE(...) [[deprecated(__VA_ARGS__)]]
#endif

// clang-format on

// Macros to disable deprecation warnings

#ifdef __clang__
#  define ARROW_SUPPRESS_DEPRECATION_WARNING \
    _Pragma("clang diagnostic push");        \
    _Pragma("clang diagnostic ignored \"-Wdeprecated-declarations\"")
#  define ARROW_UNSUPPRESS_DEPRECATION_WARNING _Pragma("clang diagnostic pop")
#elif defined(__GNUC__)
#  define ARROW_SUPPRESS_DEPRECATION_WARNING \
    _Pragma("GCC diagnostic push");          \
    _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#  define ARROW_UNSUPPRESS_DEPRECATION_WARNING _Pragma("GCC diagnostic pop")
#elif defined(_MSC_VER)
#  define ARROW_SUPPRESS_DEPRECATION_WARNING \
    __pragma(warning(push)) __pragma(warning(disable : 4996))
#  define ARROW_UNSUPPRESS_DEPRECATION_WARNING __pragma(warning(pop))
#else
#  define ARROW_SUPPRESS_DEPRECATION_WARNING
#  define ARROW_UNSUPPRESS_DEPRECATION_WARNING
#endif

// ----------------------------------------------------------------------

// Macros to disable warnings about undeclared global functions
#if defined(__GNUC__)
#  define ARROW_SUPPRESS_MISSING_DECLARATIONS_WARNING \
    _Pragma("GCC diagnostic push");                   \
    _Pragma("GCC diagnostic ignored \"-Wmissing-declarations\"")
#  define ARROW_UNSUPPRESS_MISSING_DECLARATIONS_WARNING _Pragma("GCC diagnostic pop")
#else
#  define ARROW_SUPPRESS_MISSING_DECLARATIONS_WARNING
#  define ARROW_UNSUPPRESS_MISSING_DECLARATIONS_WARNING
#endif

// ----------------------------------------------------------------------

// macros to disable padding
// these macros are portable across different compilers and platforms
//[https://github.com/google/flatbuffers/blob/master/include/flatbuffers/flatbuffers.h#L1355]
#if !defined(MANUALLY_ALIGNED_STRUCT)
#  if defined(_MSC_VER)
#    define MANUALLY_ALIGNED_STRUCT(alignment) \
      __pragma(pack(1));                       \
      struct __declspec(align(alignment))
#    define STRUCT_END(name, size) \
      __pragma(pack());            \
      static_assert(sizeof(name) == size, "compiler breaks packing rules")
#  elif defined(__GNUC__) || defined(__clang__)
#    define MANUALLY_ALIGNED_STRUCT(alignment) \
      _Pragma("pack(1)") struct __attribute__((aligned(alignment)))
#    define STRUCT_END(name, size)                          \
      _Pragma("pack()") static_assert(sizeof(name) == size, \
                                      "compiler breaks packing rules")
#  else
#    error Unknown compiler, please define structure alignment macros
#  endif
#endif  // !defined(MANUALLY_ALIGNED_STRUCT)

// ----------------------------------------------------------------------
// Convenience macro disabling a particular UBSan check in a function

#if defined(__clang__)
#  define ARROW_DISABLE_UBSAN(feature) __attribute__((no_sanitize(feature)))
#else
#  define ARROW_DISABLE_UBSAN(feature)
#endif

// ----------------------------------------------------------------------
// Machine information

#if INTPTR_MAX == INT64_MAX
#  define ARROW_BITNESS 64
#elif INTPTR_MAX == INT32_MAX
#  define ARROW_BITNESS 32
#else
#  error Unexpected INTPTR_MAX
#endif

// ----------------------------------------------------------------------
// From googletest
// (also in parquet-cpp)

// When you need to test the private or protected members of a class,
// use the FRIEND_TEST macro to declare your tests as friends of the
// class.  For example:
//
// class MyClass {
//  private:
//   void MyMethod();
//   FRIEND_TEST(MyClassTest, MyMethod);
// };
//
// class MyClassTest : public testing::Test {
//   // ...
// };
//
// TEST_F(MyClassTest, MyMethod) {
//   // Can call MyClass::MyMethod() here.
// }

#define FRIEND_TEST(test_case_name, test_name) \
  friend class test_case_name##_##test_name##_Test
