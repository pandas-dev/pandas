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

#if defined(__GNUC__)  // GCC and compatible compilers (clang, Intel ICC)
#  define ARROW_PREFETCH(addr) __builtin_prefetch(addr)
#elif defined(_MSC_VER)  // MSVC
#  if defined(ARROW_HAVE_SSE4_2) || defined(ARROW_HAVE_RUNTIME_SSE4_2)
#    include <nmmintrin.h>
#    define ARROW_PREFETCH(addr) _mm_prefetch((const char*)(addr), _MM_HINT_T0)
#  else
#    define ARROW_PREFETCH(addr)
#  endif
#else
#  define ARROW_PREFETCH(addr)
#endif
