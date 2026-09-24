// Copyright 2018 Ulf Adams
//
// The contents of this file may be used under the terms of the Apache License,
// Version 2.0.
//
//    (See accompanying file LICENSE-Apache or copy at
//     http://www.apache.org/licenses/LICENSE-2.0)
//
// Alternatively, the contents of this file may be used under the terms of
// the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE-Boost or copy at
//     https://www.boost.org/LICENSE_1_0.txt)
//
// Unless required by applicable law or agreed to in writing, this software
// is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.
#ifndef RYU_H
#define RYU_H

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>

int d2s_buffered_n(double f, char *result);
void d2s_buffered(double f, char *result);
char *d2s(double f);

int f2s_buffered_n(float f, char *result);
void f2s_buffered(float f, char *result);
char *f2s(float f);

int d2fixed_buffered_n(double d, uint32_t precision, char *result);
void d2fixed_buffered(double d, uint32_t precision, char *result);
char *d2fixed(double d, uint32_t precision);

int d2exp_buffered_n(double d, uint32_t precision, char *result);
void d2exp_buffered(double d, uint32_t precision, char *result);
char *d2exp(double d, uint32_t precision);

#ifdef __cplusplus
}
#endif

#endif // RYU_H
