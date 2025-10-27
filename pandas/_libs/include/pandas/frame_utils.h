/*
Copyright (c) 2025, Pandas Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Return whether or not the object is a local in the caller's frame
int is_local_in_caller_frame_impl(PyObject *object);
