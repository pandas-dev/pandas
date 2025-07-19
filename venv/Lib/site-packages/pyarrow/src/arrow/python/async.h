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

#include <utility>

#include "arrow/python/common.h"
#include "arrow/status.h"
#include "arrow/util/future.h"

namespace arrow::py {

/// \brief Bind a Python callback to an arrow::Future.
///
/// If the Future finishes successfully, py_wrapper is called with its
/// result value and should return a PyObject*. If py_wrapper is successful,
/// py_cb is called with its return value.
///
/// If either the Future or py_wrapper fails, py_cb is called with the
/// associated Python exception.
///
/// \param future The future to bind to.
/// \param py_cb The Python callback function. Will be passed the result of
///   py_wrapper, or a Python exception if the future failed or one was
///   raised by py_wrapper.
/// \param py_wrapper A function (likely defined in Cython) to convert the C++
///   result of the future to a Python object.
template <typename T, typename PyWrapper = PyObject* (*)(T)>
void BindFuture(Future<T> future, PyObject* py_cb, PyWrapper py_wrapper) {
  Py_INCREF(py_cb);
  OwnedRefNoGIL cb_ref(py_cb);

  auto future_cb = [cb_ref = std::move(cb_ref),
                    py_wrapper = std::move(py_wrapper)](Result<T> result) {
    SafeCallIntoPythonVoid([&]() {
      OwnedRef py_value_or_exc{WrapResult(std::move(result), std::move(py_wrapper))};
      Py_XDECREF(
          PyObject_CallFunctionObjArgs(cb_ref.obj(), py_value_or_exc.obj(), NULLPTR));
      ARROW_WARN_NOT_OK(CheckPyError(), "Internal error in async call");
    });
  };
  future.AddCallback(std::move(future_cb));
}

}  // namespace arrow::py
