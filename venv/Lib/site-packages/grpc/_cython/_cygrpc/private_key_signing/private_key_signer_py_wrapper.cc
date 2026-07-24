//
//
// Copyright 2026 gRPC authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

#include "src/python/grpcio/grpc/_cython/_cygrpc/private_key_signing/private_key_signer_py_wrapper.h"

#include <grpc/support/log.h>

#include <memory>

#include "Python.h"
#include "grpc/private_key_signer.h"
#include "absl/status/status.h"

namespace grpc_python {

std::shared_ptr<grpc_core::PrivateKeySigner> PrivateKeySignerPyWrapper::Create(
    PrivateKeySignerPyWrapper::SignWrapperForPy sign_py_wrapper,
    PyObject* py_user_sign_fn, PyObject* destroy_event) {
  PyGILState_STATE state = PyGILState_Ensure();
  Py_INCREF(py_user_sign_fn);
  Py_INCREF(destroy_event);
  PyGILState_Release(state);
  return std::make_shared<PrivateKeySignerPyWrapper>(
      sign_py_wrapper, py_user_sign_fn, destroy_event);
}

PrivateKeySignerPyWrapper::~PrivateKeySignerPyWrapper() {
  PyGILState_STATE state = PyGILState_Ensure();
  Py_DECREF(static_cast<PyObject*>(py_user_sign_fn_));
  // Python will stay alive until this event is set
  PyObject* result = PyObject_CallMethod(destroy_event_, "set", "()");
  // crash if result is nullptr? - discussing
  Py_XDECREF(result);
  PyGILState_Release(state);
}

std::variant<absl::StatusOr<std::string>,
             std::shared_ptr<grpc_core::PrivateKeySigner::AsyncSigningHandle>>
PrivateKeySignerPyWrapper::Sign(absl::string_view data_to_sign,
                                SignatureAlgorithm signature_algorithm,
                                OnSignComplete on_sign_complete) {
  auto completion_context =
      std::make_shared<CompletionContext>(std::move(on_sign_complete));

  PrivateKeySignerPyWrapperResult result = sign_py_wrapper_(
      data_to_sign, signature_algorithm, py_user_sign_fn_, completion_context);
  if (result.is_sync) {
    return result.sync_result;
  } else {
    auto handle = std::make_shared<AsyncSigningHandlePyWrapper>(
        result.async_result.cancel_wrapper,
        result.async_result.py_user_cancel_fn, std::move(completion_context));
    return handle;
  }
}

void PrivateKeySignerPyWrapper::Cancel(
    std::shared_ptr<AsyncSigningHandle> handle) {
  if (handle == nullptr) return;
  auto handle_impl =
      std::static_pointer_cast<AsyncSigningHandlePyWrapper>(handle);
  handle_impl->Cancel();
}

PrivateKeySignerPyWrapper::AsyncSigningHandlePyWrapper::
    ~AsyncSigningHandlePyWrapper() {
  PyGILState_STATE state = PyGILState_Ensure();
  Py_DECREF(py_user_cancel_fn_);
  PyGILState_Release(state);
}

void PrivateKeySignerPyWrapper::AsyncSigningHandlePyWrapper::Cancel() {
  if (cancel_py_wrapper_ != nullptr && py_user_cancel_fn_ != nullptr) {
    cancel_py_wrapper_(py_user_cancel_fn_);
  }
  completion_context_.reset();
}

std::string MakeStringForCython(const char* inp, size_t size) {
  return std::string(inp, size);
}

std::string MakeStringForCython(const char* inp) { return std::string(inp); }

}  // namespace grpc_python
