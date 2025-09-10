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

#include "_blocking_impl.h"

#if defined(_WIN32)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#define READ _read
#include <errno.h>
#include <fcntl.h>
#include <io.h>
#include <windows.h>
#else
#define READ read
#include <fcntl.h>
#include <pthread.h>
#include <unistd.h>
#endif

#include <csignal>
#include <cstring>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>

namespace pyadbc_driver_manager {

// This is somewhat derived from io_util.cc in arrow, but that implementation
// isn't easily used outside of Arrow's monolith.
namespace {
static std::once_flag kInitOnce;
// We may encounter errors below that we can't do anything about. Use this to
// print out an error, once.
static std::once_flag kWarnOnce;
// This thread reads from a pipe forever.  Whenever it reads something, it
// calls the callback below.
static std::thread kCancelThread;

static std::mutex cancel_mutex;
// This callback is registered by the Python side; basically it will call
// cancel() on an ADBC object.
static void (*cancel_callback)(void*) = nullptr;
// Callback state (a pointer to the ADBC PyObject).
static void* cancel_callback_data = nullptr;
// A nonblocking self-pipe.
static int pipe[2];
#if defined(_WIN32)
void (*old_sigint)(int);
#else
// The old signal handler (most likely Python's).
struct sigaction old_sigint;
// Our signal handler (below).
struct sigaction our_sigint;
#endif

std::string MakePipe() {
  int rc = 0;
#if defined(__linux__) && defined(__GLIBC__)
  rc = pipe2(pipe, O_CLOEXEC);
#elif defined(_WIN32)
  rc = _pipe(pipe, 4096, _O_BINARY);
#else
  rc = ::pipe(pipe);
#endif

  if (rc != 0) {
    return std::strerror(errno);
  }

#if (!defined(__linux__) || !defined(__GLIBC__)) && !defined(_WIN32)
  {
    int flags = fcntl(pipe[0], F_GETFD, 0);
    if (flags < 0) {
      return std::strerror(errno);
    }
    rc = fcntl(pipe[0], F_SETFD, flags | FD_CLOEXEC);
    if (rc < 0) {
      return std::strerror(errno);
    }

    flags = fcntl(pipe[1], F_GETFD, 0);
    if (flags < 0) {
      return std::strerror(errno);
    }
    rc = fcntl(pipe[1], F_SETFD, flags | FD_CLOEXEC);
    if (rc < 0) {
      return std::strerror(errno);
    }
  }
#endif

  // Make the write side nonblocking (the read side should stay blocking!)
#if defined(_WIN32)
  const auto handle = reinterpret_cast<HANDLE>(_get_osfhandle(pipe[1]));
  DWORD mode = PIPE_NOWAIT;
  if (!SetNamedPipeHandleState(handle, &mode, nullptr, nullptr)) {
    DWORD last_error = GetLastError();
    LPVOID message;

    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                      FORMAT_MESSAGE_IGNORE_INSERTS,
                  /*lpSource=*/nullptr, last_error,
                  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                  reinterpret_cast<LPSTR>(&message), /*nSize=*/0, /*Arguments=*/nullptr);

    std::string buffer = "(";
    buffer += std::to_string(last_error);
    buffer += ") ";
    buffer += reinterpret_cast<char*>(message);
    LocalFree(message);
    return buffer;
  }
#else
  {
    int flags = fcntl(pipe[1], F_GETFL, 0);
    if (flags < 0) {
      return std::strerror(errno);
    }
    rc = fcntl(pipe[1], F_SETFL, flags | O_NONBLOCK);
    if (rc < 0) {
      return std::strerror(errno);
    }
  }
#endif

  return "";
}

void InterruptThread() {
#if defined(__APPLE__)
  pthread_setname_np("AdbcInterrupt");
#endif

  while (true) {
    char buf = 0;
    // Anytime something is written to the pipe, attempt to call the callback
    auto bytes_read = READ(pipe[0], &buf, 1);
    if (bytes_read < 0) {
      if (errno == EINTR) continue;

      // XXX: we failed reading from the pipe
      std::string message = std::strerror(errno);
      std::call_once(kWarnOnce, [&]() {
        std::cerr << "adbc_driver_manager (native code): error handling interrupt: "
                  << message << std::endl;
      });
    } else if (bytes_read > 0) {
      // Save the callback locally instead of calling it under the lock, since
      // otherwise we may deadlock with the Python side trying to call us
      void (*local_callback)(void*) = nullptr;
      void* local_callback_data = nullptr;

      {
        std::lock_guard<std::mutex> lock(cancel_mutex);
        if (cancel_callback != nullptr) {
          local_callback = cancel_callback;
          local_callback_data = cancel_callback_data;
        }
        cancel_callback = nullptr;
        cancel_callback_data = nullptr;
      }

      if (local_callback != nullptr) {
        local_callback(local_callback_data);
      }
    }
  }
}

// We can't do much about failures here, so ignore the result.  If the pipe is
// full, that's fine; it just means the thread has fallen behind in processing
// earlier interrupts.
void SigintHandler(int) {
#if defined(_WIN32)
  (void)_write(pipe[1], "X", 1);
#else
  (void)write(pipe[1], "X", 1);
#endif
}

}  // namespace

std::string InitBlockingCallback() {
  std::string error;
  std::call_once(kInitOnce, [&]() {
    error = MakePipe();
    if (!error.empty()) {
      return;
    }

#if !defined(_WIN32)
    our_sigint.sa_handler = &SigintHandler;
    our_sigint.sa_flags = 0;
    sigemptyset(&our_sigint.sa_mask);
#endif

    kCancelThread = std::thread(InterruptThread);
#if defined(__linux__)
    pthread_setname_np(kCancelThread.native_handle(), "AdbcInterrupt");
#endif
    kCancelThread.detach();
  });
  return error;
}

std::string SetBlockingCallback(void (*callback)(void*), void* data) {
  std::lock_guard<std::mutex> lock(cancel_mutex);
  cancel_callback = callback;
  cancel_callback_data = data;

#if defined(_WIN32)
  if (old_sigint == nullptr) {
    old_sigint = signal(SIGINT, &SigintHandler);
    if (old_sigint == SIG_ERR) {
      old_sigint = nullptr;
      return std::strerror(errno);
    }
  }
#else
  // Don't set the handler again if we're somehow called twice
  if (old_sigint.sa_handler == nullptr && old_sigint.sa_sigaction == nullptr) {
    int rc = sigaction(SIGINT, &our_sigint, &old_sigint);
    if (rc != 0) {
      return std::strerror(errno);
    }
  }
#endif
  return "";
}

std::string ClearBlockingCallback() {
  std::lock_guard<std::mutex> lock(cancel_mutex);
  cancel_callback = nullptr;
  cancel_callback_data = nullptr;

#if defined(_WIN32)
  if (old_sigint != nullptr) {
    auto rc = signal(SIGINT, old_sigint);
    old_sigint = nullptr;
    if (rc == SIG_ERR) {
      return std::strerror(errno);
    }
  }
#else
  if (old_sigint.sa_handler != nullptr || old_sigint.sa_sigaction != nullptr) {
    int rc = sigaction(SIGINT, &old_sigint, nullptr);
    std::memset(&old_sigint, 0, sizeof(old_sigint));
    if (rc != 0) {
      return std::strerror(errno);
    }
  }
#endif
  return "";
}

}  // namespace pyadbc_driver_manager
