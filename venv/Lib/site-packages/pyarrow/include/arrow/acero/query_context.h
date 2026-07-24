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

#include <string_view>

#include "arrow/acero/exec_plan.h"
#include "arrow/acero/task_util.h"
#include "arrow/acero/util.h"
#include "arrow/compute/exec.h"
#include "arrow/io/interfaces.h"
#include "arrow/util/async_util.h"
#include "arrow/util/type_fwd.h"

namespace arrow {

using compute::default_exec_context;
using io::IOContext;

namespace acero {

class ARROW_ACERO_EXPORT QueryContext {
 public:
  QueryContext(QueryOptions opts = {},
               ExecContext exec_context = *default_exec_context());

  Status Init(arrow::util::AsyncTaskScheduler* scheduler);

  const ::arrow::internal::CpuInfo* cpu_info() const;
  int64_t hardware_flags() const;
  const QueryOptions& options() const { return options_; }
  MemoryPool* memory_pool() const { return exec_context_.memory_pool(); }
  ::arrow::internal::Executor* executor() const { return exec_context_.executor(); }
  ExecContext* exec_context() { return &exec_context_; }
  IOContext* io_context() { return &io_context_; }
  TaskScheduler* scheduler() { return task_scheduler_.get(); }
  arrow::util::AsyncTaskScheduler* async_scheduler() { return async_scheduler_; }

  size_t GetThreadIndex();
  size_t max_concurrency() const;

  /// \brief Start an external task
  ///
  /// This should be avoided if possible.  It is kept in for now for legacy
  /// purposes.  This should be called before the external task is started.  If
  /// a valid future is returned then it should be marked complete when the
  /// external task has finished.
  ///
  /// \param name A name to give the task for traceability and debugging
  ///
  /// \return an invalid future if the plan has already ended, otherwise this
  ///         returns a future that must be completed when the external task
  ///         finishes.
  Result<Future<>> BeginExternalTask(std::string_view name);

  /// \brief Add a single function as a task to the query's task group
  ///        on the compute threadpool.
  ///
  /// \param fn The task to run. Takes no arguments and returns a Status.
  /// \param name A name to give the task for traceability and debugging
  void ScheduleTask(std::function<Status()> fn, std::string_view name);
  /// \brief Add a single function as a task to the query's task group
  ///        on the compute threadpool.
  ///
  /// \param fn The task to run. Takes the thread index and returns a Status.
  /// \param name A name to give the task for traceability and debugging
  void ScheduleTask(std::function<Status(size_t)> fn, std::string_view name);
  /// \brief Add a single function as a task to the query's task group on
  ///        the IO thread pool
  ///
  /// \param fn The task to run. Returns a status.
  /// \param name A name to give the task for traceability and debugging
  void ScheduleIOTask(std::function<Status()> fn, std::string_view name);

  // Register/Start TaskGroup is a way of performing a "Parallel For" pattern:
  // - The task function takes the thread index and the index of the task
  // - The on_finished function takes the thread index
  // Returns an integer ID that will be used to reference the task group in
  // StartTaskGroup. At runtime, call StartTaskGroup with the ID and the number of times
  // you'd like the task to be executed. The need to register a task group before use will
  // be removed after we rewrite the scheduler.
  /// \brief Register a "parallel for" task group with the scheduler
  ///
  /// \param task The function implementing the task. Takes the thread_index and
  ///             the task index.
  /// \param on_finished The function that gets run once all tasks have been completed.
  /// Takes the thread_index.
  ///
  /// Must be called inside of ExecNode::Init.
  int RegisterTaskGroup(std::function<Status(size_t, int64_t)> task,
                        std::function<Status(size_t)> on_finished);

  /// \brief Start the task group with the specified ID. This can only
  ///        be called once per task_group_id.
  ///
  /// \param task_group_id The ID  of the task group to run
  /// \param num_tasks The number of times to run the task
  Status StartTaskGroup(int task_group_id, int64_t num_tasks);

  // This is an RAII class for keeping track of in-flight file IO. Useful for getting
  // an estimate of memory use, and how much memory we expect to be freed soon.
  // Returned by ReportTempFileIO.
  struct [[nodiscard]] TempFileIOMark {
    QueryContext* ctx_;
    size_t bytes_;

    TempFileIOMark(QueryContext* ctx, size_t bytes) : ctx_(ctx), bytes_(bytes) {
      ctx_->in_flight_bytes_to_disk_.fetch_add(bytes_, std::memory_order_acquire);
    }

    ARROW_DISALLOW_COPY_AND_ASSIGN(TempFileIOMark);

    ~TempFileIOMark() {
      ctx_->in_flight_bytes_to_disk_.fetch_sub(bytes_, std::memory_order_release);
    }
  };

  TempFileIOMark ReportTempFileIO(size_t bytes) { return {this, bytes}; }

  size_t GetCurrentTempFileIO() { return in_flight_bytes_to_disk_.load(); }

 private:
  QueryOptions options_;
  // To be replaced with Acero-specific context once scheduler is done and
  // we don't need ExecContext for kernels
  ExecContext exec_context_;
  IOContext io_context_;

  arrow::util::AsyncTaskScheduler* async_scheduler_ = NULLPTR;
  std::unique_ptr<TaskScheduler> task_scheduler_ = TaskScheduler::Make();

  ThreadIndexer thread_indexer_;

  std::atomic<size_t> in_flight_bytes_to_disk_{0};
};
}  // namespace acero
}  // namespace arrow
