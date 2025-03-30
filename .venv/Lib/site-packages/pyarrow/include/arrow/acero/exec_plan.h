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

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "arrow/acero/type_fwd.h"
#include "arrow/acero/visibility.h"
#include "arrow/compute/api_vector.h"
#include "arrow/compute/exec.h"
#include "arrow/compute/ordering.h"
#include "arrow/type_fwd.h"
#include "arrow/util/future.h"
#include "arrow/util/macros.h"
#include "arrow/util/tracing.h"
#include "arrow/util/type_fwd.h"

namespace arrow {

using compute::ExecBatch;
using compute::ExecContext;
using compute::FunctionRegistry;
using compute::GetFunctionRegistry;
using compute::Ordering;
using compute::threaded_exec_context;

namespace acero {

/// \addtogroup acero-internals
/// @{

class ARROW_ACERO_EXPORT ExecPlan : public std::enable_shared_from_this<ExecPlan> {
 public:
  // This allows operators to rely on signed 16-bit indices
  static const uint32_t kMaxBatchSize = 1 << 15;
  using NodeVector = std::vector<ExecNode*>;

  virtual ~ExecPlan() = default;

  QueryContext* query_context();

  /// \brief retrieve the nodes in the plan
  const NodeVector& nodes() const;

  /// Make an empty exec plan
  static Result<std::shared_ptr<ExecPlan>> Make(
      QueryOptions options, ExecContext exec_context = *threaded_exec_context(),
      std::shared_ptr<const KeyValueMetadata> metadata = NULLPTR);

  static Result<std::shared_ptr<ExecPlan>> Make(
      ExecContext exec_context = *threaded_exec_context(),
      std::shared_ptr<const KeyValueMetadata> metadata = NULLPTR);

  static Result<std::shared_ptr<ExecPlan>> Make(
      QueryOptions options, ExecContext* exec_context,
      std::shared_ptr<const KeyValueMetadata> metadata = NULLPTR);

  static Result<std::shared_ptr<ExecPlan>> Make(
      ExecContext* exec_context,
      std::shared_ptr<const KeyValueMetadata> metadata = NULLPTR);

  ExecNode* AddNode(std::unique_ptr<ExecNode> node);

  template <typename Node, typename... Args>
  Node* EmplaceNode(Args&&... args) {
    std::unique_ptr<Node> node{new Node{std::forward<Args>(args)...}};
    auto out = node.get();
    AddNode(std::move(node));
    return out;
  }

  Status Validate();

  /// \brief Start producing on all nodes
  ///
  /// Nodes are started in reverse topological order, such that any node
  /// is started before all of its inputs.
  void StartProducing();

  /// \brief Stop producing on all nodes
  ///
  /// Triggers all sources to stop producing new data.  In order to cleanly stop the plan
  /// will continue to run any tasks that are already in progress.  The caller should
  /// still wait for `finished` to complete before destroying the plan.
  void StopProducing();

  /// \brief A future which will be marked finished when all tasks have finished.
  Future<> finished();

  /// \brief Return whether the plan has non-empty metadata
  bool HasMetadata() const;

  /// \brief Return the plan's attached metadata
  std::shared_ptr<const KeyValueMetadata> metadata() const;

  std::string ToString() const;
};

// Acero can be extended by providing custom implementations of ExecNode.  The methods
// below are documented in detail and provide careful instruction on how to fulfill the
// ExecNode contract.  It's suggested you familiarize yourself with the Acero
// documentation in the C++ user guide.
class ARROW_ACERO_EXPORT ExecNode {
 public:
  using NodeVector = std::vector<ExecNode*>;

  virtual ~ExecNode() = default;

  virtual const char* kind_name() const = 0;

  // The number of inputs expected by this node
  int num_inputs() const { return static_cast<int>(inputs_.size()); }

  /// This node's predecessors in the exec plan
  const NodeVector& inputs() const { return inputs_; }

  /// True if the plan has no output schema (is a sink)
  bool is_sink() const { return !output_schema_; }

  /// \brief Labels identifying the function of each input.
  const std::vector<std::string>& input_labels() const { return input_labels_; }

  /// This node's successor in the exec plan
  const ExecNode* output() const { return output_; }

  /// The datatypes for batches produced by this node
  const std::shared_ptr<Schema>& output_schema() const { return output_schema_; }

  /// This node's exec plan
  ExecPlan* plan() { return plan_; }

  /// \brief An optional label, for display and debugging
  ///
  /// There is no guarantee that this value is non-empty or unique.
  const std::string& label() const { return label_; }
  void SetLabel(std::string label) { label_ = std::move(label); }

  virtual Status Validate() const;

  /// \brief the ordering of the output batches
  ///
  /// This does not guarantee the batches will be emitted by this node
  /// in order.  Instead it guarantees that the batches will have their
  /// ExecBatch::index property set in a way that respects this ordering.
  ///
  /// In other words, given the ordering {{"x", SortOrder::Ascending}} we
  /// know that all values of x in a batch with index N will be less than
  /// or equal to all values of x in a batch with index N+k (assuming k > 0).
  /// Furthermore, we also know that values will be sorted within a batch.
  /// Any row N will have a value of x that is less than the value for
  /// any row N+k.
  ///
  /// Note that an ordering can be both Ordering::Unordered and Ordering::Implicit.
  /// A node's output should be marked Ordering::Unordered if the order is
  /// non-deterministic.  For example, a hash-join has no predictable output order.
  ///
  /// If the ordering is Ordering::Implicit then there is a meaningful order but that
  /// ordering is not represented by any column in the data.  The most common case for
  /// this is when reading data from an in-memory table.  The data has an implicit "row
  /// order" which is not necessarily represented in the data set.
  ///
  /// A filter or project node will not modify the ordering.  Nothing needs to be done
  /// other than ensure the index assigned to output batches is the same as the
  /// input batch that was mapped.
  ///
  /// Other nodes may introduce order.  For example, an order-by node will emit
  /// a brand new ordering independent of the input ordering.
  ///
  /// Finally, as described above, such as a hash-join or aggregation may may
  /// destroy ordering (although these nodes could also choose to establish a
  /// new ordering based on the hash keys).
  ///
  /// Some nodes will require an ordering.  For example, a fetch node or an
  /// asof join node will only function if the input data is ordered (for fetch
  /// it is enough to be implicitly ordered.  For an asof join the ordering must
  /// be explicit and compatible with the on key.)
  ///
  /// Nodes that maintain ordering should be careful to avoid introducing gaps
  /// in the batch index.  This may require emitting empty batches in order to
  /// maintain continuity.
  virtual const Ordering& ordering() const;

  /// Upstream API:
  /// These functions are called by input nodes that want to inform this node
  /// about an updated condition (a new input batch or an impending
  /// end of stream).
  ///
  /// Implementation rules:
  /// - these may be called anytime after StartProducing() has succeeded
  ///   (and even during or after StopProducing())
  /// - these may be called concurrently
  /// - these are allowed to call back into PauseProducing(), ResumeProducing()
  ///   and StopProducing()

  /// Transfer input batch to ExecNode
  ///
  /// A node will typically perform some kind of operation on the batch
  /// and then call InputReceived on its outputs with the result.
  ///
  /// Other nodes may need to accumulate some number of inputs before any
  /// output can be produced.  These nodes will add the batch to some kind
  /// of in-memory accumulation queue and return.
  virtual Status InputReceived(ExecNode* input, ExecBatch batch) = 0;

  /// Mark the inputs finished after the given number of batches.
  ///
  /// This may be called before all inputs are received.  This simply fixes
  /// the total number of incoming batches for an input, so that the ExecNode
  /// knows when it has received all input, regardless of order.
  virtual Status InputFinished(ExecNode* input, int total_batches) = 0;

  /// \brief Perform any needed initialization
  ///
  /// This hook performs any actions in between creation of ExecPlan and the call to
  /// StartProducing. An example could be Bloom filter pushdown. The order of ExecNodes
  /// that executes this method is undefined, but the calls are made synchronously.
  ///
  /// At this point a node can rely on all inputs & outputs (and the input schemas)
  /// being well defined.
  virtual Status Init();

  /// Lifecycle API:
  /// - start / stop to initiate and terminate production
  /// - pause / resume to apply backpressure
  ///
  /// Implementation rules:
  /// - StartProducing() should not recurse into the inputs, as it is
  ///   handled by ExecPlan::StartProducing()
  /// - PauseProducing(), ResumeProducing(), StopProducing() may be called
  ///   concurrently, potentially even before the call to StartProducing
  ///   has finished.
  /// - PauseProducing(), ResumeProducing(), StopProducing() may be called
  ///   by the downstream nodes' InputReceived(), InputFinished() methods
  ///
  /// StopProducing may be called due to an error, by the user (e.g. cancel), or
  /// because a node has all the data it needs (e.g. limit, top-k on sorted data).
  /// This means the method may be called multiple times and we have the following
  /// additional rules
  /// - StopProducing() must be idempotent
  /// - StopProducing() must be forwarded to inputs (this is needed for the limit/top-k
  ///     case because we may not be stopping the entire plan)

  // Right now, since synchronous calls happen in both directions (input to
  // output and then output to input), a node must be careful to be reentrant
  // against synchronous calls from its output, *and* also concurrent calls from
  // other threads.  The most reliable solution is to update the internal state
  // first, and notify outputs only at the end.
  //
  // Concurrent calls to PauseProducing and ResumeProducing can be hard to sequence
  // as they may travel at different speeds through the plan.
  //
  // For example, consider a resume that comes quickly after a pause.  If the source
  // receives the resume before the pause the source may think the destination is full
  // and halt production which would lead to deadlock.
  //
  // To resolve this a counter is sent for all calls to pause/resume.  Only the call with
  // the highest counter value is valid.  So if a call to PauseProducing(5) comes after
  // a call to ResumeProducing(6) then the source should continue producing.

  /// \brief Start producing
  ///
  /// This must only be called once.
  ///
  /// This is typically called automatically by ExecPlan::StartProducing().
  virtual Status StartProducing() = 0;

  /// \brief Pause producing temporarily
  ///
  /// \param output Pointer to the output that is full
  /// \param counter Counter used to sequence calls to pause/resume
  ///
  /// This call is a hint that an output node is currently not willing
  /// to receive data.
  ///
  /// This may be called any number of times.
  /// However, the node is still free to produce data (which may be difficult
  /// to prevent anyway if data is produced using multiple threads).
  virtual void PauseProducing(ExecNode* output, int32_t counter) = 0;

  /// \brief Resume producing after a temporary pause
  ///
  /// \param output Pointer to the output that is now free
  /// \param counter Counter used to sequence calls to pause/resume
  ///
  /// This call is a hint that an output node is willing to receive data again.
  ///
  /// This may be called any number of times.
  virtual void ResumeProducing(ExecNode* output, int32_t counter) = 0;

  /// \brief Stop producing new data
  ///
  /// If this node is a source then the source should stop generating data
  /// as quickly as possible.  If this node is not a source then there is typically
  /// nothing that needs to be done although a node may choose to start ignoring incoming
  /// data.
  ///
  /// This method will be called when an error occurs in the plan
  /// This method may also be called by the user if they wish to end a plan early
  /// Finally, this method may be called if a node determines it no longer needs any more
  /// input (for example, a limit node).
  ///
  /// This method may be called multiple times.
  ///
  /// This is not a pause.  There will be no way to start the source again after this has
  /// been called.
  virtual Status StopProducing();

  std::string ToString(int indent = 0) const;

 protected:
  ExecNode(ExecPlan* plan, NodeVector inputs, std::vector<std::string> input_labels,
           std::shared_ptr<Schema> output_schema);

  virtual Status StopProducingImpl() = 0;

  /// Provide extra info to include in the string representation.
  virtual std::string ToStringExtra(int indent = 0) const;

  std::atomic<bool> stopped_;
  ExecPlan* plan_;
  std::string label_;

  NodeVector inputs_;
  std::vector<std::string> input_labels_;

  std::shared_ptr<Schema> output_schema_;
  ExecNode* output_ = NULLPTR;
};

/// \brief An extensible registry for factories of ExecNodes
class ARROW_ACERO_EXPORT ExecFactoryRegistry {
 public:
  using Factory = std::function<Result<ExecNode*>(ExecPlan*, std::vector<ExecNode*>,
                                                  const ExecNodeOptions&)>;

  virtual ~ExecFactoryRegistry() = default;

  /// \brief Get the named factory from this registry
  ///
  /// will raise if factory_name is not found
  virtual Result<Factory> GetFactory(const std::string& factory_name) = 0;

  /// \brief Add a factory to this registry with the provided name
  ///
  /// will raise if factory_name is already in the registry
  virtual Status AddFactory(std::string factory_name, Factory factory) = 0;
};

/// The default registry, which includes built-in factories.
ARROW_ACERO_EXPORT
ExecFactoryRegistry* default_exec_factory_registry();

/// \brief Construct an ExecNode using the named factory
inline Result<ExecNode*> MakeExecNode(
    const std::string& factory_name, ExecPlan* plan, std::vector<ExecNode*> inputs,
    const ExecNodeOptions& options,
    ExecFactoryRegistry* registry = default_exec_factory_registry()) {
  ARROW_ASSIGN_OR_RAISE(auto factory, registry->GetFactory(factory_name));
  return factory(plan, std::move(inputs), options);
}

/// @}

/// \addtogroup acero-api
/// @{

/// \brief Helper class for declaring execution nodes
///
/// A Declaration represents an unconstructed ExecNode (and potentially an entire graph
/// since its inputs may also be Declarations)
///
/// A Declaration can be converted to a plan and executed using one of the
/// DeclarationToXyz methods.
///
/// For more direct control, a Declaration can be added to an existing execution
/// plan with Declaration::AddToPlan, which will recursively construct any inputs as
/// necessary.
struct ARROW_ACERO_EXPORT Declaration {
  using Input = std::variant<ExecNode*, Declaration>;

  Declaration() {}

  /// \brief construct a declaration
  /// \param factory_name the name of the exec node to construct.  The node must have
  ///                     been added to the exec node registry with this name.
  /// \param inputs the inputs to the node, these should be other declarations
  /// \param options options that control the behavior of the node.  You must use
  ///                the appropriate subclass.  For example, if `factory_name` is
  ///                "project" then `options` should be ProjectNodeOptions.
  /// \param label a label to give the node.  Can be used to distinguish it from other
  ///              nodes of the same type in the plan.
  Declaration(std::string factory_name, std::vector<Input> inputs,
              std::shared_ptr<ExecNodeOptions> options, std::string label)
      : factory_name{std::move(factory_name)},
        inputs{std::move(inputs)},
        options{std::move(options)},
        label{std::move(label)} {}

  template <typename Options>
  Declaration(std::string factory_name, std::vector<Input> inputs, Options options,
              std::string label)
      : Declaration{std::move(factory_name), std::move(inputs),
                    std::shared_ptr<ExecNodeOptions>(
                        std::make_shared<Options>(std::move(options))),
                    std::move(label)} {}

  template <typename Options>
  Declaration(std::string factory_name, std::vector<Input> inputs, Options options)
      : Declaration{std::move(factory_name), std::move(inputs), std::move(options),
                    /*label=*/""} {}

  template <typename Options>
  Declaration(std::string factory_name, Options options)
      : Declaration{std::move(factory_name), {}, std::move(options), /*label=*/""} {}

  template <typename Options>
  Declaration(std::string factory_name, Options options, std::string label)
      : Declaration{std::move(factory_name), {}, std::move(options), std::move(label)} {}

  /// \brief Convenience factory for the common case of a simple sequence of nodes.
  ///
  /// Each of decls will be appended to the inputs of the subsequent declaration,
  /// and the final modified declaration will be returned.
  ///
  /// Without this convenience factory, constructing a sequence would require explicit,
  /// difficult-to-read nesting:
  ///
  ///     Declaration{"n3",
  ///                   {
  ///                       Declaration{"n2",
  ///                                   {
  ///                                       Declaration{"n1",
  ///                                                   {
  ///                                                       Declaration{"n0", N0Opts{}},
  ///                                                   },
  ///                                                   N1Opts{}},
  ///                                   },
  ///                                   N2Opts{}},
  ///                   },
  ///                   N3Opts{}};
  ///
  /// An equivalent Declaration can be constructed more tersely using Sequence:
  ///
  ///     Declaration::Sequence({
  ///         {"n0", N0Opts{}},
  ///         {"n1", N1Opts{}},
  ///         {"n2", N2Opts{}},
  ///         {"n3", N3Opts{}},
  ///     });
  static Declaration Sequence(std::vector<Declaration> decls);

  /// \brief add the declaration to an already created execution plan
  /// \param plan the plan to add the node to
  /// \param registry the registry to use to lookup the node factory
  ///
  /// This method will recursively call AddToPlan on all of the declaration's inputs.
  /// This method is only for advanced use when the DeclarationToXyz methods are not
  /// sufficient.
  ///
  /// \return the instantiated execution node
  Result<ExecNode*> AddToPlan(ExecPlan* plan, ExecFactoryRegistry* registry =
                                                  default_exec_factory_registry()) const;

  // Validate a declaration
  bool IsValid(ExecFactoryRegistry* registry = default_exec_factory_registry()) const;

  /// \brief the name of the factory to use when creating a node
  std::string factory_name;
  /// \brief the declarations's inputs
  std::vector<Input> inputs;
  /// \brief options to control the behavior of the node
  std::shared_ptr<ExecNodeOptions> options;
  /// \brief a label to give the node in the plan
  std::string label;
};

/// \brief How to handle unaligned buffers
enum class UnalignedBufferHandling { kWarn, kIgnore, kReallocate, kError };

/// \brief get the default behavior of unaligned buffer handling
///
/// This is configurable via the ACERO_ALIGNMENT_HANDLING environment variable which
/// can be set to "warn", "ignore", "reallocate", or "error".  If the environment
/// variable is not set, or is set to an invalid value, this will return kWarn
UnalignedBufferHandling GetDefaultUnalignedBufferHandling();

/// \brief plan-wide options that can be specified when executing an execution plan
struct ARROW_ACERO_EXPORT QueryOptions {
  /// \brief Should the plan use a legacy batching strategy
  ///
  /// This is currently in place only to support the Scanner::ToTable
  /// method.  This method relies on batch indices from the scanner
  /// remaining consistent.  This is impractical in the ExecPlan which
  /// might slice batches as needed (e.g. for a join)
  ///
  /// However, it still works for simple plans and this is the only way
  /// we have at the moment for maintaining implicit order.
  bool use_legacy_batching = false;

  /// If the output has a meaningful order then sequence the output of the plan
  ///
  /// The default behavior (std::nullopt) will sequence output batches if there
  /// is a meaningful ordering in the final node and will emit batches immediately
  /// otherwise.
  ///
  /// If explicitly set to true then plan execution will fail if there is no
  /// meaningful ordering.  This can be useful to validate a query that should
  /// be emitting ordered results.
  ///
  /// If explicitly set to false then batches will be emit immediately even if there
  /// is a meaningful ordering.  This could cause batches to be emit out of order but
  /// may offer a small decrease to latency.
  std::optional<bool> sequence_output = std::nullopt;

  /// \brief should the plan use multiple background threads for CPU-intensive work
  ///
  /// If this is false then all CPU work will be done on the calling thread.  I/O tasks
  /// will still happen on the I/O executor and may be multi-threaded (but should not use
  /// significant CPU resources).
  ///
  /// Will be ignored if custom_cpu_executor is set
  bool use_threads = true;

  /// \brief custom executor to use for CPU-intensive work
  ///
  /// Must be null or remain valid for the duration of the plan.  If this is null then
  /// a default thread pool will be chosen whose behavior will be controlled by
  /// the `use_threads` option.
  ::arrow::internal::Executor* custom_cpu_executor = NULLPTR;

  /// \brief custom executor to use for IO work
  ///
  /// Must be null or remain valid for the duration of the plan.  If this is null then
  /// the global io thread pool will be chosen whose behavior will be controlled by
  /// the "ARROW_IO_THREADS" environment.
  ::arrow::internal::Executor* custom_io_executor = NULLPTR;

  /// \brief a memory pool to use for allocations
  ///
  /// Must remain valid for the duration of the plan.
  MemoryPool* memory_pool = default_memory_pool();

  /// \brief a function registry to use for the plan
  ///
  /// Must remain valid for the duration of the plan.
  FunctionRegistry* function_registry = GetFunctionRegistry();
  /// \brief the names of the output columns
  ///
  /// If this is empty then names will be generated based on the input columns
  ///
  /// If set then the number of names must equal the number of output columns
  std::vector<std::string> field_names;

  /// \brief Policy for unaligned buffers in source data
  ///
  /// Various compute functions and acero internals will type pun array
  /// buffers from uint8_t* to some kind of value type (e.g. we might
  /// cast to int32_t* to add two int32 arrays)
  ///
  /// If the buffer is poorly aligned (e.g. an int32 array is not aligned
  /// on a 4-byte boundary) then this is technically undefined behavior in C++.
  /// However, most modern compilers and CPUs are fairly tolerant of this
  /// behavior and nothing bad (beyond a small hit to performance) is likely
  /// to happen.
  ///
  /// Note that this only applies to source buffers.  All buffers allocated internally
  /// by Acero will be suitably aligned.
  ///
  /// If this field is set to kWarn then Acero will check if any buffers are unaligned
  /// and, if they are, will emit a warning.
  ///
  /// If this field is set to kReallocate then Acero will allocate a new, suitably aligned
  /// buffer and copy the contents from the old buffer into this new buffer.
  ///
  /// If this field is set to kError then Acero will gracefully abort the plan instead.
  ///
  /// If this field is set to kIgnore then Acero will not even check if the buffers are
  /// unaligned.
  ///
  /// If this field is not set then it will be treated as kWarn unless overridden
  /// by the ACERO_ALIGNMENT_HANDLING environment variable
  std::optional<UnalignedBufferHandling> unaligned_buffer_handling;
};

/// \brief Calculate the output schema of a declaration
///
/// This does not actually execute the plan.  This operation may fail if the
/// declaration represents an invalid plan (e.g. a project node with multiple inputs)
///
/// \param declaration A declaration describing an execution plan
/// \param function_registry The function registry to use for function execution.  If null
///                          then the default function registry will be used.
///
/// \return the schema that batches would have after going through the execution plan
ARROW_ACERO_EXPORT Result<std::shared_ptr<Schema>> DeclarationToSchema(
    const Declaration& declaration, FunctionRegistry* function_registry = NULLPTR);

/// \brief Create a string representation of a plan
///
/// This representation is for debug purposes only.
///
/// Conversion to a string may fail if the declaration represents an
/// invalid plan.
///
/// Use Substrait for complete serialization of plans
///
/// \param declaration A declaration describing an execution plan
/// \param function_registry The function registry to use for function execution.  If null
///                          then the default function registry will be used.
///
/// \return a string representation of the plan suitable for debugging output
ARROW_ACERO_EXPORT Result<std::string> DeclarationToString(
    const Declaration& declaration, FunctionRegistry* function_registry = NULLPTR);

/// \brief Utility method to run a declaration and collect the results into a table
///
/// \param declaration A declaration describing the plan to run
/// \param use_threads If `use_threads` is false then all CPU work will be done on the
///                    calling thread.  I/O tasks will still happen on the I/O executor
///                    and may be multi-threaded (but should not use significant CPU
///                    resources).
/// \param memory_pool The memory pool to use for allocations made while running the plan.
/// \param function_registry The function registry to use for function execution.  If null
///                          then the default function registry will be used.
///
/// This method will add a sink node to the declaration to collect results into a
/// table.  It will then create an ExecPlan from the declaration, start the exec plan,
/// block until the plan has finished, and return the created table.
ARROW_ACERO_EXPORT Result<std::shared_ptr<Table>> DeclarationToTable(
    Declaration declaration, bool use_threads = true,
    MemoryPool* memory_pool = default_memory_pool(),
    FunctionRegistry* function_registry = NULLPTR);

ARROW_ACERO_EXPORT Result<std::shared_ptr<Table>> DeclarationToTable(
    Declaration declaration, QueryOptions query_options);

/// \brief Asynchronous version of \see DeclarationToTable
///
/// \param declaration A declaration describing the plan to run
/// \param use_threads The behavior of use_threads is slightly different than the
///                    synchronous version since we cannot run synchronously on the
///                    calling thread. Instead, if use_threads=false then a new thread
///                    pool will be created with a single thread and this will be used for
///                    all compute work.
/// \param memory_pool The memory pool to use for allocations made while running the plan.
/// \param function_registry The function registry to use for function execution. If null
///                          then the default function registry will be used.
ARROW_ACERO_EXPORT Future<std::shared_ptr<Table>> DeclarationToTableAsync(
    Declaration declaration, bool use_threads = true,
    MemoryPool* memory_pool = default_memory_pool(),
    FunctionRegistry* function_registry = NULLPTR);

/// \brief Overload of \see DeclarationToTableAsync accepting a custom exec context
///
/// The executor must be specified (cannot be null) and must be kept alive until the
/// returned future finishes.
ARROW_ACERO_EXPORT Future<std::shared_ptr<Table>> DeclarationToTableAsync(
    Declaration declaration, ExecContext custom_exec_context);

/// \brief a collection of exec batches with a common schema
struct BatchesWithCommonSchema {
  std::vector<ExecBatch> batches;
  std::shared_ptr<Schema> schema;
};

/// \brief Utility method to run a declaration and collect the results into ExecBatch
/// vector
///
/// \see DeclarationToTable for details on threading & execution
ARROW_ACERO_EXPORT Result<BatchesWithCommonSchema> DeclarationToExecBatches(
    Declaration declaration, bool use_threads = true,
    MemoryPool* memory_pool = default_memory_pool(),
    FunctionRegistry* function_registry = NULLPTR);

ARROW_ACERO_EXPORT Result<BatchesWithCommonSchema> DeclarationToExecBatches(
    Declaration declaration, QueryOptions query_options);

/// \brief Asynchronous version of \see DeclarationToExecBatches
///
/// \see DeclarationToTableAsync for details on threading & execution
ARROW_ACERO_EXPORT Future<BatchesWithCommonSchema> DeclarationToExecBatchesAsync(
    Declaration declaration, bool use_threads = true,
    MemoryPool* memory_pool = default_memory_pool(),
    FunctionRegistry* function_registry = NULLPTR);

/// \brief Overload of \see DeclarationToExecBatchesAsync accepting a custom exec context
///
/// \see DeclarationToTableAsync for details on threading & execution
ARROW_ACERO_EXPORT Future<BatchesWithCommonSchema> DeclarationToExecBatchesAsync(
    Declaration declaration, ExecContext custom_exec_context);

/// \brief Utility method to run a declaration and collect the results into a vector
///
/// \see DeclarationToTable for details on threading & execution
ARROW_ACERO_EXPORT Result<std::vector<std::shared_ptr<RecordBatch>>> DeclarationToBatches(
    Declaration declaration, bool use_threads = true,
    MemoryPool* memory_pool = default_memory_pool(),
    FunctionRegistry* function_registry = NULLPTR);

ARROW_ACERO_EXPORT Result<std::vector<std::shared_ptr<RecordBatch>>> DeclarationToBatches(
    Declaration declaration, QueryOptions query_options);

/// \brief Asynchronous version of \see DeclarationToBatches
///
/// \see DeclarationToTableAsync for details on threading & execution
ARROW_ACERO_EXPORT Future<std::vector<std::shared_ptr<RecordBatch>>>
DeclarationToBatchesAsync(Declaration declaration, bool use_threads = true,
                          MemoryPool* memory_pool = default_memory_pool(),
                          FunctionRegistry* function_registry = NULLPTR);

/// \brief Overload of \see DeclarationToBatchesAsync accepting a custom exec context
///
/// \see DeclarationToTableAsync for details on threading & execution
ARROW_ACERO_EXPORT Future<std::vector<std::shared_ptr<RecordBatch>>>
DeclarationToBatchesAsync(Declaration declaration, ExecContext exec_context);

/// \brief Utility method to run a declaration and return results as a RecordBatchReader
///
/// If an exec context is not provided then a default exec context will be used based
/// on the value of `use_threads`.  If `use_threads` is false then the CPU executor will
/// be a serial executor and all CPU work will be done on the calling thread.  I/O tasks
/// will still happen on the I/O executor and may be multi-threaded.
///
/// If `use_threads` is false then all CPU work will happen during the calls to
/// RecordBatchReader::Next and no CPU work will happen in the background.  If
/// `use_threads` is true then CPU work will happen on the CPU thread pool and tasks may
/// run in between calls to RecordBatchReader::Next.  If the returned reader is not
/// consumed quickly enough then the plan will eventually pause as the backpressure queue
/// fills up.
///
/// If a custom exec context is provided then the value of `use_threads` will be ignored.
///
/// The returned RecordBatchReader can be closed early to cancel the computation of record
/// batches. In this case, only errors encountered by the computation may be reported. In
/// particular, no cancellation error may be reported.
ARROW_ACERO_EXPORT Result<std::unique_ptr<RecordBatchReader>> DeclarationToReader(
    Declaration declaration, bool use_threads = true,
    MemoryPool* memory_pool = default_memory_pool(),
    FunctionRegistry* function_registry = NULLPTR);

ARROW_ACERO_EXPORT Result<std::unique_ptr<RecordBatchReader>> DeclarationToReader(
    Declaration declaration, QueryOptions query_options);

/// \brief Utility method to run a declaration and ignore results
///
/// This can be useful when the data are consumed as part of the plan itself, for
/// example, when the plan ends with a write node.
///
/// \see DeclarationToTable for details on threading & execution
ARROW_ACERO_EXPORT Status
DeclarationToStatus(Declaration declaration, bool use_threads = true,
                    MemoryPool* memory_pool = default_memory_pool(),
                    FunctionRegistry* function_registry = NULLPTR);

ARROW_ACERO_EXPORT Status DeclarationToStatus(Declaration declaration,
                                              QueryOptions query_options);

/// \brief Asynchronous version of \see DeclarationToStatus
///
/// This can be useful when the data are consumed as part of the plan itself, for
/// example, when the plan ends with a write node.
///
/// \see DeclarationToTableAsync for details on threading & execution
ARROW_ACERO_EXPORT Future<> DeclarationToStatusAsync(
    Declaration declaration, bool use_threads = true,
    MemoryPool* memory_pool = default_memory_pool(),
    FunctionRegistry* function_registry = NULLPTR);

/// \brief Overload of \see DeclarationToStatusAsync accepting a custom exec context
///
/// \see DeclarationToTableAsync for details on threading & execution
ARROW_ACERO_EXPORT Future<> DeclarationToStatusAsync(Declaration declaration,
                                                     ExecContext exec_context);

/// @}

/// \brief Wrap an ExecBatch generator in a RecordBatchReader.
///
/// The RecordBatchReader does not impose any ordering on emitted batches.
ARROW_ACERO_EXPORT
std::shared_ptr<RecordBatchReader> MakeGeneratorReader(
    std::shared_ptr<Schema>, std::function<Future<std::optional<ExecBatch>>()>,
    MemoryPool*);

constexpr int kDefaultBackgroundMaxQ = 32;
constexpr int kDefaultBackgroundQRestart = 16;

/// \brief Make a generator of RecordBatchReaders
///
/// Useful as a source node for an Exec plan
ARROW_ACERO_EXPORT
Result<std::function<Future<std::optional<ExecBatch>>()>> MakeReaderGenerator(
    std::shared_ptr<RecordBatchReader> reader, arrow::internal::Executor* io_executor,
    int max_q = kDefaultBackgroundMaxQ, int q_restart = kDefaultBackgroundQRestart);

}  // namespace acero
}  // namespace arrow
