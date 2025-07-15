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

#include "arrow/python/udf.h"

#include "arrow/array/array_nested.h"
#include "arrow/array/builder_base.h"
#include "arrow/buffer_builder.h"
#include "arrow/compute/api_aggregate.h"
#include "arrow/compute/api_vector.h"
#include "arrow/compute/function.h"
#include "arrow/compute/kernel.h"
#include "arrow/compute/row/grouper.h"
#include "arrow/python/common.h"
#include "arrow/python/vendored/pythoncapi_compat.h"
#include "arrow/table.h"
#include "arrow/util/checked_cast.h"
#include "arrow/util/logging.h"

namespace arrow {
using compute::ExecSpan;
using compute::Grouper;
using compute::KernelContext;
using compute::KernelState;
using internal::checked_cast;

namespace py {
namespace {

struct PythonUdfKernelState : public compute::KernelState {
  // NOTE: this KernelState constructor doesn't require the GIL.
  // If it did, the corresponding KernelInit::operator() should be wrapped
  // within SafeCallIntoPython (GH-43487).
  explicit PythonUdfKernelState(std::shared_ptr<OwnedRefNoGIL> function)
      : function(std::move(function)) {}

  std::shared_ptr<OwnedRefNoGIL> function;
};

struct PythonUdfKernelInit {
  explicit PythonUdfKernelInit(std::shared_ptr<OwnedRefNoGIL> function)
      : function(std::move(function)) {}

  Result<std::unique_ptr<compute::KernelState>> operator()(
      compute::KernelContext*, const compute::KernelInitArgs&) {
    return std::make_unique<PythonUdfKernelState>(function);
  }

  std::shared_ptr<OwnedRefNoGIL> function;
};

struct ScalarUdfAggregator : public compute::KernelState {
  virtual Status Consume(compute::KernelContext* ctx, const compute::ExecSpan& batch) = 0;
  virtual Status MergeFrom(compute::KernelContext* ctx, compute::KernelState&& src) = 0;
  virtual Status Finalize(compute::KernelContext* ctx, Datum* out) = 0;
};

struct HashUdfAggregator : public compute::KernelState {
  virtual Status Resize(KernelContext* ctx, int64_t size) = 0;
  virtual Status Consume(KernelContext* ctx, const ExecSpan& batch) = 0;
  virtual Status Merge(KernelContext* ct, KernelState&& other, const ArrayData&) = 0;
  virtual Status Finalize(KernelContext* ctx, Datum* out) = 0;
};

Status AggregateUdfConsume(compute::KernelContext* ctx, const compute::ExecSpan& batch) {
  return checked_cast<ScalarUdfAggregator*>(ctx->state())->Consume(ctx, batch);
}

Status AggregateUdfMerge(compute::KernelContext* ctx, compute::KernelState&& src,
                         compute::KernelState* dst) {
  return checked_cast<ScalarUdfAggregator*>(dst)->MergeFrom(ctx, std::move(src));
}

Status AggregateUdfFinalize(compute::KernelContext* ctx, arrow::Datum* out) {
  return checked_cast<ScalarUdfAggregator*>(ctx->state())->Finalize(ctx, out);
}

Status HashAggregateUdfResize(KernelContext* ctx, int64_t size) {
  return checked_cast<HashUdfAggregator*>(ctx->state())->Resize(ctx, size);
}

Status HashAggregateUdfConsume(KernelContext* ctx, const ExecSpan& batch) {
  return checked_cast<HashUdfAggregator*>(ctx->state())->Consume(ctx, batch);
}

Status HashAggregateUdfMerge(KernelContext* ctx, KernelState&& src,
                             const ArrayData& group_id_mapping) {
  return checked_cast<HashUdfAggregator*>(ctx->state())
      ->Merge(ctx, std::move(src), group_id_mapping);
}

Status HashAggregateUdfFinalize(KernelContext* ctx, Datum* out) {
  return checked_cast<HashUdfAggregator*>(ctx->state())->Finalize(ctx, out);
}

struct PythonTableUdfKernelInit {
  PythonTableUdfKernelInit(std::shared_ptr<OwnedRefNoGIL> function_maker,
                           UdfWrapperCallback cb)
      : function_maker(std::move(function_maker)), cb(std::move(cb)) {}

  Result<std::unique_ptr<compute::KernelState>> operator()(
      compute::KernelContext* ctx, const compute::KernelInitArgs&) {
    return SafeCallIntoPython(
        [this, ctx]() -> Result<std::unique_ptr<compute::KernelState>> {
          UdfContext udf_context{ctx->memory_pool(), /*batch_length=*/0};
          OwnedRef empty_tuple(PyTuple_New(0));
          auto function = std::make_shared<OwnedRefNoGIL>(
              cb(function_maker->obj(), udf_context, empty_tuple.obj()));
          RETURN_NOT_OK(CheckPyError());
          if (!PyCallable_Check(function->obj())) {
            return Status::TypeError("Expected a callable Python object.");
          }
          return std::make_unique<PythonUdfKernelState>(std::move(function));
        });
  }

  std::shared_ptr<OwnedRefNoGIL> function_maker;
  UdfWrapperCallback cb;
};

struct PythonUdfScalarAggregatorImpl : public ScalarUdfAggregator {
  PythonUdfScalarAggregatorImpl(std::shared_ptr<OwnedRefNoGIL> function,
                                UdfWrapperCallback cb,
                                std::vector<std::shared_ptr<DataType>> input_types,
                                std::shared_ptr<DataType> output_type)
      : function(std::move(function)),
        cb(std::move(cb)),
        output_type(std::move(output_type)) {
    std::vector<std::shared_ptr<Field>> fields;
    for (size_t i = 0; i < input_types.size(); i++) {
      fields.push_back(field("", input_types[i]));
    }
    input_schema = schema(std::move(fields));
  };

  Status Consume(compute::KernelContext* ctx, const compute::ExecSpan& batch) override {
    ARROW_ASSIGN_OR_RAISE(
        auto rb, batch.ToExecBatch().ToRecordBatch(input_schema, ctx->memory_pool()));
    values.push_back(std::move(rb));
    return Status::OK();
  }

  Status MergeFrom(compute::KernelContext* ctx, compute::KernelState&& src) override {
    auto& other_values = checked_cast<PythonUdfScalarAggregatorImpl&>(src).values;
    values.insert(values.end(), std::make_move_iterator(other_values.begin()),
                  std::make_move_iterator(other_values.end()));

    other_values.erase(other_values.begin(), other_values.end());
    return Status::OK();
  }

  Status Finalize(compute::KernelContext* ctx, Datum* out) override {
    auto state =
        arrow::internal::checked_cast<PythonUdfScalarAggregatorImpl*>(ctx->state());
    const int num_args = input_schema->num_fields();

    // Note: The way that batches are concatenated together
    // would result in using double amount of the memory.
    // This is OK for now because non decomposable aggregate
    // UDF is supposed to be used with segmented aggregation
    // where the size of the segment is more or less constant
    // so doubling that is not a big deal. This can be also
    // improved in the future to use more efficient way to
    // concatenate.
    ARROW_ASSIGN_OR_RAISE(auto table,
                          arrow::Table::FromRecordBatches(input_schema, values));
    ARROW_ASSIGN_OR_RAISE(table, table->CombineChunks(ctx->memory_pool()));
    UdfContext udf_context{ctx->memory_pool(), table->num_rows()};

    if (table->num_rows() == 0) {
      return Status::Invalid("Finalized is called with empty inputs");
    }

    RETURN_NOT_OK(SafeCallIntoPython([&] {
      std::unique_ptr<OwnedRef> result;
      OwnedRef arg_tuple(PyTuple_New(num_args));
      RETURN_NOT_OK(CheckPyError());

      for (int arg_id = 0; arg_id < num_args; arg_id++) {
        // Since we combined chunks there is only one chunk
        std::shared_ptr<Array> c_data = table->column(arg_id)->chunk(0);
        PyObject* data = wrap_array(c_data);
        PyTuple_SetItem(arg_tuple.obj(), arg_id, data);
      }
      result =
          std::make_unique<OwnedRef>(cb(function->obj(), udf_context, arg_tuple.obj()));
      RETURN_NOT_OK(CheckPyError());
      // unwrapping the output for expected output type
      if (is_scalar(result->obj())) {
        ARROW_ASSIGN_OR_RAISE(std::shared_ptr<Scalar> val, unwrap_scalar(result->obj()));
        if (*output_type != *val->type) {
          return Status::TypeError("Expected output datatype ", output_type->ToString(),
                                   ", but function returned datatype ",
                                   val->type->ToString());
        }
        out->value = std::move(val);
        return Status::OK();
      }
      return Status::TypeError("Unexpected output type: ",
                               Py_TYPE(result->obj())->tp_name, " (expected Scalar)");
    }));
    return Status::OK();
  }

  std::shared_ptr<OwnedRefNoGIL> function;
  UdfWrapperCallback cb;
  std::vector<std::shared_ptr<RecordBatch>> values;
  std::shared_ptr<Schema> input_schema;
  std::shared_ptr<DataType> output_type;
};

struct PythonUdfHashAggregatorImpl : public HashUdfAggregator {
  PythonUdfHashAggregatorImpl(std::shared_ptr<OwnedRefNoGIL> function,
                              UdfWrapperCallback cb,
                              std::vector<std::shared_ptr<DataType>> input_types,
                              std::shared_ptr<DataType> output_type)
      : function(std::move(function)),
        cb(std::move(cb)),
        output_type(std::move(output_type)) {
    std::vector<std::shared_ptr<Field>> fields;
    fields.reserve(input_types.size());
    for (size_t i = 0; i < input_types.size(); i++) {
      fields.push_back(field("", input_types[i]));
    }
    input_schema = schema(std::move(fields));
  };

  // same as ApplyGrouping in partition.cc
  // replicated the code here to avoid complicating the dependencies
  static Result<RecordBatchVector> ApplyGroupings(
      const ListArray& groupings, const std::shared_ptr<RecordBatch>& batch) {
    ARROW_ASSIGN_OR_RAISE(Datum sorted,
                          compute::Take(batch, groupings.data()->child_data[0]));

    const auto& sorted_batch = *sorted.record_batch();

    RecordBatchVector out(static_cast<size_t>(groupings.length()));
    for (size_t i = 0; i < out.size(); ++i) {
      out[i] = sorted_batch.Slice(groupings.value_offset(i), groupings.value_length(i));
    }

    return out;
  }

  Status Resize(KernelContext* ctx, int64_t new_num_groups) override {
    // We only need to change num_groups in resize
    // similar to other hash aggregate kernels
    num_groups = new_num_groups;
    return Status::OK();
  }

  Status Consume(KernelContext* ctx, const ExecSpan& batch) override {
    ARROW_ASSIGN_OR_RAISE(
        std::shared_ptr<RecordBatch> rb,
        batch.ToExecBatch().ToRecordBatch(input_schema, ctx->memory_pool()));

    // This is similar to GroupedListImpl
    // last array is the group id
    const ArraySpan& groups_array_data = batch[batch.num_values() - 1].array;
    DCHECK_EQ(groups_array_data.offset, 0);
    int64_t batch_num_values = groups_array_data.length;
    const auto* batch_groups = groups_array_data.GetValues<uint32_t>(1);
    RETURN_NOT_OK(groups.Append(batch_groups, batch_num_values));
    values.push_back(std::move(rb));
    num_values += batch_num_values;
    return Status::OK();
  }
  Status Merge(KernelContext* ctx, KernelState&& other_state,
               const ArrayData& group_id_mapping) override {
    // This is similar to GroupedListImpl
    auto& other = checked_cast<PythonUdfHashAggregatorImpl&>(other_state);
    auto& other_values = other.values;
    const uint32_t* other_raw_groups = other.groups.data();
    values.insert(values.end(), std::make_move_iterator(other_values.begin()),
                  std::make_move_iterator(other_values.end()));

    auto g = group_id_mapping.GetValues<uint32_t>(1);
    for (uint32_t other_g = 0; static_cast<int64_t>(other_g) < other.num_values;
         ++other_g) {
      // Different state can have different group_id mappings, so we
      // need to translate the ids
      RETURN_NOT_OK(groups.Append(g[other_raw_groups[other_g]]));
    }

    num_values += other.num_values;
    return Status::OK();
  }

  Status Finalize(KernelContext* ctx, Datum* out) override {
    // Exclude the last column which is the group id
    const int num_args = input_schema->num_fields() - 1;

    ARROW_ASSIGN_OR_RAISE(auto groups_buffer, groups.Finish());
    ARROW_ASSIGN_OR_RAISE(auto groupings,
                          Grouper::MakeGroupings(UInt32Array(num_values, groups_buffer),
                                                 static_cast<uint32_t>(num_groups)));

    ARROW_ASSIGN_OR_RAISE(auto table,
                          arrow::Table::FromRecordBatches(input_schema, values));
    ARROW_ASSIGN_OR_RAISE(auto rb, table->CombineChunksToBatch(ctx->memory_pool()));
    UdfContext udf_context{ctx->memory_pool(), table->num_rows()};

    if (rb->num_rows() == 0) {
      *out = Datum();
      return Status::OK();
    }

    ARROW_ASSIGN_OR_RAISE(RecordBatchVector rbs, ApplyGroupings(*groupings, rb));

    return SafeCallIntoPython([&] {
      ARROW_ASSIGN_OR_RAISE(std::unique_ptr<ArrayBuilder> builder,
                            MakeBuilder(output_type, ctx->memory_pool()));
      for (auto& group_rb : rbs) {
        std::unique_ptr<OwnedRef> result;
        OwnedRef arg_tuple(PyTuple_New(num_args));
        RETURN_NOT_OK(CheckPyError());

        for (int arg_id = 0; arg_id < num_args; arg_id++) {
          // Since we combined chunks there is only one chunk
          std::shared_ptr<Array> c_data = group_rb->column(arg_id);
          PyObject* data = wrap_array(c_data);
          PyTuple_SetItem(arg_tuple.obj(), arg_id, data);
        }

        result =
            std::make_unique<OwnedRef>(cb(function->obj(), udf_context, arg_tuple.obj()));
        RETURN_NOT_OK(CheckPyError());

        // unwrapping the output for expected output type
        if (is_scalar(result->obj())) {
          ARROW_ASSIGN_OR_RAISE(std::shared_ptr<Scalar> val,
                                unwrap_scalar(result->obj()));
          if (*output_type != *val->type) {
            return Status::TypeError("Expected output datatype ", output_type->ToString(),
                                     ", but function returned datatype ",
                                     val->type->ToString());
          }
          ARROW_RETURN_NOT_OK(builder->AppendScalar(std::move(*val)));
        } else {
          return Status::TypeError("Unexpected output type: ",
                                   Py_TYPE(result->obj())->tp_name, " (expected Scalar)");
        }
      }
      ARROW_ASSIGN_OR_RAISE(auto result, builder->Finish());
      out->value = std::move(result->data());
      return Status::OK();
    });
  }

  std::shared_ptr<OwnedRefNoGIL> function;
  UdfWrapperCallback cb;
  // Accumulated input batches
  std::vector<std::shared_ptr<RecordBatch>> values;
  // Group ids - extracted from the last column from the batch
  TypedBufferBuilder<uint32_t> groups;
  int64_t num_groups = 0;
  int64_t num_values = 0;
  std::shared_ptr<Schema> input_schema;
  std::shared_ptr<DataType> output_type;
};

struct PythonUdf : public PythonUdfKernelState {
  PythonUdf(std::shared_ptr<OwnedRefNoGIL> function, UdfWrapperCallback cb,
            std::vector<TypeHolder> input_types, compute::OutputType output_type)
      : PythonUdfKernelState(std::move(function)),
        cb(std::move(cb)),
        input_types(std::move(input_types)),
        output_type(std::move(output_type)) {}

  UdfWrapperCallback cb;
  std::vector<TypeHolder> input_types;
  compute::OutputType output_type;
  TypeHolder resolved_type;

  Result<TypeHolder> ResolveType(compute::KernelContext* ctx,
                                 const std::vector<TypeHolder>& types) {
    if (input_types == types) {
      if (!resolved_type) {
        ARROW_ASSIGN_OR_RAISE(resolved_type, output_type.Resolve(ctx, input_types));
      }
      return resolved_type;
    }
    return output_type.Resolve(ctx, types);
  }

  Status Exec(compute::KernelContext* ctx, const compute::ExecSpan& batch,
              compute::ExecResult* out) {
    auto state = arrow::internal::checked_cast<PythonUdfKernelState*>(ctx->state());
    PyObject* function = state->function->obj();
    const int num_args = batch.num_values();
    UdfContext udf_context{ctx->memory_pool(), batch.length};

    OwnedRef arg_tuple(PyTuple_New(num_args));
    RETURN_NOT_OK(CheckPyError());
    for (int arg_id = 0; arg_id < num_args; arg_id++) {
      if (batch[arg_id].is_scalar()) {
        std::shared_ptr<Scalar> c_data = batch[arg_id].scalar->GetSharedPtr();
        PyObject* data = wrap_scalar(c_data);
        PyTuple_SetItem(arg_tuple.obj(), arg_id, data);
      } else {
        std::shared_ptr<Array> c_data = batch[arg_id].array.ToArray();
        PyObject* data = wrap_array(c_data);
        PyTuple_SetItem(arg_tuple.obj(), arg_id, data);
      }
    }

    OwnedRef result(cb(function, udf_context, arg_tuple.obj()));
    RETURN_NOT_OK(CheckPyError());
    // unwrapping the output for expected output type
    if (is_array(result.obj())) {
      ARROW_ASSIGN_OR_RAISE(std::shared_ptr<Array> val, unwrap_array(result.obj()));
      ARROW_ASSIGN_OR_RAISE(TypeHolder type, ResolveType(ctx, batch.GetTypes()));
      if (type.type == NULLPTR) {
        return Status::TypeError("expected output datatype is null");
      }
      if (*type.type != *val->type()) {
        return Status::TypeError("Expected output datatype ", type.type->ToString(),
                                 ", but function returned datatype ",
                                 val->type()->ToString());
      }
      out->value = std::move(val->data());
      return Status::OK();
    } else {
      return Status::TypeError("Unexpected output type: ", Py_TYPE(result.obj())->tp_name,
                               " (expected Array)");
    }
    return Status::OK();
  }
};

Status PythonUdfExec(compute::KernelContext* ctx, const compute::ExecSpan& batch,
                     compute::ExecResult* out) {
  auto udf = static_cast<PythonUdf*>(ctx->kernel()->data.get());
  return SafeCallIntoPython([&]() -> Status { return udf->Exec(ctx, batch, out); });
}

template <class Function, class Kernel>
Status RegisterUdf(PyObject* function, compute::KernelInit kernel_init,
                   UdfWrapperCallback cb, const UdfOptions& options,
                   compute::FunctionRegistry* registry) {
  if (!PyCallable_Check(function)) {
    return Status::TypeError("Expected a callable Python object.");
  }
  auto scalar_func =
      std::make_shared<Function>(options.func_name, options.arity, options.func_doc);
  std::vector<compute::InputType> input_types;
  for (const auto& in_dtype : options.input_types) {
    input_types.emplace_back(in_dtype);
  }
  compute::OutputType output_type(options.output_type);
  // Take reference before wrapping with OwnedRefNoGIL
  Py_INCREF(function);
  auto udf_data = std::make_shared<PythonUdf>(
      std::make_shared<OwnedRefNoGIL>(function), cb,
      TypeHolder::FromTypes(options.input_types), options.output_type);
  Kernel kernel(
      compute::KernelSignature::Make(std::move(input_types), std::move(output_type),
                                     options.arity.is_varargs),
      PythonUdfExec, kernel_init);
  kernel.data = std::move(udf_data);

  kernel.mem_allocation = compute::MemAllocation::NO_PREALLOCATE;
  kernel.null_handling = compute::NullHandling::COMPUTED_NO_PREALLOCATE;
  RETURN_NOT_OK(scalar_func->AddKernel(std::move(kernel)));
  if (registry == NULLPTR) {
    registry = compute::GetFunctionRegistry();
  }
  RETURN_NOT_OK(registry->AddFunction(std::move(scalar_func)));
  return Status::OK();
}

}  // namespace

Status RegisterScalarFunction(PyObject* function, UdfWrapperCallback cb,
                              const UdfOptions& options,
                              compute::FunctionRegistry* registry) {
  return RegisterUdf<compute::ScalarFunction, compute::ScalarKernel>(
      function, PythonUdfKernelInit{std::make_shared<OwnedRefNoGIL>(function)}, cb,
      options, registry);
}

Status RegisterVectorFunction(PyObject* function, UdfWrapperCallback cb,
                              const UdfOptions& options,
                              compute::FunctionRegistry* registry) {
  return RegisterUdf<compute::VectorFunction, compute::VectorKernel>(
      function, PythonUdfKernelInit{std::make_shared<OwnedRefNoGIL>(function)}, cb,
      options, registry);
}

Status RegisterTabularFunction(PyObject* function, UdfWrapperCallback cb,
                               const UdfOptions& options,
                               compute::FunctionRegistry* registry) {
  if (options.arity.num_args != 0 || options.arity.is_varargs) {
    return Status::NotImplemented("tabular function of non-null arity");
  }
  if (options.output_type->id() != Type::type::STRUCT) {
    return Status::Invalid("tabular function with non-struct output");
  }
  return RegisterUdf<compute::ScalarFunction, compute::ScalarKernel>(
      function, PythonTableUdfKernelInit{std::make_shared<OwnedRefNoGIL>(function), cb},
      cb, options, registry);
}

Status RegisterScalarAggregateFunction(PyObject* function, UdfWrapperCallback cb,
                                       const UdfOptions& options,
                                       compute::FunctionRegistry* registry) {
  if (!PyCallable_Check(function)) {
    return Status::TypeError("Expected a callable Python object.");
  }

  if (registry == NULLPTR) {
    registry = compute::GetFunctionRegistry();
  }

  static auto default_scalar_aggregate_options =
      compute::ScalarAggregateOptions::Defaults();
  auto aggregate_func = std::make_shared<compute::ScalarAggregateFunction>(
      options.func_name, options.arity, options.func_doc,
      &default_scalar_aggregate_options);

  std::vector<compute::InputType> input_types;
  for (const auto& in_dtype : options.input_types) {
    input_types.emplace_back(in_dtype);
  }
  compute::OutputType output_type(options.output_type);

  // Take reference before wrapping with OwnedRefNoGIL
  Py_INCREF(function);
  auto function_ref = std::make_shared<OwnedRefNoGIL>(function);

  compute::KernelInit init = [cb, function_ref, options](
                                 compute::KernelContext* ctx,
                                 const compute::KernelInitArgs& args)
      -> Result<std::unique_ptr<compute::KernelState>> {
    return std::make_unique<PythonUdfScalarAggregatorImpl>(
        function_ref, cb, options.input_types, options.output_type);
  };

  auto sig = compute::KernelSignature::Make(
      std::move(input_types), std::move(output_type), options.arity.is_varargs);
  compute::ScalarAggregateKernel kernel(std::move(sig), std::move(init),
                                        AggregateUdfConsume, AggregateUdfMerge,
                                        AggregateUdfFinalize, /*ordered=*/false);
  RETURN_NOT_OK(aggregate_func->AddKernel(std::move(kernel)));
  RETURN_NOT_OK(registry->AddFunction(std::move(aggregate_func)));
  return Status::OK();
}

/// \brief Create a new UdfOptions with adjustment for hash kernel
/// \param options User provided udf options
UdfOptions AdjustForHashAggregate(const UdfOptions& options) {
  UdfOptions hash_options;
  // Append hash_ before the function name to separate from the scalar
  // version
  hash_options.func_name = "hash_" + options.func_name;
  // Extend input types with group id. Group id is appended by the group
  // aggregation node. Here we change both arity and input types
  if (options.arity.is_varargs) {
    hash_options.arity = options.arity;
  } else {
    hash_options.arity = compute::Arity(options.arity.num_args + 1, false);
  }
  // Changing the function doc shouldn't be necessarily because group id
  // is not user visible, however, this is currently needed to pass the
  // function validation. The name group_id_array is consistent with
  // hash kernels in hash_aggregate.cc
  hash_options.func_doc = options.func_doc;
  hash_options.func_doc.arg_names.emplace_back("group_id_array");
  std::vector<std::shared_ptr<DataType>> input_dtypes = options.input_types;
  input_dtypes.emplace_back(uint32());
  hash_options.input_types = std::move(input_dtypes);
  hash_options.output_type = options.output_type;
  return hash_options;
}

Status RegisterHashAggregateFunction(PyObject* function, UdfWrapperCallback cb,
                                     const UdfOptions& options,
                                     compute::FunctionRegistry* registry) {
  if (!PyCallable_Check(function)) {
    return Status::TypeError("Expected a callable Python object.");
  }

  if (registry == NULLPTR) {
    registry = compute::GetFunctionRegistry();
  }

  UdfOptions hash_options = AdjustForHashAggregate(options);

  std::vector<compute::InputType> input_types;
  for (const auto& in_dtype : hash_options.input_types) {
    input_types.emplace_back(in_dtype);
  }
  compute::OutputType output_type(hash_options.output_type);

  static auto default_hash_aggregate_options =
      compute::ScalarAggregateOptions::Defaults();
  auto hash_aggregate_func = std::make_shared<compute::HashAggregateFunction>(
      hash_options.func_name, hash_options.arity, hash_options.func_doc,
      &default_hash_aggregate_options);

  // Take reference before wrapping with OwnedRefNoGIL
  Py_INCREF(function);
  auto function_ref = std::make_shared<OwnedRefNoGIL>(function);
  compute::KernelInit init = [function_ref, cb, hash_options](
                                 compute::KernelContext* ctx,
                                 const compute::KernelInitArgs& args)
      -> Result<std::unique_ptr<compute::KernelState>> {
    return std::make_unique<PythonUdfHashAggregatorImpl>(
        function_ref, cb, hash_options.input_types, hash_options.output_type);
  };

  auto sig = compute::KernelSignature::Make(
      std::move(input_types), std::move(output_type), hash_options.arity.is_varargs);

  compute::HashAggregateKernel kernel(
      std::move(sig), std::move(init), HashAggregateUdfResize, HashAggregateUdfConsume,
      HashAggregateUdfMerge, HashAggregateUdfFinalize, /*ordered=*/false);
  RETURN_NOT_OK(hash_aggregate_func->AddKernel(std::move(kernel)));
  RETURN_NOT_OK(registry->AddFunction(std::move(hash_aggregate_func)));
  return Status::OK();
}

Status RegisterAggregateFunction(PyObject* function, UdfWrapperCallback cb,
                                 const UdfOptions& options,
                                 compute::FunctionRegistry* registry) {
  RETURN_NOT_OK(RegisterScalarAggregateFunction(function, cb, options, registry));
  RETURN_NOT_OK(RegisterHashAggregateFunction(function, cb, options, registry));

  return Status::OK();
}

Result<std::shared_ptr<RecordBatchReader>> CallTabularFunction(
    const std::string& func_name, const std::vector<Datum>& args,
    compute::FunctionRegistry* registry) {
  if (args.size() != 0) {
    return Status::NotImplemented("non-empty arguments to tabular function");
  }
  if (registry == NULLPTR) {
    registry = compute::GetFunctionRegistry();
  }
  ARROW_ASSIGN_OR_RAISE(auto func, registry->GetFunction(func_name));
  if (func->kind() != compute::Function::SCALAR) {
    return Status::Invalid("tabular function of non-scalar kind");
  }
  auto arity = func->arity();
  if (arity.num_args != 0 || arity.is_varargs) {
    return Status::NotImplemented("tabular function of non-null arity");
  }
  auto kernels =
      arrow::internal::checked_pointer_cast<compute::ScalarFunction>(func)->kernels();
  if (kernels.size() != 1) {
    return Status::NotImplemented("tabular function with non-single kernel");
  }
  const compute::ScalarKernel* kernel = kernels[0];
  auto out_type = kernel->signature->out_type();
  if (out_type.kind() != compute::OutputType::FIXED) {
    return Status::Invalid("tabular kernel of non-fixed kind");
  }
  auto datatype = out_type.type();
  if (datatype->id() != Type::type::STRUCT) {
    return Status::Invalid("tabular kernel with non-struct output");
  }
  auto struct_type = arrow::internal::checked_cast<StructType*>(datatype.get());
  auto schema = ::arrow::schema(struct_type->fields());
  std::vector<TypeHolder> in_types;
  ARROW_ASSIGN_OR_RAISE(auto func_exec,
                        GetFunctionExecutor(func_name, in_types, NULLPTR, registry));
  auto next_func = [schema, func_exec = std::move(
                                func_exec)]() -> Result<std::shared_ptr<RecordBatch>> {
    std::vector<Datum> args;
    // passed_length of -1 or 0 with args.size() of 0 leads to an empty ExecSpanIterator
    // in exec.cc and to never invoking the source function, so 1 is passed instead
    // TODO: GH-33612: Support batch size in user-defined tabular functions
    ARROW_ASSIGN_OR_RAISE(auto datum, func_exec->Execute(args, /*passed_length=*/1));
    if (!datum.is_array()) {
      return Status::Invalid("UDF result of non-array kind");
    }
    std::shared_ptr<Array> array = datum.make_array();
    if (array->length() == 0) {
      return IterationTraits<std::shared_ptr<RecordBatch>>::End();
    }
    ARROW_ASSIGN_OR_RAISE(auto batch, RecordBatch::FromStructArray(std::move(array)));
    if (!schema->Equals(batch->schema())) {
      return Status::Invalid("UDF result with shape not conforming to schema");
    }
    return std::move(batch);
  };
  return RecordBatchReader::MakeFromIterator(MakeFunctionIterator(std::move(next_func)),
                                             schema);
}

}  // namespace py
}  // namespace arrow
