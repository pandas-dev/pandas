#include <functional>
#include <optional>
#include <sstream>
#include <string_view>
#include <tuple>

#include <nanoarrow/nanoarrow.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>

using namespace nanoarrow::literals;
namespace nb = nanobind;

static auto ReleaseArrowArray(void *ptr) noexcept -> void {
  auto array = static_cast<struct ArrowArray *>(ptr);
  if (array->release != nullptr) {
    ArrowArrayRelease(array);
  }

  delete array;
}

static auto ReleaseArrowSchema(void *ptr) noexcept -> void {
  auto schema = static_cast<struct ArrowSchema *>(ptr);
  if (schema->release != nullptr) {
    ArrowSchemaRelease(schema);
  }

  delete schema;
}

template <size_t OffsetSize>
static auto CumSum(struct ArrowArrayStream *array_stream,
                   struct ArrowArray *out, bool skipna) {
  bool seen_na = false;
  std::stringstream ss{};

  nanoarrow::UniqueSchema schema{};
  NANOARROW_THROW_NOT_OK(
      ArrowArrayStreamGetSchema(array_stream, schema.get(), nullptr));

  nanoarrow::ViewArrayStream array_stream_view(array_stream);
  for (const auto &array : array_stream_view) {
    for (const auto &sv : nanoarrow::ViewArrayAsBytes<OffsetSize>(&array)) {
      if ((!sv || seen_na) && !skipna) {
        seen_na = true;
        ArrowArrayAppendNull(out, 1);
      } else {
        if (sv) {
          ss << std::string_view{(*sv).data,
                                 static_cast<size_t>((*sv).size_bytes)};
        }
        const auto str = ss.str();
        const ArrowStringView asv{str.c_str(),
                                  static_cast<int64_t>(str.size())};
        NANOARROW_THROW_NOT_OK(ArrowArrayAppendString(out, asv));
      }
    }
  }
}

// TODO: doesn't seem like all compilers in CI support this?
// template <typename T>
// concept MinOrMaxOp =
//     std::same_as<T, std::less<>> || std::same_as<T, std::greater<>>;
// template <size_t OffsetSize, auto Op>
//  requires MinOrMaxOp<decltype(Op)>
template <size_t OffsetSize, template <typename...> typename CompareOp>
static auto CumMinOrMax(struct ArrowArrayStream *array_stream,
                        struct ArrowArray *out, bool skipna) {
  bool seen_na = false;
  std::optional<std::string> current_str{};

  nanoarrow::UniqueSchema schema{};
  NANOARROW_THROW_NOT_OK(
      ArrowArrayStreamGetSchema(array_stream, schema.get(), nullptr));

  nanoarrow::ViewArrayStream array_stream_view(array_stream);
  for (const auto &array : array_stream_view) {
    for (const auto &sv : nanoarrow::ViewArrayAsBytes<OffsetSize>(&array)) {
      if ((!sv || seen_na) && !skipna) {
        seen_na = true;
        ArrowArrayAppendNull(out, 1);
      } else {
        if (sv || current_str) {
          if (sv) {
            const nb::str pyval{(*sv).data,
                                static_cast<size_t>((*sv).size_bytes)};
            if (current_str) {
              const nb::str pycurrent{current_str->data(), current_str->size()};
              if (CompareOp<const nb::str &>{}(pyval, pycurrent)) {
                current_str = std::string{
                    (*sv).data, static_cast<size_t>((*sv).size_bytes)};
              }
            } else {
              current_str = std::string{(*sv).data,
                                        static_cast<size_t>((*sv).size_bytes)};
            }
          }

          struct ArrowStringView out_sv{
              current_str->data(), static_cast<int64_t>(current_str->size())};
          NANOARROW_THROW_NOT_OK(ArrowArrayAppendString(out, out_sv));
        } else {
          ArrowArrayAppendEmpty(out, 1);
        }
      }
    }
  }
}

class ArrowStringAccumulation {
public:
  ArrowStringAccumulation(nb::object array_object, std::string accumulation,
                          bool skipna)
      : skipna_(skipna) {
    if ((accumulation == "cumsum") || (accumulation == "cummin") ||
        (accumulation == "cummax")) {
      accumulation_ = std::move(accumulation);
    } else {
      const auto error_message =
          std::string("Unsupported accumulation: ") + accumulation;
      throw nb::value_error(error_message.c_str());
    }

    const auto obj = nb::getattr(array_object, "__arrow_c_stream__")();
    const auto pycapsule_obj = nb::cast<nb::capsule>(obj);

    const auto stream = static_cast<struct ArrowArrayStream *>(
        PyCapsule_GetPointer(pycapsule_obj.ptr(), "arrow_array_stream"));
    if (stream == nullptr) {
      throw std::invalid_argument("Invalid Arrow Stream capsule provided!");
    }

    if (stream->get_schema(stream, schema_.get()) != 0) {
      std::string error_msg{stream->get_last_error(stream)};
      throw std::runtime_error("Could not read from arrow schema:" + error_msg);
    }
    struct ArrowSchemaView schema_view{};
    NANOARROW_THROW_NOT_OK(
        ArrowSchemaViewInit(&schema_view, schema_.get(), nullptr));

    switch (schema_view.type) {
    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_LARGE_STRING:
      break;
    default:
      const auto error_message =
          std::string("Expected a string-like array type, got: ") +
          ArrowTypeString(schema_view.type);
      throw std::invalid_argument(error_message);
    }

    ArrowArrayStreamMove(stream, stream_.get());
  }

  std::pair<nb::capsule, nb::capsule> Accumulate(nb::object) {
    struct ArrowSchemaView schema_view{};
    NANOARROW_THROW_NOT_OK(
        ArrowSchemaViewInit(&schema_view, schema_.get(), nullptr));
    auto uschema = nanoarrow::UniqueSchema{};
    ArrowSchemaInit(uschema.get());
    NANOARROW_THROW_NOT_OK(ArrowSchemaSetType(uschema.get(), schema_view.type));

    // TODO: even though we are reading a stream we are returning an array
    // We should return a like sized stream of data in the future
    auto uarray_out = nanoarrow::UniqueArray{};
    NANOARROW_THROW_NOT_OK(
        ArrowArrayInitFromSchema(uarray_out.get(), uschema.get(), nullptr));

    NANOARROW_THROW_NOT_OK(ArrowArrayStartAppending(uarray_out.get()));

    if (accumulation_ == "cumsum") {
      if (schema_view.type == NANOARROW_TYPE_STRING) {
        CumSum<32>(stream_.get(), uarray_out.get(), skipna_);
      } else {
        CumSum<64>(stream_.get(), uarray_out.get(), skipna_);
      }

    } else if (accumulation_ == "cummin") {
      if (schema_view.type == NANOARROW_TYPE_STRING) {
        CumMinOrMax<32, std::less>(stream_.get(), uarray_out.get(), skipna_);
      } else {
        CumMinOrMax<64, std::less>(stream_.get(), uarray_out.get(), skipna_);
      }
    } else if (accumulation_ == "cummax") {
      if (schema_view.type == NANOARROW_TYPE_STRING) {
        CumMinOrMax<32, std::greater>(stream_.get(), uarray_out.get(), skipna_);
      } else {
        CumMinOrMax<64, std::greater>(stream_.get(), uarray_out.get(), skipna_);
      }
    } else {
      throw std::runtime_error("Unexpected branch");
    }

    NANOARROW_THROW_NOT_OK(
        ArrowArrayFinishBuildingDefault(uarray_out.get(), nullptr));

    auto out_schema = new struct ArrowSchema;
    ArrowSchemaMove(uschema.get(), out_schema);
    nb::capsule schema_capsule{out_schema, "arrow_schema", &ReleaseArrowSchema};

    auto out_array = new struct ArrowArray;
    ArrowArrayMove(uarray_out.get(), out_array);
    nb::capsule array_capsule{out_array, "arrow_array", &ReleaseArrowArray};

    return std::pair<nb::capsule, nb::capsule>{schema_capsule, array_capsule};
  }

private:
  nanoarrow::UniqueArrayStream stream_;
  nanoarrow::UniqueSchema schema_;
  std::string accumulation_;
  bool skipna_;
};

NB_MODULE(arrow_string_accumulations, m) {
  nb::class_<ArrowStringAccumulation>(m, "ArrowStringAccumulation")
      .def(nb::init<nb::object, std::string, bool>())
      .def("__arrow_c_array__", &ArrowStringAccumulation::Accumulate,
           nb::arg("requested_schema") = nb::none());
}
