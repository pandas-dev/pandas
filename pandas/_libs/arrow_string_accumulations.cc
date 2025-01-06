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

static auto CumSum(const struct ArrowArrayView *array_view,
                   struct ArrowArray *out, bool skipna) {
  bool seen_na = false;
  std::stringstream ss{};

  for (int64_t i = 0; i < array_view->length; i++) {
    const bool isna = ArrowArrayViewIsNull(array_view, i);
    if (!skipna && (seen_na || isna)) {
      seen_na = true;
      ArrowArrayAppendNull(out, 1);
    } else {
      if (!isna) {
        const auto std_sv = ArrowArrayViewGetStringUnsafe(array_view, i);
        ss << std::string_view{std_sv.data,
                               static_cast<size_t>(std_sv.size_bytes)};
      }
      const auto str = ss.str();
      const ArrowStringView asv{str.c_str(), static_cast<int64_t>(str.size())};
      NANOARROW_THROW_NOT_OK(ArrowArrayAppendString(out, asv));
    }
  }
}

template <typename T>
concept MinOrMaxOp =
    std::same_as<T, std::less<>> || std::same_as<T, std::greater<>>;

template <auto Op>
  requires MinOrMaxOp<decltype(Op)>
static auto CumMinOrMax(const struct ArrowArrayView *array_view,
                        struct ArrowArray *out, bool skipna) {
  bool seen_na = false;
  std::optional<std::string> current_str{};

  for (int64_t i = 0; i < array_view->length; i++) {
    const bool isna = ArrowArrayViewIsNull(array_view, i);
    if (!skipna && (seen_na || isna)) {
      seen_na = true;
      ArrowArrayAppendNull(out, 1);
    } else {
      if (!isna || current_str) {
        if (!isna) {
          const auto asv = ArrowArrayViewGetStringUnsafe(array_view, i);
          const nb::str pyval{asv.data, static_cast<size_t>(asv.size_bytes)};

          if (current_str) {
            const nb::str pycurrent{current_str->data(), current_str->size()};
            if (Op(pyval, pycurrent)) {
              current_str =
                  std::string{asv.data, static_cast<size_t>(asv.size_bytes)};
            }
          } else {
            current_str =
                std::string{asv.data, static_cast<size_t>(asv.size_bytes)};
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
    case NANOARROW_TYPE_STRING_VIEW:
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

    nanoarrow::UniqueArray chunk{};
    int errcode{};

    while ((errcode = ArrowArrayStreamGetNext(stream_.get(), chunk.get(),
                                              nullptr) == 0) &&
           chunk->release != nullptr) {
      struct ArrowArrayView array_view{};
      NANOARROW_THROW_NOT_OK(
          ArrowArrayViewInitFromSchema(&array_view, schema_.get(), nullptr));

      NANOARROW_THROW_NOT_OK(
          ArrowArrayViewSetArray(&array_view, chunk.get(), nullptr));

      if (accumulation_ == "cumsum") {
        CumSum(&array_view, uarray_out.get(), skipna_);
      } else if (accumulation_ == "cummin") {
        CumMinOrMax<std::less{}>(&array_view, uarray_out.get(), skipna_);
      } else if (accumulation_ == "cummax") {
        CumMinOrMax<std::greater{}>(&array_view, uarray_out.get(), skipna_);
      } else {
        throw std::runtime_error("Unexpected branch");
      }

      chunk.reset();
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
