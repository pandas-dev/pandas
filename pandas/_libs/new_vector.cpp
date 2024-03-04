#include <cstring>
#include <stdint.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <string>
#include <vector>

#include "pandas/vendored/klib/cpp/khashl.hpp"

namespace nb = nanobind;

template <typename T> class PandasVector {
public:
  size_t __len__() const { return vec_.size(); }

  void Append(const T &value) { vec_.emplace_back(value); }

  void Extend(const nb::ndarray<const T, nb::ndim<1>> &values) {
    const size_t current_size = vec_.size();
    const size_t n_elems = values.shape(0);
    vec_.resize(vec_.size() + n_elems);
    memcpy(vec_.data() + current_size, values.data(), sizeof(T) * n_elems);
  }

  auto ToNdArray() const -> nb::ndarray<nb::numpy, const T, nb::ndim<1>> {
    // don't manually delete anything - RAII will take care of it
    // nb::capsule owner(*this, [](void *) noexcept {});
    return nb::ndarray<nb::numpy, const T, nb::ndim<1>>(vec_.data(),
                                                        {vec_.size()});
  }

private:
  std::vector<T> vec_;
};

template <typename T> class PandasHashTable {
public:

  explicit PandasHashTable ([[maybe_unused]] int64_t size_hint=1, bool uses_mask=false) : uses_mask_(uses_mask) { }

  auto __len__() const noexcept -> size_t{
    return 42;
  }

  auto __contains__(nb::object key) const noexcept -> bool {
    if (uses_mask_ and key.is_none()) {
      return -1 != na_position_;
    }

    auto ckey = nb::cast<T>(key);
    return hash_map_.get(ckey) != hash_map_.end();
  }

  ///
  /// Extracts the position of val from the hashtable.
  ///
  /// Parameters
  /// ----------
  /// val : Scalar
  ///     The value that is looked up in the hashtable
  ///
  /// Returns
  /// -------
  /// The position of the requested integer.
  auto GetItem(T key) -> int64_t {
    if (hash_map_.get(key) != hash_map_.end()) {
      return hash_map_.value(key);
    } else {
      throw nb::key_error(std::to_string(key).c_str());
    }
  }

  ///
  /// Extracts the position of na_value from the hashtable.
  ///
  /// Returns
  /// -------
  /// The position of the last na value.
  auto GetNA() const -> int64_t {
    // TODO: missing NotImplementedError for mask, although this should really
    // just be templated out
    if (na_position_ == -1) {
      throw nb::key_error("NA");
    }

    return na_position_;
  }

  auto SetItem(T key, ssize_t val) -> void {
    int absent;
    auto k = hash_map_.put(key, &absent);
    if (!absent) {
      hash_map_.value(k) = val;
    } else {
      throw nb::key_error(std::to_string(key).c_str());
    }
  }

  auto SetNA(ssize_t val) noexcept -> void {
    // TODO: missing NotImplementedError for mask, although this should really
    // just be templated out
    na_position_ = val;
  }

  auto MapKeysToValues(const nb::ndarray<const T, nb::ndim<1>> &keys,
                       const nb::ndarray<const T, nb::ndim<1>> &values) noexcept -> void {
    const auto keys_v = keys.view();
    const auto values_v = values.view();

    for (size_t i = 0; i < values_v.shape(0); i++) {
      hash_map_[keys_v(i)] = values_v(i);
    }
  }

  auto MapLocations(const nb::ndarray<const T, nb::ndim<1>> &values, nb::object mask_obj) -> void {
    // TODO: templating!
    if (uses_mask_ && mask_obj.is_none()) {
      throw std::invalid_argument("mask must not be None!");
    }

    // TODO: how can we release the GIL imperatively?
    // https://nanobind.readthedocs.io/en/latest/api_core.html#gil-management
    const auto values_v = values.view();
    if (uses_mask_) {
      const auto mask = nb::cast<nb::ndarray<const uint8_t, nb::ndim<1>>>(mask_obj);
      const auto mask_v = mask.view();
      size_t na_position = na_position_;  // pandas uses int8_t here - why?
      for (size_t i = 0; i < values_v.shape(0); i++) {
        if (mask_v(i)) {
          na_position = i;
        } else {
          hash_map_[values(i)] = i;
        }
      }
      na_position_ = na_position;
    } else {
      for (size_t i = 0; i < values_v.shape(0); i++) {
        hash_map_[values(i)] = i;
      }
    }
  }

  auto Lookup(const nb::ndarray<const T, nb::ndim<1>> &values, nb::object mask_obj) ->
    nb::ndarray<nb::numpy, const size_t, nb::ndim<1>> {
    // TODO: templating!
    if (uses_mask_ && mask_obj.is_none()) {
      throw std::invalid_argument("mask must not be None!");
    }

    const size_t n = values.shape(0);
    size_t *locs = new size_t[n];
    const auto values_v = values.view();

    // TODO: how can we release the GIL imperatively?
    // https://nanobind.readthedocs.io/en/latest/api_core.html#gil-management
    if (uses_mask_) {
      const auto mask = nb::cast<nb::ndarray<const uint8_t, nb::ndim<1>>>(mask_obj);
      const auto mask_v = mask.view();
      for (size_t i = 0; i < n; i++) {
        if (mask_v(i)) {
          locs[i] = na_position_;
        } else {
          const auto val = values_v(i);
          const auto position = hash_map_.get(val);
          if (position == hash_map_.end()) {
            locs[i] = -1;
          } else {
            locs[i] = hash_map_.value(val);
          }
        }
      }
    } else {
      for (size_t i = 0; i < n; i++) {
        const auto val = values_v(i);
        const auto position = hash_map_.get(val);
        if (position == hash_map_.end()) {
          locs[i] = -1;
        } else {
          locs[i] = hash_map_.value(val);
        }
      }
    }

    nb::capsule owner(locs, [](void *p) noexcept {
      delete[] (size_t *)p;
    });

    return nb::ndarray<nb::numpy, const size_t, nb::ndim<1>>(locs, {1,}, owner);
  }

  /*
  using UniqueResultT = std::tuple<nb::ndarray<nb::numpy, const T, nb::ndim<1>>,
    nb::ndarray<nb::numpy, const size_t, nb::ndim<1>>,
    nb::ndarray<nb::numpy, const bool, nb::ndim<1>>
  >;
  auto Unique(const nb::ndarray<const T, nb::ndim<1>> &values, PandasVector<T>& uniques,
              ssize_t count_prior=0, ssize_t na_sentinel=-1, nb::object na_value= nb::none(),
              bool ignore_na=false, nb::object mask_obj=nb::none(), bool return_inverse=false,
              bool use_result_mask=false) -> UniqueResultT {
    // TODO: templating!
    if (uses_mask_ && mask_obj.is_none()) {
      throw std::invalid_argument("mask must not be None!");
    }

    auto result_mask = PandasVector<uint8_t>();

  }
  */

private:
  klib::KHashMap<T, int64_t, std::hash<T>> hash_map_;
  bool uses_mask_;
  ssize_t na_position_ = -1;
};

#define BIND_VECTOR(TYPE, NAME) do { \
  nb::class_<PandasVector<TYPE>>(m, NAME) \
      .def(nb::init<>()) \
      .def("__len__", &PandasVector<TYPE>::__len__) \
      .def("append", &PandasVector<TYPE>::Append) \
      .def("extend", &PandasVector<TYPE>::Extend) \
      .def("to_array", &PandasVector<TYPE>::ToNdArray); \
  } while(0)

NB_MODULE(new_vector, m) {
  BIND_VECTOR(int8_t, "Int8Vector");
  BIND_VECTOR(int16_t, "Int16Vector");
  BIND_VECTOR(int32_t, "Int32Vector");
  BIND_VECTOR(int64_t, "Int64Vector");
  BIND_VECTOR(uint8_t, "UInt8Vector");
  BIND_VECTOR(uint16_t, "UInt16Vector");
  BIND_VECTOR(uint32_t, "UInt32Vector");
  BIND_VECTOR(uint64_t, "UInt64Vector");
  BIND_VECTOR(float, "Float32Vector");
  BIND_VECTOR(double, "Float64Vector");

  nb::class_<PandasHashTable<int64_t>>(m, "Int64HashTable")
    .def(nb::init<int64_t, bool>())
    .def("__len__", &PandasHashTable<int64_t>::__len__)
    .def("__contains__", &PandasHashTable<int64_t>::__contains__)
    .def("get_item", &PandasHashTable<int64_t>::GetItem)
    .def("set_item", &PandasHashTable<int64_t>::SetItem)
    .def("set_na", &PandasHashTable<int64_t>::SetNA)
    .def("map_keys_to_values", &PandasHashTable<int64_t>::MapKeysToValues)
    .def("map_locations", &PandasHashTable<int64_t>::MapLocations)
    .def("lookup", &PandasHashTable<int64_t>::Lookup);
}
