#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cmath>
#include <cstring>
#include <map>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/tuple.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "pandas/vendored/klib/cpp/khashl.hpp"

namespace nb = nanobind;

///
/// we typically defer to the standard library, but can create
/// specializations here to override that behavior (ex: NaN comparison)
/// and support arbitrary types
///
template <typename T> struct PandasHashFunction {
  auto operator()(const T &value) const { return std::hash<T>()(value); }
};

template <>
auto PandasHashFunction<float>::operator()(const float &value) const {
  if (std::isnan(value)) {
    return static_cast<decltype(std::hash<float>()(value))>(0);
  }

  return std::hash<float>()(value);
}

template <>
auto PandasHashFunction<double>::operator()(const double &value) const {
  if (std::isnan(value)) {
    return static_cast<decltype(std::hash<double>()(value))>(0);
  }

  return std::hash<double>()(value);
}

template <typename T> struct PandasHashEquality {
  size_t operator()(const T &lhs, const T &rhs) const { return lhs == rhs; }
};

template <>
size_t PandasHashEquality<float>::operator()(const float &lhs,
                                             const float &rhs) const {
  if (std::isnan(lhs) && std::isnan(rhs)) {
    return true;
  }

  return lhs == rhs;
}

template <>
size_t PandasHashEquality<double>::operator()(const double &lhs,
                                              const double &rhs) const {
  if (std::isnan(lhs) && std::isnan(rhs)) {
    return true;
  }

  return lhs == rhs;
}

template <typename T> class PandasVector {
public:
  explicit PandasVector<T>() : external_view_exists_(false) {}
  ~PandasVector<T>() = default;
  PandasVector<T>(PandasVector<T> const &) = delete;
  void operator=(PandasVector<T> const &) = delete;
  PandasVector<T>(PandasVector<T> &&) = default;
  PandasVector<T> &operator=(PandasVector<T> &&) = default;

  auto __len__() const { return vec_.size(); }

  auto Append(const T &value) -> void {
    if (external_view_exists_) {
      throw std::domain_error("external reference but Vector.resize() needed");
    }
    vec_.emplace_back(value);
  }

  auto Extend(const nb::ndarray<const T, nb::ndim<1>> &values) -> void {
    if (external_view_exists_) {
      throw std::domain_error("external reference but Vector.resize() needed");
    }

    const size_t current_size = vec_.size();
    const size_t n_elems = values.shape(0);
    vec_.resize(vec_.size() + n_elems);
    memcpy(vec_.data() + current_size, values.data(), sizeof(T) * n_elems);
  }

  auto ToNdArray() -> nb::object {
    const auto ndarray = nb::ndarray<nb::numpy, const T, nb::ndim<1>>(
        vec_.data(), {vec_.size()});

    external_view_exists_ = true;

    return nb::cast(ndarray, nb::rv_policy::reference_internal);
  }

private:
  std::vector<T> vec_;
  bool external_view_exists_;
};

template <typename T, bool IsMasked> class PandasHashTable {
public:
  // in English, if the return value from the hashing function is 4 bytes or less, use
  // uint32_t for the khash "int" size. Otherwise use 64 bits
  using HashValueT = typename std::conditional<sizeof(decltype(PandasHashFunction<T>()(T()))) <= 4, uint32_t, uint64_t>::type;
  explicit PandasHashTable<T, IsMasked>() = default;
  explicit PandasHashTable<T, IsMasked>(size_t new_size) {
    // TODO: C++20 std::in_range would be great to safely check cast
    hash_map_.resize(static_cast<HashValueT>(new_size));
  }

  auto __len__() const noexcept { return hash_map_.size(); }

  auto __contains__(nb::object key) const noexcept {
    if constexpr (IsMasked) {
      if (key.is_none()) {
        return -1 != na_position_;
      }
    }

    auto ckey = nb::cast<T>(key);
    return hash_map_.get(ckey) != hash_map_.end();
  }

  auto SizeOf() const noexcept {
    constexpr auto overhead = 4 * sizeof(uint32_t) + 3 * sizeof(uint32_t *);
    const auto for_flags =
      std::max(1UL, hash_map_.n_buckets() >> 5) * sizeof(uint32_t);
    const auto for_pairs =
        hash_map_.n_buckets() * (sizeof(T) + sizeof(Py_ssize_t));

    return overhead + for_flags + for_pairs;
  }

  auto GetState() const noexcept -> std::map<const char *, size_t> {
    return {{"n_buckets", hash_map_.n_buckets()}, {"size", hash_map_.size()}};
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
  auto GetItem(T key) {
    const auto k = hash_map_.get(key);
    if (k != hash_map_.end()) {
      return hash_map_.value(k);
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
  auto GetNA() const {
    // TODO: missing NotImplementedError for mask, although this should really
    // just be templated out
    if (na_position_ == -1) {
      throw nb::key_error("NA");
    }

    return na_position_;
  }

  auto SetItem(T key, Py_ssize_t val) noexcept -> void { hash_map_[key] = val; }

  auto SetNA([[maybe_unused]] Py_ssize_t val) noexcept -> void {
    if constexpr (IsMasked) {
      na_position_ = val;
    }
  }

  auto MapKeysToValues(
      const nb::ndarray<const T, nb::ndim<1>> &keys,
      const nb::ndarray<const int64_t, nb::ndim<1>> &values) noexcept -> void {
    nb::call_guard<nb::gil_scoped_release>();
    const auto keys_v = keys.view();
    const auto values_v = values.view();
    const auto n = values_v.shape(0);
    for (auto i = decltype(n){0}; i < n; i++) {
      hash_map_[keys_v(i)] = values_v(i);
    }
  }

  auto MapLocations(const nb::ndarray<const T, nb::ndim<1>> &values,
                    nb::object mask) -> void {
    if constexpr (IsMasked) {
      if (mask.is_none()) {
        throw std::invalid_argument("mask must not be None!");
      }
    }

    nb::call_guard<nb::gil_scoped_release>();
    const auto values_v = values.view();
    const auto n = values_v.shape(0);
    if constexpr (IsMasked) {
      const auto mask_base =
          nb::cast<nb::ndarray<const uint8_t, nb::ndim<1>>>(mask);
      const auto mask_v = mask_base.view();
      auto na_position = na_position_; // pandas uses int8_t here - why?
      for (auto i = decltype(n){0}; i < n; i++) {
        if (mask_v(i)) {
          na_position = i;
        } else {
          hash_map_[values_v(i)] = i;
        }
      }
      na_position_ = na_position;
    } else {
      for (auto i = decltype(n){0}; i < n; i++) {
        const auto key = values_v(i);
        hash_map_[key] = i;
      }
    }
  }

  auto Lookup(const nb::ndarray<const T, nb::ndim<1>> &values, nb::object mask)
      -> nb::ndarray<nb::numpy, const Py_ssize_t, nb::ndim<1>> {
    if constexpr (IsMasked) {
      if (mask.is_none()) {
        throw std::invalid_argument("mask must not be None!");
      }
    }

    nb::call_guard<nb::gil_scoped_release>();
    const auto n = values.shape(0);
    auto *locs = new Py_ssize_t[n];
    const auto values_v = values.view();

    // TODO: how can we release the GIL imperatively?
    // https://nanobind.readthedocs.io/en/latest/api_core.html#gil-management
    if constexpr (IsMasked) {
      const auto mask_base =
          nb::cast<nb::ndarray<const uint8_t, nb::ndim<1>>>(mask);
      const auto mask_v = mask_base.view();
      for (auto i = decltype(n){0}; i < n; i++) {
        if (mask_v(i)) {
          locs[i] = na_position_;
        } else {
          const auto val = values_v(i);
          const auto position = hash_map_.get(val);
          if (position == hash_map_.end()) {
            locs[i] = -1;
          } else {
            locs[i] = hash_map_.value(position);
          }
        }
      }
    } else {
      for (auto i = decltype(n){0}; i < n; i++) {
        const auto val = values_v(i);
        const auto position = hash_map_.get(val);
        if (position == hash_map_.end()) {
          locs[i] = -1;
        } else {
          locs[i] = hash_map_.value(position);
        }
      }
    }

    nb::gil_scoped_acquire();
    nb::capsule owner(locs, [](void *p) noexcept { delete[](size_t *) p; });

    const size_t shape[1] = {n};
    return nb::ndarray<nb::numpy, const Py_ssize_t, nb::ndim<1>>(locs, 1, shape,
                                                                 owner);
  }

  /// Mutates uniques argument
  auto Unique(const nb::ndarray<const T, nb::ndim<1>> &values,
              bool return_inverse = false, nb::object mask = nb::none())
      -> nb::object {
    PandasVector<T> uniques;
    const bool use_result_mask = mask.is_none() ? false : true;
    return UniqueInternal(values, uniques, 0, -1, nb::none(), false, mask,
                          return_inverse, use_result_mask);
  }

  auto Factorize(const nb::ndarray<const T, nb::ndim<1>> &values,
                 Py_ssize_t na_sentinel = -1, nb::object na_value = nb::none(),
                 nb::object mask = nb::none(), bool ignore_na = true)
      -> nb::object {
    PandasVector<T> uniques;
    return UniqueInternal(values, uniques, 0, na_sentinel, na_value, ignore_na,
                          mask, true, false);
  }

  auto GetLabels(const nb::ndarray<const T, nb::ndim<1>> &values,
                 PandasVector<T> &uniques, Py_ssize_t count_prior = 0,
                 Py_ssize_t na_sentinel = -1, nb::object na_value = nb::none(),
                 nb::object mask = nb::none(),
                 [[maybe_unused]] bool ignore_na = true) -> nb::object {
    const auto tup = UniqueInternal(values, uniques, count_prior, na_sentinel,
                                    na_value, true, mask, true, false);

    return tup[1];
  }

  auto GetLabelsGroupby(const nb::ndarray<const T, nb::ndim<1>> &values)
      -> std::tuple<nb::ndarray<nb::numpy, const Py_ssize_t, nb::ndim<1>>,
                    nb::object> {

    nb::call_guard<nb::gil_scoped_release>();
    const auto values_v = values.view();
    const auto n = values.shape(0);
    PandasVector<T> uniques;
    auto *labels = new Py_ssize_t[n];
    Py_ssize_t count = 0;

    for (auto i = decltype(n){0}; i < n; i++) {
      const auto val = values_v(i);

      // specific for groupby
      if (val < 0) {
        labels[i] = -1;
        continue;
      }

      auto k = hash_map_.get(val);
      if (k != hash_map_.end()) {
        labels[i] = hash_map_.value(k);
      } else {
        int dummy;
        k = hash_map_.put(val, &dummy);
        hash_map_.value(k) = count;
        uniques.Append(val);
        labels[i] = count;
        count++;
      }
    }

    nb::gil_scoped_acquire();

    const size_t shape[1] = {n};
    nb::capsule owner(labels, [](void *p) noexcept { delete[](size_t *) p; });
    const auto labels_arr =
        nb::ndarray<nb::numpy, const Py_ssize_t, nb::ndim<1>>(labels, 1, shape,
                                                              owner);
    return std::make_tuple(labels_arr, uniques.ToNdArray());
  }

private:
  auto UniqueInternal(const nb::ndarray<const T, nb::ndim<1>> &values,
                      PandasVector<T> &uniques, Py_ssize_t count_prior = 0,
                      [[maybe_unused]] Py_ssize_t na_sentinel = -1,
                      [[maybe_unused]] nb::object na_value = nb::none(),
                      bool ignore_na = false, nb::object mask_obj = nb::none(),
                      bool return_inverse = false, bool use_result_mask = false)
      -> nb::object {
    if (use_result_mask && return_inverse) {
      throw std::invalid_argument(
          "cannot supply both use_result_mask and return_inverse");
    }

    const bool use_na_value = !na_value.is_none();
    const auto na_val = use_na_value ? nb::cast<T>(na_value) : T();

    if (return_inverse) { // this is really factorize
      nb::ndarray<nb::numpy, const Py_ssize_t, nb::ndim<1>> labels;
      if (ignore_na) {
        if (use_na_value) {
          labels = UniqueWithInverse<true, true>(values, uniques, count_prior,
                                                 na_sentinel, na_val, mask_obj);
        } else {
          labels = UniqueWithInverse<true, false>(
              values, uniques, count_prior, na_sentinel, na_val, mask_obj);
        }
      } else {
        if (use_na_value) {
          labels = UniqueWithInverse<false, true>(
              values, uniques, count_prior, na_sentinel, na_val, mask_obj);
        } else {
          labels = UniqueWithInverse<false, false>(
              values, uniques, count_prior, na_sentinel, na_val, mask_obj);
        }
      }
      return nb::make_tuple(uniques.ToNdArray(), labels);
    } else if (use_result_mask) {
      PandasVector<uint8_t> mask;
      if (ignore_na) {
        if (use_na_value) {
          mask = UniqueWithResultMask<true, true>(values, uniques, mask_obj);
        } else {
          mask = UniqueWithResultMask<true, false>(values, uniques, mask_obj);
        }
      } else {
        if (use_na_value) {
          mask = UniqueWithResultMask<false, true>(values, uniques, mask_obj);
        } else {
          mask = UniqueWithResultMask<false, false>(values, uniques, mask_obj);
        }
      }

      return nb::make_tuple(uniques.ToNdArray(), mask.ToNdArray());
    } else {
      if (ignore_na) {
        if (use_na_value) {
          UniquesOnly<true, true>(values, uniques, mask_obj);
        } else {
          UniquesOnly<true, false>(values, uniques, mask_obj);
        }
      } else {
        if (use_na_value) {
          UniquesOnly<false, true>(values, uniques, mask_obj);
        } else {
          UniquesOnly<false, false>(values, uniques, mask_obj);
        }

        const auto out_array = uniques.ToNdArray();
        return nb::cast(out_array);
      }
    }

    throw std::runtime_error("should not reach this");
  }

  template <bool IgnoreNA, bool UseNAValue>
  auto UniqueWithInverse(const nb::ndarray<const T, nb::ndim<1>> &values,
                         PandasVector<T> &uniques, Py_ssize_t count_prior,
                         [[maybe_unused]] Py_ssize_t na_sentinel,
                         [[maybe_unused]] T na_value,
                         [[maybe_unused]] nb::object mask_obj = nb::none())
      -> nb::ndarray<nb::numpy, const Py_ssize_t, nb::ndim<1>> {
    if constexpr (IsMasked) {
      if (mask_obj.is_none()) {
        throw std::invalid_argument("mask must not be None!");
      }
    }

    const auto values_v = values.view();
    const auto n = values.shape(0);
    auto *labels = new Py_ssize_t[n];

    if constexpr (IsMasked) {
      using MaskT = nb::ndarray<const uint8_t, nb::ndim<1>>;
      MaskT mask;
      if (!nb::try_cast<MaskT>(mask_obj, mask, false)) {
        delete[] labels;
        throw std::invalid_argument("Could not convert mask to uint8_t array!");
      }
      nb::call_guard<nb::gil_scoped_release>();
      const auto mask_v = mask.view();

      for (auto i = decltype(n){0}; i < n; i++) {
        if constexpr (IgnoreNA) {
          if (mask_v(i)) {
            labels[i] = na_sentinel;
            continue;
          }
        }

        const auto val = values_v(i);
        auto k = hash_map_.get(val);
        if (k == hash_map_.end()) {
          int dummy;
          k = hash_map_.put(val, &dummy);
          uniques.Append(val);
          hash_map_.value(k) = count_prior;
          labels[i] = count_prior;
          count_prior++;
        } else {
          // k falls into a previous bucket
          // only relevant in case we need to construct the inverse
          labels[i] = hash_map_.value(k);
        }
      }
    } else {
      nb::call_guard<nb::gil_scoped_release>();
      for (auto i = decltype(n){0}; i < n; i++) {
        const auto val = values_v(i);

        if constexpr (IgnoreNA) {
          if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
            if (std::isnan(val)) {
              labels[i] = na_sentinel;
              continue;
            }
          } else if constexpr (UseNAValue) {
            if (val == na_value) {
              labels[i] = na_sentinel;
              continue;
            }
          }
        }

        auto k = hash_map_.get(val);
        if (k == hash_map_.end()) {
          int dummy;
          k = hash_map_.put(val, &dummy);
          uniques.Append(val);
          hash_map_.value(k) = count_prior;
          labels[i] = count_prior;
          count_prior++;
        } else {
          // k falls into a previous bucket
          // only relevant in case we need to construct the inverse
          labels[i] = hash_map_.value(k);
        }
      }
    }

    const size_t shape[1] = {n};
    nb::capsule owner(labels, [](void *p) noexcept { delete[](size_t *) p; });
    return nb::ndarray<nb::numpy, const Py_ssize_t, nb::ndim<1>>(labels, 1,
                                                                 shape, owner);
  }

  template <bool IgnoreNA, bool UseNAValue>
  auto UniqueWithResultMask(const nb::ndarray<const T, nb::ndim<1>> &values,
                            PandasVector<T> &uniques,
                            [[maybe_unused]] nb::object mask_obj = nb::none())
      -> PandasVector<uint8_t> {
    if constexpr (IsMasked) {
      if (mask_obj.is_none()) {
        throw std::invalid_argument("mask must not be None!");
      }
    }

    const auto values_v = values.view();
    const auto n = values.shape(0);
    PandasVector<uint8_t> result;

    if constexpr (IsMasked) {
      using MaskT = nb::ndarray<const uint8_t, nb::ndim<1>>;
      MaskT mask;
      if (!nb::try_cast<MaskT>(mask_obj, mask, false)) {
        throw std::invalid_argument("Could not convert mask to uint8_t array!");
      }
      nb::call_guard<nb::gil_scoped_release>();
      const auto mask_v = mask.view();

      bool seen_na = false;
      for (auto i = decltype(n){0}; i < n; i++) {
        const auto val = values_v(i);

        if constexpr (IgnoreNA) {
          // TODO: current pandas code is a bit messy here...
        } else {
          if (mask_v(i)) {
            if (seen_na) {
              continue;
            }

            seen_na = true;
            uniques.Append(val);
            result.Append(1);
            continue;
          }

          // NaN / pd.NA are treated the same? hmmm
          if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
            if (std::isnan(val)) {
              if (seen_na) {
                continue;
              }

              seen_na = true;
              uniques.Append(val);
              result.Append(1);
              continue;
            }
          }
        }

        auto k = hash_map_.get(val);
        if (k == hash_map_.end()) {
          int dummy;
          k = hash_map_.put(val, &dummy);
          uniques.Append(val);
          result.Append(0);
        }
      }
    } else {
      nb::call_guard<nb::gil_scoped_release>();
      for (auto i = decltype(n){0}; i < n; i++) {
        const auto val = values_v(i);
        auto k = hash_map_.get(val);
        if (k == hash_map_.end()) {
          int dummy;
          k = hash_map_.put(val, &dummy);
          uniques.Append(val);
          result.Append(0);
        }
      }
    }

    return result;
  }

  template <bool IgnoreNA, bool UseNAValue>
  auto UniquesOnly(const nb::ndarray<const T, nb::ndim<1>> &values,
                   PandasVector<T> &uniques,
                   [[maybe_unused]] nb::object mask_obj = nb::none()) -> void {
    if constexpr (IsMasked) {
      if (mask_obj.is_none()) {
        throw std::invalid_argument("mask must not be None!");
      }
    }

    const auto values_v = values.view();
    const auto n = values.shape(0);

    if constexpr (IsMasked) {
      using MaskT = nb::ndarray<const uint8_t, nb::ndim<1>>;
      MaskT mask;
      if (!nb::try_cast<MaskT>(mask_obj, mask, false)) {
        throw std::invalid_argument("Could not convert mask to uint8_t array!");
      }
      nb::call_guard<nb::gil_scoped_release>();

      for (auto i = decltype(n){0}; i < n; i++) {
        if constexpr (IgnoreNA) {
          // TODO: current pandas code is a bit messy here...
          // labels[i] = na_sentinel
          // continue;
        }

        const auto val = values_v(i);
        auto k = hash_map_.get(val);
        if (k == hash_map_.end()) {
          int dummy;
          k = hash_map_.put(val, &dummy);
          uniques.Append(val);
        }
      }
    } else {
      nb::call_guard<nb::gil_scoped_release>();
      for (auto i = decltype(n){0}; i < n; i++) {
        const auto val = values_v(i);
        auto k = hash_map_.get(val);
        if (k == hash_map_.end()) {
          int dummy;
          k = hash_map_.put(val, &dummy);
          uniques.Append(val);
        }
      }
    }

    return;
  }

  klib::KHashMap<T, size_t, PandasHashFunction<T>, PandasHashEquality<T>, HashValueT>
      hash_map_;
  Py_ssize_t na_position_ = -1;
};

using namespace nb::literals;

#define BIND_VECTOR(TYPE, NAME)                                                \
  do {                                                                         \
    nb::class_<PandasVector<TYPE>>(m, NAME)                                    \
        .def(nb::init<>())                                                     \
        .def("__len__", &PandasVector<TYPE>::__len__)                          \
        .def("append", &PandasVector<TYPE>::Append)                            \
        .def("extend", &PandasVector<TYPE>::Extend)                            \
        .def("to_array", &PandasVector<TYPE>::ToNdArray);                      \
  } while (0)

#define BIND_HASHTABLE(TYPE, NAME, MASKED)                                     \
  do {                                                                         \
    nb::class_<PandasHashTable<TYPE, MASKED>>(m, NAME)                         \
        .def(nb::init<>())                                                     \
        .def(nb::init<size_t>(), "size_hint"_a)                           \
        .def("__len__", &PandasHashTable<TYPE, MASKED>::__len__)               \
        .def("__contains__", &PandasHashTable<TYPE, MASKED>::__contains__)     \
        .def("sizeof", &PandasHashTable<TYPE, MASKED>::SizeOf)                 \
        .def("get_state", &PandasHashTable<TYPE, MASKED>::GetState)            \
        .def("get_item", &PandasHashTable<TYPE, MASKED>::GetItem)              \
        .def("get_na", &PandasHashTable<TYPE, MASKED>::GetNA)                  \
        .def("set_item", &PandasHashTable<TYPE, MASKED>::SetItem)              \
        .def("set_na", &PandasHashTable<TYPE, MASKED>::SetNA)                  \
        .def("map_keys_to_values",                                             \
             &PandasHashTable<TYPE, MASKED>::MapKeysToValues)                  \
        .def("map_locations", &PandasHashTable<TYPE, MASKED>::MapLocations,    \
             "values"_a, "mask"_a = nb::none())                                \
        .def("lookup", &PandasHashTable<TYPE, MASKED>::Lookup, "values"_a,     \
             "mask"_a = nb::none())                                            \
        .def("unique", &PandasHashTable<TYPE, MASKED>::Unique, "values"_a,     \
             "return_inverse"_a = false, "mask"_a = nb::none())                \
        .def("factorize", &PandasHashTable<TYPE, MASKED>::Factorize,           \
             "values"_a, "na_sentinel"_a = -1, "na_value"_a = nb::none(),      \
             "mask"_a = nb::none(), "ignore_na"_a = true)                      \
        .def("get_labels", &PandasHashTable<TYPE, MASKED>::GetLabels,          \
             "values"_a, "uniques"_a, "count_prior"_a = 0,                     \
             "na_sentinel"_a = -1, "na_value"_a = nb::none(),                  \
             "mask"_a = nb::none(), "ignore_na"_a = true)                      \
        .def("get_labels_groupby",                                             \
             &PandasHashTable<TYPE, MASKED>::GetLabelsGroupby);                \
  } while (0)

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

  BIND_HASHTABLE(int8_t, "Int8HashTable", false);
  BIND_HASHTABLE(int8_t, "Int8MaskedHashTable", true);
  BIND_HASHTABLE(int16_t, "Int16HashTable", false);
  BIND_HASHTABLE(int16_t, "Int16MaskedHashTable", true);
  BIND_HASHTABLE(int32_t, "Int32HashTable", false);
  BIND_HASHTABLE(int32_t, "Int32MaskedHashTable", true);
  BIND_HASHTABLE(int64_t, "Int64HashTable", false);
  BIND_HASHTABLE(int64_t, "Int64MaskedHashTable", true);
  BIND_HASHTABLE(uint8_t, "UInt8HashTable", false);
  BIND_HASHTABLE(uint8_t, "UInt8MaskedHashTable", true);
  BIND_HASHTABLE(uint16_t, "UInt16HashTable", false);
  BIND_HASHTABLE(uint16_t, "UInt16MaskedHashTable", true);
  BIND_HASHTABLE(uint32_t, "UInt32HashTable", false);
  BIND_HASHTABLE(uint32_t, "UInt32MaskedHashTable", true);
  BIND_HASHTABLE(uint64_t, "UInt64HashTable", false);
  BIND_HASHTABLE(uint64_t, "UInt64MaskedHashTable", true);
  BIND_HASHTABLE(float, "Float32HashTable", false);
  BIND_HASHTABLE(float, "Float32MaskedHashTable", true);
  BIND_HASHTABLE(double, "Float64HashTable", false);
  BIND_HASHTABLE(double, "Float64MaskedHashTable", true);
}
