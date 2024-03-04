#include <cstring>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <vector>

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

NB_MODULE(new_vector, m) {
  nb::class_<PandasVector<int64_t>>(m, "Int64Vector")
      .def(nb::init<>())
      .def("__len__", &PandasVector<int64_t>::__len__)
      .def("append", &PandasVector<int64_t>::Append)
      .def("extend", &PandasVector<int64_t>::Extend)
      .def("to_array", &PandasVector<int64_t>::ToNdArray);
}
