#include <cstring>
#include <vector>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

template <typename T>
class PandasVector {
public:

  size_t __len__ () const {
    return vec_.size();
  }

  void Append(const T& value) {
    vec_.emplace_back(std::forward(value));
  }

  void Extend(const nb::ndarray<const T, nb::ndim<1>>& values) {
    const size_t n_elems = values.shape(0);
    vec_.reserve(vec_.size() + n_elems);
    memcpy(vec_.data() + vec_.size(), values.data(), sizeof(T) * n_elems);
  }

private:
  std::vector<T> vec_;
};

NB_MODULE(new_vector, m) {
  nb::class_<PandasVector<int64_t>>(m, "Int64Vector")
    .def("__len__", &PandasVector<int64_t>::__len__)
    .def("append", &PandasVector<int64_t>::Append)
    .def("extend", &PandasVector<int64_t>::Extend)
    ;
}
