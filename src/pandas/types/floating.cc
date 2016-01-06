// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include "pandas/types/floating.h"

namespace pandas {

template <typename T>
PyObject* FloatingArrayImpl<T>::GetValue(size_t i) {
  return NULL;
}

template <typename T>
void FloatingArrayImpl<T>::SetValue(size_t i, PyObject* val) {
  return;
}

// Instantiate templates
template class FloatingArrayImpl<FloatType>;
template class FloatingArrayImpl<DoubleType>;

} // namespace pandas
