#define PY_SSIZE_T_CLEAN

#include "Python.h"
#include "simdjson.h"

namespace pandas {
namespace json {
using namespace simdjson;

ondemand::parser parser;

static PyObject *build_python_object(ondemand::value element);

static PyObject *object_to_dict(ondemand::object element) {
  PyObject *dict = PyDict_New();
  for (auto field : element) {
    std::string_view key = field.unescaped_key();
    PyObject *value = build_python_object(field.value());

    if (!value) {
      Py_DECREF(dict);
      return NULL;
    }

    PyObject *key_py = PyUnicode_FromStringAndSize(key.data(), key.size());
    PyDict_SetItem(dict, key_py, value);
    Py_DECREF(key_py);
    Py_DECREF(value);
  }

  return dict;
}

static PyObject *array_to_list(ondemand::array element) {
  PyObject *list = PyList_New(0);
  for (auto child : element) {
    PyObject *tmp = build_python_object(child.value());
    if (!tmp) {
      Py_DECREF(list);
      return NULL;
    }

    if (PyList_Append(list, tmp) != 0) {
      Py_DECREF(list);
      Py_DECREF(tmp);
      return NULL;
    }

    Py_DECREF(tmp);
  }
  return list;
}

static PyObject *big_int_to_pylong(ondemand::value element) {
  std::string_view s = element.raw_json_token();
  std::string null_terminated_s(s);
  return PyLong_FromString(null_terminated_s.c_str(), NULL, 10);
}

static PyObject *json_number_to_pyobject(ondemand::value element) {
  ondemand::number num = element.get_number();
  switch (num.get_number_type()) {
  case ondemand::number_type::signed_integer:
    return PyLong_FromLongLong(num.get_int64());
    break;
  case ondemand::number_type::unsigned_integer:
    return PyLong_FromUnsignedLongLong(num.get_uint64());
    break;
  case ondemand::number_type::floating_point_number:
    return PyFloat_FromDouble(num.get_double());
    break;
  case ondemand::number_type::big_integer:
    return big_int_to_pylong(element);
    break;
  }
}

static PyObject *json_str_to_pyobject(ondemand::value element) {
  std::string_view s = element.get_string(true);
  return PyUnicode_FromStringAndSize(s.data(), s.size());
}

static PyObject *build_python_object(ondemand::value element) {
  switch (element.type()) {
  case ondemand::json_type::object:
    return object_to_dict(element.get_object());
    break;
  case ondemand::json_type::array:
    return array_to_list(element.get_array());
    break;
  case ondemand::json_type::boolean:
    return element.get_bool() ? Py_True : Py_False;
    break;
  case ondemand::json_type::null:
    return Py_None;
  case ondemand::json_type::string:
    return json_str_to_pyobject(element);
    break;
  case ondemand::json_type::number:
    return json_number_to_pyobject(element);
    break;
  case ondemand::json_type::unknown:
    // TODO: improve error hadling
    PyErr_Format(PyExc_ValueError, "Some error occourred");
    break;
  }

  return NULL;
}

} // namespace json
} // namespace pandas

extern "C" {

PyObject *json_loads(PyObject *Py_UNUSED(self), PyObject *args,
                     PyObject *kwargs) {
  static const char *kwlist[] = {"obj", "precise_float", NULL};
  const char *buf;
  Py_ssize_t len;
  int *precise_float; // Unused. It's declared for compatibility with old parser
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|b", kwlist, &buf, &len,
                                   &precise_float)) {
    return NULL;
  }

  PyObject *ret;
  try {
    simdjson::padded_string padded_json(buf, len);
    simdjson::ondemand::document doc =
        pandas::json::parser.iterate(padded_json);
    switch (doc.type()) {
    case simdjson::fallback::ondemand::json_type::null:
      ret = Py_None;
      break;
    case simdjson::fallback::ondemand::json_type::boolean:
      ret = doc.get_bool() ? Py_True : Py_False;
      break;
    case simdjson::fallback::ondemand::json_type::number: {
      simdjson::ondemand::number num = doc.get_number();
      switch (num.get_number_type()) {
      case simdjson::ondemand::number_type::signed_integer:
        ret = PyLong_FromLongLong(num.get_int64());
        break;
      case simdjson::ondemand::number_type::unsigned_integer:
        ret = PyLong_FromUnsignedLongLong(num.get_uint64());
        break;
      case simdjson::ondemand::number_type::floating_point_number:
        ret = PyFloat_FromDouble(num.get_double());
        break;
      case simdjson::ondemand::number_type::big_integer:
        PyErr_Format(PyExc_ValueError, "Overflow");
        return NULL;
      }
      break;
    }
    case simdjson::fallback::ondemand::json_type::string: {
      std::string_view s = doc.get_string();
      ret = PyUnicode_FromStringAndSize(s.data(), s.size());
      break;
    }
    default:
      simdjson::ondemand::value val = doc;
      ret = pandas::json::build_python_object(val);
      break;
    }
  } catch (simdjson::simdjson_error &error) {
    Py_XDECREF(ret);
    ret = NULL;
    // TODO: get location or token where error occourred
    PyErr_Format(PyExc_ValueError, "JSON parsing error: %s", error.what());
    return NULL;
  }

  return ret;
}

static PyMethodDef json_methods[] = {
    {"simdjson_loads", (PyCFunction)(void (*)(void))json_loads,
     METH_VARARGS | METH_KEYWORDS, "Parse JSON string using simdjson"},
    {NULL, NULL, 0, NULL} /* sentinel */
};

static struct PyModuleDef json_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "pandas._libs.simdjson",
    .m_doc = "simdjson python binding",
    .m_size = 0,
    .m_methods = json_methods,
    .m_slots = NULL,
    .m_traverse = NULL,
    .m_clear = NULL,
    .m_free = NULL,
};

PyMODINIT_FUNC PyInit_simdjson(void) { return PyModuleDef_Init(&json_module); }

} // extern "C"
