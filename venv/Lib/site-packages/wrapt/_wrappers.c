/* ------------------------------------------------------------------------- */

#include "Python.h"

#include "structmember.h"

/* ------------------------------------------------------------------------- */

typedef struct
{
  PyObject_HEAD

      PyObject *dict;
  PyObject *wrapped;
  PyObject *weakreflist;
  int init_called;
} WraptObjectProxyObject;

typedef struct
{
  WraptObjectProxyObject object_proxy;

  PyObject *args;
  PyObject *kwargs;
} WraptPartialCallableObjectProxyObject;

typedef struct
{
  WraptObjectProxyObject object_proxy;

  PyObject *instance;
  PyObject *wrapper;
  PyObject *enabled;
  PyObject *binding;
  PyObject *parent;
  PyObject *owner;
} WraptFunctionWrapperObject;

/* ------------------------------------------------------------------------- */

/* Forward declaration of moduledef so module-state helpers can reference it. */

static struct PyModuleDef moduledef;

/* Per-interpreter module state. Holds the heap type objects so that methods
 * can reach them via PyType_GetModuleByDef without using static globals.
 *
 * Cached interned strings will be added here in a later stage of the
 * sub-interpreter / free-threading migration. */

typedef struct
{
  PyTypeObject *ObjectProxy_Type;
  PyTypeObject *CallableObjectProxy_Type;
  PyTypeObject *PartialCallableObjectProxy_Type;
  PyTypeObject *FunctionWrapperBase_Type;
  PyTypeObject *BoundFunctionWrapper_Type;
  PyTypeObject *FunctionWrapper_Type;

  /* Cached interned attribute / argument names. Initialized eagerly in
   * wrapt_exec, released in wrapt_clear. Per-interpreter so they remain
   * sub-interpreter safe; eager init makes them free-threading safe by
   * eliminating the lazy-init data race that the previous static-local
   * string pattern had. */

  PyObject *str_wrapped;                /* "__wrapped__" */
  PyObject *str_wrapped_factory;        /* "__wrapped_factory__" */
  PyObject *str_wrapped_get;            /* "__wrapped_get__" */
  PyObject *str_module;                 /* "__module__" */
  PyObject *str_doc;                    /* "__doc__" */
  PyObject *str_setattr_fixups;         /* "__wrapped_setattr_fixups__" */
  PyObject *str_getattr;                /* "__getattr__" */
  PyObject *str_self_;                  /* "_self_" */
  PyObject *str_callable;               /* "callable" */
  PyObject *str_bound_function_wrapper; /* "__bound_function_wrapper__" */
  PyObject *str_function;               /* "function" */
  PyObject *str_classmethod;            /* "classmethod" */
  PyObject *str_staticmethod;           /* "staticmethod" */
  PyObject *str_builtin;                /* "builtin" */
  PyObject *str_class;                  /* "class" */
  PyObject *str_instancemethod;         /* "instancemethod" */
  PyObject *str_object_proxy;           /* "__object_proxy__" */
  PyObject *str_enter;                  /* "__enter__" */
  PyObject *str_exit;                   /* "__exit__" */
  PyObject *str_aenter;                 /* "__aenter__" */
  PyObject *str_aexit;                  /* "__aexit__" */
  PyObject *str_mro_entries;            /* "__mro_entries__" */
  PyObject *str_name;                   /* "__name__" */
  PyObject *str_qualname;               /* "__qualname__" */
  PyObject *str_annotations;            /* "__annotations__" */
  PyObject *str_dunder_class;           /* "__class__" */
  PyObject *str_self;                   /* "__self__" */
  PyObject *str_set_name;               /* "__set_name__" */
  PyObject *str_self_binding;           /* "_self_binding" */
  PyObject *str_dict;                   /* "__dict__" */

  /* Cached exception type from wrapt.wrappers. Initialized eagerly in
   * wrapt_exec after type creation. The wrapt.wrappers module is guaranteed
   * to be imported before _wrappers because __wrapt__.py imports it first. */

  PyObject *WrapperNotInitializedError;
} wrapt_module_state;

static inline wrapt_module_state *wrapt_get_state(PyObject *module)
{
  return (wrapt_module_state *)PyModule_GetState(module);
}

/* Polyfill PyType_GetModuleByDef for Python < 3.11. The borrowed-reference
 * contract matches the real API so callers are version-agnostic. */

#if PY_VERSION_HEX < 0x030B0000
static PyObject *PyType_GetModuleByDef(PyTypeObject *type,
                                       struct PyModuleDef *def)
{
  PyObject *mro = type->tp_mro;
  Py_ssize_t i, n;

  if (mro == NULL)
    return NULL;

  n = PyTuple_GET_SIZE(mro);

  for (i = 0; i < n; i++)
  {
    PyTypeObject *t = (PyTypeObject *)PyTuple_GET_ITEM(mro, i);
    PyObject *m;

    if (!PyType_HasFeature(t, Py_TPFLAGS_HEAPTYPE))
      continue;

    m = PyType_GetModule(t);

    if (m != NULL && PyModule_GetDef(m) == def)
      return m;

    if (m == NULL)
    {
      /* PyType_GetModule only raises TypeError here (no module
       * association on heap type); no other failure mode is expected.
       * Guard with WriteUnraisable for completeness. */
      if (!PyErr_ExceptionMatches(PyExc_TypeError))
        PyErr_WriteUnraisable((PyObject *)t);
      PyErr_Clear();
    }
  }

  PyErr_Format(PyExc_TypeError,
               "PyType_GetModuleByDef: No superclass of '%s' has the given "
               "module",
               type->tp_name);
  return NULL;
}
#endif

/* Polyfill PyObject_GetOptionalAttrString for Python < 3.13. Matches the
 * 3.13+ semantics: returns 1 if found, 0 if not found (AttributeError
 * only), -1 on other errors with exception set. */

#if PY_VERSION_HEX < 0x030d0000
static inline int
PyObject_GetOptionalAttrString(PyObject *obj, const char *name, PyObject **result)
{
  *result = PyObject_GetAttrString(obj, name);
  if (*result)
    return 1;
  if (PyErr_ExceptionMatches(PyExc_AttributeError))
  {
    PyErr_Clear();
    return 0;
  }
  return -1;
}
#endif

/* Get module state for the wrapt module given any type whose MRO includes
 * one of our heap types. Returns NULL and sets an exception on miss. */

static wrapt_module_state *wrapt_state_from_type(PyTypeObject *type)
{
  PyObject *m = PyType_GetModuleByDef(type, &moduledef);
  return m ? wrapt_get_state(m) : NULL;
}

/* Walk type's MRO checking each tp_dict directly for the named attribute.
 * Does not invoke type's tp_getattro, so it does NOT trigger Python's lazy
 * creation of class-level attributes such as __annotations__. Used by
 * setattro to decide "is this a local attribute on the proxy type" without
 * shadowing our getset descriptors via lazy class-dict mutation. */

static int wrapt_type_has_attr(PyTypeObject *type, PyObject *name)
{
  PyObject *mro = type->tp_mro;
  Py_ssize_t i, n;

  if (mro == NULL)
    return 0;

  n = PyTuple_GET_SIZE(mro);

  for (i = 0; i < n; i++)
  {
    PyTypeObject *t = (PyTypeObject *)PyTuple_GET_ITEM(mro, i);
    PyObject *dict = t->tp_dict;
    int contains;

    if (dict == NULL)
      continue;

    contains = PyDict_Contains(dict, name);
    if (contains < 0)
    {
      /* PyDict_Contains cannot fail with string keys against a
       * string-keyed tp_dict; no failure is expected here. Guard
       * with WriteUnraisable for completeness. */
      PyErr_WriteUnraisable((PyObject *)t);
      PyErr_Clear();
      continue;
    }
    if (contains > 0)
      return 1;
  }

  return 0;
}

/* Non-error-setting predicate: is this object an instance of any wrapt
 * proxy type? Walks the MRO directly. PyType_GetModule sets an exception
 * for heap types that lack a module association (e.g. pure-Python user
 * subclasses), so we save/restore exception state to keep this predicate
 * side-effect-free for the caller. */

static int wrapt_is_proxy(PyObject *o)
{
  PyTypeObject *type = Py_TYPE(o);
  PyObject *mro = type->tp_mro;
  Py_ssize_t i, n;
  PyObject *saved_type = NULL, *saved_value = NULL, *saved_tb = NULL;
  int had_error = PyErr_Occurred() != NULL;
  int result = 0;

  if (had_error)
    PyErr_Fetch(&saved_type, &saved_value, &saved_tb);

  if (mro == NULL)
    goto done;

  n = PyTuple_GET_SIZE(mro);

  for (i = 0; i < n; i++)
  {
    PyTypeObject *t = (PyTypeObject *)PyTuple_GET_ITEM(mro, i);
    PyObject *m;

    if (!PyType_HasFeature(t, Py_TPFLAGS_HEAPTYPE))
      continue;

    m = PyType_GetModule(t);

    if (m == NULL)
    {
      /* PyType_GetModule only raises TypeError here (no module
       * association on heap type); no other failure mode is expected.
       * Guard with WriteUnraisable for completeness. */
      if (!PyErr_ExceptionMatches(PyExc_TypeError))
        PyErr_WriteUnraisable((PyObject *)t);
      PyErr_Clear();
      continue;
    }

    if (PyModule_GetDef(m) == &moduledef)
    {
      result = 1;
      break;
    }
  }

done:
  if (had_error)
    PyErr_Restore(saved_type, saved_value, saved_tb);
  else
    PyErr_Clear();
  return result;
}

/* ------------------------------------------------------------------------- */

static int raise_uninitialized_wrapper_error(WraptObjectProxyObject *object)
{
  // Before raising an exception we need to first check whether this is a lazy
  // proxy object and attempt to intialize the wrapped object using the supplied
  // callback if so. If the callback is not present then we can proceed to raise
  // the exception, but if the callback is present and returns a value, we set
  // that as the wrapped object and continue and return without raising an error.

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(object));
  if (!state)
    return -1;

  PyObject *wrapped_str = state->str_wrapped;
  PyObject *wrapped_factory_str = state->str_wrapped_factory;
  PyObject *wrapped_get_str = state->str_wrapped_get;

  PyObject *callback = NULL;
  PyObject *value = NULL;

  // Note that we use existance of __wrapped_factory__ to gate whether we can
  // attempt to initialize the wrapped object lazily, but it is __wrapped_get__
  // that we actually call to do the initialization. This is so that we can
  // handle multithreading correctly by having __wrapped_get__ use a lock to
  // protect against multiple threads trying to initialize the wrapped object
  // at the same time.

  callback = PyObject_GenericGetAttr((PyObject *)object, wrapped_factory_str);

  if (callback)
  {
    if (callback != Py_None)
    {
      Py_DECREF(callback);

      callback = PyObject_GenericGetAttr((PyObject *)object, wrapped_get_str);

      if (!callback)
        return -1;

      value = PyObject_CallObject(callback, NULL);

      Py_DECREF(callback);

      if (value)
      {
        // We use setattr so that special dunder methods will be properly set.

        if (PyObject_SetAttr((PyObject *)object, wrapped_str, value) == -1)
        {
          Py_DECREF(value);
          return -1;
        }

        Py_DECREF(value);

        return 0;
      }
      else
      {
        return -1;
      }
    }
    else
    {
      Py_DECREF(callback);
    }
  }
  else
  {
    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
      return -1;
    PyErr_Clear();
  }

  if (object->init_called)
  {
    PyErr_SetString(state->WrapperNotInitializedError,
                    "wrapper is in an inconsistent state: __wrapped__ is not set");
  }
  else
  {
    PyErr_Format(PyExc_AttributeError,
                 "'%.100s' object has no attribute '__wrapped__'",
                 Py_TYPE(object)->tp_name);
  }

  return -1;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_new(PyTypeObject *type, PyObject *args,
                                      PyObject *kwds)
{
  WraptObjectProxyObject *self;

  self = (WraptObjectProxyObject *)type->tp_alloc(type, 0);

  if (!self)
    return NULL;

  self->dict = PyDict_New();

  if (!self->dict)
  {
    Py_DECREF(self);
    return NULL;
  }

  self->wrapped = NULL;
  self->weakreflist = NULL;
  self->init_called = 0;

  return (PyObject *)self;
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_raw_init(WraptObjectProxyObject *self,
                                     PyObject *wrapped)
{
  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  PyObject *module_str = state->str_module;
  PyObject *doc_str = state->str_doc;
  PyObject *wrapped_factory_str = state->str_wrapped_factory;

  PyObject *object = NULL;

  // If wrapped is Py_None and we have a __wrapped_factory__ attribute
  // then we defer initialization of the wrapped object until it is first needed.

  if (wrapped == Py_None)
  {
    PyObject *callback = NULL;

    callback = PyObject_GenericGetAttr((PyObject *)self, wrapped_factory_str);

    if (callback)
    {
      if (callback != Py_None)
      {
        Py_DECREF(callback);
        self->init_called = 1;
        return 0;
      }
      else
      {
        Py_DECREF(callback);
      }
    }
    else
    {
      if (!PyErr_ExceptionMatches(PyExc_AttributeError))
        return -1;
      PyErr_Clear();
    }
  }

  Py_INCREF(wrapped);
  Py_XSETREF(self->wrapped, wrapped);

  self->init_called = 1;

  object = PyObject_GetAttr(wrapped, module_str);

  if (object)
  {
    if (PyDict_SetItem(self->dict, module_str, object) == -1)
    {
      Py_DECREF(object);
      return -1;
    }
    Py_DECREF(object);
  }
  else
  {
    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
      return -1;
    PyErr_Clear();
  }

  object = PyObject_GetAttr(wrapped, doc_str);

  if (object)
  {
    if (PyDict_SetItem(self->dict, doc_str, object) == -1)
    {
      Py_DECREF(object);
      return -1;
    }
    Py_DECREF(object);
  }
  else
  {
    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
      return -1;
    PyErr_Clear();
  }

  return 0;
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_init(WraptObjectProxyObject *self, PyObject *args,
                                 PyObject *kwds)
{
  PyObject *wrapped = NULL;

  char *const kwlist[] = {"wrapped", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:ObjectProxy", kwlist,
                                   &wrapped))
  {
    return -1;
  }

  return WraptObjectProxy_raw_init(self, wrapped);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_traverse(WraptObjectProxyObject *self,
                                     visitproc visit, void *arg)
{
  Py_VISIT(Py_TYPE(self));
  Py_VISIT(self->dict);
  Py_VISIT(self->wrapped);

  return 0;
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_clear(WraptObjectProxyObject *self)
{
  Py_CLEAR(self->dict);
  Py_CLEAR(self->wrapped);

  return 0;
}

/* ------------------------------------------------------------------------- */

static void WraptObjectProxy_dealloc(WraptObjectProxyObject *self)
{
  PyTypeObject *tp = Py_TYPE(self);

  PyObject_GC_UnTrack(self);

  if (self->weakreflist != NULL)
    PyObject_ClearWeakRefs((PyObject *)self);

  WraptObjectProxy_clear(self);

  tp->tp_free(self);

#if PY_VERSION_HEX >= 0x030C0000
  PyObject *exc = PyErr_GetRaisedException();
  Py_DECREF(tp);
  PyErr_SetRaisedException(exc);
#else
  PyObject *exc_type, *exc_value, *exc_tb;
  PyErr_Fetch(&exc_type, &exc_value, &exc_tb);
  Py_DECREF(tp);
  PyErr_Restore(exc_type, exc_value, exc_tb);
#endif
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_repr(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  char self_addr[24], wrapped_addr[24];
  snprintf(self_addr, sizeof(self_addr), "0x%zx", (Py_uintptr_t)self);
  snprintf(wrapped_addr, sizeof(wrapped_addr), "0x%zx",
           (Py_uintptr_t)self->wrapped);
  return PyUnicode_FromFormat("<%s at %s for %s at %s>",
                              Py_TYPE(self)->tp_name, self_addr,
                              Py_TYPE(self->wrapped)->tp_name, wrapped_addr);
}

/* ------------------------------------------------------------------------- */

static Py_hash_t WraptObjectProxy_hash(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  return PyObject_Hash(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_str(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_Str(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_add(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Add(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_subtract(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Subtract(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_multiply(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Multiply(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_remainder(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Remainder(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_divmod(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Divmod(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_power(PyObject *o1, PyObject *o2,
                                        PyObject *modulo)
{
  /*
   * If modulo is a proxy we cannot call PyNumber_Power with it, otherwise
   * we would recurse back into this slot via CPython's ternary_op fallback
   * and blow the C stack. Return NotImplemented so a TypeError is raised,
   * matching the pure Python implementation which does not unwrap modulo.
   */

  if (wrapt_is_proxy(modulo))
  {
    Py_RETURN_NOTIMPLEMENTED;
  }

  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Power(o1, o2, modulo);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_negative(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyNumber_Negative(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_positive(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyNumber_Positive(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_absolute(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyNumber_Absolute(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_bool(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  return PyObject_IsTrue(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_invert(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyNumber_Invert(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_lshift(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Lshift(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_rshift(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Rshift(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_and(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_And(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_xor(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Xor(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_or(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_Or(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_long(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyNumber_Long(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_float(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyNumber_Float(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_add(WraptObjectProxyObject *self,
                                              PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__iadd__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceAdd(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Add(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_subtract(WraptObjectProxyObject *self,
                                                   PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__isub__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceSubtract(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Subtract(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_multiply(WraptObjectProxyObject *self,
                                                   PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__imul__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceMultiply(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Multiply(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptObjectProxy_inplace_remainder(WraptObjectProxyObject *self,
                                   PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__imod__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceRemainder(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Remainder(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_power(WraptObjectProxyObject *self,
                                                PyObject *other,
                                                PyObject *modulo)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__ipow__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlacePower(self->wrapped, other, modulo);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Power(self->wrapped, other, modulo);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_lshift(WraptObjectProxyObject *self,
                                                 PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__ilshift__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceLshift(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Lshift(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_rshift(WraptObjectProxyObject *self,
                                                 PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__irshift__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceRshift(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Rshift(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_and(WraptObjectProxyObject *self,
                                              PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__iand__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceAnd(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_And(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_xor(WraptObjectProxyObject *self,
                                              PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__ixor__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceXor(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Xor(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_or(WraptObjectProxyObject *self,
                                             PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__ior__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceOr(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_Or(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_floor_divide(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_FloorDivide(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_true_divide(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_TrueDivide(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptObjectProxy_inplace_floor_divide(WraptObjectProxyObject *self,
                                      PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__ifloordiv__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceFloorDivide(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_FloorDivide(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptObjectProxy_inplace_true_divide(WraptObjectProxyObject *self,
                                     PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__itruediv__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceTrueDivide(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_TrueDivide(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_index(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyNumber_Index(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_matrix_multiply(PyObject *o1, PyObject *o2)
{
  if (wrapt_is_proxy(o1))
  {
    if (!((WraptObjectProxyObject *)o1)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o1) == -1)
        return NULL;
    }

    o1 = ((WraptObjectProxyObject *)o1)->wrapped;
  }

  if (wrapt_is_proxy(o2))
  {
    if (!((WraptObjectProxyObject *)o2)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)o2) == -1)
        return NULL;
    }

    o2 = ((WraptObjectProxyObject *)o2)->wrapped;
  }

  return PyNumber_MatrixMultiply(o1, o2);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_inplace_matrix_multiply(
    WraptObjectProxyObject *self, PyObject *other)
{
  PyObject *object = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (wrapt_is_proxy(other))
  {
    if (!((WraptObjectProxyObject *)other)->wrapped)
    {
      if (raise_uninitialized_wrapper_error((WraptObjectProxyObject *)other) == -1)
        return NULL;
    }

    other = ((WraptObjectProxyObject *)other)->wrapped;
  }

  PyObject *attr;
  int rc = PyObject_GetOptionalAttrString(self->wrapped, "__imatmul__", &attr);
  Py_XDECREF(attr);
  if (rc < 0)
    return NULL;
  if (rc)
  {
    object = PyNumber_InPlaceMatrixMultiply(self->wrapped, other);

    if (!object)
      return NULL;

    Py_SETREF(self->wrapped, object);

    Py_INCREF(self);
    return (PyObject *)self;
  }
  else
  {
    PyObject *result = PyNumber_MatrixMultiply(self->wrapped, other);

    if (!result)
      return NULL;

    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_type = PyObject_GetAttr((PyObject *)self, state->str_object_proxy);

    if (!proxy_type)
    {
      Py_DECREF(result);
      return NULL;
    }

    PyObject *proxy_args = PyTuple_Pack(1, result);

    Py_DECREF(result);

    if (!proxy_args)
    {
      Py_DECREF(proxy_type);
      return NULL;
    }

    PyObject *proxy_instance = PyObject_Call(proxy_type, proxy_args, NULL);

    Py_DECREF(proxy_type);
    Py_DECREF(proxy_args);

    return proxy_instance;
  }
}

/* ------------------------------------------------------------------------- */

static Py_ssize_t WraptObjectProxy_length(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  return PyObject_Length(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_contains(WraptObjectProxyObject *self,
                                     PyObject *value)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  return PySequence_Contains(self->wrapped, value);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_getitem(WraptObjectProxyObject *self,
                                          PyObject *key)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_GetItem(self->wrapped, key);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_setitem(WraptObjectProxyObject *self, PyObject *key,
                                    PyObject *value)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  if (value == NULL)
    return PyObject_DelItem(self->wrapped, key);
  else
    return PyObject_SetItem(self->wrapped, key, value);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_self_setattr(WraptObjectProxyObject *self,
                                               PyObject *args)
{
  PyObject *name = NULL;
  PyObject *value = NULL;

  if (!PyArg_ParseTuple(args, "UO:__self_setattr__", &name, &value))
    return NULL;

  if (PyObject_GenericSetAttr((PyObject *)self, name, value) != 0)
  {
    return NULL;
  }

  Py_RETURN_NONE;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_dir(WraptObjectProxyObject *self,
                                      PyObject *args)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_Dir(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_enter(WraptObjectProxyObject *self,
                                        PyObject *args, PyObject *kwds)
{
  PyObject *method = NULL;
  PyObject *result = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  method = PyObject_GetAttr(self->wrapped, state->str_enter);

  if (!method)
    return NULL;

  result = PyObject_Call(method, args, kwds);

  Py_DECREF(method);

  return result;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_exit(WraptObjectProxyObject *self,
                                       PyObject *args, PyObject *kwds)
{
  PyObject *method = NULL;
  PyObject *result = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  method = PyObject_GetAttr(self->wrapped, state->str_exit);

  if (!method)
    return NULL;

  result = PyObject_Call(method, args, kwds);

  Py_DECREF(method);

  return result;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_aenter(WraptObjectProxyObject *self,
                                         PyObject *args, PyObject *kwds)
{
  PyObject *method = NULL;
  PyObject *result = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  method = PyObject_GetAttr(self->wrapped, state->str_aenter);

  if (!method)
    return NULL;

  result = PyObject_Call(method, args, kwds);

  Py_DECREF(method);

  return result;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_aexit(WraptObjectProxyObject *self,
                                        PyObject *args, PyObject *kwds)
{
  PyObject *method = NULL;
  PyObject *result = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  method = PyObject_GetAttr(self->wrapped, state->str_aexit);

  if (!method)
    return NULL;

  result = PyObject_Call(method, args, kwds);

  Py_DECREF(method);

  return result;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_copy(WraptObjectProxyObject *self,
                                       PyObject *Py_UNUSED(ignored))
{
  PyErr_SetString(PyExc_NotImplementedError,
                  "object proxy must define __copy__()");

  return NULL;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_deepcopy(WraptObjectProxyObject *self,
                                           PyObject *args, PyObject *kwds)
{
  PyErr_SetString(PyExc_NotImplementedError,
                  "object proxy must define __deepcopy__()");

  return NULL;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_reduce(WraptObjectProxyObject *self,
                                         PyObject *Py_UNUSED(ignored))
{
  PyErr_SetString(PyExc_NotImplementedError,
                  "object proxy must define __reduce__()");

  return NULL;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_bytes(WraptObjectProxyObject *self,
                                        PyObject *args)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_Bytes(self->wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_format(WraptObjectProxyObject *self,
                                         PyObject *args)
{
  PyObject *format_spec = NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (!PyArg_ParseTuple(args, "O:format", &format_spec))
    return NULL;

  return PyObject_Format(self->wrapped, format_spec);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_reversed(WraptObjectProxyObject *self,
                                           PyObject *args)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_CallFunctionObjArgs((PyObject *)&PyReversed_Type,
                                      self->wrapped, NULL);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_round(WraptObjectProxyObject *self,
                                        PyObject *args, PyObject *kwds)
{
  PyObject *ndigits = NULL;

  PyObject *module = NULL;
  PyObject *round = NULL;

  PyObject *result = NULL;

  char *const kwlist[] = {"ndigits", NULL};

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O:ObjectProxy", kwlist,
                                   &ndigits))
  {
    return NULL;
  }

  module = PyImport_ImportModule("builtins");

  if (!module)
    return NULL;

  round = PyObject_GetAttrString(module, "round");

  if (!round)
  {
    Py_DECREF(module);
    return NULL;
  }

  Py_DECREF(module);

  result = PyObject_CallFunctionObjArgs(round, self->wrapped, ndigits, NULL);

  Py_DECREF(round);

  return result;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_complex(WraptObjectProxyObject *self,
                                          PyObject *args)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_CallFunctionObjArgs((PyObject *)&PyComplex_Type,
                                      self->wrapped, NULL);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_mro_entries(WraptObjectProxyObject *self,
                                              PyObject *args, PyObject *kwds)
{
  PyObject *wrapped = NULL;
  PyObject *mro_entries_method = NULL;
  PyObject *result = NULL;
  int is_type = 0;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapped = self->wrapped;

  // Check if wrapped is a type (class).

  is_type = PyType_Check(wrapped);

  // If wrapped is not a type and has __mro_entries__, forward to it.

  if (!is_type)
  {
    wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
    if (!state)
      return NULL;

    mro_entries_method = PyObject_GetAttr(wrapped, state->str_mro_entries);

    if (mro_entries_method)
    {
      // Call wrapped.__mro_entries__(bases).

      result = PyObject_Call(mro_entries_method, args, kwds);

      Py_DECREF(mro_entries_method);

      return result;
    }
    else
    {
      if (!PyErr_ExceptionMatches(PyExc_AttributeError))
        return NULL;
      PyErr_Clear();
    }
  }

  return Py_BuildValue("(O)", wrapped);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_name(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  return PyObject_GetAttr(self->wrapped, state->str_name);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_set_name(WraptObjectProxyObject *self,
                                     PyObject *value)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  return PyObject_SetAttr(self->wrapped, state->str_name, value);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_qualname(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  return PyObject_GetAttr(self->wrapped, state->str_qualname);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_set_qualname(WraptObjectProxyObject *self,
                                         PyObject *value)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  return PyObject_SetAttr(self->wrapped, state->str_qualname, value);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_module(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  return PyObject_GetAttr(self->wrapped, state->str_module);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_set_module(WraptObjectProxyObject *self,
                                       PyObject *value)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  if (PyObject_SetAttr(self->wrapped, state->str_module, value) == -1)
    return -1;

  return PyDict_SetItemString(self->dict, "__module__", value);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_doc(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  return PyObject_GetAttr(self->wrapped, state->str_doc);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_set_doc(WraptObjectProxyObject *self,
                                    PyObject *value)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  if (PyObject_SetAttr(self->wrapped, state->str_doc, value) == -1)
    return -1;

  return PyDict_SetItemString(self->dict, "__doc__", value);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_class(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  return PyObject_GetAttr(self->wrapped, state->str_dunder_class);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_set_class(WraptObjectProxyObject *self,
                                      PyObject *value)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  return PyObject_SetAttr(self->wrapped, state->str_dunder_class, value);
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptObjectProxy_get_annotations(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  return PyObject_GetAttr(self->wrapped, state->str_annotations);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_set_annotations(WraptObjectProxyObject *self,
                                            PyObject *value)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  return PyObject_SetAttr(self->wrapped, state->str_annotations, value);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_dict(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  return PyObject_GetAttr(self->wrapped, state->str_dict);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_set_dict(WraptObjectProxyObject *self,
                                     PyObject *value)
{
  if (!value)
  {
    PyErr_SetString(PyExc_AttributeError,
                    "property '__dict__' of 'ObjectProxy' object has no deleter");
    return -1;
  }

  PyErr_SetString(PyExc_AttributeError,
                  "property '__dict__' of 'ObjectProxy' object has no setter");
  return -1;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_self_dict(WraptObjectProxyObject *self)
{
  if (!self->dict)
  {
    self->dict = PyDict_New();
    if (!self->dict)
      return NULL;
  }

  Py_INCREF(self->dict);
  return self->dict;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_wrapped(WraptObjectProxyObject *self)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  Py_INCREF(self->wrapped);
  return self->wrapped;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_get_object_proxy(WraptObjectProxyObject *self)
{
  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;
  Py_INCREF(state->ObjectProxy_Type);
  return (PyObject *)state->ObjectProxy_Type;
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_set_wrapped(WraptObjectProxyObject *self,
                                        PyObject *value)
{
  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  PyObject *fixups = NULL;

  if (!value)
  {
    PyErr_SetString(PyExc_TypeError, "can't delete __wrapped__ attribute");
    return -1;
  }

  Py_INCREF(value);
  Py_XSETREF(self->wrapped, value);

  fixups = PyObject_GetAttr((PyObject *)self, state->str_setattr_fixups);

  if (fixups)
  {
    PyObject *result = NULL;

    result = PyObject_CallObject(fixups, NULL);
    Py_DECREF(fixups);

    if (!result)
      return -1;

    Py_DECREF(result);
  }
  else
  {
    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
      return -1;
    PyErr_Clear();
  }

  return 0;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_getattro(WraptObjectProxyObject *self,
                                           PyObject *name)
{
  PyObject *object = NULL;
  PyObject *result = NULL;

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  /* Intercept __module__ and __doc__ so they proxy to the wrapped object.
   * These cannot be getset descriptors on heap types because that would
   * place a descriptor object in tp_dict, shadowing the type-level
   * __module__/__doc__ strings that type.__module__ reads via raw dict
   * lookup. */

  if (name == state->str_module)
    return WraptObjectProxy_get_module(self);
  if (name == state->str_doc)
    return WraptObjectProxy_get_doc(self);

  object = PyObject_GenericGetAttr((PyObject *)self, name);

  if (object)
    return object;

  if (!PyErr_ExceptionMatches(PyExc_AttributeError))
    return NULL;

  PyErr_Clear();

  object = PyObject_GenericGetAttr((PyObject *)self, state->str_getattr);

  if (!object)
    return NULL;

  result = PyObject_CallFunctionObjArgs(object, name, NULL);

  Py_DECREF(object);

  return result;
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_getattr(WraptObjectProxyObject *self,
                                          PyObject *args)
{
  PyObject *name = NULL;

  if (!PyArg_ParseTuple(args, "U:__getattr__", &name))
    return NULL;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_GetAttr(self->wrapped, name);
}

/* ------------------------------------------------------------------------- */

static int WraptObjectProxy_setattro(WraptObjectProxyObject *self,
                                     PyObject *name, PyObject *value)
{
  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  int match = PyUnicode_Tailmatch(name, state->str_self_, 0, PY_SSIZE_T_MAX, -1);

  if (match < 0)
    return -1;

  if (match)
    return PyObject_GenericSetAttr((PyObject *)self, name, value);

  /* Intercept __module__ and __doc__ so they forward to the wrapped
   * object rather than being stored locally. Without this, the string
   * values in the type dict (from PyType_FromModuleAndSpec) would cause
   * wrapt_type_has_attr to match and GenericSetAttr to store locally. */

  if (name == state->str_module)
    return WraptObjectProxy_set_module(self, value);
  if (name == state->str_doc)
    return WraptObjectProxy_set_doc(self, value);

  if (wrapt_type_has_attr(Py_TYPE(self), name))
    return PyObject_GenericSetAttr((PyObject *)self, name, value);

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return -1;
  }

  return PyObject_SetAttr(self->wrapped, name, value);
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptObjectProxy_richcompare(WraptObjectProxyObject *self,
                                              PyObject *other, int opcode)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_RichCompare(self->wrapped, other, opcode);
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptObjectProxy_instancecheck(WraptObjectProxyObject *self,
                               PyObject *instance)
{
  int check = 0;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  check = PyObject_IsInstance(instance, self->wrapped);

  if (check < 0)
  {
    return NULL;
  }

  return PyBool_FromLong(check);
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptObjectProxy_subclasscheck(WraptObjectProxyObject *self,
                               PyObject *args)
{
  PyObject *subclass = NULL;
  PyObject *object = NULL;

  int check = 0;

  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  if (!PyArg_ParseTuple(args, "O", &subclass))
    return NULL;

  wrapt_module_state *state =
      wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  object = PyObject_GetAttr(subclass, state->str_wrapped);

  if (!object)
  {
    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
      return NULL;
    PyErr_Clear();
  }

  check = PyObject_IsSubclass(object ? object : subclass,
                              self->wrapped);

  Py_XDECREF(object);

  if (check == -1)
    return NULL;

  return PyBool_FromLong(check);
}

/* ------------------------------------------------------------------------- */

static PyMethodDef WraptObjectProxy_methods[] = {
    {"__self_setattr__", (PyCFunction)WraptObjectProxy_self_setattr,
     METH_VARARGS, 0},
    {"__dir__", (PyCFunction)WraptObjectProxy_dir, METH_NOARGS, 0},
    {"__enter__", (PyCFunction)WraptObjectProxy_enter,
     METH_VARARGS | METH_KEYWORDS, 0},
    {"__exit__", (PyCFunction)WraptObjectProxy_exit,
     METH_VARARGS | METH_KEYWORDS, 0},
    {"__aenter__", (PyCFunction)WraptObjectProxy_aenter,
     METH_VARARGS | METH_KEYWORDS, 0},
    {"__aexit__", (PyCFunction)WraptObjectProxy_aexit,
     METH_VARARGS | METH_KEYWORDS, 0},
    {"__copy__", (PyCFunction)WraptObjectProxy_copy, METH_NOARGS, 0},
    {"__deepcopy__", (PyCFunction)WraptObjectProxy_deepcopy,
     METH_VARARGS | METH_KEYWORDS, 0},
    {"__reduce__", (PyCFunction)WraptObjectProxy_reduce, METH_NOARGS, 0},
    {"__getattr__", (PyCFunction)WraptObjectProxy_getattr, METH_VARARGS, 0},
    {"__bytes__", (PyCFunction)WraptObjectProxy_bytes, METH_NOARGS, 0},
    {"__format__", (PyCFunction)WraptObjectProxy_format, METH_VARARGS, 0},
    {"__reversed__", (PyCFunction)WraptObjectProxy_reversed, METH_NOARGS, 0},
    {"__round__", (PyCFunction)WraptObjectProxy_round,
     METH_VARARGS | METH_KEYWORDS, 0},
    {"__complex__", (PyCFunction)WraptObjectProxy_complex, METH_NOARGS, 0},
    {"__mro_entries__", (PyCFunction)WraptObjectProxy_mro_entries,
     METH_VARARGS | METH_KEYWORDS, 0},
    {"__instancecheck__", (PyCFunction)WraptObjectProxy_instancecheck,
     METH_O, 0},
    {"__subclasscheck__", (PyCFunction)WraptObjectProxy_subclasscheck,
     METH_VARARGS, 0},
    {"__class_getitem__", Py_GenericAlias, METH_O | METH_CLASS, 0},
    {NULL, NULL},
};

static PyGetSetDef WraptObjectProxy_getset[] = {
    {"__name__", (getter)WraptObjectProxy_get_name,
     (setter)WraptObjectProxy_set_name, 0},
    {"__qualname__", (getter)WraptObjectProxy_get_qualname,
     (setter)WraptObjectProxy_set_qualname, 0},
    {"__class__", (getter)WraptObjectProxy_get_class,
     (setter)WraptObjectProxy_set_class, 0},
    {"__annotations__", (getter)WraptObjectProxy_get_annotations,
     (setter)WraptObjectProxy_set_annotations, 0},
    {"__dict__", (getter)WraptObjectProxy_get_dict,
     (setter)WraptObjectProxy_set_dict, 0},
    {"__self_dict__", (getter)WraptObjectProxy_get_self_dict, 0, 0},
    {"__wrapped__", (getter)WraptObjectProxy_get_wrapped,
     (setter)WraptObjectProxy_set_wrapped, 0},
    {"__object_proxy__", (getter)WraptObjectProxy_get_object_proxy, 0, 0},
    {NULL},
};

static PyMemberDef WraptObjectProxy_members[] = {
    {"__dictoffset__", T_PYSSIZET,
     offsetof(WraptObjectProxyObject, dict), READONLY, NULL},
    {"__weaklistoffset__", T_PYSSIZET,
     offsetof(WraptObjectProxyObject, weakreflist), READONLY, NULL},
    {NULL},
};

static PyType_Slot WraptObjectProxy_slots[] = {
    {Py_tp_dealloc, WraptObjectProxy_dealloc},
    {Py_tp_repr, WraptObjectProxy_repr},
    {Py_tp_hash, WraptObjectProxy_hash},
    {Py_tp_str, WraptObjectProxy_str},
    {Py_tp_getattro, WraptObjectProxy_getattro},
    {Py_tp_setattro, WraptObjectProxy_setattro},
    {Py_tp_traverse, WraptObjectProxy_traverse},
    {Py_tp_clear, WraptObjectProxy_clear},
    {Py_tp_richcompare, WraptObjectProxy_richcompare},
    {Py_tp_methods, WraptObjectProxy_methods},
    {Py_tp_members, WraptObjectProxy_members},
    {Py_tp_getset, WraptObjectProxy_getset},
    {Py_tp_init, WraptObjectProxy_init},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_new, WraptObjectProxy_new},
    {Py_tp_free, PyObject_GC_Del},
    /* PyNumberMethods */
    {Py_nb_add, WraptObjectProxy_add},
    {Py_nb_subtract, WraptObjectProxy_subtract},
    {Py_nb_multiply, WraptObjectProxy_multiply},
    {Py_nb_remainder, WraptObjectProxy_remainder},
    {Py_nb_divmod, WraptObjectProxy_divmod},
    {Py_nb_power, WraptObjectProxy_power},
    {Py_nb_negative, WraptObjectProxy_negative},
    {Py_nb_positive, WraptObjectProxy_positive},
    {Py_nb_absolute, WraptObjectProxy_absolute},
    {Py_nb_bool, WraptObjectProxy_bool},
    {Py_nb_invert, WraptObjectProxy_invert},
    {Py_nb_lshift, WraptObjectProxy_lshift},
    {Py_nb_rshift, WraptObjectProxy_rshift},
    {Py_nb_and, WraptObjectProxy_and},
    {Py_nb_xor, WraptObjectProxy_xor},
    {Py_nb_or, WraptObjectProxy_or},
    {Py_nb_int, WraptObjectProxy_long},
    {Py_nb_float, WraptObjectProxy_float},
    {Py_nb_inplace_add, WraptObjectProxy_inplace_add},
    {Py_nb_inplace_subtract, WraptObjectProxy_inplace_subtract},
    {Py_nb_inplace_multiply, WraptObjectProxy_inplace_multiply},
    {Py_nb_inplace_remainder, WraptObjectProxy_inplace_remainder},
    {Py_nb_inplace_power, WraptObjectProxy_inplace_power},
    {Py_nb_inplace_lshift, WraptObjectProxy_inplace_lshift},
    {Py_nb_inplace_rshift, WraptObjectProxy_inplace_rshift},
    {Py_nb_inplace_and, WraptObjectProxy_inplace_and},
    {Py_nb_inplace_xor, WraptObjectProxy_inplace_xor},
    {Py_nb_inplace_or, WraptObjectProxy_inplace_or},
    {Py_nb_floor_divide, WraptObjectProxy_floor_divide},
    {Py_nb_true_divide, WraptObjectProxy_true_divide},
    {Py_nb_inplace_floor_divide, WraptObjectProxy_inplace_floor_divide},
    {Py_nb_inplace_true_divide, WraptObjectProxy_inplace_true_divide},
    {Py_nb_index, WraptObjectProxy_index},
    {Py_nb_matrix_multiply, WraptObjectProxy_matrix_multiply},
    {Py_nb_inplace_matrix_multiply, WraptObjectProxy_inplace_matrix_multiply},
    /* PySequenceMethods */
    {Py_sq_length, WraptObjectProxy_length},
    {Py_sq_contains, WraptObjectProxy_contains},
    /* PyMappingMethods */
    {Py_mp_length, WraptObjectProxy_length},
    {Py_mp_subscript, WraptObjectProxy_getitem},
    {Py_mp_ass_subscript, WraptObjectProxy_setitem},
    {0, NULL},
};

static PyType_Spec WraptObjectProxy_spec = {
    .name = "_wrappers.ObjectProxy",
    .basicsize = sizeof(WraptObjectProxyObject),
    .itemsize = 0,
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = WraptObjectProxy_slots,
};

/* ------------------------------------------------------------------------- */

static PyObject *WraptCallableObjectProxy_call(WraptObjectProxyObject *self,
                                               PyObject *args, PyObject *kwds)
{
  if (!self->wrapped)
  {
    if (raise_uninitialized_wrapper_error(self) == -1)
      return NULL;
  }

  return PyObject_Call(self->wrapped, args, kwds);
}

/* ------------------------------------------------------------------------- */;

static PyType_Slot WraptCallableObjectProxy_slots[] = {
    {Py_tp_dealloc, WraptObjectProxy_dealloc},
    {Py_tp_traverse, WraptObjectProxy_traverse},
    {Py_tp_clear, WraptObjectProxy_clear},
    {Py_tp_call, WraptCallableObjectProxy_call},
    {Py_tp_init, WraptObjectProxy_init},
    {0, NULL},
};

static PyType_Spec WraptCallableObjectProxy_spec = {
    .name = "_wrappers.CallableObjectProxy",
    .basicsize = sizeof(WraptObjectProxyObject),
    .itemsize = 0,
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = WraptCallableObjectProxy_slots,
};

/* ------------------------------------------------------------------------- */

static PyObject *WraptPartialCallableObjectProxy_new(PyTypeObject *type,
                                                     PyObject *args,
                                                     PyObject *kwds)
{
  WraptPartialCallableObjectProxyObject *self;

  self = (WraptPartialCallableObjectProxyObject *)WraptObjectProxy_new(
      type, args, kwds);

  if (!self)
    return NULL;

  self->args = NULL;
  self->kwargs = NULL;

  return (PyObject *)self;
}

/* ------------------------------------------------------------------------- */

static int WraptPartialCallableObjectProxy_raw_init(
    WraptPartialCallableObjectProxyObject *self, PyObject *wrapped,
    PyObject *args, PyObject *kwargs)
{
  int result = 0;

  result = WraptObjectProxy_raw_init((WraptObjectProxyObject *)self, wrapped);

  if (result == 0)
  {
    Py_INCREF(args);
    Py_XSETREF(self->args, args);

    Py_XINCREF(kwargs);
    Py_XSETREF(self->kwargs, kwargs);
  }

  return result;
}

/* ------------------------------------------------------------------------- */

static int WraptPartialCallableObjectProxy_init(
    WraptPartialCallableObjectProxyObject *self, PyObject *args,
    PyObject *kwds)
{
  PyObject *wrapped = NULL;
  PyObject *fnargs = NULL;

  int result = 0;

  if (PyTuple_GET_SIZE(args) < 1)
  {
    PyErr_SetString(PyExc_TypeError, "__init__ of partial needs an argument");
    return -1;
  }

  wrapped = PyTuple_GetItem(args, 0);

  if (!PyCallable_Check(wrapped))
  {
    PyErr_SetString(PyExc_TypeError, "the first argument must be callable");
    return -1;
  }

  fnargs = PyTuple_GetSlice(args, 1, PyTuple_GET_SIZE(args));

  if (!fnargs)
    return -1;

  result =
      WraptPartialCallableObjectProxy_raw_init(self, wrapped, fnargs, kwds);

  Py_DECREF(fnargs);

  return result;
}

/* ------------------------------------------------------------------------- */

static int WraptPartialCallableObjectProxy_traverse(
    WraptPartialCallableObjectProxyObject *self, visitproc visit, void *arg)
{
  int err = WraptObjectProxy_traverse((WraptObjectProxyObject *)self, visit, arg);
  if (err)
    return err;

  Py_VISIT(self->args);
  Py_VISIT(self->kwargs);

  return 0;
}

/* ------------------------------------------------------------------------- */

static int WraptPartialCallableObjectProxy_clear(
    WraptPartialCallableObjectProxyObject *self)
{
  WraptObjectProxy_clear((WraptObjectProxyObject *)self);

  Py_CLEAR(self->args);
  Py_CLEAR(self->kwargs);

  return 0;
}

/* ------------------------------------------------------------------------- */

static void WraptPartialCallableObjectProxy_dealloc(
    WraptPartialCallableObjectProxyObject *self)
{
  PyTypeObject *tp = Py_TYPE(self);

  PyObject_GC_UnTrack(self);

  if (self->object_proxy.weakreflist != NULL)
    PyObject_ClearWeakRefs((PyObject *)self);

  WraptPartialCallableObjectProxy_clear(self);

  tp->tp_free(self);

#if PY_VERSION_HEX >= 0x030C0000
  PyObject *exc = PyErr_GetRaisedException();
  Py_DECREF(tp);
  PyErr_SetRaisedException(exc);
#else
  PyObject *exc_type, *exc_value, *exc_tb;
  PyErr_Fetch(&exc_type, &exc_value, &exc_tb);
  Py_DECREF(tp);
  PyErr_Restore(exc_type, exc_value, exc_tb);
#endif
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptPartialCallableObjectProxy_call(
    WraptPartialCallableObjectProxyObject *self, PyObject *args,
    PyObject *kwds)
{
  PyObject *fnargs = NULL;
  PyObject *fnkwargs = NULL;

  PyObject *result = NULL;

  long i;
  long offset;

  if (!self->object_proxy.wrapped)
  {
    if (raise_uninitialized_wrapper_error(&self->object_proxy) == -1)
      return NULL;
  }

  fnargs = PyTuple_New(PyTuple_GET_SIZE(self->args) + PyTuple_GET_SIZE(args));

  if (!fnargs)
    return NULL;

  for (i = 0; i < PyTuple_GET_SIZE(self->args); i++)
  {
    PyObject *item;
    item = PyTuple_GetItem(self->args, i);
    Py_INCREF(item);
    PyTuple_SetItem(fnargs, i, item);
  }

  offset = PyTuple_GET_SIZE(self->args);

  for (i = 0; i < PyTuple_GET_SIZE(args); i++)
  {
    PyObject *item;
    item = PyTuple_GetItem(args, i);
    Py_INCREF(item);
    PyTuple_SetItem(fnargs, offset + i, item);
  }

  fnkwargs = PyDict_New();

  if (!fnkwargs)
  {
    Py_DECREF(fnargs);
    return NULL;
  }

  if (self->kwargs && PyDict_Update(fnkwargs, self->kwargs) == -1)
  {
    Py_DECREF(fnargs);
    Py_DECREF(fnkwargs);
    return NULL;
  }

  if (kwds && PyDict_Update(fnkwargs, kwds) == -1)
  {
    Py_DECREF(fnargs);
    Py_DECREF(fnkwargs);
    return NULL;
  }

  result = PyObject_Call(self->object_proxy.wrapped, fnargs, fnkwargs);

  Py_DECREF(fnargs);
  Py_DECREF(fnkwargs);

  return result;
}

/* ------------------------------------------------------------------------- */;

static PyType_Slot WraptPartialCallableObjectProxy_slots[] = {
    {Py_tp_dealloc, WraptPartialCallableObjectProxy_dealloc},
    {Py_tp_call, WraptPartialCallableObjectProxy_call},
    {Py_tp_traverse, WraptPartialCallableObjectProxy_traverse},
    {Py_tp_clear, WraptPartialCallableObjectProxy_clear},
    {Py_tp_init, WraptPartialCallableObjectProxy_init},
    {Py_tp_new, WraptPartialCallableObjectProxy_new},
    {0, NULL},
};

static PyType_Spec WraptPartialCallableObjectProxy_spec = {
    .name = "_wrappers.PartialCallableObjectProxy",
    .basicsize = sizeof(WraptPartialCallableObjectProxyObject),
    .itemsize = 0,
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = WraptPartialCallableObjectProxy_slots,
};

/* ------------------------------------------------------------------------- */

static PyObject *WraptFunctionWrapperBase_new(PyTypeObject *type,
                                              PyObject *args, PyObject *kwds)
{
  WraptFunctionWrapperObject *self;

  self = (WraptFunctionWrapperObject *)WraptObjectProxy_new(type, args, kwds);

  if (!self)
    return NULL;

  self->instance = NULL;
  self->wrapper = NULL;
  self->enabled = NULL;
  self->binding = NULL;
  self->parent = NULL;
  self->owner = NULL;

  return (PyObject *)self;
}

/* ------------------------------------------------------------------------- */

static int WraptFunctionWrapperBase_raw_init(
    WraptFunctionWrapperObject *self, PyObject *wrapped, PyObject *instance,
    PyObject *wrapper, PyObject *enabled, PyObject *binding, PyObject *parent,
    PyObject *owner)
{
  int result = 0;

  result = WraptObjectProxy_raw_init((WraptObjectProxyObject *)self, wrapped);

  if (result == 0)
  {
    Py_INCREF(instance);
    Py_XSETREF(self->instance, instance);

    Py_INCREF(wrapper);
    Py_XSETREF(self->wrapper, wrapper);

    Py_INCREF(enabled);
    Py_XSETREF(self->enabled, enabled);

    Py_INCREF(binding);
    Py_XSETREF(self->binding, binding);

    Py_INCREF(parent);
    Py_XSETREF(self->parent, parent);

    Py_INCREF(owner);
    Py_XSETREF(self->owner, owner);
  }

  return result;
}

/* ------------------------------------------------------------------------- */

static int WraptFunctionWrapperBase_init(WraptFunctionWrapperObject *self,
                                         PyObject *args, PyObject *kwds)
{
  PyObject *wrapped = NULL;
  PyObject *instance = NULL;
  PyObject *wrapper = NULL;
  PyObject *enabled = Py_None;
  PyObject *binding = NULL;
  PyObject *parent = Py_None;
  PyObject *owner = Py_None;

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  char *const kwlist[] = {"wrapped", "instance", "wrapper", "enabled",
                          "binding", "parent", "owner", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|OOOO:FunctionWrapperBase",
                                   kwlist, &wrapped, &instance, &wrapper,
                                   &enabled, &binding, &parent, &owner))
  {
    return -1;
  }

  if (!binding)
    binding = state->str_callable;

  return WraptFunctionWrapperBase_raw_init(self, wrapped, instance, wrapper,
                                           enabled, binding, parent, owner);
}

/* ------------------------------------------------------------------------- */

static int WraptFunctionWrapperBase_traverse(WraptFunctionWrapperObject *self,
                                             visitproc visit, void *arg)
{
  int err = WraptObjectProxy_traverse((WraptObjectProxyObject *)self, visit, arg);
  if (err)
    return err;

  Py_VISIT(self->instance);
  Py_VISIT(self->wrapper);
  Py_VISIT(self->enabled);
  Py_VISIT(self->binding);
  Py_VISIT(self->parent);
  Py_VISIT(self->owner);

  return 0;
}

/* ------------------------------------------------------------------------- */

static int WraptFunctionWrapperBase_clear(WraptFunctionWrapperObject *self)
{
  WraptObjectProxy_clear((WraptObjectProxyObject *)self);

  Py_CLEAR(self->instance);
  Py_CLEAR(self->wrapper);
  Py_CLEAR(self->enabled);
  Py_CLEAR(self->binding);
  Py_CLEAR(self->parent);
  Py_CLEAR(self->owner);

  return 0;
}

/* ------------------------------------------------------------------------- */

static void WraptFunctionWrapperBase_dealloc(WraptFunctionWrapperObject *self)
{
  PyTypeObject *tp = Py_TYPE(self);

  PyObject_GC_UnTrack(self);

  if (self->object_proxy.weakreflist != NULL)
    PyObject_ClearWeakRefs((PyObject *)self);

  WraptFunctionWrapperBase_clear(self);

  tp->tp_free(self);

#if PY_VERSION_HEX >= 0x030C0000
  PyObject *exc = PyErr_GetRaisedException();
  Py_DECREF(tp);
  PyErr_SetRaisedException(exc);
#else
  PyObject *exc_type, *exc_value, *exc_tb;
  PyErr_Fetch(&exc_type, &exc_value, &exc_tb);
  Py_DECREF(tp);
  PyErr_Restore(exc_type, exc_value, exc_tb);
#endif
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptFunctionWrapperBase_call(WraptFunctionWrapperObject *self,
                                               PyObject *args, PyObject *kwds)
{
  PyObject *param_kwds = NULL;

  PyObject *result = NULL;

  if (!self->object_proxy.wrapped)
  {
    if (raise_uninitialized_wrapper_error(&self->object_proxy) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  if (self->enabled != Py_None)
  {
    if (PyCallable_Check(self->enabled))
    {
      PyObject *object = NULL;
      int is_false;

      object = PyObject_CallFunctionObjArgs(self->enabled, NULL);

      if (!object)
        return NULL;

      is_false = PyObject_Not(object);
      Py_DECREF(object);

      if (is_false < 0)
        return NULL;

      if (is_false)
        return PyObject_Call(self->object_proxy.wrapped, args, kwds);
    }
    else
    {
      int is_false = PyObject_Not(self->enabled);

      if (is_false < 0)
        return NULL;

      if (is_false)
        return PyObject_Call(self->object_proxy.wrapped, args, kwds);
    }
  }

  if (!kwds)
  {
    param_kwds = PyDict_New();
    if (!param_kwds)
      return NULL;
    kwds = param_kwds;
  }

  if (self->instance == Py_None)
  {
    int matched =
        PyUnicode_CompareWithASCIIString(self->binding, "function") == 0 ||
        PyUnicode_CompareWithASCIIString(self->binding, "instancemethod") == 0 ||
        PyUnicode_CompareWithASCIIString(self->binding, "callable") == 0 ||
        PyUnicode_CompareWithASCIIString(self->binding, "classmethod") == 0;

    if (matched)
    {
      PyObject *instance = NULL;

      instance =
          PyObject_GetAttr(self->object_proxy.wrapped, state->str_self);

      if (instance)
      {
        result = PyObject_CallFunctionObjArgs(self->wrapper,
                                              self->object_proxy.wrapped,
                                              instance, args, kwds, NULL);

        Py_XDECREF(param_kwds);

        Py_DECREF(instance);

        return result;
      }
      else
      {
        if (!PyErr_ExceptionMatches(PyExc_AttributeError))
        {
          Py_XDECREF(param_kwds);
          return NULL;
        }
        PyErr_Clear();
      }
    }
  }

  result =
      PyObject_CallFunctionObjArgs(self->wrapper, self->object_proxy.wrapped,
                                   self->instance, args, kwds, NULL);

  Py_XDECREF(param_kwds);

  return result;
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptFunctionWrapperBase_descr_get(WraptFunctionWrapperObject *self,
                                   PyObject *obj, PyObject *type)
{
  PyObject *bound_type = NULL;
  PyObject *descriptor = NULL;
  PyObject *result = NULL;

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  PyObject *bound_type_str = state->str_bound_function_wrapper;

  if (!self->object_proxy.wrapped)
  {
    if (raise_uninitialized_wrapper_error(&self->object_proxy) == -1)
      return NULL;
  }

  if (self->parent == Py_None)
  {
    if (PyUnicode_CompareWithASCIIString(self->binding, "builtin") == 0)
    {
      Py_INCREF(self);
      return (PyObject *)self;
    }

    if (PyUnicode_CompareWithASCIIString(self->binding, "class") == 0)
    {
      Py_INCREF(self);
      return (PyObject *)self;
    }

    if (Py_TYPE(self->object_proxy.wrapped)->tp_descr_get == NULL)
    {
      Py_INCREF(self);
      return (PyObject *)self;
    }

    descriptor = (Py_TYPE(self->object_proxy.wrapped)->tp_descr_get)(
        self->object_proxy.wrapped, obj, type);

    if (!descriptor)
      return NULL;

    if (Py_TYPE(self) != state->FunctionWrapper_Type)
    {
      bound_type = PyObject_GenericGetAttr((PyObject *)self, bound_type_str);

      if (!bound_type)
      {
        if (!PyErr_ExceptionMatches(PyExc_AttributeError))
        {
          Py_DECREF(descriptor);
          return NULL;
        }
        PyErr_Clear();
      }
    }

    if (obj == NULL)
      obj = Py_None;

    result = PyObject_CallFunctionObjArgs(
        bound_type ? bound_type : (PyObject *)state->BoundFunctionWrapper_Type,
        descriptor, obj, self->wrapper, self->enabled, self->binding, self,
        type, NULL);

    Py_XDECREF(bound_type);
    Py_DECREF(descriptor);

    return result;
  }

  if (self->instance == Py_None)
  {
    int matched =
        PyUnicode_CompareWithASCIIString(self->binding, "function") == 0 ||
        PyUnicode_CompareWithASCIIString(self->binding, "instancemethod") == 0 ||
        PyUnicode_CompareWithASCIIString(self->binding, "callable") == 0;

    if (matched)
    {
      PyObject *wrapped = NULL;

      if (PyObject_TypeCheck(self->parent, state->ObjectProxy_Type))
      {
        WraptObjectProxyObject *parent_proxy =
            (WraptObjectProxyObject *)self->parent;

        if (!parent_proxy->wrapped)
        {
          if (raise_uninitialized_wrapper_error(parent_proxy) == -1)
            return NULL;
        }

        wrapped = parent_proxy->wrapped;
        Py_INCREF(wrapped);
      }
      else
      {
        /* Fallback for the unusual case where parent is not a wrapt proxy. */
        wrapped = PyObject_GetAttr((PyObject *)self->parent, state->str_wrapped);

        if (!wrapped)
          return NULL;
      }

      if (Py_TYPE(wrapped)->tp_descr_get == NULL)
      {
        PyErr_Format(PyExc_AttributeError,
                     "'%s' object has no attribute '__get__'",
                     Py_TYPE(wrapped)->tp_name);
        Py_DECREF(wrapped);
        return NULL;
      }

      descriptor = (Py_TYPE(wrapped)->tp_descr_get)(wrapped, obj, type);

      Py_DECREF(wrapped);

      if (!descriptor)
        return NULL;

      if (Py_TYPE(self->parent) != state->FunctionWrapper_Type)
      {
        bound_type =
            PyObject_GenericGetAttr((PyObject *)self->parent, bound_type_str);

        if (!bound_type)
        {
          if (!PyErr_ExceptionMatches(PyExc_AttributeError))
          {
            Py_DECREF(descriptor);
            return NULL;
          }
          PyErr_Clear();
        }
      }

      if (obj == NULL)
        obj = Py_None;

      result = PyObject_CallFunctionObjArgs(
          bound_type ? bound_type : (PyObject *)state->BoundFunctionWrapper_Type,
          descriptor, obj, self->wrapper, self->enabled, self->binding,
          self->parent, type, NULL);

      Py_XDECREF(bound_type);
      Py_DECREF(descriptor);

      return result;
    }
  }

  Py_INCREF(self);
  return (PyObject *)self;
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptFunctionWrapperBase_set_name(WraptFunctionWrapperObject *self,
                                  PyObject *args, PyObject *kwds)
{
  PyObject *method = NULL;
  PyObject *result = NULL;

  if (!self->object_proxy.wrapped)
  {
    if (raise_uninitialized_wrapper_error(&self->object_proxy) == -1)
      return NULL;
  }

  wrapt_module_state *state =
      wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  method = PyObject_GetAttr(self->object_proxy.wrapped, state->str_set_name);

  if (!method)
  {
    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
      return NULL;
    PyErr_Clear();
    Py_RETURN_NONE;
  }

  result = PyObject_Call(method, args, kwds);

  Py_DECREF(method);

  return result;
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptFunctionWrapperBase_get_self_instance(WraptFunctionWrapperObject *self,
                                           void *closure)
{
  if (!self->instance)
  {
    Py_RETURN_NONE;
  }

  Py_INCREF(self->instance);
  return self->instance;
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptFunctionWrapperBase_get_self_wrapper(WraptFunctionWrapperObject *self,
                                          void *closure)
{
  if (!self->wrapper)
  {
    Py_RETURN_NONE;
  }

  Py_INCREF(self->wrapper);
  return self->wrapper;
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptFunctionWrapperBase_get_self_enabled(WraptFunctionWrapperObject *self,
                                          void *closure)
{
  if (!self->enabled)
  {
    Py_RETURN_NONE;
  }

  Py_INCREF(self->enabled);
  return self->enabled;
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptFunctionWrapperBase_get_self_binding(WraptFunctionWrapperObject *self,
                                          void *closure)
{
  if (!self->binding)
  {
    Py_RETURN_NONE;
  }

  Py_INCREF(self->binding);
  return self->binding;
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptFunctionWrapperBase_get_self_parent(WraptFunctionWrapperObject *self,
                                         void *closure)
{
  if (!self->parent)
  {
    Py_RETURN_NONE;
  }

  Py_INCREF(self->parent);
  return self->parent;
}

/* ------------------------------------------------------------------------- */

static PyObject *
WraptFunctionWrapperBase_get_self_owner(WraptFunctionWrapperObject *self,
                                        void *closure)
{
  if (!self->owner)
  {
    Py_RETURN_NONE;
  }

  Py_INCREF(self->owner);
  return self->owner;
}

/* ------------------------------------------------------------------------- */;

static PyMethodDef WraptFunctionWrapperBase_methods[] = {
    {"__set_name__", (PyCFunction)WraptFunctionWrapperBase_set_name,
     METH_VARARGS | METH_KEYWORDS, 0},
    {NULL, NULL},
};

/* ------------------------------------------------------------------------- */;

static PyGetSetDef WraptFunctionWrapperBase_getset[] = {
    {"_self_instance", (getter)WraptFunctionWrapperBase_get_self_instance, NULL,
     0},
    {"_self_wrapper", (getter)WraptFunctionWrapperBase_get_self_wrapper, NULL,
     0},
    {"_self_enabled", (getter)WraptFunctionWrapperBase_get_self_enabled, NULL,
     0},
    {"_self_binding", (getter)WraptFunctionWrapperBase_get_self_binding, NULL,
     0},
    {"_self_parent", (getter)WraptFunctionWrapperBase_get_self_parent, NULL, 0},
    {"_self_owner", (getter)WraptFunctionWrapperBase_get_self_owner, NULL, 0},
    {NULL},
};

static PyType_Slot WraptFunctionWrapperBase_slots[] = {
    {Py_tp_dealloc, WraptFunctionWrapperBase_dealloc},
    {Py_tp_call, WraptFunctionWrapperBase_call},
    {Py_tp_traverse, WraptFunctionWrapperBase_traverse},
    {Py_tp_clear, WraptFunctionWrapperBase_clear},
    {Py_tp_methods, WraptFunctionWrapperBase_methods},
    {Py_tp_getset, WraptFunctionWrapperBase_getset},
    {Py_tp_descr_get, WraptFunctionWrapperBase_descr_get},
    {Py_tp_init, WraptFunctionWrapperBase_init},
    {Py_tp_new, WraptFunctionWrapperBase_new},
    {0, NULL},
};

static PyType_Spec WraptFunctionWrapperBase_spec = {
    .name = "_wrappers._FunctionWrapperBase",
    .basicsize = sizeof(WraptFunctionWrapperObject),
    .itemsize = 0,
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = WraptFunctionWrapperBase_slots,
};

/* ------------------------------------------------------------------------- */

static PyObject *
WraptBoundFunctionWrapper_call(WraptFunctionWrapperObject *self, PyObject *args,
                               PyObject *kwds)
{
  PyObject *param_args = NULL;
  PyObject *param_kwds = NULL;

  PyObject *wrapped = NULL;
  PyObject *instance = NULL;

  PyObject *result = NULL;

  int matched;

  if (!self->object_proxy.wrapped)
  {
    if (raise_uninitialized_wrapper_error(&self->object_proxy) == -1)
      return NULL;
  }

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return NULL;

  if (self->enabled != Py_None)
  {
    if (PyCallable_Check(self->enabled))
    {
      PyObject *object = NULL;
      int is_false;

      object = PyObject_CallFunctionObjArgs(self->enabled, NULL);

      if (!object)
        return NULL;

      is_false = PyObject_Not(object);
      Py_DECREF(object);

      if (is_false < 0)
        return NULL;

      if (is_false)
        return PyObject_Call(self->object_proxy.wrapped, args, kwds);
    }
    else
    {
      int is_false = PyObject_Not(self->enabled);

      if (is_false < 0)
        return NULL;

      if (is_false)
        return PyObject_Call(self->object_proxy.wrapped, args, kwds);
    }
  }

  /*
   * We need to do things different depending on whether we are likely
   * wrapping an instance method vs a static method or class method.
   */

  matched =
      PyUnicode_CompareWithASCIIString(self->binding, "function") == 0 ||
      PyUnicode_CompareWithASCIIString(self->binding, "callable") == 0;

  if (matched)
  {

    if (self->instance == Py_None && PyTuple_GET_SIZE(args) != 0)
    {
      /*
       * This situation can occur where someone is calling the
       * instancemethod via the class type and passing the
       * instance as the first argument. We need to shift the args
       * before making the call to the wrapper and effectively
       * bind the instance to the wrapped function using a partial
       * so the wrapper doesn't see anything as being different.
       */

      instance = PyTuple_GetItem(args, 0);

      if (!instance)
        return NULL;

      int check = PyObject_IsInstance(instance, self->owner);

      if (check < 0)
        return NULL;

      if (check)
      {
        wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
        if (!state)
          return NULL;
        wrapped = PyObject_CallFunctionObjArgs(
            (PyObject *)state->PartialCallableObjectProxy_Type,
            self->object_proxy.wrapped, instance, NULL);

        if (!wrapped)
          return NULL;

        param_args = PyTuple_GetSlice(args, 1, PyTuple_GET_SIZE(args));

        if (!param_args)
        {
          Py_DECREF(wrapped);
          return NULL;
        }

        args = param_args;
      }
      else
      {
        instance = self->instance;
      }
    }
    else
    {
      instance = self->instance;
    }

    if (!wrapped)
    {
      Py_INCREF(self->object_proxy.wrapped);
      wrapped = self->object_proxy.wrapped;
    }

    if (!kwds)
    {
      param_kwds = PyDict_New();
      if (!param_kwds)
      {
        Py_XDECREF(param_args);
        Py_DECREF(wrapped);
        return NULL;
      }
      kwds = param_kwds;
    }

    result = PyObject_CallFunctionObjArgs(self->wrapper, wrapped, instance,
                                          args, kwds, NULL);

    Py_XDECREF(param_args);
    Py_XDECREF(param_kwds);
    Py_DECREF(wrapped);

    return result;
  }
  else
  {
    /*
     * As in this case we would be dealing with a classmethod or
     * staticmethod, then _self_instance will only tell us whether
     * when calling the classmethod or staticmethod they did it via
     * an instance of the class it is bound to and not the case
     * where done by the class type itself. We thus ignore
     * _self_instance and use the __self__ attribute of the bound
     * function instead. For a classmethod, this means instance will
     * be the class type and for a staticmethod it will be None.
     * This is probably the more useful thing we can pass through
     * even though we loose knowledge of whether they were called on
     * the instance vs the class type, as it reflects what they have
     * available in the decoratored function.
     */

    instance = PyObject_GetAttr(self->object_proxy.wrapped, state->str_self);

    if (!instance)
    {
      if (!PyErr_ExceptionMatches(PyExc_AttributeError))
        return NULL;
      PyErr_Clear();
      Py_INCREF(Py_None);
      instance = Py_None;
    }

    if (!kwds)
    {
      param_kwds = PyDict_New();
      if (!param_kwds)
      {
        Py_DECREF(instance);
        return NULL;
      }
      kwds = param_kwds;
    }

    result = PyObject_CallFunctionObjArgs(
        self->wrapper, self->object_proxy.wrapped, instance, args, kwds, NULL);

    Py_XDECREF(param_kwds);

    Py_DECREF(instance);

    return result;
  }
}

/* ------------------------------------------------------------------------- */

static PyObject *WraptBoundFunctionWrapper_getattr(
    WraptFunctionWrapperObject *self, PyObject *args)
{
  PyObject *name = NULL;

  if (!PyArg_ParseTuple(args, "U:__getattr__", &name))
    return NULL;

  if (self->parent && self->parent != Py_None)
  {
    PyObject *result = PyObject_GetAttr(self->parent, name);

    if (result)
      return result;

    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
      return NULL;

    PyErr_Clear();
  }

  return WraptObjectProxy_getattr((WraptObjectProxyObject *)self, args);
}

/* ------------------------------------------------------------------------- */

static int WraptBoundFunctionWrapper_setattro(
    WraptFunctionWrapperObject *self, PyObject *name, PyObject *value)
{
  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  int match = PyUnicode_Tailmatch(name, state->str_self_, 0, PY_SSIZE_T_MAX, -1);

  if (match < 0)
    return -1;

  if (match)
  {
    /* Check if this is an internal slot (has a getset descriptor on the type) */
    if (wrapt_type_has_attr(Py_TYPE(self), name))
      return PyObject_GenericSetAttr((PyObject *)self, name, value);

    /* User-defined _self_ attribute — delegate to parent */
    if (self->parent && self->parent != Py_None)
      return PyObject_GenericSetAttr(self->parent, name, value);
  }

  /* Fall through to base ObjectProxy setattro for everything else */
  return WraptObjectProxy_setattro((WraptObjectProxyObject *)self, name, value);
}

/* ------------------------------------------------------------------------- */

static PyMethodDef WraptBoundFunctionWrapper_methods[] = {
    {"__getattr__", (PyCFunction)WraptBoundFunctionWrapper_getattr, METH_VARARGS,
     0},
    {NULL, NULL},
};

/* ------------------------------------------------------------------------- */

static PyType_Slot WraptBoundFunctionWrapper_slots[] = {
    {Py_tp_dealloc, WraptFunctionWrapperBase_dealloc},
    {Py_tp_traverse, WraptFunctionWrapperBase_traverse},
    {Py_tp_clear, WraptFunctionWrapperBase_clear},
    {Py_tp_call, WraptBoundFunctionWrapper_call},
    {Py_tp_setattro, WraptBoundFunctionWrapper_setattro},
    {Py_tp_methods, WraptBoundFunctionWrapper_methods},
    {0, NULL},
};

static PyType_Spec WraptBoundFunctionWrapper_spec = {
    .name = "_wrappers.BoundFunctionWrapper",
    .basicsize = sizeof(WraptFunctionWrapperObject),
    .itemsize = 0,
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = WraptBoundFunctionWrapper_slots,
};

/* ------------------------------------------------------------------------- */

static int WraptFunctionWrapper_init(WraptFunctionWrapperObject *self,
                                     PyObject *args, PyObject *kwds)
{
  PyObject *wrapped = NULL;
  PyObject *wrapper = NULL;
  PyObject *enabled = Py_None;
  PyObject *binding = NULL;
  PyObject *binding_owned = NULL;
  PyObject *instance = NULL;

  wrapt_module_state *state = wrapt_state_from_type(Py_TYPE(self));
  if (!state)
    return -1;

  PyObject *function_str = state->str_function;
  PyObject *classmethod_str = state->str_classmethod;
  PyObject *staticmethod_str = state->str_staticmethod;
  PyObject *callable_str = state->str_callable;
  PyObject *builtin_str = state->str_builtin;
  PyObject *class_str = state->str_class;
  PyObject *instancemethod_str = state->str_instancemethod;

  int result = 0;

  char *const kwlist[] = {"wrapped", "wrapper", "enabled", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|O:FunctionWrapper", kwlist,
                                   &wrapped, &wrapper, &enabled))
  {
    return -1;
  }

  if (PyObject_TypeCheck(wrapped, state->FunctionWrapperBase_Type))
  {
    binding_owned = PyObject_GetAttr(wrapped, state->str_self_binding);
    if (!binding_owned)
    {
      if (!PyErr_ExceptionMatches(PyExc_AttributeError))
        return -1;
      PyErr_Clear();
    }
    else
    {
      binding = binding_owned;
    }
  }

  if (!binding)
  {
    int check;

    if (PyCFunction_Check(wrapped))
    {
      binding = builtin_str;
    }
    else
    {
      check = PyObject_IsInstance(wrapped, (PyObject *)&PyFunction_Type);
      if (check < 0)
        goto error;
      if (check)
        binding = function_str;
    }

    if (!binding)
    {
      check = PyObject_IsInstance(wrapped, (PyObject *)&PyType_Type);
      if (check < 0)
        goto error;
      if (check)
        binding = class_str;
    }

    if (!binding)
    {
      check = PyObject_IsInstance(wrapped, (PyObject *)&PyClassMethod_Type);
      if (check < 0)
        goto error;
      if (check)
        binding = classmethod_str;
    }

    if (!binding)
    {
      check = PyObject_IsInstance(wrapped, (PyObject *)&PyStaticMethod_Type);
      if (check < 0)
        goto error;
      if (check)
        binding = staticmethod_str;
    }

    if (!binding)
    {
      instance = PyObject_GetAttr(wrapped, state->str_self);
      if (instance != NULL)
      {
        check = PyObject_IsInstance(instance, (PyObject *)&PyType_Type);
        if (check < 0)
        {
          Py_DECREF(instance);
          goto error;
        }
        if (check)
        {
          binding = classmethod_str;
        }
        else
        {
          check = PyObject_IsInstance(wrapped, (PyObject *)&PyMethod_Type);
          if (check < 0)
          {
            Py_DECREF(instance);
            goto error;
          }
          if (check)
            binding = instancemethod_str;
          else
            binding = callable_str;
        }

        Py_DECREF(instance);
      }
      else
      {
        if (!PyErr_ExceptionMatches(PyExc_AttributeError))
          goto error;
        PyErr_Clear();

        binding = callable_str;
      }
    }
  }

  result = WraptFunctionWrapperBase_raw_init(
      self, wrapped, Py_None, wrapper, enabled, binding, Py_None, Py_None);

  Py_XDECREF(binding_owned);

  return result;

error:
  Py_XDECREF(binding_owned);
  return -1;
}

/* ------------------------------------------------------------------------- */

static PyType_Slot WraptFunctionWrapper_slots[] = {
    {Py_tp_dealloc, WraptFunctionWrapperBase_dealloc},
    {Py_tp_traverse, WraptFunctionWrapperBase_traverse},
    {Py_tp_clear, WraptFunctionWrapperBase_clear},
    {Py_tp_init, WraptFunctionWrapper_init},
    {0, NULL},
};

static PyType_Spec WraptFunctionWrapper_spec = {
    .name = "_wrappers.FunctionWrapper",
    .basicsize = sizeof(WraptFunctionWrapperObject),
    .itemsize = 0,
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = WraptFunctionWrapper_slots,
};

/* ------------------------------------------------------------------------- */

/* PyModule_AddObjectRef polyfill for Python 3.9. */
#if PY_VERSION_HEX >= 0x030A0000
#define WRAPT_ADD_TYPE(mod, name, type)                               \
  do                                                                  \
  {                                                                   \
    if (PyModule_AddObjectRef((mod), (name), (PyObject *)(type)) < 0) \
      return -1;                                                      \
  } while (0)
#else
#define WRAPT_ADD_TYPE(mod, name, type)                            \
  do                                                               \
  {                                                                \
    Py_INCREF((PyObject *)(type));                                 \
    if (PyModule_AddObject((mod), (name), (PyObject *)(type)) < 0) \
    {                                                              \
      Py_DECREF((PyObject *)(type));                               \
      return -1;                                                   \
    }                                                              \
  } while (0)
#endif

/* Helper to create a heap type via PyType_FromModuleAndSpecWithBases and
 * store it on the module. Returns 0 on success, -1 on failure. */

static int wrapt_create_type(PyObject *module, PyTypeObject **out,
                             PyType_Spec *spec, PyObject *bases,
                             const char *name)
{
  PyObject *type;

  type = PyType_FromModuleAndSpec(module, spec, bases);
  if (type == NULL)
    return -1;

  /* PyType_FromModuleAndSpec returns a new reference. Add the type to
   * the module dict first (WRAPT_ADD_TYPE adds another ref via
   * PyModule_AddObjectRef on 3.10+, or via manual incref +
   * PyModule_AddObject on 3.9). Only then store the owned reference
   * in module state so that *out is never populated on failure. */

  WRAPT_ADD_TYPE(module, name, type);

  *out = (PyTypeObject *)type;

  return 0;
}

/* Helper: intern a C string into the given module-state slot. Returns 0
 * on success, -1 on failure (with PyErr set). */

static int wrapt_intern_string(PyObject **slot, const char *s)
{
  PyObject *o = PyUnicode_InternFromString(s);
  if (o == NULL)
    return -1;
  *slot = o;
  return 0;
}

static int wrapt_init_strings(wrapt_module_state *state)
{
  if (wrapt_intern_string(&state->str_wrapped, "__wrapped__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_wrapped_factory,
                          "__wrapped_factory__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_wrapped_get, "__wrapped_get__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_module, "__module__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_doc, "__doc__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_setattr_fixups,
                          "__wrapped_setattr_fixups__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_getattr, "__getattr__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_self_, "_self_") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_callable, "callable") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_bound_function_wrapper,
                          "__bound_function_wrapper__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_function, "function") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_classmethod, "classmethod") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_staticmethod, "staticmethod") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_builtin, "builtin") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_class, "class") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_instancemethod, "instancemethod") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_object_proxy, "__object_proxy__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_enter, "__enter__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_exit, "__exit__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_aenter, "__aenter__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_aexit, "__aexit__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_mro_entries, "__mro_entries__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_name, "__name__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_qualname, "__qualname__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_annotations, "__annotations__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_dunder_class, "__class__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_self, "__self__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_set_name, "__set_name__") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_self_binding, "_self_binding") < 0)
    return -1;
  if (wrapt_intern_string(&state->str_dict, "__dict__") < 0)
    return -1;
  return 0;
}

static int wrapt_exec(PyObject *module)
{
  wrapt_module_state *state = wrapt_get_state(module);
  PyObject *bases = NULL;

  /* Eagerly intern all cached strings. m_clear releases anything that
   * was set before a failure, so we can return -1 on the first error. */
  if (wrapt_init_strings(state) < 0)
    return -1;

  /* ObjectProxy: base = object (default). */
  if (wrapt_create_type(module, &state->ObjectProxy_Type,
                        &WraptObjectProxy_spec, NULL, "ObjectProxy") < 0)
    return -1;

  /* CallableObjectProxy: base = ObjectProxy. */
  bases = PyTuple_Pack(1, (PyObject *)state->ObjectProxy_Type);
  if (!bases)
    return -1;
  if (wrapt_create_type(module, &state->CallableObjectProxy_Type,
                        &WraptCallableObjectProxy_spec, bases,
                        "CallableObjectProxy") < 0)
  {
    Py_DECREF(bases);
    return -1;
  }
  Py_DECREF(bases);

  /* PartialCallableObjectProxy: base = ObjectProxy. */
  bases = PyTuple_Pack(1, (PyObject *)state->ObjectProxy_Type);
  if (!bases)
    return -1;
  if (wrapt_create_type(module, &state->PartialCallableObjectProxy_Type,
                        &WraptPartialCallableObjectProxy_spec, bases,
                        "PartialCallableObjectProxy") < 0)
  {
    Py_DECREF(bases);
    return -1;
  }
  Py_DECREF(bases);

  /* _FunctionWrapperBase: base = ObjectProxy. */
  bases = PyTuple_Pack(1, (PyObject *)state->ObjectProxy_Type);
  if (!bases)
    return -1;
  if (wrapt_create_type(module, &state->FunctionWrapperBase_Type,
                        &WraptFunctionWrapperBase_spec, bases,
                        "_FunctionWrapperBase") < 0)
  {
    Py_DECREF(bases);
    return -1;
  }
  Py_DECREF(bases);

  /* BoundFunctionWrapper: base = _FunctionWrapperBase. */
  bases = PyTuple_Pack(1, (PyObject *)state->FunctionWrapperBase_Type);
  if (!bases)
    return -1;
  if (wrapt_create_type(module, &state->BoundFunctionWrapper_Type,
                        &WraptBoundFunctionWrapper_spec, bases,
                        "BoundFunctionWrapper") < 0)
  {
    Py_DECREF(bases);
    return -1;
  }
  Py_DECREF(bases);

  /* FunctionWrapper: base = _FunctionWrapperBase. */
  bases = PyTuple_Pack(1, (PyObject *)state->FunctionWrapperBase_Type);
  if (!bases)
    return -1;
  if (wrapt_create_type(module, &state->FunctionWrapper_Type,
                        &WraptFunctionWrapper_spec, bases,
                        "FunctionWrapper") < 0)
  {
    Py_DECREF(bases);
    return -1;
  }
  Py_DECREF(bases);

  /* Cache WrapperNotInitializedError from wrapt.wrappers. The module is
   * already in sys.modules because __wrapt__.py imports it before us. */

  {
    PyObject *wrapt_wrappers_module = PyImport_ImportModule("wrapt.wrappers");
    if (!wrapt_wrappers_module)
      return -1;

    state->WrapperNotInitializedError = PyObject_GetAttrString(
        wrapt_wrappers_module, "WrapperNotInitializedError");

    Py_DECREF(wrapt_wrappers_module);

    if (!state->WrapperNotInitializedError)
      return -1;
  }

  return 0;
}

static int wrapt_traverse(PyObject *module, visitproc visit, void *arg)
{
  wrapt_module_state *state = wrapt_get_state(module);
  if (state == NULL)
    return 0;
  Py_VISIT(state->ObjectProxy_Type);
  Py_VISIT(state->CallableObjectProxy_Type);
  Py_VISIT(state->PartialCallableObjectProxy_Type);
  Py_VISIT(state->FunctionWrapperBase_Type);
  Py_VISIT(state->BoundFunctionWrapper_Type);
  Py_VISIT(state->FunctionWrapper_Type);
  Py_VISIT(state->str_wrapped);
  Py_VISIT(state->str_wrapped_factory);
  Py_VISIT(state->str_wrapped_get);
  Py_VISIT(state->str_module);
  Py_VISIT(state->str_doc);
  Py_VISIT(state->str_setattr_fixups);
  Py_VISIT(state->str_getattr);
  Py_VISIT(state->str_self_);
  Py_VISIT(state->str_callable);
  Py_VISIT(state->str_bound_function_wrapper);
  Py_VISIT(state->str_function);
  Py_VISIT(state->str_classmethod);
  Py_VISIT(state->str_staticmethod);
  Py_VISIT(state->str_builtin);
  Py_VISIT(state->str_class);
  Py_VISIT(state->str_instancemethod);
  Py_VISIT(state->str_object_proxy);
  Py_VISIT(state->str_enter);
  Py_VISIT(state->str_exit);
  Py_VISIT(state->str_aenter);
  Py_VISIT(state->str_aexit);
  Py_VISIT(state->str_mro_entries);
  Py_VISIT(state->str_name);
  Py_VISIT(state->str_qualname);
  Py_VISIT(state->str_annotations);
  Py_VISIT(state->str_dunder_class);
  Py_VISIT(state->str_self);
  Py_VISIT(state->str_set_name);
  Py_VISIT(state->str_self_binding);
  Py_VISIT(state->str_dict);
  Py_VISIT(state->WrapperNotInitializedError);
  return 0;
}

static int wrapt_clear(PyObject *module)
{
  wrapt_module_state *state = wrapt_get_state(module);
  if (state == NULL)
    return 0;
  Py_CLEAR(state->ObjectProxy_Type);
  Py_CLEAR(state->CallableObjectProxy_Type);
  Py_CLEAR(state->PartialCallableObjectProxy_Type);
  Py_CLEAR(state->FunctionWrapperBase_Type);
  Py_CLEAR(state->BoundFunctionWrapper_Type);
  Py_CLEAR(state->FunctionWrapper_Type);
  Py_CLEAR(state->str_wrapped);
  Py_CLEAR(state->str_wrapped_factory);
  Py_CLEAR(state->str_wrapped_get);
  Py_CLEAR(state->str_module);
  Py_CLEAR(state->str_doc);
  Py_CLEAR(state->str_setattr_fixups);
  Py_CLEAR(state->str_getattr);
  Py_CLEAR(state->str_self_);
  Py_CLEAR(state->str_callable);
  Py_CLEAR(state->str_bound_function_wrapper);
  Py_CLEAR(state->str_function);
  Py_CLEAR(state->str_classmethod);
  Py_CLEAR(state->str_staticmethod);
  Py_CLEAR(state->str_builtin);
  Py_CLEAR(state->str_class);
  Py_CLEAR(state->str_instancemethod);
  Py_CLEAR(state->str_object_proxy);
  Py_CLEAR(state->str_enter);
  Py_CLEAR(state->str_exit);
  Py_CLEAR(state->str_aenter);
  Py_CLEAR(state->str_aexit);
  Py_CLEAR(state->str_mro_entries);
  Py_CLEAR(state->str_name);
  Py_CLEAR(state->str_qualname);
  Py_CLEAR(state->str_annotations);
  Py_CLEAR(state->str_dunder_class);
  Py_CLEAR(state->str_self);
  Py_CLEAR(state->str_set_name);
  Py_CLEAR(state->str_self_binding);
  Py_CLEAR(state->str_dict);
  Py_CLEAR(state->WrapperNotInitializedError);
  return 0;
}

static void wrapt_free(void *module)
{
  (void)wrapt_clear((PyObject *)module);
}

static PyModuleDef_Slot wrapt_slots[] = {
    {Py_mod_exec, wrapt_exec},
#if PY_VERSION_HEX >= 0x030C0000
    /* Per-interpreter GIL support (PEP 684, 3.12+). All process-wide
     * mutable state has been removed: types are heap types created per
     * interpreter via PyType_FromModuleAndSpec, cached strings live in
     * per-interpreter module state, and there are no static PyObject *
     * globals. */
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if PY_VERSION_HEX >= 0x030D0000
    /* Free-threading support (3.13+). The previous Py_MOD_GIL_NOT_USED
     * declaration via PyUnstable_Module_SetGIL was unsound because the
     * 19 cached interned strings were lazily initialized via racy
     * static locals. With eager init in wrapt_exec the declaration is
     * now actually correct. */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_wrappers",
    .m_doc = NULL,
    .m_size = sizeof(wrapt_module_state),
    .m_methods = NULL,
    .m_slots = wrapt_slots,
    .m_traverse = wrapt_traverse,
    .m_clear = wrapt_clear,
    .m_free = wrapt_free,
};

PyMODINIT_FUNC PyInit__wrappers(void)
{
  return PyModuleDef_Init(&moduledef);
}

/* ------------------------------------------------------------------------- */
