/*
 * General object operations and protocol implementations,
 * including their specialisations for certain builtins.
 *
 * Optional optimisations for builtins are in Optimize.c.
 *
 * Required replacements of builtins are in Builtins.c.
 */

/////////////// RaiseNoneIterError.proto ///////////////

static CYTHON_INLINE void __Pyx_RaiseNoneNotIterableError(void);

/////////////// RaiseNoneIterError ///////////////

static CYTHON_INLINE void __Pyx_RaiseNoneNotIterableError(void) {
    PyErr_SetString(PyExc_TypeError, "'NoneType' object is not iterable");
}

/////////////// RaiseTooManyValuesToUnpack.proto ///////////////

static CYTHON_INLINE void __Pyx_RaiseTooManyValuesError(Py_ssize_t expected);

/////////////// RaiseTooManyValuesToUnpack ///////////////

static CYTHON_INLINE void __Pyx_RaiseTooManyValuesError(Py_ssize_t expected) {
    PyErr_Format(PyExc_ValueError,
                 "too many values to unpack (expected %" CYTHON_FORMAT_SSIZE_T "d)", expected);
}

/////////////// RaiseNeedMoreValuesToUnpack.proto ///////////////

static CYTHON_INLINE void __Pyx_RaiseNeedMoreValuesError(Py_ssize_t index);

/////////////// RaiseNeedMoreValuesToUnpack ///////////////

static CYTHON_INLINE void __Pyx_RaiseNeedMoreValuesError(Py_ssize_t index) {
    PyErr_Format(PyExc_ValueError,
                 "need more than %" CYTHON_FORMAT_SSIZE_T "d value%.1s to unpack",
                 index, (index == 1) ? "" : "s");
}

/////////////// UnpackTupleError.proto ///////////////

static void __Pyx_UnpackTupleError(PyObject *, Py_ssize_t index); /*proto*/

/////////////// UnpackTupleError ///////////////
//@requires: RaiseNoneIterError
//@requires: RaiseNeedMoreValuesToUnpack
//@requires: RaiseTooManyValuesToUnpack

static void __Pyx_UnpackTupleError(PyObject *t, Py_ssize_t index) {
    if (t == Py_None) {
      __Pyx_RaiseNoneNotIterableError();
    } else {
      Py_ssize_t size = __Pyx_PyTuple_GET_SIZE(t);
 #if !CYTHON_ASSUME_SAFE_SIZE
      if (unlikely(size < 0)) return;
 #endif
      if (size < index) {
        __Pyx_RaiseNeedMoreValuesError(size);
      } else {
        __Pyx_RaiseTooManyValuesError(index);
      }
    }
}

/////////////// UnpackItemEndCheck.proto ///////////////

static int __Pyx_IternextUnpackEndCheck(PyObject *retval, Py_ssize_t expected); /*proto*/

/////////////// UnpackItemEndCheck ///////////////
//@requires: RaiseTooManyValuesToUnpack
//@requires: IterFinish

static int __Pyx_IternextUnpackEndCheck(PyObject *retval, Py_ssize_t expected) {
    if (unlikely(retval)) {
        Py_DECREF(retval);
        __Pyx_RaiseTooManyValuesError(expected);
        return -1;
    }

    return __Pyx_IterFinish();
}

/////////////// UnpackTuple2.proto ///////////////

static CYTHON_INLINE int __Pyx_unpack_tuple2(
    PyObject* tuple, PyObject** value1, PyObject** value2, int is_tuple, int has_known_size, int decref_tuple);

static CYTHON_INLINE int __Pyx_unpack_tuple2_exact(
    PyObject* tuple, PyObject** value1, PyObject** value2, int decref_tuple);
static int __Pyx_unpack_tuple2_generic(
    PyObject* tuple, PyObject** value1, PyObject** value2, int has_known_size, int decref_tuple);

/////////////// UnpackTuple2 ///////////////
//@requires: UnpackItemEndCheck
//@requires: UnpackTupleError
//@requires: RaiseNeedMoreValuesToUnpack

static CYTHON_INLINE int __Pyx_unpack_tuple2(
        PyObject* tuple, PyObject** value1, PyObject** value2, int is_tuple, int has_known_size, int decref_tuple) {
    if (likely(is_tuple || PyTuple_Check(tuple))) {
        Py_ssize_t size;
        if (has_known_size) {
            return __Pyx_unpack_tuple2_exact(tuple, value1, value2, decref_tuple);
        }
        size = __Pyx_PyTuple_GET_SIZE(tuple);
        if (likely(size == 2)) {
            return __Pyx_unpack_tuple2_exact(tuple, value1, value2, decref_tuple);
        }
        if (size >= 0) {
            // "size == -1" indicates an error already.
            __Pyx_UnpackTupleError(tuple, 2);
        }
        return -1;
    } else {
        return __Pyx_unpack_tuple2_generic(tuple, value1, value2, has_known_size, decref_tuple);
    }
}

static CYTHON_INLINE int __Pyx_unpack_tuple2_exact(
        PyObject* tuple, PyObject** pvalue1, PyObject** pvalue2, int decref_tuple) {
    PyObject *value1 = NULL, *value2 = NULL;
#if CYTHON_AVOID_BORROWED_REFS || !CYTHON_ASSUME_SAFE_MACROS
    value1 = __Pyx_PySequence_ITEM(tuple, 0);  if (unlikely(!value1)) goto bad;
    value2 = __Pyx_PySequence_ITEM(tuple, 1);  if (unlikely(!value2)) goto bad;
#else
    value1 = PyTuple_GET_ITEM(tuple, 0);  Py_INCREF(value1);
    value2 = PyTuple_GET_ITEM(tuple, 1);  Py_INCREF(value2);
#endif
    if (decref_tuple) {
        Py_DECREF(tuple);
    }

    *pvalue1 = value1;
    *pvalue2 = value2;
    return 0;
#if CYTHON_AVOID_BORROWED_REFS || !CYTHON_ASSUME_SAFE_MACROS
bad:
    Py_XDECREF(value1);
    Py_XDECREF(value2);
    if (decref_tuple) { Py_XDECREF(tuple); }
    return -1;
#endif
}

static int __Pyx_unpack_tuple2_generic(PyObject* tuple, PyObject** pvalue1, PyObject** pvalue2,
                                       int has_known_size, int decref_tuple) {
    Py_ssize_t index;
    PyObject *value1 = NULL, *value2 = NULL, *iter = NULL;
    iternextfunc iternext;

    iter = PyObject_GetIter(tuple);
    if (unlikely(!iter)) goto bad;
    if (decref_tuple) { Py_DECREF(tuple); tuple = NULL; }

    iternext = __Pyx_PyObject_GetIterNextFunc(iter);
    value1 = iternext(iter); if (unlikely(!value1)) { index = 0; goto unpacking_failed; }
    value2 = iternext(iter); if (unlikely(!value2)) { index = 1; goto unpacking_failed; }
    if (!has_known_size && unlikely(__Pyx_IternextUnpackEndCheck(iternext(iter), 2))) goto bad;

    Py_DECREF(iter);
    *pvalue1 = value1;
    *pvalue2 = value2;
    return 0;

unpacking_failed:
    if (!has_known_size && __Pyx_IterFinish() == 0)
        __Pyx_RaiseNeedMoreValuesError(index);
bad:
    Py_XDECREF(iter);
    Py_XDECREF(value1);
    Py_XDECREF(value2);
    if (decref_tuple) { Py_XDECREF(tuple); }
    return -1;
}


/////////////// IterNextPlain.proto ///////////////

static CYTHON_INLINE PyObject *__Pyx_PyIter_Next_Plain(PyObject *iterator);
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
static PyObject *__Pyx_GetBuiltinNext_LimitedAPI(void);
#endif

/////////////// IterNextPlain.module_state_decls ///////////////

#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
PyObject *__Pyx_GetBuiltinNext_LimitedAPI_cache;
#endif

/////////////// IterNextPlain ///////////////
//@requires: GetBuiltinName

// There is no replacement for "Py_TYPE(obj)->iternext" in the C-API.
// PyIter_Next() discards the StopIteration, unlike Python's "next()".
// PyObject_GetSlot() rejects non-heap types in CPython < 3.10 (and PyPy).

#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
static PyObject *__Pyx_GetBuiltinNext_LimitedAPI(void) {
    if (unlikely(!CGLOBAL(__Pyx_GetBuiltinNext_LimitedAPI_cache)))
        CGLOBAL(__Pyx_GetBuiltinNext_LimitedAPI_cache) = __Pyx_GetBuiltinName(PYIDENT("next"));
    // Returns the globally owned reference, not a new reference!
    return CGLOBAL(__Pyx_GetBuiltinNext_LimitedAPI_cache);
}
#endif

static CYTHON_INLINE PyObject *__Pyx_PyIter_Next_Plain(PyObject *iterator) {
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
    PyObject *result;
    PyObject *next = __Pyx_GetBuiltinNext_LimitedAPI();
    if (unlikely(!next)) return NULL;
    result = PyObject_CallFunctionObjArgs(next, iterator, NULL);
    return result;
#else
    (void)__Pyx_GetBuiltinName; // only for early limited API
    iternextfunc iternext = __Pyx_PyObject_GetIterNextFunc(iterator);
    assert(iternext);
    return iternext(iterator);
#endif
}


/////////////// IterNext.proto ///////////////

#define __Pyx_PyIter_Next(obj) __Pyx_PyIter_Next2(obj, NULL)
static CYTHON_INLINE PyObject *__Pyx_PyIter_Next2(PyObject *, PyObject *); /*proto*/

/////////////// IterNext ///////////////
//@requires: Exceptions.c::PyThreadStateGet
//@requires: Exceptions.c::PyErrFetchRestore
//@requires: GetBuiltinName
//@requires: IterNextPlain

static PyObject *__Pyx_PyIter_Next2Default(PyObject* defval) {
    PyObject* exc_type;
    __Pyx_PyThreadState_declare
    __Pyx_PyThreadState_assign
    exc_type = __Pyx_PyErr_CurrentExceptionType();
    if (unlikely(exc_type)) {
        if (!defval || unlikely(!__Pyx_PyErr_GivenExceptionMatches(exc_type, PyExc_StopIteration)))
            return NULL;
        __Pyx_PyErr_Clear();
        Py_INCREF(defval);
        return defval;
    }
    if (defval) {
        Py_INCREF(defval);
        return defval;
    }
    __Pyx_PyErr_SetNone(PyExc_StopIteration);
    return NULL;
}

static void __Pyx_PyIter_Next_ErrorNoIterator(PyObject *iterator) {
    __Pyx_TypeName iterator_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(iterator));
    PyErr_Format(PyExc_TypeError,
        __Pyx_FMT_TYPENAME " object is not an iterator", iterator_type_name);
    __Pyx_DECREF_TypeName(iterator_type_name);
}

// originally copied from Py3's builtin_next()
static CYTHON_INLINE PyObject *__Pyx_PyIter_Next2(PyObject* iterator, PyObject* defval) {
    PyObject* next;
#if !CYTHON_COMPILING_IN_LIMITED_API
    // We always do a quick slot check because calling PyIter_Check() is so wasteful.
    iternextfunc iternext = __Pyx_PyObject_TryGetSlot(iterator, tp_iternext, iternextfunc);
    if (likely(iternext)) {
        next = iternext(iterator);
        if (likely(next))
            return next;
    #if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX < 0x030d0000
        // If the reason for failure was "not implemented", we're done raising the appropriate exception.
        if (unlikely(iternext == &_PyObject_NextNotImplemented))
            return NULL;
    #endif
    } else if (CYTHON_USE_TYPE_SLOTS) {
        // If CYTHON_USE_TYPE_SLOTS, then the slot was really not set and we don't have an iterable.
        __Pyx_PyIter_Next_ErrorNoIterator(iterator);
        return NULL;
    } else
    // Without CYTHON_USE_TYPE_SLOTS, don't trust "tp_iternext" and rely on PyIter_Check().
#endif
    if (unlikely(!PyIter_Check(iterator))) {
        __Pyx_PyIter_Next_ErrorNoIterator(iterator);
        return NULL;
    } else {
        // We'll probably only end up here in the Limited API, where we'd call Python's "next()" now.
        // If we have a default value, we don't need the StopIteration and can use PyIter_Next().
        next = defval ? PyIter_Next(iterator) : __Pyx_PyIter_Next_Plain(iterator);
        if (likely(next))
            return next;
    }
    return __Pyx_PyIter_Next2Default(defval);
}


/////////////// IterFinish.proto ///////////////

static CYTHON_INLINE int __Pyx_IterFinish(void); /*proto*/

/////////////// IterFinish ///////////////
//@requires: Exceptions.c::PyThreadStateGet
//@requires: Exceptions.c::PyErrFetchRestore

// When PyIter_Next(iter) has returned NULL in order to signal termination,
// this function does the right cleanup and returns 0 on success.  If it
// detects an error that occurred in the iterator, it returns -1.

static CYTHON_INLINE int __Pyx_IterFinish(void) {
    PyObject* exc_type;
    __Pyx_PyThreadState_declare
    __Pyx_PyThreadState_assign
    exc_type = __Pyx_PyErr_CurrentExceptionType();
    if (unlikely(exc_type)) {
        if (unlikely(!__Pyx_PyErr_GivenExceptionMatches(exc_type, PyExc_StopIteration)))
            return -1;
        __Pyx_PyErr_Clear();
        return 0;
    }
    return 0;
}


/////////////// ObjectGetItem.proto ///////////////

#if CYTHON_USE_TYPE_SLOTS
static CYTHON_INLINE PyObject *__Pyx_PyObject_GetItem(PyObject *obj, PyObject *key);/*proto*/
#else
#define __Pyx_PyObject_GetItem(obj, key)  PyObject_GetItem(obj, key)
#endif

/////////////// ObjectGetItem ///////////////
// //@requires: GetItemInt - added in IndexNode as it uses templating.
//@requires: PyObjectGetAttrStrNoError
//@requires: PyObjectCallOneArg

#if CYTHON_USE_TYPE_SLOTS
static PyObject *__Pyx_PyObject_GetIndex(PyObject *obj, PyObject *index) {
    // Get element from sequence object `obj` at index `index`.
    PyObject *runerr = NULL;
    Py_ssize_t key_value;
    key_value = __Pyx_PyIndex_AsSsize_t(index);
    if (likely(key_value != -1 || !(runerr = PyErr_Occurred()))) {
        return __Pyx_GetItemInt_Fast(obj, key_value, 0, 1, 1, 1);
    }

    // Error handling code -- only manage OverflowError differently.
    if (PyErr_GivenExceptionMatches(runerr, PyExc_OverflowError)) {
        __Pyx_TypeName index_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(index));
        PyErr_Clear();
        PyErr_Format(PyExc_IndexError,
            "cannot fit '" __Pyx_FMT_TYPENAME "' into an index-sized integer", index_type_name);
        __Pyx_DECREF_TypeName(index_type_name);
    }
    return NULL;
}

static PyObject *__Pyx_PyObject_GetItem_Slow(PyObject *obj, PyObject *key) {
    __Pyx_TypeName obj_type_name;
    // Handles less common slow-path checks for GetItem
    if (likely(PyType_Check(obj))) {
        PyObject *meth = __Pyx_PyObject_GetAttrStrNoError(obj, PYIDENT("__class_getitem__"));
        if (!meth) {
            PyErr_Clear();
        } else {
            PyObject *result = __Pyx_PyObject_CallOneArg(meth, key);
            Py_DECREF(meth);
            return result;
        }
    }

    obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
    PyErr_Format(PyExc_TypeError,
        "'" __Pyx_FMT_TYPENAME "' object is not subscriptable", obj_type_name);
    __Pyx_DECREF_TypeName(obj_type_name);
    return NULL;
}

static PyObject *__Pyx_PyObject_GetItem(PyObject *obj, PyObject *key) {
    PyTypeObject *tp = Py_TYPE(obj);
    PyMappingMethods *mm = tp->tp_as_mapping;
    PySequenceMethods *sm = tp->tp_as_sequence;

    if (likely(mm && mm->mp_subscript)) {
        return mm->mp_subscript(obj, key);
    }
    if (likely(sm && sm->sq_item)) {
        return __Pyx_PyObject_GetIndex(obj, key);
    }
    return __Pyx_PyObject_GetItem_Slow(obj, key);
}
#endif


/////////////// DictGetItem.proto ///////////////

#if !CYTHON_COMPILING_IN_PYPY
static PyObject *__Pyx_PyDict_GetItem(PyObject *d, PyObject* key);/*proto*/

#define __Pyx_PyObject_Dict_GetItem(obj, name) \
    (likely(PyDict_CheckExact(obj)) ? \
     __Pyx_PyDict_GetItem(obj, name) : PyObject_GetItem(obj, name))

#else
#define __Pyx_PyDict_GetItem(d, key) PyObject_GetItem(d, key)
#define __Pyx_PyObject_Dict_GetItem(obj, name)  PyObject_GetItem(obj, name)
#endif

/////////////// DictGetItem ///////////////

#if !CYTHON_COMPILING_IN_PYPY
static PyObject *__Pyx_PyDict_GetItem(PyObject *d, PyObject* key) {
    PyObject *value;
    if (unlikely(__Pyx_PyDict_GetItemRef(d, key, &value) == 0)) { // no value, no error
        if (unlikely(PyTuple_Check(key))) {
            // CPython interprets tuples as separate arguments => must wrap them in another tuple.
            PyObject* args = PyTuple_Pack(1, key);
            if (likely(args)) {
                PyErr_SetObject(PyExc_KeyError, args);
                Py_DECREF(args);
            }
        } else {
            // Avoid tuple packing if possible.
            PyErr_SetObject(PyExc_KeyError, key);
        }
    }
    return value;
}
#endif

/////////////// GetItemInt.proto ///////////////
//@substitute: tempita

#define __Pyx_GetItemInt(o, i, type, is_signed, to_py_func, is_list, wraparound, boundscheck, has_gil, unsafe_shared) \
    (__Pyx_fits_Py_ssize_t(i, type, is_signed) ? \
    __Pyx_GetItemInt_Fast(o, (Py_ssize_t)i, is_list, wraparound, boundscheck, unsafe_shared) : \
    (is_list ? (PyErr_SetString(PyExc_IndexError, "list index out of range"), (PyObject*)NULL) : \
               __Pyx_GetItemInt_Generic(o, to_py_func(i))))

{{for type in ['List', 'Tuple']}}
#define __Pyx_GetItemInt_{{type}}(o, i, type, is_signed, to_py_func, is_list, wraparound, boundscheck, has_gil, unsafe_shared) \
    (__Pyx_fits_Py_ssize_t(i, type, is_signed) ? \
    __Pyx_GetItemInt_{{type}}_Fast(o, (Py_ssize_t)i, wraparound, boundscheck, unsafe_shared) : \
    (PyErr_SetString(PyExc_IndexError, "{{ type.lower() }} index out of range"), (PyObject*)NULL))

static CYTHON_INLINE PyObject *__Pyx_GetItemInt_{{type}}_Fast(PyObject *o, Py_ssize_t i,
                                                              int wraparound, int boundscheck, int unsafe_shared);
{{endfor}}

static PyObject *__Pyx_GetItemInt_Generic(PyObject *o, PyObject* j);
static CYTHON_INLINE PyObject *__Pyx_GetItemInt_Fast(PyObject *o, Py_ssize_t i,
                                                     int is_list, int wraparound, int boundscheck, int unsafe_shared);

/////////////// GetItemInt ///////////////
//@substitute: tempita

static PyObject *__Pyx_GetItemInt_Generic(PyObject *o, PyObject* j) {
    PyObject *r;
    if (unlikely(!j)) return NULL;
    r = PyObject_GetItem(o, j);
    Py_DECREF(j);
    return r;
}

{{for type in ['List', 'Tuple']}}
static CYTHON_INLINE PyObject *__Pyx_GetItemInt_{{type}}_Fast(PyObject *o, Py_ssize_t i,
                                                              int wraparound, int boundscheck, int unsafe_shared) {
    CYTHON_MAYBE_UNUSED_VAR(unsafe_shared);
#if CYTHON_ASSUME_SAFE_SIZE{{if type == 'Tuple'}} && CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS{{endif}}
    Py_ssize_t wrapped_i = i;
    if (wraparound & unlikely(i < 0)) {
        wrapped_i += Py{{type}}_GET_SIZE(o);
    }
    {{if type == 'List'}}
    if ((CYTHON_AVOID_BORROWED_REFS || CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS || !CYTHON_ASSUME_SAFE_MACROS)) {
        // Note that there's a minor thread-safety issue where the size for wraparound may change before we use it.
        // CPython doesn't worry about it and serious violations will be caught by boundschecking anyway.
        return __Pyx_PyList_GetItemRefFast(o, wrapped_i, unsafe_shared);
    } else
    {{endif}}
    if ((!boundscheck) || likely(__Pyx_is_valid_index(wrapped_i, Py{{type}}_GET_SIZE(o)))) {
        return __Pyx_NewRef(Py{{type}}_GET_ITEM(o, wrapped_i));
    }
    return __Pyx_GetItemInt_Generic(o, PyLong_FromSsize_t(i));
#else
    (void)wraparound;
    (void)boundscheck;
    return PySequence_GetItem(o, i);
#endif
}
{{endfor}}

static CYTHON_INLINE PyObject *__Pyx_GetItemInt_Fast(PyObject *o, Py_ssize_t i, int is_list,
                                                     int wraparound, int boundscheck, int unsafe_shared) {
    CYTHON_MAYBE_UNUSED_VAR(unsafe_shared);
#if CYTHON_ASSUME_SAFE_MACROS && CYTHON_ASSUME_SAFE_SIZE
    if (is_list || PyList_CheckExact(o)) {
        // We specifically aren't worrying about thread-safety issues from the delay between getting the size
        // and using the size because Python itself doesn't either. If we have boundschecking on they'll
        // be caught eventually.
        Py_ssize_t n = ((!wraparound) | likely(i >= 0)) ? i : i + PyList_GET_SIZE(o);
        if ((CYTHON_AVOID_BORROWED_REFS || CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS)) {
            return __Pyx_PyList_GetItemRefFast(o, n, unsafe_shared);
        } else if ((!boundscheck) || (likely(__Pyx_is_valid_index(n, PyList_GET_SIZE(o))))) {
            return __Pyx_NewRef(PyList_GET_ITEM(o, n));
        }
    } else
    #if !CYTHON_AVOID_BORROWED_REFS
    if (PyTuple_CheckExact(o)) {
        Py_ssize_t n = ((!wraparound) | likely(i >= 0)) ? i : i + PyTuple_GET_SIZE(o);
        if ((!boundscheck) || likely(__Pyx_is_valid_index(n, PyTuple_GET_SIZE(o)))) {
            return __Pyx_NewRef(PyTuple_GET_ITEM(o, n));
        }
    } else
    #endif
#endif
#if CYTHON_USE_TYPE_SLOTS && !CYTHON_COMPILING_IN_PYPY
    {
        // inlined PySequence_GetItem() + special cased length overflow
        PyMappingMethods *mm = Py_TYPE(o)->tp_as_mapping;
        PySequenceMethods *sm = Py_TYPE(o)->tp_as_sequence;
        if (!is_list && mm && mm->mp_subscript) {
            PyObject *r, *key = PyLong_FromSsize_t(i);
            if (unlikely(!key)) return NULL;
            r = mm->mp_subscript(o, key);
            Py_DECREF(key);
            return r;
        }
        if (is_list || likely(sm && sm->sq_item)) {
            if (wraparound && unlikely(i < 0) && likely(sm->sq_length)) {
                Py_ssize_t l = sm->sq_length(o);
                if (likely(l >= 0)) {
                    i += l;
                } else {
                    // if length > max(Py_ssize_t), maybe the object can wrap around itself?
                    if (!PyErr_ExceptionMatches(PyExc_OverflowError))
                        return NULL;
                    PyErr_Clear();
                }
            }
            return sm->sq_item(o, i);
        }
    }
#else
    // PySequence_GetItem behaves differently to PyObject_GetItem for i<0
    // and possibly some other cases so can't generally be substituted
    if (is_list || !PyMapping_Check(o)) {
        return PySequence_GetItem(o, i);
    }
#endif
    (void)wraparound;
    (void)boundscheck;
    return __Pyx_GetItemInt_Generic(o, PyLong_FromSsize_t(i));
}

/////////////// SetItemInt.proto ///////////////

#define __Pyx_SetItemInt(o, i, v, type, is_signed, to_py_func, is_list, wraparound, boundscheck, has_gil, unsafe_shared) \
    (__Pyx_fits_Py_ssize_t(i, type, is_signed) ? \
    __Pyx_SetItemInt_Fast(o, (Py_ssize_t)i, v, is_list, wraparound, boundscheck, unsafe_shared) : \
    (is_list ? (PyErr_SetString(PyExc_IndexError, "list assignment index out of range"), -1) : \
               __Pyx_SetItemInt_Generic(o, to_py_func(i), v)))

static int __Pyx_SetItemInt_Generic(PyObject *o, PyObject *j, PyObject *v);
static CYTHON_INLINE int __Pyx_SetItemInt_Fast(PyObject *o, Py_ssize_t i, PyObject *v,
                                               int is_list, int wraparound, int boundscheck, int unsafe_shared);

/////////////// SetItemInt ///////////////

static int __Pyx_SetItemInt_Generic(PyObject *o, PyObject *j, PyObject *v) {
    int r;
    if (unlikely(!j)) return -1;
    r = PyObject_SetItem(o, j, v);
    Py_DECREF(j);
    return r;
}

static CYTHON_INLINE int __Pyx_SetItemInt_Fast(PyObject *o, Py_ssize_t i, PyObject *v, int is_list,
                                               int wraparound, int boundscheck, int unsafe_shared) {
    CYTHON_MAYBE_UNUSED_VAR(unsafe_shared);
#if CYTHON_ASSUME_SAFE_MACROS && CYTHON_ASSUME_SAFE_SIZE && !CYTHON_AVOID_BORROWED_REFS
    if (is_list || PyList_CheckExact(o)) {
        Py_ssize_t n = (!wraparound) ? i : ((likely(i >= 0)) ? i : i + PyList_GET_SIZE(o));
        if ((CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS && !__Pyx_IS_UNIQUELY_REFERENCED(o, unsafe_shared))) {
            Py_INCREF(v);
            return PyList_SetItem(o, n, v);
        } else if ((!boundscheck) || likely(__Pyx_is_valid_index(n, PyList_GET_SIZE(o)))) {
            PyObject* old;
            Py_INCREF(v);
            old = PyList_GET_ITEM(o, n);
            PyList_SET_ITEM(o, n, v);
            Py_DECREF(old);
            return 0;
        }
    } else
#endif
#if CYTHON_USE_TYPE_SLOTS && !CYTHON_COMPILING_IN_PYPY
    {
        // inlined PySequence_SetItem() + special cased length overflow
        PyMappingMethods *mm = Py_TYPE(o)->tp_as_mapping;
        PySequenceMethods *sm = Py_TYPE(o)->tp_as_sequence;
        if (!is_list && mm && mm->mp_ass_subscript) {
            int r;
            PyObject *key = PyLong_FromSsize_t(i);
            if (unlikely(!key)) return -1;
            r = mm->mp_ass_subscript(o, key, v);
            Py_DECREF(key);
            return r;
        }
        if (is_list || likely(sm && sm->sq_ass_item)) {
            if (wraparound && unlikely(i < 0) && likely(sm->sq_length)) {
                Py_ssize_t l = sm->sq_length(o);
                if (likely(l >= 0)) {
                    i += l;
                } else {
                    // if length > max(Py_ssize_t), maybe the object can wrap around itself?
                    if (!PyErr_ExceptionMatches(PyExc_OverflowError))
                        return -1;
                    PyErr_Clear();
                }
            }
            return sm->sq_ass_item(o, i, v);
        }
    }
#else
    // PySequence_SetItem behaves differently to PyObject_SetItem for i<0
    // and possibly some other cases, so can't generally be substituted.
    if (is_list || !PyMapping_Check(o)) {
        return PySequence_SetItem(o, i, v);
    }
#endif
    (void)wraparound;
    (void)boundscheck;
    return __Pyx_SetItemInt_Generic(o, PyLong_FromSsize_t(i), v);
}


/////////////// DelItemInt.proto ///////////////

#define __Pyx_DelItemInt(o, i, type, is_signed, to_py_func, is_list, wraparound, boundscheck, has_gil, unsafe_shared) \
    (__Pyx_fits_Py_ssize_t(i, type, is_signed) ? \
    __Pyx_DelItemInt_Fast(o, (Py_ssize_t)i, is_list, wraparound) : \
    (is_list ? (PyErr_SetString(PyExc_IndexError, "list assignment index out of range"), -1) : \
               __Pyx_DelItem_Generic(o, to_py_func(i))))

static int __Pyx_DelItem_Generic(PyObject *o, PyObject *j);
static CYTHON_INLINE int __Pyx_DelItemInt_Fast(PyObject *o, Py_ssize_t i,
                                               int is_list, int wraparound);

/////////////// DelItemInt ///////////////

static int __Pyx_DelItem_Generic(PyObject *o, PyObject *j) {
    int r;
    if (unlikely(!j)) return -1;
    r = PyObject_DelItem(o, j);
    Py_DECREF(j);
    return r;
}

static CYTHON_INLINE int __Pyx_DelItemInt_Fast(PyObject *o, Py_ssize_t i,
                                               int is_list, CYTHON_NCP_UNUSED int wraparound) {
#if !CYTHON_USE_TYPE_SLOTS
    // PySequence_DelItem behaves differently to PyObject_DelItem for i<0
    // and possibly some other cases so can't generally be substituted
    if (is_list || !PyMapping_Check(o)) {
        return PySequence_DelItem(o, i);
    }
#else
    // inlined PySequence_DelItem() + special cased length overflow
    PyMappingMethods *mm = Py_TYPE(o)->tp_as_mapping;
    PySequenceMethods *sm = Py_TYPE(o)->tp_as_sequence;
    if ((!is_list) && mm && mm->mp_ass_subscript) {
        PyObject *key = PyLong_FromSsize_t(i);
        return likely(key) ? mm->mp_ass_subscript(o, key, (PyObject *)NULL) : -1;
    }
    if (likely(sm && sm->sq_ass_item)) {
        if (wraparound && unlikely(i < 0) && likely(sm->sq_length)) {
            Py_ssize_t l = sm->sq_length(o);
            if (likely(l >= 0)) {
                i += l;
            } else {
                // if length > max(Py_ssize_t), maybe the object can wrap around itself?
                if (!PyErr_ExceptionMatches(PyExc_OverflowError))
                    return -1;
                PyErr_Clear();
            }
        }
        return sm->sq_ass_item(o, i, (PyObject *)NULL);
    }
#endif
    return __Pyx_DelItem_Generic(o, PyLong_FromSsize_t(i));
}


/////////////// SliceObject.proto ///////////////

// we pass pointer addresses to show the C compiler what is NULL and what isn't
{{if access == 'Get'}}
static CYTHON_INLINE PyObject* __Pyx_PyObject_GetSlice(
        PyObject* obj, Py_ssize_t cstart, Py_ssize_t cstop,
        PyObject** py_start, PyObject** py_stop, PyObject** py_slice,
        int has_cstart, int has_cstop, int wraparound);
{{else}}
#define __Pyx_PyObject_DelSlice(obj, cstart, cstop, py_start, py_stop, py_slice, has_cstart, has_cstop, wraparound) \
    __Pyx_PyObject_SetSlice(obj, (PyObject*)NULL, cstart, cstop, py_start, py_stop, py_slice, has_cstart, has_cstop, wraparound)

// we pass pointer addresses to show the C compiler what is NULL and what isn't
static CYTHON_INLINE int __Pyx_PyObject_SetSlice(
        PyObject* obj, PyObject* value, Py_ssize_t cstart, Py_ssize_t cstop,
        PyObject** py_start, PyObject** py_stop, PyObject** py_slice,
        int has_cstart, int has_cstop, int wraparound);
{{endif}}

/////////////// SliceObject ///////////////

{{if access == 'Get'}}
static CYTHON_INLINE PyObject* __Pyx_PyObject_GetSlice(PyObject* obj,
{{else}}
static CYTHON_INLINE int __Pyx_PyObject_SetSlice(PyObject* obj, PyObject* value,
{{endif}}
        Py_ssize_t cstart, Py_ssize_t cstop,
        PyObject** _py_start, PyObject** _py_stop, PyObject** _py_slice,
        int has_cstart, int has_cstop, CYTHON_UNUSED int wraparound) {
    __Pyx_TypeName obj_type_name;
#if CYTHON_USE_TYPE_SLOTS
    PyMappingMethods* mp = Py_TYPE(obj)->tp_as_mapping;
{{if access == 'Get'}}
    if (likely(mp && mp->mp_subscript))
{{else}}
    if (likely(mp && mp->mp_ass_subscript))
{{endif}}
#endif
    {
        {{if access == 'Get'}}PyObject*{{else}}int{{endif}} result;
        PyObject *py_slice, *py_start, *py_stop;
        if (_py_slice) {
            py_slice = *_py_slice;
        } else {
            PyObject* owned_start = NULL;
            PyObject* owned_stop = NULL;
            if (_py_start) {
                py_start = *_py_start;
            } else {
                if (has_cstart) {
                    owned_start = py_start = PyLong_FromSsize_t(cstart);
                    if (unlikely(!py_start)) goto bad;
                } else
                    py_start = Py_None;
            }
            if (_py_stop) {
                py_stop = *_py_stop;
            } else {
                if (has_cstop) {
                    owned_stop = py_stop = PyLong_FromSsize_t(cstop);
                    if (unlikely(!py_stop)) {
                        Py_XDECREF(owned_start);
                        goto bad;
                    }
                } else
                    py_stop = Py_None;
            }
            py_slice = PySlice_New(py_start, py_stop, Py_None);
            Py_XDECREF(owned_start);
            Py_XDECREF(owned_stop);
            if (unlikely(!py_slice)) goto bad;
        }
#if CYTHON_USE_TYPE_SLOTS
{{if access == 'Get'}}
        result = mp->mp_subscript(obj, py_slice);
#else
        result = PyObject_GetItem(obj, py_slice);
{{else}}
        result = mp->mp_ass_subscript(obj, py_slice, value);
#else
        result = value ? PyObject_SetItem(obj, py_slice, value) : PyObject_DelItem(obj, py_slice);
{{endif}}
#endif
        if (!_py_slice) {
            Py_DECREF(py_slice);
        }
        return result;
    }
    obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
    PyErr_Format(PyExc_TypeError,
{{if access == 'Get'}}
        "'" __Pyx_FMT_TYPENAME "' object is unsliceable", obj_type_name);
{{else}}
        "'" __Pyx_FMT_TYPENAME "' object does not support slice %.10s",
        obj_type_name, value ? "assignment" : "deletion");
{{endif}}
    __Pyx_DECREF_TypeName(obj_type_name);

bad:
    return {{if access == 'Get'}}NULL{{else}}-1{{endif}};
}


/////////////// TupleAndListFromArray.proto ///////////////

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_PyList_FromArray(PyObject *const *src, Py_ssize_t n);
#endif
#if CYTHON_COMPILING_IN_CPYTHON || CYTHON_METH_FASTCALL
static CYTHON_INLINE PyObject* __Pyx_PyTuple_FromArray(PyObject *const *src, Py_ssize_t n);
#endif


/////////////// TupleAndListFromArray ///////////////

#if !CYTHON_COMPILING_IN_CPYTHON && CYTHON_METH_FASTCALL
static CYTHON_INLINE PyObject *
__Pyx_PyTuple_FromArray(PyObject *const *src, Py_ssize_t n)
{
    PyObject *res;
    Py_ssize_t i;
    if (n <= 0) {
        return __Pyx_NewRef(EMPTY(tuple));
    }
    res = PyTuple_New(n);
    if (unlikely(res == NULL)) return NULL;
    for (i = 0; i < n; i++) {
        Py_INCREF(src[i]);
        if (unlikely(__Pyx_PyTuple_SET_ITEM(res, i, src[i]) < (0))) {
            Py_DECREF(res);
            return NULL;
        }
    }
    return res;
}
#elif CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE void __Pyx_copy_object_array(PyObject *const *CYTHON_RESTRICT src, PyObject** CYTHON_RESTRICT dest, Py_ssize_t length) {
    PyObject *v;
    Py_ssize_t i;
    for (i = 0; i < length; i++) {
        v = dest[i] = src[i];
        Py_INCREF(v);
    }
}

static CYTHON_INLINE PyObject *
__Pyx_PyTuple_FromArray(PyObject *const *src, Py_ssize_t n)
{
    PyObject *res;
    if (n <= 0) {
        return __Pyx_NewRef(EMPTY(tuple));
    }
    res = PyTuple_New(n);
    if (unlikely(res == NULL)) return NULL;
    __Pyx_copy_object_array(src, ((PyTupleObject*)res)->ob_item, n);
    return res;
}

static CYTHON_INLINE PyObject *
__Pyx_PyList_FromArray(PyObject *const *src, Py_ssize_t n)
{
    PyObject *res;
    if (n <= 0) {
        return PyList_New(0);
    }
    res = PyList_New(n);
    if (unlikely(res == NULL)) return NULL;
    __Pyx_copy_object_array(src, ((PyListObject*)res)->ob_item, n);
    return res;
}
#endif


/////////////// SliceTupleAndList.proto ///////////////

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_PyList_GetSlice(PyObject* src, Py_ssize_t start, Py_ssize_t stop);
static CYTHON_INLINE PyObject* __Pyx_PyTuple_GetSlice(PyObject* src, Py_ssize_t start, Py_ssize_t stop);
#else
#define __Pyx_PyList_GetSlice(seq, start, stop)   PySequence_GetSlice(seq, start, stop)
#define __Pyx_PyTuple_GetSlice(seq, start, stop)  PySequence_GetSlice(seq, start, stop)
#endif

/////////////// SliceTupleAndList ///////////////
//@requires: TupleAndListFromArray
//@requires: Synchronization.c::CriticalSections

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE void __Pyx_crop_slice(Py_ssize_t* _start, Py_ssize_t* _stop, Py_ssize_t* _length) {
    Py_ssize_t start = *_start, stop = *_stop, length = *_length;
    if (start < 0) {
        start += length;
        if (start < 0)
            start = 0;
    }

    if (stop < 0)
        stop += length;
    else if (stop > length)
        stop = length;

    *_length = stop - start;
    *_start = start;
    *_stop = stop;
}

static CYTHON_INLINE PyObject* __Pyx_PyTuple_GetSlice(
            PyObject* src, Py_ssize_t start, Py_ssize_t stop) {
    Py_ssize_t length = PyTuple_GET_SIZE(src);
    __Pyx_crop_slice(&start, &stop, &length);
    return __Pyx_PyTuple_FromArray(((PyTupleObject*)src)->ob_item + start, length);
}

static CYTHON_INLINE PyObject* __Pyx_PyList_GetSlice_locked(
            PyObject* src, Py_ssize_t start, Py_ssize_t stop) {
    Py_ssize_t length = PyList_GET_SIZE(src);
    __Pyx_crop_slice(&start, &stop, &length);
    if (length <= 0) {
        // Avoid undefined behaviour when accessing `ob_item` of an empty list.
        return PyList_New(0);
    }
    return __Pyx_PyList_FromArray(((PyListObject*)src)->ob_item + start, length);
}

static CYTHON_INLINE PyObject* __Pyx_PyList_GetSlice(
            PyObject* src, Py_ssize_t start, Py_ssize_t stop) {
    PyObject *result;
    __Pyx_BEGIN_CRITICAL_SECTION(src);
    result = __Pyx_PyList_GetSlice_locked(src, start, stop);
    __Pyx_END_CRITICAL_SECTION();
    return result;
}
#endif // CYTHON_COMPILING_IN_CPYTHON


/////////////// CalculateMetaclass.proto ///////////////

static PyObject *__Pyx_CalculateMetaclass(PyTypeObject *metaclass, PyObject *bases);

/////////////// CalculateMetaclass ///////////////

static PyObject *__Pyx_CalculateMetaclass(PyTypeObject *metaclass, PyObject *bases) {
    Py_ssize_t i, nbases;
#if CYTHON_ASSUME_SAFE_SIZE
    nbases = PyTuple_GET_SIZE(bases);
#else
    nbases = PyTuple_Size(bases);
    if (nbases < 0) return NULL;
#endif
    for (i=0; i < nbases; i++) {
        PyTypeObject *tmptype;
#if CYTHON_ASSUME_SAFE_MACROS
        PyObject *tmp = PyTuple_GET_ITEM(bases, i);
#else
        PyObject *tmp = PyTuple_GetItem(bases, i);
        if (!tmp) return NULL;
#endif
        tmptype = Py_TYPE(tmp);
        if (!metaclass) {
            metaclass = tmptype;
            continue;
        }
        if (PyType_IsSubtype(metaclass, tmptype))
            continue;
        if (PyType_IsSubtype(tmptype, metaclass)) {
            metaclass = tmptype;
            continue;
        }
        // else:
        PyErr_SetString(PyExc_TypeError,
                        "metaclass conflict: "
                        "the metaclass of a derived class "
                        "must be a (non-strict) subclass "
                        "of the metaclasses of all its bases");
        return NULL;
    }
    if (!metaclass) {
        metaclass = &PyType_Type;
    }
    // make owned reference
    Py_INCREF((PyObject*) metaclass);
    return (PyObject*) metaclass;
}


/////////////// FindInheritedMetaclass.proto ///////////////

static PyObject *__Pyx_FindInheritedMetaclass(PyObject *bases); /*proto*/

/////////////// FindInheritedMetaclass ///////////////
//@requires: PyObjectGetAttrStr
//@requires: CalculateMetaclass

static PyObject *__Pyx_FindInheritedMetaclass(PyObject *bases) {
    PyObject *metaclass;
    #if CYTHON_ASSUME_SAFE_SIZE
    if (PyTuple_Check(bases) && PyTuple_GET_SIZE(bases) > 0)
    #else
    Py_ssize_t tuple_size = PyTuple_Check(bases) ? PyTuple_Size(bases) : 0;
    if (unlikely(tuple_size < 0)) return NULL;
    if (tuple_size > 0)
    #endif
    {
        PyTypeObject *metatype;
#if CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS
        PyObject *base = PyTuple_GET_ITEM(bases, 0);
#else
        PyObject *base = __Pyx_PySequence_ITEM(bases, 0);
#endif
        metatype = Py_TYPE(base);
        metaclass = __Pyx_CalculateMetaclass(metatype, bases);
#if !(CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS)
        Py_DECREF(base);
#endif
    } else {
        // no bases => use default metaclass
        metaclass = (PyObject *) &PyType_Type;
        Py_INCREF(metaclass);
    }
    return metaclass;
}

/////////////// Py3MetaclassGet.proto ///////////////

static PyObject *__Pyx_Py3MetaclassGet(PyObject *bases, PyObject *mkw); /*proto*/

/////////////// Py3MetaclassGet ///////////////
//@requires: FindInheritedMetaclass
//@requires: CalculateMetaclass

static PyObject *__Pyx_Py3MetaclassGet(PyObject *bases, PyObject *mkw) {
    PyObject *metaclass = mkw ? __Pyx_PyDict_GetItemStr(mkw, PYIDENT("metaclass")) : NULL;
    if (metaclass) {
        Py_INCREF(metaclass);
        if (PyDict_DelItem(mkw, PYIDENT("metaclass")) < 0) {
            Py_DECREF(metaclass);
            return NULL;
        }
        if (PyType_Check(metaclass)) {
            PyObject* orig = metaclass;
            metaclass = __Pyx_CalculateMetaclass((PyTypeObject*) metaclass, bases);
            Py_DECREF(orig);
        }
        return metaclass;
    }
    return __Pyx_FindInheritedMetaclass(bases);
}

/////////////// CreateClass.proto ///////////////

static PyObject *__Pyx_CreateClass(PyObject *bases, PyObject *dict, PyObject *name,
                                   PyObject *qualname, PyObject *modname); /*proto*/

/////////////// CreateClass ///////////////
//@requires: FindInheritedMetaclass
//@requires: CalculateMetaclass

static PyObject *__Pyx_CreateClass(PyObject *bases, PyObject *dict, PyObject *name,
                                   PyObject *qualname, PyObject *modname) {
    PyObject *result;
    PyObject *metaclass;

    if (unlikely(PyDict_SetItem(dict, PYIDENT("__module__"), modname) < 0))
        return NULL;
    if (unlikely(PyDict_SetItem(dict, PYIDENT("__qualname__"), qualname) < 0))
        return NULL;

    /* Python2 __metaclass__ */
    metaclass = __Pyx_PyDict_GetItemStr(dict, PYIDENT("__metaclass__"));
    if (metaclass) {
        Py_INCREF(metaclass);
        if (PyType_Check(metaclass)) {
            PyObject* orig = metaclass;
            metaclass = __Pyx_CalculateMetaclass((PyTypeObject*) metaclass, bases);
            Py_DECREF(orig);
        }
    } else {
        metaclass = __Pyx_FindInheritedMetaclass(bases);
    }
    if (unlikely(!metaclass))
        return NULL;
    result = PyObject_CallFunctionObjArgs(metaclass, name, bases, dict, NULL);
    Py_DECREF(metaclass);
    return result;
}

/////////////// Py3UpdateBases.proto ///////////////

static PyObject* __Pyx_PEP560_update_bases(PyObject *bases); /* proto */

/////////////// Py3UpdateBases /////////////////////
//@requires: PyObjectCallOneArg
//@requires: PyObjectGetAttrStrNoError

/* Shamelessly adapted from cpython/bltinmodule.c update_bases */
static PyObject*
__Pyx_PEP560_update_bases(PyObject *bases)
{
    Py_ssize_t i, j, size_bases;
    PyObject *base = NULL, *meth, *new_base, *result, *new_bases = NULL;
    /*assert(PyTuple_Check(bases));*/

#if CYTHON_ASSUME_SAFE_SIZE
    size_bases = PyTuple_GET_SIZE(bases);
#else
    size_bases = PyTuple_Size(bases);
    if (size_bases < 0) return NULL;
#endif
    for (i = 0; i < size_bases; i++) {
        // original code in CPython: base  = args[i];
#if CYTHON_AVOID_BORROWED_REFS
        Py_CLEAR(base);
#endif
#if CYTHON_ASSUME_SAFE_MACROS
        base = PyTuple_GET_ITEM(bases, i);
#else
        base = PyTuple_GetItem(bases, i);
        if (!base) goto error;
#endif
#if CYTHON_AVOID_BORROWED_REFS
        Py_INCREF(base);
#endif
        if (PyType_Check(base)) {
            if (new_bases) {
                // If we already have made a replacement, then we append every normal base,
                // otherwise just skip it.
                if (PyList_Append(new_bases, base) < 0) {
                    goto error;
                }
            }
            continue;
        }
        // original code in CPython:
        // if (_PyObject_LookupAttrId(base, &PyId___mro_entries__, &meth) < 0) {
        meth = __Pyx_PyObject_GetAttrStrNoError(base, PYIDENT("__mro_entries__"));
        if (!meth && PyErr_Occurred()) {
            goto error;
        }
        if (!meth) {
            if (new_bases) {
                if (PyList_Append(new_bases, base) < 0) {
                    goto error;
                }
            }
            continue;
        }
        new_base = __Pyx_PyObject_CallOneArg(meth, bases);
        Py_DECREF(meth);
        if (!new_base) {
            goto error;
        }
        if (!PyTuple_Check(new_base)) {
            PyErr_SetString(PyExc_TypeError,
                            "__mro_entries__ must return a tuple");
            Py_DECREF(new_base);
            goto error;
        }
        if (!new_bases) {
            // If this is a first successful replacement, create new_bases list and
            // copy previously encountered bases.
            if (!(new_bases = PyList_New(i))) {
                goto error;
            }
            for (j = 0; j < i; j++) {
                // original code in CPython: base = args[j];
                // We use a different name here to keep refcounting separate from base
                PyObject *base_from_list;
#if CYTHON_ASSUME_SAFE_MACROS
                base_from_list = PyTuple_GET_ITEM(bases, j);
                PyList_SET_ITEM(new_bases, j, base_from_list);
                Py_INCREF(base_from_list);
#else
                base_from_list = PyTuple_GetItem(bases, j);
                if (!base_from_list) goto error;
                Py_INCREF(base_from_list);
                if (PyList_SetItem(new_bases, j, base_from_list) < 0) goto error;
#endif
            }
        }
#if CYTHON_ASSUME_SAFE_SIZE
        j = PyList_GET_SIZE(new_bases);
#else
        j = PyList_Size(new_bases);
        if (j < 0) goto error;
#endif
        if (PyList_SetSlice(new_bases, j, j, new_base) < 0) {
            goto error;
        }
        Py_DECREF(new_base);
    }
    if (!new_bases) {
        // unlike the CPython implementation, always return a new reference
        Py_INCREF(bases);
        return bases;
    }
    result = PyList_AsTuple(new_bases);
    Py_DECREF(new_bases);
#if CYTHON_AVOID_BORROWED_REFS
    Py_XDECREF(base);
#endif
    return result;

error:
    Py_XDECREF(new_bases);
#if CYTHON_AVOID_BORROWED_REFS
    Py_XDECREF(base);
#endif
    return NULL;
}

/////////////// Py3ClassCreate.proto ///////////////

static PyObject *__Pyx_Py3MetaclassPrepare(PyObject *metaclass, PyObject *bases, PyObject *name, PyObject *qualname,
                                           PyObject *mkw, PyObject *modname, PyObject *doc); /*proto*/
static PyObject *__Pyx_Py3ClassCreate(PyObject *metaclass, PyObject *name, PyObject *bases, PyObject *dict,
                                      PyObject *mkw, int calculate_metaclass, int allow_py2_metaclass); /*proto*/

/////////////// Py3ClassCreate ///////////////
//@requires: PyObjectGetAttrStrNoError
//@requires: CalculateMetaclass
//@requires: PyObjectFastCall
//@requires: PyObjectCall2Args
//@requires: PyObjectLookupSpecial
// only in fallback code:

static PyObject *__Pyx_Py3MetaclassPrepare(PyObject *metaclass, PyObject *bases, PyObject *name,
                                           PyObject *qualname, PyObject *mkw, PyObject *modname, PyObject *doc) {
    PyObject *ns;
    if (metaclass) {
        PyObject *prep = __Pyx_PyObject_GetAttrStrNoError(metaclass, PYIDENT("__prepare__"));
        if (prep) {
            PyObject *pargs[3] = {NULL, name, bases};
            ns = __Pyx_PyObject_FastCallDict(prep, pargs+1, 2 | __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET, mkw);
            Py_DECREF(prep);
        } else {
            if (unlikely(PyErr_Occurred()))
                return NULL;
            ns = PyDict_New();
        }
    } else {
        ns = PyDict_New();
    }

    if (unlikely(!ns))
        return NULL;

    /* Required here to emulate assignment order */
    if (unlikely(PyObject_SetItem(ns, PYIDENT("__module__"), modname) < 0)) goto bad;
    if (unlikely(PyObject_SetItem(ns, PYIDENT("__qualname__"), qualname) < 0)) goto bad;
    if (unlikely(doc && PyObject_SetItem(ns, PYIDENT("__doc__"), doc) < 0)) goto bad;
    return ns;
bad:
    Py_DECREF(ns);
    return NULL;
}

static PyObject *__Pyx_Py3ClassCreate(PyObject *metaclass, PyObject *name, PyObject *bases,
                                      PyObject *dict, PyObject *mkw,
                                      int calculate_metaclass, int allow_py2_metaclass) {
    PyObject *result;
    PyObject *owned_metaclass = NULL;
    PyObject *margs[4] = {NULL, name, bases, dict};
    if (allow_py2_metaclass) {
        /* honour Python2 __metaclass__ for backward compatibility */
        owned_metaclass = PyObject_GetItem(dict, PYIDENT("__metaclass__"));
        if (owned_metaclass) {
            metaclass = owned_metaclass;
        } else if (likely(PyErr_ExceptionMatches(PyExc_KeyError))) {
            PyErr_Clear();
        } else {
            return NULL;
        }
    }
    if (calculate_metaclass && (!metaclass || PyType_Check(metaclass))) {
        metaclass = __Pyx_CalculateMetaclass((PyTypeObject*) metaclass, bases);
        Py_XDECREF(owned_metaclass);
        if (unlikely(!metaclass))
            return NULL;
        owned_metaclass = metaclass;
    }
    result = __Pyx_PyObject_FastCallDict(metaclass, margs+1, 3 | __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET, mkw);
    Py_XDECREF(owned_metaclass);
    return result;
}

/////////////// ExtTypeTest.proto ///////////////

static CYTHON_INLINE int __Pyx_TypeTest(PyObject *obj, PyTypeObject *type); /*proto*/

/////////////// ExtTypeTest ///////////////

static CYTHON_INLINE int __Pyx_TypeTest(PyObject *obj, PyTypeObject *type) {
    __Pyx_TypeName obj_type_name;
    __Pyx_TypeName type_name;
    if (unlikely(!type)) {
        PyErr_SetString(PyExc_SystemError, "Missing type object");
        return 0;
    }
    if (likely(__Pyx_TypeCheck(obj, type)))
        return 1;
    obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
    type_name = __Pyx_PyType_GetFullyQualifiedName(type);
    PyErr_Format(PyExc_TypeError,
                 "Cannot convert " __Pyx_FMT_TYPENAME " to " __Pyx_FMT_TYPENAME,
                 obj_type_name, type_name);
    __Pyx_DECREF_TypeName(obj_type_name);
    __Pyx_DECREF_TypeName(type_name);
    return 0;
}

/////////////// CallableCheck.proto ///////////////

// PyPy does not set tp_call on classes with metaclass.
// See https://github.com/cython/cython/issues/6892
#if CYTHON_USE_TYPE_SLOTS && !CYTHON_COMPILING_IN_PYPY
#define __Pyx_PyCallable_Check(obj)   (Py_TYPE(obj)->tp_call != NULL)
#else
#define __Pyx_PyCallable_Check(obj)   PyCallable_Check(obj)
#endif

/////////////// PyDictContains.proto ///////////////

static CYTHON_INLINE int __Pyx_PyDict_ContainsTF(PyObject* item, PyObject* dict, int eq) {
    int result = PyDict_Contains(dict, item);
    return unlikely(result < 0) ? result : (result == (eq == Py_EQ));
}

/////////////// PySetContains.proto ///////////////

static CYTHON_INLINE int __Pyx_PySet_ContainsTF(PyObject* key, PyObject* set, int eq); /* proto */

/////////////// PySetContains ///////////////
//@requires: Builtins.c::pyfrozenset_new

static int __Pyx_PySet_ContainsUnhashable(PyObject *set, PyObject *key) {
    int result = -1;
    if (PySet_Check(key) && PyErr_ExceptionMatches(PyExc_TypeError)) {
        /* Convert key to frozenset */
        PyObject *tmpkey;
        PyErr_Clear();
        tmpkey = __Pyx_PyFrozenSet_New(key);
        if (tmpkey != NULL) {
            result = PySet_Contains(set, tmpkey);
            Py_DECREF(tmpkey);
        }
    }
    return result;
}

static CYTHON_INLINE int __Pyx_PySet_ContainsTF(PyObject* key, PyObject* set, int eq) {
    int result = PySet_Contains(set, key);

    if (unlikely(result < 0)) {
        result = __Pyx_PySet_ContainsUnhashable(set, key);
    }
    return unlikely(result < 0) ? result : (result == (eq == Py_EQ));
}

/////////////// PySequenceContains.proto ///////////////

static CYTHON_INLINE int __Pyx_PySequence_ContainsTF(PyObject* item, PyObject* seq, int eq) {
    int result = PySequence_Contains(seq, item);
    return unlikely(result < 0) ? result : (result == (eq == Py_EQ));
}

/////////////// PyBoolOrNullFromLong.proto ///////////////

static CYTHON_INLINE PyObject* __Pyx_PyBoolOrNull_FromLong(long b) {
    return unlikely(b < 0) ? NULL : __Pyx_PyBool_FromLong(b);
}

/////////////// GetBuiltinName.proto ///////////////

static PyObject *__Pyx_GetBuiltinName(PyObject *name); /*proto*/

/////////////// GetBuiltinName ///////////////
//@requires: PyObjectGetAttrStrNoError

static PyObject *__Pyx_GetBuiltinName(PyObject *name) {
    PyObject* result = __Pyx_PyObject_GetAttrStrNoError(NAMED_CGLOBAL(builtins_cname), name);
    if (unlikely(!result) && !PyErr_Occurred()) {
        PyErr_Format(PyExc_NameError,
            "name '%U' is not defined", name);
    }
    return result;
}

/////////////// GetNameInClass.proto ///////////////

#define __Pyx_GetNameInClass(var, nmspace, name)  (var) = __Pyx__GetNameInClass(nmspace, name)
static PyObject *__Pyx__GetNameInClass(PyObject *nmspace, PyObject *name); /*proto*/

/////////////// GetNameInClass ///////////////
//@requires: GetModuleGlobalName

static PyObject *__Pyx__GetNameInClass(PyObject *nmspace, PyObject *name) {
    PyObject *result;
    PyObject *dict;
    assert(PyType_Check(nmspace));
#if CYTHON_USE_TYPE_SLOTS
    dict = ((PyTypeObject*)nmspace)->tp_dict;
    Py_XINCREF(dict);
#else
    dict = PyObject_GetAttr(nmspace, PYIDENT("__dict__"));
#endif
    if (likely(dict)) {
        result = PyObject_GetItem(dict, name);
        Py_DECREF(dict);
        if (result) {
            return result;
        }
    }
    PyErr_Clear();
    __Pyx_GetModuleGlobalNameUncached(result, name);
    return result;
}


/////////////// SetNameInClass.proto ///////////////

#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX < 0x030d0000
// Identifier names are always interned and have a pre-calculated hash value.
#define __Pyx_SetNameInClass(ns, name, value) \
    (likely(PyDict_CheckExact(ns)) ? _PyDict_SetItem_KnownHash(ns, name, value, ((PyASCIIObject *) name)->hash) : PyObject_SetItem(ns, name, value))
#elif CYTHON_COMPILING_IN_CPYTHON
#define __Pyx_SetNameInClass(ns, name, value) \
    (likely(PyDict_CheckExact(ns)) ? PyDict_SetItem(ns, name, value) : PyObject_SetItem(ns, name, value))
#else
#define __Pyx_SetNameInClass(ns, name, value)  PyObject_SetItem(ns, name, value)
#endif

/////////////// SetNewInClass.proto ///////////////

static int __Pyx_SetNewInClass(PyObject *ns, PyObject *name, PyObject *value);

/////////////// SetNewInClass ///////////////
//@requires: SetNameInClass

// Special-case setting __new__: if it's a Cython function, wrap it in a
// staticmethod. This is similar to what Python does for a Python function
// called __new__.
static int __Pyx_SetNewInClass(PyObject *ns, PyObject *name, PyObject *value) {
#ifdef __Pyx_CyFunction_USED
    int ret;
    if (__Pyx_CyFunction_Check(value)) {
        PyObject *staticnew;
#if !CYTHON_COMPILING_IN_LIMITED_API
        staticnew = PyStaticMethod_New(value);
#else
        PyObject *builtins, *staticmethod_str, *staticmethod;
        builtins = PyEval_GetBuiltins(); // borrowed
        if (!builtins) return -1;
        staticmethod_str = PyUnicode_FromStringAndSize("staticmethod", 12);
        if (!staticmethod_str) return -1;
        staticmethod = PyObject_GetItem(builtins, staticmethod_str);
        Py_DECREF(staticmethod_str);
        if (!staticmethod) return -1;
        staticnew = PyObject_CallFunctionObjArgs(staticmethod, value, NULL);
        Py_DECREF(staticmethod);
#endif
        if (unlikely(!staticnew)) return -1;
        ret = __Pyx_SetNameInClass(ns, name, staticnew);
        Py_DECREF(staticnew);
        return ret;
    }
#endif
    return __Pyx_SetNameInClass(ns, name, value);
}


/////////////// GetModuleGlobalName.proto ///////////////
//@requires: PyDictVersioning

#if CYTHON_USE_DICT_VERSIONS
#define __Pyx_GetModuleGlobalName(var, name)  do { \
    static PY_UINT64_T __pyx_dict_version = 0; \
    static PyObject *__pyx_dict_cached_value = NULL; \
    (var) = (likely(__pyx_dict_version == __PYX_GET_DICT_VERSION(NAMED_CGLOBAL(moddict_cname)))) ? \
        (likely(__pyx_dict_cached_value) ? __Pyx_NewRef(__pyx_dict_cached_value) : __Pyx_GetBuiltinName(name)) : \
        __Pyx__GetModuleGlobalName(name, &__pyx_dict_version, &__pyx_dict_cached_value); \
} while(0)
#define __Pyx_GetModuleGlobalNameUncached(var, name)  do { \
    PY_UINT64_T __pyx_dict_version; \
    PyObject *__pyx_dict_cached_value; \
    (var) = __Pyx__GetModuleGlobalName(name, &__pyx_dict_version, &__pyx_dict_cached_value); \
} while(0)
static PyObject *__Pyx__GetModuleGlobalName(PyObject *name, PY_UINT64_T *dict_version, PyObject **dict_cached_value); /*proto*/
#else
#define __Pyx_GetModuleGlobalName(var, name)  (var) = __Pyx__GetModuleGlobalName(name)
#define __Pyx_GetModuleGlobalNameUncached(var, name)  (var) = __Pyx__GetModuleGlobalName(name)
static CYTHON_INLINE PyObject *__Pyx__GetModuleGlobalName(PyObject *name); /*proto*/
#endif


/////////////// GetModuleGlobalName ///////////////
//@requires: GetBuiltinName
//@substitute: naming

#if CYTHON_USE_DICT_VERSIONS
static PyObject *__Pyx__GetModuleGlobalName(PyObject *name, PY_UINT64_T *dict_version, PyObject **dict_cached_value)
#else
static CYTHON_INLINE PyObject *__Pyx__GetModuleGlobalName(PyObject *name)
#endif
{
    PyObject *result;
#if CYTHON_COMPILING_IN_LIMITED_API
    if (unlikely(!$module_cname)) {
        if (!PyErr_Occurred())
            PyErr_SetNone(PyExc_NameError);
        return NULL;
    }
    result = PyObject_GetAttr($module_cname, name);
    if (likely(result)) {
        return result;
    }
    PyErr_Clear();
#elif CYTHON_AVOID_BORROWED_REFS || CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
    if (unlikely(__Pyx_PyDict_GetItemRef(NAMED_CGLOBAL(moddict_cname), name, &result) == -1)) PyErr_Clear();
    __PYX_UPDATE_DICT_CACHE(NAMED_CGLOBAL(moddict_cname), result, *dict_cached_value, *dict_version)
    if (likely(result)) {
        return result;
    }
#else
    // Identifier names are always interned and have a pre-calculated hash value.
    result = _PyDict_GetItem_KnownHash(NAMED_CGLOBAL(moddict_cname), name, ((PyASCIIObject *) name)->hash);
    __PYX_UPDATE_DICT_CACHE(NAMED_CGLOBAL(moddict_cname), result, *dict_cached_value, *dict_version)
    if (likely(result)) {
        return __Pyx_NewRef(result);
    }
    PyErr_Clear();
#endif
    return __Pyx_GetBuiltinName(name);
}

//////////////////// GetAttr.proto ////////////////////

static CYTHON_INLINE PyObject *__Pyx_GetAttr(PyObject *, PyObject *); /*proto*/

//////////////////// GetAttr ////////////////////
//@requires: PyObjectGetAttrStr

static CYTHON_INLINE PyObject *__Pyx_GetAttr(PyObject *o, PyObject *n) {
#if CYTHON_USE_TYPE_SLOTS
    if (likely(PyUnicode_Check(n)))
        return __Pyx_PyObject_GetAttrStr(o, n);
#endif
    return PyObject_GetAttr(o, n);
}


/////////////// PyObjectLookupSpecial.proto ///////////////

#if CYTHON_USE_PYTYPE_LOOKUP && CYTHON_USE_TYPE_SLOTS
#define __Pyx_PyObject_LookupSpecialNoError(obj, attr_name)  __Pyx__PyObject_LookupSpecial(obj, attr_name, 0)
#define __Pyx_PyObject_LookupSpecial(obj, attr_name)  __Pyx__PyObject_LookupSpecial(obj, attr_name, 1)

static CYTHON_INLINE PyObject* __Pyx__PyObject_LookupSpecial(PyObject* obj, PyObject* attr_name, int with_error); /*proto*/

#else
#define __Pyx_PyObject_LookupSpecialNoError(o,n) __Pyx_PyObject_GetAttrStrNoError(o,n)
#define __Pyx_PyObject_LookupSpecial(o,n) __Pyx_PyObject_GetAttrStr(o,n)
#endif

/////////////// PyObjectLookupSpecial ///////////////
//@requires: PyObjectGetAttrStr
//@requires: PyObjectGetAttrStrNoError

#if CYTHON_USE_PYTYPE_LOOKUP && CYTHON_USE_TYPE_SLOTS
static CYTHON_INLINE PyObject* __Pyx__PyObject_LookupSpecial(PyObject* obj, PyObject* attr_name, int with_error) {
    PyObject *res;
    PyTypeObject *tp = Py_TYPE(obj);
    // adapted from CPython's special_lookup() in ceval.c
    res = _PyType_Lookup(tp, attr_name);
    if (likely(res)) {
        descrgetfunc f = Py_TYPE(res)->tp_descr_get;
        if (!f) {
            Py_INCREF(res);
        } else {
            res = f(res, obj, (PyObject *)tp);
        }
    } else if (with_error) {
        PyErr_SetObject(PyExc_AttributeError, attr_name);
    }
    return res;
}
#endif

/////////////// PyObjectGetAttrStrNoError.proto ///////////////

static CYTHON_INLINE PyObject* __Pyx_PyObject_GetAttrStrNoError(PyObject* obj, PyObject* attr_name);/*proto*/

/////////////// PyObjectGetAttrStrNoError ///////////////
//@requires: PyObjectGetAttrStr
//@requires: Exceptions.c::PyThreadStateGet
//@requires: Exceptions.c::PyErrFetchRestore
//@requires: Exceptions.c::PyErrExceptionMatches

#if __PYX_LIMITED_VERSION_HEX < 0x030d0000
static void __Pyx_PyObject_GetAttrStr_ClearAttributeError(void) {
    __Pyx_PyThreadState_declare
    __Pyx_PyThreadState_assign
    if (likely(__Pyx_PyErr_ExceptionMatches(PyExc_AttributeError)))
        __Pyx_PyErr_Clear();
}
#endif

static CYTHON_INLINE PyObject* __Pyx_PyObject_GetAttrStrNoError(PyObject* obj, PyObject* attr_name) {
    PyObject *result;
#if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
    (void) PyObject_GetOptionalAttr(obj, attr_name, &result);
    return result;
#else
#if CYTHON_COMPILING_IN_CPYTHON && CYTHON_USE_TYPE_SLOTS
    // _PyObject_GenericGetAttrWithDict() in CPython 3.7+ can avoid raising the AttributeError.
    // See https://bugs.python.org/issue32544
    PyTypeObject* tp = Py_TYPE(obj);
    if (likely(tp->tp_getattro == PyObject_GenericGetAttr)) {
        return _PyObject_GenericGetAttrWithDict(obj, attr_name, NULL, 1);
    }
#endif
    result = __Pyx_PyObject_GetAttrStr(obj, attr_name);
    if (unlikely(!result)) {
        __Pyx_PyObject_GetAttrStr_ClearAttributeError();
    }
    return result;
#endif
}


/////////////// PyObjectGetAttrStr.proto ///////////////

#if CYTHON_USE_TYPE_SLOTS
static CYTHON_INLINE PyObject* __Pyx_PyObject_GetAttrStr(PyObject* obj, PyObject* attr_name);/*proto*/
#else
#define __Pyx_PyObject_GetAttrStr(o,n) PyObject_GetAttr(o,n)
#endif

/////////////// PyObjectGetAttrStr ///////////////

#if CYTHON_USE_TYPE_SLOTS
static CYTHON_INLINE PyObject* __Pyx_PyObject_GetAttrStr(PyObject* obj, PyObject* attr_name) {
    PyTypeObject* tp = Py_TYPE(obj);
    if (likely(tp->tp_getattro))
        return tp->tp_getattro(obj, attr_name);
    return PyObject_GetAttr(obj, attr_name);
}
#endif

/////////////// PyObjectDelAttr.proto //////////////////

#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
#define __Pyx_PyObject_DelAttr(o, n) PyObject_SetAttr(o, n, NULL)
#else
#define __Pyx_PyObject_DelAttr(o, n) PyObject_DelAttr(o, n)
#endif

/////////////// PyObjectSetAttrStr.proto ///////////////
//@requires: PyObjectDelAttr

#if CYTHON_USE_TYPE_SLOTS
#define __Pyx_PyObject_DelAttrStr(o,n) __Pyx_PyObject_SetAttrStr(o, n, NULL)
static CYTHON_INLINE int __Pyx_PyObject_SetAttrStr(PyObject* obj, PyObject* attr_name, PyObject* value);/*proto*/
#else
#define __Pyx_PyObject_DelAttrStr(o,n)   __Pyx_PyObject_DelAttr(o,n)
#define __Pyx_PyObject_SetAttrStr(o,n,v) PyObject_SetAttr(o,n,v)
#endif

/////////////// PyObjectSetAttrStr ///////////////

#if CYTHON_USE_TYPE_SLOTS
static CYTHON_INLINE int __Pyx_PyObject_SetAttrStr(PyObject* obj, PyObject* attr_name, PyObject* value) {
    PyTypeObject* tp = Py_TYPE(obj);
    if (likely(tp->tp_setattro))
        return tp->tp_setattro(obj, attr_name, value);
    return PyObject_SetAttr(obj, attr_name, value);
}
#endif


/////////////// PyObjectGetMethod.proto ///////////////

#if !(CYTHON_VECTORCALL && (__PYX_LIMITED_VERSION_HEX >= 0x030C0000 || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x03090000)))
static int __Pyx_PyObject_GetMethod(PyObject *obj, PyObject *name, PyObject **method);/*proto*/
#endif

/////////////// PyObjectGetMethod ///////////////
//@requires: PyObjectGetAttrStr

#if !(CYTHON_VECTORCALL && (__PYX_LIMITED_VERSION_HEX >= 0x030C0000 || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x03090000)))
static int __Pyx_PyObject_GetMethod(PyObject *obj, PyObject *name, PyObject **method) {
    PyObject *attr;
#if CYTHON_UNPACK_METHODS && CYTHON_COMPILING_IN_CPYTHON && CYTHON_USE_PYTYPE_LOOKUP
    __Pyx_TypeName type_name;
    // Copied from _PyObject_GetMethod() in CPython 3.7
    PyTypeObject *tp = Py_TYPE(obj);
    PyObject *descr;
    descrgetfunc f = NULL;
    PyObject **dictptr, *dict;
    int meth_found = 0;

    assert (*method == NULL);

    if (unlikely(tp->tp_getattro != PyObject_GenericGetAttr)) {
        attr = __Pyx_PyObject_GetAttrStr(obj, name);
        goto try_unpack;
    }
    if (unlikely(tp->tp_dict == NULL) && unlikely(PyType_Ready(tp) < 0)) {
        return 0;
    }

    descr = _PyType_Lookup(tp, name);
    if (likely(descr != NULL)) {
        Py_INCREF(descr);
#if defined(Py_TPFLAGS_METHOD_DESCRIPTOR) && Py_TPFLAGS_METHOD_DESCRIPTOR
        if (__Pyx_PyType_HasFeature(Py_TYPE(descr), Py_TPFLAGS_METHOD_DESCRIPTOR))
#else
        // Repeating the condition below accommodates for MSVC's inability to test macros inside of macro expansions.
        #ifdef __Pyx_CyFunction_USED
        if (likely(PyFunction_Check(descr) || __Pyx_IS_TYPE(descr, &PyMethodDescr_Type) || __Pyx_CyFunction_Check(descr)))
        #else
        if (likely(PyFunction_Check(descr) || __Pyx_IS_TYPE(descr, &PyMethodDescr_Type)))
        #endif
#endif
        {
            meth_found = 1;
        } else {
            f = Py_TYPE(descr)->tp_descr_get;
            if (f != NULL && PyDescr_IsData(descr)) {
                attr = f(descr, obj, (PyObject *)Py_TYPE(obj));
                Py_DECREF(descr);
                goto try_unpack;
            }
        }
    }

    dictptr = _PyObject_GetDictPtr(obj);
    if (dictptr != NULL && (dict = *dictptr) != NULL) {
        Py_INCREF(dict);
        attr = __Pyx_PyDict_GetItemStr(dict, name);
        if (attr != NULL) {
            Py_INCREF(attr);
            Py_DECREF(dict);
            Py_XDECREF(descr);
            goto try_unpack;
        }
        Py_DECREF(dict);
    }

    if (meth_found) {
        *method = descr;
        return 1;
    }

    if (f != NULL) {
        attr = f(descr, obj, (PyObject *)Py_TYPE(obj));
        Py_DECREF(descr);
        goto try_unpack;
    }

    if (likely(descr != NULL)) {
        *method = descr;
        return 0;
    }

    type_name = __Pyx_PyType_GetFullyQualifiedName(tp);
    PyErr_Format(PyExc_AttributeError,
                 "'" __Pyx_FMT_TYPENAME "' object has no attribute '%U'",
                 type_name, name);
    __Pyx_DECREF_TypeName(type_name);
    return 0;

// Generic fallback implementation using normal attribute lookup.
#else
    attr = __Pyx_PyObject_GetAttrStr(obj, name);
    goto try_unpack;
#endif

try_unpack:
#if CYTHON_UNPACK_METHODS
    // Even if we failed to avoid creating a bound method object, it's still worth unpacking it now, if possible.
    if (likely(attr) && PyMethod_Check(attr) && likely(PyMethod_GET_SELF(attr) == obj)) {
        PyObject *function = PyMethod_GET_FUNCTION(attr);
        Py_INCREF(function);
        Py_DECREF(attr);
        *method = function;
        return 1;
    }
#endif
    *method = attr;
    return 0;
}
#endif


/////////////// UnpackUnboundCMethod.proto ///////////////
//@requires:  Synchronization.c::Atomics

typedef struct {
    PyObject *type;
    PyObject **method_name;
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
    // 0 for uninitialized, 1 for initializing, 2 for initialized
    __pyx_atomic_int_type initialized;
#endif
    // "func" is set on first access (direct C function pointer)
    PyCFunction func;
    // "method" is set on first access (fallback)
    PyObject *method;
    int flag;
} __Pyx_CachedCFunction;

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
static CYTHON_INLINE int __Pyx_CachedCFunction_GetAndSetInitializing(__Pyx_CachedCFunction *cfunc) {
#if !CYTHON_ATOMICS
    // always say "we're initializing". This'll always go an inefficient route,
    // but in reality we always expect to have atomics.
    return 1;
#else
    __pyx_nonatomic_int_type expected = 0;
    if (__pyx_atomic_int_cmp_exchange(&cfunc->initialized, &expected, 1)) {
        // we were uninitialized
        return 0;
    }
    return expected;
#endif
}

static CYTHON_INLINE void __Pyx_CachedCFunction_SetFinishedInitializing(__Pyx_CachedCFunction *cfunc) {
#if CYTHON_ATOMICS
    __pyx_atomic_store(&cfunc->initialized, 2);
#endif
}
#else
#define __Pyx_CachedCFunction_GetAndSetInitializing(cfunc) 2
#define __Pyx_CachedCFunction_SetFinishedInitializing(cfunc)
#endif

/////////////// UnpackUnboundCMethod ///////////////
//@requires: PyObjectGetAttrStr

#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030C0000
static PyObject *__Pyx_SelflessCall(PyObject *method, PyObject *args, PyObject *kwargs) {
    PyObject *result;
    PyObject *selfless_args = PyTuple_GetSlice(args, 1, PyTuple_Size(args));
    if (unlikely(!selfless_args)) return NULL;

    result = PyObject_Call(method, selfless_args, kwargs);
    Py_DECREF(selfless_args);
    return result;
}
#elif CYTHON_COMPILING_IN_PYPY && PY_VERSION_HEX < 0x03090000
// PyPy 3.8 does not define 'args' as 'const'.
// See https://github.com/cython/cython/issues/6867
static PyObject *__Pyx_SelflessCall(PyObject *method, PyObject **args, Py_ssize_t nargs, PyObject *kwnames) {
        return _PyObject_Vectorcall
            (method, args ? args+1 : NULL, nargs ? nargs-1 : 0, kwnames);
}
#else
static PyObject *__Pyx_SelflessCall(PyObject *method, PyObject *const *args, Py_ssize_t nargs, PyObject *kwnames) {
    return
#if PY_VERSION_HEX < 0x03090000
    _PyObject_Vectorcall
#else
    PyObject_Vectorcall
#endif
        (method, args ? args+1 : NULL, nargs ? (size_t) nargs-1 : 0, kwnames);
}
#endif

static PyMethodDef __Pyx_UnboundCMethod_Def = {
    /* .ml_name  = */ "CythonUnboundCMethod",
    /* .ml_meth  = */ __PYX_REINTERPRET_FUNCION(PyCFunction, __Pyx_SelflessCall),
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030C0000
    /* .ml_flags = */ METH_VARARGS | METH_KEYWORDS,
#else
    /* .ml_flags = */ METH_FASTCALL | METH_KEYWORDS,
#endif
    /* .ml_doc   = */ NULL
};

static int __Pyx_TryUnpackUnboundCMethod(__Pyx_CachedCFunction* target) {
    PyObject *method, *result=NULL;
    method = __Pyx_PyObject_GetAttrStr(target->type, *target->method_name);
    if (unlikely(!method))
        return -1;
    result = method;
// FIXME: use functionality from CythonFunction.c/ClassMethod
#if CYTHON_COMPILING_IN_CPYTHON
    if (likely(__Pyx_TypeCheck(method, &PyMethodDescr_Type)))
    {
        PyMethodDescrObject *descr = (PyMethodDescrObject*) method;
        target->func = descr->d_method->ml_meth;
        target->flag = descr->d_method->ml_flags & ~(METH_CLASS | METH_STATIC | METH_COEXIST | METH_STACKLESS);
    } else
#endif
    // bound classmethods need special treatment
#if CYTHON_COMPILING_IN_PYPY
    // In PyPy, functions are regular methods, so just do the self check.
#else
    if (PyCFunction_Check(method))
#endif
    {
        PyObject *self;
        int self_found;
#if CYTHON_COMPILING_IN_LIMITED_API || CYTHON_COMPILING_IN_PYPY
        self = PyObject_GetAttrString(method, "__self__");
        if (!self) {
            PyErr_Clear();
        }
#else
        self = PyCFunction_GET_SELF(method);
#endif
        self_found = (self && self != Py_None);
#if CYTHON_COMPILING_IN_LIMITED_API || CYTHON_COMPILING_IN_PYPY
        Py_XDECREF(self);
#endif
        if (self_found) {
            PyObject *unbound_method = PyCFunction_New(&__Pyx_UnboundCMethod_Def, method);
            if (unlikely(!unbound_method)) return -1;
            // New PyCFunction will own method reference, thus decref __Pyx_PyObject_GetAttrStr
            Py_DECREF(method);
            result = unbound_method;
        }
    }
#if !CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    // Unlikely corner-case where another thread has initialized target->method first.
    // Can't happen in free-threading because the whole thing is locked.
    if (unlikely(target->method)) {
        Py_DECREF(result);
    } else
#endif
    target->method = result;
    return 0;
}

/////////////// CallCFunction.proto ///////////////
#define __Pyx_CallCFunction(cfunc, self, args) \
    ((PyCFunction)(void(*)(void))(cfunc)->func)(self, args)
#define __Pyx_CallCFunctionWithKeywords(cfunc, self, args, kwargs) \
    ((PyCFunctionWithKeywords)(void(*)(void))(cfunc)->func)(self, args, kwargs)
#define __Pyx_CallCFunctionFast(cfunc, self, args, nargs) \
    ((__Pyx_PyCFunctionFast)(void(*)(void))(PyCFunction)(cfunc)->func)(self, args, nargs)
#define __Pyx_CallCFunctionFastWithKeywords(cfunc, self, args, nargs, kwnames) \
    ((__Pyx_PyCFunctionFastWithKeywords)(void(*)(void))(PyCFunction)(cfunc)->func)(self, args, nargs, kwnames)

/////////////// CallUnboundCMethod0.proto ///////////////

CYTHON_UNUSED
static PyObject* __Pyx__CallUnboundCMethod0(__Pyx_CachedCFunction* cfunc, PyObject* self); /*proto*/

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_CallUnboundCMethod0(__Pyx_CachedCFunction* cfunc, PyObject* self); /* proto */
#else
#define __Pyx_CallUnboundCMethod0(cfunc, self)  __Pyx__CallUnboundCMethod0(cfunc, self)
#endif

/////////////// CallUnboundCMethod0 ///////////////
//@requires: CallCFunction
//@requires: UnpackUnboundCMethod
//@requires: PyObjectCallOneArg

#if CYTHON_COMPILING_IN_CPYTHON
// FASTCALL methods receive "&empty_tuple" as simple "PyObject[0]*"
static CYTHON_INLINE PyObject* __Pyx_CallUnboundCMethod0(__Pyx_CachedCFunction* cfunc, PyObject* self) {
    int was_initialized = __Pyx_CachedCFunction_GetAndSetInitializing(cfunc);

    if (likely(was_initialized == 2 && cfunc->func)) {
        if (likely(cfunc->flag == METH_NOARGS))
            return __Pyx_CallCFunction(cfunc, self, NULL);
        if (likely(cfunc->flag == METH_FASTCALL))
            return __Pyx_CallCFunctionFast(cfunc, self, NULL, 0);
        if (cfunc->flag == (METH_FASTCALL | METH_KEYWORDS))
            return __Pyx_CallCFunctionFastWithKeywords(cfunc, self, NULL, 0, NULL);
        if (likely(cfunc->flag == (METH_VARARGS | METH_KEYWORDS)))
            return __Pyx_CallCFunctionWithKeywords(cfunc, self, EMPTY(tuple), NULL);
        if (cfunc->flag == METH_VARARGS)
            return __Pyx_CallCFunction(cfunc, self, EMPTY(tuple));
        return __Pyx__CallUnboundCMethod0(cfunc, self);
    }
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    else if (unlikely(was_initialized == 1)) {
        // If it's being simultaneously initialized, just work on a temp
        __Pyx_CachedCFunction tmp_cfunc = {
#ifndef __cplusplus
            0
#endif
        };
        tmp_cfunc.type = cfunc->type;
        tmp_cfunc.method_name = cfunc->method_name;
        return __Pyx__CallUnboundCMethod0(&tmp_cfunc, self);
    }
#endif
    PyObject *result = __Pyx__CallUnboundCMethod0(cfunc, self);
    __Pyx_CachedCFunction_SetFinishedInitializing(cfunc);
    return result;
}
#endif

static PyObject* __Pyx__CallUnboundCMethod0(__Pyx_CachedCFunction* cfunc, PyObject* self) {
    PyObject *result;
    if (unlikely(!cfunc->method) && unlikely(__Pyx_TryUnpackUnboundCMethod(cfunc) < 0)) return NULL;
    result = __Pyx_PyObject_CallOneArg(cfunc->method, self);
    return result;
}


/////////////// CallUnboundCMethod1.proto ///////////////

CYTHON_UNUSED
static PyObject* __Pyx__CallUnboundCMethod1(__Pyx_CachedCFunction* cfunc, PyObject* self, PyObject* arg);/*proto*/

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_CallUnboundCMethod1(__Pyx_CachedCFunction* cfunc, PyObject* self, PyObject* arg);/*proto*/
#else
#define __Pyx_CallUnboundCMethod1(cfunc, self, arg)  __Pyx__CallUnboundCMethod1(cfunc, self, arg)
#endif

/////////////// CallUnboundCMethod1 ///////////////
//@requires: CallCFunction
//@requires: UnpackUnboundCMethod
//@requires: PyObjectCall2Args

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_CallUnboundCMethod1(__Pyx_CachedCFunction* cfunc, PyObject* self, PyObject* arg) {
    int was_initialized =  __Pyx_CachedCFunction_GetAndSetInitializing(cfunc);

    if (likely(was_initialized == 2 && cfunc->func)) {
        int flag = cfunc->flag;
        if (flag == METH_O) {
            return __Pyx_CallCFunction(cfunc, self, arg);
        } else if (flag == METH_FASTCALL) {
            return __Pyx_CallCFunctionFast(cfunc, self, &arg, 1);
        } else if (flag == (METH_FASTCALL | METH_KEYWORDS)) {
            return __Pyx_CallCFunctionFastWithKeywords(cfunc, self, &arg, 1, NULL);
        }
    }
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    else if (unlikely(was_initialized == 1)) {
        // Race to initialize - use a temp
        __Pyx_CachedCFunction tmp_cfunc = {
#ifndef __cplusplus
            0
#endif
        };
        tmp_cfunc.type = cfunc->type;
        tmp_cfunc.method_name = cfunc->method_name;
        return __Pyx__CallUnboundCMethod1(&tmp_cfunc, self, arg);
    }
#endif
    PyObject* result = __Pyx__CallUnboundCMethod1(cfunc, self, arg);
    __Pyx_CachedCFunction_SetFinishedInitializing(cfunc);
    return result;
}
#endif

static PyObject* __Pyx__CallUnboundCMethod1(__Pyx_CachedCFunction* cfunc, PyObject* self, PyObject* arg){
    PyObject *result = NULL;
    if (unlikely(!cfunc->func && !cfunc->method) && unlikely(__Pyx_TryUnpackUnboundCMethod(cfunc) < 0)) return NULL;
#if CYTHON_COMPILING_IN_CPYTHON
    if (cfunc->func && (cfunc->flag & METH_VARARGS)) {
        PyObject *args = PyTuple_New(1);
        if (unlikely(!args)) return NULL;
        Py_INCREF(arg);
        PyTuple_SET_ITEM(args, 0, arg);
        if (cfunc->flag & METH_KEYWORDS)
            result = __Pyx_CallCFunctionWithKeywords(cfunc, self, args, NULL);
        else
            result = __Pyx_CallCFunction(cfunc, self, args);
        Py_DECREF(args);
    } else
#endif
    {
        result = __Pyx_PyObject_Call2Args(cfunc->method, self, arg);
    }
    return result;
}


/////////////// CallUnboundCMethod2.proto ///////////////

CYTHON_UNUSED
static PyObject* __Pyx__CallUnboundCMethod2(__Pyx_CachedCFunction* cfunc, PyObject* self, PyObject* arg1, PyObject* arg2); /*proto*/

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject *__Pyx_CallUnboundCMethod2(__Pyx_CachedCFunction *cfunc, PyObject *self, PyObject *arg1, PyObject *arg2); /*proto*/
#else
#define __Pyx_CallUnboundCMethod2(cfunc, self, arg1, arg2)  __Pyx__CallUnboundCMethod2(cfunc, self, arg1, arg2)
#endif

/////////////// CallUnboundCMethod2 ///////////////
//@requires: CallCFunction
//@requires: UnpackUnboundCMethod
//@requires: PyObjectFastCall

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject *__Pyx_CallUnboundCMethod2(__Pyx_CachedCFunction *cfunc, PyObject *self, PyObject *arg1, PyObject *arg2) {
    int was_initialized = __Pyx_CachedCFunction_GetAndSetInitializing(cfunc);

    if (likely(was_initialized == 2 && cfunc->func)) {
        PyObject *args[2] = {arg1, arg2};
        if (cfunc->flag == METH_FASTCALL) {
            return __Pyx_CallCFunctionFast(cfunc, self, args, 2);
        }
        if (cfunc->flag == (METH_FASTCALL | METH_KEYWORDS))
            return __Pyx_CallCFunctionFastWithKeywords(cfunc, self, args, 2, NULL);
    }
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    else if (unlikely(was_initialized == 1)) {
        // Race to initialize - run this on a temp function instead.
        __Pyx_CachedCFunction tmp_cfunc = {
#ifndef __cplusplus
            0
#endif
        };
        tmp_cfunc.type = cfunc->type;
        tmp_cfunc.method_name = cfunc->method_name;
        return __Pyx__CallUnboundCMethod2(&tmp_cfunc, self, arg1, arg2);
    }
#endif
    PyObject *result = __Pyx__CallUnboundCMethod2(cfunc, self, arg1, arg2);
    __Pyx_CachedCFunction_SetFinishedInitializing(cfunc);
    return result;
}
#endif

static PyObject* __Pyx__CallUnboundCMethod2(__Pyx_CachedCFunction* cfunc, PyObject* self, PyObject* arg1, PyObject* arg2){
    if (unlikely(!cfunc->func && !cfunc->method) && unlikely(__Pyx_TryUnpackUnboundCMethod(cfunc) < 0)) return NULL;
#if CYTHON_COMPILING_IN_CPYTHON
    if (cfunc->func && (cfunc->flag & METH_VARARGS)) {
        PyObject *result = NULL;
        PyObject *args = PyTuple_New(2);
        if (unlikely(!args)) return NULL;
        Py_INCREF(arg1);
        PyTuple_SET_ITEM(args, 0, arg1);
        Py_INCREF(arg2);
        PyTuple_SET_ITEM(args, 1, arg2);
        if (cfunc->flag & METH_KEYWORDS)
            result = __Pyx_CallCFunctionWithKeywords(cfunc, self, args, NULL);
        else
            result = __Pyx_CallCFunction(cfunc, self, args);
        Py_DECREF(args);
        return result;
    }
#endif
    {
        PyObject *args[4] = {NULL, self, arg1, arg2};
        return __Pyx_PyObject_FastCall(cfunc->method, args+1, 3 | __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET);
    }
}


/////////////// PyObjectFastCall.proto ///////////////

#define __Pyx_PyObject_FastCall(func, args, nargs)  __Pyx_PyObject_FastCallDict(func, args, (size_t)(nargs), NULL)
static CYTHON_INLINE PyObject* __Pyx_PyObject_FastCallDict(PyObject *func, PyObject * const*args, size_t nargs, PyObject *kwargs); /*proto*/

/////////////// PyObjectFastCall ///////////////
//@requires: PyObjectCall
//@requires: PyObjectCallMethO

#if PY_VERSION_HEX < 0x03090000 || CYTHON_COMPILING_IN_LIMITED_API
static PyObject* __Pyx_PyObject_FastCall_fallback(PyObject *func, PyObject * const*args, size_t nargs, PyObject *kwargs) {
    PyObject *argstuple;
    PyObject *result = 0;
    size_t i;

    argstuple = PyTuple_New((Py_ssize_t)nargs);
    if (unlikely(!argstuple)) return NULL;
    for (i = 0; i < nargs; i++) {
        Py_INCREF(args[i]);
        if (__Pyx_PyTuple_SET_ITEM(argstuple, (Py_ssize_t)i, args[i]) != (0)) goto bad;
    }
    result = __Pyx_PyObject_Call(func, argstuple, kwargs);
  bad:
    Py_DECREF(argstuple);
    return result;
}
#endif

#if CYTHON_VECTORCALL && !CYTHON_COMPILING_IN_LIMITED_API
  #if PY_VERSION_HEX < 0x03090000
    #define __Pyx_PyVectorcall_Function(callable) _PyVectorcall_Function(callable)
  #elif CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE vectorcallfunc __Pyx_PyVectorcall_Function(PyObject *callable) {
    PyTypeObject *tp = Py_TYPE(callable);

    #if defined(__Pyx_CyFunction_USED)
    if (__Pyx_CyFunction_CheckExact(callable)) {
        return __Pyx_CyFunction_func_vectorcall(callable);
    }
    #endif

    if (!PyType_HasFeature(tp, Py_TPFLAGS_HAVE_VECTORCALL)) {
        return NULL;
    }
    assert(PyCallable_Check(callable));

    Py_ssize_t offset = tp->tp_vectorcall_offset;
    assert(offset > 0);

    vectorcallfunc ptr;
    memcpy(&ptr, (char *) callable + offset, sizeof(ptr));
    return ptr;
}
  #else
    #define __Pyx_PyVectorcall_Function(callable) PyVectorcall_Function(callable)
  #endif
#endif

static CYTHON_INLINE PyObject* __Pyx_PyObject_FastCallDict(PyObject *func, PyObject *const *args, size_t _nargs, PyObject *kwargs) {
    // Special fast paths for 0 and 1 arguments
    // NOTE: in many cases, this is called with a constant value for nargs
    // which is known at compile-time. So the branches below will typically
    // be optimized away.
    Py_ssize_t nargs = __Pyx_PyVectorcall_NARGS(_nargs);
#if CYTHON_COMPILING_IN_CPYTHON
    if (nargs == 0 && kwargs == NULL) {
        if (__Pyx_CyOrPyCFunction_Check(func) && likely( __Pyx_CyOrPyCFunction_GET_FLAGS(func) & METH_NOARGS))
            return __Pyx_PyObject_CallMethO(func, NULL);
    }
    else if (nargs == 1 && kwargs == NULL) {
        if (__Pyx_CyOrPyCFunction_Check(func) && likely( __Pyx_CyOrPyCFunction_GET_FLAGS(func) & METH_O))
            return __Pyx_PyObject_CallMethO(func, args[0]);
    }
#endif

    if (kwargs == NULL) {
        #if CYTHON_VECTORCALL
          #if CYTHON_COMPILING_IN_LIMITED_API
            return PyObject_Vectorcall(func, args, _nargs, NULL);
          #else
            vectorcallfunc f = __Pyx_PyVectorcall_Function(func);
            if (f) {
                return f(func, args, _nargs, NULL);
            }
          #endif
        #endif
    }

    if (nargs == 0) {
        return __Pyx_PyObject_Call(func, EMPTY(tuple), kwargs);
    }
    #if PY_VERSION_HEX >= 0x03090000 && !CYTHON_COMPILING_IN_LIMITED_API
    return PyObject_VectorcallDict(func, args, (size_t)nargs, kwargs);
    #else
    return __Pyx_PyObject_FastCall_fallback(func, args, (size_t)nargs, kwargs);
    #endif
}

/////////////// PyObjectFastCallMethod.proto ///////////////
//@requires: PyObjectFastCall

#if CYTHON_VECTORCALL && PY_VERSION_HEX >= 0x03090000
#define __Pyx_PyObject_FastCallMethod(name, args, nargsf) PyObject_VectorcallMethod(name, args, nargsf, NULL)
#else
static PyObject *__Pyx_PyObject_FastCallMethod(PyObject *name, PyObject *const *args, size_t nargsf); /* proto */
#endif

/////////////// PyObjectFastCallMethod ///////////////

#if !CYTHON_VECTORCALL || PY_VERSION_HEX < 0x03090000
static PyObject *__Pyx_PyObject_FastCallMethod(PyObject *name, PyObject *const *args, size_t nargsf) {
    PyObject *result;
    PyObject *attr = PyObject_GetAttr(args[0], name);
    if (unlikely(!attr))
        return NULL;
    result = __Pyx_PyObject_FastCall(attr, args+1, nargsf - 1);
    Py_DECREF(attr);
    return result;
}
#endif

/////////////// PyObjectVectorCallKwBuilder.proto ////////////////
//@requires: PyObjectFastCall
// For versions that define PyObject_Vectorcall, use PyObject_Vectorcall and define functions to build a kwnames tuple and add arguments to args.
// For versions that don't, use __Pyx_PyObject_FastCallDict and functions to build a keyword dictionary

CYTHON_UNUSED static int __Pyx_VectorcallBuilder_AddArg_Check(PyObject *key, PyObject *value, PyObject *builder, PyObject **args, int n); /* proto */

#if CYTHON_VECTORCALL
#if PY_VERSION_HEX >= 0x03090000
#define __Pyx_Object_Vectorcall_CallFromBuilder PyObject_Vectorcall
#else
#define __Pyx_Object_Vectorcall_CallFromBuilder _PyObject_Vectorcall
#endif

#define __Pyx_MakeVectorcallBuilderKwds(n) PyTuple_New(n)

static int __Pyx_VectorcallBuilder_AddArg(PyObject *key, PyObject *value, PyObject *builder, PyObject **args, int n); /* proto */
static int __Pyx_VectorcallBuilder_AddArgStr(const char *key, PyObject *value, PyObject *builder, PyObject **args, int n); /* proto */
#else
#define __Pyx_Object_Vectorcall_CallFromBuilder __Pyx_PyObject_FastCallDict

#define __Pyx_MakeVectorcallBuilderKwds(n) __Pyx_PyDict_NewPresized(n)

#define __Pyx_VectorcallBuilder_AddArg(key, value, builder, args, n) PyDict_SetItem(builder, key, value)
#define __Pyx_VectorcallBuilder_AddArgStr(key, value, builder, args, n) PyDict_SetItemString(builder, key, value)
#endif

/////////////// PyObjectVectorCallKwBuilder ////////////////

#if CYTHON_VECTORCALL
static int __Pyx_VectorcallBuilder_AddArg(PyObject *key, PyObject *value, PyObject *builder, PyObject **args, int n) {
    (void)__Pyx_PyObject_FastCallDict;

    Py_INCREF(key);
    if (__Pyx_PyTuple_SET_ITEM(builder, n, key) != (0)) return -1;
    args[n] = value;
    return 0;
}

CYTHON_UNUSED static int __Pyx_VectorcallBuilder_AddArg_Check(PyObject *key, PyObject *value, PyObject *builder, PyObject **args, int n) {
    (void)__Pyx_VectorcallBuilder_AddArgStr;
    if (unlikely(!PyUnicode_Check(key))) {
        PyErr_SetString(PyExc_TypeError, "keywords must be strings");
        return -1;
    }
    return __Pyx_VectorcallBuilder_AddArg(key, value, builder, args, n);
}

static int __Pyx_VectorcallBuilder_AddArgStr(const char *key, PyObject *value, PyObject *builder, PyObject **args, int n) {
    PyObject *pyKey = PyUnicode_FromString(key);
    if (!pyKey) return -1;
    return __Pyx_VectorcallBuilder_AddArg(pyKey, value, builder, args, n);
}
#else // CYTHON_VECTORCALL
CYTHON_UNUSED static int __Pyx_VectorcallBuilder_AddArg_Check(PyObject *key, PyObject *value, PyObject *builder, CYTHON_UNUSED PyObject **args, CYTHON_UNUSED int n) {
    if (unlikely(!PyUnicode_Check(key))) {
        PyErr_SetString(PyExc_TypeError, "keywords must be strings");
        return -1;
    }
    return PyDict_SetItem(builder, key, value);
}
#endif

/////////////// PyObjectVectorCallMethodKwBuilder.proto ////////////////

#if CYTHON_VECTORCALL && PY_VERSION_HEX >= 0x03090000
#define __Pyx_Object_VectorcallMethod_CallFromBuilder PyObject_VectorcallMethod
#else
static PyObject *__Pyx_Object_VectorcallMethod_CallFromBuilder(PyObject *name, PyObject *const *args, size_t nargsf, PyObject *kwnames); /* proto */
#endif

/////////////// PyObjectVectorCallMethodKwBuilder ////////////////
//@requires: PyObjectVectorCallKwBuilder

#if !CYTHON_VECTORCALL || PY_VERSION_HEX < 0x03090000
static PyObject *__Pyx_Object_VectorcallMethod_CallFromBuilder(PyObject *name, PyObject *const *args, size_t nargsf, PyObject *kwnames) {
    PyObject *result;
    PyObject *obj = PyObject_GetAttr(args[0], name);
    if (unlikely(!obj))
        return NULL;
    result = __Pyx_Object_Vectorcall_CallFromBuilder(obj, args+1, nargsf-1, kwnames);
    Py_DECREF(obj);
    return result;
}
#endif

/////////////// PyObjectCallMethod0.proto ///////////////

static PyObject* __Pyx_PyObject_CallMethod0(PyObject* obj, PyObject* method_name); /*proto*/

/////////////// PyObjectCallMethod0 ///////////////
//@requires: PyObjectGetMethod
//@requires: PyObjectCallOneArg
//@requires: PyObjectCallNoArg

static PyObject* __Pyx_PyObject_CallMethod0(PyObject* obj, PyObject* method_name) {
#if CYTHON_VECTORCALL && (__PYX_LIMITED_VERSION_HEX >= 0x030C0000 || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x03090000))
    PyObject *args[1] = {obj};
    // avoid unused functions
    (void) __Pyx_PyObject_CallOneArg;
    (void) __Pyx_PyObject_CallNoArg;
    return PyObject_VectorcallMethod(method_name, args, 1 | PY_VECTORCALL_ARGUMENTS_OFFSET, NULL);
#else
    PyObject *method = NULL, *result = NULL;
    int is_method = __Pyx_PyObject_GetMethod(obj, method_name, &method);
    if (likely(is_method)) {
        result = __Pyx_PyObject_CallOneArg(method, obj);
        Py_DECREF(method);
        return result;
    }
    if (unlikely(!method)) goto bad;
    result = __Pyx_PyObject_CallNoArg(method);
    Py_DECREF(method);
bad:
    return result;
#endif
}


/////////////// PyObjectCallMethod1.proto ///////////////

static PyObject* __Pyx_PyObject_CallMethod1(PyObject* obj, PyObject* method_name, PyObject* arg); /*proto*/

/////////////// PyObjectCallMethod1 ///////////////
//@requires: PyObjectGetMethod
//@requires: PyObjectCallOneArg
//@requires: PyObjectCall2Args

#if !(CYTHON_VECTORCALL && (__PYX_LIMITED_VERSION_HEX >= 0x030C0000 || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x03090000)))
static PyObject* __Pyx__PyObject_CallMethod1(PyObject* method, PyObject* arg) {
    // Separate function to avoid excessive inlining.
    PyObject *result = __Pyx_PyObject_CallOneArg(method, arg);
    Py_DECREF(method);
    return result;
}
#endif

static PyObject* __Pyx_PyObject_CallMethod1(PyObject* obj, PyObject* method_name, PyObject* arg) {
#if CYTHON_VECTORCALL && (__PYX_LIMITED_VERSION_HEX >= 0x030C0000 || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x03090000))
    PyObject *args[2] = {obj, arg};
    // avoid unused functions
    (void) __Pyx_PyObject_CallOneArg;
    (void) __Pyx_PyObject_Call2Args;
    return PyObject_VectorcallMethod(method_name, args, 2 | PY_VECTORCALL_ARGUMENTS_OFFSET, NULL);
#else
    PyObject *method = NULL, *result;
    int is_method = __Pyx_PyObject_GetMethod(obj, method_name, &method);
    if (likely(is_method)) {
        result = __Pyx_PyObject_Call2Args(method, obj, arg);
        Py_DECREF(method);
        return result;
    }
    if (unlikely(!method)) return NULL;
    return __Pyx__PyObject_CallMethod1(method, arg);
#endif
}


/////////////// tp_new.proto ///////////////

#define __Pyx_tp_new(type_obj, args) __Pyx_tp_new_kwargs(type_obj, args, NULL)
static CYTHON_INLINE PyObject* __Pyx_tp_new_kwargs(PyObject* type_obj, PyObject* args, PyObject* kwargs); /*proto*/

/////////////// tp_new ///////////////

#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
static PyObject* __Pyx_tp_new_fallback(PyObject* type_obj, PyObject* args, PyObject* kwargs) {
    PyObject *new_func = NULL, *new_args = NULL, *obj;
    Py_ssize_t i, nargs = PyTuple_Size(args);

    if (unlikely(nargs < 0)) goto bad;
    new_args = PyTuple_New(nargs + 1);
    if (unlikely(!new_args)) goto bad;
    for (i = 0; i < nargs + 1; i++) {
        PyObject *item = (i == 0) ? type_obj : PyTuple_GetItem(args, i - 1);
        if (unlikely(!item)) goto bad;
        Py_INCREF(item);
        if (unlikely(PyTuple_SetItem(new_args, i, item)) < 0) goto bad;
    }

    new_func = PyObject_GetAttrString(type_obj, "__new__");
    if (unlikely(!new_func)) goto bad;
    obj = PyObject_Call(new_func, new_args, kwargs);
    Py_DECREF(new_func);
    Py_DECREF(new_args);
    return obj;

 bad:
    Py_XDECREF(new_func);
    Py_XDECREF(new_args);
    return NULL;
}
#endif

static PyObject* __Pyx_tp_new_kwargs(PyObject* type_obj, PyObject* args, PyObject* kwargs) {
    newfunc tp_new = __Pyx_PyType_TryGetSlot((PyTypeObject*)type_obj, tp_new, newfunc);
#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
    if (!tp_new) return __Pyx_tp_new_fallback(type_obj, args, kwargs);
#else
    assert(tp_new != NULL);
#endif
    return tp_new((PyTypeObject*)type_obj, args, kwargs);
}


/////////////// PyObjectCall.proto ///////////////

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_PyObject_Call(PyObject *func, PyObject *arg, PyObject *kw); /*proto*/
#else
#define __Pyx_PyObject_Call(func, arg, kw) PyObject_Call(func, arg, kw)
#endif

/////////////// PyObjectCall ///////////////

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_PyObject_Call(PyObject *func, PyObject *arg, PyObject *kw) {
    PyObject *result;
    ternaryfunc call = Py_TYPE(func)->tp_call;

    if (unlikely(!call))
        return PyObject_Call(func, arg, kw);
    if (unlikely(Py_EnterRecursiveCall(" while calling a Python object")))
        return NULL;
    result = (*call)(func, arg, kw);
    Py_LeaveRecursiveCall();
    if (unlikely(!result) && unlikely(!PyErr_Occurred())) {
        PyErr_SetString(
            PyExc_SystemError,
            "NULL result without error in PyObject_Call");
    }
    return result;
}
#endif


/////////////// PyObjectCallMethO.proto ///////////////

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_PyObject_CallMethO(PyObject *func, PyObject *arg); /*proto*/
#endif

/////////////// PyObjectCallMethO ///////////////

#if CYTHON_COMPILING_IN_CPYTHON
static CYTHON_INLINE PyObject* __Pyx_PyObject_CallMethO(PyObject *func, PyObject *arg) {
    PyObject *self, *result;
    PyCFunction cfunc;
    cfunc = __Pyx_CyOrPyCFunction_GET_FUNCTION(func);
    self = __Pyx_CyOrPyCFunction_GET_SELF(func);

    if (unlikely(Py_EnterRecursiveCall(" while calling a Python object")))
        return NULL;
    result = cfunc(self, arg);
    Py_LeaveRecursiveCall();
    if (unlikely(!result) && unlikely(!PyErr_Occurred())) {
        PyErr_SetString(
            PyExc_SystemError,
            "NULL result without error in PyObject_Call");
    }
    return result;
}
#endif


/////////////// PyObjectCall2Args.proto ///////////////

static CYTHON_INLINE PyObject* __Pyx_PyObject_Call2Args(PyObject* function, PyObject* arg1, PyObject* arg2); /*proto*/

/////////////// PyObjectCall2Args ///////////////
//@requires: PyObjectFastCall

static CYTHON_INLINE PyObject* __Pyx_PyObject_Call2Args(PyObject* function, PyObject* arg1, PyObject* arg2) {
    PyObject *args[3] = {NULL, arg1, arg2};
    return __Pyx_PyObject_FastCall(function, args+1, 2 | __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET);
}


/////////////// PyObjectCallOneArg.proto ///////////////

static CYTHON_INLINE PyObject* __Pyx_PyObject_CallOneArg(PyObject *func, PyObject *arg); /*proto*/

/////////////// PyObjectCallOneArg ///////////////
//@requires: PyObjectFastCall

static CYTHON_INLINE PyObject* __Pyx_PyObject_CallOneArg(PyObject *func, PyObject *arg) {
    PyObject *args[2] = {NULL, arg};
    return __Pyx_PyObject_FastCall(func, args+1, 1 | __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET);
}


/////////////// PyObjectCallNoArg.proto ///////////////

static CYTHON_INLINE PyObject* __Pyx_PyObject_CallNoArg(PyObject *func); /*proto*/

/////////////// PyObjectCallNoArg ///////////////
//@requires: PyObjectFastCall

static CYTHON_INLINE PyObject* __Pyx_PyObject_CallNoArg(PyObject *func) {
    PyObject *arg[2] = {NULL, NULL};
    return __Pyx_PyObject_FastCall(func, arg + 1, 0 | __Pyx_PY_VECTORCALL_ARGUMENTS_OFFSET);
}


/////////////// PyVectorcallFastCallDict.proto ///////////////

#if CYTHON_METH_FASTCALL && CYTHON_VECTORCALL
static CYTHON_INLINE PyObject *__Pyx_PyVectorcall_FastCallDict(PyObject *func, __pyx_vectorcallfunc vc, PyObject *const *args, size_t nargs, PyObject *kw);
#endif

/////////////// PyVectorcallFastCallDict ///////////////
//@requires: OwnedDictNext

#if CYTHON_METH_FASTCALL && CYTHON_VECTORCALL
// Slow path when kw is non-empty
static PyObject *__Pyx_PyVectorcall_FastCallDict_kw(PyObject *func, __pyx_vectorcallfunc vc, PyObject *const *args, size_t nargs, PyObject *kw)
{
    // Code based on _PyObject_FastCallDict() and _PyStack_UnpackDict() from CPython
    PyObject *res = NULL;
    PyObject *kwnames;
    PyObject **newargs;
    PyObject **kwvalues;
    Py_ssize_t i;
    #if CYTHON_AVOID_BORROWED_REFS
    PyObject *pos;
    #else
    Py_ssize_t pos;
    #endif
    size_t j;
    PyObject *key, *value;
    unsigned long keys_are_strings;
    #if !CYTHON_ASSUME_SAFE_SIZE
    Py_ssize_t nkw = PyDict_Size(kw);
    if (unlikely(nkw == -1)) return NULL;
    #else
    Py_ssize_t nkw = PyDict_GET_SIZE(kw);
    #endif

    // Copy positional arguments
    newargs = (PyObject **)PyMem_Malloc((nargs + (size_t)nkw) * sizeof(args[0]));
    if (unlikely(newargs == NULL)) {
        PyErr_NoMemory();
        return NULL;
    }
    for (j = 0; j < nargs; j++) newargs[j] = args[j];

    // Copy keyword arguments
    kwnames = PyTuple_New(nkw);
    if (unlikely(kwnames == NULL)) {
        PyMem_Free(newargs);
        return NULL;
    }
    kwvalues = newargs + nargs;
    pos = 0;
    i = 0;
    keys_are_strings = Py_TPFLAGS_UNICODE_SUBCLASS;
    while (__Pyx_PyDict_NextRef(kw, &pos, &key, &value)) {
        keys_are_strings &=
        #if CYTHON_COMPILING_IN_LIMITED_API
            PyType_GetFlags(Py_TYPE(key));
        #else
            Py_TYPE(key)->tp_flags;
        #endif
        #if !CYTHON_ASSUME_SAFE_MACROS
        if (unlikely(PyTuple_SetItem(kwnames, i, key) < 0)) goto cleanup;
        #else
        PyTuple_SET_ITEM(kwnames, i, key);
        #endif
        kwvalues[i] = value;
        i++;
    }
    if (unlikely(!keys_are_strings)) {
        PyErr_SetString(PyExc_TypeError, "keywords must be strings");
        goto cleanup;
    }

    // The actual call
    res = vc(func, newargs, nargs, kwnames);

cleanup:
    #if CYTHON_AVOID_BORROWED_REFS
    Py_DECREF(pos);
    #endif
    Py_DECREF(kwnames);
    for (i = 0; i < nkw; i++)
        Py_DECREF(kwvalues[i]);
    PyMem_Free(newargs);
    return res;
}

static CYTHON_INLINE PyObject *__Pyx_PyVectorcall_FastCallDict(PyObject *func, __pyx_vectorcallfunc vc, PyObject *const *args, size_t nargs, PyObject *kw)
{
    Py_ssize_t kw_size =
        likely(kw == NULL) ?
        0 :
#if !CYTHON_ASSUME_SAFE_SIZE
        PyDict_Size(kw);
#else
        PyDict_GET_SIZE(kw);
#endif

    if (kw_size == 0) {
        return vc(func, args, nargs, NULL);
    }
#if !CYTHON_ASSUME_SAFE_SIZE
    else if (unlikely(kw_size == -1)) {
        return NULL;
    }
#endif
    return __Pyx_PyVectorcall_FastCallDict_kw(func, vc, args, nargs, kw);
}
#endif


/////////////// MatrixMultiply.proto ///////////////

// TODO: remove
#define __Pyx_PyNumber_MatrixMultiply(x,y)         PyNumber_MatrixMultiply(x,y)
#define __Pyx_PyNumber_InPlaceMatrixMultiply(x,y)  PyNumber_InPlaceMatrixMultiply(x,y)


/////////////// PyDictVersioning.proto ///////////////

#if CYTHON_USE_DICT_VERSIONS && CYTHON_USE_TYPE_SLOTS
#define __PYX_DICT_VERSION_INIT  ((PY_UINT64_T) -1)
#define __PYX_GET_DICT_VERSION(dict)  (((PyDictObject*)(dict))->ma_version_tag)
#define __PYX_UPDATE_DICT_CACHE(dict, value, cache_var, version_var) \
    (version_var) = __PYX_GET_DICT_VERSION(dict); \
    (cache_var) = (value);

// The lookup function should return a new reference.
#define __PYX_PY_DICT_LOOKUP_IF_MODIFIED(VAR, DICT, LOOKUP) { \
    static PY_UINT64_T __pyx_dict_version = 0; \
    static PyObject *__pyx_dict_cached_value = NULL; \
    if (likely(__PYX_GET_DICT_VERSION(DICT) == __pyx_dict_version)) { \
        (VAR) = __Pyx_XNewRef(__pyx_dict_cached_value); \
    } else { \
        (VAR) = __pyx_dict_cached_value = (LOOKUP); \
        __pyx_dict_version = __PYX_GET_DICT_VERSION(DICT); \
    } \
}

static CYTHON_INLINE PY_UINT64_T __Pyx_get_tp_dict_version(PyObject *obj); /*proto*/
static CYTHON_INLINE PY_UINT64_T __Pyx_get_object_dict_version(PyObject *obj); /*proto*/
static CYTHON_INLINE int __Pyx_object_dict_version_matches(PyObject* obj, PY_UINT64_T tp_dict_version, PY_UINT64_T obj_dict_version); /*proto*/

#else
#define __PYX_GET_DICT_VERSION(dict)  (0)
#define __PYX_UPDATE_DICT_CACHE(dict, value, cache_var, version_var)
#define __PYX_PY_DICT_LOOKUP_IF_MODIFIED(VAR, DICT, LOOKUP)  (VAR) = (LOOKUP);
#endif

/////////////// PyDictVersioning ///////////////

#if CYTHON_USE_DICT_VERSIONS && CYTHON_USE_TYPE_SLOTS
static CYTHON_INLINE PY_UINT64_T __Pyx_get_tp_dict_version(PyObject *obj) {
    PyObject *dict = Py_TYPE(obj)->tp_dict;
    return likely(dict) ? __PYX_GET_DICT_VERSION(dict) : 0;
}

static CYTHON_INLINE PY_UINT64_T __Pyx_get_object_dict_version(PyObject *obj) {
    PyObject **dictptr = NULL;
    Py_ssize_t offset = Py_TYPE(obj)->tp_dictoffset;
    if (offset) {
#if CYTHON_COMPILING_IN_CPYTHON
        dictptr = (likely(offset > 0)) ? (PyObject **) ((char *)obj + offset) : _PyObject_GetDictPtr(obj);
#else
        dictptr = _PyObject_GetDictPtr(obj);
#endif
    }
    return (dictptr && *dictptr) ? __PYX_GET_DICT_VERSION(*dictptr) : 0;
}

static CYTHON_INLINE int __Pyx_object_dict_version_matches(PyObject* obj, PY_UINT64_T tp_dict_version, PY_UINT64_T obj_dict_version) {
    PyObject *dict = Py_TYPE(obj)->tp_dict;
    if (unlikely(!dict) || unlikely(tp_dict_version != __PYX_GET_DICT_VERSION(dict)))
        return 0;
    return obj_dict_version == __Pyx_get_object_dict_version(obj);
}
#endif

////////////// CachedMethodType.module_state_decls ////////////

#if CYTHON_COMPILING_IN_LIMITED_API
PyObject *__Pyx_CachedMethodType;
#endif

////////////// CachedMethodType.init //////////////

#if CYTHON_COMPILING_IN_LIMITED_API
{
    PyObject *typesModule=NULL;
    typesModule = PyImport_ImportModule("types");
    if (typesModule) {
        CGLOBAL(__Pyx_CachedMethodType) = PyObject_GetAttrString(typesModule, "MethodType");
        Py_DECREF(typesModule);
    }
} // error handling follows
#endif

/////////////// CachedMethodType.cleanup ////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
Py_CLEAR(CGLOBAL(__Pyx_CachedMethodType));
#endif

/////////////// PyMethodNew.proto ///////////////

static PyObject *__Pyx_PyMethod_New(PyObject *func, PyObject *self, PyObject *typ); /* proto */

/////////////// PyMethodNew ///////////////
//@requires: CachedMethodType

#if CYTHON_COMPILING_IN_LIMITED_API
static PyObject *__Pyx_PyMethod_New(PyObject *func, PyObject *self, PyObject *typ) {
    PyObject *result;
    CYTHON_UNUSED_VAR(typ);
    if (!self)
        return __Pyx_NewRef(func);
    #if __PYX_LIMITED_VERSION_HEX >= 0x030C0000
    {
        PyObject *args[] = {func, self};
        result = PyObject_Vectorcall(CGLOBAL(__Pyx_CachedMethodType), args, 2, NULL);
    }
    #else
    result = PyObject_CallFunctionObjArgs(CGLOBAL(__Pyx_CachedMethodType), func, self, NULL);
    #endif
    return result;
}
#else
// This should be an actual function (not a macro), such that we can put it
// directly in a tp_descr_get slot.
static PyObject *__Pyx_PyMethod_New(PyObject *func, PyObject *self, PyObject *typ) {
    CYTHON_UNUSED_VAR(typ);
    if (!self)
        return __Pyx_NewRef(func);
    return PyMethod_New(func, self);
}
#endif

///////////// PyMethodNew2Arg.proto /////////////
//@requires: CachedMethodType

// TODO: remove
// Another wrapping of PyMethod_New that matches the Python3 signature
#if CYTHON_COMPILING_IN_LIMITED_API
static PyObject *__Pyx_PyMethod_New2Arg(PyObject *func, PyObject *self); /* proto */
#else
#define __Pyx_PyMethod_New2Arg PyMethod_New
#endif

///////////// PyMethodNew2Arg /////////////
//@requires: CachedMethodType

#if CYTHON_COMPILING_IN_LIMITED_API
static PyObject *__Pyx_PyMethod_New2Arg(PyObject *func, PyObject *self) {
    PyObject *result;
    #if __PYX_LIMITED_VERSION_HEX >= 0x030C0000
    {
        PyObject *args[] = {func, self};
        result = PyObject_Vectorcall(CGLOBAL(__Pyx_CachedMethodType), args, 2, NULL);
    }
    #else
    result = PyObject_CallFunctionObjArgs(CGLOBAL(__Pyx_CachedMethodType), func, self, NULL);
    #endif
    return result;
}
#endif

/////////////// UnicodeConcatInPlace.proto ////////////////

# if CYTHON_COMPILING_IN_CPYTHON
// __Pyx_PyUnicode_ConcatInPlace may modify the first argument 'left'
// However, unlike `PyUnicode_Append` it will never NULL it.
// It behaves like a regular function - returns a new reference and NULL on error
    #if CYTHON_REFNANNY
        #define __Pyx_PyUnicode_ConcatInPlace(left, right, unsafe_shared) __Pyx_PyUnicode_ConcatInPlaceImpl(&left, right, unsafe_shared, __pyx_refnanny)
    #else
        #define __Pyx_PyUnicode_ConcatInPlace(left, right, unsafe_shared) __Pyx_PyUnicode_ConcatInPlaceImpl(&left, right, unsafe_shared)
    #endif
    #define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_DefinitelyUniqueInPlace(left, right) __Pyx_PyUnicode_ConcatInPlace(left, right, __Pyx_ReferenceSharing_DefinitelyUnique)
    #define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_OwnStrongReferenceInPlace(left, right) __Pyx_PyUnicode_ConcatInPlace(left, right, __Pyx_ReferenceSharing_OwnStrongReference)
    #define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_FunctionArgumentInPlace(left, right) __Pyx_PyUnicode_ConcatInPlace(left, right, __Pyx_ReferenceSharing_DefinitelyUnique)
    #define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_SharedReferenceInPlace(left, right) __Pyx_PyUnicode_ConcatInPlace(left, right, __Pyx_ReferenceSharing_SharedReference)
    // __Pyx_PyUnicode_ConcatInPlace is slightly odd because it has the potential to modify the input
    // argument (but only in cases where no user should notice). Therefore, it needs to keep Cython's
    // refnanny informed.
    static CYTHON_INLINE PyObject *__Pyx_PyUnicode_ConcatInPlaceImpl(PyObject **p_left, PyObject *right, int unsafe_shared
        #if CYTHON_REFNANNY
        , void* __pyx_refnanny
        #endif
    ); /* proto */
#else
#define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_DefinitelyUniqueInPlace __Pyx_PyUnicode_Concat
#define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_OwnStrongReferenceInPlace __Pyx_PyUnicode_Concat
#define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_FunctionArgumentInPlace __Pyx_PyUnicode_Concat
#define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_SharedReferenceInPlace __Pyx_PyUnicode_Concat
#endif
#define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_DefinitelyUniqueInPlaceSafe(left, right) \
    ((unlikely((left) == Py_None) || unlikely((right) == Py_None)) ? \
    PyNumber_InPlaceAdd(left, right) : __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_DefinitelyUniqueInPlace(left, right))
#define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_OwnStrongReferenceInPlaceSafe(left, right) \
    ((unlikely((left) == Py_None) || unlikely((right) == Py_None)) ? \
    PyNumber_InPlaceAdd(left, right) : __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_OwnStrongReferenceInPlace(left, right))
#define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_FunctionArgumentInPlaceSafe(left, right) \
    ((unlikely((left) == Py_None) || unlikely((right) == Py_None)) ? \
    PyNumber_InPlaceAdd(left, right) : __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_FunctionArgumentInPlace(left, right))
#define __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_SharedReferenceInPlaceSafe(left, right) \
    ((unlikely((left) == Py_None) || unlikely((right) == Py_None)) ? \
    PyNumber_InPlaceAdd(left, right) : __Pyx_PyUnicode_Concat__Pyx_ReferenceSharing_SharedReferenceInPlace(left, right))

/////////////// UnicodeConcatInPlace ////////////////

# if CYTHON_COMPILING_IN_CPYTHON
// copied directly from unicode_object.c "unicode_modifiable"
// removing _PyUnicode_HASH if we don't have access to it
//  - this is OK because trying PyUnicode_Resize on a non-modifyable
//  object will still work, it just won't happen in place
static int
__Pyx_unicode_modifiable(PyObject *unicode, int unsafe_shared)
{
    if (!__Pyx_IS_UNIQUELY_REFERENCED(unicode, unsafe_shared))
        return 0;
#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX > 0x030F0000
    if (PyUnstable_Unicode_GET_CACHED_HASH(unicode) != -1)
        return 0;
#endif
    if (!PyUnicode_CheckExact(unicode))
        return 0;
    if (PyUnicode_CHECK_INTERNED(unicode))
        return 0;
    return 1;
}

static CYTHON_INLINE PyObject *__Pyx_PyUnicode_ConcatInPlaceImpl(PyObject **p_left, PyObject *right, int unsafe_shared
        #if CYTHON_REFNANNY
        , void* __pyx_refnanny
        #endif
    ) {
    // heavily based on PyUnicode_Append
    PyObject *left = *p_left;
    Py_ssize_t left_len, right_len, new_len;

    if (unlikely(__Pyx_PyUnicode_READY(left) == -1))
        return NULL;
    if (unlikely(__Pyx_PyUnicode_READY(right) == -1))
        return NULL;

    // Shortcuts
    left_len = PyUnicode_GET_LENGTH(left);
    if (left_len == 0) {
        Py_INCREF(right);
        return right;
    }
    right_len = PyUnicode_GET_LENGTH(right);
    if (right_len == 0) {
        Py_INCREF(left);
        return left;
    }
    if (unlikely(left_len > PY_SSIZE_T_MAX - right_len)) {
        PyErr_SetString(PyExc_OverflowError,
                        "strings are too large to concat");
        return NULL;
    }
    new_len = left_len + right_len;

    if (__Pyx_unicode_modifiable(left, unsafe_shared)
            && PyUnicode_CheckExact(right)
            && PyUnicode_KIND(right) <= PyUnicode_KIND(left)
            // Don't resize for ascii += latin1. Convert ascii to latin1 requires
            //   to change the structure size, but characters are stored just after
            //   the structure, and so it requires to move all characters which is
            //   not so different than duplicating the string.
            && !(PyUnicode_IS_ASCII(left) && !PyUnicode_IS_ASCII(right))) {

        int ret;
        // GIVEREF/GOTREF since we expect *p_left to change (although it won't change on failures).
        __Pyx_GIVEREF(*p_left);
        ret = PyUnicode_Resize(p_left, new_len);
        __Pyx_GOTREF(*p_left);
        if (unlikely(ret != 0))
            return NULL;

        // copy 'right' into the newly allocated area of 'left'
        #if PY_VERSION_HEX >= 0x030d0000
        if (unlikely(PyUnicode_CopyCharacters(*p_left, left_len, right, 0, right_len) < 0)) return NULL;
        #else
        _PyUnicode_FastCopyCharacters(*p_left, left_len, right, 0, right_len);
        #endif
        __Pyx_INCREF(*p_left);
        __Pyx_GIVEREF(*p_left);
        return *p_left;
    } else {
        return __Pyx_PyUnicode_Concat(left, right);
    }
  }
#endif


/////////////// PySequenceMultiply.proto ///////////////

#define __Pyx_PySequence_Multiply_Left(mul, seq)  __Pyx_PySequence_Multiply(seq, mul)
#if !CYTHON_USE_TYPE_SLOTS
#define  __Pyx_PySequence_Multiply PySequence_Repeat
#else
static CYTHON_INLINE PyObject* __Pyx_PySequence_Multiply(PyObject *seq, Py_ssize_t mul);
#endif

/////////////// PySequenceMultiply ///////////////

#if CYTHON_USE_TYPE_SLOTS
static PyObject* __Pyx_PySequence_Multiply_Generic(PyObject *seq, Py_ssize_t mul) {
    PyObject *result, *pymul = PyLong_FromSsize_t(mul);
    if (unlikely(!pymul))
        return NULL;
    result = PyNumber_Multiply(seq, pymul);
    Py_DECREF(pymul);
    return result;
}

static CYTHON_INLINE PyObject* __Pyx_PySequence_Multiply(PyObject *seq, Py_ssize_t mul) {
    PyTypeObject *type = Py_TYPE(seq);
    if (likely(type->tp_as_sequence && type->tp_as_sequence->sq_repeat)) {
        return type->tp_as_sequence->sq_repeat(seq, mul);
    } else {
        return __Pyx_PySequence_Multiply_Generic(seq, mul);
    }
}
#endif

/////////////// BuiltinSequenceMultiply.proto ///////////////

static CYTHON_INLINE PyObject* __Pyx_{{typeobj}}_Multiply(PyObject *seq, Py_ssize_t mul);

/////////////// BuiltinSequenceMultiply ////////////////

static CYTHON_INLINE PyObject* __Pyx_{{typeobj}}_Multiply(PyObject *seq, Py_ssize_t mul) {
    // It's important that this function always calls the exact typeobj slot, because it may
    // be used with a subclass deliberately to access the base class slot.
    ssizeargfunc slot = __Pyx_PyType_TryGetSubSlot(&{{typeobj}}, tp_as_sequence, sq_repeat, ssizeargfunc);
    #if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
    if (unlikely(!slot)) {
        return PyObject_CallMethod((PyObject*)&{{typeobj}}, "__mul__", "On", seq, mul);
    }
    #endif
    return slot(seq, mul);
}

/////////////// FormatTypeName.proto ///////////////

#if CYTHON_COMPILING_IN_LIMITED_API
typedef PyObject *__Pyx_TypeName;
#define __Pyx_FMT_TYPENAME "%U"
#define __Pyx_DECREF_TypeName(obj) Py_XDECREF(obj)

#if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
#define __Pyx_PyType_GetFullyQualifiedName PyType_GetFullyQualifiedName
#else
static __Pyx_TypeName __Pyx_PyType_GetFullyQualifiedName(PyTypeObject* tp); /*proto*/
#endif
#else  // !LIMITED_API
typedef const char *__Pyx_TypeName;
#define __Pyx_FMT_TYPENAME "%.200s"
#define __Pyx_PyType_GetFullyQualifiedName(tp) ((tp)->tp_name)
#define __Pyx_DECREF_TypeName(obj)
#endif

/////////////// FormatTypeName ///////////////

#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
static __Pyx_TypeName
__Pyx_PyType_GetFullyQualifiedName(PyTypeObject* tp)
{
    PyObject *module = NULL, *name = NULL, *result = NULL;

    #if __PYX_LIMITED_VERSION_HEX < 0x030b0000
    name = __Pyx_PyObject_GetAttrStr((PyObject *)tp,
                                               PYIDENT("__qualname__"));
    #else
    name = PyType_GetQualName(tp);
    #endif
    if (unlikely(name == NULL) || unlikely(!PyUnicode_Check(name))) goto bad;
    module = __Pyx_PyObject_GetAttrStr((PyObject *)tp,
                                               PYIDENT("__module__"));
    if (unlikely(module == NULL) || unlikely(!PyUnicode_Check(module))) goto bad;

    if (PyUnicode_CompareWithASCIIString(module, "builtins") == 0) {
        result = name;
        name = NULL;
        goto done;
    }

    result = PyUnicode_FromFormat("%U.%U", module, name);
    if (unlikely(result == NULL)) goto bad;

  done:
    Py_XDECREF(name);
    Py_XDECREF(module);
    return result;

  bad:
    PyErr_Clear();
    if (name) {
        // Use this as second-best fallback
        result = name;
        name = NULL;
    } else {
        result = __Pyx_NewRef(PYUNICODE("?"));
    }
    goto done;
}
#endif


/////////////// RaiseUnexpectedTypeError.proto ///////////////

static int __Pyx_RaiseUnexpectedTypeError(const char *expected, PyObject *obj); /*proto*/

/////////////// RaiseUnexpectedTypeError ///////////////

static int
__Pyx_RaiseUnexpectedTypeError(const char *expected, PyObject *obj)
{
    __Pyx_TypeName obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
    PyErr_Format(PyExc_TypeError, "Expected %s, got " __Pyx_FMT_TYPENAME,
                 expected, obj_type_name);
    __Pyx_DECREF_TypeName(obj_type_name);
    return 0;
}


/////////////// RaiseUnboundLocalError.proto ///////////////
static void __Pyx_RaiseUnboundLocalError(const char *varname);/*proto*/

/////////////// RaiseUnboundLocalError ///////////////
static void __Pyx_RaiseUnboundLocalError(const char *varname) {
    PyErr_Format(PyExc_UnboundLocalError, "local variable '%s' referenced before assignment", varname);
}

/////////////// RaiseUnboundLocalErrorNogil.proto ///////////////
static void __Pyx_RaiseUnboundLocalErrorNogil(const char *varname);/*proto*/

/////////////// RaiseUnboundLocalErrorNogil ///////////////
//@requires: RaiseUnboundLocalError

// Don't inline the function, it should really never be called in production
static void __Pyx_RaiseUnboundLocalErrorNogil(const char *varname) {
    PyGILState_STATE gilstate = PyGILState_Ensure();
    __Pyx_RaiseUnboundLocalError(varname);
    PyGILState_Release(gilstate);
}


/////////////// RaiseClosureNameError.proto ///////////////
static void __Pyx_RaiseClosureNameError(const char *varname);/*proto*/

/////////////// RaiseClosureNameError ///////////////
static void __Pyx_RaiseClosureNameError(const char *varname) {
    PyErr_Format(PyExc_NameError, "free variable '%s' referenced before assignment in enclosing scope", varname);
}

/////////////// RaiseClosureNameErrorNogil.proto ///////////////
static void __Pyx_RaiseClosureNameErrorNogil(const char *varname);/*proto*/

/////////////// RaiseClosureNameErrorNogil ///////////////
//@requires: RaiseClosureNameError

static void __Pyx_RaiseClosureNameErrorNogil(const char *varname) {
    PyGILState_STATE gilstate = PyGILState_Ensure();
    __Pyx_RaiseClosureNameError(varname);
    PyGILState_Release(gilstate);
}

//////////////// RaiseCppGlobalNameError.proto ///////////////////////
static void __Pyx_RaiseCppGlobalNameError(const char *varname); /*proto*/

/////////////// RaiseCppGlobalNameError //////////////////////////////
static void __Pyx_RaiseCppGlobalNameError(const char *varname) {
    PyErr_Format(PyExc_NameError, "C++ global '%s' is not initialized", varname);
}

//////////////// RaiseCppGlobalNameErrorNogil.proto ///////////////////////
static void __Pyx_RaiseCppGlobalNameErrorNogil(const char *varname); /*proto*/

/////////////// RaiseCppGlobalNameErrorNogil //////////////////////////////
//@requires: RaiseCppGlobalNameError

static void __Pyx_RaiseCppGlobalNameErrorNogil(const char *varname) {
    PyGILState_STATE gilstate = PyGILState_Ensure();
    __Pyx_RaiseCppGlobalNameError(varname);
    PyGILState_Release(gilstate);
}

//////////////// RaiseCppAttributeError.proto ///////////////////////
static void __Pyx_RaiseCppAttributeError(const char *varname); /*proto*/

/////////////// RaiseCppAttributeError //////////////////////////////
static void __Pyx_RaiseCppAttributeError(const char *varname) {
    PyErr_Format(PyExc_AttributeError, "C++ attribute '%s' is not initialized", varname);
}

//////////////// RaiseCppAttributeErrorNogil.proto ///////////////////////
static void __Pyx_RaiseCppAttributeErrorNogil(const char *varname); /*proto*/

/////////////// RaiseCppAttributeErrorNogil //////////////////////////////
//@requires: RaiseCppGlobalNameError

static void __Pyx_RaiseCppAttributeErrorNogil(const char *varname) {
    PyGILState_STATE gilstate = PyGILState_Ensure();
    __Pyx_RaiseCppAttributeError(varname);
    PyGILState_Release(gilstate);
}

/////////////// ListPack.proto //////////////////////////////////////

// Equivalent to PyTuple_Pack for sections where we want to reduce the code size
static PyObject *__Pyx_PyList_Pack(Py_ssize_t n, ...); /* proto */

/////////////// ListPack //////////////////////////////////////

static PyObject *__Pyx_PyList_Pack(Py_ssize_t n, ...) {
    va_list va;
    PyObject *l = PyList_New(n);
    va_start(va, n);
    if (unlikely(!l)) goto end;

    for (Py_ssize_t i=0; i<n; ++i) {
        PyObject *arg = va_arg(va, PyObject*);
        Py_INCREF(arg);
        if (__Pyx_PyList_SET_ITEM(l, i, arg) != (0)) {
            Py_CLEAR(l);
            goto end;
        }
    }

    end:
    va_end(va);
    return l;
}

/////////////// LengthHint.proto //////////////////////////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
// We could reimplement it fully, however since it's only an optimization
// the simpler version is to make it the default.
#define __Pyx_PyObject_LengthHint(o, defaultval)  (defaultval)
#else
#define __Pyx_PyObject_LengthHint(o, defaultval)  PyObject_LengthHint(o, defaultval)
#endif

//////////////////////// OwnedDictNext.proto ///////////////////////////////

// Like PyDict_Next but pkey and pvalue are new references (and thus must be decrefed).
// pkey and pvalue may be NULL but this must be the same for all subsequent calls.
// If it encounters an exception it prints it as unraisable and stops.
#if CYTHON_AVOID_BORROWED_REFS
static int __Pyx_PyDict_NextRef(PyObject *p, PyObject **ppos, PyObject **pkey, PyObject **pvalue); /* proto */
#else
CYTHON_INLINE
static int __Pyx_PyDict_NextRef(PyObject *p, Py_ssize_t *ppos, PyObject **pkey, PyObject **pvalue); /* proto */
#endif


//////////////////////// OwnedDictNext //////////////////////////////////////
//@requires: Builtins.c::py_dict_values
//@requires: Builtins.c::py_dict_items

#if CYTHON_AVOID_BORROWED_REFS
static int __Pyx_PyDict_NextRef(PyObject *p, PyObject **ppos, PyObject **pkey, PyObject **pvalue) {
    PyObject *next = NULL;
    if (!*ppos) {
        // Intentionally using lazy iterators instead of PyDict_Items/Keys/Values
        // at the cost of a function call.
        if (pvalue) {
            PyObject *dictview = pkey ? __Pyx_PyDict_Items(p) : __Pyx_PyDict_Values(p);
            if (unlikely(!dictview)) goto bad;
            *ppos = PyObject_GetIter(dictview);
            Py_DECREF(dictview);
        } else {
            *ppos = PyObject_GetIter(p);
        }
        if (unlikely(!*ppos)) goto bad;
    }
    next = PyIter_Next(*ppos);
    if (!next) {
        if (PyErr_Occurred()) goto bad;
        return 0;
    }
    if (pkey && pvalue) {
        *pkey = __Pyx_PySequence_ITEM(next, 0);
        if (unlikely(*pkey)) goto bad;
        *pvalue = __Pyx_PySequence_ITEM(next, 1);
        if (unlikely(*pvalue)) goto bad;
        Py_DECREF(next);
    } else if (pkey) {
        *pkey = next;
    } else {
        assert(pvalue);
        *pvalue = next;
    }
    return 1;

  bad:
    Py_XDECREF(next);
    // PyDict_Next can't fail, so neither can this
#if !CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX >= 0x030d0000
    PyErr_FormatUnraisable("Exception ignored in __Pyx_PyDict_NextRef");
#else
    PyErr_WriteUnraisable(PYIDENT("__Pyx_PyDict_NextRef"));
#endif
    if (pkey) *pkey = NULL;
    if (pvalue) *pvalue = NULL;
    return 0;
}

#else // !CYTHON_AVOID_BORROWED_REFS
static int __Pyx_PyDict_NextRef(PyObject *p, Py_ssize_t *ppos, PyObject **pkey, PyObject **pvalue) {
    int result = PyDict_Next(p, ppos, pkey, pvalue);
    if (likely(result == 1)) {
        if (pkey) Py_INCREF(*pkey);
        if (pvalue) Py_INCREF(*pvalue);
    }
    return result;
}
#endif
