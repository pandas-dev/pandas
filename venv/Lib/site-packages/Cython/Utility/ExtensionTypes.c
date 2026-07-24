/////////////// FixUpExtensionType.proto ///////////////

static CYTHON_INLINE int __Pyx_fix_up_extension_type_from_spec(PyType_Spec *spec, PyTypeObject *type); /*proto*/

/////////////// FixUpExtensionType ///////////////
//@requires:ModuleSetupCode.c::IncludeStructmemberH
//@requires:StringTools.c::IncludeStringH
//@requires:SetItemOnTypeDict

static int __Pyx_fix_up_extension_type_from_spec(PyType_Spec *spec, PyTypeObject *type) {
#if __PYX_LIMITED_VERSION_HEX > 0x030900B1
    CYTHON_UNUSED_VAR(spec);
    CYTHON_UNUSED_VAR(type);
    CYTHON_UNUSED_VAR(__Pyx__SetItemOnTypeDict);
#else
    // Set tp_weakreflist, tp_dictoffset, tp_vectorcalloffset
    // Copied and adapted from https://bugs.python.org/issue38140
    const PyType_Slot *slot = spec->slots;
    int changed = 0;
#if !CYTHON_COMPILING_IN_LIMITED_API
    while (slot && slot->slot && slot->slot != Py_tp_members)
        slot++;
    if (slot && slot->slot == Py_tp_members) {
#if !CYTHON_COMPILING_IN_CPYTHON
        const
#endif  // !CYTHON_COMPILING_IN_CPYTHON)
            PyMemberDef *memb = (PyMemberDef*) slot->pfunc;
        while (memb && memb->name) {
            if (memb->name[0] == '_' && memb->name[1] == '_') {
                if (strcmp(memb->name, "__weaklistoffset__") == 0) {
                    // The PyMemberDef must be a Py_ssize_t and readonly.
                    assert(memb->type == T_PYSSIZET);
                    assert(memb->flags == READONLY);
                    type->tp_weaklistoffset = memb->offset;
                    // FIXME: is it even worth calling PyType_Modified() here?
                    changed = 1;
                }
                else if (strcmp(memb->name, "__dictoffset__") == 0) {
                    // The PyMemberDef must be a Py_ssize_t and readonly.
                    assert(memb->type == T_PYSSIZET);
                    assert(memb->flags == READONLY);
                    type->tp_dictoffset = memb->offset;
                    // FIXME: is it even worth calling PyType_Modified() here?
                    changed = 1;
                }
#if CYTHON_METH_FASTCALL
                else if (strcmp(memb->name, "__vectorcalloffset__") == 0) {
                    // The PyMemberDef must be a Py_ssize_t and readonly.
                    assert(memb->type == T_PYSSIZET);
                    assert(memb->flags == READONLY);
                    type->tp_vectorcall_offset = memb->offset;
                    // FIXME: is it even worth calling PyType_Modified() here?
                    changed = 1;
                }
#endif  // CYTHON_METH_FASTCALL
#if !CYTHON_COMPILING_IN_PYPY
                else if (strcmp(memb->name, "__module__") == 0) {
                    // PyType_FromSpec() in CPython <= 3.9b1 overwrites this field with a constant string.
                    // See https://bugs.python.org/issue40703
                    PyObject *descr;
                    // The PyMemberDef must be an object and normally readable, possibly writable.
                    assert(memb->type == T_OBJECT);
                    assert(memb->flags == 0 || memb->flags == READONLY);
                    descr = PyDescr_NewMember(type, memb);
                    if (unlikely(!descr))
                        return -1;
                    int set_item_result = PyDict_SetItem(type->tp_dict, PyDescr_NAME(descr), descr);
                    Py_DECREF(descr);
                    if (unlikely(set_item_result < 0)) {
                        return -1;
                    }
                    changed = 1;
                }
#endif  // !CYTHON_COMPILING_IN_PYPY
            }
            memb++;
        }
    }
#endif  // !CYTHON_COMPILING_IN_LIMITED_API
#if !CYTHON_COMPILING_IN_PYPY
    slot = spec->slots;
    while (slot && slot->slot && slot->slot != Py_tp_getset)
        slot++;
    if (slot && slot->slot == Py_tp_getset) {
        PyGetSetDef *getset = (PyGetSetDef*) slot->pfunc;
        while (getset && getset->name) {
            if (getset->name[0] == '_' && getset->name[1] == '_' && strcmp(getset->name, "__module__") == 0) {
                PyObject *descr = PyDescr_NewGetSet(type, getset);
                if (unlikely(!descr))
                    return -1;
                #if CYTHON_COMPILING_IN_LIMITED_API
                PyObject *pyname = PyUnicode_FromString(getset->name);
                if (unlikely(!pyname)) {
                    Py_DECREF(descr);
                    return -1;
                }
                int set_item_result = __Pyx_SetItemOnTypeDict(type, pyname, descr);
                Py_DECREF(pyname);
                #else
                CYTHON_UNUSED_VAR(__Pyx__SetItemOnTypeDict);
                int set_item_result = PyDict_SetItem(type->tp_dict, PyDescr_NAME(descr), descr);
                #endif
                Py_DECREF(descr);
                if (unlikely(set_item_result < 0)) {
                    return -1;
                }
                changed = 1;
            }
            ++getset;
        }
    }
#else
    CYTHON_UNUSED_VAR(__Pyx__SetItemOnTypeDict);
#endif  // !CYTHON_COMPILING_IN_PYPY
    if (changed)
        PyType_Modified(type);
#endif  // PY_VERSION_HEX > 0x030900B1
    return 0;
}


/////////////// ValidateBasesTuple.proto ///////////////

#if CYTHON_COMPILING_IN_CPYTHON || CYTHON_COMPILING_IN_LIMITED_API || CYTHON_USE_TYPE_SPECS
static int __Pyx_validate_bases_tuple(const char *type_name, Py_ssize_t dictoffset, PyObject *bases); /*proto*/
#endif

/////////////// ValidateBasesTuple ///////////////

#if CYTHON_COMPILING_IN_CPYTHON || CYTHON_COMPILING_IN_LIMITED_API || CYTHON_USE_TYPE_SPECS
static int __Pyx_validate_bases_tuple(const char *type_name, Py_ssize_t dictoffset, PyObject *bases) {
    // Loop over all bases (except the first) and check that those
    // really are heap types. Otherwise, it would not be safe to
    // subclass them.
    //
    // We also check tp_dictoffset: it is unsafe to inherit
    // tp_dictoffset from a base class because the object structures
    // would not be compatible. So, if our extension type doesn't set
    // tp_dictoffset (i.e. there is no __dict__ attribute in the object
    // structure), we need to check that none of the base classes sets
    // it either.
    Py_ssize_t i, n;
#if CYTHON_ASSUME_SAFE_SIZE
    n = PyTuple_GET_SIZE(bases);
#else
    n = PyTuple_Size(bases);
    if (unlikely(n < 0)) return -1;
#endif
    for (i = 1; i < n; i++)  /* Skip first base */
    {
        PyTypeObject *b;
#if CYTHON_AVOID_BORROWED_REFS
        PyObject *b0 = PySequence_GetItem(bases, i);
        if (!b0) return -1;
#elif CYTHON_ASSUME_SAFE_MACROS
        PyObject *b0 = PyTuple_GET_ITEM(bases, i);
#else
        PyObject *b0 = PyTuple_GetItem(bases, i);
        if (!b0) return -1;
#endif
        b = (PyTypeObject*) b0;
        if (!__Pyx_PyType_HasFeature(b, Py_TPFLAGS_HEAPTYPE))
        {
            __Pyx_TypeName b_name = __Pyx_PyType_GetFullyQualifiedName(b);
            PyErr_Format(PyExc_TypeError,
                "base class '" __Pyx_FMT_TYPENAME "' is not a heap type", b_name);
            __Pyx_DECREF_TypeName(b_name);
#if CYTHON_AVOID_BORROWED_REFS
            Py_DECREF(b0);
#endif
            return -1;
        }
        if (dictoffset == 0)
        {
            Py_ssize_t b_dictoffset = 0;
#if CYTHON_USE_TYPE_SLOTS
            b_dictoffset = b->tp_dictoffset;
#else
            PyObject *py_b_dictoffset = PyObject_GetAttrString((PyObject*)b, "__dictoffset__");
            if (!py_b_dictoffset) goto dictoffset_return;
            b_dictoffset = PyLong_AsSsize_t(py_b_dictoffset);
            Py_DECREF(py_b_dictoffset);
            if (b_dictoffset == -1 && PyErr_Occurred()) goto dictoffset_return;
#endif
            if (b_dictoffset) {
                {
                    __Pyx_TypeName b_name = __Pyx_PyType_GetFullyQualifiedName(b);
                    PyErr_Format(PyExc_TypeError,
                        "extension type '%.200s' has no __dict__ slot, "
                        "but base type '" __Pyx_FMT_TYPENAME "' has: "
                        "either add 'cdef dict __dict__' to the extension type "
                        "or add '__slots__ = [...]' to the base type",
                        type_name, b_name);
                    __Pyx_DECREF_TypeName(b_name);
                }
#if !CYTHON_USE_TYPE_SLOTS
              dictoffset_return:
#endif
#if CYTHON_AVOID_BORROWED_REFS
                Py_DECREF(b0);
#endif
                return -1;
            }
        }
#if CYTHON_AVOID_BORROWED_REFS
        Py_DECREF(b0);
#endif
    }
    return 0;
}
#endif


/////////////// PyType_Ready.proto ///////////////

// unused when using type specs
CYTHON_UNUSED static int __Pyx_PyType_Ready(PyTypeObject *t);/*proto*/

/////////////// PyType_Ready ///////////////
//@requires: ObjectHandling.c::PyObjectCallMethod0
//@requires: ValidateBasesTuple

CYTHON_UNUSED static int __Pyx_PyType_HasMultipleInheritance(PyTypeObject *t) {
    while (t) {
        PyObject *bases = __Pyx_PyType_GetSlot(t, tp_bases, PyObject*);
        if (bases) {
            return 1;
        }
        t = __Pyx_PyType_GetSlot(t, tp_base, PyTypeObject*);
    }
    return 0;
}

// Wrapper around PyType_Ready() with some runtime checks and fixes
// to deal with multiple inheritance.
static int __Pyx_PyType_Ready(PyTypeObject *t) {

#if CYTHON_USE_TYPE_SPECS || !CYTHON_COMPILING_IN_CPYTHON || defined(PYSTON_MAJOR_VERSION)
    // avoid C warning about unused helper function
    (void)__Pyx_PyObject_CallMethod0;
#if CYTHON_USE_TYPE_SPECS
    (void)__Pyx_validate_bases_tuple;
#endif

    return PyType_Ready(t);

#else
    int r;

    if (!__Pyx_PyType_HasMultipleInheritance(t)) {
        // shortcut - if none of the base classes do multiple inheritance then we don't need to
        // (and shouldn't) mess around with faking heaptypes.
        return PyType_Ready(t);
    }
    PyObject *bases = __Pyx_PyType_GetSlot(t, tp_bases, PyObject*);
    if (bases && unlikely(__Pyx_validate_bases_tuple(t->tp_name, t->tp_dictoffset, bases) == -1))
        return -1;

#if !defined(PYSTON_MAJOR_VERSION)
    {
        // Make sure GC does not pick up our non-heap type as heap type with this hack!
        // For details, see https://github.com/cython/cython/issues/3603
        int gc_was_enabled;
    #if PY_VERSION_HEX >= 0x030A00b1
        // finally added in Py3.10 :)
        gc_was_enabled = PyGC_Disable();
        (void)__Pyx_PyObject_CallMethod0;

    #else
        // Call gc.disable() as a backwards compatible fallback, but only if needed.
        PyObject *ret, *py_status;
        PyObject *gc = NULL;
        #if (!CYTHON_COMPILING_IN_PYPY || PYPY_VERSION_NUM+0 >= 0x07030400) && \
                !CYTHON_COMPILING_IN_GRAAL
        // https://foss.heptapod.net/pypy/pypy/-/issues/3385
        gc = PyImport_GetModule(PYUNICODE("gc"));
        #endif
        if (unlikely(!gc)) gc = PyImport_Import(PYUNICODE("gc"));
        if (unlikely(!gc)) return -1;
        py_status = __Pyx_PyObject_CallMethod0(gc, PYUNICODE("isenabled"));
        if (unlikely(!py_status)) {
            Py_DECREF(gc);
            return -1;
        }
        gc_was_enabled = __Pyx_PyObject_IsTrue(py_status);
        Py_DECREF(py_status);
        if (gc_was_enabled > 0) {
            ret = __Pyx_PyObject_CallMethod0(gc, PYUNICODE("disable"));
            if (unlikely(!ret)) {
                Py_DECREF(gc);
                return -1;
            }
            Py_DECREF(ret);
        } else if (unlikely(gc_was_enabled == -1)) {
            Py_DECREF(gc);
            return -1;
        }
    #endif

        // As of https://github.com/python/cpython/issues/66277
        // PyType_Ready enforces that all bases of a non-heap type are
        // non-heap. We know that this is the case for the solid base but
        // other bases are heap allocated and are kept alive through the
        // tp_bases reference.
        // Other than this check, the Py_TPFLAGS_HEAPTYPE flag is unused
        // in PyType_Ready().
        t->tp_flags |= Py_TPFLAGS_HEAPTYPE;
#if PY_VERSION_HEX >= 0x030A0000
        // As of https://github.com/python/cpython/pull/25520
        // PyType_Ready marks types as immutable if they are static types
        // and requires the Py_TPFLAGS_IMMUTABLETYPE flag to mark types as
        // immutable
        // Manually set the Py_TPFLAGS_IMMUTABLETYPE flag, since the type
        // is immutable
        t->tp_flags |= Py_TPFLAGS_IMMUTABLETYPE;
#endif
#else
        // avoid C warning about unused helper function
        (void)__Pyx_PyObject_CallMethod0;
#endif

    r = PyType_Ready(t);

#if !defined(PYSTON_MAJOR_VERSION)
        t->tp_flags &= ~Py_TPFLAGS_HEAPTYPE;

    #if PY_VERSION_HEX >= 0x030A00b1
        if (gc_was_enabled)
            PyGC_Enable();
    #else
        if (gc_was_enabled) {
            PyObject *tp, *v, *tb;
            PyErr_Fetch(&tp, &v, &tb);
            ret = __Pyx_PyObject_CallMethod0(gc, PYUNICODE("enable"));
            if (likely(ret || r == -1)) {
                Py_XDECREF(ret);
                // do not overwrite exceptions raised by PyType_Ready() above
                PyErr_Restore(tp, v, tb);
            } else {
                // PyType_Ready() succeeded, but gc.enable() failed.
                Py_XDECREF(tp);
                Py_XDECREF(v);
                Py_XDECREF(tb);
                r = -1;
            }
        }
        Py_DECREF(gc);
    #endif
    }
#endif

    return r;
#endif
}


/////////////// PyTrashcan.proto ///////////////

// These macros are taken from https://github.com/python/cpython/pull/11841
// Unlike the Py_TRASHCAN_SAFE_BEGIN/Py_TRASHCAN_SAFE_END macros, they
// allow dealing correctly with subclasses.

#if CYTHON_COMPILING_IN_CPYTHON
// https://github.com/python/cpython/pull/11841 merged so Cython reimplementation
// is no longer necessary
#define __Pyx_TRASHCAN_BEGIN Py_TRASHCAN_BEGIN
#define __Pyx_TRASHCAN_END Py_TRASHCAN_END

#elif CYTHON_COMPILING_IN_CPYTHON

#define __Pyx_TRASHCAN_BEGIN_CONDITION(op, cond) \
    do { \
        PyThreadState *_tstate = NULL; \
        // If "cond" is false, then _tstate remains NULL and the deallocator
        // is run normally without involving the trashcan
        if (cond) { \
            _tstate = PyThreadState_GET(); \
            if (_tstate->trash_delete_nesting >= PyTrash_UNWIND_LEVEL) { \
                // Store the object (to be deallocated later) and jump past
                // Py_TRASHCAN_END, skipping the body of the deallocator
                _PyTrash_thread_deposit_object((PyObject*)(op)); \
                break; \
            } \
            ++_tstate->trash_delete_nesting; \
        }
        // The body of the deallocator is here.
#define __Pyx_TRASHCAN_END \
        if (_tstate) { \
            --_tstate->trash_delete_nesting; \
            if (_tstate->trash_delete_later && _tstate->trash_delete_nesting <= 0) \
                _PyTrash_thread_destroy_chain(); \
        } \
    } while (0);

#define __Pyx_TRASHCAN_BEGIN(op, dealloc) __Pyx_TRASHCAN_BEGIN_CONDITION(op, \
        __Pyx_PyObject_GetSlot(op, tp_dealloc, destructor) == (destructor)(dealloc))

#else
// The trashcan is a no-op on other Python implementations
// or old CPython versions
#define __Pyx_TRASHCAN_BEGIN(op, dealloc)
#define __Pyx_TRASHCAN_END
#endif

/////////////// DeallocKeepAlive.proto ///////////////

// Keep `o` alive across a call into user __dealloc__ code, to prevent the user
// code from triggering recursive deallocation, then drop the temporary reference.
//
// On the GIL build this is a plain refcount bump.  On free-threading we must NOT
// use Py_SET_REFCNT(o, 1): the object reaching tp_dealloc is unowned and merged
// (ob_ref_shared == _Py_REF_MERGED), so Py_SET_REFCNT takes the unowned branch and
// writes ob_ref_shared = _Py_REF_SHARED(1, _Py_REF_MERGED), i.e. "merged, shared
// count 1" -- indistinguishable from a live object.  A concurrent
// PyUnstable_TryIncRef() would then succeed and resurrect the dying object
// (use-after-free).  Instead we mirror CPython's _PyObject_ResurrectStart/End:
// re-take ownership and keep the *shared* refcount at the refuse sentinel (0), so
// the bump is invisible to other threads and TryIncRef keeps refusing.

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
#define __Pyx_DeallocKeepAliveBegin(o) do {                                  \
        _Py_atomic_store_uintptr_relaxed(&(o)->ob_tid, _Py_ThreadId());      \
        _Py_atomic_store_uint32_relaxed(&(o)->ob_ref_local, 1);              \
        _Py_atomic_store_ssize_relaxed(&(o)->ob_ref_shared, 0);             \
    } while (0)
#define __Pyx_DeallocKeepAliveEnd(o) \
        _Py_atomic_store_uint32_relaxed(&(o)->ob_ref_local, 0)
#else
#define __Pyx_DeallocKeepAliveBegin(o) __Pyx_SET_REFCNT(o, Py_REFCNT(o) + 1)
#define __Pyx_DeallocKeepAliveEnd(o)   __Pyx_SET_REFCNT(o, Py_REFCNT(o) - 1)
#endif

/////////////// CallNextTpDealloc.proto ///////////////

static void __Pyx_call_next_tp_dealloc(PyObject* obj, destructor current_tp_dealloc);

/////////////// CallNextTpDealloc ///////////////

static void __Pyx_call_next_tp_dealloc(PyObject* obj, destructor current_tp_dealloc) {
    PyTypeObject* type = Py_TYPE(obj);
    destructor tp_dealloc = NULL;
    /* try to find the first parent type that has a different tp_dealloc() function */
    while (type && __Pyx_PyType_GetSlot(type, tp_dealloc, destructor) != current_tp_dealloc)
        type = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
    while (type && (tp_dealloc = __Pyx_PyType_GetSlot(type, tp_dealloc, destructor)) == current_tp_dealloc)
        type = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
    if (type)
        tp_dealloc(obj);
}

/////////////// CallNextTpTraverse.proto ///////////////

static int __Pyx_call_next_tp_traverse(PyObject* obj, visitproc v, void *a, traverseproc current_tp_traverse);

/////////////// CallNextTpTraverse ///////////////

static int __Pyx_call_next_tp_traverse(PyObject* obj, visitproc v, void *a, traverseproc current_tp_traverse) {
    PyTypeObject* type = Py_TYPE(obj);
    traverseproc tp_traverse = NULL;
    /* try to find the first parent type that has a different tp_traverse() function */
    while (type && __Pyx_PyType_GetSlot(type, tp_traverse, traverseproc) != current_tp_traverse)
        type = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
    while (type && (tp_traverse = __Pyx_PyType_GetSlot(type, tp_traverse, traverseproc)) == current_tp_traverse)
        type = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
    if (type && tp_traverse)
        return tp_traverse(obj, v, a);
    // FIXME: really ignore?
    return 0;
}

/////////////// CallNextTpClear.proto ///////////////

static void __Pyx_call_next_tp_clear(PyObject* obj, inquiry current_tp_clear);

/////////////// CallNextTpClear ///////////////

static void __Pyx_call_next_tp_clear(PyObject* obj, inquiry current_tp_clear) {
    PyTypeObject* type = Py_TYPE(obj);
    inquiry tp_clear = NULL;
    /* try to find the first parent type that has a different tp_clear() function */
    while (type && __Pyx_PyType_GetSlot(type, tp_clear, inquiry) != current_tp_clear)
        type = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
    while (type && (tp_clear = __Pyx_PyType_GetSlot(type, tp_clear, inquiry)) == current_tp_clear)
        type = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
    if (type && tp_clear)
        tp_clear(obj);
}


/////////////// SetupReduce.proto ///////////////

static int __Pyx_setup_reduce(PyObject* type_obj);

/////////////// SetupReduce ///////////////
//@requires: ObjectHandling.c::PyObjectGetAttrStrNoError
//@requires: ObjectHandling.c::PyObjectGetAttrStr
//@requires: SetItemOnTypeDict
//@requires: DelItemOnTypeDict

static int __Pyx_setup_reduce_is_named(PyObject* meth, PyObject* name) {
  int ret;
  PyObject *name_attr;

  name_attr = __Pyx_PyObject_GetAttrStrNoError(meth, PYIDENT("__name__"));
  if (likely(name_attr)) {
      ret = PyObject_RichCompareBool(name_attr, name, Py_EQ);
  } else {
      ret = -1;
  }

  if (unlikely(ret < 0)) {
      PyErr_Clear();
      ret = 0;
  }

  Py_XDECREF(name_attr);
  return ret;
}

static int __Pyx_setup_reduce(PyObject* type_obj) {
    int ret = 0;
    PyObject *object_reduce = NULL;
    PyObject *object_getstate = NULL;
    PyObject *object_reduce_ex = NULL;
    PyObject *reduce = NULL;
    PyObject *reduce_ex = NULL;
    PyObject *reduce_cython = NULL;
    PyObject *setstate = NULL;
    PyObject *setstate_cython = NULL;
    PyObject *getstate = NULL;

#if CYTHON_USE_PYTYPE_LOOKUP
    getstate = _PyType_Lookup((PyTypeObject*)type_obj, PYIDENT("__getstate__"));
#else
    getstate = __Pyx_PyObject_GetAttrStrNoError(type_obj, PYIDENT("__getstate__"));
    if (!getstate && PyErr_Occurred()) {
        goto __PYX_BAD;
    }
#endif
    if (getstate) {
        // Python 3.11 introduces object.__getstate__. Because it's version-specific failure to find it should not be an error
#if CYTHON_USE_PYTYPE_LOOKUP
        object_getstate = _PyType_Lookup(&PyBaseObject_Type, PYIDENT("__getstate__"));
#else
        object_getstate = __Pyx_PyObject_GetAttrStrNoError((PyObject*)&PyBaseObject_Type, PYIDENT("__getstate__"));
        if (!object_getstate && PyErr_Occurred()) {
            goto __PYX_BAD;
        }
#endif
        if (object_getstate != getstate) {
            goto __PYX_GOOD;
        }
    }

#if CYTHON_USE_PYTYPE_LOOKUP
    object_reduce_ex = _PyType_Lookup(&PyBaseObject_Type, PYIDENT("__reduce_ex__")); if (!object_reduce_ex) goto __PYX_BAD;
#else
    object_reduce_ex = __Pyx_PyObject_GetAttrStr((PyObject*)&PyBaseObject_Type, PYIDENT("__reduce_ex__")); if (!object_reduce_ex) goto __PYX_BAD;
#endif

    reduce_ex = __Pyx_PyObject_GetAttrStr(type_obj, PYIDENT("__reduce_ex__")); if (unlikely(!reduce_ex)) goto __PYX_BAD;
    if (reduce_ex == object_reduce_ex) {

#if CYTHON_USE_PYTYPE_LOOKUP
        object_reduce = _PyType_Lookup(&PyBaseObject_Type, PYIDENT("__reduce__")); if (!object_reduce) goto __PYX_BAD;
#else
        object_reduce = __Pyx_PyObject_GetAttrStr((PyObject*)&PyBaseObject_Type, PYIDENT("__reduce__")); if (!object_reduce) goto __PYX_BAD;
#endif
        reduce = __Pyx_PyObject_GetAttrStr(type_obj, PYIDENT("__reduce__")); if (unlikely(!reduce)) goto __PYX_BAD;

        if (reduce == object_reduce || __Pyx_setup_reduce_is_named(reduce, PYIDENT("__reduce_cython__"))) {
            reduce_cython = __Pyx_PyObject_GetAttrStrNoError(type_obj, PYIDENT("__reduce_cython__"));
            if (likely(reduce_cython)) {
                ret = __Pyx_SetItemOnTypeDict((PyTypeObject*)type_obj, PYIDENT("__reduce__"), reduce_cython); if (unlikely(ret < 0)) goto __PYX_BAD;
                ret = __Pyx_DelItemOnTypeDict((PyTypeObject*)type_obj, PYIDENT("__reduce_cython__")); if (unlikely(ret < 0)) goto __PYX_BAD;
            } else if (reduce == object_reduce || PyErr_Occurred()) {
                // Ignore if we're done, i.e. if 'reduce' already has the right name and the original is gone.
                // Otherwise: error.
                goto __PYX_BAD;
            }

            setstate = __Pyx_PyObject_GetAttrStrNoError(type_obj, PYIDENT("__setstate__"));
            if (!setstate) PyErr_Clear();
            if (!setstate || __Pyx_setup_reduce_is_named(setstate, PYIDENT("__setstate_cython__"))) {
                setstate_cython = __Pyx_PyObject_GetAttrStrNoError(type_obj, PYIDENT("__setstate_cython__"));
                if (likely(setstate_cython)) {
                    ret = __Pyx_SetItemOnTypeDict((PyTypeObject*)type_obj, PYIDENT("__setstate__"), setstate_cython); if (unlikely(ret < 0)) goto __PYX_BAD;
                    ret = __Pyx_DelItemOnTypeDict((PyTypeObject*)type_obj, PYIDENT("__setstate_cython__")); if (unlikely(ret < 0)) goto __PYX_BAD;
                } else if (!setstate || PyErr_Occurred()) {
                    // Ignore if we're done, i.e. if 'setstate' already has the right name and the original is gone.
                    // Otherwise: error.
                    goto __PYX_BAD;
                }
            }
            PyType_Modified((PyTypeObject*)type_obj);
        }
    }
    goto __PYX_GOOD;

__PYX_BAD:
    if (!PyErr_Occurred()) {
        __Pyx_TypeName type_obj_name =
            __Pyx_PyType_GetFullyQualifiedName((PyTypeObject*)type_obj);
        PyErr_Format(PyExc_RuntimeError,
            "Unable to initialize pickling for " __Pyx_FMT_TYPENAME, type_obj_name);
        __Pyx_DECREF_TypeName(type_obj_name);
    }
    ret = -1;
__PYX_GOOD:
#if !CYTHON_USE_PYTYPE_LOOKUP
    Py_XDECREF(object_reduce);
    Py_XDECREF(object_reduce_ex);
    Py_XDECREF(object_getstate);
    Py_XDECREF(getstate);
#endif
    Py_XDECREF(reduce);
    Py_XDECREF(reduce_ex);
    Py_XDECREF(reduce_cython);
    Py_XDECREF(setstate);
    Py_XDECREF(setstate_cython);
    return ret;
}


/////////////// CheckUnpickleChecksum.proto ///////////////

static CYTHON_INLINE int __Pyx_CheckUnpickleChecksum(long checksum, long checksum1, long checksum2, long checksum3, const char *members); /*proto*/

/////////////// CheckUnpickleChecksum ///////////////

static void __Pyx_RaiseUnpickleChecksumError(long checksum, long checksum1, long checksum2, long checksum3, const char *members) {
    // This function always raises some kind of error, either the expected one or a different one.
    PyObject *pickle_module = PyImport_ImportModule("pickle");
    if (unlikely(!pickle_module)) return;

    PyObject *pickle_error = PyObject_GetAttrString(pickle_module, "PickleError");
    Py_DECREF(pickle_module);
    if (unlikely(!pickle_error)) return;

    if (checksum2 == checksum1) {
        PyErr_Format(pickle_error, "Incompatible checksums (0x%x vs (0x%x) = (%s))",
            checksum, checksum1, members);
    } else if (checksum3 == checksum2) {
        PyErr_Format(pickle_error, "Incompatible checksums (0x%x vs (0x%x, 0x%x) = (%s))",
            checksum, checksum1, checksum2, members);
    } else {
        PyErr_Format(pickle_error, "Incompatible checksums (0x%x vs (0x%x, 0x%x, 0x%x) = (%s))",
            checksum, checksum1, checksum2, checksum3, members);
    }
    Py_DECREF(pickle_error);
}

static int __Pyx_CheckUnpickleChecksum(long checksum, long checksum1, long checksum2, long checksum3, const char *members) {
    int found = 0;
    found |= checksum1 == checksum;
    found |= checksum2 == checksum;
    found |= checksum3 == checksum;
    if (likely(found))
        return 0;

    __Pyx_RaiseUnpickleChecksumError(checksum, checksum1, checksum2, checksum3, members);
    return -1;
}


/////////////// UpdateUnpickledDict.proto ///////////////

static int __Pyx_UpdateUnpickledDict(PyObject *obj, PyObject *state, Py_ssize_t index); /*proto*/

/////////////// UpdateUnpickledDict ///////////////
//@requires: ObjectHandling.c::PyObjectCallMethod1

static int __Pyx__UpdateUnpickledDict(PyObject *obj, PyObject *state, Py_ssize_t index) {
    PyObject *state_dict = __Pyx_PySequence_ITEM(state, index);
    if (unlikely(!state_dict)) {
        return -1;
    }

    int non_empty = PyObject_IsTrue(state_dict);
    if (non_empty == 0) {
        // Nothing to do.
        Py_DECREF(state_dict);
        return 0;
    } else if (unlikely(non_empty == -1)) {
        return -1;
    }

    // Since we received a non-empty dict from pickling,
    // we assume that we're unpickling an object that has one, too.
    PyObject *dict;
    #if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
    dict = PyObject_GetAttrString(obj, "__dict__");
    #else
    dict = PyObject_GenericGetDict(obj, NULL);
    #endif
    if (unlikely(!dict)) {
        // It is debatable if it is a fatal error if we cannot reassign the state of the dict because
        // the unpickled object does not have a '__dict__'. But the user should probably know.
        Py_DECREF(state_dict);
        return -1;
    }

    int result;
    if (likely(PyDict_CheckExact(dict))) {
        result = PyDict_Update(dict, state_dict);
    } else {
        PyObject *obj_result = __Pyx_PyObject_CallMethod1(dict, PYIDENT("update"), state_dict);
        if (likely(obj_result)) {
            Py_DECREF(obj_result);
            result = 0;
        } else {
            result = -1;
        }
    }
    Py_DECREF(state_dict);
    Py_DECREF(dict);
    return result;
}

static int __Pyx_UpdateUnpickledDict(PyObject *obj, PyObject *state, Py_ssize_t index) {
    Py_ssize_t state_size = __Pyx_PyTuple_GET_SIZE(state);
    #if !CYTHON_ASSUME_SAFE_SIZE
    if (unlikely(state_size == -1)) return -1;
    #endif
    if (state_size <= index) {
        // No dict from pickling.
        return 0;
    }

    return __Pyx__UpdateUnpickledDict(obj, state, index);
}


/////////////// BinopSlot ///////////////

static CYTHON_INLINE PyObject *{{func_name}}_maybe_call_slot(PyTypeObject* type, PyObject *left, PyObject *right {{extra_arg_decl}}) {
    {{slot_type}} slot;
#if CYTHON_USE_TYPE_SLOTS
    slot = type->tp_as_number ? type->tp_as_number->{{slot_name}} : NULL;
#else
    slot = ({{slot_type}}) PyType_GetSlot(type, Py_{{slot_name}});
#endif
    return slot ? slot(left, right {{extra_arg}}) : __Pyx_NewRef(Py_NotImplemented);
}

static PyObject *{{func_name}}(PyObject *left, PyObject *right {{extra_arg_decl}}) {
    int maybe_self_is_left, maybe_self_is_right = 0;
    maybe_self_is_left = Py_TYPE(left) == Py_TYPE(right)
#if CYTHON_USE_TYPE_SLOTS
            || (Py_TYPE(left)->tp_as_number && Py_TYPE(left)->tp_as_number->{{slot_name}} == &{{func_name}})
#endif
            || __Pyx_TypeCheck(left, {{type_cname}});

    // Optimize for the common case where the left operation is defined (and successful).
    {{if not overloads_left}}
    maybe_self_is_right = Py_TYPE(left) == Py_TYPE(right)
#if CYTHON_USE_TYPE_SLOTS
            || (Py_TYPE(right)->tp_as_number && Py_TYPE(right)->tp_as_number->{{slot_name}} == &{{func_name}})
#endif
            || __Pyx_TypeCheck(right, {{type_cname}});
    {{endif}}

    if (maybe_self_is_left) {
        PyObject *res;

        {{if overloads_right and not overloads_left}}
        if (maybe_self_is_right) {
            res = {{call_right}};
            if (res != Py_NotImplemented) return res;
            Py_DECREF(res);
            // Don't bother calling it again.
            maybe_self_is_right = 0;
        }
        {{endif}}

        res = {{call_left}};
        if (res != Py_NotImplemented) return res;
        Py_DECREF(res);
    }

    {{if overloads_left}}
    maybe_self_is_right = Py_TYPE(left) == Py_TYPE(right)
#if CYTHON_USE_TYPE_SLOTS
            || (Py_TYPE(right)->tp_as_number && Py_TYPE(right)->tp_as_number->{{slot_name}} == &{{func_name}})
#endif
            || PyType_IsSubtype(Py_TYPE(right), {{type_cname}});
    {{endif}}

    if (maybe_self_is_right) {
        return {{call_right}};
    }
    return __Pyx_NewRef(Py_NotImplemented);
}

/////////////// ValidateExternBase.proto ///////////////

static int __Pyx_validate_extern_base(PyTypeObject *base); /* proto */

/////////////// ValidateExternBase ///////////////
//@requires: ObjectHandling.c::FormatTypeName

static int __Pyx_validate_extern_base(PyTypeObject *base) {
    Py_ssize_t itemsize;
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *py_itemsize;
#endif
#if !CYTHON_COMPILING_IN_LIMITED_API
    itemsize = ((PyTypeObject *)base)->tp_itemsize;
#else
    py_itemsize = PyObject_GetAttrString((PyObject*)base, "__itemsize__");
    if (!py_itemsize)
        return -1;
    itemsize = PyLong_AsSsize_t(py_itemsize);
    Py_DECREF(py_itemsize);
    py_itemsize = 0;
    if (itemsize == (Py_ssize_t)-1 && PyErr_Occurred())
        return -1;
#endif
    if (itemsize) {
        __Pyx_TypeName b_name = __Pyx_PyType_GetFullyQualifiedName(base);
        PyErr_Format(PyExc_TypeError,
                "inheritance from PyVarObject types like '" __Pyx_FMT_TYPENAME "' not currently supported", b_name);
        __Pyx_DECREF_TypeName(b_name);
        return -1;
    }
    return 0;
}

/////////////// CallTypeTraverse.proto /////////////////////

#if !CYTHON_USE_TYPE_SPECS || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX < 0x03090000)
// Without type specs, we're never responsible for this.
// If we *know* the Python version is less that 3.9 we're also not responsible
#define __Pyx_call_type_traverse(o, always_call, visit, arg) 0
#else
static int __Pyx_call_type_traverse(PyObject *o, int always_call, visitproc visit, void *arg); /* proto */
#endif

/////////////// CallTypeTraverse ////////////////////////////

#if !CYTHON_USE_TYPE_SPECS || (!CYTHON_COMPILING_IN_LIMITED_API && PY_VERSION_HEX < 0x03090000)
// nothing to do
#else
static int __Pyx_call_type_traverse(PyObject *o, int always_call, visitproc visit, void *arg) {
    #if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x03090000
    // We have to work out whether to call traverse based on the runtime version.
    // __Pyx_get_runtime_version is always available so no need to require it.
    if (__Pyx_get_runtime_version() < 0x03090000) return 0;
    #endif
    if (!always_call) {
        // Written with reference to https://docs.python.org/3/howto/isolating-extensions.html
        PyTypeObject *base = __Pyx_PyObject_GetSlot(o, tp_base, PyTypeObject*);
        unsigned long flags = PyType_GetFlags(base);
        if (flags & Py_TPFLAGS_HEAPTYPE) {
            // The base class should have handled it
            return 0;
        }
    }
    Py_VISIT((PyObject*)Py_TYPE(o));
    return 0;
}
#endif


////////////////// LimitedApiGetTypeDict.proto //////////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
// This is a little hacky - the Limited API works quite hard to stop us getting
// the dict of a type object. But apparently not hard enough...
//
// In future we should prefer to work with mutable types, and then make them immutable
// once we're done (pending C API support for this).
static PyObject *__Pyx_GetTypeDict(PyTypeObject *tp); /* proto */
#endif

////////////////// LimitedApiGetTypeDict //////////////////////

#if CYTHON_COMPILING_IN_LIMITED_API
static Py_ssize_t __Pyx_GetTypeDictOffset(void) {
    PyObject *tp_dictoffset_o;
    Py_ssize_t tp_dictoffset;
    tp_dictoffset_o = PyObject_GetAttrString((PyObject*)(&PyType_Type), "__dictoffset__");
    if (unlikely(!tp_dictoffset_o)) return -1;
    tp_dictoffset = PyLong_AsSsize_t(tp_dictoffset_o);
    Py_DECREF(tp_dictoffset_o);

    if (unlikely(tp_dictoffset == 0)) {
        PyErr_SetString(
            PyExc_TypeError,
            "'type' doesn't have a dictoffset");
        return -1;
    } else if (unlikely(tp_dictoffset < 0)) {
        // This isn't completely future proof. dictoffset can be
        // negative, but isn't in Python <=3.13 (current at time
        // of writing).  It's awkward to calculate in the limited
        // API because we need to know the object size.  For now
        // just raise an error and fix it if it every changes.
        PyErr_SetString(
            PyExc_TypeError,
            "'type' has an unexpected negative dictoffset. "
            "Please report this as Cython bug");
        return -1;
    }

    return tp_dictoffset;
}

static PyObject *__Pyx_GetTypeDict(PyTypeObject *tp) {
    // TODO - if we ever support custom metatypes for extension types then
    // we have to modify this caching.
    static Py_ssize_t tp_dictoffset = 0;
    if (unlikely(tp_dictoffset == 0)) {
        tp_dictoffset = __Pyx_GetTypeDictOffset();
        // Note that negative dictoffsets are definitely allowed.
        // A dictoffset of -1 seems unlikely but isn't obviously forbidden.
        if (unlikely(tp_dictoffset == -1 && PyErr_Occurred())) {
            tp_dictoffset = 0; // try again next time?
            return NULL;
        }
    }
    return *(PyObject**)((char*)tp + tp_dictoffset);
}
#endif


////////////////// SetItemOnTypeDict.proto //////////////////////////
//@requires: LimitedApiGetTypeDict

static int __Pyx__SetItemOnTypeDict(PyTypeObject *tp, PyObject *k, PyObject *v); /* proto */

#define __Pyx_SetItemOnTypeDict(tp, k, v) __Pyx__SetItemOnTypeDict((PyTypeObject*)tp, k, v)

////////////////// SetItemOnTypeDict //////////////////////////

static int __Pyx__SetItemOnTypeDict(PyTypeObject *tp, PyObject *k, PyObject *v) {
    int result;
    PyObject *tp_dict;
#if CYTHON_COMPILING_IN_LIMITED_API
    tp_dict = __Pyx_GetTypeDict(tp);
    if (unlikely(!tp_dict)) return -1;
#else
    tp_dict = tp->tp_dict;
#endif
    result = PyDict_SetItem(tp_dict, k, v);
    if (likely(!result)) {
        PyType_Modified(tp);
        if (unlikely(PyObject_HasAttr(v, PYIDENT("__set_name__")))) {
            PyObject *setNameResult = PyObject_CallMethodObjArgs(v, PYIDENT("__set_name__"),  (PyObject *) tp, k, NULL);
            if (!setNameResult) return -1;
            Py_DECREF(setNameResult);
        }
    }
    return result;
}

////////////////// DelItemOnTypeDict.proto //////////////////////////

static int __Pyx__DelItemOnTypeDict(PyTypeObject *tp, PyObject *k); /* proto */

#define __Pyx_DelItemOnTypeDict(tp, k) __Pyx__DelItemOnTypeDict((PyTypeObject*)tp, k)

////////////////// DelItemOnTypeDict //////////////////////////
//@requires: LimitedApiGetTypeDict

static int __Pyx__DelItemOnTypeDict(PyTypeObject *tp, PyObject *k) {
    int result;
    PyObject *tp_dict;
#if CYTHON_COMPILING_IN_LIMITED_API
    tp_dict = __Pyx_GetTypeDict(tp);
    if (unlikely(!tp_dict)) return -1;
#else
    tp_dict = tp->tp_dict;
#endif
    result = PyDict_DelItem(tp_dict, k);
    if (likely(!result)) PyType_Modified(tp);
    return result;
}

////////////////// AllocateExtensionType.proto ///////////////////////

static PyObject *__Pyx_AllocateExtensionType(PyTypeObject *t, int is_final); /* proto */

////////////////// AllocateExtensionType ////////////////////////////

static PyObject *__Pyx_AllocateExtensionType(PyTypeObject *t, int is_final) {
    if (is_final || likely(!__Pyx_PyType_HasFeature(t, Py_TPFLAGS_IS_ABSTRACT))) {
        allocfunc alloc_func = __Pyx_PyType_GetSlot(t, tp_alloc, allocfunc);
        return alloc_func(t, 0);
    } else {
        // Call PyBaseObject_Type.tp_new. This is is expected to fail, generating an appropriate
        // error message about allocating an abstract type.
        newfunc tp_new = __Pyx_PyType_TryGetSlot(&PyBaseObject_Type, tp_new, newfunc);
    #if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000
        if (!tp_new) {
            PyObject *new_str = PyUnicode_FromString("__new__");
            if (likely(new_str)) {
                PyObject *o = PyObject_CallMethodObjArgs((PyObject *)&PyBaseObject_Type, new_str, t, NULL);
                Py_DECREF(new_str);
                return o;
            } else
                return NULL;
        } else
    #endif
        return tp_new(t, EMPTY(tuple), 0);
    }
}

///////////////////// CheckTypeForFreelists.proto ////////////////////////

#if CYTHON_USE_FREELISTS
#if CYTHON_USE_TYPE_SPECS
// with CYTHON_USE_TYPE_SPECS we can only reasonably use freelists for an exact type match,
// because it isn't easy to look up basicsize is the limited API, and because everything fails the heap-type check.
#define __PYX_CHECK_FINAL_TYPE_FOR_FREELISTS(t, expected_tp, expected_size) ((int) ((t) == (expected_tp)))
#define __PYX_CHECK_TYPE_FOR_FREELIST_FLAGS  Py_TPFLAGS_IS_ABSTRACT
#else
#define __PYX_CHECK_FINAL_TYPE_FOR_FREELISTS(t, expected_tp, expected_size) ((int) ((t)->tp_basicsize == (expected_size)))
#define __PYX_CHECK_TYPE_FOR_FREELIST_FLAGS  (Py_TPFLAGS_IS_ABSTRACT | Py_TPFLAGS_HEAPTYPE)
#endif

#define __PYX_CHECK_TYPE_FOR_FREELISTS(t, expected_tp, expected_size) \
    (__PYX_CHECK_FINAL_TYPE_FOR_FREELISTS((t), (expected_tp), (expected_size)) & \
     (int) (!__Pyx_PyType_HasFeature((t), __PYX_CHECK_TYPE_FOR_FREELIST_FLAGS)))
#endif