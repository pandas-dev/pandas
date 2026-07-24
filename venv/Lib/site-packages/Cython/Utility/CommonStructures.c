/////////////// FetchSharedCythonModule.proto ///////

static PyObject *__Pyx_FetchSharedCythonABIModule(void);

/////////////// FetchSharedCythonModule ////////////
//@requires: ModuleSetupCode.c::AddModuleRef

static PyObject *__Pyx_FetchSharedCythonABIModule(void) {
    return __Pyx_PyImport_AddModuleRef(__PYX_ABI_MODULE_NAME);
}

/////////////// FetchCommonType.proto ///////////////

static PyTypeObject* __Pyx_FetchCommonTypeFromSpec(PyTypeObject *metaclass, PyObject *module, PyType_Spec *spec, PyObject *bases);

/////////////// FetchCommonType ///////////////
//@requires:ExtensionTypes.c::FixUpExtensionType
//@requires: FetchSharedCythonModule
//@requires:StringTools.c::IncludeStringH
//@requires:Builtins.c::dict_setdefault

#if __PYX_LIMITED_VERSION_HEX < 0x030C0000
static PyObject* __Pyx_PyType_FromMetaclass(PyTypeObject *metaclass, PyObject *module, PyType_Spec *spec, PyObject *bases) {
    PyObject *result = __Pyx_PyType_FromModuleAndSpec(module, spec, bases);
    if (result && metaclass) {
        PyObject *old_tp = (PyObject*)Py_TYPE(result);
    Py_INCREF((PyObject*)metaclass);
#if __PYX_LIMITED_VERSION_HEX >= 0x03090000
        Py_SET_TYPE(result, metaclass);
#else
        result->ob_type = metaclass;
#endif
        Py_DECREF(old_tp);
    }
    return result;
}
#else
#define __Pyx_PyType_FromMetaclass(me, mo, s, b) PyType_FromMetaclass(me, mo, s, b)
#endif

static int __Pyx_VerifyCachedType(PyObject *cached_type,
                               const char *name,
                               Py_ssize_t expected_basicsize) {
    Py_ssize_t basicsize;

    if (!PyType_Check(cached_type)) {
        PyErr_Format(PyExc_TypeError,
            "Shared Cython type %.200s is not a type object", name);
        return -1;
    }

    if (expected_basicsize == 0) {
        return 0; // size is inherited, nothing useful to check
    }

#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *py_basicsize;
    py_basicsize = PyObject_GetAttrString(cached_type, "__basicsize__");
    if (unlikely(!py_basicsize)) return -1;
    basicsize = PyLong_AsSsize_t(py_basicsize);
    Py_DECREF(py_basicsize);
    py_basicsize = NULL;
    if (unlikely(basicsize == (Py_ssize_t)-1) && PyErr_Occurred()) return -1;
#else
    basicsize = ((PyTypeObject*) cached_type)->tp_basicsize;
#endif

    if (basicsize != expected_basicsize) {
        PyErr_Format(PyExc_TypeError,
            "Shared Cython type %.200s has the wrong size, try recompiling",
            name);
        return -1;
    }
    return 0;
}

static PyTypeObject *__Pyx_FetchCommonTypeFromSpec(PyTypeObject *metaclass, PyObject *module, PyType_Spec *spec, PyObject *bases) {
    PyObject *abi_module = NULL, *cached_type = NULL, *abi_module_dict, *new_cached_type, *py_object_name;
    int get_item_ref_result;
    // get the final part of the object name (after the last dot)
    const char* object_name = strrchr(spec->name, '.');
    object_name = object_name ? object_name+1 : spec->name;

    py_object_name = PyUnicode_FromString(object_name);
    if (!py_object_name) return NULL;

    abi_module = __Pyx_FetchSharedCythonABIModule();
    if (!abi_module) goto done;
    abi_module_dict = PyModule_GetDict(abi_module);
    if (!abi_module_dict) goto done;


    get_item_ref_result = __Pyx_PyDict_GetItemRef(abi_module_dict, py_object_name, &cached_type);
    if (get_item_ref_result == 1) {
        if (__Pyx_VerifyCachedType(
              cached_type,
              object_name,
              spec->basicsize) < 0) {
            goto bad;
        }
        goto done;
    } else if (unlikely(get_item_ref_result == -1)) {
        goto bad;
    }

    // Without module-state, pass the ABI module reference to avoid keeping the user module alive by foreign type usages.
    // With module-state it's important to keep the user module alive though.
    cached_type = __Pyx_PyType_FromMetaclass(
        metaclass,
        CYTHON_USE_MODULE_STATE ? module : abi_module,
        spec, bases);
    if (unlikely(!cached_type)) goto bad;
    if (unlikely(__Pyx_fix_up_extension_type_from_spec(spec, (PyTypeObject *) cached_type) < 0)) goto bad;

    new_cached_type = __Pyx_PyDict_SetDefault(abi_module_dict, py_object_name, cached_type);
    if (unlikely(new_cached_type != cached_type)) {
        if (unlikely(!new_cached_type)) goto bad;
        // race to initialize it - use the value that's already been set.
        Py_DECREF(cached_type);
        cached_type = new_cached_type;

        if (__Pyx_VerifyCachedType(
                cached_type,
                object_name,
                spec->basicsize) < 0) {
            goto bad;
        }
        goto done;
    } else {
        Py_DECREF(new_cached_type);
    }

done:
    Py_XDECREF(abi_module);
    Py_DECREF(py_object_name);
    // NOTE: always returns owned reference, or NULL on error
    assert(cached_type == NULL || PyType_Check(cached_type));
    return (PyTypeObject *) cached_type;

bad:
    Py_XDECREF(cached_type);
    cached_type = NULL;
    goto done;
}

//////////////////// CommonTypesMetaclass.module_state_decls ////////////////////

PyTypeObject *__pyx_CommonTypesMetaclassType;

//////////////////// CommonTypesMetaclass.module_state_traverse ///////////////////

Py_VISIT(traverse_module_state->__pyx_CommonTypesMetaclassType);

//////////////////// CommonTypesMetaclass.module_state_clear ///////////////////

Py_CLEAR(clear_module_state->__pyx_CommonTypesMetaclassType);

//////////////////// CommonTypesMetaclass.init //////////////////
//@substitute: naming

if (likely(__pyx_CommonTypesMetaclass_init($module_cname) == 0)); else

/////////////////////////// CommonTypesMetaclass.proto ////////////////////////

static int __pyx_CommonTypesMetaclass_init(PyObject *module); /* proto */

#define __Pyx_CommonTypesMetaclass_USED

//////////////////// CommonTypesMetaclass ////////////////////
//@requires: FetchCommonType
//@substitute: naming

static PyObject* __pyx_CommonTypesMetaclass_get_module(CYTHON_UNUSED PyObject *self, CYTHON_UNUSED void* context) {
    return PyUnicode_FromString(__PYX_ABI_MODULE_NAME);
}

#if __PYX_LIMITED_VERSION_HEX < 0x030A0000
static PyObject* __pyx_CommonTypesMetaclass_call(CYTHON_UNUSED PyObject *self, CYTHON_UNUSED PyObject *args, CYTHON_UNUSED PyObject *kwds) {
    PyErr_SetString(PyExc_TypeError, "Cannot instantiate Cython internal types");
    return NULL;
}

static int __pyx_CommonTypesMetaclass_setattr(CYTHON_UNUSED PyObject *self, CYTHON_UNUSED PyObject *attr, CYTHON_UNUSED PyObject *value) {
    PyErr_SetString(PyExc_TypeError, "Cython internal types are immutable");
    return -1;
}
#endif

static PyGetSetDef __pyx_CommonTypesMetaclass_getset[] = {
    {"__module__", __pyx_CommonTypesMetaclass_get_module, NULL, NULL, NULL},
    {0, 0, 0, 0, 0}
};

static PyType_Slot __pyx_CommonTypesMetaclass_slots[] = {
    {Py_tp_getset, (void *)__pyx_CommonTypesMetaclass_getset},
    #if __PYX_LIMITED_VERSION_HEX < 0x030A0000
    {Py_tp_call, (void*)__pyx_CommonTypesMetaclass_call},
    {Py_tp_new, (void*)__pyx_CommonTypesMetaclass_call},
    {Py_tp_setattro, (void*)__pyx_CommonTypesMetaclass_setattr},
    #endif
    {0, 0}
};

static PyType_Spec __pyx_CommonTypesMetaclass_spec = {
    __PYX_TYPE_MODULE_PREFIX "_common_types_metatype",
    0,
    0,
    Py_TPFLAGS_IMMUTABLETYPE |
    Py_TPFLAGS_DISALLOW_INSTANTIATION |
    Py_TPFLAGS_DEFAULT, /*tp_flags*/
    __pyx_CommonTypesMetaclass_slots
};

static int __pyx_CommonTypesMetaclass_init(PyObject *module) {
    $modulestatetype_cname *mstate = __Pyx_PyModule_GetState(module);
    PyObject *bases = PyTuple_Pack(1, &PyType_Type);
    if (unlikely(!bases)) {
        return -1;
    }
    mstate->__pyx_CommonTypesMetaclassType = __Pyx_FetchCommonTypeFromSpec(NULL, module, &__pyx_CommonTypesMetaclass_spec, bases);
    Py_DECREF(bases);
    if (unlikely(mstate->__pyx_CommonTypesMetaclassType == NULL)) {
        return -1;
    }
    return 0;
}
