/////////////// Import.proto ///////////////////

static CYTHON_INLINE PyObject *__Pyx_Import(PyObject *name, PyObject *const *imported_names, Py_ssize_t len_imported_names, PyObject *qualname, int level); /*proto*/

/////////////// Import /////////////////////////
//@requires: ImportImpl

static PyObject *__Pyx_Import(PyObject *name, PyObject *const *imported_names, Py_ssize_t len_imported_names, PyObject *qualname, int level) {
    // Forward to the shared utility code with an extra moddict argument (because that's not possible to look up in shared utility code)
    return __Pyx__Import(name, imported_names, len_imported_names, qualname, NAMED_CGLOBAL(moddict_cname), level);
}

/////////////// ImportImpl.export ///////////////

static PyObject *__Pyx__Import(PyObject *name, PyObject *const *imported_names, Py_ssize_t len_imported_names, PyObject *qualname, PyObject *moddict, int level); /*proto*/

/////////////// ImportImpl ///////////////
//@requires: StringTools.c::IncludeStringH
//@requires: Builtins.c::HasAttr
//@requires: ObjectHandling.c::TupleAndListFromArray
//@requires: ObjectHandling.c::PyObjectCallOneArg

static int __Pyx__Import_GetModule(PyObject *qualname, PyObject **module) {
    PyObject *imported_module = PyImport_GetModule(qualname);
    if (unlikely(!imported_module)) {
        *module = NULL;
        if (PyErr_Occurred()) {
            return -1;
        }
        return 0;
    }
    *module = imported_module;
    return 1;
}

static int __Pyx__Import_Lookup(PyObject *qualname, PyObject *const *imported_names, Py_ssize_t len_imported_names, PyObject **module) {
    PyObject *imported_module;
    PyObject *top_level_package_name;
    Py_ssize_t i;
    int status, module_found;
    Py_ssize_t dot_index;

    module_found = __Pyx__Import_GetModule(qualname, &imported_module);
    if (unlikely(!module_found || module_found == -1)) {
        *module = NULL;
        return module_found;
    }

    if (imported_names) {
        for (i = 0; i < len_imported_names; i++) {
            PyObject *imported_name = imported_names[i];
#if __PYX_LIMITED_VERSION_HEX < 0x030d0000
            int has_imported_attribute = PyObject_HasAttr(imported_module, imported_name);
#else
            int has_imported_attribute = PyObject_HasAttrWithError(imported_module, imported_name);
            if (unlikely(has_imported_attribute == -1)) goto error;
#endif
            if (!has_imported_attribute) {
                goto not_found;
            }
        }
        *module = imported_module;
        return 1;
    }

    // Get top-level module
    dot_index = PyUnicode_FindChar(qualname, '.', 0, PY_SSIZE_T_MAX, 1);
    if (dot_index == -1) {
        // No dot, so we already have the top-level module
        *module = imported_module;
        return 1;
    }
    if (unlikely(dot_index == -2)) goto error;
    top_level_package_name = PyUnicode_Substring(qualname, 0, dot_index);
    if (unlikely(!top_level_package_name)) goto error;

    Py_DECREF(imported_module);

    status = __Pyx__Import_GetModule(top_level_package_name, module);
    Py_DECREF(top_level_package_name);
    return status;

error:
    Py_DECREF(imported_module);
    *module = NULL;
    return -1;

not_found:
    Py_DECREF(imported_module);
    *module = NULL;
    return 0;
}

static PyObject *__Pyx__Import(PyObject *name, PyObject *const *imported_names, Py_ssize_t len_imported_names, PyObject *qualname, PyObject *moddict, int level) {
    PyObject *module = 0;
    PyObject *empty_dict = 0;
    PyObject *from_list = 0;
    int module_found;

    if (!qualname) {
        qualname = name;
    }
    module_found = __Pyx__Import_Lookup(qualname, imported_names, len_imported_names, &module);
    if (likely(module_found == 1)) {
        return module;
    } else if (unlikely(module_found == -1)) {
        return NULL;
    }
    empty_dict = PyDict_New();
    if (unlikely(!empty_dict))
        goto bad;
    if (imported_names) {
#if CYTHON_COMPILING_IN_CPYTHON
        from_list = __Pyx_PyList_FromArray(imported_names, len_imported_names);
        if (unlikely(!from_list))
            goto bad;
#else
        from_list = PyList_New(len_imported_names);
        if (unlikely(!from_list)) goto bad;
        for (Py_ssize_t i=0; i<len_imported_names; ++i) {
            if (PyList_SetItem(from_list, i, __Pyx_NewRef(imported_names[i])) < 0) goto bad;
        }
#endif
    }
    if (level == -1) {
        const char* package_sep = strchr(__Pyx_MODULE_NAME, '.');
        if (package_sep != (0)) {
            /* try package relative import first */
            module = PyImport_ImportModuleLevelObject(
                name, moddict, empty_dict, from_list, 1);
            if (unlikely(!module)) {
                if (unlikely(!PyErr_ExceptionMatches(PyExc_ImportError)))
                    goto bad;
                PyErr_Clear();
            }
        }
        level = 0; /* try absolute import on failure */
    }
    if (!module) {
        module = PyImport_ImportModuleLevelObject(
            name, moddict, empty_dict, from_list, level);
    }
bad:
    Py_XDECREF(from_list);
    Py_XDECREF(empty_dict);
    return module;
}


/////////////// ImportFrom.proto ///////////////

static PyObject* __Pyx_ImportFrom(PyObject* module, PyObject* name); /*proto*/

/////////////// ImportFrom ///////////////
//@requires: ObjectHandling.c::PyObjectGetAttrStr

static PyObject* __Pyx_ImportFrom(PyObject* module, PyObject* name) {
    PyObject* value = __Pyx_PyObject_GetAttrStr(module, name);
    if (unlikely(!value) && PyErr_ExceptionMatches(PyExc_AttributeError)) {
        // 'name' may refer to a (sub-)module which has not finished initialization
        // yet, and may not be assigned as an attribute to its parent, so try
        // finding it by full name.
        const char* module_name_str = 0;
        PyObject* module_name = 0;
        PyObject* module_dot = 0;
        PyObject* full_name = 0;
        PyErr_Clear();
        module_name_str = PyModule_GetName(module);
        if (unlikely(!module_name_str)) { goto modbad; }
        module_name = PyUnicode_FromString(module_name_str);
        if (unlikely(!module_name)) { goto modbad; }
        module_dot = PyUnicode_Concat(module_name, PYUNICODE("."));
        if (unlikely(!module_dot)) { goto modbad; }
        full_name = PyUnicode_Concat(module_dot, name);
        if (unlikely(!full_name)) { goto modbad; }
        #if (CYTHON_COMPILING_IN_PYPY && PYPY_VERSION_NUM  < 0x07030400) || \
                CYTHON_COMPILING_IN_GRAAL
        {
            PyObject *modules = PyImport_GetModuleDict();
            if (unlikely(!modules))
                goto modbad;
            value = PyObject_GetItem(modules, full_name);
        }
        #else
        value = PyImport_GetModule(full_name);
        #endif

      modbad:
        Py_XDECREF(full_name);
        Py_XDECREF(module_dot);
        Py_XDECREF(module_name);
    }
    if (unlikely(!value)) {
        PyErr_Format(PyExc_ImportError, "cannot import name %S", name);
    }
    return value;
}


/////////////// ImportStar ///////////////
//@substitute: naming

/* import_all_from is an unexposed function from ceval.c */

static int
__Pyx_import_all_from(PyObject *locals, PyObject *v)
{
    PyObject *all = PyObject_GetAttrString(v, "__all__");
    PyObject *dict, *name, *value;
    int skip_leading_underscores = 0;
    int pos, err;

    if (all == NULL) {
        if (!PyErr_ExceptionMatches(PyExc_AttributeError))
            return -1; /* Unexpected error */
        PyErr_Clear();
        dict = PyObject_GetAttrString(v, "__dict__");
        if (dict == NULL) {
            if (!PyErr_ExceptionMatches(PyExc_AttributeError))
                return -1;
            PyErr_SetString(PyExc_ImportError,
            "from-import-* object has no __dict__ and no __all__");
            return -1;
        }
        all = PyMapping_Keys(dict);
        Py_DECREF(dict);
        if (all == NULL)
            return -1;
        skip_leading_underscores = 1;
    }

    for (pos = 0, err = 0; ; pos++) {
        name = PySequence_GetItem(all, pos);
        if (name == NULL) {
            if (!PyErr_ExceptionMatches(PyExc_IndexError))
                err = -1;
            else
                PyErr_Clear();
            break;
        }
        if (skip_leading_underscores && likely(PyUnicode_Check(name))) {
            Py_ssize_t length = __Pyx_PyUnicode_GET_LENGTH(name);
            #if !CYTHON_ASSUME_SAFE_SIZE
            if (unlikely(length < 0)) {
                Py_DECREF(name);
                return -1;
            }
            #endif
            if (likely(length) && __Pyx_PyUnicode_READ_CHAR(name, 0) == '_') {
                Py_DECREF(name);
                continue;
            }
        }
        value = PyObject_GetAttr(v, name);
        if (value == NULL)
            err = -1;
        else if (PyDict_CheckExact(locals))
            err = PyDict_SetItem(locals, name, value);
        else
            err = PyObject_SetItem(locals, name, value);
        Py_DECREF(name);
        Py_XDECREF(value);
        if (err != 0)
            break;
    }
    Py_DECREF(all);
    return err;
}


static int ${import_star}(PyObject* m) {

    int i;
    int ret = -1;
    const char* s;
    PyObject *locals = 0;
    PyObject *list = 0;
    PyObject *utf8_name = 0;
    PyObject *name;
    PyObject *item;
    PyObject *import_obj;
    Py_ssize_t size;
    ${modulestatetype_cname} *mstate = __Pyx_PyModule_GetState(m);

    locals = PyDict_New();              if (!locals) goto bad;
    if (__Pyx_import_all_from(locals, m) < 0) goto bad;
    list = PyDict_Items(locals);        if (!list) goto bad;

    size = __Pyx_PyList_GET_SIZE(list);
    #if !CYTHON_ASSUME_SAFE_SIZE
    if (size < 0) goto bad;
    #endif
    for(i=0; i<size; i++) {
        import_obj = __Pyx_PyList_GET_ITEM(list, i); if (!import_obj) goto bad;
        name = __Pyx_PyTuple_GET_ITEM(import_obj, 0); if (!name) goto bad;
        item = __Pyx_PyTuple_GET_ITEM(import_obj, 1); if (!item) goto bad;

        utf8_name = PyUnicode_AsUTF8String(name);
        if (!utf8_name) goto bad;
        s = __Pyx_PyBytes_AsString(utf8_name); if (!s) goto bad;
        if (${import_star_set}(mstate, item, name, s) < 0) goto bad;
        Py_DECREF(utf8_name); utf8_name = 0;
    }
    ret = 0;

bad:
    Py_XDECREF(locals);
    Py_XDECREF(list);
    Py_XDECREF(utf8_name);
    return ret;
}


/////////////// SetPackagePathFromImportLib.proto ///////////////

#if !CYTHON_PEP489_MULTI_PHASE_INIT
static int __Pyx_SetPackagePathFromImportLib(PyObject *module_name);
#else
#define __Pyx_SetPackagePathFromImportLib(a) 0
#endif

/////////////// SetPackagePathFromImportLib ///////////////
//@substitute: naming

#if !CYTHON_PEP489_MULTI_PHASE_INIT
static int __Pyx_SetPackagePathFromImportLib(PyObject *module_name) {
    PyObject *importlib, *osmod, *ossep, *parts, *package_path;
    PyObject *file_path = NULL;
    PyObject *item;
    int result;
    PyObject *spec;
    // package_path = [importlib.util.find_spec(module_name).origin.rsplit(os.sep, 1)[0]]
    importlib = PyImport_ImportModule("importlib.util");
    if (unlikely(!importlib))
        goto bad;
    spec = PyObject_CallMethod(importlib, "find_spec", "(O)", module_name);
    Py_DECREF(importlib);
    if (unlikely(!spec))
        goto bad;
    file_path = PyObject_GetAttrString(spec, "origin");
    Py_DECREF(spec);
    if (unlikely(!file_path))
        goto bad;

    if (unlikely(PyObject_SetAttrString($module_cname, "__file__", file_path) < 0))
        goto bad;

    osmod = PyImport_ImportModule("os");
    if (unlikely(!osmod))
        goto bad;
    ossep = PyObject_GetAttrString(osmod, "sep");
    Py_DECREF(osmod);
    if (unlikely(!ossep))
        goto bad;
    parts = PyObject_CallMethod(file_path, "rsplit", "(Oi)", ossep, 1);
    Py_DECREF(file_path); file_path = NULL;
    Py_DECREF(ossep);
    if (unlikely(!parts))
        goto bad;
#if CYTHON_ASSUME_SAFE_MACROS
    package_path = Py_BuildValue("[O]", PyList_GET_ITEM(parts, 0));
#else
    item = PyList_GetItem(parts, 0);
    if (unlikely(!item))
        goto bad;
    package_path = Py_BuildValue("[O]", item);
#endif
    Py_DECREF(parts);
    if (unlikely(!package_path))
        goto bad;
    goto set_path;

bad:
    PyErr_WriteUnraisable(module_name);
    Py_XDECREF(file_path);

    // set an empty path list on failure
    PyErr_Clear();
    package_path = PyList_New(0);
    if (unlikely(!package_path))
        return -1;

set_path:
    result = PyObject_SetAttrString($module_cname, "__path__", package_path);
    Py_DECREF(package_path);
    return result;
}
#endif


/////////////// TypeImport.proto ///////////////
//@substitute: naming

#ifndef __PYX_HAVE_RT_ImportType_proto_$cyversion
#define __PYX_HAVE_RT_ImportType_proto_$cyversion

#if defined (__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#include <stdalign.h>
#endif

#if (defined (__STDC_VERSION__) && __STDC_VERSION__ >= 201112L) || __cplusplus >= 201103L
#define __PYX_GET_STRUCT_ALIGNMENT_$cyversion(s) alignof(s)
#else
// best guess at what the alignment could be since we can't measure it
#define __PYX_GET_STRUCT_ALIGNMENT_$cyversion(s) sizeof(void*)
#endif

enum __Pyx_ImportType_CheckSize_$cyversion {
   __Pyx_ImportType_CheckSize_Error_$cyversion = 0,
   __Pyx_ImportType_CheckSize_Warn_$cyversion = 1,
   __Pyx_ImportType_CheckSize_Ignore_$cyversion = 2
};

static PyTypeObject *__Pyx_ImportType_$cyversion(PyObject* module, const char *module_name, const char *class_name, size_t size, size_t alignment, enum __Pyx_ImportType_CheckSize_$cyversion check_size);  /*proto*/

#endif

/////////////// TypeImport ///////////////
//@substitute: naming

// Note that this goes into headers so CYTHON_COMPILING_IN_LIMITED_API may not be available.

#ifndef __PYX_HAVE_RT_ImportType_$cyversion
#define __PYX_HAVE_RT_ImportType_$cyversion
static PyTypeObject *__Pyx_ImportType_$cyversion(PyObject *module, const char *module_name, const char *class_name,
    size_t size, size_t alignment, enum __Pyx_ImportType_CheckSize_$cyversion check_size)
{
    PyObject *result = 0;
    Py_ssize_t basicsize;
    Py_ssize_t itemsize;
#if defined(Py_LIMITED_API) || (defined(CYTHON_COMPILING_IN_LIMITED_API) && CYTHON_COMPILING_IN_LIMITED_API)
    PyObject *py_basicsize;
    PyObject *py_itemsize;
#endif

    result = PyObject_GetAttrString(module, class_name);
    if (!result)
        goto bad;
    if (!PyType_Check(result)) {
        PyErr_Format(PyExc_TypeError,
            "%.200s.%.200s is not a type object",
            module_name, class_name);
        goto bad;
    }
#if !( defined(Py_LIMITED_API) || (defined(CYTHON_COMPILING_IN_LIMITED_API) && CYTHON_COMPILING_IN_LIMITED_API) )
    basicsize = ((PyTypeObject *)result)->tp_basicsize;
    itemsize = ((PyTypeObject *)result)->tp_itemsize;
#else
    if (size == 0) {
        return (PyTypeObject *)result;
    }
    py_basicsize = PyObject_GetAttrString(result, "__basicsize__");
    if (!py_basicsize)
        goto bad;
    basicsize = PyLong_AsSsize_t(py_basicsize);
    Py_DECREF(py_basicsize);
    py_basicsize = 0;
    if (basicsize == (Py_ssize_t)-1 && PyErr_Occurred())
        goto bad;
    py_itemsize = PyObject_GetAttrString(result, "__itemsize__");
    if (!py_itemsize)
        goto bad;
    itemsize = PyLong_AsSsize_t(py_itemsize);
    Py_DECREF(py_itemsize);
    py_itemsize = 0;
    if (itemsize == (Py_ssize_t)-1 && PyErr_Occurred())
        goto bad;
#endif
    if (itemsize) {
        // If itemsize is smaller than the alignment the struct can end up with some extra
        // padding at the end. In this case we need to work out the maximum size that
        // the padding could be when calculating the range of valid struct sizes.
        if (size % alignment) {
            // if this is true we've probably calculated the alignment wrongly
            // (most likely because alignof isn't available)
            alignment = size % alignment;
        }
        if (itemsize < (Py_ssize_t)alignment)
            itemsize = (Py_ssize_t)alignment;
    }
    if ((size_t)(basicsize + itemsize) < size) {
        PyErr_Format(PyExc_ValueError,
            "%.200s.%.200s size changed, may indicate binary incompatibility. "
            "Expected %zd from C header, got %zd from PyObject",
            module_name, class_name, size, basicsize+itemsize);
        goto bad;
    }
    // varobjects almost have structs  between basicsize and basicsize + itemsize
    // but the struct isn't always one of the two limiting values
    if (check_size == __Pyx_ImportType_CheckSize_Error_$cyversion &&
            ((size_t)basicsize > size || (size_t)(basicsize + itemsize) < size)) {
        PyErr_Format(PyExc_ValueError,
            "%.200s.%.200s size changed, may indicate binary incompatibility. "
            "Expected %zd from C header, got %zd-%zd from PyObject",
            module_name, class_name, size, basicsize, basicsize+itemsize);
        goto bad;
    }
    else if (check_size == __Pyx_ImportType_CheckSize_Warn_$cyversion && (size_t)basicsize > size) {
        if (PyErr_WarnFormat(NULL, 0,
                "%.200s.%.200s size changed, may indicate binary incompatibility. "
                "Expected %zd from C header, got %zd from PyObject",
                module_name, class_name, size, basicsize) < 0) {
            goto bad;
        }
    }
    /* check_size == __Pyx_ImportType_CheckSize_Ignore does not warn nor error */
    return (PyTypeObject *)result;
bad:
    Py_XDECREF(result);
    return NULL;
}
#endif

/////////////// FunctionImport.proto ///////////////
//@substitute: naming

static int __Pyx_ImportFunction_$cyversion(PyObject *module, const char *funcname, void (**f)(void), const char *sig); /*proto*/

/////////////// FunctionImport ///////////////
//@substitute: naming
//@requires: PxdImportShared

#ifndef __PYX_HAVE_RT_ImportFunction_$cyversion
#define __PYX_HAVE_RT_ImportFunction_$cyversion
static int __Pyx_ImportFunction_$cyversion(PyObject *module, const char *funcname, void (**f)(void), const char *sig) {
    union {
        void (*fp)(void);
        void *p;
    } tmp;
    int result = __Pyx_ImportFromPxd_$cyversion(module, funcname, &tmp.p, sig, "function");
    if (result == 0) {
        *f = tmp.fp;
    }
    return result;
}
#endif


/////////////// GetApiDict.proto ///////////////

static PyObject *__Pyx_ApiExport_GetApiDict(void); /*proto*/

/////////////// GetApiDict ///////////////
//@substitute: naming

static PyObject *__Pyx_ApiExport_GetApiDict(void) {
    PyObject *d;
    if (__Pyx_PyDict_GetItemRef(NAMED_CGLOBAL(moddict_cname), PYIDENT("$api_name"), &d) == -1)
        return NULL;
    if (!d) {
        d = PyDict_New();
        if (!d)
            goto bad;
        if (PyObject_SetAttr($module_cname, PYIDENT("$api_name"), d) < 0)
            goto bad;
    }
    return d;

bad:
    Py_XDECREF(d);
    return NULL;
}


/////////////// FunctionExport.proto ///////////////

static int __Pyx_ExportFunction(PyObject *api_dict, const char *name, void (*f)(void), const char *sig); /*proto*/

/////////////// FunctionExport ///////////////

static int __Pyx_ExportFunction(PyObject *api_dict, const char *name, void (*f)(void), const char *sig) {
    PyObject *cobj;
    union {
        void (*fp)(void);
        void *p;
    } tmp;

    tmp.fp = f;
    cobj = PyCapsule_New(tmp.p, sig, 0);
    if (!cobj)
        goto bad;
    if (PyDict_SetItemString(api_dict, name, cobj) < 0)
        goto bad;
    Py_DECREF(cobj);
    return 0;
bad:
    Py_XDECREF(cobj);
    return -1;
}

/////////////// PxdImportShared ///////////////
//@substitute: naming

#ifndef __PYX_HAVE_RT_ImportFromPxd_$cyversion
#define __PYX_HAVE_RT_ImportFromPxd_$cyversion
static int __Pyx_ImportFromPxd_$cyversion(PyObject *module, const char *name, void **p, const char *sig, const char *what) {
    PyObject *d = 0;
    PyObject *cobj = 0;

    d = PyObject_GetAttrString(module, "$api_name");
    if (!d)
        goto bad;
// potentially defined in headers so we can't rely on __PYX_LIMITED_VERSION_HEX
#if (defined(Py_LIMITED_API) && Py_LIMITED_API >= 0x030d0000) || (!defined(Py_LIMITED_API) && PY_VERSION_HEX >= 0x030d0000)
    PyDict_GetItemStringRef(d, name, &cobj);
#else
    cobj = PyDict_GetItemString(d, name);
    Py_XINCREF(cobj);
#endif
    if (!cobj) {
        PyErr_Format(PyExc_ImportError,
            "%.200s does not export expected C %.8s %.200s",
                PyModule_GetName(module), what, name);
        goto bad;
    }
    if (!PyCapsule_IsValid(cobj, sig)) {
        PyErr_Format(PyExc_TypeError,
            "C %.8s %.200s.%.200s has wrong signature (expected %.500s, got %.500s)",
             what, PyModule_GetName(module), name, sig, PyCapsule_GetName(cobj));
        goto bad;
    }
    *p = PyCapsule_GetPointer(cobj, sig);
    if (!(*p))
        goto bad;
    Py_DECREF(d);
    Py_DECREF(cobj);
    return 0;
bad:
    Py_XDECREF(d);
    Py_XDECREF(cobj);
    return -1;
}
#endif

/////////////// VoidPtrImport.proto ///////////////
//@substitute: naming

static int __Pyx_ImportVoidPtr_$cyversion(PyObject *module, const char *name, void **p, const char *sig); /*proto*/

/////////////// VoidPtrImport ///////////////
//@substitute: naming
//@requires: PxdImportShared

#ifndef __PYX_HAVE_RT_ImportVoidPtr_$cyversion
#define __PYX_HAVE_RT_ImportVoidPtr_$cyversion
static int __Pyx_ImportVoidPtr_$cyversion(PyObject *module, const char *name, void **p, const char *sig) {
    return __Pyx_ImportFromPxd_$cyversion(module, name, p, sig, "variable");
}
#endif


/////////////// VoidPtrExport.proto ///////////////

static int __Pyx_ExportVoidPtr(PyObject *api_dict, const char *name, void *p, const char *sig); /*proto*/

/////////////// VoidPtrExport ///////////////

static int __Pyx_ExportVoidPtr(PyObject *api_dict, const char *name, void *p, const char *sig) {
    PyObject *cobj;

    cobj = PyCapsule_New(p, sig, 0);
    if (!cobj)
        goto bad;
    if (PyDict_SetItemString(api_dict, name, cobj) < 0)
        goto bad;
    Py_DECREF(cobj);
    return 0;
bad:
    Py_XDECREF(cobj);
    return -1;
}


/////////////// SetVTable.proto ///////////////

static int __Pyx_SetVtable(PyTypeObject* typeptr , void* vtable); /*proto*/

/////////////// SetVTable ///////////////

static int __Pyx_SetVtable(PyTypeObject *type, void *vtable) {
    PyObject *ob = PyCapsule_New(vtable, 0, 0);
    if (unlikely(!ob))
        goto bad;
#if CYTHON_COMPILING_IN_LIMITED_API
    if (unlikely(PyObject_SetAttr((PyObject *) type, PYIDENT("__pyx_vtable__"), ob) < 0))
#else
    if (unlikely(PyDict_SetItem(type->tp_dict, PYIDENT("__pyx_vtable__"), ob) < 0))
#endif
        goto bad;
    Py_DECREF(ob);
    return 0;
bad:
    Py_XDECREF(ob);
    return -1;
}


/////////////// GetVTable.proto ///////////////

static void* __Pyx_GetVtable(PyTypeObject *type); /*proto*/

/////////////// GetVTable ///////////////

static void* __Pyx_GetVtable(PyTypeObject *type) {
    void* ptr;
#if CYTHON_COMPILING_IN_LIMITED_API
    PyObject *ob = PyObject_GetAttr((PyObject *)type, PYIDENT("__pyx_vtable__"));
#else
    PyObject *ob = PyObject_GetItem(type->tp_dict, PYIDENT("__pyx_vtable__"));
#endif
    if (!ob)
        goto bad;
    ptr = PyCapsule_GetPointer(ob, 0);
    if (!ptr && !PyErr_Occurred())
        PyErr_SetString(PyExc_RuntimeError, "invalid vtable found for imported type");
    Py_DECREF(ob);
    return ptr;
bad:
    Py_XDECREF(ob);
    return NULL;
}


/////////////// MergeVTables.proto ///////////////
//@requires: GetVTable

static int __Pyx_MergeVtables(PyTypeObject *type); /*proto*/

/////////////// MergeVTables ///////////////

static int __Pyx_MergeVtables(PyTypeObject *type) {
    int i=0;
    Py_ssize_t size;
    void** base_vtables;
    __Pyx_TypeName tp_base_name = NULL;
    __Pyx_TypeName base_name = NULL;
    void* unknown = (void*)-1;
    PyObject* bases = __Pyx_PyType_GetSlot(type, tp_bases, PyObject*);
    int base_depth = 0;
    {
        PyTypeObject* base = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
        while (base) {
            base_depth += 1;
            base = __Pyx_PyType_GetSlot(base, tp_base, PyTypeObject*);
        }
    }
    base_vtables = (void**) PyMem_Malloc(sizeof(void*) * (size_t)(base_depth + 1));
    base_vtables[0] = unknown;
    // Could do MRO resolution of individual methods in the future, assuming
    // compatible vtables, but for now simply require a common vtable base.
    // Note that if the vtables of various bases are extended separately,
    // resolution isn't possible and we must reject it just as when the
    // instance struct is so extended.  (It would be good to also do this
    // check when a multiple-base class is created in pure Python as well.)
#if CYTHON_COMPILING_IN_LIMITED_API
    size = PyTuple_Size(bases);
    if (size < 0) goto other_failure;
#else
    size = PyTuple_GET_SIZE(bases);
#endif
    for (i = 1; i < size; i++) {
        PyObject *basei;
        void* base_vtable;
#if CYTHON_AVOID_BORROWED_REFS
        basei = PySequence_GetItem(bases, i);
        if (unlikely(!basei)) goto other_failure;
#elif !CYTHON_ASSUME_SAFE_MACROS
        basei = PyTuple_GetItem(bases, i);
        if (unlikely(!basei)) goto other_failure;
#else
        basei = PyTuple_GET_ITEM(bases, i);
#endif
        base_vtable = __Pyx_GetVtable((PyTypeObject*)basei);
#if CYTHON_AVOID_BORROWED_REFS
        Py_DECREF(basei);
#endif
        if (base_vtable != NULL) {
            int j;
            PyTypeObject* base = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
            for (j = 0; j < base_depth; j++) {
                if (base_vtables[j] == unknown) {
                    base_vtables[j] = __Pyx_GetVtable(base);
                    base_vtables[j + 1] = unknown;
                }
                if (base_vtables[j] == base_vtable) {
                    break;
                } else if (base_vtables[j] == NULL) {
                    // No more potential matching bases (with vtables).
                    goto bad;
                }
                base = __Pyx_PyType_GetSlot(base, tp_base, PyTypeObject*);
            }
        }
    }
    PyErr_Clear();
    PyMem_Free(base_vtables);
    return 0;
bad:
    {
        PyTypeObject* basei = NULL;
        PyTypeObject* tp_base = __Pyx_PyType_GetSlot(type, tp_base, PyTypeObject*);
        tp_base_name = __Pyx_PyType_GetFullyQualifiedName(tp_base);
#if CYTHON_AVOID_BORROWED_REFS
        basei = (PyTypeObject*)PySequence_GetItem(bases, i);
        if (unlikely(!basei)) goto really_bad;
#elif !CYTHON_ASSUME_SAFE_MACROS
        basei = (PyTypeObject*)PyTuple_GetItem(bases, i);
        if (unlikely(!basei)) goto really_bad;
#else
        basei = (PyTypeObject*)PyTuple_GET_ITEM(bases, i);
#endif
        base_name = __Pyx_PyType_GetFullyQualifiedName(basei);
#if CYTHON_AVOID_BORROWED_REFS
        Py_DECREF(basei);
#endif
    }
    PyErr_Format(PyExc_TypeError,
        "multiple bases have vtable conflict: '" __Pyx_FMT_TYPENAME "' and '" __Pyx_FMT_TYPENAME "'", tp_base_name, base_name);
#if CYTHON_AVOID_BORROWED_REFS || !CYTHON_ASSUME_SAFE_MACROS
really_bad: // bad has failed!
#endif
    __Pyx_DECREF_TypeName(tp_base_name);
    __Pyx_DECREF_TypeName(base_name);
#if CYTHON_COMPILING_IN_LIMITED_API || CYTHON_AVOID_BORROWED_REFS || !CYTHON_ASSUME_SAFE_MACROS
other_failure:
#endif
    PyMem_Free(base_vtables);
    return -1;
}

/////////////// ImportNumPyArray.module_state_decls //////////////
//@requires: Synchronization.c::Atomics

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
__pyx_atomic_ptr_type __pyx_numpy_ndarray;
#else
// If freethreading but not atomics, then this is guarded by ndarray_mutex
// in __Pyx__ImportNumPyArrayTypeIfAvailable
PyObject *__pyx_numpy_ndarray;
#endif

/////////////// ImportNumPyArray.proto ///////////////

static PyObject* __Pyx_ImportNumPyArrayTypeIfAvailable(void); /*proto*/

/////////////// ImportNumPyArray.cleanup ///////////////

// I'm not actually worried about thread-safety in the cleanup function.
// The CYTHON_ATOMICS code is only for the type-casting.
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
{
    PyObject* old = (PyObject*)__pyx_atomic_pointer_exchange(&CGLOBAL(__pyx_numpy_ndarray), NULL);
    Py_XDECREF(old);
}
#else
Py_CLEAR(CGLOBAL(__pyx_numpy_ndarray));
#endif

/////////////// ImportNumPyArray ///////////////
//@requires: ImportExport.c::Import

static PyObject* __Pyx__ImportNumPyArray(void) {
    PyObject *numpy_module, *ndarray_object = NULL;
    numpy_module = __Pyx_Import(PYIDENT("numpy"), NULL, 0, NULL, 0);
    if (likely(numpy_module)) {
        ndarray_object = PyObject_GetAttrString(numpy_module, "ndarray");
        Py_DECREF(numpy_module);
    }
    if (unlikely(!ndarray_object)) {
        // ImportError, AttributeError, ...
        PyErr_Clear();
    }
    if (unlikely(!ndarray_object || !PyObject_TypeCheck(ndarray_object, &PyType_Type))) {
        Py_XDECREF(ndarray_object);
        Py_INCREF(Py_None);
        ndarray_object = Py_None;
    }
    return ndarray_object;
}

// Returns a borrowed reference, and encapsulates all thread safety stuff needed
static CYTHON_INLINE PyObject* __Pyx__ImportNumPyArrayTypeIfAvailable(void) {
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    static PyMutex ndarray_mutex = {0};
#endif
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
    PyObject *result = (PyObject*)__pyx_atomic_pointer_load_relaxed(&CGLOBAL(__pyx_numpy_ndarray));
    if (unlikely(!result)) {
        PyMutex_Lock(&ndarray_mutex);
        // Now we've got the lock and know that no-one else is modifying it, check again
        // that it hasn't already been set.
        result = (PyObject*)__pyx_atomic_pointer_load_acquire(&CGLOBAL(__pyx_numpy_ndarray));
        if (!result) {
            result = __Pyx__ImportNumPyArray();
            __pyx_atomic_pointer_exchange(&CGLOBAL(__pyx_numpy_ndarray), result);
        }
        PyMutex_Unlock(&ndarray_mutex);
    }
    return result;
#else
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING // but not atomics
    PyMutex_Lock(&ndarray_mutex);
#endif
    if (unlikely(!CGLOBAL(__pyx_numpy_ndarray)))
    {
        CGLOBAL(__pyx_numpy_ndarray) = __Pyx__ImportNumPyArray();
    }
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING // but not atomics
    PyMutex_Unlock(&ndarray_mutex);
#endif
    return CGLOBAL(__pyx_numpy_ndarray);
#endif
}

static CYTHON_INLINE PyObject* __Pyx_ImportNumPyArrayTypeIfAvailable(void) {
    PyObject *result = __Pyx__ImportNumPyArrayTypeIfAvailable();
    Py_INCREF(result);
    return result;
}
