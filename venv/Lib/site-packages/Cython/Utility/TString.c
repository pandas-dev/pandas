//////////////////////////// InitializeTemplateLib.module_state_decls /////////////////////
//@requires: Synchronization.c::Atomics

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
__pyx_atomic_ptr_type __pyx_templatelib_Template;
__pyx_atomic_ptr_type __pyx_templatelib_Interpolation;
#else
// If freethreading but not atomics then this is unguarded
PyObject *__pyx_templatelib_Template;
PyObject *__pyx_templatelib_Interpolation;
#endif

//////////////////////////// InitializeTemplateLib.module_state_traverse ////////////////////////////

Py_VISIT((PyObject*)traverse_module_state->__pyx_templatelib_Template);
Py_VISIT((PyObject*)traverse_module_state->__pyx_templatelib_Interpolation);

//////////////////////////// InitializeTemplateLib.module_state_clear ////////////////////////////

Py_XDECREF((PyObject*)clear_module_state->__pyx_templatelib_Template);
clear_module_state->__pyx_templatelib_Template = 0;
Py_XDECREF((PyObject*)clear_module_state->__pyx_templatelib_Interpolation);
clear_module_state->__pyx_templatelib_Interpolation = 0;

//////////////////////////// InitializeTemplateLib.proto ///////////////////////////

// Returns Template if template, else Interpolation
static PyObject* __Pyx__GetObjectFromTemplateLib(int is_template); /* proto */

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && !CYTHON_ATOMICS
static PyObject* __Pyx_GetObjectFromTemplateLib(int is_template); /* proto */
#else
#define __Pyx_GetObjectFromTemplateLib __Pyx__GetObjectFromTemplateLib
#endif


//////////////////////////// InitializeTemplateLib ///////////////////////////

#if __PYX_LIMITED_VERSION_HEX < 0x030E0000
static PyObject *__Pyx_TemplateLibFallback(void) {
    PyErr_Clear();

    // The assumption here is that Interpolation and Template are fairly simple classes
    // and the cost of compiling them with Cython (for all Python versions) is probably
    // higher than the cost of using a plain-Python fallback. 
    const char code_str[] = CSTRING("""
class Interpolation:
    __module__ = 'string.templatelib'
    __slots__ = __match_args__ = ('value', 'expression', 'conversion', 'format_spec')
    def __setattr__(self, attr, value):
        raise AttributeError('Interpolation is immutable')
    def __new__(cls, value, expression='', conversion=None, format_spec=''):
        obj = super().__new__(cls)
        super().__setattr__(obj, 'value', value)
        super().__setattr__(obj, 'expression', expression)
        super().__setattr__(obj, 'conversion', conversion)
        super().__setattr__(obj, 'format_spec', format_spec)
        return obj
    def __repr__(self):
        return f'Interpolation({self.value!r}, {self.expression!r}, {self.conversion!r}, {self.format_spec!r})'
    def __reduce__(self):
        # This probably won't work unless a t-string has already been created from Cython
        return (type(self), (self.value, self.expression, self.conversion, self.format_spec))
    def __init_subclass__(cls, **kwds):
        raise TypeError('Interpolation is not an acceptable base type')
class Template:
    __module__ = 'string.templatelib'
    __slots__ = ('strings', 'interpolations')
    def __setattr__(self, attr, value):
        raise AttributeError('Template is immutable')
    def __new__(cls, *args, strings=None, interpolations=None):
        if strings is None and interpolations is None:
            strings = []
            interpolations = []
            last_string = ''
            for arg in args:
                if isinstance(arg, str):
                    last_string += arg
                elif isinstance(arg, Interpolation):
                    strings.append(last_string)
                    last_string = ''
                    interpolations.append(arg)
                else:
                    raise TypeError('Unexpected argument to Template')
            strings.append(last_string)
        else:
            if args:
                raise ValueError("'strings' or 'interpolations' should not be passed with positional arguments")
            if strings is None: strings = ()
            if interpolations is None: interpolations = ()
        obj = super().__new__(cls)
        super().__setattr__(obj, 'strings', tuple(strings))
        super().__setattr__(obj, 'interpolations', tuple(interpolations))
        return obj
    def __repr__(self):
        return f'Template(strings={self.strings!r}, interpolations={self.interpolations!r})'
    def __reduce__(self):
        # This probably won't work unless a t-string has already been created from Cython.
        # It also doesn't quite match how CPython pickles them.
        values = tuple(iter(self))
        return (type(self), values)
    def __iter__(self):
        for n in range(len(self.interpolations)):
            if (s := self.strings[n]):
               yield s
            yield self.interpolations[n]
        if (s := self.strings[-1]):
           yield s
    def __add__(self, other):
        if not (isinstance(self, Template) and isinstance(other, Template)):
            raise TypeError('can only concatenate Template to Template')
        interpolations = self.interpolations + other.interpolations
        middle_string = self.strings[-1] + other.strings[0]
        strings = self.strings[:-1] + (middle_string,) + other.strings[1:]
        return Template(strings=strings, interpolations=interpolations)
    @property
    def values(self):
        return tuple(i.value for i in self.interpolations)
    def __init_subclass__(cls, **kwds):
        raise TypeError('Template is not an acceptable base type')
""");

    PyObject *code=NULL, *eval_result=NULL, *module=NULL, *module_dict=NULL;
    PyObject *dict = PyDict_New();
    if (unlikely(!dict)) return NULL;
#if __PYX_LIMITED_VERSION_HEX < 0x030A0000
    {
        PyObject *builtins = PyEval_GetBuiltins();
        if (unlikely(!builtins)) goto end;
        if (unlikely(PyDict_SetItemString(dict, "__builtins__", builtins) < 0)) goto end;
    }
#endif
    
    code = Py_CompileString(code_str, "<cython string.templatelib fallback>", Py_file_input);
    if (unlikely(!code)) goto end;
    eval_result = PyEval_EvalCode(code, dict, NULL);
    Py_DECREF(code);
    if (unlikely(!eval_result)) goto end;
    Py_DECREF(eval_result);

    module = __Pyx_PyImport_AddModuleRef("string.templatelib");
    if (!module) goto end;
    module_dict = PyModule_GetDict(module);
    if (!module_dict) goto bad;
    if (unlikely(PyDict_Merge(module_dict, dict, 0) < 0)) goto bad;

  end:
    Py_XDECREF(dict);
    return module;
  bad:
    Py_CLEAR(module);
    goto end;
}
#endif

static int __Pyx_InitializeTemplateLib(void) {
    // Even in earlier versions of Python, still try the import. We're happy
    // to use what's there if someone's patched it with something compatible.
    PyObject *templatelib = PyImport_ImportModule("string.templatelib");
    if (!templatelib) {
#if __PYX_LIMITED_VERSION_HEX < 0x030E0000
        templatelib = __Pyx_TemplateLibFallback();
        if (!templatelib)
#endif
        return -1;
    }
    PyObject *template_=NULL, *interpolation=NULL;
    int result = -1;
    template_ = PyObject_GetAttrString(templatelib, "Template");
    if (unlikely(!template_)) goto end;
    interpolation = PyObject_GetAttrString(templatelib, "Interpolation");
    if (unlikely(!interpolation)) goto end;

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
    {
        __pyx_nonatomic_ptr_type expected = 0;
        if (!__pyx_atomic_pointer_cmp_exchange(&CGLOBAL(__pyx_templatelib_Template), &expected, (__pyx_nonatomic_ptr_type)template_)) {
            // Already written - that's fine.
            Py_DECREF(template_);
        }
        expected = 0;
        if (!__pyx_atomic_pointer_cmp_exchange(&CGLOBAL(__pyx_templatelib_Interpolation), &expected, (__pyx_nonatomic_ptr_type)interpolation)) {
            // Already written - that's fine.
            Py_DECREF(interpolation);
        }
    }
#else
    if (unlikely(CGLOBAL(__pyx_templatelib_Template))) {
        Py_DECREF(template_);
    } else {
        CGLOBAL(__pyx_templatelib_Template) = template_;
    }
    if (unlikely(CGLOBAL(__pyx_templatelib_Interpolation))) {
        Py_DECREF(interpolation);
    } else {
        CGLOBAL(__pyx_templatelib_Interpolation) = interpolation;
    }
#endif
    result = 0;

  end:
    Py_DECREF(templatelib);
    Py_XDECREF(template_);
    Py_XDECREF(interpolation);
    return result;
}

static PyObject* __Pyx__GetObjectFromTemplateLib(int is_template) {
    PyObject *lookup;
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
    __pyx_atomic_ptr_type* ptr;
#else
    PyObject **ptr;
#endif
    ptr = is_template ? &CGLOBAL(__pyx_templatelib_Template) : &CGLOBAL(__pyx_templatelib_Interpolation);
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
    lookup = (PyObject*)__pyx_atomic_pointer_load_relaxed(ptr);
#else
    lookup = *ptr;
#endif
    if (unlikely(!lookup)) {
        if (unlikely(__Pyx_InitializeTemplateLib()) < 0) return NULL;
#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && CYTHON_ATOMICS
        lookup = (PyObject*)__pyx_atomic_pointer_load_acquire(ptr);
#else
        lookup = *ptr;
#endif
    }
    Py_XINCREF(lookup);
    return lookup;
}

#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING && !CYTHON_ATOMICS
static PyObject* __Pyx_GetObjectFromTemplateLib(int is_template) {
    static PyMutex mutex = {0};
    PyMutex_Lock(&mutex);
    PyObject *result = __Pyx__GetObjectFromTemplateLib(is_template);
    PyMutex_Unlock(&mutex);
    return result;
}
#endif




//////////////////////////// MakeTemplateLibInterpolation.proto //////////////////////

static PyObject* __Pyx_MakeTemplateLibInterpolation(PyObject *value, PyObject *expression, PyObject *conversion_char, PyObject *format_spec); /* proto */

//////////////////////////// MakeTemplateLibInterpolation ////////////////////////
//@requires: InitializeTemplateLib
//@requires: ObjectHandling.c::PyObjectFastCall

static PyObject* __Pyx_MakeTemplateLibInterpolation(PyObject *value, PyObject *expression, PyObject *conversion_char, PyObject *format_spec) {
    PyObject *tp = __Pyx_GetObjectFromTemplateLib(0);
    if (unlikely(!tp)) return NULL;
    PyObject *args[] = {value, expression, conversion_char, format_spec};

    PyObject *result = __Pyx_PyObject_FastCallDict(tp, args, 4, NULL);
    Py_DECREF(tp);

    return result;
}

//////////////////////////// MakeTemplateLibTemplate.proto //////////////////////

static PyObject* __Pyx_MakeTemplateLibTemplate(PyObject *strings, PyObject *interpolations); /* proto */

//////////////////////////// MakeTemplateLibTemplate ////////////////////////
//@requires: InitializeTemplateLib
//@requires: ObjectHandling.c::PyObjectVectorCallKwBuilder

#if PY_VERSION_HEX >= 0x030E0000 && CYTHON_COMPILING_IN_CPYTHON
#ifndef Py_BUILD_CORE
#define Py_BUILD_CORE
#endif
#include "internal/pycore_template.h"

static PyObject* __Pyx_MakeTemplateLibTemplate(PyObject *strings, PyObject *interpolations) {
    (void)__Pyx_GetObjectFromTemplateLib;
    return _PyTemplate_Build(strings, interpolations);
}
#else

static PyObject* __Pyx_MakeTemplateLibTemplate(PyObject *strings, PyObject *interpolations) {
    PyObject *tp = __Pyx_GetObjectFromTemplateLib(1);
    PyObject *result = NULL, *zipped_tuple = NULL;
    Py_ssize_t zipped_index = 0;
    if (unlikely(!tp)) return NULL;
    
#if __PYX_LIMITED_VERSION_HEX < 0x030E0000
    if (__Pyx_get_runtime_version() < 0x030E0000) {
        // There's a high chance (but not certain) that we're using our internal
        // fallback version of template. In this case we can try to use a better
        // constructor.
        PyObject *args[] = { NULL, NULL };
        PyObject *kwargs_builder = __Pyx_MakeVectorcallBuilderKwds(2);
        if (unlikely(!kwargs_builder)) goto failed_shortcut; 
        if (unlikely(__Pyx_VectorcallBuilder_AddArg(PYIDENT("strings"), strings, kwargs_builder, args, 0)<0))
            goto failed_shortcut;
        if (unlikely(__Pyx_VectorcallBuilder_AddArg(PYIDENT("interpolations"), interpolations, kwargs_builder, args, 1)<0))
            goto failed_shortcut;
        result = __Pyx_Object_Vectorcall_CallFromBuilder(tp, args, 0, kwargs_builder);
        if (result) {
            Py_DECREF(kwargs_builder);
            Py_DECREF(tp);
            return result;
        }

      failed_shortcut:
        Py_CLEAR(kwargs_builder);
        PyErr_Clear();
    }
#endif

    // Slightly frustratingly, we have to go to the trouble of zipping together
    // the strings and the interpolations. And then internally the Template does
    // the exact opposite. 
    Py_ssize_t strings_len, interpolations_len;
    strings_len = __Pyx_PyTuple_GET_SIZE(strings);
#if !CYTHON_ASSUME_SAFE_SIZE
    if (unlikely(strings_len < 0)) goto end;
#endif
    interpolations_len = __Pyx_PyTuple_GET_SIZE(interpolations);
#if !CYTHON_ASSUME_SAFE_SIZE
    if (unlikely(interpolations_len < 0)) goto end;
#endif
    zipped_tuple = PyTuple_New(strings_len + interpolations_len);
    if (!zipped_tuple) goto end;
    for (Py_ssize_t i=0; (i<interpolations_len || i<strings_len); ++i) {
        if (i < strings_len) {
            PyObject *s = __Pyx_PyTuple_GET_ITEM(strings, i);
#if !CYTHON_ASSUME_SAFE_MACROS
            if (unlikely(!s)) goto end;
#endif
            Py_INCREF(s);
            if (unlikely(__Pyx_PyTuple_SET_ITEM(zipped_tuple, zipped_index, s) < 0)) goto end;
            ++zipped_index;
        }
        if (i < interpolations_len) {
            PyObject *interpolation = __Pyx_PyTuple_GET_ITEM(interpolations, i);
#if !CYTHON_ASSUME_SAFE_MACROS
            if (unlikely(!interpolation)) goto end;
#endif
            Py_INCREF(interpolation);
            if (unlikely(__Pyx_PyTuple_SET_ITEM(zipped_tuple, zipped_index, interpolation) < 0)) goto end;
            ++zipped_index;
        }
    }

    result = PyObject_CallObject(tp, zipped_tuple);

  end:
    Py_XDECREF(zipped_tuple);
    Py_DECREF(tp);
    return result;
}
#endif
