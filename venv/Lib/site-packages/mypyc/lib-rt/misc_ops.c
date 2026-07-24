#include "pythoncapi_compat.h"

// Misc primitive operations + C helpers
//
// These are registered in mypyc.primitives.misc_ops.

#include <Python.h>
#include <patchlevel.h>
#include "CPy.h"

PyObject *CPy_GetCoro(PyObject *obj)
{
    // If the type has an __await__ method, call it,
    // otherwise, fallback to calling __iter__.
    PyAsyncMethods* async_struct = Py_TYPE(obj)->tp_as_async;
    if (async_struct != NULL && async_struct->am_await != NULL) {
        return (async_struct->am_await)(obj);
    } else {
        // TODO: We should check that the type is a generator decorated with
        // asyncio.coroutine
        return PyObject_GetIter(obj);
    }
}

PyObject *CPyIter_Send(PyObject *iter, PyObject *val)
{
    // Do a send, or a next if second arg is None.
    // (This behavior is to match the PEP 380 spec for yield from.)
    if (Py_IsNone(val)) {
        return CPyIter_Next(iter);
    } else {
        return PyObject_CallMethodOneArg(iter, mypyc_interned_str.send, val);
    }
}

// A somewhat hairy implementation of specifically most of the error handling
// in `yield from` error handling. The point here is to reduce code size.
//
// This implements most of the bodies of the `except` blocks in the
// pseudocode in PEP 380.
//
// Returns true (1) if a StopIteration was received and we should return.
// Returns false (0) if a value should be yielded.
// In both cases the value is stored in outp.
// Signals an error (2) if the an exception should be propagated.
int CPy_YieldFromErrorHandle(PyObject *iter, PyObject **outp)
{
    PyObject *exc_type = (PyObject *)Py_TYPE(CPy_ExcState()->exc_value);
    PyObject *type, *value, *traceback;
    PyObject *_m;
    PyObject *res;
    *outp = NULL;

    if (PyErr_GivenExceptionMatches(exc_type, PyExc_GeneratorExit)) {
        _m = PyObject_GetAttr(iter, mypyc_interned_str.close_);
        if (_m) {
            res = PyObject_CallNoArgs(_m);
            Py_DECREF(_m);
            if (!res)
                return 2;
            Py_DECREF(res);
        } else if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
            PyErr_Clear();
        } else {
            return 2;
        }
    } else {
        _m = PyObject_GetAttr(iter, mypyc_interned_str.throw_);
        if (_m) {
            _CPy_GetExcInfo(&type, &value, &traceback);
            res = PyObject_CallFunctionObjArgs(_m, type, value, traceback, NULL);
            Py_DECREF(type);
            Py_DECREF(value);
            Py_DECREF(traceback);
            Py_DECREF(_m);
            if (res) {
                *outp = res;
                return 0;
            } else {
                res = CPy_FetchStopIterationValue();
                if (res) {
                    *outp = res;
                    return 1;
                } else {
                    return 2;
                }
            }
        } else if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
            PyErr_Clear();
        } else {
            return 2;
        }
    }

    CPy_Reraise();
    return 2;
}

PyObject *CPy_FetchStopIterationValue(void)
{
    PyObject *val = NULL;
    _PyGen_FetchStopIterationValue(&val);
    return val;
}

static bool _CPy_IsSafeMetaClass(PyTypeObject *metaclass) {
    // mypyc classes can't work with metaclasses in
    // general. Through some various nasty hacks we *do*
    // manage to work with TypingMeta and its friends.
    if (metaclass == &PyType_Type)
        return true;
    PyObject *module = PyObject_GetAttr((PyObject *)metaclass, mypyc_interned_str.__module__);
    if (!module) {
        PyErr_Clear();
        return false;
    }

    bool matches = false;
    if (PyUnicode_CompareWithASCIIString(module, "typing") == 0 &&
            (strcmp(metaclass->tp_name, "TypingMeta") == 0
             || strcmp(metaclass->tp_name, "GenericMeta") == 0
             || strcmp(metaclass->tp_name, "_ProtocolMeta") == 0)) {
        matches = true;
    } else if (PyUnicode_CompareWithASCIIString(module, "typing_extensions") == 0 &&
               strcmp(metaclass->tp_name, "_ProtocolMeta") == 0) {
        matches = true;
    } else if (PyUnicode_CompareWithASCIIString(module, "abc") == 0 &&
               strcmp(metaclass->tp_name, "ABCMeta") == 0) {
        matches = true;
    }
    Py_DECREF(module);
    return matches;
}

#if CPY_3_13_FEATURES

// Adapted from CPython 3.13.0b3
/* Determine the most derived metatype. */
PyObject *CPy_CalculateMetaclass(PyObject *metatype, PyObject *bases)
{
    Py_ssize_t i, nbases;
    PyTypeObject *winner;
    PyObject *tmp;
    PyTypeObject *tmptype;

    /* Determine the proper metatype to deal with this,
       and check for metatype conflicts while we're at it.
       Note that if some other metatype wins to contract,
       it's possible that its instances are not types. */

    nbases = PyTuple_GET_SIZE(bases);
    winner = (PyTypeObject *)metatype;
    for (i = 0; i < nbases; i++) {
        tmp = PyTuple_GET_ITEM(bases, i);
        tmptype = Py_TYPE(tmp);
        if (PyType_IsSubtype(winner, tmptype))
            continue;
        if (PyType_IsSubtype(tmptype, winner)) {
            winner = tmptype;
            continue;
        }
        /* else: */
        PyErr_SetString(PyExc_TypeError,
                        "metaclass conflict: "
                        "the metaclass of a derived class "
                        "must be a (non-strict) subclass "
                        "of the metaclasses of all its bases");
        return NULL;
    }
    return (PyObject *)winner;
}

#else

PyObject *CPy_CalculateMetaclass(PyObject *metatype, PyObject *bases) {
    return (PyObject *)_PyType_CalculateMetaclass((PyTypeObject *)metatype, bases);
}

#endif

// Create a heap type based on a template non-heap type.
// This is super hacky and maybe we should suck it up and use PyType_FromSpec instead.
// We allow bases to be NULL to represent just inheriting from object.
// We don't support NULL bases and a non-type metaclass.
PyObject *CPyType_FromTemplate(PyObject *template,
                               PyObject *orig_bases,
                               PyObject *modname) {
    PyTypeObject *template_ = (PyTypeObject *)template;
    PyHeapTypeObject *t = NULL;
    PyTypeObject *dummy_class = NULL;
    PyObject *name = NULL;
    PyObject *bases = NULL;
    PyObject *slots;

    // If the type of the class (the metaclass) is NULL, we default it
    // to being type.  (This allows us to avoid needing to initialize
    // it explicitly on windows.)
    if (!Py_TYPE(template_)) {
        Py_SET_TYPE(template_, &PyType_Type);
    }
    PyTypeObject *metaclass = Py_TYPE(template_);

    if (orig_bases) {
        bases = update_bases(orig_bases);
        // update_bases doesn't increment the refcount if nothing changes,
        // so we do it to make sure we have distinct "references" to both
        if (bases == orig_bases)
            Py_INCREF(bases);

        // Find the appropriate metaclass from our base classes. We
        // care about this because Generic uses a metaclass prior to
        // Python 3.7.
        metaclass = (PyTypeObject *)CPy_CalculateMetaclass((PyObject *)metaclass, bases);
        if (!metaclass)
            goto error;

        if (!_CPy_IsSafeMetaClass(metaclass)) {
            PyErr_SetString(PyExc_TypeError, "mypyc classes can't have a metaclass");
            goto error;
        }
    }

    name = PyUnicode_FromString(template_->tp_name);
    if (!name)
        goto error;

    if (template_->tp_doc) {
        // cpython expects tp_doc to be heap-allocated so convert it here to
        // avoid segfaults on deallocation.
        Py_ssize_t size = strlen(template_->tp_doc) + 1;
        char *doc = (char *)PyMem_Malloc(size);
        if (!doc)
            goto error;
        memcpy(doc, template_->tp_doc, size);
        template_->tp_doc = doc;
    }

    // Allocate the type and then copy the main stuff in.
    t = (PyHeapTypeObject*)PyType_GenericAlloc(&PyType_Type, 0);
    if (!t)
        goto error;
    memcpy((char *)t + sizeof(PyVarObject),
           (char *)template_ + sizeof(PyVarObject),
           sizeof(PyTypeObject) - sizeof(PyVarObject));

    if (bases != orig_bases) {
        if (PyObject_SetAttr((PyObject *)t, mypyc_interned_str.__orig_bases__, orig_bases) < 0)
            goto error;
    }

    // Having tp_base set is I think required for stuff to get
    // inherited in PyType_Ready, which we needed for subclassing
    // BaseException. XXX: Taking the first element is wrong I think though.
    if (bases) {
        t->ht_type.tp_base = (PyTypeObject *)PyTuple_GET_ITEM(bases, 0);
        Py_INCREF((PyObject *)t->ht_type.tp_base);
    }

    t->ht_name = name;
    Py_INCREF(name);
    t->ht_qualname = name;
    t->ht_type.tp_bases = bases;
    // references stolen so NULL these out
    bases = name = NULL;

    if (PyType_Ready((PyTypeObject *)t) < 0)
        goto error;

    assert(t->ht_type.tp_base != NULL);

    // XXX: This is a terrible hack to work around a cpython check on
    // the mro. It was needed for mypy.stats. I need to investigate
    // what is actually going on here.
    Py_INCREF(metaclass);
    Py_SET_TYPE(t, metaclass);

    if (dummy_class) {
        if (PyDict_Merge(t->ht_type.tp_dict, dummy_class->tp_dict, 0) != 0)
            goto error;
        // This is the *really* tasteless bit. GenericMeta's __new__
        // in certain versions of typing sets _gorg to point back to
        // the class. We need to override it to keep it from pointing
        // to the proxy.
        if (PyDict_SetItemString(t->ht_type.tp_dict, "_gorg", (PyObject *)t) < 0)
            goto error;
    }

    // Reject anything that would give us a nontrivial __slots__,
    // because the layout will conflict
    slots = PyObject_GetAttr((PyObject *)t, mypyc_interned_str.__slots__);
    if (slots) {
        // don't fail on an empty __slots__
        int is_true = PyObject_IsTrue(slots);
        Py_DECREF(slots);
        if (is_true > 0)
            PyErr_SetString(PyExc_TypeError, "mypyc classes can't have __slots__");
        if (is_true != 0)
            goto error;
    } else {
        PyErr_Clear();
    }

    if (PyObject_SetAttr((PyObject *)t, mypyc_interned_str.__module__, modname) < 0)
        goto error;

    Py_XDECREF(dummy_class);

    // Unlike the tp_doc slots of most other object, a heap type's tp_doc
    // must be heap allocated.
    if (template_->tp_doc) {
        // Silently truncate the docstring if it contains a null byte
        Py_ssize_t size = strlen(template_->tp_doc) + 1;
        char *tp_doc = (char *)PyMem_Malloc(size);
        if (tp_doc == NULL) {
            PyErr_NoMemory();
            goto error;
        }

        memcpy(tp_doc, template_->tp_doc, size);
        t->ht_type.tp_doc = tp_doc;
    }

#if PY_MINOR_VERSION == 11
    // This is a hack. Python 3.11 doesn't include good public APIs to work with managed
    // dicts, which are the default for heap types. So we try to opt-out until Python 3.12.
    t->ht_type.tp_flags &= ~Py_TPFLAGS_MANAGED_DICT;
#endif
    return (PyObject *)t;

error:
    Py_XDECREF(t);
    Py_XDECREF(bases);
    Py_XDECREF(dummy_class);
    Py_XDECREF(name);
    return NULL;
}

// Call __init_subclass__ on the appropriate base class of type.
// This is separated from CPyType_FromTemplate so that class attributes
// can be set before __init_subclass__ is called.
bool CPy_InitSubclass(PyObject *type) {
    if (init_subclass((PyTypeObject *)type, NULL)) {
        return false;
    }
    return true;
}

static int _CPy_UpdateObjFromDict(PyObject *obj, PyObject *dict)
{
    Py_ssize_t pos = 0;
    PyObject *key, *value;
    while (PyDict_Next(dict, &pos, &key, &value)) {
        if (PyObject_SetAttr(obj, key, value) != 0) {
            return -1;
        }
    }
    return 0;
}

/* Support for our partial built-in support for dataclasses.
 *
 * Take a class we want to make a dataclass, remove any descriptors
 * for annotated attributes, swap in the actual values of the class
 * variables invoke dataclass, and then restore all of the
 * descriptors.
 *
 * The purpose of all this is that dataclasses uses the values of
 * class variables to drive which attributes are required and what the
 * default values/factories are for optional attributes. This means
 * that the class dict needs to contain those values instead of getset
 * descriptors for the attributes when we invoke dataclass.
 *
 * We need to remove descriptors for attributes even when there is no
 * default value for them, or else dataclass will think the descriptor
 * is the default value. We remove only the attributes, since we don't
 * want dataclasses to try generating functions when they are already
 * implemented.
 *
 * Args:
 *   dataclass_dec: The decorator to apply
 *   tp: The class we are making a dataclass
 *   dict: The dictionary containing values that dataclasses needs
 *   annotations: The type annotation dictionary
 *   dataclass_type: A str object with the return value of util.py:dataclass_type()
 */
int
CPyDataclass_SleightOfHand(PyObject *dataclass_dec, PyObject *tp,
                           PyObject *dict, PyObject *annotations,
                           PyObject *dataclass_type) {
    PyTypeObject *ttp = (PyTypeObject *)tp;
    Py_ssize_t pos;
    PyObject *res = NULL;

    /* Make a copy of the original class __dict__ */
    PyObject *orig_dict = PyDict_Copy(ttp->tp_dict);
    if (!orig_dict) {
        goto fail;
    }

    /* Delete anything that had an annotation */
    pos = 0;
    PyObject *key;
    while (PyDict_Next(annotations, &pos, &key, NULL)) {
        // Check and delete key. Key may be absent from tp for InitVar variables.
        if (PyObject_HasAttr(tp, key) == 1 && PyObject_DelAttr(tp, key) != 0) {
            goto fail;
        }
    }

    /* Copy in all the attributes that we want dataclass to see */
    if (_CPy_UpdateObjFromDict(tp, dict) != 0) {
        goto fail;
    }

    /* Run the @dataclass descriptor */
    res = PyObject_CallOneArg(dataclass_dec, tp);
    if (!res) {
        goto fail;
    }
    const char *dataclass_type_ptr = PyUnicode_AsUTF8(dataclass_type);
    if (dataclass_type_ptr == NULL) {
        goto fail;
    }
    if (strcmp(dataclass_type_ptr, "attr") == 0 ||
        strcmp(dataclass_type_ptr, "attr-auto") == 0) {
        // These attributes are added or modified by @attr.s(slots=True).
        const char * const keys[] = {"__attrs_attrs__", "__attrs_own_setattr__", "__init__", ""};
        for (const char * const *key_iter = keys; **key_iter != '\0'; key_iter++) {
            PyObject *value = NULL;
            int rv = PyObject_GetOptionalAttrString(res, *key_iter, &value);
            if (rv == 1) {
                PyObject_SetAttrString(tp, *key_iter, value);
                Py_DECREF(value);
            } else if (rv == -1) {
              goto fail;
            }
        }
    }

    /* Copy back the original contents of the dict */
    if (_CPy_UpdateObjFromDict(tp, orig_dict) != 0) {
        goto fail;
    }

    Py_DECREF(res);
    Py_DECREF(orig_dict);
    return 1;

fail:
    Py_XDECREF(res);
    Py_XDECREF(orig_dict);
    return 0;
}

// Support for pickling; reusable getstate and setstate functions
PyObject *
CPyPickle_SetState(PyObject *obj, PyObject *state)
{
    if (_CPy_UpdateObjFromDict(obj, state) != 0) {
        return NULL;
    }
    Py_RETURN_NONE;
}

PyObject *
CPyPickle_GetState(PyObject *obj)
{
    PyObject *attrs = NULL, *state = NULL;

    attrs = PyObject_GetAttr((PyObject *)Py_TYPE(obj), mypyc_interned_str.__mypyc_attrs__);
    if (!attrs) {
        goto fail;
    }
    if (!PyTuple_Check(attrs)) {
        PyErr_SetString(PyExc_TypeError, "__mypyc_attrs__ is not a tuple");
        goto fail;
    }
    state = PyDict_New();
    if (!state) {
        goto fail;
    }

    // Collect all the values of attributes in __mypyc_attrs__
    // Attributes that are missing we just ignore
    int i;
    for (i = 0; i < PyTuple_GET_SIZE(attrs); i++) {
        PyObject *key = PyTuple_GET_ITEM(attrs, i);
        PyObject *value = PyObject_GetAttr(obj, key);
        if (!value) {
            if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
                PyErr_Clear();
                continue;
            }
            goto fail;
        }
        int result = PyDict_SetItem(state, key, value);
        Py_DECREF(value);
        if (result != 0) {
            goto fail;
        }
    }

    Py_DECREF(attrs);

    return state;
fail:
    Py_XDECREF(attrs);
    Py_XDECREF(state);
    return NULL;
}

CPyTagged CPyTagged_Id(PyObject *o) {
    return CPyTagged_FromVoidPtr(o);
}

#define MAX_INT_CHARS 22
#define _PyUnicode_LENGTH(op)                           \
    (((PyASCIIObject *)(op))->length)

// using snprintf or PyUnicode_FromFormat was way slower than
// boxing the int and calling PyObject_Str on it, so we implement our own
static int fmt_ssize_t(char *out, Py_ssize_t n) {
	bool neg = n < 0;
	if (neg) n = -n;

	// buf gets filled backward and then we copy it forward
	char buf[MAX_INT_CHARS];
	int i = 0;
	do {
		buf[i] = (n % 10) + '0';
		n /= 10;
		i++;
	} while (n);


	int len = i;
	int j = 0;
	if (neg) {
		out[j++] = '-';
		len++;
	}

	for (; j < len; j++, i--) {
		out[j] = buf[i-1];
	}
	out[j] = '\0';

	return len;
}

static PyObject *CPyTagged_ShortToStr(Py_ssize_t n) {
    PyObject *obj = PyUnicode_New(MAX_INT_CHARS, 127);
    if (!obj) return NULL;
    int len = fmt_ssize_t((char *)PyUnicode_1BYTE_DATA(obj), n);
    _PyUnicode_LENGTH(obj) = len;
    return obj;
}

PyObject *CPyTagged_Str(CPyTagged n) {
    if (CPyTagged_CheckShort(n)) {
        return CPyTagged_ShortToStr(CPyTagged_ShortAsSsize_t(n));
    } else {
        return PyObject_Str(CPyTagged_AsObject(n));
    }
}

static PyObject *CPyTagged_ShortToAsciiBytes(Py_ssize_t n) {
    PyObject *obj = PyBytes_FromStringAndSize(NULL, MAX_INT_CHARS);
    if (!obj) return NULL;
    int len = fmt_ssize_t(PyBytes_AsString(obj), n);
    Py_SET_SIZE(obj, len);
    return obj;
}

PyObject *CPyTagged_AsciiBytes(CPyTagged n) {
    if (CPyTagged_CheckShort(n)) {
        return CPyTagged_ShortToAsciiBytes(CPyTagged_ShortAsSsize_t(n));
    }
    PyObject *str = PyObject_Str(CPyTagged_AsObject(n));
    PyObject *bytes = PyUnicode_AsASCIIString(str);
    CPy_DECREF(str);
    return bytes;
}

void CPyDebug_Print(const char *msg) {
    printf("%s\n", msg);
    fflush(stdout);
}

void CPyDebug_PrintObject(PyObject *obj) {
    // Printing can cause errors. We don't want this to affect any existing
    // state so we'll save any existing error and restore it at the end.
    PyObject *exc_type, *exc_value, *exc_traceback;
    PyErr_Fetch(&exc_type, &exc_value, &exc_traceback);

    if (PyObject_Print(obj, stderr, 0) == -1) {
        PyErr_Print();
    } else {
        fprintf(stderr, "\n");
    }
    fflush(stderr);

    PyErr_Restore(exc_type, exc_value, exc_traceback);
}

int CPySequence_CheckUnpackCount(PyObject *sequence, Py_ssize_t expected) {
    Py_ssize_t actual = Py_SIZE(sequence);
    if (unlikely(actual != expected)) {
        if (actual < expected) {
            PyErr_Format(PyExc_ValueError, "not enough values to unpack (expected %zd, got %zd)",
                         expected, actual);
        } else {
            PyErr_Format(PyExc_ValueError, "too many values to unpack (expected %zd)", expected);
        }
        return -1;
    }
    return 0;
}

// Parse an integer (size_t) encoded as a variable-length binary sequence.
static const char *parse_int(const char *s, size_t *len) {
    Py_ssize_t n = 0;
    while ((unsigned char)*s >= 0x80) {
        n = (n << 7) + (*s & 0x7f);
        s++;
    }
    n = (n << 7) | *s++;
    *len = n;
    return s;
}

// Initialize static constant array of literal values
int CPyStatics_Initialize(PyObject **statics,
                          const char * const *strings,
                          const char * const *bytestrings,
                          const char * const *ints,
                          const double *floats,
                          const double *complex_numbers,
                          const int *tuples,
                          const int *frozensets) {
    PyObject **result = statics;
    // Start with some hard-coded values
    *result++ = Py_None;
    Py_INCREF(Py_None);
    *result++ = Py_False;
    Py_INCREF(Py_False);
    *result++ = Py_True;
    Py_INCREF(Py_True);
    if (strings) {
        for (; **strings != '\0'; strings++) {
            size_t num;
            const char *data = *strings;
            data = parse_int(data, &num);
            while (num-- > 0) {
                size_t len;
                data = parse_int(data, &len);
                PyObject *obj = PyUnicode_DecodeUTF8(data, len, "surrogatepass");
                if (obj == NULL) {
                    return -1;
                }
                PyUnicode_InternInPlace(&obj);
                *result++ = obj;
                data += len;
            }
        }
    }
    if (bytestrings) {
        for (; **bytestrings != '\0'; bytestrings++) {
            size_t num;
            const char *data = *bytestrings;
            data = parse_int(data, &num);
            while (num-- > 0) {
                size_t len;
                data = parse_int(data, &len);
                PyObject *obj = PyBytes_FromStringAndSize(data, len);
                if (obj == NULL) {
                    return -1;
                }
                *result++ = obj;
                data += len;
            }
        }
    }
    if (ints) {
        for (; **ints != '\0'; ints++) {
            size_t num;
            const char *data = *ints;
            data = parse_int(data, &num);
            while (num-- > 0) {
                char *end;
                PyObject *obj = PyLong_FromString(data, &end, 10);
                if (obj == NULL) {
                    return -1;
                }
                data = end;
                data++;
                *result++ = obj;
            }
        }
    }
    if (floats) {
        size_t num = (size_t)*floats++;
        while (num-- > 0) {
            PyObject *obj = PyFloat_FromDouble(*floats++);
            if (obj == NULL) {
                return -1;
            }
            *result++ = obj;
        }
    }
    if (complex_numbers) {
        size_t num = (size_t)*complex_numbers++;
        while (num-- > 0) {
            double real = *complex_numbers++;
            double imag = *complex_numbers++;
            PyObject *obj = PyComplex_FromDoubles(real, imag);
            if (obj == NULL) {
                return -1;
            }
            *result++ = obj;
        }
    }
    if (tuples) {
        int num = *tuples++;
        while (num-- > 0) {
            int num_items = *tuples++;
            PyObject *obj = PyTuple_New(num_items);
            if (obj == NULL) {
                return -1;
            }
            int i;
            for (i = 0; i < num_items; i++) {
                PyObject *item = statics[*tuples++];
                Py_INCREF(item);
                PyTuple_SET_ITEM(obj, i, item);
            }
            *result++ = obj;
        }
    }
    if (frozensets) {
        int num = *frozensets++;
        while (num-- > 0) {
            int num_items = *frozensets++;
            PyObject *obj = PyFrozenSet_New(NULL);
            if (obj == NULL) {
                return -1;
            }
            for (int i = 0; i < num_items; i++) {
                PyObject *item = statics[*frozensets++];
                Py_INCREF(item);
                if (PySet_Add(obj, item) == -1) {
                    return -1;
                }
            }
            *result++ = obj;
        }
    }
    return 0;
}

// Call super(type(self), self)
PyObject *
CPy_Super(PyObject *builtins, PyObject *self) {
    PyObject *super_type = PyObject_GetAttr(builtins, mypyc_interned_str.super);
    if (!super_type)
        return NULL;
    PyObject *result = PyObject_CallFunctionObjArgs(
        super_type, (PyObject*)Py_TYPE(self), self, NULL);
    Py_DECREF(super_type);
    return result;
}

static bool import_single(PyObject *mod_id, PyObject **mod_static,
                          PyObject *globals_id, PyObject *globals_name, PyObject *globals) {
    if (Py_IsNone(*mod_static)) {
        CPyModule *mod = PyImport_Import(mod_id);
        if (mod == NULL) {
            return false;
        }
        *mod_static = mod;
    }

    PyObject *mod_dict = PyImport_GetModuleDict();
    CPyModule *globals_mod = CPyDict_GetItem(mod_dict, globals_id);
    if (globals_mod == NULL) {
        return false;
    }
    int ret = CPyDict_SetItem(globals, globals_name, globals_mod);
    Py_DECREF(globals_mod);
    if (ret < 0) {
        return false;
    }

    return true;
}

// Table-driven import helper. See transform_import() in irbuild for the details.
bool CPyImport_ImportMany(PyObject *modules, CPyModule **statics[], PyObject *globals,
                          PyObject *tb_path, PyObject *tb_function, Py_ssize_t *tb_lines) {
    for (Py_ssize_t i = 0; i < PyTuple_GET_SIZE(modules); i++) {
        PyObject *module = PyTuple_GET_ITEM(modules, i);
        PyObject *mod_id = PyTuple_GET_ITEM(module, 0);
        PyObject *globals_id = PyTuple_GET_ITEM(module, 1);
        PyObject *globals_name = PyTuple_GET_ITEM(module, 2);

        if (!import_single(mod_id, statics[i], globals_id, globals_name, globals)) {
            assert(PyErr_Occurred() && "error indicator should be set on bad import!");
            PyObject *typ, *val, *tb;
            PyErr_Fetch(&typ, &val, &tb);
            const char *path = PyUnicode_AsUTF8(tb_path);
            if (path == NULL) {
                path = "<unable to display>";
            }
            const char *function = PyUnicode_AsUTF8(tb_function);
            if (function == NULL) {
                function = "<unable to display>";
            }
            PyErr_Restore(typ, val, tb);
            CPy_AddTraceback(path, function, tb_lines[i], globals);
            return false;
        }
    }
    return true;
}

// This helper function is a simplification of cpython/ceval.c/import_from()
static PyObject *CPyImport_ImportFrom(PyObject *module, PyObject *package_name,
                                      PyObject *import_name, PyObject *as_name) {
    // check if the imported module has an attribute by that name
    PyObject *x = PyObject_GetAttr(module, import_name);
    if (x == NULL) {
        // Attribute lookup failed. The name may still be a submodule that's
        // been imported already; look it up directly in sys.modules.
        PyObject *fullmodname = PyUnicode_FromFormat("%U.%U", package_name, import_name);
        if (fullmodname == NULL) {
            goto fail;
        }
        PyErr_Clear();
        x = PyImport_GetModule(fullmodname);
        Py_DECREF(fullmodname);
        if (x == NULL) {
            goto fail;
        }
    }
    return x;

fail:
    PyErr_Clear();
    PyObject *package_path = PyModule_GetFilenameObject(module);
    PyObject *path_for_msg = package_path != NULL ? package_path : Py_None;
    PyObject *errmsg = PyUnicode_FromFormat("cannot import name %R from %R (%S)",
                                            import_name, package_name, path_for_msg);
    // NULL checks for errmsg and package_name done by PyErr_SetImportError.
    PyErr_SetImportError(errmsg, package_name, package_path);
    Py_XDECREF(package_path);
    Py_XDECREF(errmsg);
    return NULL;
}

PyObject *CPyImport_ImportFromMany(PyObject *mod_id, PyObject *names, PyObject *as_names,
                                   PyObject *globals) {
    PyObject *mod = PyImport_ImportModuleLevelObject(mod_id, globals, 0, names, 0);
    if (mod == NULL) {
        return NULL;
    }

    for (Py_ssize_t i = 0; i < PyTuple_GET_SIZE(names); i++) {
        PyObject *name = PyTuple_GET_ITEM(names, i);
        PyObject *as_name = PyTuple_GET_ITEM(as_names, i);
        PyObject *obj = CPyImport_ImportFrom(mod, mod_id, name, as_name);
        if (obj == NULL) {
            Py_DECREF(mod);
            return NULL;
        }
        int ret = CPyDict_SetItem(globals, as_name, obj);
        Py_DECREF(obj);
        if (ret < 0) {
            Py_DECREF(mod);
            return NULL;
        }
    }
    return mod;
}

// Import attributes from an already-imported native module and store them
// in the globals dict.  Returns the module on success, NULL on error.
PyObject *CPyImport_GetNativeAttrs(PyObject *mod_id, PyObject *names,
                                   PyObject *as_names, PyObject *globals) {
    PyObject *mod = PyImport_GetModule(mod_id);
    if (mod == NULL) {
        if (!PyErr_Occurred()) {
            PyErr_Format(PyExc_ImportError, "module '%U' is not in sys.modules", mod_id);
        }
        return NULL;
    }
    for (Py_ssize_t i = 0; i < PyTuple_GET_SIZE(names); i++) {
        PyObject *name = PyTuple_GET_ITEM(names, i);
        PyObject *as_name = PyTuple_GET_ITEM(as_names, i);
        PyObject *obj = PyObject_GetAttr(mod, name);
        if (obj == NULL) {
            Py_DECREF(mod);
            return NULL;
        }
        int ret = CPyDict_SetItem(globals, as_name, obj);
        Py_DECREF(obj);
        if (ret < 0) {
            Py_DECREF(mod);
            return NULL;
        }
    }
    return mod;
}

// From CPython
static PyObject *
CPy_BinopTypeError(PyObject *left, PyObject *right, const char *op) {
    PyErr_Format(PyExc_TypeError,
                 "unsupported operand type(s) for %.100s: "
                 "'%.100s' and '%.100s'",
                 op,
                 Py_TYPE(left)->tp_name,
                 Py_TYPE(right)->tp_name);
    return NULL;
}

PyObject *
CPy_CallReverseOpMethod(PyObject *left,
                        PyObject *right,
                        const char *op,
                        PyObject *method) {
    // Look up reverse method
    PyObject *m = PyObject_GetAttr(right, method);
    if (m == NULL) {
        // If reverse method not defined, generate TypeError instead AttributeError
        if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
            CPy_BinopTypeError(left, right, op);
        }
        return NULL;
    }
    // Call reverse method
    PyObject *result = PyObject_CallOneArg(m, left);
    Py_DECREF(m);
    return result;
}

PyObject *CPySingledispatch_RegisterFunction(PyObject *singledispatch_func,
                                             PyObject *cls,
                                             PyObject *func) {
    PyObject *registry = PyObject_GetAttr(singledispatch_func, mypyc_interned_str.registry);
    PyObject *register_func = NULL;
    PyObject *typing = NULL;
    PyObject *get_type_hints = NULL;
    PyObject *type_hints = NULL;

    if (registry == NULL) goto fail;
    if (func == NULL) {
        // one argument case
        if (PyType_Check(cls)) {
            // passed a class
            // bind cls to the first argument so that register gets called again with both the
            // class and the function
            register_func = PyObject_GetAttr(singledispatch_func, mypyc_interned_str.register_);
            if (register_func == NULL) goto fail;
            return PyMethod_New(register_func, cls);
        }
        // passed a function
        PyObject *annotations = PyFunction_GetAnnotations(cls);
        const char *invalid_first_arg_msg =
            "Invalid first argument to `register()`: %R. "
            "Use either `@register(some_class)` or plain `@register` "
            "on an annotated function.";

        if (annotations == NULL) {
            PyErr_Format(PyExc_TypeError, invalid_first_arg_msg, cls);
            goto fail;
        }

        Py_INCREF(annotations);

        func = cls;
        typing = PyImport_ImportModule("typing");
        if (typing == NULL) goto fail;
        get_type_hints = PyObject_GetAttr(typing, mypyc_interned_str.get_type_hints);

        type_hints = PyObject_CallOneArg(get_type_hints, func);
        PyObject *argname;
        Py_ssize_t pos = 0;
        if (!PyDict_Next(type_hints, &pos, &argname, &cls)) {
            // the functools implementation raises the same type error if annotations is an empty dict
            PyErr_Format(PyExc_TypeError, invalid_first_arg_msg, cls);
            goto fail;
        }
        if (!PyType_Check(cls)) {
            const char *invalid_annotation_msg = "Invalid annotation for %R. %R is not a class.";
            PyErr_Format(PyExc_TypeError, invalid_annotation_msg, argname, cls);
            goto fail;
        }
    }
    if (PyDict_SetItem(registry, cls, func) == -1) {
        goto fail;
    }

    // clear the cache so we consider the newly added function when dispatching
    PyObject *dispatch_cache = PyObject_GetAttr(singledispatch_func, mypyc_interned_str.dispatch_cache);
    if (dispatch_cache == NULL) goto fail;
    PyDict_Clear(dispatch_cache);

    Py_INCREF(func);
    return func;

fail:
    Py_XDECREF(registry);
    Py_XDECREF(register_func);
    Py_XDECREF(typing);
    Py_XDECREF(get_type_hints);
    Py_XDECREF(type_hints);
    return NULL;

}

// Adapted from ceval.c GET_AITER
PyObject *CPy_GetAIter(PyObject *obj)
{
    unaryfunc getter = NULL;
    PyTypeObject *type = Py_TYPE(obj);

    if (type->tp_as_async != NULL) {
        getter = type->tp_as_async->am_aiter;
    }

    if (getter == NULL) {
        PyErr_Format(PyExc_TypeError,
                     "'async for' requires an object with "
                     "__aiter__ method, got %.100s",
                     type->tp_name);
        Py_DECREF(obj);
        return NULL;
    }

    PyObject *iter = (*getter)(obj);
    if (!iter) {
        return NULL;
    }

    if (Py_TYPE(iter)->tp_as_async == NULL ||
        Py_TYPE(iter)->tp_as_async->am_anext == NULL) {

        PyErr_Format(PyExc_TypeError,
                     "'async for' received an object from __aiter__ "
                     "that does not implement __anext__: %.100s",
                     Py_TYPE(iter)->tp_name);
        Py_DECREF(iter);
        return NULL;
    }

    return iter;
}

// Adapted from ceval.c GET_ANEXT
PyObject *CPy_GetANext(PyObject *aiter)
{
    unaryfunc getter = NULL;
    PyObject *next_iter = NULL;
    PyObject *awaitable = NULL;
    PyTypeObject *type = Py_TYPE(aiter);

    if (PyAsyncGen_CheckExact(aiter)) {
        awaitable = type->tp_as_async->am_anext(aiter);
        if (awaitable == NULL) {
            goto error;
        }
    } else {
        if (type->tp_as_async != NULL){
            getter = type->tp_as_async->am_anext;
        }

        if (getter != NULL) {
            next_iter = (*getter)(aiter);
            if (next_iter == NULL) {
                goto error;
            }
        }
        else {
            PyErr_Format(PyExc_TypeError,
                         "'async for' requires an iterator with "
                         "__anext__ method, got %.100s",
                         type->tp_name);
            goto error;
        }

        awaitable = CPyCoro_GetAwaitableIter(next_iter);
        if (awaitable == NULL) {
            _PyErr_FormatFromCause(
                PyExc_TypeError,
                "'async for' received an invalid object "
                "from __anext__: %.100s",
                Py_TYPE(next_iter)->tp_name);

            Py_DECREF(next_iter);
            goto error;
        } else {
            Py_DECREF(next_iter);
        }
    }

    return awaitable;
error:
    return NULL;
}

#if CPY_3_11_FEATURES

// Return obj.__name__ (specialized to type objects, which are the most common target).
PyObject *CPy_GetName(PyObject *obj) {
    if (PyType_Check(obj)) {
        return PyType_GetName((PyTypeObject *)obj);
    }
    return PyObject_GetAttr(obj, mypyc_interned_str.__name__);
}

#endif

#ifdef MYPYC_LOG_TRACE

// This is only compiled in if trace logging is enabled by user

static int TraceCounter = 0;
static const int TRACE_EVERY_NTH = 1009;  // Should be a prime number
#define TRACE_LOG_FILE_NAME "mypyc_trace.txt"
static FILE *TraceLogFile = NULL;

// Log a tracing event on every Nth call
void CPyTrace_LogEvent(const char *location, const char *line, const char *op, const char *details) {
    if (TraceLogFile == NULL) {
        if ((TraceLogFile = fopen(TRACE_LOG_FILE_NAME, "w")) == NULL) {
            fprintf(stderr, "error: Could not open trace file %s\n", TRACE_LOG_FILE_NAME);
            abort();
        }
    }
    if (TraceCounter == 0) {
        fprintf(TraceLogFile, "%s:%s:%s:%s\n", location, line, op, details);
    }
    TraceCounter++;
    if (TraceCounter == TRACE_EVERY_NTH) {
        TraceCounter = 0;
    }
}

#endif

#if CPY_3_12_FEATURES

// Copied from Python 3.12.3, since this struct is internal to CPython. It defines
// the structure of typing.TypeAliasType objects. We need it since compute_value is
// not part of the public API, and we need to set it to match Python runtime semantics.
//
// IMPORTANT: This needs to be kept in sync with CPython!
typedef struct {
    PyObject_HEAD
    PyObject *name;
#if CPY_3_15_FEATURES
    PyObject *qualname;
#endif
    PyObject *type_params;
    PyObject *compute_value;
    PyObject *value;
    PyObject *module;
} typealiasobject;

void CPy_SetTypeAliasTypeComputeFunction(PyObject *alias, PyObject *compute_value) {
    typealiasobject *obj = (typealiasobject *)alias;
    if (obj->value != NULL) {
        Py_DECREF(obj->value);
    }
    obj->value = NULL;
    Py_INCREF(compute_value);
    if (obj->compute_value != NULL) {
        Py_DECREF(obj->compute_value);
    }
    obj->compute_value = compute_value;
}

#endif

#ifdef _WIN32
#define SEP "\\"
#else
#define SEP "/"
#endif

// Cached class references for __spec__ / __loader__ construction.
static PyObject *CPyImport_ModuleSpecClass = NULL;
static PyObject *CPyImport_ExtFileLoaderClass = NULL;
static PyObject *CPyImport_SpecKwnames = NULL;  // ("origin", "is_package")

// Initialize cached references for ModuleSpec and ExtensionFileLoader.
// Returns 0 on success, -1 on error.
static int CPyImport_InitSpecClasses(void) {
    if (CPyImport_ModuleSpecClass != NULL) {
        return 0;
    }
    PyObject *machinery = PyImport_ImportModule("importlib.machinery");
    if (machinery == NULL) {
        return -1;
    }
    CPyImport_ModuleSpecClass = PyObject_GetAttrString(machinery, "ModuleSpec");
    CPyImport_ExtFileLoaderClass = PyObject_GetAttrString(machinery, "ExtensionFileLoader");
    Py_DECREF(machinery);
    if (CPyImport_ModuleSpecClass == NULL || CPyImport_ExtFileLoaderClass == NULL) {
        Py_CLEAR(CPyImport_ModuleSpecClass);
        Py_CLEAR(CPyImport_ExtFileLoaderClass);
        return -1;
    }
    PyObject *origin_str = PyUnicode_InternFromString("origin");
    PyObject *is_package_str = PyUnicode_InternFromString("is_package");
    if (origin_str == NULL || is_package_str == NULL) {
        CPyError_OutOfMemory();
    }
    CPyImport_SpecKwnames = PyTuple_Pack(2, origin_str, is_package_str);
    Py_DECREF(origin_str);
    Py_DECREF(is_package_str);
    if (CPyImport_SpecKwnames == NULL) {
        CPyError_OutOfMemory();
    }
    return 0;
}

// Set __package__ before executing the module body so it is available
// during module initialization. For a package, __package__ is the module
// name itself. For a non-package submodule "a.b.c", it is "a.b". For a
// top-level non-package module, it is "".
static int CPyImport_SetModulePackage(PyObject *modobj, PyObject *module_name,
                                      Py_ssize_t is_package) {
    PyObject *pkg = NULL;
    int rc = PyObject_GetOptionalAttrString(modobj, "__package__", &pkg);
    if (rc < 0) {
        return -1;
    }
    if (pkg != NULL && pkg != Py_None) {
        Py_DECREF(pkg);
        return 0;
    }
    Py_XDECREF(pkg);

    PyObject *package_name = NULL;
    if (is_package) {
        package_name = module_name;
        Py_INCREF(package_name);
    } else {
        Py_ssize_t name_len = PyUnicode_GetLength(module_name);
        if (name_len < 0) {
            return -1;
        }
        Py_ssize_t dot = PyUnicode_FindChar(module_name, '.', 0, name_len, -1);
        if (dot >= 0) {
            package_name = PyUnicode_Substring(module_name, 0, dot);
        } else {
            package_name = PyUnicode_FromString("");
        }
    }
    if (package_name == NULL) {
        return -1;
    }
    rc = PyObject_SetAttrString(modobj, "__package__", package_name);
    Py_DECREF(package_name);
    return rc;
}

// Derive and set __file__ on modobj from the shared library path, module name,
// and extension suffix. Returns 0 on success, -1 on error.
static int CPyImport_SetModuleFile(PyObject *modobj, PyObject *module_name,
                                    PyObject *shared_lib_file, PyObject *ext_suffix,
                                    Py_ssize_t is_package) {
    PyObject *file = NULL;
    int rc = PyObject_GetOptionalAttrString(modobj, "__file__", &file);
    if (rc < 0) {
        return -1;
    }
    if (file != NULL) {
        // __file__ already set, nothing to do.
        Py_DECREF(file);
        return 0;
    }
    // Derive __file__ from the shared lib's directory, the module
    // name, and the extension suffix. Two layouts:
    //
    //  Monolithic: one shared lib above the package tree holds many
    //    modules, so append the full dotted module path.
    //  separate=True: each module has its own "<segment>__mypyc.so"
    //    next to the module, so dirname(shared_lib) is already inside
    //    the parent package. Append only the last segment.
    //
    // Detect the separate=True case by matching the shared lib's
    // basename against "<last_segment>__mypyc<ext>".
    PyObject *derived_file = NULL;
    if (shared_lib_file != NULL && shared_lib_file != Py_None &&
            PyUnicode_Check(shared_lib_file)) {
        Py_ssize_t sf_len = PyUnicode_GetLength(shared_lib_file);
        // Find the last path separator, checking both '/' and '\\'
        // for cross-platform support.
        Py_ssize_t sep = PyUnicode_FindChar(shared_lib_file, '/', 0, sf_len, -1);
        Py_ssize_t bsep = PyUnicode_FindChar(shared_lib_file, '\\', 0, sf_len, -1);
        if (bsep > sep) {
            sep = bsep;
        }
        // Use the same separator character found in the path, or
        // the platform default if no separator was found.
        Py_UCS4 sep_char = sep >= 0
            ? PyUnicode_ReadChar(shared_lib_file, sep)
            : (Py_UCS4)SEP[0];
        PyObject *dot_str = PyUnicode_FromString(".");
        PyObject *sep_str = PyUnicode_FromOrdinal(sep_char);
        if (dot_str == NULL || sep_str == NULL) {
            CPyError_OutOfMemory();
        }
        PyObject *module_path = PyUnicode_Replace(module_name, dot_str, sep_str, -1);
        Py_DECREF(dot_str);
        Py_DECREF(sep_str);
        if (module_path == NULL) {
            return -1;
        }

        // Compute the module's last dotted segment for the separate=True check.
        Py_ssize_t name_len = PyUnicode_GetLength(module_name);
        Py_ssize_t last_dot = PyUnicode_FindChar(module_name, '.', 0, name_len, -1);
        PyObject *last_segment;
        if (last_dot >= 0) {
            last_segment = PyUnicode_Substring(module_name, last_dot + 1, name_len);
        } else {
            last_segment = module_name;
            Py_INCREF(last_segment);
        }
        if (last_segment == NULL) {
            Py_DECREF(module_path);
            return -1;
        }
        // Compare shared_lib_file basename against "<last_segment>__mypyc<ext>".
        PyObject *expected_basename = PyUnicode_FromFormat(
            "%U__mypyc%U", last_segment, ext_suffix);
        PyObject *actual_basename;
        if (sep >= 0) {
            actual_basename = PyUnicode_Substring(shared_lib_file, sep + 1, sf_len);
        } else {
            actual_basename = shared_lib_file;
            Py_INCREF(actual_basename);
        }
        int is_per_module_lib = 0;
        if (expected_basename != NULL && actual_basename != NULL) {
            is_per_module_lib =
                (PyUnicode_Compare(expected_basename, actual_basename) == 0);
        }
        Py_XDECREF(expected_basename);
        Py_XDECREF(actual_basename);

        // For packages, __file__ should point to __init__<ext>,
        // e.g. "a/b/__init__.cpython-312-x86_64-linux-gnu.so".
        PyObject *file_path = is_per_module_lib ? last_segment : module_path;
        if (sep >= 0) {
            PyObject *dir = PyUnicode_Substring(shared_lib_file, 0, sep);
            if (dir != NULL) {
                if (is_package) {
                    derived_file = PyUnicode_FromFormat(
                        "%U%c%U%c__init__%U", dir, (int)sep_char,
                        file_path, (int)sep_char, ext_suffix);
                } else {
                    derived_file = PyUnicode_FromFormat(
                        "%U%c%U%U", dir, (int)sep_char,
                        file_path, ext_suffix);
                }
                Py_DECREF(dir);
            }
        } else {
            if (is_package) {
                derived_file = PyUnicode_FromFormat(
                    "%U%c__init__%U", file_path, (int)SEP[0], ext_suffix);
            } else {
                derived_file = PyUnicode_FromFormat("%U%U", file_path, ext_suffix);
            }
        }
        Py_DECREF(last_segment);
        Py_DECREF(module_path);
    }
    if (derived_file == NULL && !PyErr_Occurred()) {
        derived_file = module_name;
        Py_INCREF(derived_file);
    }
    if (derived_file == NULL ||
            PyObject_SetAttrString(modobj, "__file__", derived_file) < 0) {
        Py_XDECREF(derived_file);
        return -1;
    }
    Py_DECREF(derived_file);
    return 0;
}

// Set __path__ to [dirname(__file__)] on a package module if not already set.
// Returns 0 on success, -1 on error.
static int CPyImport_SetModulePath(PyObject *modobj) {
    PyObject *existing_path = NULL;
    int rc = PyObject_GetOptionalAttrString(modobj, "__path__", &existing_path);
    if (rc < 0) {
        return -1;
    }
    if (existing_path != NULL) {
        Py_DECREF(existing_path);
        return 0;
    }
    PyObject *file = NULL;
    rc = PyObject_GetOptionalAttrString(modobj, "__file__", &file);
    if (rc <= 0) {
        return rc;
    }
    PyObject *os_path = PyImport_ImportModule("os.path");
    if (os_path == NULL) {
        Py_DECREF(file);
        return -1;
    }
    PyObject *dir = PyObject_CallMethod(os_path, "dirname", "O", file);
    Py_DECREF(os_path);
    Py_DECREF(file);
    if (dir == NULL) {
        return -1;
    }
    PyObject *path_list = PyList_New(1);
    if (path_list == NULL) {
        CPyError_OutOfMemory();
    }
    PyList_SET_ITEM(path_list, 0, dir);  // steals ref to dir
    int ret = PyObject_SetAttrString(modobj, "__path__", path_list);
    Py_DECREF(path_list);
    return ret;
}

// Set __spec__ and __loader__ on modobj if not already set.
// Returns 0 on success, -1 on error.
static int CPyImport_SetModuleSpec(PyObject *modobj, PyObject *module_name,
                                    Py_ssize_t is_package) {
    PyObject *spec = NULL;
    int rc = PyObject_GetOptionalAttrString(modobj, "__spec__", &spec);
    if (rc < 0) {
        return -1;
    }
    if (spec != NULL && spec != Py_None) {
        // __spec__ already set.
        Py_DECREF(spec);
        return 0;
    }
    Py_XDECREF(spec);
    if (CPyImport_InitSpecClasses() < 0) {
        return -1;
    }
    PyObject *file = NULL;
    if (PyObject_GetOptionalAttrString(modobj, "__file__", &file) < 0) {
        return -1;
    }
    if (file == NULL) {
        file = Py_None;
        Py_INCREF(file);
    }
    // ExtensionFileLoader(name, path)
    PyObject *loader = PyObject_CallFunctionObjArgs(
        CPyImport_ExtFileLoaderClass, module_name, file, NULL);
    if (loader == NULL) {
        Py_DECREF(file);
        return -1;
    }
    // ModuleSpec(name, loader, *, origin=file, is_package=is_package)
    PyObject *is_pkg_obj = is_package ? Py_True : Py_False;
    PyObject *spec_args[] = {NULL, module_name, loader, file, is_pkg_obj};
    spec = PyObject_Vectorcall(
        CPyImport_ModuleSpecClass,
        spec_args + 1,
        2 | PY_VECTORCALL_ARGUMENTS_OFFSET,
        CPyImport_SpecKwnames);
    Py_DECREF(file);
    if (spec == NULL) {
        Py_DECREF(loader);
        return -1;
    }
    if (is_package) {
        // Set submodule_search_locations from __path__
        PyObject *path = NULL;
        if (PyObject_GetOptionalAttrString(modobj, "__path__", &path) < 0) {
            Py_DECREF(spec);
            Py_DECREF(loader);
            return -1;
        }
        if (path != NULL) {
            if (PyObject_SetAttrString(spec, "submodule_search_locations", path) < 0) {
                Py_DECREF(path);
                Py_DECREF(spec);
                Py_DECREF(loader);
                return -1;
            }
            Py_DECREF(path);
        }
    }
    if (PyObject_SetAttrString(modobj, "__spec__", spec) < 0 ||
            PyObject_SetAttrString(modobj, "__loader__", loader) < 0) {
        Py_DECREF(spec);
        Py_DECREF(loader);
        return -1;
    }
    Py_DECREF(spec);
    Py_DECREF(loader);
    return 0;
}

PyObject *CPyImport_ImportNative(PyObject *module_name,
                                 PyObject *(*init_only_fn)(void),
                                 int (*exec_fn)(PyObject *),
                                 CPyModule **module_static,
                                 PyObject *shared_lib_file, PyObject *ext_suffix,
                                 Py_ssize_t is_package) {
    PyObject *parent_module = NULL;
    PyObject *child_name = NULL;
    PyObject *exc_type, *exc_val, *exc_tb;
    Py_ssize_t name_len = PyUnicode_GetLength(module_name);
    if (name_len < 0) {
        return NULL;
    }
    Py_ssize_t dot = PyUnicode_FindChar(module_name, '.', 0, name_len, -1);
    if (dot >= 0) {
        // Import the parent package first to preserve import ordering semantics.
        PyObject *parent_name = PyUnicode_Substring(module_name, 0, dot);
        if (parent_name == NULL) {
            CPyError_OutOfMemory();
        }
        child_name = PyUnicode_Substring(module_name, dot + 1, name_len);
        if (child_name == NULL) {
            CPyError_OutOfMemory();
        }
        parent_module = PyImport_Import(parent_name);
        Py_DECREF(parent_name);
        if (parent_module == NULL) {
            Py_DECREF(child_name);
            return NULL;
        }
    }

    // Create the module object without executing the module body.
    // CPyInitOnly_* uses an internal static to cache the module object.
    // We then check sys.modules to determine whether the module body
    // has already been executed (or is being executed in a circular import).
    PyObject *module_dict = PyImport_GetModuleDict();
    if (module_dict == NULL) {
        Py_XDECREF(parent_module);
        Py_XDECREF(child_name);
        return NULL;
    }

    PyObject *existing = PyDict_GetItemWithError(module_dict, module_name);
    if (existing != NULL) {
        if (*module_static != NULL) {
            if (existing == (PyObject *)*module_static) {
                Py_INCREF(existing);
                Py_XDECREF(parent_module);
                Py_XDECREF(child_name);
                return existing;
            }
            PyErr_Format(PyExc_ImportError,
                         "native module '%U' in sys.modules was replaced after initialization",
                         module_name);
            Py_XDECREF(parent_module);
            Py_XDECREF(child_name);
            return NULL;
        }
    }
    if (PyErr_Occurred()) {
        Py_XDECREF(parent_module);
        Py_XDECREF(child_name);
        return NULL;
    }

    PyObject *modobj = init_only_fn();
    if (modobj == NULL) {
        Py_XDECREF(parent_module);
        Py_XDECREF(child_name);
        return NULL;
    }

    if (PyObject_SetItem(module_dict, module_name, modobj) < 0) {
        goto fail;
    }

    if (*module_static != (CPyModule *)modobj) {
        PyErr_Format(PyExc_ImportError,
                     "native module '%U' was initialized inconsistently",
                     module_name);
        goto fail;
    }

    if (CPyImport_SetDunderAttrs(modobj, module_name, shared_lib_file, ext_suffix, is_package) < 0) {
        goto fail;
    }

    // Now execute the module body, with __file__ and __package__ already set.
    if (exec_fn(modobj) != 0) {
        goto fail;
    }

    // Match CPython import semantics: publish parent.child only after the
    // child module finished executing successfully.
    if (parent_module != NULL && PyObject_SetAttr(parent_module, child_name, modobj) < 0) {
        goto fail;
    }

    Py_XDECREF(parent_module);
    Py_XDECREF(child_name);
    return modobj;

fail:
    // Clean up on failure so that a subsequent import attempt will retry
    // initialization.
    PyErr_Fetch(&exc_type, &exc_val, &exc_tb);
    PyObject_DelItem(module_dict, module_name);
    PyErr_Clear();
    PyErr_Restore(exc_type, exc_val, exc_tb);
    Py_XDECREF(parent_module);
    Py_XDECREF(child_name);
    Py_CLEAR(*module_static);
    return NULL;
}

int CPyImport_SetDunderAttrs(PyObject *module, PyObject *module_name, PyObject *shared_lib_file,
                             PyObject *ext_suffix, Py_ssize_t is_package)
{
    int res = CPyImport_SetModulePackage(module, module_name, is_package);
    if (res < 0) {
        return res;
    }

    res = CPyImport_SetModuleFile(module, module_name, shared_lib_file, ext_suffix,
                                  is_package);
    if (res < 0) {
        return res;
    }

    if (is_package) {
        res = CPyImport_SetModulePath(module);
        if (res < 0) {
            return res;
        }
    }

    return CPyImport_SetModuleSpec(module, module_name, is_package);
}

#if CPY_3_14_FEATURES

#include "internal/pycore_object.h"

void CPy_SetImmortal(PyObject *obj) {
    _Py_SetImmortal(obj);
}

#endif
