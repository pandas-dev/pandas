//////////////////// ArgTypeTest.proto ////////////////////
//@requires: ArgTypeTestFunc

// Exact is 0 (False), 1 (True) or 2 (True and from annotation)
// The latter gives a small amount of extra error diagnostics
#define __Pyx_ArgTypeTest(obj, type, none_allowed, name, exact) \
    ((likely(__Pyx_IS_TYPE(obj, type) | (none_allowed && (obj == Py_None)))) ? 1 : \
        __Pyx__ArgTypeTest(obj, type, name, exact))


//////////////////// ArgTypeTestFunc.export ////////////////////

static int __Pyx__ArgTypeTest(PyObject *obj, PyTypeObject *type, const char *name, int exact); /*proto*/

//////////////////// ArgTypeTestFunc ////////////////////

static int __Pyx__ArgTypeTest(PyObject *obj, PyTypeObject *type, const char *name, int exact)
{
    __Pyx_TypeName type_name;
    __Pyx_TypeName obj_type_name;
    PyObject *extra_info = EMPTY(unicode);
    int from_annotation_subclass = 0;
    if (unlikely(!type)) {
        PyErr_SetString(PyExc_SystemError, "Missing type object");
        return 0;
    }
    else if (!exact) {
        if (likely(__Pyx_TypeCheck(obj, type))) return 1;
    } else if (exact == 2) {
        // type from annotation
        if (__Pyx_TypeCheck(obj, type)) {
            from_annotation_subclass = 1;
            extra_info = PYUNICODE("Note that Cython is deliberately stricter than PEP-484 and rejects subclasses of builtin types. If you need to pass subclasses then set the 'annotation_typing' directive to False.");
        }
    }
    type_name = __Pyx_PyType_GetFullyQualifiedName(type);
    obj_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(obj));
    PyErr_Format(PyExc_TypeError,
        "Argument '%.200s' has incorrect type (expected " __Pyx_FMT_TYPENAME
        ", got " __Pyx_FMT_TYPENAME ")"
#if __PYX_LIMITED_VERSION_HEX < 0x030C0000
        "%s%U"
#endif
        , name, type_name, obj_type_name
#if __PYX_LIMITED_VERSION_HEX < 0x030C0000
        , (from_annotation_subclass ? ". " : ""), extra_info
#endif
        );
#if __PYX_LIMITED_VERSION_HEX >= 0x030C0000
    // Set the extra_info as a note instead. In principle it'd be possible to do this
    // from Python 3.11 up, but PyErr_GetRaisedException makes it much easier so do it
    // from Python 3.12 instead.
    if (exact == 2 && from_annotation_subclass) {
        PyObject *res;
        PyObject *vargs[2];
        vargs[0] = PyErr_GetRaisedException();
        vargs[1] = extra_info;
        res = PyObject_VectorcallMethod(PYUNICODE("add_note"), vargs, 2, NULL);
        Py_XDECREF(res);
        PyErr_SetRaisedException(vargs[0]);
    }
#endif
    __Pyx_DECREF_TypeName(type_name);
    __Pyx_DECREF_TypeName(obj_type_name);
    return 0;
}


//////////////////// RaiseArgTupleInvalid.proto ////////////////////

static void __Pyx_RaiseArgtupleInvalid(const char* func_name, int exact,
    Py_ssize_t num_min, Py_ssize_t num_max, Py_ssize_t num_found); /*proto*/

//////////////////// RaiseArgTupleInvalid ////////////////////

//  __Pyx_RaiseArgtupleInvalid raises the correct exception when too
//  many or too few positional arguments were found.  This handles
//  Py_ssize_t formatting correctly.

static void __Pyx_RaiseArgtupleInvalid(
    const char* func_name,
    int exact,
    Py_ssize_t num_min,
    Py_ssize_t num_max,
    Py_ssize_t num_found)
{
    Py_ssize_t num_expected;
    const char *more_or_less;

    if (num_found < num_min) {
        num_expected = num_min;
        more_or_less = "at least";
    } else {
        num_expected = num_max;
        more_or_less = "at most";
    }
    if (exact) {
        more_or_less = "exactly";
    }
    PyErr_Format(PyExc_TypeError,
                 "%.200s() takes %.8s %" CYTHON_FORMAT_SSIZE_T "d positional argument%.1s (%" CYTHON_FORMAT_SSIZE_T "d given)",
                 func_name, more_or_less, num_expected,
                 (num_expected == 1) ? "" : "s", num_found);
}


//////////////////// RaiseKeywordRequired.proto ////////////////////

static void __Pyx_RaiseKeywordRequired(const char* func_name, PyObject* kw_name); /*proto*/

//////////////////// RaiseKeywordRequired ////////////////////

static void __Pyx_RaiseKeywordRequired(const char* func_name, PyObject* kw_name) {
    PyErr_Format(PyExc_TypeError,
        "%s() needs keyword-only argument %U", func_name, kw_name);
}


//////////////////// RaiseDoubleKeywords.proto ////////////////////

static void __Pyx_RaiseDoubleKeywordsError(const char* func_name, PyObject* kw_name); /*proto*/

//////////////////// RaiseDoubleKeywords ////////////////////

static void __Pyx_RaiseDoubleKeywordsError(
    const char* func_name,
    PyObject* kw_name)
{
    PyErr_Format(PyExc_TypeError,
        "%s() got multiple values for keyword argument '%U'", func_name, kw_name);
}


//////////////////// RaiseMappingExpected.proto ////////////////////

static void __Pyx_RaiseMappingExpectedError(PyObject* arg); /*proto*/

//////////////////// RaiseMappingExpected ////////////////////

static void __Pyx_RaiseMappingExpectedError(PyObject* arg) {
    __Pyx_TypeName arg_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(arg));
    PyErr_Format(PyExc_TypeError,
        "'" __Pyx_FMT_TYPENAME "' object is not a mapping", arg_type_name);
    __Pyx_DECREF_TypeName(arg_type_name);
}


//////////////////// KeywordStringCheck.proto ////////////////////

static CYTHON_INLINE int __Pyx_CheckKeywordStrings(const char* function_name, PyObject *kw); /*proto*/

//////////////////// KeywordStringCheck ////////////////////

// __Pyx_CheckKeywordStrings raises an error if non-string keywords
// were passed to a function.
//
// The "kw" argument is either a dict (for METH_VARARGS) or a tuple
// (for METH_FASTCALL), both non-empty.

static int __Pyx_CheckKeywordStrings(
    const char* function_name,
    PyObject *kw)
{
    // PyPy appears to check keyword types at call time, not at unpacking time.
#if CYTHON_COMPILING_IN_PYPY && !defined(PyArg_ValidateKeywordArguments)
    CYTHON_UNUSED_VAR(function_name);
    CYTHON_UNUSED_VAR(kw);
    return 0;
#else

    // Validate keyword types.
    if (CYTHON_METH_FASTCALL && likely(PyTuple_Check(kw))) {
        // On CPython >= 3.9, the FASTCALL protocol guarantees that keyword
        // names are strings (see https://github.com/python/cpython/issues/81721)
#if PY_VERSION_HEX >= 0x03090000
        CYTHON_UNUSED_VAR(function_name);
#else

        Py_ssize_t kwsize;
        #if CYTHON_ASSUME_SAFE_SIZE
        kwsize = PyTuple_GET_SIZE(kw);
        #else
        kwsize = PyTuple_Size(kw);
        if (unlikely(kwsize < 0)) return -1;
        #endif

        for (Py_ssize_t pos = 0; pos < kwsize; pos++) {
            PyObject* key = NULL;
            #if CYTHON_ASSUME_SAFE_MACROS
            key = PyTuple_GET_ITEM(kw, pos);
            #else
            key = PyTuple_GetItem(kw, pos);
            if (unlikely(!key)) return -1;
            #endif

            if (unlikely(!PyUnicode_Check(key))) {
                PyErr_Format(PyExc_TypeError,
                    "%.200s() keywords must be strings", function_name);
                return -1;
            }
        }
#endif
    } else {
        // Otherwise, 'kw' is a dict: check if it's unicode-keys-only and let Python set the error otherwise.
        if (unlikely(!PyArg_ValidateKeywordArguments(kw))) return -1;
    }

    return 0;
#endif
}


//////////////////// RejectKeywords.export ////////////////////

static void __Pyx_RejectKeywords(const char* function_name, PyObject *kwds); /*proto*/

//////////////////// RejectKeywords ////////////////////
//@requires: ObjectHandling.c::OwnedDictNext

static void __Pyx_RejectKeywords(const char* function_name, PyObject *kwds) {
    // Get the first keyword argument (there is at least one) and raise a TypeError for it.
    PyObject *key = NULL;
    if (CYTHON_METH_FASTCALL && likely(PyTuple_Check(kwds))) {
        key = __Pyx_PySequence_ITEM(kwds, 0);
    } else {
#if CYTHON_AVOID_BORROWED_REFS
        PyObject *pos = NULL;
#else
        Py_ssize_t pos = 0;
#endif
#if !CYTHON_COMPILING_IN_PYPY || defined(PyArg_ValidateKeywordArguments)
        // Check if dict is unicode-keys-only and let Python set the error otherwise.
        if (unlikely(!PyArg_ValidateKeywordArguments(kwds))) return;
#endif
        // Read first key.
        __Pyx_PyDict_NextRef(kwds, &pos, &key, NULL);
#if CYTHON_AVOID_BORROWED_REFS
        Py_XDECREF(pos);
#endif
    }

    if (likely(key)) {
        PyErr_Format(PyExc_TypeError,
            "%s() got an unexpected keyword argument '%U'",
            function_name, key);
        Py_DECREF(key);
    }
}


//////////////////// ParseKeywords.proto ////////////////////
//@requires: ParseKeywordsImpl

static CYTHON_INLINE int __Pyx_ParseKeywords(
    PyObject *kwds, PyObject *const *kwvalues, PyObject ** const argnames[],
    PyObject *kwds2, PyObject *values[],
    Py_ssize_t num_pos_args, Py_ssize_t num_kwargs,
    const char* function_name,
    int ignore_unknown_kwargs
); /*proto*/

//////////////////// ParseKeywords ////////////////////
//  __Pyx_ParseOptionalKeywords copies the optional/unknown keyword
//  arguments from kwds into the dict kwds2.  If kwds2 is NULL, unknown
//  keywords will raise an invalid keyword error.
//
//  When not using METH_FASTCALL, kwds is a dict and kwvalues is NULL.
//  Otherwise, kwds is a tuple with keyword names and kwvalues is a C
//  array with the corresponding values.
//
//  Three kinds of errors are checked: 1) non-string keywords, 2)
//  unexpected keywords and 3) overlap with positional arguments.
//
//  If num_posargs is greater 0, it denotes the number of positional
//  arguments that were passed and that must therefore not appear
//  amongst the keywords as well.
//
//  This method does not check for required keyword arguments.

static int __Pyx_ParseKeywords(
    PyObject *kwds,
    PyObject * const *kwvalues,
    PyObject ** const argnames[],
    PyObject *kwds2,
    PyObject *values[],
    Py_ssize_t num_pos_args,
    Py_ssize_t num_kwargs,
    const char* function_name,
    int ignore_unknown_kwargs)
{
    // Only called if kwds contains at least one optional keyword argument.
    if (CYTHON_METH_FASTCALL && likely(PyTuple_Check(kwds)))
        return __Pyx_ParseKeywordsTuple(kwds, kwvalues, argnames, kwds2, values, num_pos_args, num_kwargs, function_name, ignore_unknown_kwargs);
    else if (kwds2)
        return __Pyx_ParseKeywordDictToDict(kwds, argnames, kwds2, values, num_pos_args, function_name);
    else
        return __Pyx_ParseKeywordDict(kwds, argnames, values, num_pos_args, num_kwargs, function_name, ignore_unknown_kwargs);
}


//////////////////// ParseKeywordsImpl.export ////////////////////
//@requires: RaiseDoubleKeywords
//@requires: Synchronization.c::CriticalSections
//@requires: ObjectHandling.c::OwnedDictNext


static int __Pyx_ParseKeywordsTuple(
    PyObject *kwds,
    PyObject * const *kwvalues,
    PyObject ** const argnames[],
    PyObject *kwds2,
    PyObject *values[],
    Py_ssize_t num_pos_args,
    Py_ssize_t num_kwargs,
    const char* function_name,
    int ignore_unknown_kwargs
);

static int __Pyx_ParseKeywordDictToDict(
    PyObject *kwds,
    PyObject ** const argnames[],
    PyObject *kwds2,
    PyObject *values[],
    Py_ssize_t num_pos_args,
    const char* function_name
);

static int __Pyx_ParseKeywordDict(
    PyObject *kwds,
    PyObject ** const argnames[],
    PyObject *values[],
    Py_ssize_t num_pos_args,
    Py_ssize_t num_kwargs,
    const char* function_name,
    int ignore_unknown_kwargs
);

//////////////////// ParseKeywordsImpl ////////////////////

static int __Pyx_ValidateDuplicatePosArgs(
    PyObject *kwds,
    PyObject ** const argnames[],
    PyObject ** const *first_kw_arg,
    const char* function_name)
{
    PyObject ** const *name = argnames;
    while (name != first_kw_arg) {
        PyObject *key = **name;
        int found = PyDict_Contains(kwds, key);
        if (unlikely(found)) {
            if (found == 1) __Pyx_RaiseDoubleKeywordsError(function_name, key);
            goto bad;
        }
        name++;
    }
    return 0;

bad:
    return -1;
}

#if CYTHON_USE_UNICODE_INTERNALS
static CYTHON_INLINE int __Pyx_UnicodeKeywordsEqual(PyObject *s1, PyObject *s2) {
    int kind;
    Py_ssize_t len = PyUnicode_GET_LENGTH(s1);
    if (len != PyUnicode_GET_LENGTH(s2)) return 0;

    kind = PyUnicode_KIND(s1);
    if (kind != PyUnicode_KIND(s2)) return 0;

    const void *data1 = PyUnicode_DATA(s1);
    const void *data2 = PyUnicode_DATA(s2);
    return (memcmp(data1, data2, (size_t) len * (size_t) kind) == 0);
}
#endif

static int __Pyx_MatchKeywordArg_str(
    PyObject *key,
    PyObject ** const argnames[],
    PyObject ** const *first_kw_arg,
    size_t *index_found,
    const char *function_name)
{
    PyObject ** const *name;
    #if CYTHON_USE_UNICODE_INTERNALS
    // The key hash is probably pre-calculated.
    Py_hash_t key_hash = ((PyASCIIObject*)key)->hash;
    if (unlikely(key_hash == -1)) {
        key_hash = PyObject_Hash(key);
        if (unlikely(key_hash == -1))
            goto bad;
    }
    #endif

    // Compare strings for non-interned matches.
    name = first_kw_arg;
    while (*name) {
        PyObject *name_str = **name;

        #if CYTHON_USE_UNICODE_INTERNALS
        // The hash value of our interned argument name is definitely pre-calculated.
        if (key_hash == ((PyASCIIObject*)name_str)->hash && __Pyx_UnicodeKeywordsEqual(name_str, key)) {
            *index_found = (size_t) (name - argnames);
            return 1;
        }
        #else

        #if CYTHON_ASSUME_SAFE_SIZE
        if (PyUnicode_GET_LENGTH(name_str) == PyUnicode_GET_LENGTH(key))
        #endif
        {
            int cmp = PyUnicode_Compare(name_str, key);
            if (cmp < 0 && unlikely(PyErr_Occurred())) goto bad;
            if (cmp == 0) {
                *index_found = (size_t) (name - argnames);
                return 1;
            }
        }
        #endif
        name++;
    }

    // Not found in keyword parameters, check for (unlikely) duplicate positional argument.
    name = argnames;
    while (name != first_kw_arg) {
        PyObject *name_str = **name;

        #if CYTHON_USE_UNICODE_INTERNALS
        if (unlikely(key_hash == ((PyASCIIObject*)name_str)->hash)) {
            if (__Pyx_UnicodeKeywordsEqual(name_str, key))
                goto arg_passed_twice;
        }

        #else

        #if CYTHON_ASSUME_SAFE_SIZE
        if (PyUnicode_GET_LENGTH(name_str) == PyUnicode_GET_LENGTH(key))
        #endif
        {
            if (unlikely(name_str == key)) goto arg_passed_twice;
            int cmp = PyUnicode_Compare(name_str, key);
            if (cmp < 0 && unlikely(PyErr_Occurred())) goto bad;
            if (cmp == 0) goto arg_passed_twice;
        }

        #endif
        name++;
    }

    return 0;

arg_passed_twice:
    __Pyx_RaiseDoubleKeywordsError(function_name, key);
    goto bad;
bad:
    return -1;
}

static int __Pyx_MatchKeywordArg_nostr(
    PyObject *key,
    PyObject ** const argnames[],
    PyObject ** const *first_kw_arg,
    size_t *index_found,
    const char *function_name)
{
    // Conservatively handle str subclasses.
    PyObject ** const *name;

    if (unlikely(!PyUnicode_Check(key))) goto invalid_keyword_type;

    // Match keyword argument names.
    name = first_kw_arg;
    while (*name) {
        int cmp = PyObject_RichCompareBool(**name, key, Py_EQ);
        if (cmp == 1) {
            *index_found = (size_t) (name - argnames);
            return 1;
        }
        if (unlikely(cmp == -1)) goto bad;
        name++;
    }
    // Reject collisions with positional arguments.
    name = argnames;
    while (name != first_kw_arg) {
        int cmp = PyObject_RichCompareBool(**name, key, Py_EQ);
        if (unlikely(cmp != 0)) {
            if (cmp == 1) goto arg_passed_twice;
            else goto bad;
        }
        name++;
    }
    return 0;

arg_passed_twice:
    __Pyx_RaiseDoubleKeywordsError(function_name, key);
    goto bad;
invalid_keyword_type:
    PyErr_Format(PyExc_TypeError,
        "%.200s() keywords must be strings", function_name);
    goto bad;
bad:
    return -1;
}

static CYTHON_INLINE int __Pyx_MatchKeywordArg(
    PyObject *key,
    PyObject ** const argnames[],
    PyObject ** const *first_kw_arg,
    size_t *index_found,
    const char *function_name)
{
    // Optimise for plain str behaviour.
    return likely(PyUnicode_CheckExact(key)) ?
        __Pyx_MatchKeywordArg_str(key, argnames, first_kw_arg, index_found, function_name) :
        __Pyx_MatchKeywordArg_nostr(key, argnames, first_kw_arg, index_found, function_name);
}

static void __Pyx_RejectUnknownKeyword(
    PyObject *kwds,
    PyObject ** const argnames[],
    PyObject ** const *first_kw_arg,
    const char *function_name)
{
    // Find the first unknown keyword and raise an error. There must be at least one.
    #if CYTHON_AVOID_BORROWED_REFS
    PyObject *pos = NULL;
    #else
    Py_ssize_t pos = 0;
    #endif

    PyObject *key = NULL;

    __Pyx_BEGIN_CRITICAL_SECTION(kwds);
    while (
        #if CYTHON_AVOID_BORROWED_REFS
        __Pyx_PyDict_NextRef(kwds, &pos, &key, NULL)
        #else
        PyDict_Next(kwds, &pos, &key, NULL)
        #endif
    ) {
        // Quickly exclude the 'obviously' valid/known keyword arguments (exact pointer match).
        PyObject** const *name = first_kw_arg;
        while (*name && (**name != key)) name++;

        if (!*name) {
            // No exact match found:
            // compare against positional (always reject) and keyword (reject unknown) names.

            size_t index_found = 0;
            int cmp = __Pyx_MatchKeywordArg(key, argnames, first_kw_arg, &index_found, function_name);
            if (cmp != 1) {
                if (cmp == 0) {
                    PyErr_Format(PyExc_TypeError,
                        "%s() got an unexpected keyword argument '%U'",
                        function_name, key);
                }
                #if CYTHON_AVOID_BORROWED_REFS
                Py_DECREF(key);
                #endif

                break;
            }
        }
        #if CYTHON_AVOID_BORROWED_REFS
        Py_DECREF(key);
        #endif
    }
    __Pyx_END_CRITICAL_SECTION();
    #if CYTHON_AVOID_BORROWED_REFS
    Py_XDECREF(pos);
    #endif

    assert(PyErr_Occurred());
}

static int __Pyx_ParseKeywordDict(
    PyObject *kwds,
    PyObject ** const argnames[],
    PyObject *values[],
    Py_ssize_t num_pos_args,
    Py_ssize_t num_kwargs,
    const char* function_name,
    int ignore_unknown_kwargs)
{
    PyObject** const *name;
    PyObject** const *first_kw_arg = argnames + num_pos_args;
    Py_ssize_t extracted = 0;

#if !CYTHON_COMPILING_IN_PYPY || defined(PyArg_ValidateKeywordArguments)
    // Check if dict is unicode-keys-only and let Python set the error otherwise.
    if (unlikely(!PyArg_ValidateKeywordArguments(kwds))) return -1;
#endif

    // Extract declared keyword arguments.
    name = first_kw_arg;
    while (*name && num_kwargs > extracted) {
        PyObject * key = **name;
        PyObject *value;
        int found = 0;

        #if __PYX_LIMITED_VERSION_HEX >= 0x030d0000
        found = PyDict_GetItemRef(kwds, key, &value);
        #else
        value = PyDict_GetItemWithError(kwds, key);
        if (value) {
            Py_INCREF(value);
            found = 1;
        } else {
            if (unlikely(PyErr_Occurred())) goto bad;
        }
        #endif

        if (found) {
            if (unlikely(found < 0)) goto bad;
            values[name-argnames] = value;
            extracted++;
        }

        name++;
    }

    if (num_kwargs > extracted) {
        if (ignore_unknown_kwargs) {
            // Make sure the remaining kwargs are not duplicate posargs.
            if (unlikely(__Pyx_ValidateDuplicatePosArgs(kwds, argnames, first_kw_arg, function_name) == -1))
                goto bad;
        } else {
            // Any remaining kwarg is an error.
            __Pyx_RejectUnknownKeyword(kwds, argnames, first_kw_arg, function_name);
            goto bad;
        }
    }
    return 0;

bad:
    return -1;
}

static int __Pyx_ParseKeywordDictToDict(
    PyObject *kwds,
    PyObject ** const argnames[],
    PyObject *kwds2,
    PyObject *values[],
    Py_ssize_t num_pos_args,
    const char* function_name)
{
    // Validate and parse keyword arguments from kwds dict.
    PyObject** const *name;
    PyObject** const *first_kw_arg = argnames + num_pos_args;
    Py_ssize_t len;

#if !CYTHON_COMPILING_IN_PYPY || defined(PyArg_ValidateKeywordArguments)
    // Check if dict is unicode-keys-only and let Python set the error otherwise.
    if (unlikely(!PyArg_ValidateKeywordArguments(kwds))) return -1;
#endif

    // Fast copy of all kwargs.
    if (PyDict_Update(kwds2, kwds) < 0) goto bad;

    // Extract declared keyword arguments (if any).
    name = first_kw_arg;
    while (*name) {
        PyObject *key = **name;
        PyObject *value;

#if !CYTHON_COMPILING_IN_LIMITED_API && (PY_VERSION_HEX >= 0x030d00A2 || defined(PyDict_Pop))
        int found = PyDict_Pop(kwds2, key, &value);
        if (found) {
            if (unlikely(found < 0)) goto bad;
            values[name-argnames] = value;
        }
#elif __PYX_LIMITED_VERSION_HEX >= 0x030d0000
        int found = PyDict_GetItemRef(kwds2, key, &value);
        if (found) {
            if (unlikely(found < 0)) goto bad;
            values[name-argnames] = value;
            if (unlikely(PyDict_DelItem(kwds2, key) < 0)) goto bad;
        }
#else
    // We use 'kwds2' as sentinel value to dict.pop() to avoid an exception on missing key.
    #if CYTHON_COMPILING_IN_CPYTHON
        value = _PyDict_Pop(kwds2, key, kwds2);
    #else
        value = CALL_UNBOUND_METHOD(PyDict_Type, "pop", kwds2, key, kwds2);
    #endif
        if (value == kwds2) {
            // Not found.
            Py_DECREF(value);
        } else {
            if (unlikely(!value)) goto bad;
            values[name-argnames] = value;
        }
#endif
        name++;
    }

    // If unmatched keywords remain, check for duplicates of positional arguments.
    len = PyDict_Size(kwds2);
    if (len > 0) {
        return __Pyx_ValidateDuplicatePosArgs(kwds, argnames, first_kw_arg, function_name);
    } else if (unlikely(len == -1)) {
        goto bad;
    }

    return 0;

bad:
    return -1;
}

static int __Pyx_ParseKeywordsTuple(
    PyObject *kwds,
    PyObject * const *kwvalues,
    PyObject ** const argnames[],
    PyObject *kwds2,
    PyObject *values[],
    Py_ssize_t num_pos_args,
    Py_ssize_t num_kwargs,
    const char* function_name,
    int ignore_unknown_kwargs)
{
    PyObject *key = NULL;
    PyObject** const * name;
    PyObject** const *first_kw_arg = argnames + num_pos_args;

    for (Py_ssize_t pos = 0; pos < num_kwargs; pos++) {
#if CYTHON_AVOID_BORROWED_REFS
        key = __Pyx_PySequence_ITEM(kwds, pos);
#else
        key = __Pyx_PyTuple_GET_ITEM(kwds, pos);
#endif
#if !CYTHON_ASSUME_SAFE_MACROS
        if (unlikely(!key)) goto bad;
#endif

        // Quick pointer search for interned parameter matches (will usually succeed).
        name = first_kw_arg;
        while (*name && (**name != key)) name++;
        if (*name) {
            // Declared keyword: **name == key
            PyObject *value = kwvalues[pos];
            values[name-argnames] = __Pyx_NewRef(value);
        } else {
            size_t index_found = 0;
            int cmp = __Pyx_MatchKeywordArg(key, argnames, first_kw_arg, &index_found, function_name);

            if (cmp == 1) {
                // Found in declared keywords => assign value.
                PyObject *value = kwvalues[pos];
                values[index_found] = __Pyx_NewRef(value);
            } else {
                if (unlikely(cmp == -1)) goto bad;
                if (kwds2) {
                    PyObject *value = kwvalues[pos];
                    if (unlikely(PyDict_SetItem(kwds2, key, value))) goto bad;
                } else if (!ignore_unknown_kwargs) {
                    goto invalid_keyword;
                }
            }
        }

        #if CYTHON_AVOID_BORROWED_REFS
        Py_DECREF(key);
        key = NULL;
        #endif
    }
    return 0;

invalid_keyword:
    PyErr_Format(PyExc_TypeError,
        "%s() got an unexpected keyword argument '%U'",
        function_name, key);
    goto bad;
bad:
    #if CYTHON_AVOID_BORROWED_REFS
    Py_XDECREF(key);
    #endif
    return -1;
}


//////////////////// MergeKeywords.proto ////////////////////

static int __Pyx_MergeKeywords(PyObject *kwdict, PyObject *source_mapping); /*proto*/

//////////////////// MergeKeywords ////////////////////
//@requires: RaiseDoubleKeywords
//@requires: Optimize.c::dict_iter
//@requires: ObjectHandling.c::OwnedDictNext

static int __Pyx_MergeKeywords_dict(PyObject *kwdict, PyObject *source_dict) {
    Py_ssize_t len1, len2;

    len2 = PyDict_Size(source_dict);
    if (unlikely(len2 == -1)) return -1;
    if (len2 == 0) {
        // There's nothing to copy from an empty dict.
        return 0;
    }

    len1 = PyDict_Size(kwdict);
    if (unlikely(len1 == -1)) return -1;

    if (len1 > 0) {
        PyObject *key, *smaller_dict, *larger_dict;
        #if CYTHON_AVOID_BORROWED_REFS
        PyObject *ppos = NULL;
        #else
        Py_ssize_t ppos = 0;
        #endif
        int duplicates_found = 0;

        if (len1 <= len2) {
            smaller_dict = kwdict;
            larger_dict = source_dict;
        } else {
            smaller_dict = source_dict;
            larger_dict = kwdict;
        }

        __Pyx_BEGIN_CRITICAL_SECTION(smaller_dict);
        while (
            #if CYTHON_AVOID_BORROWED_REFS || CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
            __Pyx_PyDict_NextRef(smaller_dict, &ppos, &key, NULL)
            #else
            PyDict_Next(smaller_dict, &ppos, &key, NULL)
            #endif
        ) {
            if (unlikely(PyDict_Contains(larger_dict, key))) {
                __Pyx_RaiseDoubleKeywordsError("function", key);
                #if CYTHON_AVOID_BORROWED_REFS || CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
                Py_DECREF(key);
                #endif
                duplicates_found = 1;
                break;
            };
            #if CYTHON_AVOID_BORROWED_REFS || CYTHON_AVOID_THREAD_UNSAFE_BORROWED_REFS
            Py_DECREF(key);
            #endif
        }
        __Pyx_END_CRITICAL_SECTION();
        #if CYTHON_AVOID_BORROWED_REFS
        Py_XDECREF(ppos);
        #endif

        if (unlikely(duplicates_found))
            return -1;
    }

    return PyDict_Update(kwdict, source_dict);
}

static int __Pyx_MergeKeywords_any(PyObject *kwdict, PyObject *source_mapping) {
    PyObject *iter, *key = NULL, *value = NULL;
    int source_is_dict, result;
    Py_ssize_t orig_length, ppos = 0;

    iter = __Pyx_dict_iterator(source_mapping, 0, PYIDENT("items"), &orig_length, &source_is_dict);
    if (unlikely(!iter)) {
        // slow fallback: try converting to dict, then iterate
        PyObject *args;
        if (unlikely(!PyErr_ExceptionMatches(PyExc_AttributeError))) goto bad;
        PyErr_Clear();
        args = PyTuple_Pack(1, source_mapping);
        if (likely(args)) {
            PyObject *fallback = PyObject_Call((PyObject*)&PyDict_Type, args, NULL);
            Py_DECREF(args);
            if (likely(fallback)) {
                result = __Pyx_MergeKeywords_dict(kwdict, fallback);
                Py_DECREF(fallback);
                return result;
            }
        }
        goto bad;
    }

    while (1) {
        result = __Pyx_dict_iter_next(iter, orig_length, &ppos, &key, &value, NULL, source_is_dict);
        if (unlikely(result < 0)) goto bad;
        if (!result) break;

    #if PY_VERSION_HEX >= 0x030d0000 && !CYTHON_COMPILING_IN_LIMITED_API
        {
            int inserted = PyDict_SetDefaultRef(kwdict, key, value, NULL);
            if (unlikely(inserted != 0)) {
                if (inserted == 1) __Pyx_RaiseDoubleKeywordsError("function", key);
                result = -1;
            }
        }
    #else
        if (unlikely(PyDict_Contains(kwdict, key))) {
            __Pyx_RaiseDoubleKeywordsError("function", key);
            result = -1;
        } else {
            result = PyDict_SetItem(kwdict, key, value);
        }
    #endif

        Py_DECREF(key);
        Py_DECREF(value);
        if (unlikely(result < 0)) goto bad;
    }
    Py_XDECREF(iter);
    return 0;

bad:
    Py_XDECREF(iter);
    return -1;
}

static CYTHON_INLINE int __Pyx_MergeKeywords(PyObject *kwdict, PyObject *source_mapping) {
    assert(PyDict_Check(kwdict));
    if (likely(PyDict_Check(source_mapping))) {
        return __Pyx_MergeKeywords_dict(kwdict, source_mapping);
    } else {
        return __Pyx_MergeKeywords_any(kwdict, source_mapping);
    }
}


/////////////// fastcall.proto ///////////////

// We define various functions and macros with two variants:
//..._FASTCALL and ..._VARARGS

// The first is used when METH_FASTCALL is enabled and the second is used
// otherwise. If the Python implementation does not support METH_FASTCALL
// (because it's an old version of CPython or it's not CPython at all),
// then the ..._FASTCALL macros simply alias ..._VARARGS

#if CYTHON_AVOID_BORROWED_REFS
    #define __Pyx_ArgRef_VARARGS(args, i) __Pyx_PySequence_ITEM(args, i)
#elif CYTHON_ASSUME_SAFE_MACROS
    #define __Pyx_ArgRef_VARARGS(args, i) __Pyx_NewRef(__Pyx_PyTuple_GET_ITEM(args, i))
#else
    #define __Pyx_ArgRef_VARARGS(args, i) __Pyx_XNewRef(PyTuple_GetItem(args, i))
#endif

#define __Pyx_NumKwargs_VARARGS(kwds) PyDict_Size(kwds)
#define __Pyx_KwValues_VARARGS(args, nargs) NULL
#define __Pyx_GetKwValue_VARARGS(kw, kwvalues, s) __Pyx_PyDict_GetItemStrWithError(kw, s)
#define __Pyx_KwargsAsDict_VARARGS(kw, kwvalues) PyDict_Copy(kw)
#if CYTHON_METH_FASTCALL
    #define __Pyx_ArgRef_FASTCALL(args, i) __Pyx_NewRef(args[i])
    #define __Pyx_NumKwargs_FASTCALL(kwds) __Pyx_PyTuple_GET_SIZE(kwds)
    #define __Pyx_KwValues_FASTCALL(args, nargs) ((args) + (nargs))
    static CYTHON_INLINE PyObject * __Pyx_GetKwValue_FASTCALL(PyObject *kwnames, PyObject *const *kwvalues, PyObject *s);
  #if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX >= 0x030d0000 || CYTHON_COMPILING_IN_LIMITED_API
    CYTHON_UNUSED static PyObject *__Pyx_KwargsAsDict_FASTCALL(PyObject *kwnames, PyObject *const *kwvalues);/*proto*/
  #else
    #define __Pyx_KwargsAsDict_FASTCALL(kw, kwvalues) _PyStack_AsDict(kwvalues, kw)
  #endif
#else
    #define __Pyx_ArgRef_FASTCALL __Pyx_ArgRef_VARARGS
    #define __Pyx_NumKwargs_FASTCALL __Pyx_NumKwargs_VARARGS
    #define __Pyx_KwValues_FASTCALL __Pyx_KwValues_VARARGS
    #define __Pyx_GetKwValue_FASTCALL __Pyx_GetKwValue_VARARGS
    #define __Pyx_KwargsAsDict_FASTCALL __Pyx_KwargsAsDict_VARARGS
#endif

#define __Pyx_ArgsSlice_VARARGS(args, start, stop) PyTuple_GetSlice(args, start, stop)

#if CYTHON_METH_FASTCALL || (CYTHON_COMPILING_IN_CPYTHON && CYTHON_ASSUME_SAFE_MACROS && !CYTHON_AVOID_BORROWED_REFS)
#define __Pyx_ArgsSlice_FASTCALL(args, start, stop) __Pyx_PyTuple_FromArray(args + start, stop - start)
#else
#define __Pyx_ArgsSlice_FASTCALL(args, start, stop) PyTuple_GetSlice(args, start, stop)
#endif

/////////////// fastcall ///////////////
//@requires: ObjectHandling.c::TupleAndListFromArray
//@requires: StringTools.c::UnicodeEquals

#if CYTHON_METH_FASTCALL
// kwnames: tuple with names of keyword arguments
// kwvalues: C array with values of keyword arguments
// s: str with the keyword name to look for
static CYTHON_INLINE PyObject * __Pyx_GetKwValue_FASTCALL(PyObject *kwnames, PyObject *const *kwvalues, PyObject *s)
{
    // Search the kwnames array for s and return the corresponding value.
    // We do two loops: a first one to compare pointers (which will find a
    // match if the name in kwnames is interned, given that s is interned
    // by Cython). A second loop compares the actual strings.
    Py_ssize_t i, n = __Pyx_PyTuple_GET_SIZE(kwnames);
    #if !CYTHON_ASSUME_SAFE_SIZE
    if (unlikely(n == -1)) return NULL;
    #endif
    for (i = 0; i < n; i++)
    {
        PyObject *namei = __Pyx_PyTuple_GET_ITEM(kwnames, i);
        #if !CYTHON_ASSUME_SAFE_MACROS
        if (unlikely(!namei)) return NULL;
        #endif
        if (s == namei) return kwvalues[i];
    }
    for (i = 0; i < n; i++)
    {
        PyObject *namei = __Pyx_PyTuple_GET_ITEM(kwnames, i);
        #if !CYTHON_ASSUME_SAFE_MACROS
        if (unlikely(!namei)) return NULL;
        #endif
        int eq = __Pyx_PyUnicode_Equals(s, namei, Py_EQ);
        if (unlikely(eq != 0)) {
            if (unlikely(eq < 0)) return NULL;  /* error */
            return kwvalues[i];
        }
    }
    return NULL;  /* not found (no exception set) */
}

#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX >= 0x030d0000 || CYTHON_COMPILING_IN_LIMITED_API
CYTHON_UNUSED static PyObject *__Pyx_KwargsAsDict_FASTCALL(PyObject *kwnames, PyObject *const *kwvalues) {
    Py_ssize_t i, nkwargs;
    PyObject *dict;
#if !CYTHON_ASSUME_SAFE_SIZE
    nkwargs = PyTuple_Size(kwnames);
    if (unlikely(nkwargs < 0)) return NULL;
#else
    nkwargs = PyTuple_GET_SIZE(kwnames);
#endif

    dict = PyDict_New();
    if (unlikely(!dict))
        return NULL;

    for (i=0; i<nkwargs; i++) {
#if !CYTHON_ASSUME_SAFE_MACROS
        PyObject *key = PyTuple_GetItem(kwnames, i);
        if (!key) goto bad;
#else
        PyObject *key = PyTuple_GET_ITEM(kwnames, i);
#endif
        if (unlikely(PyDict_SetItem(dict, key, kwvalues[i]) < 0))
            goto bad;
    }
    return dict;

bad:
    Py_DECREF(dict);
    return NULL;
}
#endif

#endif
