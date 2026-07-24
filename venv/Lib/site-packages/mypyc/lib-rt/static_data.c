#ifndef STATIC_DATA
#define STATIC_DATA

#include "static_data.h"

// Adopted from numpy 2.4.0: numpy/_core/src/multiarry/npy_static_data.c

mypyc_interned_str_struct mypyc_interned_str;

#define INTERN_STRING(struct_member, string) \
    assert(mypyc_interned_str.struct_member == NULL); \
    mypyc_interned_str.struct_member = PyUnicode_InternFromString(string); \
    if (mypyc_interned_str.struct_member == NULL) { \
        return -1; \
    }

int
intern_strings(void) {
    if (mypyc_interned_str.values != NULL) {
        // Already interned.
        return 0;
    }
    INTERN_STRING(__init_subclass__, "__init_subclass__");
    INTERN_STRING(__module__, "__module__");
    INTERN_STRING(__mro_entries__, "__mro_entries__");
    INTERN_STRING(__mypyc_attrs__, "__mypyc_attrs__");
    INTERN_STRING(__name__, "__name__");
    INTERN_STRING(__orig_bases__, "__orig_bases__");
    INTERN_STRING(__qualname__, "__qualname__");
    INTERN_STRING(__slots__, "__slots__");
    INTERN_STRING(__radd__, "__radd__");
    INTERN_STRING(__rsub__, "__rsub__");
    INTERN_STRING(__rmul__, "__rmul__");
    INTERN_STRING(__rtruediv__, "__rtruediv__");
    INTERN_STRING(__rmod__, "__rmod__");
    INTERN_STRING(__rdivmod__, "__rdivmod__");
    INTERN_STRING(__rfloordiv__, "__rfloordiv__");
    INTERN_STRING(__rpow__, "__rpow__");
    INTERN_STRING(__rmatmul__, "__rmatmul__");
    INTERN_STRING(__rand__, "__rand__");
    INTERN_STRING(__ror__, "__ror__");
    INTERN_STRING(__rxor__, "__rxor__");
    INTERN_STRING(__rlshift__, "__rlshift__");
    INTERN_STRING(__rrshift__, "__rrshift__");
    INTERN_STRING(__eq__, "__eq__");
    INTERN_STRING(__ne__, "__ne__");
    INTERN_STRING(__gt__, "__gt__");
    INTERN_STRING(__le__, "__le__");
    INTERN_STRING(__lt__, "__lt__");
    INTERN_STRING(__ge__, "__ge__");
    INTERN_STRING(clear, "clear");
    INTERN_STRING(close_, "close");
    INTERN_STRING(copy, "copy");
    INTERN_STRING(dispatch_cache, "dispatch_cache");
    INTERN_STRING(endswith, "endswith");
    INTERN_STRING(get_type_hints, "get_type_hints");
    INTERN_STRING(keys, "keys");
    INTERN_STRING(lower, "lower");
    INTERN_STRING(items, "items");
    INTERN_STRING(join, "join");
    INTERN_STRING(register_, "register");
    INTERN_STRING(registry, "registry");
    INTERN_STRING(send, "send");
    INTERN_STRING(setdefault, "setdefault");
    INTERN_STRING(startswith, "startswith");
    INTERN_STRING(super, "super");
    INTERN_STRING(throw_, "throw");
    INTERN_STRING(translate, "translate");
    INTERN_STRING(update, "update");
    INTERN_STRING(upper, "upper");
    INTERN_STRING(values, "values");
    return 0;
}

#endif
