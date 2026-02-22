"""
Itanium CXX ABI Mangler

Reference: https://itanium-cxx-abi.github.io/cxx-abi/abi.html

The basics of the mangling scheme.

We are hijacking the CXX mangling scheme for our use.  We map Python modules
into CXX namespace.  A `module1.submodule2.foo` is mapped to
`module1::submodule2::foo`.   For parameterized numba types, we treat them as
templated types; for example, `array(int64, 1d, C)` becomes an
`array<int64, 1, C>`.

All mangled names are prefixed with "_Z".  It is followed by the name of the
entity.  A name contains one or more identifiers.  Each identifier is encoded
as "<num of char><name>".   If the name is namespaced and, therefore,
has multiple identifiers, the entire name is encoded as "N<name>E".

For functions, arguments types follow.  There are condensed encodings for basic
built-in types; e.g. "i" for int, "f" for float.  For other types, the
previously mentioned name encoding should be used.

For templated types, the template parameters are encoded immediately after the
name.  If it is namespaced, it should be within the 'N' 'E' marker.  Template
parameters are encoded in "I<params>E", where each parameter is encoded using
the mentioned name encoding scheme.  Template parameters can contain literal
values like the '1' in the array type shown earlier.  There is special encoding
scheme for them to avoid leading digits.
"""


import re

from numba.core import types, config


# According the scheme, valid characters for mangled names are [a-zA-Z0-9_].
# We borrow the '_' as the escape character to encode invalid char into
# '_xx' where 'xx' is the hex codepoint.
_re_invalid_char = re.compile(r'[^a-z0-9_]', re.I)

PREFIX = "_Z"

# Numba types to mangled type code. These correspond with the codes listed in
# https://itanium-cxx-abi.github.io/cxx-abi/abi.html#mangling-builtin
if config.USE_LEGACY_TYPE_SYSTEM: # Old type system
    N2CODE = {
        types.void: 'v',
        types.boolean: 'b',
        types.uint8: 'h',
        types.int8: 'a',
        types.uint16: 't',
        types.int16: 's',
        types.uint32: 'j',
        types.int32: 'i',
        types.uint64: 'y',
        types.int64: 'x',
        types.float16: 'Dh',
        types.float32: 'f',
        types.float64: 'd'
    }
else:
    N2CODE = {
        types.void: 'v',
        types.py_bool: 'b',
        types.py_int: 'x',
        types.py_float: 'd',
        types.np_bool_: 'b',
        types.np_uint8: 'h',
        types.np_int8: 'a',
        types.np_uint16: 't',
        types.np_int16: 's',
        types.np_uint32: 'j',
        types.np_int32: 'i',
        types.np_uint64: 'y',
        types.np_int64: 'x',
        types.np_float16: 'Dh',
        types.np_float32: 'f',
        types.np_float64: 'd'
    }


def _escape_string(text):
    """Escape the given string so that it only contains ASCII characters
    of [a-zA-Z0-9_$].

    The dollar symbol ($) and other invalid characters are escaped into
    the string sequence of "$xx" where "xx" is the hex codepoint of the char.

    Multibyte characters are encoded into utf8 and converted into the above
    hex format.
    """

    def repl(m):
        return ''.join(('_%02x' % ch)
                       for ch in m.group(0).encode('utf8'))
    ret = re.sub(_re_invalid_char, repl, text)
    # Return str if we got a unicode (for py2)
    if not isinstance(ret, str):
        return ret.encode('ascii')
    return ret


def _fix_lead_digit(text):
    """
    Fix text with leading digit
    """
    if text and text[0].isdigit():
        return '_' + text
    else:
        return text


def _len_encoded(string):
    """
    Prefix string with digit indicating the length.
    Add underscore if string is prefixed with digits.
    """
    string = _fix_lead_digit(string)
    return '%u%s' % (len(string), string)


def mangle_abi_tag(abi_tag: str) -> str:
    return "B" + _len_encoded(_escape_string(abi_tag))


def mangle_identifier(ident, template_params='', *, abi_tags=(), uid=None):
    """
    Mangle the identifier with optional template parameters and abi_tags.

    Note:

    This treats '.' as '::' in C++.
    """
    if uid is not None:
        # Add uid to abi-tags
        abi_tags = (f"v{uid}", *abi_tags)
    parts = [_len_encoded(_escape_string(x)) for x in ident.split('.')]
    enc_abi_tags = list(map(mangle_abi_tag, abi_tags))
    extras = template_params + ''.join(enc_abi_tags)
    if len(parts) > 1:
        return 'N%s%sE' % (''.join(parts), extras)
    else:
        return '%s%s' % (parts[0], extras)


def mangle_type_or_value(typ):
    """
    Mangle type parameter and arbitrary value.
    """
    # Handle numba types
    if isinstance(typ, types.Type):
        if typ in N2CODE:
            return N2CODE[typ]
        else:
            return mangle_templated_ident(*typ.mangling_args)
    # Handle integer literal
    elif isinstance(typ, int):
        return 'Li%dE' % typ
    # Handle str as identifier
    elif isinstance(typ, str):
        return mangle_identifier(typ)
    # Otherwise
    else:
        enc = _escape_string(str(typ))
        return _len_encoded(enc)


# Alias
mangle_type = mangle_type_or_value
mangle_value = mangle_type_or_value


def mangle_templated_ident(identifier, parameters):
    """
    Mangle templated identifier.
    """
    template_params = ('I%sE' % ''.join(map(mangle_type_or_value, parameters))
                       if parameters else '')
    return mangle_identifier(identifier, template_params)


def mangle_args(argtys):
    """
    Mangle sequence of Numba type objects and arbitrary values.
    """
    return ''.join([mangle_type_or_value(t) for t in argtys])


def mangle(ident, argtys, *, abi_tags=(), uid=None):
    """
    Mangle identifier with Numba type objects and abi-tags.
    """
    return ''.join([PREFIX,
                    mangle_identifier(ident, abi_tags=abi_tags, uid=uid),
                    mangle_args(argtys)])


def prepend_namespace(mangled, ns):
    """
    Prepend namespace to mangled name.
    """
    if not mangled.startswith(PREFIX):
        raise ValueError('input is not a mangled name')
    elif mangled.startswith(PREFIX + 'N'):
        # nested
        remaining = mangled[3:]
        ret = PREFIX + 'N' + mangle_identifier(ns) + remaining
    else:
        # non-nested
        remaining = mangled[2:]
        head, tail = _split_mangled_ident(remaining)
        ret = PREFIX + 'N' + mangle_identifier(ns) + head + 'E' + tail
    return ret


def _split_mangled_ident(mangled):
    """
    Returns `(head, tail)` where `head` is the `<len> + <name>` encoded
    identifier and `tail` is the remaining.
    """
    ct = int(mangled)
    ctlen = len(str(ct))
    at = ctlen + ct
    return mangled[:at], mangled[at:]
