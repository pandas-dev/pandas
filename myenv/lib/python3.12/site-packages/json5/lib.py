# Copyright 2015 Google Inc. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import math
import re
from typing import (
    Any,
    Callable,
    IO,
    Iterable,
    Mapping,
    Optional,
    Set,
    Tuple,
    Union,
)
import unicodedata

from .parser import Parser


def load(
    fp: IO,
    *,
    encoding: Optional[str] = None,
    cls: None = None,
    object_hook: Optional[Callable[[Mapping[str, Any]], Any]] = None,
    parse_float: Optional[Callable[[str], Any]] = None,
    parse_int: Optional[Callable[[str], Any]] = None,
    parse_constant: Optional[Callable[[str], Any]] = None,
    object_pairs_hook: Optional[
        Callable[[Iterable[Tuple[str, Any]]], Any]
    ] = None,
    allow_duplicate_keys: bool = True,
) -> Any:
    """Deserialize ``fp`` (a ``.read()``-supporting file-like object
    containing a JSON document) to a Python object.

    Supports almost the same arguments as ``json.load()`` except that:
        - the `cls` keyword is ignored.
        - an extra `allow_duplicate_keys` parameter supports checking for
          duplicate keys in a object; by default, this is True for
          compatibility with ``json.load()``, but if set to False and
          the object contains duplicate keys, a ValueError will be raised.
    """

    s = fp.read()
    return loads(
        s,
        encoding=encoding,
        cls=cls,
        object_hook=object_hook,
        parse_float=parse_float,
        parse_int=parse_int,
        parse_constant=parse_constant,
        object_pairs_hook=object_pairs_hook,
        allow_duplicate_keys=allow_duplicate_keys,
    )


def loads(
    s: str,
    *,
    encoding: Optional[str] = None,
    cls: None = None,
    object_hook: Optional[Callable[[Mapping[str, Any]], Any]] = None,
    parse_float: Optional[Callable[[str], Any]] = None,
    parse_int: Optional[Callable[[str], Any]] = None,
    parse_constant: Optional[Callable[[str], Any]] = None,
    object_pairs_hook: Optional[
        Callable[[Iterable[Tuple[str, Any]]], Any]
    ] = None,
    allow_duplicate_keys: bool = True,
):
    """Deserialize ``s`` (a string containing a JSON5 document) to a Python
    object.

    Supports the same arguments as ``json.load()`` except that:
        - the `cls` keyword is ignored.
        - an extra `allow_duplicate_keys` parameter supports checking for
          duplicate keys in a object; by default, this is True for
          compatibility with ``json.load()``, but if set to False and
          the object contains duplicate keys, a ValueError will be raised.
    """

    assert cls is None, 'Custom decoders are not supported'

    if isinstance(s, bytes):
        encoding = encoding or 'utf-8'
        s = s.decode(encoding)

    if not s:
        raise ValueError('Empty strings are not legal JSON5')
    parser = Parser(s, '<string>')
    ast, err, _ = parser.parse()
    if err:
        raise ValueError(err)

    def _fp_constant_parser(s):
        return float(s.replace('Infinity', 'inf').replace('NaN', 'nan'))

    if object_pairs_hook:
        dictify = object_pairs_hook
    elif object_hook:

        def dictify(pairs):
            return object_hook(dict(pairs))
    else:
        dictify = dict

    if not allow_duplicate_keys:
        _orig_dictify = dictify

        def dictify(pairs):  # pylint: disable=function-redefined
            return _reject_duplicate_keys(pairs, _orig_dictify)

    parse_float = parse_float or float
    parse_int = parse_int or int
    parse_constant = parse_constant or _fp_constant_parser

    return _walk_ast(ast, dictify, parse_float, parse_int, parse_constant)


def _reject_duplicate_keys(pairs, dictify):
    keys = set()
    for key, _ in pairs:
        if key in keys:
            raise ValueError(f'Duplicate key "{key}" found in object')
        keys.add(key)
    return dictify(pairs)


def _walk_ast(
    el,
    dictify: Callable[[Iterable[Tuple[str, Any]]], Any],
    parse_float,
    parse_int,
    parse_constant,
):
    if el == 'None':
        return None
    if el == 'True':
        return True
    if el == 'False':
        return False
    ty, v = el
    if ty == 'number':
        if v.startswith('0x') or v.startswith('0X'):
            return parse_int(v, base=16)
        if '.' in v or 'e' in v or 'E' in v:
            return parse_float(v)
        if 'Infinity' in v or 'NaN' in v:
            return parse_constant(v)
        return parse_int(v)
    if ty == 'string':
        return v
    if ty == 'object':
        pairs = []
        for key, val_expr in v:
            val = _walk_ast(
                val_expr, dictify, parse_float, parse_int, parse_constant
            )
            pairs.append((key, val))
        return dictify(pairs)
    if ty == 'array':
        return [
            _walk_ast(el, dictify, parse_float, parse_int, parse_constant)
            for el in v
        ]
    raise ValueError('unknown el: ' + el)  # pragma: no cover


def dump(
    obj: Any,
    fp: IO,
    *,
    skipkeys: bool = False,
    ensure_ascii: bool = True,
    check_circular: bool = True,
    allow_nan: bool = True,
    cls: None = None,
    indent: Optional[Union[int, str]] = None,
    separators: Optional[Tuple[str, str]] = None,
    default: Optional[Callable[[Any], Any]] = None,
    sort_keys: bool = False,
    quote_keys: bool = False,
    trailing_commas: bool = True,
    allow_duplicate_keys: bool = True,
    **kwargs,
):
    """Serialize ``obj`` to a JSON5-formatted stream to ``fp``,
    a ``.write()``-supporting file-like object.

    Supports the same arguments as ``json.dump()``, except that:

    - The ``cls`` keyword is not supported.
    - The ``encoding`` keyword is ignored; Unicode strings are always
      written.
    - By default, object keys that are legal identifiers are not quoted;
      if you pass ``quote_keys=True``, they will be.
    - By default, if lists and objects span multiple lines of output (i.e.,
      when ``indent`` >=0), the last item will have a trailing comma
      after it. If you pass ``trailing_commas=False``, it will not.
    - If you use a number, a boolean, or ``None`` as a key value in a dict,
      it will be converted to the corresponding JSON string value, e.g.
      "1", "true", or "null". By default, ``dump()`` will match the `json`
      modules behavior and produce malformed JSON if you mix keys of
      different types that have the same converted value; e.g.,
      ``{1: "foo", "1": "bar"}`` produces '{"1": "foo", "1": "bar"}', an
      object with duplicated keys. If you pass
      ``allow_duplicate_keys=False``, an exception will be raised instead.
    - If `quote_keys` is true, then keys of objects will be enclosed in
      quotes, as in regular JSON. Otherwise, keys will not be enclosed in
      quotes unless they contain whitespace.
    - If `trailing_commas` is false, then commas will not be inserted after
      the final elements of objects and arrays, as in regular JSON.
      Otherwise, such commas will be inserted.
    - If `allow_duplicate_keys` is false, then only the last entry with a
      given key will be written. Otherwise, all entries with the same key
      will be written.

    Calling ``dump(obj, fp, quote_keys=True, trailing_commas=False, \
                   allow_duplicate_keys=True)``
    should produce exactly the same output as ``json.dump(obj, fp).``
    """

    del kwargs
    fp.write(
        dumps(
            obj=obj,
            skipkeys=skipkeys,
            ensure_ascii=ensure_ascii,
            check_circular=check_circular,
            allow_nan=allow_nan,
            cls=cls,
            indent=indent,
            separators=separators,
            default=default,
            sort_keys=sort_keys,
            quote_keys=quote_keys,
            trailing_commas=trailing_commas,
            allow_duplicate_keys=allow_duplicate_keys,
        )
    )


def dumps(
    obj: Any,
    *,
    skipkeys: bool = False,
    ensure_ascii: bool = True,
    check_circular: bool = True,
    allow_nan: bool = True,
    cls: None = None,
    indent: Optional[Union[int, str]] = None,
    separators: Optional[Tuple[str, str]] = None,
    default: Optional[Callable[[Any], Any]] = None,
    sort_keys: bool = False,
    quote_keys: bool = False,
    trailing_commas: bool = True,
    allow_duplicate_keys: bool = True,
    **kwargs,
):
    """Serialize ``obj`` to a JSON5-formatted string.

    Supports the same arguments as ``json.dumps()``, except that:

    - The ``cls`` keyword is not supported.
    - The ``encoding`` keyword is ignored; Unicode strings are always
      written.
    - By default, object keys that are legal identifiers are not quoted;
      if you pass ``quote_keys=True``, they will be.
    - By default, if lists and objects span multiple lines of output (i.e.,
      when ``indent`` >=0), the last item will have a trailing comma
      after it. If you pass ``trailing_commas=False``, it will not.
    - If you use a number, a boolean, or ``None`` as a key value in a dict,
      it will be converted to the corresponding JSON string value, e.g.
      "1", "true", or "null". By default, ``dump()`` will match the `json`
      modules behavior and produce malformed JSON if you mix keys of
      different types that have the same converted value; e.g.,
      ``{1: "foo", "1": "bar"}`` produces '{"1": "foo", "1": "bar"}', an
      object with duplicated keys. If you pass
      ``allow_duplicate_keys=False``, an exception will be raised instead.
    - If `quote_keys` is true, then keys of objects will be enclosed
      in quotes, as in regular JSON. Otheriwse, keys will not be enclosed
      in quotes unless they contain whitespace.
    - If `trailing_commas` is false, then commas will not be inserted after
      the final elements of objects and arrays, as in regular JSON.
      Otherwise, such commas will be inserted.
    - If `allow_duplicate_keys` is false, then only the last entry with a
      given key will be written. Otherwise, all entries with the same key
      will be written.

    Calling ``dumps(obj, quote_keys=True, trailing_commas=False, \
                    allow_duplicate_keys=True)``
    should produce exactly the same output as ``json.dumps(obj).``
    """

    assert kwargs.get('cls', None) is None, 'Custom encoders are not supported'
    del cls

    if separators is None:
        if indent is None:
            separators = (', ', ': ')
        else:
            separators = (',', ': ')

    default = default or _raise_type_error

    if check_circular:
        seen: Optional[Set[int]] = set()
    else:
        seen = None

    level = 1
    is_key = False

    _, v = _dumps(
        obj,
        skipkeys,
        ensure_ascii,
        check_circular,
        allow_nan,
        indent,
        separators,
        default,
        sort_keys,
        quote_keys,
        trailing_commas,
        allow_duplicate_keys,
        seen,
        level,
        is_key,
    )
    return v


def _dumps(
    obj,
    skipkeys,
    ensure_ascii,
    check_circular,
    allow_nan,
    indent,
    separators,
    default,
    sort_keys,
    quote_keys,
    trailing_commas,
    allow_duplicate_keys,
    seen: Optional[Set[int]],
    level: int,
    is_key: bool,
):
    # pylint: disable=too-many-statements
    if obj is True:
        s = 'true'
    elif obj is False:
        s = 'false'
    elif obj is None:
        s = 'null'
    elif obj == float('inf'):
        if allow_nan:
            s = 'Infinity'
        else:
            raise ValueError()
    elif obj == float('-inf'):
        if allow_nan:
            s = '-Infinity'
        else:
            raise ValueError()
    elif isinstance(obj, float) and math.isnan(obj):
        if allow_nan:
            s = 'NaN'
        else:
            raise ValueError()
    elif isinstance(obj, str):
        if (
            is_key
            and _is_ident(obj)
            and not quote_keys
            and not _is_reserved_word(obj)
        ):
            return True, obj
        return True, _dump_str(obj, ensure_ascii)
    elif isinstance(obj, int):
        # Subclasses of `int` and `float` may have custom
        # __repr__ or __str__ methods, but the `JSON` library
        # ignores them in order to ensure that the representation
        # are just bare numbers. In order to match JSON's behavior
        # we call the methods of the `float` and `int` class directly.
        s = int.__repr__(obj)
    elif isinstance(obj, float):
        # See comment above for int
        s = float.__repr__(obj)
    else:
        s = None

    if is_key:
        if s is not None:
            return True, f'"{s}"'
        if skipkeys:
            return False, None
        raise TypeError(f'invalid key {repr(obj)}')

    if s is not None:
        return True, s

    if indent is not None:
        end_str = ''
        if trailing_commas:
            end_str = ','
        if isinstance(indent, int):
            if indent > 0:
                indent_str = '\n' + ' ' * indent * level
                end_str += '\n' + ' ' * indent * (level - 1)
            else:
                indent_str = '\n'
                end_str += '\n'
        else:
            indent_str = '\n' + indent * level
            end_str += '\n' + indent * (level - 1)
    else:
        indent_str = ''
        end_str = ''

    item_sep, kv_sep = separators
    item_sep += indent_str

    if seen is not None:
        i = id(obj)
        if i in seen:
            raise ValueError('Circular reference detected.')
        seen.add(i)

    # Ideally we'd use collections.abc.Mapping and collections.abc.Sequence
    # here, but for backwards-compatibility with potential old callers,
    # we only check for the two attributes we need in each case.
    if hasattr(obj, 'keys') and hasattr(obj, '__getitem__'):
        s = _dump_dict(
            obj,
            skipkeys,
            ensure_ascii,
            check_circular,
            allow_nan,
            indent,
            separators,
            default,
            sort_keys,
            quote_keys,
            trailing_commas,
            allow_duplicate_keys,
            seen,
            level + 1,
            item_sep,
            kv_sep,
            indent_str,
            end_str,
        )
    elif hasattr(obj, '__getitem__') and hasattr(obj, '__iter__'):
        s = _dump_array(
            obj,
            skipkeys,
            ensure_ascii,
            check_circular,
            allow_nan,
            indent,
            separators,
            default,
            sort_keys,
            quote_keys,
            trailing_commas,
            allow_duplicate_keys,
            seen,
            level + 1,
            item_sep,
            indent_str,
            end_str,
        )
    else:
        s = _dumps(
            default(obj),
            skipkeys,
            ensure_ascii,
            check_circular,
            allow_nan,
            indent,
            separators,
            default,
            sort_keys,
            quote_keys,
            trailing_commas,
            allow_duplicate_keys,
            seen,
            level,
            is_key,
        )[1]

    if seen is not None:
        seen.remove(i)
    return False, s


def _dump_dict(
    obj,
    skipkeys,
    ensure_ascii,
    check_circular,
    allow_nan,
    indent,
    separators,
    default,
    sort_keys,
    quote_keys,
    trailing_commas,
    allow_duplicate_keys,
    seen,
    level,
    item_sep,
    kv_sep,
    indent_str,
    end_str,
):
    if not obj:
        return '{}'

    if sort_keys:
        keys = sorted(obj.keys())
    else:
        keys = obj.keys()

    s = '{' + indent_str

    num_items_added = 0
    new_keys = set()
    for key in keys:
        valid_key, key_str = _dumps(
            key,
            skipkeys,
            ensure_ascii,
            check_circular,
            allow_nan,
            indent,
            separators,
            default,
            sort_keys,
            quote_keys,
            trailing_commas,
            allow_duplicate_keys,
            seen,
            level,
            is_key=True,
        )

        if skipkeys and not valid_key:
            continue

        if not allow_duplicate_keys:
            if key_str in new_keys:
                raise ValueError(f'duplicate key {repr(key)}')
            new_keys.add(key_str)

        if num_items_added:
            s += item_sep

        s += (
            key_str
            + kv_sep
            + _dumps(
                obj[key],
                skipkeys,
                ensure_ascii,
                check_circular,
                allow_nan,
                indent,
                separators,
                default,
                sort_keys,
                quote_keys,
                trailing_commas,
                allow_duplicate_keys,
                seen,
                level,
                is_key=False,
            )[1]
        )
        num_items_added += 1

    s += end_str + '}'
    return s


def _dump_array(
    obj,
    skipkeys,
    ensure_ascii,
    check_circular,
    allow_nan,
    indent,
    separators,
    default,
    sort_keys,
    quote_keys,
    trailing_commas,
    allow_duplicate_keys,
    seen,
    level,
    item_sep,
    indent_str,
    end_str,
):
    if not obj:
        return '[]'
    return (
        '['
        + indent_str
        + item_sep.join(
            [
                _dumps(
                    el,
                    skipkeys,
                    ensure_ascii,
                    check_circular,
                    allow_nan,
                    indent,
                    separators,
                    default,
                    sort_keys,
                    quote_keys,
                    trailing_commas,
                    allow_duplicate_keys,
                    seen,
                    level,
                    False,
                )[1]
                for el in obj
            ]
        )
        + end_str
        + ']'
    )


def _dump_str(obj, ensure_ascii):
    ret = ['"']
    for ch in obj:
        if ch == '\\':
            ret.append('\\\\')
        elif ch == '"':
            ret.append('\\"')
        elif ch == '\u2028':
            ret.append('\\u2028')
        elif ch == '\u2029':
            ret.append('\\u2029')
        elif ch == '\n':
            ret.append('\\n')
        elif ch == '\r':
            ret.append('\\r')
        elif ch == '\b':
            ret.append('\\b')
        elif ch == '\f':
            ret.append('\\f')
        elif ch == '\t':
            ret.append('\\t')
        elif ch == '\v':
            ret.append('\\v')
        elif ch == '\0':
            ret.append('\\0')
        elif not ensure_ascii:
            ret.append(ch)
        else:
            o = ord(ch)
            if 32 <= o < 128:
                ret.append(ch)
            elif o < 65536:
                ret.append(f'\\u{o:04x}')
            else:
                val = o - 0x10000
                high = 0xD800 + (val >> 10)
                low = 0xDC00 + (val & 0x3FF)
                ret.append(f'\\u{high:04x}\\u{low:04x}')
    return ''.join(ret) + '"'


def _is_ident(k):
    if not k or not _is_id_start(k[0]) and k[0] not in ('$', '_'):
        return False
    for ch in k[1:]:
        if not _is_id_continue(ch) and ch not in ('$', '_'):
            return False
    return True


def _is_id_start(ch):
    return unicodedata.category(ch) in (
        'Lu',
        'Ll',
        'Li',
        'Lt',
        'Lm',
        'Lo',
        'Nl',
    )


def _is_id_continue(ch):
    return unicodedata.category(ch) in (
        'Lu',
        'Ll',
        'Li',
        'Lt',
        'Lm',
        'Lo',
        'Nl',
        'Nd',
        'Mn',
        'Mc',
        'Pc',
    )


_reserved_word_re = None


def _is_reserved_word(k):
    global _reserved_word_re

    if _reserved_word_re is None:
        # List taken from section 7.6.1 of ECMA-262.
        _reserved_word_re = re.compile(
            '('
            + '|'.join(
                [
                    'break',
                    'case',
                    'catch',
                    'class',
                    'const',
                    'continue',
                    'debugger',
                    'default',
                    'delete',
                    'do',
                    'else',
                    'enum',
                    'export',
                    'extends',
                    'false',
                    'finally',
                    'for',
                    'function',
                    'if',
                    'import',
                    'in',
                    'instanceof',
                    'new',
                    'null',
                    'return',
                    'super',
                    'switch',
                    'this',
                    'throw',
                    'true',
                    'try',
                    'typeof',
                    'var',
                    'void',
                    'while',
                    'with',
                ]
            )
            + ')$'
        )
    return _reserved_word_re.match(k) is not None


def _raise_type_error(obj):
    raise TypeError(f'{repr(obj)} is not JSON5 serializable')
