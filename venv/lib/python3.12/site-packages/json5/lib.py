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

import enum
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
    Type,
    Union,
)
import unicodedata

from json5.parser import Parser


# Used when encoding keys, below.
_reserved_word_re: Optional[re.Pattern] = None


class QuoteStyle(enum.Enum):
    """Controls how strings will be quoted during encoding.

    By default, for compatibility with the `json` module and older versions of
    `json5`, strings (not being used as keys and that are legal identifiers)
    will always be double-quoted, and any double quotes in the string will be
    escaped. This is `QuoteStyle.ALWAYS_DOUBLE`.  If you pass
    `QuoteStyle.ALWAYS_SINGLE`, then strings will always be single-quoted, and
    any single quotes in the string will be escaped.  If you pass
    `QuoteStyle.PREFER_DOUBLE`, then the behavior is the same as ALWAYS_DOUBLE
    and strings will be double-quoted *unless* the string contains more double
    quotes than single quotes, in which case the string will be single-quoted
    and single quotes will be escaped. If you pass `QuoteStyle.PREFER_SINGLE`,
    then the behavior is the same as ALWAYS_SINGLE and strings will be
    single-quoted *unless* the string contains more single quotes than double
    quotes, in which case the string will be double-quoted and any double
    quotes will be escaped.

    *Note:* PREFER_DOUBLE and PREFER_SINGLE can impact performance, since in
    order to know which encoding to use you have to iterate over the entire
    string to count the number of single and double quotes. The codes guesses
    at an encoding while doing so, but if it guess wrong, the entire string has
    to be re-encoded, which will slow things down. If you are very concerned
    about performance (a) you probably shouldn't be using this library in the
    first place, because it just isn't very fast, and (b) you should use
    ALWAYS_DOUBLE or ALWAYS_SINGLE, which won't have this issue.
    """

    ALWAYS_DOUBLE = 'always_double'
    ALWAYS_SINGLE = 'always_single'
    PREFER_DOUBLE = 'prefer_double'
    PREFER_SINGLE = 'prefer_single'


def load(
    fp: IO,
    *,
    encoding: Optional[str] = None,
    cls: Any = None,
    object_hook: Optional[Callable[[Mapping[str, Any]], Any]] = None,
    parse_float: Optional[Callable[[str], Any]] = None,
    parse_int: Optional[Callable[[str], Any]] = None,
    parse_constant: Optional[Callable[[str], Any]] = None,
    strict: bool = True,
    object_pairs_hook: Optional[
        Callable[[Iterable[Tuple[str, Any]]], Any]
    ] = None,
    allow_duplicate_keys: bool = True,
    consume_trailing: bool = True,
    start: Optional[int] = None,
) -> Any:
    """Deserialize ``fp`` (a ``.read()``-supporting file-like object
    containing a JSON document) to a Python object.

    Supports almost the same arguments as ``json.load()`` except that:
        - the `cls` keyword is ignored.
        - an extra `allow_duplicate_keys` parameter supports checking for
          duplicate keys in a object; by default, this is True for
          compatibility with ``json.load()``, but if set to False and
          the object contains duplicate keys, a ValueError will be raised.
        - an extra `consume_trailing` parameter specifies whether to
          consume any trailing characters after a valid object has been
          parsed. By default, this value is True and the only legal
          trailing characters are whitespace. If this value is set to False,
          parsing will stop when a valid object has been parsed and any
          trailing characters in the string will be ignored.
        - an extra `start` parameter specifies the zero-based offset into the
          file to start parsing at. If `start` is None, parsing will
          start at the current position in the file, and line number
          and column values will be reported as if starting from the
          beginning of the file; If `start` is not None,
          `load` will seek to zero and then read (and discard) the
          appropriate number of characters before beginning parsing;
          the file must be seekable for this to work correctly.

    You can use `load(..., consume_trailing=False)` to repeatedly read
    values from a file. However, in the current implementation `load` does
    this by reading the entire file into memory before doing anything, so
    it is not very efficient.

    Raises
        - `ValueError` if given an invalid document. This is different
          from the `json` module, which raises `json.JSONDecodeError`.
        - `UnicodeDecodeError` if given a byte string that is not a
          legal UTF-8 document (or the equivalent, if using a different
          `encoding`). This matches the `json` module.
    """

    s = fp.read()
    val, err, _ = parse(
        s,
        encoding=encoding,
        cls=cls,
        object_hook=object_hook,
        parse_float=parse_float,
        parse_int=parse_int,
        parse_constant=parse_constant,
        strict=strict,
        object_pairs_hook=object_pairs_hook,
        allow_duplicate_keys=allow_duplicate_keys,
        consume_trailing=consume_trailing,
        start=start,
    )
    if err:
        raise ValueError(err)
    return val


def loads(
    s: str,
    *,
    encoding: Optional[str] = None,
    cls: Any = None,
    object_hook: Optional[Callable[[Mapping[str, Any]], Any]] = None,
    parse_float: Optional[Callable[[str], Any]] = None,
    parse_int: Optional[Callable[[str], Any]] = None,
    parse_constant: Optional[Callable[[str], Any]] = None,
    strict: bool = True,
    object_pairs_hook: Optional[
        Callable[[Iterable[Tuple[str, Any]]], Any]
    ] = None,
    allow_duplicate_keys: bool = True,
    consume_trailing: bool = True,
    start: Optional[int] = None,
):
    """Deserialize ``s`` (a string containing a JSON5 document) to a Python
    object.

    Supports the same arguments as ``json.load()`` except that:
        - the `cls` keyword is ignored.
        - an extra `allow_duplicate_keys` parameter supports checking for
          duplicate keys in a object; by default, this is True for
          compatibility with ``json.load()``, but if set to False and
          the object contains duplicate keys, a ValueError will be raised.
        - an extra `consume_trailing` parameter specifies whether to
          consume any trailing characters after a valid object has been
          parsed. By default, this value is True and the only legal
          trailing characters are whitespace. If this value is set to False,
          parsing will stop when a valid object has been parsed and any
          trailing characters in the string will be ignored.
        - an extra `start` parameter specifies the zero-based offset into the
          string to start parsing at.

    Raises
        - `ValueError` if given an invalid document. This is different
          from the `json` module, which raises `json.JSONDecodeError`.
        - `UnicodeDecodeError` if given a byte string that is not a
          legal UTF-8 document (or the equivalent, if using a different
          `encoding`). This matches the `json` module.
    """

    val, err, _ = parse(
        s=s,
        encoding=encoding,
        cls=cls,
        object_hook=object_hook,
        parse_float=parse_float,
        parse_int=parse_int,
        parse_constant=parse_constant,
        strict=strict,
        object_pairs_hook=object_pairs_hook,
        allow_duplicate_keys=allow_duplicate_keys,
        consume_trailing=consume_trailing,
        start=start,
    )
    if err:
        raise ValueError(err)
    return val


def parse(
    s: str,
    *,
    encoding: Optional[str] = None,
    cls: Any = None,
    object_hook: Optional[Callable[[Mapping[str, Any]], Any]] = None,
    parse_float: Optional[Callable[[str], Any]] = None,
    parse_int: Optional[Callable[[str], Any]] = None,
    parse_constant: Optional[Callable[[str], Any]] = None,
    strict: bool = True,
    object_pairs_hook: Optional[
        Callable[[Iterable[Tuple[str, Any]]], Any]
    ] = None,
    allow_duplicate_keys: bool = True,
    consume_trailing: bool = True,
    start: Optional[int] = None,
):
    """Parse ```s``, returning positional information along with a value.

    This works exactly like `loads()`, except that (a) it returns the
    position in the string where the parsing stopped (either due to
    hitting an error or parsing a valid value) and any error as a string,
    (b) it takes an optional `consume_trailing` parameter that says whether
    to keep parsing the string after a valid value has been parsed; if True
    (the default), any trailing characters must be whitespace. If False,
    parsing stops when a valid value has been reached, (c) it takes an
    optional `start` parameter that specifies a zero-based offset to start
    parsing from in the string, and (d) the return value is different, as
    described below.

    `parse()` is useful if you have a string that might contain multiple
    values and you need to extract all of them; you can do so by repeatedly
    calling `parse`, setting `start` to the value returned in `position`
    from the previous call.

    Returns a tuple of (value, error_string, position). If the string
        was a legal value, `value` will be the deserialized value,
        `error_string` will be `None`, and `position` will be one
        past the zero-based offset where the parser stopped reading.
        If the string was not a legal value,
        `value` will be `None`, `error_string` will be the string value
        of the exception that would've been raised, and `position` will
        be the zero-based farthest offset into the string where the parser
        hit an error.

    Raises:
        - `UnicodeDecodeError` if given a byte string that is not a
          legal UTF-8 document (or the equivalent, if using a different
          `encoding`). This matches the `json` module.

    Note that this does *not* raise a `ValueError`; instead any error is
    returned as the second value in the tuple.

    You can use this method to read in a series of values from a string
    `s` as follows:

    >>> import json5
    >>> s = '1 2 3 4'
    >>> values = []
    >>> start = 0
    >>> while True:
    ...     v, err, pos = json5.parse(s, start=start, consume_trailing=False)
    ...     if v:
    ...         values.append(v)
    ...         start = pos
    ...         if start == len(s) or s[start:].isspace():
    ...             # Reached the end of the string (ignoring trailing
    ...             # whitespace
    ...             break
    ...         continue
    ...     raise ValueError(err)
    >>> values
    [1, 2, 3, 4]

    """
    assert cls is None, 'Custom decoders are not supported'

    if isinstance(s, bytes):
        encoding = encoding or 'utf-8'
        s = s.decode(encoding)

    if not s:
        raise ValueError('Empty strings are not legal JSON5')
    start = start or 0
    parser = Parser(s, '<string>', pos=start)
    ast, err, pos = parser.parse(
        global_vars={'_strict': strict, '_consume_trailing': consume_trailing}
    )
    if err:
        return None, err, pos

    try:
        value = _convert(
            ast,
            object_hook=object_hook,
            parse_float=parse_float,
            parse_int=parse_int,
            parse_constant=parse_constant,
            object_pairs_hook=object_pairs_hook,
            allow_duplicate_keys=allow_duplicate_keys,
        )
        return value, None, pos
    except ValueError as e:
        return None, str(e), pos


def _convert(
    ast,
    object_hook,
    parse_float,
    parse_int,
    parse_constant,
    object_pairs_hook,
    allow_duplicate_keys,
):
    def _fp_constant_parser(s):
        return float(s.replace('Infinity', 'inf').replace('NaN', 'nan'))

    def _dictify(pairs):
        if not allow_duplicate_keys:
            keys = set()
            for key, _ in pairs:
                if key in keys:
                    raise ValueError(f'Duplicate key "{key}" found in object')
                keys.add(key)

        if object_pairs_hook:
            return object_pairs_hook(pairs)
        if object_hook:
            return object_hook(dict(pairs))
        return dict(pairs)

    parse_float = parse_float or float
    parse_int = parse_int or int
    parse_constant = parse_constant or _fp_constant_parser

    return _walk_ast(ast, _dictify, parse_float, parse_int, parse_constant)


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
    cls: Optional[Type['JSON5Encoder']] = None,
    indent: Optional[Union[int, str]] = None,
    separators: Optional[Tuple[str, str]] = None,
    default: Optional[Callable[[Any], Any]] = None,
    sort_keys: bool = False,
    quote_keys: bool = False,
    trailing_commas: bool = True,
    allow_duplicate_keys: bool = True,
    quote_style: QuoteStyle = QuoteStyle.ALWAYS_DOUBLE,
    **kw,
):
    """Serialize ``obj`` to a JSON5-formatted stream to ``fp``,
    a ``.write()``-supporting file-like object.

    Supports the same arguments as ``dumps()``, below.

    Calling ``dump(obj, fp, quote_keys=True, trailing_commas=False, \
                   allow_duplicate_keys=True)``
    should produce exactly the same output as ``json.dump(obj, fp).``
    """

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
            quote_style=quote_style,
            **kw,
        )
    )


def dumps(
    obj: Any,
    *,
    skipkeys: bool = False,
    ensure_ascii: bool = True,
    check_circular: bool = True,
    allow_nan: bool = True,
    cls: Optional[Type['JSON5Encoder']] = None,
    indent: Optional[Union[int, str]] = None,
    separators: Optional[Tuple[str, str]] = None,
    default: Optional[Callable[[Any], Any]] = None,
    sort_keys: bool = False,
    quote_keys: bool = False,
    trailing_commas: bool = True,
    allow_duplicate_keys: bool = True,
    quote_style: QuoteStyle = QuoteStyle.ALWAYS_DOUBLE,
    **kw: Any,
):
    """Serialize ``obj`` to a JSON5-formatted string.

    Supports the same arguments as ``json.dumps()``, except that:

    - The ``encoding`` keyword is ignored; Unicode strings are always written.
    - By default, object keys that are legal identifiers are not quoted; if you
      pass ``quote_keys=True``, they will be.
    - By default, if lists and objects span multiple lines of output (i.e.,
      when ``indent`` >=0), the last item will have a trailing comma after it.
      If you pass ``trailing_commas=False``, it will not.
    - If you use a number, a boolean, or ``None`` as a key value in a dict, it
      will be converted to the corresponding JSON string value, e.g.  "1",
      "true", or "null". By default, ``dump()`` will match the `json` modules
      behavior and produce malformed JSON if you mix keys of different types
      that have the same converted value; e.g., ``{1: "foo", "1": "bar"}``
      produces '{"1": "foo", "1": "bar"}', an object with duplicated keys. If
      you pass ``allow_duplicate_keys=False``, an exception will be raised
      instead.
    - If `quote_keys` is true, then keys of objects will be enclosed in quotes,
      as in regular JSON. Otheriwse, keys will not be enclosed in quotes unless
      they contain whitespace.
    - If `trailing_commas` is false, then commas will not be inserted after the
      final elements of objects and arrays, as in regular JSON.  Otherwise,
      such commas will be inserted.
    - If `allow_duplicate_keys` is false, then only the last entry with a given
      key will be written. Otherwise, all entries with the same key will be
      written.
    - `quote_style` controls how strings are encoded. See the documentation
      for the `QuoteStyle` class, above, for how this is used.

      *Note*: Strings that are being used as unquoted keys are not affected
      by this parameter and remain unquoted.

      *`quote_style` was added in version 0.10.0*.

    Other keyword arguments are allowed and will be passed to the
    encoder so custom encoders can get them, but otherwise they will
    be ignored in an attempt to provide some amount of forward-compatibility.

    *Note:* the standard JSON module explicitly calls `int.__repr(obj)__`
    and `float.__repr(obj)__` to encode ints and floats, thereby bypassing
    any custom representations you might have for objects that are subclasses
    of ints and floats, and, for compatibility, JSON5 does the same thing.
    To override this behavior, create a subclass of JSON5Encoder
    that overrides `encode()` and handles your custom representation.

    For example:

    ```
    >>> import json5
    >>> from typing import Any, Set
    >>>
    >>> class Hex(int):
    ...    def __repr__(self):
    ...        return hex(self)
    >>>
    >>> class CustomEncoder(json5.JSON5Encoder):
    ...    def encode(
    ...        self, obj: Any, seen: Set, level: int, *, as_key: bool
    ...    ) -> str:
    ...        if isinstance(obj, Hex):
    ...            return repr(obj)
    ...        return super().encode(obj, seen, level, as_key=as_key)
    ...
    >>> json5.dumps([20, Hex(20)], cls=CustomEncoder)
    '[20, 0x14]'

    ```

    *Note:* calling ``dumps(obj, quote_keys=True, trailing_commas=False, \
                            allow_duplicate_keys=True)``
    should produce exactly the same output as ``json.dumps(obj).``
    """

    cls = cls or JSON5Encoder
    enc = cls(
        skipkeys=skipkeys,
        ensure_ascii=ensure_ascii,
        check_circular=check_circular,
        allow_nan=allow_nan,
        indent=indent,
        separators=separators,
        default=default,
        sort_keys=sort_keys,
        quote_keys=quote_keys,
        trailing_commas=trailing_commas,
        allow_duplicate_keys=allow_duplicate_keys,
        quote_style=quote_style,
        **kw,
    )
    return enc.encode(obj, seen=set(), level=0, as_key=False)


class JSON5Encoder:
    def __init__(
        self,
        *,
        skipkeys: bool = False,
        ensure_ascii: bool = True,
        check_circular: bool = True,
        allow_nan: bool = True,
        indent: Optional[Union[int, str]] = None,
        separators: Optional[Tuple[str, str]] = None,
        default: Optional[Callable[[Any], Any]] = None,
        sort_keys: bool = False,
        quote_keys: bool = False,
        trailing_commas: bool = True,
        allow_duplicate_keys: bool = True,
        quote_style: QuoteStyle = QuoteStyle.ALWAYS_DOUBLE,
        **kw,
    ):
        """Provides a class that may be overridden to customize the behavior
        of `dumps()`. The keyword args are the same as for that function.
        *Added in version 0.10.0"""
        # Ignore unrecognized keyword arguments in the hope of providing
        # some level of backwards- and forwards-compatibility.
        del kw

        self.skipkeys = skipkeys
        self.ensure_ascii = ensure_ascii
        self.check_circular = check_circular
        self.allow_nan = allow_nan
        self.indent = indent
        self.separators = separators
        if separators is None:
            separators = (', ', ': ') if indent is None else (',', ': ')
        self.item_separator, self.kv_separator = separators
        self.default_fn = default or _raise_type_error
        self.sort_keys = sort_keys
        self.quote_keys = quote_keys
        self.trailing_commas = trailing_commas
        self.allow_duplicate_keys = allow_duplicate_keys
        self.quote_style = quote_style

    def default(self, obj: Any) -> Any:
        """Provides a last-ditch option to encode a value that the encoder
        doesn't otherwise recognize, by converting `obj` to a value that
        *can* (and will) be serialized by the other methods in the class.

        Note: this must not return a serialized value (i.e., string)
        directly, as that'll result in a doubly-encoded value."""
        return self.default_fn(obj)

    def encode(
        self,
        obj: Any,
        seen: Set,
        level: int,
        *,
        as_key: bool,
    ) -> str:
        """Returns an JSON5-encoded version of an arbitrary object. This can
        be used to provide customized serialization of objects. Overridden
        methods of this class should handle their custom objects and then
        fall back to super.encode() if they've been passed a normal object.

        `seen` is used for duplicate object tracking when `check_circular`
        is True.

        `level` represents the current indentation level, which increases
        by one for each recursive invocation of encode (i.e., whenever
        we're encoding the values of a dict or a list).

        May raise `TypeError` if the object is the wrong type to be
        encoded (i.e., your custom routine can't handle it either), and
        `ValueError` if there's something wrong with the value, e.g.
        a float value of NaN when `allow_nan` is false.

        If `as_key` is true, the return value should be a double-quoted string
        representation of the object, unless obj is a string that can be an
        identifier (and quote_keys is false and obj isn't a reserved word).
        If the object should not be used as a key, `TypeError` should be
        raised; that allows the base implementation to implement `skipkeys`
        properly.
        """
        seen = seen or set()
        s = self._encode_basic_type(obj, as_key=as_key)
        if s is not None:
            return s

        if as_key:
            raise TypeError(f'Invalid key f{obj}')
        return self._encode_non_basic_type(obj, seen, level)

    def _encode_basic_type(self, obj: Any, *, as_key: bool) -> Optional[str]:
        """Returns None if the object is not a basic type."""

        if isinstance(obj, str):
            return self._encode_str(obj, as_key=as_key)

        # Check for True/False before ints because True and False are
        # also considered ints and so would be represented as 1 and 0
        # if we did ints first.
        if obj is True:
            return '"true"' if as_key else 'true'
        if obj is False:
            return '"false"' if as_key else 'false'
        if obj is None:
            return '"null"' if as_key else 'null'

        if isinstance(obj, int):
            return self._encode_int(obj, as_key=as_key)

        if isinstance(obj, float):
            return self._encode_float(obj, as_key=as_key)

        return None

    def _encode_int(self, obj: int, *, as_key: bool) -> str:
        s = int.__repr__(obj)
        return f'"{s}"' if as_key else s

    def _encode_float(self, obj: float, *, as_key: bool) -> str:
        if obj == float('inf'):
            allowed = self.allow_nan
            s = 'Infinity'
        elif obj == float('-inf'):
            allowed = self.allow_nan
            s = '-Infinity'
        elif math.isnan(obj):
            allowed = self.allow_nan
            s = 'NaN'
        else:
            allowed = True
            s = float.__repr__(obj)

        if not allowed:
            raise ValueError('Illegal JSON5 value: f{obj}')
        return f'"{s}"' if as_key else s

    def _encode_str(self, obj: str, *, as_key: bool) -> str:
        if (
            as_key
            and self.is_identifier(obj)
            and not self.quote_keys
            and not self.is_reserved_word(obj)
        ):
            return obj

        return self._encode_quoted_str(obj, self.quote_style)

    def _encode_quoted_str(self, obj: str, quote_style: QuoteStyle) -> str:
        """Returns a quoted string with a minimal number of escaped quotes."""
        ret = []
        double_quotes_seen = 0
        single_quotes_seen = 0
        sq = "'"
        dq = '"'
        for ch in obj:
            if ch == dq:
                # At first we will guess at which quotes to escape. If
                # we guess wrong, we reencode the string below.
                double_quotes_seen += 1
                if quote_style in (
                    QuoteStyle.ALWAYS_DOUBLE,
                    QuoteStyle.PREFER_DOUBLE,
                ):
                    encoded_ch = self._escape_ch(dq)
                else:
                    encoded_ch = dq
            elif ch == sq:
                single_quotes_seen += 1
                if quote_style in (
                    QuoteStyle.ALWAYS_SINGLE,
                    QuoteStyle.PREFER_SINGLE,
                ):
                    encoded_ch = self._escape_ch(sq)
                else:
                    encoded_ch = sq
            elif ch == '\\':
                encoded_ch = self._escape_ch(ch)
            else:
                o = ord(ch)
                if o < 32:
                    encoded_ch = self._escape_ch(ch)
                elif o < 128:
                    encoded_ch = ch
                elif not self.ensure_ascii and ch not in ('\u2028', '\u2029'):
                    encoded_ch = ch
                else:
                    encoded_ch = self._escape_ch(ch)
            ret.append(encoded_ch)

        # We may have guessed wrong and need to reencode the string.
        if (
            double_quotes_seen > single_quotes_seen
            and quote_style == QuoteStyle.PREFER_DOUBLE
        ):
            return self._encode_quoted_str(obj, QuoteStyle.ALWAYS_SINGLE)
        if (
            single_quotes_seen > double_quotes_seen
            and quote_style == QuoteStyle.PREFER_SINGLE
        ):
            return self._encode_quoted_str(obj, QuoteStyle.ALWAYS_DOUBLE)

        if quote_style in (QuoteStyle.ALWAYS_DOUBLE, QuoteStyle.PREFER_DOUBLE):
            return '"' + ''.join(ret) + '"'
        return "'" + ''.join(ret) + "'"

    def _escape_ch(self, ch: str) -> str:
        """Returns the backslash-escaped representation of the char."""
        if ch == '\\':
            return '\\\\'
        if ch == "'":
            return r'\''
        if ch == '"':
            return r'\"'
        if ch == '\n':
            return r'\n'
        if ch == '\r':
            return r'\r'
        if ch == '\t':
            return r'\t'
        if ch == '\b':
            return r'\b'
        if ch == '\f':
            return r'\f'
        if ch == '\v':
            return r'\v'
        if ch == '\0':
            return r'\0'

        o = ord(ch)
        if o < 65536:
            return rf'\u{o:04x}'

        val = o - 0x10000
        high = 0xD800 + (val >> 10)
        low = 0xDC00 + (val & 0x3FF)
        return rf'\u{high:04x}\u{low:04x}'

    def _encode_non_basic_type(self, obj, seen: Set, level: int) -> str:
        # Basic types can't be recursive so we only check for circularity
        # on non-basic types. If for some reason the caller was using a
        # subclass of a basic type and wanted to check circularity on it,
        # it'd have to do so directly in a subclass of JSON5Encoder.
        if self.check_circular:
            i = id(obj)
            if i in seen:
                raise ValueError('Circular reference detected.')
            seen.add(i)

        # Ideally we'd use collections.abc.Mapping and collections.abc.Sequence
        # here, but for backwards-compatibility with potential old callers,
        # we only check for the two attributes we need in each case.
        if hasattr(obj, 'keys') and hasattr(obj, '__getitem__'):
            s = self._encode_dict(obj, seen, level + 1)
        elif hasattr(obj, '__getitem__') and hasattr(obj, '__iter__'):
            s = self._encode_array(obj, seen, level + 1)
        else:
            s = self.encode(self.default(obj), seen, level, as_key=False)
            assert s is not None

        if self.check_circular:
            seen.remove(i)
        return s

    def _encode_dict(self, obj: Any, seen: set, level: int) -> str:
        if not obj:
            return '{}'

        indent_str, end_str = self._spacers(level)
        item_sep = self.item_separator + indent_str
        kv_sep = self.kv_separator

        if self.sort_keys:
            keys = sorted(obj.keys())
        else:
            keys = obj.keys()

        s = '{' + indent_str

        first_key = True
        new_keys = set()
        for key in keys:
            try:
                key_str = self.encode(key, seen, level, as_key=True)
            except TypeError:
                if self.skipkeys:
                    continue
                raise

            if not self.allow_duplicate_keys:
                if key_str in new_keys:
                    raise ValueError(f'duplicate key {repr(key)}')
                new_keys.add(key_str)

            if first_key:
                first_key = False
            else:
                s += item_sep

            val_str = self.encode(obj[key], seen, level, as_key=False)
            s += key_str + kv_sep + val_str

        s += end_str + '}'
        return s

    def _encode_array(self, obj: Any, seen: Set, level: int) -> str:
        if not obj:
            return '[]'

        indent_str, end_str = self._spacers(level)
        item_sep = self.item_separator + indent_str
        return (
            '['
            + indent_str
            + item_sep.join(
                self.encode(el, seen, level, as_key=False) for el in obj
            )
            + end_str
            + ']'
        )

    def _spacers(self, level: int) -> Tuple[str, str]:
        if self.indent is not None:
            end_str = ''
            if self.trailing_commas:
                end_str = ','
            if isinstance(self.indent, int):
                if self.indent > 0:
                    indent_str = '\n' + ' ' * self.indent * level
                    end_str += '\n' + ' ' * self.indent * (level - 1)
                else:
                    indent_str = '\n'
                    end_str += '\n'
            else:
                indent_str = '\n' + self.indent * level
                end_str += '\n' + self.indent * (level - 1)
        else:
            indent_str = ''
            end_str = ''
        return indent_str, end_str

    def is_identifier(self, key: str) -> bool:
        """Returns whether the string could be used as a legal
        EcmaScript/JavaScript identifier.

        There should normally be no reason to override this, unless
        the definition of identifiers change in later versions of the
        JSON5 spec and this implementation hasn't been updated to handle
        the changes yet."""
        if (
            not key
            or not self._is_id_start(key[0])
            and key[0] not in ('$', '_')
        ):
            return False
        for ch in key[1:]:
            if not self._is_id_continue(ch) and ch not in ('$', '_'):
                return False
        return True

    def _is_id_start(self, ch: str) -> bool:
        return unicodedata.category(ch) in (
            'Lu',
            'Ll',
            'Li',
            'Lt',
            'Lm',
            'Lo',
            'Nl',
        )

    def _is_id_continue(self, ch: str) -> bool:
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

    def is_reserved_word(self, key: str) -> bool:
        """Returns whether the key is a reserved word.

        There should normally be no need to override this, unless there
        have been reserved words added in later versions of the JSON5
        spec and this implementation has not yet been updated to handle
        the changes yet."""
        global _reserved_word_re
        if _reserved_word_re is None:
            # List taken from section 7.6.1 of ECMA-262, version 5.1.
            # https://262.ecma-international.org/5.1/#sec-7.6.1.
            # This includes currently reserved words, words reserved
            # for future use (both as of 5.1), null, true, and false.
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
                        'implements',
                        'import',
                        'in',
                        'instanceof',
                        'interface',
                        'let',
                        'new',
                        'null',
                        'package',
                        'private',
                        'protected',
                        'public',
                        'return',
                        'static',
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
                        'yield',
                    ]
                )
                + ')$'
            )
        return _reserved_word_re.match(key) is not None


def _raise_type_error(obj) -> Any:
    raise TypeError(f'{repr(obj)} is not JSON5 serializable')
