cimport cython
from cython cimport Py_ssize_t

import numpy as np

from cpython cimport (
    PyAnySet_Check,
    PyBool_Check,
    PyBytes_Check,
    PyBytes_GET_SIZE,
    PyDict_Check,
    PyFloat_CheckExact,
    PyList_Check,
    PyList_GET_ITEM,
    PyLong_CheckExact,
    PyTuple_Check,
    PyUnicode_CheckExact,
    PyUnicode_FromStringAndSize,
    PyUnicode_GET_LENGTH,
)
from cpython.bytes cimport (
    PyBytes_AS_STRING,
    PyBytes_FromStringAndSize,
    _PyBytes_Resize,
)
from cpython.mem cimport (
    PyMem_Free,
    PyMem_Malloc,
)
from cpython.ref cimport (
    Py_DECREF,
    Py_INCREF,
    PyObject,
)
from libc.math cimport isfinite
from libc.stdint cimport INT64_MAX
from libc.stdio cimport snprintf
from libc.string cimport (
    memchr,
    memcmp,
    memcpy,
)
from numpy cimport (
    int64_t,
    ndarray,
    uint8_t,
)

from decimal import Decimal

from pandas._libs.missing cimport C_NA
from pandas._libs.tslibs.nattype cimport c_NaT
from pandas._libs.tslibs.np_datetime cimport (
    NPY_DATETIME_NAT,
    NPY_DATETIMEUNIT,
    NPY_FR_ms,
    NPY_FR_ns,
    NPY_FR_s,
    NPY_FR_us,
    import_pandas_datetime,
    npy_datetimestruct,
    pandas_datetime_to_datetimestruct,
    pandas_timedelta_to_timedeltastruct,
    pandas_timedeltastruct,
)
from pandas._libs.tslibs.timedeltas cimport _Timedelta
from pandas._libs.tslibs.timestamps cimport _Timestamp

import_pandas_datetime()

ctypedef fused pandas_string:
    str
    bytes


@cython.boundscheck(False)
@cython.wraparound(False)
def write_csv_rows(
    list data,
    ndarray data_index,
    Py_ssize_t nlevels,
    ndarray cols,
    object writer
) -> None:
    """
    Write the given data to the writer object, pre-allocating where possible
    for performance improvements.

    Parameters
    ----------
    data : list[ArrayLike]
    data_index : ndarray
    nlevels : int
    cols : ndarray
    writer : _csv.writer
    """
    # In crude testing, N>100 yields little marginal improvement
    cdef:
        Py_ssize_t i, j = 0, k = len(data_index), N = 100, ncols = len(cols)
        list rows

    # pre-allocate rows
    rows = [[None] * (nlevels + ncols) for _ in range(N)]

    if nlevels == 1:
        for j in range(k):
            row = rows[j % N]
            row[0] = data_index[j]
            for i in range(ncols):
                row[1 + i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)
    elif nlevels > 1:
        for j in range(k):
            row = rows[j % N]
            row[:nlevels] = list(data_index[j])
            for i in range(ncols):
                row[nlevels + i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)
    else:
        for j in range(k):
            row = rows[j % N]
            for i in range(ncols):
                row[i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)

    if j >= 0 and (j < N - 1 or (j % N) != N - 1):
        writer.writerows(rows[:((j + 1) % N)])


@cython.boundscheck(False)
@cython.wraparound(False)
def convert_json_to_lines(arr: str) -> str:
    """
    replace comma separated json with line feeds, paying special attention
    to quotes & brackets
    """
    cdef:
        Py_ssize_t i = 0, num_open_brackets_seen = 0, length
        bint in_quotes = False, is_escaping = False
        ndarray[uint8_t, ndim=1] narr
        unsigned char val, newline, comma, left_bracket, right_bracket, quote
        unsigned char backslash

    newline = ord("\n")
    comma = ord(",")
    left_bracket = ord("{")
    right_bracket = ord("}")
    quote = ord('"')
    backslash = ord("\\")

    narr = np.frombuffer(arr.encode("utf-8"), dtype="u1").copy()
    length = narr.shape[0]
    for i in range(length):
        val = narr[i]
        if val == quote and i > 0 and not is_escaping:
            in_quotes = ~in_quotes
        if val == backslash or is_escaping:
            is_escaping = ~is_escaping
        if val == comma:  # commas that should be \n
            if num_open_brackets_seen == 0 and not in_quotes:
                narr[i] = newline
        elif val == left_bracket:
            if not in_quotes:
                num_open_brackets_seen += 1
        elif val == right_bracket:
            if not in_quotes:
                num_open_brackets_seen -= 1

    return narr.tobytes().decode("utf-8") + "\n"  # GH:36888


@cython.boundscheck(False)
@cython.wraparound(False)
def build_json_rows(list columns, Py_ssize_t nrows, list keys=None) -> list:
    """
    Transpose per-column lists into per-row containers for ``to_json``.

    Parameters
    ----------
    columns : list of list
        One list per column, each of length ``nrows``.
    nrows : int
    keys : list of str, optional
        Column labels. When given, each row is a dict mapping label to value
        (``orient="records"``); otherwise each row is a list
        (``orient="values"``).

    Returns
    -------
    list
    """
    cdef:
        Py_ssize_t i, j, ncols = len(columns)
        list out, row, col
        dict record

    for col in columns:
        if len(col) != nrows:
            raise ValueError("all columns must have length nrows")
    if keys is not None and len(keys) != ncols:
        raise ValueError("keys must have one entry per column")

    out = [None] * nrows
    if keys is None:
        for i in range(nrows):
            row = [None] * ncols
            for j in range(ncols):
                col = <list>PyList_GET_ITEM(columns, j)
                row[j] = <object>PyList_GET_ITEM(col, i)
            out[i] = row
    else:
        for i in range(nrows):
            record = {}
            for j in range(ncols):
                col = <list>PyList_GET_ITEM(columns, j)
                record[<object>PyList_GET_ITEM(keys, j)] = (
                    <object>PyList_GET_ITEM(col, i)
                )
            out[i] = record
    return out


# ----------------------------------------------------------------------
# datetime-like formatting for to_json

cdef NPY_DATETIMEUNIT _json_unit(str unit) except? NPY_FR_ns:
    if unit == "s":
        return NPY_FR_s
    if unit == "ms":
        return NPY_FR_ms
    if unit == "us":
        return NPY_FR_us
    if unit == "ns":
        return NPY_FR_ns
    raise ValueError(f"Invalid value '{unit}' for option 'date_unit'")


cdef inline int64_t _unit_per_second(NPY_DATETIMEUNIT unit) noexcept nogil:
    if unit == NPY_FR_s:
        return 1
    if unit == NPY_FR_ms:
        return 1_000
    if unit == NPY_FR_us:
        return 1_000_000
    return 1_000_000_000


cdef int _scale_epoch(
    int64_t value, NPY_DATETIMEUNIT from_unit, NPY_DATETIMEUNIT to_unit,
    int64_t* out
) noexcept nogil:
    """
    Rescale an epoch value between s/ms/us/ns, rounding down. Returns -1 on
    overflow, 0 otherwise.
    """
    cdef int64_t a = _unit_per_second(from_unit), b = _unit_per_second(to_unit)
    if a == b:
        out[0] = value
    elif a > b:
        out[0] = value // (a // b)
    else:
        if value > INT64_MAX // (b // a) or value < -(INT64_MAX // (b // a)):
            return -1
        out[0] = value * (b // a)
    return 0


cdef inline int _write_null(char* buf, bint as_string) noexcept nogil:
    if as_string:
        memcpy(buf, b'"null"', 6)
        return 6
    memcpy(buf, b"null", 4)
    return 4


cdef inline void _write_digits(char* buf, int64_t value, int width) noexcept nogil:
    """Write ``value`` zero-padded to ``width`` digits."""
    cdef int i
    for i in range(width - 1, -1, -1):
        buf[i] = <char>(48 + value % 10)
        value //= 10


cdef int _write_iso_datetime(
    char* buf, int64_t value, NPY_DATETIMEUNIT value_unit,
    NPY_DATETIMEUNIT out_unit, bint utc
) noexcept nogil:
    """Write ``value`` as an ISO 8601 string to ``buf``; return its length."""
    cdef:
        npy_datetimestruct dts
        int n

    pandas_datetime_to_datetimestruct(value, value_unit, &dts)
    if 0 <= dts.year <= 9999:
        _write_digits(buf, dts.year, 4)
        n = 4
    else:
        n = snprintf(buf, 24, "%lld", <long long>dts.year)
    buf[n] = b"-"
    _write_digits(buf + n + 1, dts.month, 2)
    buf[n + 3] = b"-"
    _write_digits(buf + n + 4, dts.day, 2)
    buf[n + 6] = b"T"
    _write_digits(buf + n + 7, dts.hour, 2)
    buf[n + 9] = b":"
    _write_digits(buf + n + 10, dts.min, 2)
    buf[n + 12] = b":"
    _write_digits(buf + n + 13, dts.sec, 2)
    n += 15
    if out_unit == NPY_FR_ms:
        buf[n] = b"."
        _write_digits(buf + n + 1, dts.us // 1000, 3)
        n += 4
    elif out_unit == NPY_FR_us:
        buf[n] = b"."
        _write_digits(buf + n + 1, dts.us, 6)
        n += 7
    elif out_unit == NPY_FR_ns:
        buf[n] = b"."
        _write_digits(buf + n + 1, dts.us, 6)
        _write_digits(buf + n + 7, dts.ps // 1000, 3)
        n += 10
    if utc:
        buf[n] = b"Z"
        n += 1
    return n


cdef int _write_iso_duration(
    char* buf, int64_t value, NPY_DATETIMEUNIT unit
) noexcept nogil:
    """
    Write ``value`` as an ISO 8601 duration to ``buf``; return its length.

    The fractional seconds are written in groups of three digits down to the
    smallest non-zero unit.
    """
    cdef:
        pandas_timedeltastruct tds
        int n

    pandas_timedelta_to_timedeltastruct(value, unit, &tds)
    n = snprintf(
        buf, 48, "P%lldDT%dH%dM%d", <long long>tds.days, tds.hrs, tds.min, tds.sec
    )
    if tds.ns:
        n += snprintf(buf + n, 12, ".%03d%03d%03d", tds.ms, tds.us, tds.ns)
    elif tds.us:
        n += snprintf(buf + n, 12, ".%03d%03d", tds.ms, tds.us)
    elif tds.ms:
        n += snprintf(buf + n, 12, ".%03d", tds.ms)
    buf[n] = b"S"
    return n + 1


@cython.boundscheck(False)
@cython.wraparound(False)
def datetime64_to_json(
    const int64_t[:] values, str value_unit, str out_unit, bint utc,
    bint nat_as_string=False,
) -> bytes:
    """
    A JSON array of ISO 8601 strings for datetime64 ``values`` in
    ``value_unit``, written with the precision of ``out_unit``; NaT becomes
    ``null`` (or the string ``"null"`` with ``nat_as_string``, for labels).
    """
    cdef:
        NPY_DATETIMEUNIT vunit = _json_unit(value_unit)
        NPY_DATETIMEUNIT ounit = _json_unit(out_unit)
        Py_ssize_t i, n = values.shape[0], pos = 1
        PyObject* raw
        char* buf

    raw = _new_bytes(n * 48 + 2, &buf)
    buf[0] = b"["
    with nogil:
        for i in range(n):
            if i:
                buf[pos] = b","
                pos += 1
            if values[i] == NPY_DATETIME_NAT:
                pos += _write_null(buf + pos, nat_as_string)
                continue
            buf[pos] = b'"'
            pos += 1
            pos += _write_iso_datetime(buf + pos, values[i], vunit, ounit, utc)
            buf[pos] = b'"'
            pos += 1
    buf[pos] = b"]"
    return _finish_bytes(raw, pos + 1)


@cython.boundscheck(False)
@cython.wraparound(False)
def timedelta64_to_json(
    const int64_t[:] values, str unit, bint nat_as_string=False
) -> bytes:
    """
    A JSON array of ISO 8601 duration strings for timedelta64 ``values`` in
    ``unit``; NaT becomes ``null``.
    """
    cdef:
        NPY_DATETIMEUNIT vunit = _json_unit(unit)
        Py_ssize_t i, n = values.shape[0], pos = 1
        PyObject* raw
        char* buf

    raw = _new_bytes(n * 64 + 2, &buf)
    buf[0] = b"["
    with nogil:
        for i in range(n):
            if i:
                buf[pos] = b","
                pos += 1
            if values[i] == NPY_DATETIME_NAT:
                pos += _write_null(buf + pos, nat_as_string)
                continue
            buf[pos] = b'"'
            pos += 1
            pos += _write_iso_duration(buf + pos, values[i], vunit)
            buf[pos] = b'"'
            pos += 1
    buf[pos] = b"]"
    return _finish_bytes(raw, pos + 1)


cdef Py_ssize_t _scan_json(
    const unsigned char* s, Py_ssize_t n, int* sign, Py_ssize_t* nan_count,
    bint* needs_stdlib,
) noexcept nogil:
    """
    Scan JSON text ``s`` outside of its string literals.

    Returns the position of the first integer outside of
    ``[-2**63, 2**64 - 1]`` with ``sign`` 1 (too big) or -1 (too small), or
    -1. ``nan_count`` receives the number of ``NaN`` literals; ``needs_stdlib``
    is set when the text holds ``Infinity`` or a float that may overflow,
    which only the stdlib decoder reads.
    """
    cdef:
        Py_ssize_t i = 0, start, dstart, ndig, exp
        unsigned char c
        bint neg, is_float

    while i < n:
        c = s[i]
        if c == b'"':
            i = _skip_json_string(s, i, n)
            continue
        if c == b"N":
            if i + 3 <= n and memcmp(s + i, b"NaN", 3) == 0:
                nan_count[0] += 1
                i += 3
            else:
                i += 1
            continue
        if c == b"I":
            if i + 8 <= n and memcmp(s + i, b"Infinity", 8) == 0:
                needs_stdlib[0] = True
                i += 8
            else:
                i += 1
            continue
        if c == b"-" or (c >= b"0" and c <= b"9"):
            start = i
            neg = c == b"-"
            if neg:
                i += 1
            dstart = i
            while i < n and s[i] >= b"0" and s[i] <= b"9":
                i += 1
            ndig = i - dstart
            is_float = i < n and (s[i] == b"." or s[i] == b"e" or s[i] == b"E")
            if ndig >= 19 and not is_float:
                if neg:
                    if ndig > 19 or memcmp(s + dstart, b"9223372036854775808", 19) > 0:
                        sign[0] = -1
                        return start
                elif ndig > 20 or (
                    ndig == 20
                    and memcmp(s + dstart, b"18446744073709551615", 20) > 0
                ):
                    sign[0] = 1
                    return start
            if is_float:
                if s[i] == b".":
                    i += 1
                    while i < n and s[i] >= b"0" and s[i] <= b"9":
                        i += 1
                exp = 0
                if i < n and (s[i] == b"e" or s[i] == b"E"):
                    i += 1
                    if i < n and s[i] == b"+":
                        i += 1
                    elif i < n and s[i] == b"-":
                        # a negative exponent cannot overflow
                        i += 1
                        exp = -100000
                    while i < n and s[i] >= b"0" and s[i] <= b"9":
                        if exp < 1000:
                            exp = exp * 10 + (s[i] - 48)
                        i += 1
                if exp >= 0 and ndig + exp >= 309:
                    needs_stdlib[0] = True
            continue
        i += 1
    return -1


@cython.boundscheck(False)
@cython.wraparound(False)
def json_scan(
    const unsigned char[::1] text
) -> tuple[Py_ssize_t, int, Py_ssize_t, bool]:
    """
    ``(position, sign, nan_count, needs_stdlib)`` for JSON ``text``: the
    first integer outside of the 64-bit range (sign 1 for too big, -1 for
    too small, position -1 if none), the number of bare ``NaN`` literals, and
    whether the text needs the stdlib decoder (``Infinity`` or a float that
    may overflow).
    """
    cdef:
        int sign = 0
        Py_ssize_t pos, nan_count = 0
        bint needs_stdlib = False

    if text.shape[0] == 0:
        return -1, 0, 0, False
    with nogil:
        pos = _scan_json(&text[0], text.shape[0], &sign, &nan_count, &needs_stdlib)
    return pos, sign, nan_count, needs_stdlib


@cython.boundscheck(False)
@cython.wraparound(False)
def json_nan_to_null(const unsigned char[::1] text, Py_ssize_t nan_count) -> bytes:
    """
    A copy of JSON ``text`` with its ``nan_count`` bare ``NaN`` literals
    replaced by ``null``.
    """
    cdef:
        Py_ssize_t i = 0, pos = 0, n = text.shape[0], j
        const unsigned char* s
        PyObject* raw
        char* buf

    if n == 0:
        return b""
    s = &text[0]
    raw = _new_bytes(n + nan_count, &buf)
    with nogil:
        while i < n:
            if s[i] == b'"':
                j = _skip_json_string(s, i, n)
                memcpy(buf + pos, s + i, j - i)
                pos += j - i
                i = j
                continue
            if s[i] == b"N" and i + 3 <= n and memcmp(s + i, b"NaN", 3) == 0:
                memcpy(buf + pos, b"null", 4)
                pos += 4
                i += 3
                continue
            buf[pos] = s[i]
            pos += 1
            i += 1
    return _finish_bytes(raw, pos)


@cython.boundscheck(False)
@cython.wraparound(False)
def int64_to_json(const int64_t[:] values, bint nat_as_string=False) -> bytes:
    """A JSON array of integers; the NaT sentinel becomes ``null``."""
    cdef:
        Py_ssize_t i, n = values.shape[0], pos = 1
        PyObject* raw
        char* buf

    raw = _new_bytes(n * 22 + 2, &buf)
    buf[0] = b"["
    with nogil:
        for i in range(n):
            if i:
                buf[pos] = b","
                pos += 1
            if values[i] == NPY_DATETIME_NAT:
                pos += _write_null(buf + pos, nat_as_string)
                continue
            pos += snprintf(buf + pos, 22, "%lld", <long long>values[i])
    buf[pos] = b"]"
    return _finish_bytes(raw, pos + 1)


# ----------------------------------------------------------------------
# object normalization for to_json

cdef class _ObjectNormalizer:
    """
    Convert arbitrary Python objects to JSON-native ones, with fast paths for
    the types that commonly appear in object columns. Anything else goes
    through ``fallback`` (the Python normalizer in pandas.io.json._encode).
    """
    cdef:
        object fallback
        bint nan_to_none, iso_dates
        NPY_DATETIMEUNIT out_unit

    def __init__(self, fallback, bint nan_to_none, bint iso_dates, str date_unit):
        self.fallback = fallback
        self.nan_to_none = nan_to_none
        self.iso_dates = iso_dates
        self.out_unit = _json_unit(date_unit)

    cdef object normalize(self, object val):
        cdef:
            char buf[64]
            int n
            int64_t scaled
            object key

        if val is None or val is c_NaT or val is C_NA:
            return None
        if (
            PyUnicode_CheckExact(val)
            or PyLong_CheckExact(val)
            or PyBool_Check(val)
        ):
            return val
        if PyFloat_CheckExact(val):
            if self.nan_to_none and not isfinite(<double>val):
                return None
            return val
        if isinstance(val, _Timestamp):
            if self.iso_dates:
                n = _write_iso_datetime(
                    buf, (<_Timestamp>val)._value, (<_Timestamp>val)._creso,
                    self.out_unit, val.tzinfo is not None,
                )
                return PyUnicode_FromStringAndSize(buf, n)
            if _scale_epoch(
                (<_Timestamp>val)._value, (<_Timestamp>val)._creso,
                self.out_unit, &scaled
            ) == 0:
                return scaled
            return self.fallback(val)
        if isinstance(val, _Timedelta):
            if self.iso_dates:
                n = _write_iso_duration(
                    buf, (<_Timedelta>val)._value, (<_Timedelta>val)._creso
                )
                return PyUnicode_FromStringAndSize(buf, n)
            if _scale_epoch(
                (<_Timedelta>val)._value, (<_Timedelta>val)._creso,
                self.out_unit, &scaled
            ) == 0:
                return scaled
            return self.fallback(val)
        if PyDict_Check(val):
            if self.is_native_dict(<dict>val):
                return val
            return {
                self.normalize_key(key): self.normalize(item)
                for key, item in (<dict>val).items()
            }
        if PyList_Check(val):
            if self.is_native_sequence(val):
                return val
            return [self.normalize(item) for item in val]
        if PyTuple_Check(val) or PyAnySet_Check(val):
            return [self.normalize(item) for item in val]
        if type(val) is Decimal:
            if val.is_nan():
                return None
            return format(val, "f")
        return self.fallback(val)

    cdef bint is_native(self, object val) noexcept:
        """Whether ``val`` can be handed to the engine as is."""
        if val is None or PyUnicode_CheckExact(val) or PyLong_CheckExact(val):
            return True
        if PyBool_Check(val):
            return True
        if PyFloat_CheckExact(val):
            return not self.nan_to_none or isfinite(<double>val)
        if PyList_Check(val):
            return self.is_native_sequence(val)
        if PyDict_Check(val):
            return self.is_native_dict(<dict>val)
        return False

    cdef bint is_native_sequence(self, object seq) noexcept:
        for item in seq:
            if not self.is_native(item):
                return False
        return True

    cdef bint is_native_dict(self, dict d) noexcept:
        for key, item in d.items():
            if not PyUnicode_CheckExact(key) or not self.is_native(item):
                return False
        return True

    cdef object normalize_key(self, object key):
        if PyUnicode_CheckExact(key):
            return key
        if PyBytes_Check(key):
            return (<bytes>key).decode("utf-8")
        return str(key)


@cython.boundscheck(False)
@cython.wraparound(False)
def normalize_json_objects(
    ndarray[object] arr,
    object fallback,
    bint nan_to_none,
    bint iso_dates,
    str date_unit,
) -> list:
    """
    Convert an object array to a list of JSON-native values for ``to_json``.

    str, int, bool, None, float, Timestamp, Timedelta, Decimal, dict, list,
    tuple and set are handled here (recursively for containers); everything
    else goes through ``fallback``.
    """
    cdef:
        Py_ssize_t i, n = len(arr)
        list out = [None] * n
        _ObjectNormalizer normalizer = _ObjectNormalizer(
            fallback, nan_to_none, iso_dates, date_unit
        )

    for i in range(n):
        out[i] = normalizer.normalize(arr[i])
    return out


# ----------------------------------------------------------------------
# JSON text splicing for to_json
#
# Building a dict per row (orient="records" / "index") or a dict per column
# (orient="columns") is the dominant cost of to_json for large frames. When
# the engine can splice pre-encoded text (orjson.Fragment), labels and values
# are encoded as flat JSON arrays instead, and zipped into objects here by
# scanning the text.


cdef inline Py_ssize_t _skip_json_string(
    const unsigned char* s, Py_ssize_t i, Py_ssize_t end
) noexcept nogil:
    """``s[i]`` is an opening quote; return the index after the closing quote."""
    cdef:
        const unsigned char* q
        Py_ssize_t j, backslashes

    i += 1
    while i < end:
        q = <const unsigned char*>memchr(s + i, b'"', end - i)
        if q == NULL:
            return end
        j = q - s
        # a quote preceded by an odd number of backslashes is escaped
        backslashes = 0
        while j - 1 - backslashes >= i and s[j - 1 - backslashes] == b"\\":
            backslashes += 1
        if backslashes % 2 == 0:
            return j + 1
        i = j + 1
    return end


cdef inline bint _next_json_element(
    const unsigned char* s, Py_ssize_t* pos, Py_ssize_t end,
    Py_ssize_t* begin, Py_ssize_t* stop
) noexcept nogil:
    """
    Find the next top-level element of a JSON array body in ``s[*pos:end]``
    (the text between the outer brackets). Sets ``begin``/``stop`` to the
    element's span and advances ``pos`` past it; returns False when the body
    is exhausted.
    """
    cdef:
        Py_ssize_t i = pos[0], depth = 0
        const unsigned char* q
        unsigned char c

    while i < end:
        c = s[i]
        if c == b" " or c == b"\n" or c == b"\t" or c == b"\r" or c == b",":
            i += 1
        else:
            break
    if i >= end:
        return False
    begin[0] = i
    c = s[i]
    if c == b'"':
        i = _skip_json_string(s, i, end)
    elif c != b"[" and c != b"{":
        # number, true, false or null: runs until the next comma
        q = <const unsigned char*>memchr(s + i, b",", end - i)
        i = end if q == NULL else q - s
    else:
        while i < end:
            c = s[i]
            if c == b'"':
                i = _skip_json_string(s, i, end)
                continue
            if c == b"[" or c == b"{":
                depth += 1
            elif c == b"]" or c == b"}":
                depth -= 1
            elif c == b"," and depth == 0:
                break
            i += 1
    stop[0] = i
    pos[0] = i
    return True


cdef inline PyObject* _new_bytes(Py_ssize_t size, char** buf) except NULL:
    """A bytes object of ``size`` owned by the caller (refcount 1)."""
    cdef bytes tmp = PyBytes_FromStringAndSize(NULL, size)
    cdef PyObject* raw = <PyObject*>tmp
    Py_INCREF(tmp)
    buf[0] = PyBytes_AS_STRING(tmp)
    return raw


cdef inline bytes _finish_bytes(PyObject* raw, Py_ssize_t size):
    """Shrink the bytes object from ``_new_bytes`` to ``size`` and return it."""
    _PyBytes_Resize(&raw, size)
    out = <bytes>raw
    Py_DECREF(out)
    return out


cdef inline void _check_json_array(const unsigned char* s, Py_ssize_t n):
    if n < 2 or s[0] != b"[" or s[n - 1] != b"]":
        raise ValueError("expected a JSON array")


cdef inline Py_ssize_t _count_json_elements(
    const unsigned char* s, Py_ssize_t n
) noexcept nogil:
    cdef Py_ssize_t pos = 1, beg, end, count = 0
    while _next_json_element(s, &pos, n - 1, &beg, &end):
        count += 1
    return count


@cython.boundscheck(False)
@cython.wraparound(False)
def json_zip_object(bytes keys, bytes values) -> bytes:
    """
    Combine two JSON arrays of equal length into a JSON object.

    ``keys`` must be an array of strings: ``["a","b"]`` and ``[1,2]`` give
    ``{"a":1,"b":2}``.
    """
    cdef:
        const unsigned char* k = <const unsigned char*>PyBytes_AS_STRING(keys)
        const unsigned char* v = <const unsigned char*>PyBytes_AS_STRING(values)
        Py_ssize_t nk = len(keys), nv = len(values), pos
        Py_ssize_t kpos, vpos, kbeg, kend, vbeg, vend
        bint has_key, has_val, first = True
        PyObject* raw
        char* buf

    _check_json_array(k, nk)
    _check_json_array(v, nv)
    raw = _new_bytes(nk + nv + 2, &buf)
    buf[0] = b"{"
    pos = 1
    kpos = 1
    vpos = 1
    with nogil:
        while True:
            has_key = _next_json_element(k, &kpos, nk - 1, &kbeg, &kend)
            has_val = _next_json_element(v, &vpos, nv - 1, &vbeg, &vend)
            if not has_key or not has_val:
                break
            if not first:
                buf[pos] = b","
                pos += 1
            first = False
            memcpy(buf + pos, k + kbeg, kend - kbeg)
            pos += kend - kbeg
            buf[pos] = b":"
            pos += 1
            memcpy(buf + pos, v + vbeg, vend - vbeg)
            pos += vend - vbeg
    if has_key != has_val:
        Py_DECREF(<object>raw)
        raise ValueError("keys and values must have the same length")
    buf[pos] = b"}"
    pos += 1
    return _finish_bytes(raw, pos)


@cython.boundscheck(False)
@cython.wraparound(False)
def json_zip_records(bytes keys, bytes rows, bytes labels=None) -> bytes:
    """
    Turn a JSON array of arrays into an array of objects.

    ``keys`` is a JSON array of strings with one entry per inner array:
    ``["a","b"]`` and ``[[1,2],[3,4]]`` give ``[{"a":1,"b":2},{"a":3,"b":4}]``.

    With ``labels``, a JSON array of strings with one entry per row, the
    result is an object keyed by label instead:
    ``{"r":{"a":1,"b":2},"s":{"a":3,"b":4}}``.
    """
    cdef:
        const unsigned char* k = <const unsigned char*>PyBytes_AS_STRING(keys)
        const unsigned char* r = <const unsigned char*>PyBytes_AS_STRING(rows)
        const unsigned char* lab = NULL
        Py_ssize_t nk = len(keys), nr = len(rows), nlab = 0
        Py_ssize_t pos, j, ncols, nrows
        Py_ssize_t kpos, rpos, vpos, lpos = 1, beg, end, rbeg, rend, vbeg, vend
        Py_ssize_t* kbeg
        Py_ssize_t* kend
        bint bad_row = False, bad_labels = False
        PyObject* raw
        char* buf

    _check_json_array(k, nk)
    _check_json_array(r, nr)
    if labels is not None:
        lab = <const unsigned char*>PyBytes_AS_STRING(labels)
        nlab = len(labels)
        _check_json_array(lab, nlab)

    ncols = _count_json_elements(k, nk)
    nrows = _count_json_elements(r, nr)
    kbeg = <Py_ssize_t*>PyMem_Malloc(max(ncols, 1) * sizeof(Py_ssize_t))
    kend = <Py_ssize_t*>PyMem_Malloc(max(ncols, 1) * sizeof(Py_ssize_t))
    if kbeg == NULL or kend == NULL:
        PyMem_Free(kbeg)
        PyMem_Free(kend)
        raise MemoryError()
    kpos = 1
    j = 0
    while _next_json_element(k, &kpos, nk - 1, &beg, &end):
        kbeg[j] = beg
        kend[j] = end
        j += 1

    # every row repeats the keys, plus ":" per value; the row text itself
    # already holds the brackets and commas
    raw = _new_bytes(nr + nlab + nrows * (nk + ncols + 1) + 2, &buf)
    buf[0] = b"{" if labels is not None else b"["
    pos = 1
    rpos = 1
    with nogil:
        while _next_json_element(r, &rpos, nr - 1, &rbeg, &rend):
            if pos > 1:
                buf[pos] = b","
                pos += 1
            if lab != NULL:
                if not _next_json_element(lab, &lpos, nlab - 1, &beg, &end):
                    bad_labels = True
                    break
                memcpy(buf + pos, lab + beg, end - beg)
                pos += end - beg
                buf[pos] = b":"
                pos += 1
            if rend - rbeg < 2 or r[rbeg] != b"[" or r[rend - 1] != b"]":
                bad_row = True
                break
            buf[pos] = b"{"
            pos += 1
            vpos = rbeg + 1
            for j in range(ncols):
                if not _next_json_element(r, &vpos, rend - 1, &vbeg, &vend):
                    bad_row = True
                    break
                if j:
                    buf[pos] = b","
                    pos += 1
                memcpy(buf + pos, k + kbeg[j], kend[j] - kbeg[j])
                pos += kend[j] - kbeg[j]
                buf[pos] = b":"
                pos += 1
                memcpy(buf + pos, r + vbeg, vend - vbeg)
                pos += vend - vbeg
            if bad_row:
                break
            if _next_json_element(r, &vpos, rend - 1, &vbeg, &vend):
                bad_row = True
                break
            buf[pos] = b"}"
            pos += 1
        if lab != NULL and not bad_row and not bad_labels:
            if _next_json_element(lab, &lpos, nlab - 1, &beg, &end):
                bad_labels = True
    PyMem_Free(kbeg)
    PyMem_Free(kend)
    if bad_row or bad_labels:
        Py_DECREF(<object>raw)
        if bad_row:
            raise ValueError("every row must be an array with one value per key")
        raise ValueError("labels and rows must have the same length")
    buf[pos] = b"}" if labels is not None else b"]"
    pos += 1
    return _finish_bytes(raw, pos)


@cython.boundscheck(False)
@cython.wraparound(False)
def json_quote_elements(bytes values) -> bytes:
    """
    Wrap every element of a flat JSON array of numbers in double quotes:
    ``[1,2]`` gives ``["1","2"]``.
    """
    cdef:
        const unsigned char* v = <const unsigned char*>PyBytes_AS_STRING(values)
        Py_ssize_t nv = len(values), pos = 1, vpos = 1, vbeg, vend
        PyObject* raw
        char* buf

    _check_json_array(v, nv)
    raw = _new_bytes(3 * nv + 2, &buf)
    buf[0] = b"["
    with nogil:
        while _next_json_element(v, &vpos, nv - 1, &vbeg, &vend):
            if pos > 1:
                buf[pos] = b","
                pos += 1
            buf[pos] = b'"'
            pos += 1
            memcpy(buf + pos, v + vbeg, vend - vbeg)
            pos += vend - vbeg
            buf[pos] = b'"'
            pos += 1
    buf[pos] = b"]"
    pos += 1
    return _finish_bytes(raw, pos)


@cython.boundscheck(False)
@cython.wraparound(False)
def json_object_of_fragments(bytes keys, list values) -> bytes:
    """
    Build a JSON object from a JSON array of string keys and a list of
    already-encoded JSON values (bytes), one per key.
    """
    cdef:
        const unsigned char* k = <const unsigned char*>PyBytes_AS_STRING(keys)
        Py_ssize_t nk = len(keys), n = len(values), i, pos = 1, kpos = 1
        Py_ssize_t kbeg, kend, total = nk + 2 + n
        bytes value
        PyObject* raw
        char* buf

    _check_json_array(k, nk)
    for value in values:
        total += len(value)
    raw = _new_bytes(total, &buf)
    buf[0] = b"{"
    for i in range(n):
        if not _next_json_element(k, &kpos, nk - 1, &kbeg, &kend):
            Py_DECREF(<object>raw)
            raise ValueError("keys and values must have the same length")
        if i:
            buf[pos] = b","
            pos += 1
        memcpy(buf + pos, k + kbeg, kend - kbeg)
        pos += kend - kbeg
        buf[pos] = b":"
        pos += 1
        value = <bytes>PyList_GET_ITEM(values, i)
        memcpy(buf + pos, PyBytes_AS_STRING(value), len(value))
        pos += len(value)
    if _next_json_element(k, &kpos, nk - 1, &kbeg, &kend):
        Py_DECREF(<object>raw)
        raise ValueError("keys and values must have the same length")
    buf[pos] = b"}"
    pos += 1
    return _finish_bytes(raw, pos)


@cython.boundscheck(False)
@cython.wraparound(False)
def json_transpose(list columns, bytes keys=None, bytes labels=None) -> bytes:
    """
    Turn per-column JSON arrays into per-row output.

    ``columns`` is a list of bytes, each a JSON array with one element per
    row. Without ``keys`` the result is an array of arrays; with ``keys`` (a
    JSON array of strings, one per column) an array of objects; with
    ``labels`` as well (a JSON array of strings, one per row) an object keyed
    by label.
    """
    cdef:
        Py_ssize_t ncols = len(columns), i, j, pos, nrows, total, nk = 0, nlab = 0
        Py_ssize_t beg, end, lbeg, lend, lpos = 1, kpos = 1
        const unsigned char** cols
        Py_ssize_t* ncol
        Py_ssize_t* cpos
        Py_ssize_t* kbeg
        Py_ssize_t* kend
        const unsigned char* k = NULL
        const unsigned char* lab = NULL
        bytes col
        bint short_column = False
        PyObject* raw
        char* buf

    if keys is not None:
        k = <const unsigned char*>PyBytes_AS_STRING(keys)
        nk = len(keys)
        _check_json_array(k, nk)
        if _count_json_elements(k, nk) != ncols:
            raise ValueError("keys must have one entry per column")
    if labels is not None:
        if keys is None:
            raise ValueError("labels require keys")
        lab = <const unsigned char*>PyBytes_AS_STRING(labels)
        nlab = len(labels)
        _check_json_array(lab, nlab)

    cols = <const unsigned char**>PyMem_Malloc(max(ncols, 1) * sizeof(void*))
    ncol = <Py_ssize_t*>PyMem_Malloc(max(ncols, 1) * 4 * sizeof(Py_ssize_t))
    if cols == NULL or ncol == NULL:
        PyMem_Free(cols)
        PyMem_Free(ncol)
        raise MemoryError()
    cpos = ncol + ncols
    kbeg = cpos + ncols
    kend = kbeg + ncols
    total = 2
    for i in range(ncols):
        col = <bytes>PyList_GET_ITEM(columns, i)
        cols[i] = <const unsigned char*>PyBytes_AS_STRING(col)
        ncol[i] = len(col)
        _check_json_array(cols[i], ncol[i])
        cpos[i] = 1
        total += ncol[i]
    if keys is not None:
        for j in range(ncols):
            _next_json_element(k, &kpos, nk - 1, &kbeg[j], &kend[j])
    nrows = _count_json_elements(cols[0], ncol[0]) if ncols else 0
    if labels is not None and _count_json_elements(lab, nlab) != nrows:
        PyMem_Free(cols)
        PyMem_Free(ncol)
        raise ValueError("labels must have one entry per row")
    # per row: brackets/braces, separators, keys with ":" and the label
    total += nrows * (3 + ncols + nk + ncols) + nlab

    raw = _new_bytes(total, &buf)
    buf[0] = b"{" if labels is not None else b"["
    pos = 1
    with nogil:
        for i in range(nrows):
            if i:
                buf[pos] = b","
                pos += 1
            if lab != NULL:
                _next_json_element(lab, &lpos, nlab - 1, &lbeg, &lend)
                memcpy(buf + pos, lab + lbeg, lend - lbeg)
                pos += lend - lbeg
                buf[pos] = b":"
                pos += 1
            buf[pos] = b"{" if k != NULL else b"["
            pos += 1
            for j in range(ncols):
                if not _next_json_element(cols[j], &cpos[j], ncol[j] - 1, &beg, &end):
                    short_column = True
                    break
                if j:
                    buf[pos] = b","
                    pos += 1
                if k != NULL:
                    memcpy(buf + pos, k + kbeg[j], kend[j] - kbeg[j])
                    pos += kend[j] - kbeg[j]
                    buf[pos] = b":"
                    pos += 1
                memcpy(buf + pos, cols[j] + beg, end - beg)
                pos += end - beg
            if short_column:
                break
            buf[pos] = b"}" if k != NULL else b"]"
            pos += 1
        if not short_column:
            for j in range(ncols):
                if _next_json_element(cols[j], &cpos[j], ncol[j] - 1, &beg, &end):
                    short_column = True
                    break
    PyMem_Free(cols)
    PyMem_Free(ncol)
    if short_column:
        Py_DECREF(<object>raw)
        raise ValueError("all columns must have the same length")
    buf[pos] = b"}" if labels is not None else b"]"
    return _finish_bytes(raw, pos + 1)


# stata, pytables
@cython.boundscheck(False)
@cython.wraparound(False)
def max_len_string_array(pandas_string[:] arr) -> Py_ssize_t:
    """
    Return the maximum size of elements in a 1-dim string array.
    """
    cdef:
        Py_ssize_t i, m = 0, wlen = 0, length = arr.shape[0]
        pandas_string val

    for i in range(length):
        val = arr[i]
        wlen = word_len(val)

        if wlen > m:
            m = wlen

    return m


cpdef inline Py_ssize_t word_len(object val):
    """
    Return the maximum length of a string or bytes value.
    """
    cdef:
        Py_ssize_t wlen = 0

    if isinstance(val, str):
        wlen = PyUnicode_GET_LENGTH(val)
    elif isinstance(val, bytes):
        wlen = PyBytes_GET_SIZE(val)

    return wlen

# ------------------------------------------------------------------
# PyTables Helpers


@cython.boundscheck(False)
@cython.wraparound(False)
def string_array_replace_from_nan_rep(
    ndarray[object, ndim=1] arr,
    object nan_rep,
) -> None:
    """
    Replace the values in the array with np.nan if they are nan_rep.
    """
    cdef:
        Py_ssize_t length = len(arr), i = 0

    for i in range(length):
        if arr[i] == nan_rep:
            arr[i] = np.nan
