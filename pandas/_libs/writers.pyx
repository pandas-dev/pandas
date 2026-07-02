cimport cython
from cython cimport Py_ssize_t

import numpy as np

cimport numpy as cnp
from cpython cimport (
    PyBytes_GET_SIZE,
    PyUnicode_GET_LENGTH,
)
from cpython.mem cimport (
    PyMem_Free,
    PyMem_Malloc,
    PyMem_Realloc,
)
from cpython.object cimport PyObject
from cpython.unicode cimport PyUnicode_DecodeUTF8
from libc.string cimport (
    memchr,
    memcpy,
    strchr,
)
from numpy cimport (
    int64_t,
    ndarray,
    uint8_t,
    uint64_t,
)

cnp.import_array()

from pandas._libs.tslibs.nattype cimport NPY_NAT
from pandas._libs.tslibs.np_datetime cimport (
    NPY_DATETIMEUNIT,
    import_pandas_datetime,
    npy_datetimestruct,
    pandas_datetime_to_datetimestruct,
)

import_pandas_datetime()

cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8AndSize(object unicode, Py_ssize_t* size) except NULL

cdef extern from "pandas/double_repr.h":
    int pd_double_repr(double val, char *out)
    int pd_float32_repr(float val, char *out)
    int pd_float32_repr_available()

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


# ------------------------------------------------------------------
# Fast whole-chunk CSV rendering (bypasses csv.writer)

# Column kind codes understood by write_csv_chunk; mirrored in
# pandas.io.formats.csvs. CSV_KIND_OBJ columns hold pre-converted cells
# (via get_values_for_csv); the other kinds are rendered natively here.
CSV_KIND_OBJ = 0
CSV_KIND_FLOAT64 = 1
CSV_KIND_FLOAT32 = 2
CSV_KIND_INT64 = 3
CSV_KIND_UINT64 = 4
CSV_KIND_BOOL = 5
CSV_KIND_DT64 = 6

# float32 rendering needs floating-point std::to_chars (no fallback)
FLOAT32_NATIVE = bool(pd_float32_repr_available())

cdef enum:
    KIND_OBJ = 0
    KIND_FLOAT64 = 1
    KIND_FLOAT32 = 2
    KIND_INT64 = 3
    KIND_UINT64 = 4
    KIND_BOOL = 5
    KIND_DT64 = 6

# All bytes a natively-rendered field can contain, per kind. If the
# delimiter/quotechar cannot occur in the rendered text, the per-cell
# quoting scan is skipped ('\r'/'\n' can never occur).
cdef const char* KIND_CHARSET_FLOAT = "0123456789.+-einfa"
cdef const char* KIND_CHARSET_INT = "-0123456789"
cdef const char* KIND_CHARSET_UINT = "0123456789"
cdef const char* KIND_CHARSET_BOOL = "TrueFals"
cdef const char* KIND_CHARSET_DT64 = "0123456789-: ."


cdef struct CsvBuf:
    char* data
    size_t size
    size_t cap


cdef int buf_grow(CsvBuf* buf, size_t needed) except -1:
    cdef:
        size_t new_cap = buf.cap
        char* new_data
    while new_cap < buf.size + needed:
        new_cap += new_cap // 2 + 64
    new_data = <char*>PyMem_Realloc(buf.data, new_cap)
    if new_data == NULL:
        raise MemoryError()
    buf.data = new_data
    buf.cap = new_cap
    return 0


cdef inline int buf_reserve(CsvBuf* buf, size_t needed) except -1:
    if buf.cap - buf.size < needed:
        buf_grow(buf, needed)
    return 0


cdef inline void buf_append_raw(CsvBuf* buf, const char* src, size_t n) noexcept:
    memcpy(buf.data + buf.size, src, n)
    buf.size += n


cdef int append_field(CsvBuf* buf, const char* field, Py_ssize_t n,
                      char sep, char quote, bint lone_field) except -1:
    """
    Append one field with csv.writer QUOTE_MINIMAL semantics: quote (doubling
    embedded quotechars) iff the field contains the delimiter, the quotechar,
    or a literal CR/LF; quote an empty field iff it is the row's only field.
    """
    cdef:
        Py_ssize_t i
        bint need_quote

    if n == 0:
        if lone_field:
            buf_reserve(buf, 2)
            buf.data[buf.size] = quote
            buf.data[buf.size + 1] = quote
            buf.size += 2
        return 0

    need_quote = (
        memchr(field, sep, n) != NULL
        or memchr(field, quote, n) != NULL
        or memchr(field, c"\r", n) != NULL
        or memchr(field, c"\n", n) != NULL
    )
    if not need_quote:
        buf_reserve(buf, n)
        buf_append_raw(buf, field, n)
    else:
        buf_reserve(buf, 2 * n + 2)
        buf.data[buf.size] = quote
        buf.size += 1
        for i in range(n):
            if field[i] == quote:
                buf.data[buf.size] = quote
                buf.size += 1
            buf.data[buf.size] = field[i]
            buf.size += 1
        buf.data[buf.size] = quote
        buf.size += 1
    return 0


cdef int append_obj_cell(CsvBuf* buf, object cell, char sep, char quote,
                         bint lone_field) except -1:
    """
    Append a pre-converted cell like csv.writer would: str used as-is,
    None as empty, anything else via str().
    """
    cdef:
        const char* cstr
        Py_ssize_t length
        object tmp

    if cell is None:
        return append_field(buf, NULL, 0, sep, quote, lone_field)
    if isinstance(cell, str):
        cstr = PyUnicode_AsUTF8AndSize(cell, &length)
        return append_field(buf, cstr, length, sep, quote, lone_field)
    tmp = str(cell)
    cstr = PyUnicode_AsUTF8AndSize(tmp, &length)
    return append_field(buf, cstr, length, sep, quote, lone_field)


cdef inline int u64toa(uint64_t val, char* out) noexcept:
    cdef:
        char tmp[20]
        int pos = 0, i = 0

    if val == 0:
        out[0] = c"0"
        return 1
    while val != 0:
        tmp[pos] = <char>(c"0" + val % 10)
        val //= 10
        pos += 1
    while pos > 0:
        pos -= 1
        out[i] = tmp[pos]
        i += 1
    return i


cdef inline int i64toa(int64_t val, char* out) noexcept:
    cdef uint64_t uval
    if val >= 0:
        return u64toa(<uint64_t>val, out)
    out[0] = c"-"
    # two's complement negation handles INT64_MIN
    uval = ~(<uint64_t>val) + 1
    return 1 + u64toa(uval, out + 1)


cdef inline void write2(char* out, int val) noexcept:
    out[0] = <char>(c"0" + val // 10)
    out[1] = <char>(c"0" + val % 10)


cdef inline void write_padded(char* out, int64_t val, int width) noexcept:
    cdef int i
    for i in range(width - 1, -1, -1):
        out[i] = <char>(c"0" + val % 10)
        val //= 10


cdef int render_dt64(int64_t val, NPY_DATETIMEUNIT reso, int frac_digits,
                     bint dates_only, char* out) noexcept:
    """
    Render a datetime64 the way tslib.format_array_from_datetime does for
    tz-naive values: "YYYY-MM-DD[ HH:MM:SS[.fff[fff[fff]]]]" with an
    unpadded (possibly negative) year and 0/3/6/9 fractional digits.
    """
    cdef:
        npy_datetimestruct dts
        int64_t frac
        int pos

    pandas_datetime_to_datetimestruct(val, reso, &dts)
    pos = i64toa(dts.year, out)
    out[pos] = c"-"
    write2(out + pos + 1, dts.month)
    out[pos + 3] = c"-"
    write2(out + pos + 4, dts.day)
    pos += 6
    if dates_only:
        return pos
    out[pos] = c" "
    write2(out + pos + 1, dts.hour)
    out[pos + 3] = c":"
    write2(out + pos + 4, dts.min)
    out[pos + 6] = c":"
    write2(out + pos + 7, dts.sec)
    pos += 9
    if frac_digits == 0:
        return pos
    if frac_digits == 9:
        frac = dts.us * 1000 + dts.ps // 1000
    elif frac_digits == 6:
        frac = dts.us
    else:
        frac = dts.us // 1000
    out[pos] = c"."
    write_padded(out + pos + 1, frac, frac_digits)
    return pos + 1 + frac_digits


cdef inline bint kind_needs_scan(int kind, char sep, char quote) noexcept:
    """
    Whether natively-rendered fields of this kind can contain the delimiter
    or quotechar and therefore need the per-cell quoting scan.
    """
    cdef const char* charset
    if kind == KIND_FLOAT64 or kind == KIND_FLOAT32:
        charset = KIND_CHARSET_FLOAT
    elif kind == KIND_INT64:
        charset = KIND_CHARSET_INT
    elif kind == KIND_UINT64:
        charset = KIND_CHARSET_UINT
    elif kind == KIND_BOOL:
        charset = KIND_CHARSET_BOOL
    else:
        charset = KIND_CHARSET_DT64
    return strchr(charset, sep) != NULL or strchr(charset, quote) != NULL


@cython.boundscheck(False)
@cython.wraparound(False)
def write_csv_chunk(
    list cols,
    Py_ssize_t nrows,
    str sep,
    str quotechar,
    str lineterminator,
    str na_rep,
) -> str:
    """
    Render a chunk of to_csv output to a single str, replicating
    csv.writer(..., quoting=QUOTE_MINIMAL, doublequote=True) byte for byte.

    Parameters
    ----------
    cols : list of (kind, ndarray, creso, frac_digits, dates_only) tuples,
        one per output field (index column, if any, first). kind is one of
        the CSV_KIND_* constants; creso/frac_digits/dates_only are only
        meaningful for CSV_KIND_DT64.
    nrows : int
    sep : str
        Single ASCII character.
    quotechar : str
        Single ASCII character.
    lineterminator : str
        One or two ASCII characters.
    na_rep : str
        Written for NaN/NaT in natively-rendered columns.

    Returns
    -------
    str

    Raises
    ------
    UnicodeEncodeError
        For strings that cannot be encoded to UTF-8 (lone surrogates);
        the caller falls back to the legacy writer for the chunk.
    """
    cdef:
        Py_ssize_t nfields = len(cols), j, k
        bint lone_field = nfields == 1
        char csep = <char>ord(sep)
        char cquote = <char>ord(quotechar)
        char lt0 = 0, lt1 = 0
        int lt_len, kind, slen
        bytes na_field_b
        const char* na_ptr
        Py_ssize_t na_len
        int* kinds = NULL
        char** data_ptrs = NULL
        Py_ssize_t* strides = NULL
        int* resos = NULL
        int* fracs = NULL
        char* donly = NULL
        char* scans = NULL
        CsvBuf buf
        char scratch[64]
        double dval
        float fval
        int64_t ival
        ndarray arr
        object cell, result

    lt_bytes = lineterminator.encode("ascii")
    lt_len = len(lt_bytes)
    if lt_len < 1 or lt_len > 2:
        raise ValueError(f"invalid lineterminator: {lineterminator!r}")
    lt0 = (<const char*>lt_bytes)[0]
    if lt_len == 2:
        lt1 = (<const char*>lt_bytes)[1]

    # pre-render the field emitted for missing values, quoting rules applied
    if na_rep == "":
        na_field_b = (quotechar * 2).encode("utf-8") if lone_field else b""
    elif (
        sep in na_rep
        or quotechar in na_rep
        or "\r" in na_rep
        or "\n" in na_rep
    ):
        na_field_b = (
            quotechar + na_rep.replace(quotechar, quotechar * 2) + quotechar
        ).encode("utf-8")
    else:
        na_field_b = na_rep.encode("utf-8")
    na_ptr = <const char*>na_field_b
    na_len = len(na_field_b)

    buf.data = NULL
    buf.size = 0
    buf.cap = 0

    kinds = <int*>PyMem_Malloc(nfields * sizeof(int))
    data_ptrs = <char**>PyMem_Malloc(nfields * sizeof(char*))
    strides = <Py_ssize_t*>PyMem_Malloc(nfields * sizeof(Py_ssize_t))
    resos = <int*>PyMem_Malloc(nfields * sizeof(int))
    fracs = <int*>PyMem_Malloc(nfields * sizeof(int))
    donly = <char*>PyMem_Malloc(nfields * sizeof(char))
    scans = <char*>PyMem_Malloc(nfields * sizeof(char))
    if (
        (nfields > 0)
        and (
            kinds == NULL or data_ptrs == NULL or strides == NULL or resos == NULL
            or fracs == NULL or donly == NULL or scans == NULL
        )
    ):
        PyMem_Free(kinds)
        PyMem_Free(data_ptrs)
        PyMem_Free(strides)
        PyMem_Free(resos)
        PyMem_Free(fracs)
        PyMem_Free(donly)
        PyMem_Free(scans)
        raise MemoryError()

    try:
        for k in range(nfields):
            kind, arr, creso, frac_digits, dates_only = cols[k]
            if arr.ndim != 1 or nrows > arr.shape[0]:
                raise ValueError("invalid column array")
            if kind == KIND_OBJ:
                if arr.dtype != object:
                    raise ValueError("CSV_KIND_OBJ requires object dtype")
            elif kind == KIND_FLOAT64:
                if arr.dtype != np.float64:
                    raise ValueError("CSV_KIND_FLOAT64 requires float64")
            elif kind == KIND_FLOAT32:
                if arr.dtype != np.float32:
                    raise ValueError("CSV_KIND_FLOAT32 requires float32")
                if not FLOAT32_NATIVE:
                    raise ValueError("float32 rendering unavailable")
            elif kind == KIND_INT64:
                if arr.dtype != np.int64:
                    raise ValueError("CSV_KIND_INT64 requires int64")
            elif kind == KIND_UINT64:
                if arr.dtype != np.uint64:
                    raise ValueError("CSV_KIND_UINT64 requires uint64")
            elif kind == KIND_BOOL:
                if arr.dtype != np.bool_:
                    raise ValueError("CSV_KIND_BOOL requires bool")
            elif kind == KIND_DT64:
                if arr.dtype != np.int64:
                    raise ValueError("CSV_KIND_DT64 requires int64 view")
            else:
                raise ValueError(f"unknown column kind {kind}")
            kinds[k] = kind
            data_ptrs[k] = <char*>cnp.PyArray_DATA(arr)
            strides[k] = cnp.PyArray_STRIDE(arr, 0)
            resos[k] = creso
            fracs[k] = frac_digits
            donly[k] = dates_only
            scans[k] = kind != KIND_OBJ and kind_needs_scan(kind, csep, cquote)

        buf.cap = <size_t>nrows * (<size_t>nfields * 8 + 8) + 64
        buf.data = <char*>PyMem_Malloc(buf.cap)
        if buf.data == NULL:
            raise MemoryError()

        for j in range(nrows):
            for k in range(nfields):
                if k > 0:
                    buf_reserve(&buf, 1)
                    buf.data[buf.size] = csep
                    buf.size += 1
                kind = kinds[k]
                if kind == KIND_OBJ:
                    cell = <object>(
                        (<PyObject**>(data_ptrs[k] + j * strides[k]))[0]
                    )
                    append_obj_cell(&buf, cell, csep, cquote, lone_field)
                    continue

                slen = -1
                if kind == KIND_FLOAT64:
                    dval = (<double*>(data_ptrs[k] + j * strides[k]))[0]
                    if dval != dval:
                        slen = 0
                    else:
                        slen = pd_double_repr(dval, scratch)
                        if slen < 0:
                            raise ValueError("float formatting failed")
                elif kind == KIND_FLOAT32:
                    fval = (<float*>(data_ptrs[k] + j * strides[k]))[0]
                    if fval != fval:
                        slen = 0
                    else:
                        slen = pd_float32_repr(fval, scratch)
                        if slen < 0:
                            raise ValueError("float formatting failed")
                elif kind == KIND_INT64:
                    ival = (<int64_t*>(data_ptrs[k] + j * strides[k]))[0]
                    slen = i64toa(ival, scratch)
                elif kind == KIND_UINT64:
                    slen = u64toa(
                        (<uint64_t*>(data_ptrs[k] + j * strides[k]))[0], scratch
                    )
                elif kind == KIND_BOOL:
                    if (<uint8_t*>(data_ptrs[k] + j * strides[k]))[0]:
                        memcpy(scratch, b"True", 4)
                        slen = 4
                    else:
                        memcpy(scratch, b"False", 5)
                        slen = 5
                else:  # KIND_DT64
                    ival = (<int64_t*>(data_ptrs[k] + j * strides[k]))[0]
                    if ival == NPY_NAT:
                        slen = 0
                    else:
                        slen = render_dt64(
                            ival,
                            <NPY_DATETIMEUNIT>resos[k],
                            fracs[k],
                            donly[k],
                            scratch,
                        )

                if slen == 0:
                    # missing value: na_field_b already has quoting (and the
                    # lone-empty-field rule) applied
                    if na_len > 0:
                        buf_reserve(&buf, na_len)
                        buf_append_raw(&buf, na_ptr, na_len)
                elif scans[k]:
                    append_field(&buf, scratch, slen, csep, cquote, lone_field)
                else:
                    buf_reserve(&buf, slen)
                    buf_append_raw(&buf, scratch, slen)

            buf_reserve(&buf, 2)
            buf.data[buf.size] = lt0
            buf.size += 1
            if lt_len == 2:
                buf.data[buf.size] = lt1
                buf.size += 1

        result = PyUnicode_DecodeUTF8(buf.data, buf.size, NULL)
    finally:
        PyMem_Free(buf.data)
        PyMem_Free(kinds)
        PyMem_Free(data_ptrs)
        PyMem_Free(strides)
        PyMem_Free(resos)
        PyMem_Free(fracs)
        PyMem_Free(donly)
        PyMem_Free(scans)

    return result


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
