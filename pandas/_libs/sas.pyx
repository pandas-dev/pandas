# cython: language_level=3, initializedcheck=False
# cython: warn.maybe_uninitialized=True, warn.unused=True
cimport cython
from cython cimport Py_ssize_t
from libc.stddef cimport size_t
from libc.stdint cimport (
    int64_t,
    uint8_t,
    uint16_t,
    uint32_t,
    uint64_t,
)
from libc.stdlib cimport (
    calloc,
    free,
    realloc,
)
from libc.string cimport memcpy

import numpy as np

import pandas.io.sas.sas_constants as const


cdef object np_nan = np.nan


cdef struct Buffer:
    # Convenience wrapper for uint8_t data to allow fast and safe reads and writes.
    # We use this as a replacement for np.array(..., dtype=np.uint8) because it's
    # much slower to create NumPy arrays and we create Buffer instances many times
    # when reading a SAS7BDAT file (roughly once per row that is being read).
    uint8_t *data
    size_t length


cdef uint8_t buf_get(Buffer buf, size_t offset) except? 255:
    assert offset < buf.length, "Out of bounds read"
    return buf.data[offset]


cdef bint buf_set(Buffer buf, size_t offset, uint8_t value) except 0:
    assert offset < buf.length, "Out of bounds write"
    buf.data[offset] = value
    return True


cdef bytes buf_as_bytes(Buffer buf, size_t offset, size_t length):
    assert offset + length <= buf.length, "Out of bounds read"
    return buf.data[offset:offset+length]


cdef Buffer buf_new(size_t length) except *:
    cdef uint8_t *data = <uint8_t *>calloc(length, sizeof(uint8_t))
    if data is NULL:
        raise MemoryError(f"Failed to allocate {length} bytes")
    return Buffer(data, length)


cdef buf_free(Buffer buf):
    if buf.data != NULL:
        free(buf.data)


cdef struct StrColumn:
    # Growable utf-8 value buffer backing one string column.
    #
    # Grown by doubling rather than preallocated to the declared column width
    # times nrows: SAS pads character fields out to a fixed width, so the
    # worst case badly overstates how many bytes the data actually needs.
    uint8_t *data
    int64_t size
    int64_t capacity


cdef bint strcol_reserve(StrColumn *col, int64_t extra) except 0:
    """Ensure `col` has room for `extra` more bytes, growing by doubling."""
    cdef:
        int64_t needed = col.size + extra
        int64_t capacity = col.capacity
        uint8_t *data

    if needed <= capacity:
        return True
    if capacity == 0:
        capacity = 4096
    while capacity < needed:
        capacity *= 2
    data = <uint8_t *>realloc(col.data, capacity)
    if data is NULL:
        raise MemoryError(f"Failed to allocate {capacity} bytes")
    col.data = data
    col.capacity = capacity
    return True

# rle_decompress decompresses data using a Run Length Encoding
# algorithm.  It is partially documented here:
#
# https://cran.r-project.org/package=sas7bdat/vignettes/sas7bdat.pdf
# Licence at LICENSES/SAS7BDAT_LICENSE
cdef int rle_decompress(Buffer inbuff, Buffer outbuff) except? 0:

    cdef:
        uint8_t control_byte, x
        int rpos = 0
        int i, nbytes, end_of_first_byte
        size_t ipos = 0
        Py_ssize_t _

    while ipos < inbuff.length:
        control_byte = buf_get(inbuff, ipos) & 0xF0
        end_of_first_byte = <int>(buf_get(inbuff, ipos) & 0x0F)
        ipos += 1

        if control_byte == 0x00:
            nbytes = <int>(buf_get(inbuff, ipos)) + 64 + end_of_first_byte * 256
            ipos += 1
            for _ in range(nbytes):
                buf_set(outbuff, rpos, buf_get(inbuff, ipos))
                rpos += 1
                ipos += 1
        elif control_byte == 0x40:
            # not documented
            nbytes = <int>(buf_get(inbuff, ipos)) + 18 + end_of_first_byte * 256
            ipos += 1
            for _ in range(nbytes):
                buf_set(outbuff, rpos, buf_get(inbuff, ipos))
                rpos += 1
            ipos += 1
        elif control_byte == 0x60:
            nbytes = end_of_first_byte * 256 + <int>(buf_get(inbuff, ipos)) + 17
            ipos += 1
            for _ in range(nbytes):
                buf_set(outbuff, rpos, 0x20)
                rpos += 1
        elif control_byte == 0x70:
            nbytes = end_of_first_byte * 256 + <int>(buf_get(inbuff, ipos)) + 17
            ipos += 1
            for _ in range(nbytes):
                buf_set(outbuff, rpos, 0x00)
                rpos += 1
        elif control_byte == 0x80:
            nbytes = end_of_first_byte + 1
            for i in range(nbytes):
                buf_set(outbuff, rpos, buf_get(inbuff, ipos + i))
                rpos += 1
            ipos += nbytes
        elif control_byte == 0x90:
            nbytes = end_of_first_byte + 17
            for i in range(nbytes):
                buf_set(outbuff, rpos, buf_get(inbuff, ipos + i))
                rpos += 1
            ipos += nbytes
        elif control_byte == 0xA0:
            nbytes = end_of_first_byte + 33
            for i in range(nbytes):
                buf_set(outbuff, rpos, buf_get(inbuff, ipos + i))
                rpos += 1
            ipos += nbytes
        elif control_byte == 0xB0:
            nbytes = end_of_first_byte + 49
            for i in range(nbytes):
                buf_set(outbuff, rpos, buf_get(inbuff, ipos + i))
                rpos += 1
            ipos += nbytes
        elif control_byte == 0xC0:
            nbytes = end_of_first_byte + 3
            x = buf_get(inbuff, ipos)
            ipos += 1
            for _ in range(nbytes):
                buf_set(outbuff, rpos, x)
                rpos += 1
        elif control_byte == 0xD0:
            nbytes = end_of_first_byte + 2
            for _ in range(nbytes):
                buf_set(outbuff, rpos, 0x40)
                rpos += 1
        elif control_byte == 0xE0:
            nbytes = end_of_first_byte + 2
            for _ in range(nbytes):
                buf_set(outbuff, rpos, 0x20)
                rpos += 1
        elif control_byte == 0xF0:
            nbytes = end_of_first_byte + 2
            for _ in range(nbytes):
                buf_set(outbuff, rpos, 0x00)
                rpos += 1
        else:
            raise ValueError(f"unknown control byte: {control_byte}")

    return rpos


# rdc_decompress decompresses data using the Ross Data Compression algorithm:
#
# http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/1992/9210/ross/ross.htm
cdef int rdc_decompress(Buffer inbuff, Buffer outbuff) except? 0:

    cdef:
        uint8_t cmd
        uint16_t ctrl_bits = 0, ctrl_mask = 0, ofs, cnt
        int rpos = 0, k, ii
        size_t ipos = 0

    ii = -1

    while ipos < inbuff.length:
        ii += 1
        ctrl_mask = ctrl_mask >> 1
        if ctrl_mask == 0:
            ctrl_bits = ((<uint16_t>buf_get(inbuff, ipos) << 8) +
                         <uint16_t>buf_get(inbuff, ipos + 1))
            ipos += 2
            ctrl_mask = 0x8000

        if ctrl_bits & ctrl_mask == 0:
            buf_set(outbuff, rpos, buf_get(inbuff, ipos))
            ipos += 1
            rpos += 1
            continue

        cmd = (buf_get(inbuff, ipos) >> 4) & 0x0F
        cnt = <uint16_t>(buf_get(inbuff, ipos) & 0x0F)
        ipos += 1

        # short RLE
        if cmd == 0:
            cnt += 3
            for k in range(cnt):
                buf_set(outbuff, rpos + k, buf_get(inbuff, ipos))
            rpos += cnt
            ipos += 1

        # long RLE
        elif cmd == 1:
            cnt += <uint16_t>buf_get(inbuff, ipos) << 4
            cnt += 19
            ipos += 1
            for k in range(cnt):
                buf_set(outbuff, rpos + k, buf_get(inbuff, ipos))
            rpos += cnt
            ipos += 1

        # long pattern
        elif cmd == 2:
            ofs = cnt + 3
            ofs += <uint16_t>buf_get(inbuff, ipos) << 4
            ipos += 1
            cnt = <uint16_t>buf_get(inbuff, ipos)
            ipos += 1
            cnt += 16
            for k in range(cnt):
                buf_set(outbuff, rpos + k, buf_get(outbuff, rpos - <int>ofs + k))
            rpos += cnt

        # short pattern
        else:
            ofs = cnt + 3
            ofs += <uint16_t>buf_get(inbuff, ipos) << 4
            ipos += 1
            for k in range(cmd):
                buf_set(outbuff, rpos + k, buf_get(outbuff, rpos - <int>ofs + k))
            rpos += cmd

    return rpos


cdef enum ColumnTypes:
    column_type_decimal = 1
    column_type_string = 2


# Const aliases
assert len(const.page_meta_types) == 2
cdef:
    int page_meta_types_0 = const.page_meta_types[0]
    int page_meta_types_1 = const.page_meta_types[1]
    int page_mix_type = const.page_mix_type
    int page_data_type = const.page_data_type
    int subheader_pointers_offset = const.subheader_pointers_offset

    # See sas_constants.py for what each string mode means. string_mode_table
    # needs no alias here: it is whatever the other two are not.
    int string_mode_object = const.string_mode_object
    int string_mode_utf8 = const.string_mode_utf8
    uint8_t undefined_byte = const.undefined_byte

    # Copy of subheader_signature_to_index that allows for much faster lookups.
    # Lookups are done in get_subheader_index. The C structures are initialized
    # in _init_subheader_signatures().
    uint32_t subheader_signatures_32bit[13]
    int subheader_indices_32bit[13]
    uint64_t subheader_signatures_64bit[17]
    int subheader_indices_64bit[17]
    int data_subheader_index = const.SASIndex.data_subheader_index

    int c_page_type_offset = const.page_type_offset
    int c_page_type_mask2 = const.page_type_mask2
    int c_block_count_offset = const.block_count_offset
    int c_subheader_count_offset = const.subheader_count_offset


def _init_subheader_signatures():
    subheaders_32bit = [
        (sig, idx)
        for sig, idx in const.subheader_signature_to_index.items()
        if len(sig) == 4
    ]
    subheaders_64bit = [
        (sig, idx)
        for sig, idx in const.subheader_signature_to_index.items()
        if len(sig) == 8
    ]
    assert len(subheaders_32bit) == 13
    assert len(subheaders_64bit) == 17
    assert len(const.subheader_signature_to_index) == 13 + 17
    for i, (signature, idx) in enumerate(subheaders_32bit):
        subheader_signatures_32bit[i] = (<uint32_t *><char *>signature)[0]
        subheader_indices_32bit[i] = idx
    for i, (signature, idx) in enumerate(subheaders_64bit):
        subheader_signatures_64bit[i] = (<uint64_t *><char *>signature)[0]
        subheader_indices_64bit[i] = idx


_init_subheader_signatures()


def get_subheader_index(bytes signature):
    """Fast version of 'subheader_signature_to_index.get(signature)'."""
    cdef:
        uint32_t sig32
        uint64_t sig64
        Py_ssize_t i
    assert len(signature) in (4, 8)
    if len(signature) == 4:
        sig32 = (<uint32_t *><char *>signature)[0]
        for i in range(len(subheader_signatures_32bit)):
            if subheader_signatures_32bit[i] == sig32:
                return subheader_indices_32bit[i]
    else:
        sig64 = (<uint64_t *><char *>signature)[0]
        for i in range(len(subheader_signatures_64bit)):
            if subheader_signatures_64bit[i] == sig64:
                return subheader_indices_64bit[i]

    return data_subheader_index


cdef class Parser:

    cdef:
        int column_count
        int64_t[:] lengths
        int64_t[:] offsets
        int64_t[:] column_types
        uint8_t[:, :] byte_chunk
        object[:, :] string_chunk
        # Output for the pyarrow string fast path (string_mode_utf8/_table):
        # str_cols[js] accumulates column js's utf-8 bytes, str_offsets[js, r+1]
        # records the end of row r within it, and str_valid[js, r] is 0 for a
        # null cell. Unused when str_mode is string_mode_object.
        int str_mode
        StrColumn *str_cols
        int n_str_cols
        int64_t[:, ::1] str_offsets
        uint8_t[:, ::1] str_valid
        # str_table[b * str_table_width : ...] holds the utf-8 for source byte
        # b, str_table_len[b] its length (undefined_byte if b is not a valid
        # character in str_encoding). Only set for string_mode_table. const
        # because these are shared, read-only arrays cached across readers.
        const uint8_t[::1] str_table
        const uint8_t[::1] str_table_len
        int str_table_width
        bint str_ascii_identity
        object str_encoding
        uint8_t *cached_page
        int cached_page_len
        int current_row_on_page_index
        int current_page_block_count
        int current_page_data_subheader_pointers_len
        int current_page_subheaders_count
        int current_row_in_chunk_index
        int current_row_in_file_index
        bint blank_missing
        int header_length
        int row_length
        int bit_offset
        int subheader_pointer_length
        int current_page_type
        bint is_little_endian
        bint need_byteswap
        int (*decompress)(Buffer, Buffer) except? 0
        Buffer decompressed_buf
        object parser

    def __cinit__(self, object parser):
        self.decompress = NULL
        self.decompressed_buf = Buffer(NULL, 0)
        self.str_cols = NULL
        self.n_str_cols = 0

    @cython.wraparound(False)
    @cython.boundscheck(False)
    def __init__(self, object parser):
        cdef:
            int j
            char[:] column_types

        self.parser = parser
        self.blank_missing = parser.blank_missing
        self.header_length = self.parser.header_length
        self.column_count = parser.column_count
        self.lengths = parser.column_data_lengths()
        self.offsets = parser.column_data_offsets()
        self.byte_chunk = parser._byte_chunk
        self.string_chunk = parser._string_chunk
        # sas7bdat.py always supplies these, using empty arrays for whichever
        # of the two string outputs is inactive, so the memoryview assignments
        # below never see None.
        self.str_mode = parser._str_mode
        self.str_offsets = parser._str_offsets
        self.str_valid = parser._str_valid
        self.str_table = parser._str_table
        self.str_table_len = parser._str_table_len
        self.str_table_width = parser._str_table_width
        self.str_ascii_identity = parser._str_ascii_identity
        self.str_encoding = parser._str_encoding
        if self.str_mode != string_mode_object:
            self.n_str_cols = self.str_offsets.shape[0]
            self.str_cols = <StrColumn *>calloc(self.n_str_cols, sizeof(StrColumn))
            if self.str_cols is NULL:
                raise MemoryError("Failed to allocate string column buffers")
        self.row_length = parser.row_length
        self.bit_offset = self.parser._page_bit_offset
        self.subheader_pointer_length = self.parser._subheader_pointer_length
        self.is_little_endian = parser.byte_order == "<"
        self.need_byteswap = parser.need_byteswap
        self.column_types = np.empty(self.column_count, dtype="int64")

        # page indicators
        self.update_next_page()

        column_types = parser.column_types()

        # map column types
        for j in range(self.column_count):
            if column_types[j] == b"d":
                self.column_types[j] = column_type_decimal
            elif column_types[j] == b"s":
                self.column_types[j] = column_type_string
            else:
                raise ValueError(f"unknown column type: {self.parser.columns[j].ctype}")

        # compression
        if parser.compression == const.rle_compression:
            self.decompress = rle_decompress
        elif parser.compression == const.rdc_compression:
            self.decompress = rdc_decompress
        else:
            self.decompress = NULL

        if self.decompress != NULL:
            self.decompressed_buf = buf_new(self.row_length)

        # update to current state of the parser
        self.current_row_in_chunk_index = parser._current_row_in_chunk_index
        self.current_row_in_file_index = parser._current_row_in_file_index
        self.current_row_on_page_index = parser._current_row_on_page_index

    def __dealloc__(self):
        cdef Py_ssize_t js
        buf_free(self.decompressed_buf)
        if self.str_cols is not NULL:
            for js in range(self.n_str_cols):
                free(self.str_cols[js].data)
            free(self.str_cols)

    def string_values(self):
        """
        Return one exact-size uint8 array of utf-8 bytes per string column.

        Copies out of the growable parse buffers so the returned arrays do not
        pin their (up to 2x) spare capacity for the lifetime of the result.
        """
        cdef:
            Py_ssize_t js
            StrColumn *col
            uint8_t[::1] out

        result = []
        for js in range(self.n_str_cols):
            col = &self.str_cols[js]
            values = np.empty(col.size, dtype=np.uint8)
            if col.size > 0:
                out = values
                memcpy(&out[0], col.data, col.size)
            result.append(values)
        return result

    def read(self, int nrows):
        cdef:
            bint done
            Py_ssize_t _

        for _ in range(nrows):
            done = self.readline()
            if done:
                break

        # update the parser
        self.parser._current_row_on_page_index = self.current_row_on_page_index
        self.parser._current_row_in_chunk_index = self.current_row_in_chunk_index
        self.parser._current_row_in_file_index = self.current_row_in_file_index

    cdef uint16_t read_uint16(self, int offset):
        cdef uint16_t val = 0
        assert offset + 2 <= self.cached_page_len, "Out of bounds read"
        memcpy(&val, &self.cached_page[offset], sizeof(uint16_t))
        if self.need_byteswap:
            val = ((val >> 8) | (val << 8)) & 0xFFFF
        return val

    cdef void _parse_page_header(self):
        # Read page header fields directly from the cached page buffer in C.
        self.cached_page = <uint8_t *>self.parser._cached_page
        self.cached_page_len = len(self.parser._cached_page)
        self.current_row_on_page_index = 0
        self.current_page_type = (
            self.read_uint16(c_page_type_offset + self.bit_offset)
            & c_page_type_mask2
        )
        self.current_page_block_count = self.read_uint16(
            c_block_count_offset + self.bit_offset
        )
        self.current_page_subheaders_count = self.read_uint16(
            c_subheader_count_offset + self.bit_offset
        )

    cdef bint read_next_page(self) except? True:
        cdef bint done

        while True:
            done = self.parser._read_page_data()
            if done:
                self.cached_page = NULL
                return True

            self._parse_page_header()

            if (self.current_page_type == page_meta_types_0
                    or self.current_page_type == page_meta_types_1):
                # Set Python attr needed by _process_page_metadata
                self.parser._current_page_subheaders_count = (
                    self.current_page_subheaders_count
                )
                self.parser._process_page_metadata()
                self.current_page_data_subheader_pointers_len = len(
                    self.parser._current_page_data_subheader_pointers
                )
                return False
            elif (self.current_page_type == page_data_type
                    or self.current_page_type == page_mix_type):
                self.current_page_data_subheader_pointers_len = 0
                return False
            # else: unsupported page type (e.g. AMD), skip to next page

    cdef update_next_page(self):
        # Called once from __init__ to sync with the page the Python parser
        # left off on after _parse_metadata.
        self._parse_page_header()
        self.current_page_data_subheader_pointers_len = len(
            self.parser._current_page_data_subheader_pointers
        )

    cdef bint readline(self) except? True:

        cdef:
            int offset, length, bit_offset, align_correction
            int subheader_pointer_length, mn
            bint done, flag

        bit_offset = self.bit_offset
        subheader_pointer_length = self.subheader_pointer_length

        # If there is no page, go to the end of the header and read a page.
        if self.cached_page == NULL:
            self.parser._path_or_buf.seek(self.header_length)
            done = self.read_next_page()
            if done:
                return True

        # Loop until a data row is read
        while True:
            if self.current_page_type in (page_meta_types_0, page_meta_types_1):
                flag = self.current_row_on_page_index >=\
                    self.current_page_data_subheader_pointers_len
                if flag:
                    done = self.read_next_page()
                    if done:
                        return True
                    continue
                offset, length = self.parser._current_page_data_subheader_pointers[
                    self.current_row_on_page_index
                ]
                self.process_byte_array_with_data(offset, length)
                return False
            elif self.current_page_type == page_mix_type:
                align_correction = (
                    bit_offset
                    + subheader_pointers_offset
                    + self.current_page_subheaders_count * subheader_pointer_length
                )
                align_correction = align_correction % 8
                offset = bit_offset + align_correction
                offset += subheader_pointers_offset
                offset += self.current_page_subheaders_count * subheader_pointer_length
                offset += self.current_row_on_page_index * self.row_length
                self.process_byte_array_with_data(offset, self.row_length)
                mn = min(self.parser.row_count, self.parser._mix_page_row_count)
                if self.current_row_on_page_index == mn:
                    done = self.read_next_page()
                    if done:
                        return True
                return False
            elif self.current_page_type == page_data_type:
                self.process_byte_array_with_data(
                    bit_offset
                    + subheader_pointers_offset
                    + self.current_row_on_page_index * self.row_length,
                    self.row_length,
                )
                flag = self.current_row_on_page_index == self.current_page_block_count
                if flag:
                    done = self.read_next_page()
                    if done:
                        return True
                return False
            else:
                raise ValueError(f"unknown page type: {self.current_page_type}")

    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef void process_byte_array_with_data(self, int offset, int length) except *:

        cdef:
            Py_ssize_t j
            int s, k, m, jb, js, current_row, rpos
            int bi, blen, table_width
            int64_t lngt, start, ct
            uint8_t bval
            bint ascii_identity, mode_object, mode_utf8
            Buffer source
            StrColumn *col
            int64_t[:] column_types
            int64_t[:] lengths
            int64_t[:] offsets
            uint8_t[:, :] byte_chunk
            object[:, :] string_chunk
            StrColumn *str_cols
            int64_t[:, ::1] str_offsets
            uint8_t[:, ::1] str_valid
            const uint8_t[::1] str_table
            const uint8_t[::1] str_table_len
            bint compressed

        assert offset + length <= self.cached_page_len, "Out of bounds read"
        source = Buffer(&self.cached_page[offset], length)

        compressed = self.decompress != NULL and length < self.row_length
        if compressed:
            rpos = self.decompress(source, self.decompressed_buf)
            if rpos != self.row_length:
                raise ValueError(
                    f"Expected decompressed line of length {self.row_length} bytes "
                    f"but decompressed {rpos} bytes"
                )
            source = self.decompressed_buf

        current_row = self.current_row_in_chunk_index
        column_types = self.column_types
        lengths = self.lengths
        offsets = self.offsets
        byte_chunk = self.byte_chunk
        string_chunk = self.string_chunk
        # Resolved once per row so that the per-cell branches below test a
        # local rather than reloading the module-level mode constants.
        mode_object = self.str_mode == string_mode_object
        mode_utf8 = self.str_mode == string_mode_utf8
        str_cols = self.str_cols
        str_offsets = self.str_offsets
        str_valid = self.str_valid
        str_table = self.str_table
        str_table_len = self.str_table_len
        table_width = self.str_table_width
        ascii_identity = self.str_ascii_identity
        s = 8 * self.current_row_in_chunk_index
        js = 0
        jb = 0
        for j in range(self.column_count):
            lngt = lengths[j]
            if lngt == 0:
                break
            start = offsets[j]
            ct = column_types[j]
            if ct == column_type_decimal:
                # decimal — copy lngt (3..8) bytes from source to byte_chunk.
                # Fast path for the common case of full-precision SAS doubles
                # (lngt == 8): a single 8-byte load/store, ~2x faster than
                # the per-byte fallback.
                if self.is_little_endian:
                    m = s + 8 - lngt
                else:
                    m = s
                if lngt == 8:
                    memcpy(&byte_chunk[jb, m], &source.data[start], sizeof(uint64_t))
                else:
                    for k in range(lngt):
                        byte_chunk[jb, m + k] = buf_get(source, start + k)
                jb += 1
            elif column_types[j] == column_type_string:
                # string
                # Skip trailing whitespace. This is equivalent to calling
                # .rstrip(b"\x00 ") but without Python call overhead.
                while lngt > 0 and buf_get(source, start + lngt - 1) in b"\x00 ":
                    lngt -= 1
                if mode_object:
                    if lngt == 0 and self.blank_missing:
                        string_chunk[js, current_row] = np_nan
                    else:
                        string_chunk[js, current_row] = buf_as_bytes(
                            source, start, lngt
                        )
                else:
                    col = &str_cols[js]
                    if lngt > 0 or not self.blank_missing:
                        # Spell out the bound the raw-pointer reads below rely
                        # on; the rstrip loop's buf_get established it.
                        assert start + lngt <= <int64_t>source.length, \
                            "Out of bounds read"
                        # Cells start null (str_valid is zero-initialized), so
                        # only non-null ones need marking.
                        str_valid[js, current_row] = 1
                        strcol_reserve(col, lngt * table_width)
                        if mode_utf8:
                            memcpy(col.data + col.size, &source.data[start], lngt)
                            col.size += lngt
                        else:
                            k = 0
                            if ascii_identity:
                                # Bytes below 0x80 are their own utf-8 in this
                                # encoding, so an all-ascii cell (the common
                                # case) costs one memcpy instead of a byte loop.
                                while k < lngt and source.data[start + k] < 0x80:
                                    k += 1
                                memcpy(col.data + col.size, &source.data[start], k)
                                col.size += k
                            while k < lngt:
                                bval = source.data[start + k]
                                blen = str_table_len[bval]
                                if blen == undefined_byte:
                                    # Re-decode the cell so the error matches
                                    # the UnicodeDecodeError that the per-cell
                                    # strict decode used to raise.
                                    buf_as_bytes(source, start, lngt).decode(
                                        self.str_encoding
                                    )
                                    raise AssertionError(
                                        "decoding an undefined byte should have raised"
                                    )
                                for bi in range(blen):
                                    col.data[col.size + bi] = (
                                        str_table[bval * table_width + bi]
                                    )
                                col.size += blen
                                k += 1
                    str_offsets[js, current_row + 1] = col.size
                js += 1

        self.current_row_on_page_index += 1
        self.current_row_in_chunk_index += 1
        self.current_row_in_file_index += 1
