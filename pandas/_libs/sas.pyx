# cython: language_level=3, initializedcheck=False
# cython: warn.maybe_uninitialized=True, warn.unused=True
cimport cython
from cpython.unicode cimport (
    PyUnicode_DecodeASCII,
    PyUnicode_DecodeLatin1,
    PyUnicode_FromStringAndSize,
)
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
        # Buffers for the pyarrow large_string fast path; populated when
        # use_pyarrow_strings is set. str_values is one flat utf-8 buffer;
        # col js writes into [str_col_starts[js], str_col_starts[js+1]).
        # str_offsets[js, :] holds large_string offsets (length nrows+1);
        # str_valid[js, r] is 1 for non-null, 0 for null.
        # str_encoding_mode: 0 = utf-8 (copy verbatim), 1 = latin-1 (expand
        # inline), 2 = other single-byte encoding via str_byte_table.
        bint use_pyarrow_strings
        int str_encoding_mode
        uint8_t[::1] str_values
        int64_t[:, ::1] str_offsets
        uint8_t[:, ::1] str_valid
        int64_t[::1] str_col_starts
        # Single-byte -> utf-8 translation (str_encoding_mode == 2):
        # str_byte_table[b, :str_byte_len[b]] is the utf-8 for source byte b;
        # str_byte_len[b] == 0xFF marks a byte undefined in str_table_encoding.
        uint8_t[:, ::1] str_byte_table
        uint8_t[::1] str_byte_len
        object str_table_encoding
        # Object-mode dispatch when use_pyarrow_strings is unset:
        # 0 = raw bytes, 1 = utf-8, 2 = latin-1, 3 = ascii.
        int str_object_mode
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
        self.use_pyarrow_strings = parser._use_pyarrow_strings
        self.str_object_mode = parser._str_object_mode
        if self.use_pyarrow_strings:
            self.str_encoding_mode = parser._str_encoding_mode
            self.str_values = parser._str_values_buf
            self.str_offsets = parser._str_offsets
            self.str_valid = parser._str_valid
            self.str_col_starts = parser._str_col_starts
            # None unless str_encoding_mode == 1 (single-byte table).
            self.str_byte_table = parser._str_byte_table
            self.str_byte_len = parser._str_byte_len
            self.str_table_encoding = parser._str_table_encoding
        else:
            self.string_chunk = parser._string_chunk
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
        buf_free(self.decompressed_buf)

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
            int bi, blen
            int64_t lngt, start, ct
            int64_t write_pos, col_start
            uint8_t bval
            Buffer source
            int64_t[:] column_types
            int64_t[:] lengths
            int64_t[:] offsets
            uint8_t[:, :] byte_chunk
            object[:, :] string_chunk
            uint8_t[::1] str_values
            int64_t[:, ::1] str_offsets
            uint8_t[:, ::1] str_valid
            int64_t[::1] str_col_starts
            uint8_t[:, ::1] str_byte_table
            uint8_t[::1] str_byte_len
            bint compressed
            bint use_pa
            int enc_mode, obj_mode

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
        # Assign every string-output local up front. The inactive path's
        # memoryviews are None and never dereferenced; doing it unconditionally
        # keeps these out of the maybe-uninitialized analysis.
        use_pa = self.use_pyarrow_strings
        enc_mode = self.str_encoding_mode
        obj_mode = self.str_object_mode
        string_chunk = self.string_chunk
        str_values = self.str_values
        str_offsets = self.str_offsets
        str_valid = self.str_valid
        str_col_starts = self.str_col_starts
        str_byte_table = self.str_byte_table
        str_byte_len = self.str_byte_len
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
                    (<uint64_t *>&byte_chunk[jb, m])[0] = (
                        (<uint64_t *>&source.data[start])[0]
                    )
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
                if use_pa:
                    write_pos = str_offsets[js, current_row]
                    col_start = str_col_starts[js]
                    if lngt == 0 and self.blank_missing:
                        str_valid[js, current_row] = 0
                    elif enc_mode == 0:
                        # utf-8: copy bytes verbatim (validated after the
                        # large_string array is built).
                        for k in range(lngt):
                            str_values[col_start + write_pos + k] = (
                                buf_get(source, start + k)
                            )
                        write_pos += lngt
                    elif enc_mode == 1:
                        # latin-1 -> utf-8: 1 byte expands to 1 or 2 bytes
                        for k in range(lngt):
                            bval = buf_get(source, start + k)
                            if bval < 0x80:
                                str_values[col_start + write_pos] = bval
                                write_pos += 1
                            else:
                                str_values[col_start + write_pos] = (
                                    0xC0 | (bval >> 6)
                                )
                                str_values[col_start + write_pos + 1] = (
                                    0x80 | (bval & 0x3F)
                                )
                                write_pos += 2
                    else:
                        # other single-byte encoding -> utf-8 via lookup table
                        for k in range(lngt):
                            bval = buf_get(source, start + k)
                            blen = str_byte_len[bval]
                            if blen == 0xFF:
                                # byte undefined in this encoding; raise the
                                # same error a strict per-cell decode would.
                                bytes([bval]).decode(self.str_table_encoding)
                            for bi in range(blen):
                                str_values[col_start + write_pos + bi] = (
                                    str_byte_table[bval, bi]
                                )
                            write_pos += blen
                    str_offsets[js, current_row + 1] = write_pos
                else:
                    if lngt == 0 and self.blank_missing:
                        string_chunk[js, current_row] = np_nan
                    elif obj_mode == 1:
                        # utf-8: PyUnicode_FromStringAndSize interprets bytes as
                        # utf-8 and raises on invalid sequences.
                        string_chunk[js, current_row] = (
                            PyUnicode_FromStringAndSize(
                                <const char *>&source.data[start], lngt
                            )
                        )
                    elif obj_mode == 2:
                        # latin-1: every byte is a valid codepoint, this call
                        # cannot fail.
                        string_chunk[js, current_row] = PyUnicode_DecodeLatin1(
                            <const char *>&source.data[start], lngt, NULL
                        )
                    elif obj_mode == 3:
                        # ascii: raises on any byte >= 0x80.
                        string_chunk[js, current_row] = PyUnicode_DecodeASCII(
                            <const char *>&source.data[start], lngt, NULL
                        )
                    else:
                        string_chunk[js, current_row] = buf_as_bytes(
                            source, start, lngt
                        )
                js += 1

        self.current_row_on_page_index += 1
        self.current_row_in_chunk_index += 1
        self.current_row_in_file_index += 1
