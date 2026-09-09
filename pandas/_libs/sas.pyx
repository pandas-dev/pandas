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
from libc.string cimport (
    memcpy,
    memset,
)

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


cdef bint buf_copy(
    Buffer dst, size_t dst_offset, Buffer src, size_t src_offset, size_t length
) except 0:
    # Copy a run of bytes, checking the bounds once per run rather than once per
    # byte. `src` and `dst` must not overlap.
    assert src_offset + length <= src.length, "Out of bounds read"
    assert dst_offset + length <= dst.length, "Out of bounds write"
    memcpy(&dst.data[dst_offset], &src.data[src_offset], length)
    return True


cdef bint buf_fill(Buffer buf, size_t offset, size_t length, uint8_t value) except 0:
    # Write a run of `length` copies of the same value, bounds-checked once.
    assert offset + length <= buf.length, "Out of bounds write"
    memset(&buf.data[offset], value, length)
    return True


cdef bint buf_copy_backref(
    Buffer buf, int dst_offset, int src_offset, int length
) except 0:
    # Copy within a single buffer, where the source is allowed to overlap the
    # destination -- that overlap is how an RDC back-reference expands a short
    # pattern -- so this has to stay a forward byte-at-a-time copy. Only the
    # bounds checks are hoisted out of the loop.
    cdef int idx
    assert src_offset >= 0, "Out of bounds read"
    assert <size_t>(dst_offset + length) <= buf.length, "Out of bounds write"
    assert <size_t>(src_offset + length) <= buf.length, "Out of bounds read"
    for idx in range(length):
        buf.data[dst_offset + idx] = buf.data[src_offset + idx]
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

    # A zero capacity means data is still NULL, so allocate even when nothing is
    # needed: an empty first cell would otherwise leave callers memcpy-ing into
    # NULL, which is undefined behavior even for a length of zero.
    if needed <= capacity and capacity > 0:
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
        int nbytes, end_of_first_byte
        size_t ipos = 0

    while ipos < inbuff.length:
        control_byte = buf_get(inbuff, ipos) & 0xF0
        end_of_first_byte = <int>(buf_get(inbuff, ipos) & 0x0F)
        ipos += 1

        if control_byte == 0x00:
            nbytes = <int>(buf_get(inbuff, ipos)) + 64 + end_of_first_byte * 256
            ipos += 1
            buf_copy(outbuff, rpos, inbuff, ipos, nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0x40:
            # not documented
            nbytes = <int>(buf_get(inbuff, ipos)) + 18 + end_of_first_byte * 256
            ipos += 1
            buf_fill(outbuff, rpos, nbytes, buf_get(inbuff, ipos))
            rpos += nbytes
            ipos += 1
        elif control_byte == 0x60:
            nbytes = end_of_first_byte * 256 + <int>(buf_get(inbuff, ipos)) + 17
            ipos += 1
            buf_fill(outbuff, rpos, nbytes, 0x20)
            rpos += nbytes
        elif control_byte == 0x70:
            nbytes = end_of_first_byte * 256 + <int>(buf_get(inbuff, ipos)) + 17
            ipos += 1
            buf_fill(outbuff, rpos, nbytes, 0x00)
            rpos += nbytes
        elif control_byte == 0x80:
            nbytes = end_of_first_byte + 1
            buf_copy(outbuff, rpos, inbuff, ipos, nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0x90:
            nbytes = end_of_first_byte + 17
            buf_copy(outbuff, rpos, inbuff, ipos, nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0xA0:
            nbytes = end_of_first_byte + 33
            buf_copy(outbuff, rpos, inbuff, ipos, nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0xB0:
            nbytes = end_of_first_byte + 49
            buf_copy(outbuff, rpos, inbuff, ipos, nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0xC0:
            nbytes = end_of_first_byte + 3
            x = buf_get(inbuff, ipos)
            ipos += 1
            buf_fill(outbuff, rpos, nbytes, x)
            rpos += nbytes
        elif control_byte == 0xD0:
            nbytes = end_of_first_byte + 2
            buf_fill(outbuff, rpos, nbytes, 0x40)
            rpos += nbytes
        elif control_byte == 0xE0:
            nbytes = end_of_first_byte + 2
            buf_fill(outbuff, rpos, nbytes, 0x20)
            rpos += nbytes
        elif control_byte == 0xF0:
            nbytes = end_of_first_byte + 2
            buf_fill(outbuff, rpos, nbytes, 0x00)
            rpos += nbytes
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
        uint8_t *in_data = inbuff.data
        uint8_t *out_data = outbuff.data
        size_t in_len = inbuff.length
        size_t out_len = outbuff.length
        int rpos = 0
        size_t ipos = 0

    # Unlike RLE, RDC streams are dominated by single-byte literals, so the
    # accessor overhead on each input byte is the cost here rather than the runs;
    # hence the raw reads with the bounds check written out once per command.
    while ipos < in_len:
        ctrl_mask = ctrl_mask >> 1
        if ctrl_mask == 0:
            assert ipos + 2 <= in_len, "Out of bounds read"
            ctrl_bits = ((<uint16_t>in_data[ipos] << 8) +
                         <uint16_t>in_data[ipos + 1])
            ipos += 2
            ctrl_mask = 0x8000

        if ctrl_bits & ctrl_mask == 0:
            assert ipos < in_len, "Out of bounds read"
            assert <size_t>rpos < out_len, "Out of bounds write"
            out_data[rpos] = in_data[ipos]
            ipos += 1
            rpos += 1
            continue

        assert ipos < in_len, "Out of bounds read"
        cmd = (in_data[ipos] >> 4) & 0x0F
        cnt = <uint16_t>(in_data[ipos] & 0x0F)
        ipos += 1

        # short RLE
        if cmd == 0:
            assert ipos < in_len, "Out of bounds read"
            cnt += 3
            buf_fill(outbuff, rpos, cnt, in_data[ipos])
            rpos += cnt
            ipos += 1

        # long RLE
        elif cmd == 1:
            assert ipos + 2 <= in_len, "Out of bounds read"
            cnt += <uint16_t>in_data[ipos] << 4
            cnt += 19
            buf_fill(outbuff, rpos, cnt, in_data[ipos + 1])
            rpos += cnt
            ipos += 2

        # long pattern
        elif cmd == 2:
            assert ipos + 2 <= in_len, "Out of bounds read"
            ofs = cnt + 3
            ofs += <uint16_t>in_data[ipos] << 4
            cnt = <uint16_t>in_data[ipos + 1] + 16
            ipos += 2
            buf_copy_backref(outbuff, rpos, rpos - <int>ofs, cnt)
            rpos += cnt

        # short pattern
        else:
            assert ipos < in_len, "Out of bounds read"
            ofs = cnt + 3
            ofs += <uint16_t>in_data[ipos] << 4
            ipos += 1
            buf_copy_backref(outbuff, rpos, rpos - <int>ofs, cmd)
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
    # Lookups are done in lookup_subheader_index_32/64. The C structures are
    # initialized in _init_subheader_signatures().
    uint32_t subheader_signatures_32bit[13]
    int subheader_indices_32bit[13]
    uint64_t subheader_signatures_64bit[17]
    int subheader_indices_64bit[17]
    int data_subheader_index = const.SASIndex.data_subheader_index

    int c_page_type_offset = const.page_type_offset
    int c_page_type_mask2 = const.page_type_mask2
    int c_block_count_offset = const.block_count_offset
    int c_subheader_count_offset = const.subheader_count_offset
    int c_truncated_subheader_id = const.truncated_subheader_id
    int c_compressed_subheader_id = const.compressed_subheader_id
    int c_compressed_subheader_type = const.compressed_subheader_type


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


cdef int lookup_subheader_index_32(uint32_t sig) noexcept:
    """Fast version of 'subheader_signature_to_index.get(signature)', 4-byte."""
    cdef Py_ssize_t i
    for i in range(len(subheader_signatures_32bit)):
        if subheader_signatures_32bit[i] == sig:
            return subheader_indices_32bit[i]
    return data_subheader_index


cdef int lookup_subheader_index_64(uint64_t sig) noexcept:
    """Fast version of 'subheader_signature_to_index.get(signature)', 8-byte."""
    cdef Py_ssize_t i
    for i in range(len(subheader_signatures_64bit)):
        if subheader_signatures_64bit[i] == sig:
            return subheader_indices_64bit[i]
    return data_subheader_index


cdef uint32_t bswap32(uint32_t val) noexcept:
    return (
        (val >> 24)
        | ((val >> 8) & 0x0000FF00)
        | ((val << 8) & 0x00FF0000)
        | (val << 24)
    )


cdef uint64_t bswap64(uint64_t val) noexcept:
    return (
        (val >> 56)
        | ((val >> 40) & <uint64_t>0xFF00)
        | ((val >> 24) & <uint64_t>0xFF0000)
        | ((val >> 8) & <uint64_t>0xFF000000)
        | ((val << 8) & <uint64_t>0xFF00000000)
        | ((val << 24) & <uint64_t>0xFF0000000000)
        | ((val << 40) & <uint64_t>0xFF000000000000)
        | (val << 56)
    )


@cython.wraparound(False)
@cython.boundscheck(False)
def collect_page_subheaders(parser):
    """
    Walk the subheader pointers on the parser's current cached page.

    Data subheaders (compressed rows) are written into
    parser._data_subheader_offsets / _data_subheader_lengths; every other
    subheader is dispatched to the matching entry of
    parser._subheader_processors. Returns the number of data subheaders found.

    On a compressed file nearly every subheader on a page is a data subheader,
    so keeping this walk in C avoids building a Python list of tuples -- and
    indexing it, with a tuple unpack and two unboxes -- once per row.
    """
    cdef:
        bytes cached_page = parser._cached_page
        const uint8_t *page = <const uint8_t *>cached_page
        size_t page_len = len(cached_page)
        int n_subheaders = parser._current_page_subheaders_count
        int subheader_pointer_length = parser._subheader_pointer_length
        int int_length = parser._int_length
        int bit_offset = parser._page_bit_offset
        bint need_byteswap = parser.need_byteswap
        Py_ssize_t i
        int total_offset
        uint64_t subheader_offset = 0, subheader_length = 0
        uint32_t off32 = 0, len32 = 0, sig32 = 0
        uint64_t sig64 = 0
        uint8_t subheader_compression, subheader_type
        int subheader_index
        int64_t[::1] data_offsets = parser._data_subheader_offsets
        int64_t[::1] data_lengths = parser._data_subheader_lengths
        Py_ssize_t data_capacity = data_offsets.shape[0]
        int count = 0
        int row_start
        list processors = parser._subheader_processors
        object processor
        bint compressed = bool(parser.compression)

    row_start = subheader_pointers_offset + bit_offset

    for i in range(n_subheaders):
        total_offset = row_start + subheader_pointer_length * <int>i
        # These three page-bounds checks raise rather than assert so that a
        # corrupt file keeps reporting as bad input, the way it did when these
        # reads went through _read_bytes, instead of looking like a pandas bug.
        if <size_t>(total_offset + 2 * int_length + 2) > page_len:
            raise ValueError("The cached page is too small.")

        if int_length == 8:
            memcpy(&subheader_offset, &page[total_offset], 8)
            memcpy(&subheader_length, &page[total_offset + 8], 8)
            if need_byteswap:
                subheader_offset = bswap64(subheader_offset)
                subheader_length = bswap64(subheader_length)
        else:
            memcpy(&off32, &page[total_offset], 4)
            memcpy(&len32, &page[total_offset + 4], 4)
            if need_byteswap:
                off32 = bswap32(off32)
                len32 = bswap32(len32)
            subheader_offset = off32
            subheader_length = len32

        subheader_compression = page[total_offset + 2 * int_length]
        subheader_type = page[total_offset + 2 * int_length + 1]

        if subheader_length == 0 or subheader_compression == c_truncated_subheader_id:
            continue

        # Read the int_length-byte signature the pointer points at. Both pointer
        # fields are raw file data, so these subtract rather than add -- adding
        # would wrap in uint64 and let a crafted offset through. The check above
        # guarantees page_len > int_length.
        if subheader_offset > page_len - int_length:
            raise ValueError("The cached page is too small.")
        # Both branches below narrow the length to a C int -- readline() for data
        # subheaders, the processor call for everything else -- so a 64-bit field
        # would otherwise truncate to a wrong but perfectly in-bounds value.
        if subheader_length > page_len - subheader_offset:
            raise ValueError("The cached page is too small.")
        if int_length == 8:
            memcpy(&sig64, &page[subheader_offset], 8)
            subheader_index = lookup_subheader_index_64(sig64)
        else:
            memcpy(&sig32, &page[subheader_offset], 4)
            subheader_index = lookup_subheader_index_32(sig32)

        if subheader_index == data_subheader_index:
            # A compressed data row.
            if not (
                compressed
                and (
                    subheader_compression == c_compressed_subheader_id
                    or subheader_compression == 0
                )
                and subheader_type == c_compressed_subheader_type
            ):
                raise ValueError(
                    f"Unknown subheader signature "
                    f"{page[subheader_offset:subheader_offset + int_length]}"
                )
            if count >= data_capacity:
                raise ValueError(
                    f"Page holds more data subheaders ({count + 1}) than the "
                    f"page can address ({data_capacity})"
                )
            data_offsets[count] = <int64_t>subheader_offset
            data_lengths[count] = <int64_t>subheader_length
            count += 1
        else:
            processor = processors[subheader_index]
            processor(<int>subheader_offset, <int>subheader_length)
            # _process_columntext_subheader can turn compression on partway
            # through a page, so re-read it rather than hoisting once.
            compressed = bool(parser.compression)

    return count


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
        # (offset, length) of each compressed data subheader on the current
        # page, filled in by collect_page_subheaders and indexed by
        # current_row_on_page_index in readline().
        int64_t[::1] data_subheader_offsets
        int64_t[::1] data_subheader_lengths
        uint8_t *cached_page
        Py_ssize_t cached_page_len
        int current_row_on_page_index
        # How many rows a mix page hands out before the parser moves on. Held
        # here rather than read off the reader per row so that a metadata page
        # reached partway through a chunk cannot move the boundary for the rest
        # of it -- sas7bdat.py rejects the file once the chunk is done. An int,
        # like the row index above that it is compared against: a wider one
        # would take a count no page can hold, never reach it, and hand out the
        # rest of the page as rows.
        int mix_page_row_count
        int current_page_block_count
        int n_data_subheaders
        int current_page_subheaders_count
        int current_row_in_chunk_index
        int current_row_in_file_index
        bint blank_missing
        int header_length
        Py_ssize_t row_length
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
        self.data_subheader_offsets = parser._data_subheader_offsets
        self.data_subheader_lengths = parser._data_subheader_lengths
        self.row_length = parser.row_length
        self.mix_page_row_count = min(parser.row_count, parser._mix_page_row_count)
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
                # Set Python attr needed by collect_page_subheaders
                self.parser._current_page_subheaders_count = (
                    self.current_page_subheaders_count
                )
                try:
                    self.n_data_subheaders = collect_page_subheaders(self.parser)
                except Exception:
                    self.parser.close()
                    raise
                self.parser._data_subheader_count = self.n_data_subheaders
                return False
            elif (self.current_page_type == page_data_type
                    or self.current_page_type == page_mix_type):
                self.n_data_subheaders = 0
                self.parser._data_subheader_count = 0
                return False
            # else: unsupported page type (e.g. AMD), skip to next page

    cdef update_next_page(self):
        # Called once from __init__ to sync with the page the Python parser
        # left off on after _parse_metadata.
        self._parse_page_header()
        self.n_data_subheaders = self.parser._data_subheader_count

    cdef bint readline(self) except? True:

        cdef:
            Py_ssize_t offset, length
            int bit_offset, align_correction
            int subheader_pointer_length
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
                flag = self.current_row_on_page_index >= self.n_data_subheaders
                if flag:
                    done = self.read_next_page()
                    if done:
                        return True
                    continue
                offset = <Py_ssize_t>self.data_subheader_offsets[
                    self.current_row_on_page_index
                ]
                length = <Py_ssize_t>self.data_subheader_lengths[
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
                if self.current_row_on_page_index == self.mix_page_row_count:
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
    cdef void process_byte_array_with_data(
        self, Py_ssize_t offset, Py_ssize_t length
    ) except *:

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
