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
)
from libc.string cimport (
    memcpy,
    memset,
)

import numpy as np

import pandas.io.sas.sas_constants as const


cdef object np_nan = np.nan


cdef struct Buffer:
    # Lightweight pointer + length pair. Used for both compressed-page slices
    # (where data is borrowed from the cached page) and the per-row decompression
    # output buffer (heap-allocated by buf_new). Hot paths access .data directly
    # and assert their own bounds rather than going through wrappers — we create
    # Buffer instances roughly once per row, so wrapper overhead is significant.
    uint8_t *data
    size_t length


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
@cython.wraparound(False)
@cython.boundscheck(False)
cdef int rle_decompress(Buffer inbuff, Buffer outbuff) except? 0:

    cdef:
        uint8_t control_byte
        uint8_t *in_data = inbuff.data
        uint8_t *out_data = outbuff.data
        size_t in_len = inbuff.length
        size_t out_len = outbuff.length
        int rpos = 0
        int nbytes, end_of_first_byte
        size_t ipos = 0

    while ipos < in_len:
        control_byte = in_data[ipos] & 0xF0
        end_of_first_byte = <int>(in_data[ipos] & 0x0F)
        ipos += 1

        if control_byte == 0x00:
            assert ipos < in_len, "Out of bounds read"
            nbytes = <int>in_data[ipos] + 64 + end_of_first_byte * 256
            ipos += 1
            assert ipos + nbytes <= in_len, "Out of bounds read"
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memcpy(&out_data[rpos], &in_data[ipos], nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0x40:
            # not documented: 1-byte literal repeated nbytes times
            assert ipos + 1 < in_len, "Out of bounds read"
            nbytes = <int>in_data[ipos] + 18 + end_of_first_byte * 256
            ipos += 1
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], in_data[ipos], nbytes)
            rpos += nbytes
            ipos += 1
        elif control_byte == 0x60:
            assert ipos < in_len, "Out of bounds read"
            nbytes = end_of_first_byte * 256 + <int>in_data[ipos] + 17
            ipos += 1
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], 0x20, nbytes)
            rpos += nbytes
        elif control_byte == 0x70:
            assert ipos < in_len, "Out of bounds read"
            nbytes = end_of_first_byte * 256 + <int>in_data[ipos] + 17
            ipos += 1
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], 0x00, nbytes)
            rpos += nbytes
        elif control_byte == 0x80:
            nbytes = end_of_first_byte + 1
            assert ipos + nbytes <= in_len, "Out of bounds read"
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memcpy(&out_data[rpos], &in_data[ipos], nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0x90:
            nbytes = end_of_first_byte + 17
            assert ipos + nbytes <= in_len, "Out of bounds read"
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memcpy(&out_data[rpos], &in_data[ipos], nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0xA0:
            nbytes = end_of_first_byte + 33
            assert ipos + nbytes <= in_len, "Out of bounds read"
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memcpy(&out_data[rpos], &in_data[ipos], nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0xB0:
            nbytes = end_of_first_byte + 49
            assert ipos + nbytes <= in_len, "Out of bounds read"
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memcpy(&out_data[rpos], &in_data[ipos], nbytes)
            rpos += nbytes
            ipos += nbytes
        elif control_byte == 0xC0:
            assert ipos < in_len, "Out of bounds read"
            nbytes = end_of_first_byte + 3
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], in_data[ipos], nbytes)
            rpos += nbytes
            ipos += 1
        elif control_byte == 0xD0:
            nbytes = end_of_first_byte + 2
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], 0x40, nbytes)
            rpos += nbytes
        elif control_byte == 0xE0:
            nbytes = end_of_first_byte + 2
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], 0x20, nbytes)
            rpos += nbytes
        elif control_byte == 0xF0:
            nbytes = end_of_first_byte + 2
            assert rpos + nbytes <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], 0x00, nbytes)
            rpos += nbytes
        else:
            raise ValueError(f"unknown control byte: {control_byte}")

    return rpos


# rdc_decompress decompresses data using the Ross Data Compression algorithm:
#
# http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/1992/9210/ross/ross.htm
@cython.wraparound(False)
@cython.boundscheck(False)
cdef int rdc_decompress(Buffer inbuff, Buffer outbuff) except? 0:

    cdef:
        uint8_t cmd
        uint16_t ctrl_bits = 0, ctrl_mask = 0, ofs, cnt
        uint8_t *in_data = inbuff.data
        uint8_t *out_data = outbuff.data
        size_t in_len = inbuff.length
        size_t out_len = outbuff.length
        int rpos = 0, k
        size_t ipos = 0

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
            assert rpos < <int>out_len, "Out of bounds write"
            out_data[rpos] = in_data[ipos]
            ipos += 1
            rpos += 1
            continue

        assert ipos < in_len, "Out of bounds read"
        cmd = (in_data[ipos] >> 4) & 0x0F
        cnt = <uint16_t>(in_data[ipos] & 0x0F)
        ipos += 1

        # short RLE: repeat single byte (cnt+3) times
        if cmd == 0:
            assert ipos < in_len, "Out of bounds read"
            cnt += 3
            assert rpos + cnt <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], in_data[ipos], cnt)
            rpos += cnt
            ipos += 1

        # long RLE: repeat single byte (cnt + nextbyte<<4 + 19) times
        elif cmd == 1:
            assert ipos + 1 < in_len, "Out of bounds read"
            cnt += <uint16_t>in_data[ipos] << 4
            cnt += 19
            ipos += 1
            assert rpos + cnt <= <int>out_len, "Out of bounds write"
            memset(&out_data[rpos], in_data[ipos], cnt)
            rpos += cnt
            ipos += 1

        # long pattern: copy cnt bytes from earlier in output
        elif cmd == 2:
            assert ipos + 1 < in_len, "Out of bounds read"
            ofs = cnt + 3
            ofs += <uint16_t>in_data[ipos] << 4
            ipos += 1
            cnt = <uint16_t>in_data[ipos]
            ipos += 1
            cnt += 16
            assert rpos + cnt <= <int>out_len, "Out of bounds write"
            assert rpos >= <int>ofs, "Out of bounds read"
            # Source and destination may overlap when ofs < cnt — the byte loop
            # is intentional so each copied byte feeds the next read (this is
            # how the format encodes repeating patterns).
            if <int>ofs >= <int>cnt:
                memcpy(&out_data[rpos], &out_data[rpos - <int>ofs], cnt)
            else:
                for k in range(cnt):
                    out_data[rpos + k] = out_data[rpos - <int>ofs + k]
            rpos += cnt

        # short pattern: copy cmd bytes from earlier in output
        else:
            assert ipos < in_len, "Out of bounds read"
            ofs = cnt + 3
            ofs += <uint16_t>in_data[ipos] << 4
            ipos += 1
            assert rpos + cmd <= <int>out_len, "Out of bounds write"
            assert rpos >= <int>ofs, "Out of bounds read"
            if <int>ofs >= <int>cmd:
                memcpy(&out_data[rpos], &out_data[rpos - <int>ofs], cmd)
            else:
                for k in range(cmd):
                    out_data[rpos + k] = out_data[rpos - <int>ofs + k]
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


cdef inline int _lookup_subheader_index_32(uint32_t sig) noexcept nogil:
    cdef Py_ssize_t i
    for i in range(13):
        if subheader_signatures_32bit[i] == sig:
            return subheader_indices_32bit[i]
    return data_subheader_index


cdef inline int _lookup_subheader_index_64(uint64_t sig) noexcept nogil:
    cdef Py_ssize_t i
    for i in range(17):
        if subheader_signatures_64bit[i] == sig:
            return subheader_indices_64bit[i]
    return data_subheader_index


cdef inline uint16_t _bswap16(uint16_t val) noexcept nogil:
    return ((val >> 8) | (val << 8)) & 0xFFFF


cdef inline uint32_t _bswap32(uint32_t val) noexcept nogil:
    return (
        (val >> 24)
        | ((val >> 8) & 0x0000FF00)
        | ((val << 8) & 0x00FF0000)
        | (val << 24)
    )


cdef inline uint64_t _bswap64(uint64_t val) noexcept nogil:
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
    Walk the subheader pointers on the parser's current cached page, populating
    parser._data_subheader_offsets and parser._data_subheader_lengths with the
    offset/length of each compressed data subheader.

    For non-data subheaders, calls back into Python to invoke the corresponding
    entry of parser._subheader_processors. This is the cold path on
    data-bearing meta pages of compressed files (which are dominated by data
    subheaders) but the warm path during initial metadata parsing.

    Returns the number of data subheaders collected on this page.
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
        uint64_t subheader_offset, subheader_length
        uint32_t off32, len32, sig32
        uint64_t sig64
        uint8_t subheader_compression, subheader_type
        int subheader_index
        int64_t[::1] data_offsets = parser._data_subheader_offsets
        int64_t[::1] data_lengths = parser._data_subheader_lengths
        int64_t data_capacity = data_offsets.shape[0]
        int count = 0
        int row_start
        list processors = parser._subheader_processors
        object processor

    row_start = subheader_pointers_offset + bit_offset

    for i in range(n_subheaders):
        total_offset = row_start + subheader_pointer_length * i
        assert total_offset + 2 * int_length + 2 <= <int>page_len, (
            "subheader pointer out of bounds"
        )

        if int_length == 8:
            memcpy(&subheader_offset, &page[total_offset], 8)
            memcpy(&subheader_length, &page[total_offset + 8], 8)
            if need_byteswap:
                subheader_offset = _bswap64(subheader_offset)
                subheader_length = _bswap64(subheader_length)
        else:
            memcpy(&off32, &page[total_offset], 4)
            memcpy(&len32, &page[total_offset + 4], 4)
            if need_byteswap:
                off32 = _bswap32(off32)
                len32 = _bswap32(len32)
            subheader_offset = off32
            subheader_length = len32

        subheader_compression = page[total_offset + 2 * int_length]
        subheader_type = page[total_offset + 2 * int_length + 1]

        if subheader_length == 0 or subheader_compression == c_truncated_subheader_id:
            continue

        # Read the int_length-byte signature at subheader_offset.
        assert subheader_offset + int_length <= page_len, "signature out of bounds"
        if int_length == 4:
            memcpy(&sig32, &page[subheader_offset], 4)
            subheader_index = _lookup_subheader_index_32(sig32)
        else:
            memcpy(&sig64, &page[subheader_offset], 8)
            subheader_index = _lookup_subheader_index_64(sig64)

        if subheader_index == data_subheader_index:
            # Data subheader (compressed row). Validate compression markers
            # and append into the int64 arrays. parser.compression is read
            # dynamically because it can be set to a non-empty literal by
            # _process_columntext_subheader earlier on this same page.
            if not (
                parser.compression
                and (
                    subheader_compression == c_compressed_subheader_id
                    or subheader_compression == 0
                )
                and subheader_type == c_compressed_subheader_type
            ):
                raise ValueError(
                    f"Unknown subheader signature "
                    f"{bytes(page[subheader_offset:subheader_offset + int_length])}"
                )
            if count >= data_capacity:
                raise ValueError(
                    f"More data subheaders ({count + 1}) than preallocated "
                    f"capacity ({data_capacity}) — this is a pandas bug"
                )
            data_offsets[count] = <int64_t>subheader_offset
            data_lengths[count] = <int64_t>subheader_length
            count += 1
        else:
            # Cold path: dispatch to the Python processor for this index.
            processor = processors[subheader_index]
            if processor is None:
                raise ValueError(
                    f"Unknown subheader signature "
                    f"{bytes(page[subheader_offset:subheader_offset + int_length])}"
                )
            processor(<int>subheader_offset, <int>subheader_length)

    return count


cdef class Parser:

    cdef:
        int column_count
        int64_t[:] lengths
        int64_t[:] offsets
        int64_t[:] column_types
        uint8_t[:, :] byte_chunk
        object[:, :] string_chunk
        # Data subheader (offset, length) pairs found on the current page.
        # Populated in Cython by collect_page_subheaders; indexed by
        # current_row_on_page_index in readline().
        int64_t[::1] data_subheader_offsets
        int64_t[::1] data_subheader_lengths
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
        self.string_chunk = parser._string_chunk
        self.data_subheader_offsets = parser._data_subheader_offsets
        self.data_subheader_lengths = parser._data_subheader_lengths
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
                self.parser._current_page_subheaders_count = (
                    self.current_page_subheaders_count
                )
                self.current_page_data_subheader_pointers_len = (
                    collect_page_subheaders(self.parser)
                )
                self.parser._data_subheader_count = (
                    self.current_page_data_subheader_pointers_len
                )
                return False
            elif (self.current_page_type == page_data_type
                    or self.current_page_type == page_mix_type):
                self.current_page_data_subheader_pointers_len = 0
                self.parser._data_subheader_count = 0
                return False
            # else: unsupported page type (e.g. AMD), skip to next page

    cdef update_next_page(self):
        # Called once from __init__ to sync with the page the Python parser
        # left off on after _parse_metadata.
        self._parse_page_header()
        self.current_page_data_subheader_pointers_len = (
            self.parser._data_subheader_count
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
                offset = <int>self.data_subheader_offsets[
                    self.current_row_on_page_index
                ]
                length = <int>self.data_subheader_lengths[
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
            int s, m, jb, js, current_row, rpos
            int64_t lngt, start, ct
            Buffer source
            uint8_t *src_data
            uint8_t last_byte
            int64_t[:] column_types
            int64_t[:] lengths
            int64_t[:] offsets
            uint8_t[:, :] byte_chunk
            object[:, :] string_chunk
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

        src_data = source.data
        current_row = self.current_row_in_chunk_index
        column_types = self.column_types
        lengths = self.lengths
        offsets = self.offsets
        byte_chunk = self.byte_chunk
        string_chunk = self.string_chunk
        s = 8 * self.current_row_in_chunk_index
        js = 0
        jb = 0
        for j in range(self.column_count):
            lngt = lengths[j]
            if lngt == 0:
                break
            start = offsets[j]
            ct = column_types[j]
            assert start + lngt <= <int64_t>source.length, "Out of bounds read"
            if ct == column_type_decimal:
                # decimal: copy lngt bytes into the right-aligned slot for
                # little-endian, left-aligned for big-endian. memcpy with a
                # small dynamic length compiles to a branch-free load/store
                # on common CPUs and generalizes the lngt==8 fast path.
                if self.is_little_endian:
                    m = s + 8 - <int>lngt
                else:
                    m = s
                memcpy(&byte_chunk[jb, m], &src_data[start], lngt)
                jb += 1
            elif ct == column_type_string:
                # string: skip trailing 0x00 / 0x20 (equivalent to
                # .rstrip(b"\x00 ") but without Python call overhead).
                while lngt > 0:
                    last_byte = src_data[start + lngt - 1]
                    if last_byte != 0x00 and last_byte != 0x20:
                        break
                    lngt -= 1
                if lngt == 0 and self.blank_missing:
                    string_chunk[js, current_row] = np_nan
                else:
                    string_chunk[js, current_row] = src_data[start:start + lngt]
                js += 1

        self.current_row_on_page_index += 1
        self.current_row_in_chunk_index += 1
        self.current_row_in_file_index += 1
