"""
Read SAS7BDAT files

Based on code written by Jared Hobbs:
  https://bitbucket.org/jaredhobbs/sas7bdat

See also:
  https://github.com/BioStatMatt/sas7bdat

Partial documentation of the file format:
  https://cran.r-project.org/web/packages/sas7bdat/vignettes/sas7bdat.pdf

Reference for binary data compression:
  http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/1992/9210/ross/ross.htm
"""

import pandas as pd
from pandas import compat
from pandas.io.common import get_filepath_or_buffer, BaseIterator
import numpy as np
import struct
from .saslib import (_rle_decompress, _rdc_decompress,
                     process_byte_array_with_data)

_magic = (b"\x00\x00\x00\x00\x00\x00\x00\x00" +
          b"\x00\x00\x00\x00\xc2\xea\x81\x60" +
          b"\xb3\x14\x11\xcf\xbd\x92\x08\x00" +
          b"\x09\xc7\x31\x8c\x18\x1f\x10\x11")

_align_1_checker_value = b'3'
_align_1_offset = 32
_align_1_length = 1
_align_1_value = 4
_u64_byte_checker_value = b'3'
_align_2_offset = 35
_align_2_length = 1
_align_2_value = 4
_endianness_offset = 37
_endianness_length = 1
_platform_offset = 39
_platform_length = 1
_encoding_offset = 70
_encoding_length = 1
_dataset_offset = 92
_dataset_length = 64
_file_type_offset = 156
_file_type_length = 8
_date_created_offset = 164
_date_created_length = 8
_date_modified_offset = 172
_date_modified_length = 8
_header_size_offset = 196
_header_size_length = 4
_page_size_offset = 200
_page_size_length = 4
_page_count_offset = 204
_page_count_length = 4
_sas_release_offset = 216
_sas_release_length = 8
_sas_server_type_offset = 224
_sas_server_type_length = 16
_os_version_number_offset = 240
_os_version_number_length = 16
_os_maker_offset = 256
_os_maker_length = 16
_os_name_offset = 272
_os_name_length = 16
_page_bit_offset_x86 = 16
_page_bit_offset_x64 = 32
_subheader_pointer_length_x86 = 12
_subheader_pointer_length_x64 = 24
_page_type_offset = 0
_page_type_length = 2
_block_count_offset = 2
_block_count_length = 2
_subheader_count_offset = 4
_subheader_count_length = 2
_page_meta_type = 0
_page_data_type = 256
_page_amd_type = 1024
_page_metc_type = 16384
_page_comp_type = -28672
_page_mix_types = [512, 640]
_subheader_pointers_offset = 8
_truncated_subheader_id = 1
_compressed_subheader_id = 4
_compressed_subheader_type = 1
_text_block_size_length = 2
_row_length_offset_multiplier = 5
_row_count_offset_multiplier = 6
_col_count_p1_multiplier = 9
_col_count_p2_multiplier = 10
_row_count_on_mix_page_offset_multiplier = 15
_column_name_pointer_length = 8
_column_name_text_subheader_offset = 0
_column_name_text_subheader_length = 2
_column_name_offset_offset = 2
_column_name_offset_length = 2
_column_name_length_offset = 4
_column_name_length_length = 2
_column_data_offset_offset = 8
_column_data_length_offset = 8
_column_data_length_length = 4
_column_type_offset = 14
_column_type_length = 1
_column_format_text_subheader_index_offset = 22
_column_format_text_subheader_index_length = 2
_column_format_offset_offset = 24
_column_format_offset_length = 2
_column_format_length_offset = 26
_column_format_length_length = 2
_column_label_text_subheader_index_offset = 28
_column_label_text_subheader_index_length = 2
_column_label_offset_offset = 30
_column_label_offset_length = 2
_column_label_length_offset = 32
_column_label_length_length = 2
_rle_compression = 'SASYZCRL'
_rdc_compression = 'SASYZCR2'

_compression_literals = [_rle_compression, _rdc_compression]

# Incomplete list of encodings
_encoding_names = {29: "latin1", 20: "utf-8", 33: "cyrillic", 60: "wlatin2",
                   61: "wcyrillic", 62: "wlatin1", 90: "ebcdic870"}

# Should be enum


class _index:
    rowSizeIndex = 0
    columnSizeIndex = 1
    subheaderCountsIndex = 2
    columnTextIndex = 3
    columnNameIndex = 4
    columnAttributesIndex = 5
    formatAndLabelIndex = 6
    columnListIndex = 7
    dataSubheaderIndex = 8


_subheader_signature_to_index = {
    b"\xF7\xF7\xF7\xF7": _index.rowSizeIndex,
    b"\x00\x00\x00\x00\xF7\xF7\xF7\xF7": _index.rowSizeIndex,
    b"\xF7\xF7\xF7\xF7\x00\x00\x00\x00": _index.rowSizeIndex,
    b"\xF7\xF7\xF7\xF7\xFF\xFF\xFB\xFE": _index.rowSizeIndex,
    b"\xF6\xF6\xF6\xF6": _index.columnSizeIndex,
    b"\x00\x00\x00\x00\xF6\xF6\xF6\xF6": _index.columnSizeIndex,
    b"\xF6\xF6\xF6\xF6\x00\x00\x00\x00": _index.columnSizeIndex,
    b"\xF6\xF6\xF6\xF6\xFF\xFF\xFB\xFE": _index.columnSizeIndex,
    b"\x00\xFC\xFF\xFF": _index.subheaderCountsIndex,
    b"\xFF\xFF\xFC\x00": _index.subheaderCountsIndex,
    b"\x00\xFC\xFF\xFF\xFF\xFF\xFF\xFF": _index.subheaderCountsIndex,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFC\x00": _index.subheaderCountsIndex,
    b"\xFD\xFF\xFF\xFF": _index.columnTextIndex,
    b"\xFF\xFF\xFF\xFD": _index.columnTextIndex,
    b"\xFD\xFF\xFF\xFF\xFF\xFF\xFF\xFF": _index.columnTextIndex,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFD": _index.columnTextIndex,
    b"\xFF\xFF\xFF\xFF": _index.columnNameIndex,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF": _index.columnNameIndex,
    b"\xFC\xFF\xFF\xFF": _index.columnAttributesIndex,
    b"\xFF\xFF\xFF\xFC": _index.columnAttributesIndex,
    b"\xFC\xFF\xFF\xFF\xFF\xFF\xFF\xFF": _index.columnAttributesIndex,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFC": _index.columnAttributesIndex,
    b"\xFE\xFB\xFF\xFF": _index.formatAndLabelIndex,
    b"\xFF\xFF\xFB\xFE": _index.formatAndLabelIndex,
    b"\xFE\xFB\xFF\xFF\xFF\xFF\xFF\xFF": _index.formatAndLabelIndex,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFB\xFE": _index.formatAndLabelIndex,
    b"\xFE\xFF\xFF\xFF": _index.columnListIndex,
    b"\xFF\xFF\xFF\xFE": _index.columnListIndex,
    b"\xFE\xFF\xFF\xFF\xFF\xFF\xFF\xFF": _index.columnListIndex,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFE": _index.columnListIndex}


class _subheader_pointer(object):
    pass


class _column(object):
    pass


# SAS7BDAT represents a SAS data file in SAS7BDAT format.
class SAS7BDATReader(BaseIterator):
    """
    Read SAS files in SAS7BDAT format.

    Parameters
    ----------
    path_or_buf : path name or buffer
        Name of SAS file or file-like object pointing to SAS file
        contents.
    index : column identifier, defaults to None
        Column to use as index.
    convert_dates : boolean, defaults to True
        Attempt to convert dates to Pandas datetime values.  Note all
        SAS date formats are supported.
    blank_missing : boolean, defaults to True
        Convert empty strings to missing values (SAS uses blanks to
        indicate missing character variables).
    chunksize : int, defaults to None
        Return SAS7BDATReader object for iterations, returns chunks
        with given number of lines.
    encoding : string, defaults to None
        String encoding.  If None, text variables are left as raw bytes.
    """

    def __init__(self, path_or_buf, index=None, convert_dates=True,
                 blank_missing=True, chunksize=None, encoding=None):

        self.index = index
        self.convert_dates = convert_dates
        self.blank_missing = blank_missing
        self.chunksize = chunksize
        self.encoding = encoding

        self.compression = ""
        self.column_names_strings = []
        self.column_names = []
        self.column_types = []
        self.column_formats = []
        self.columns = []

        self._current_page_data_subheader_pointers = []
        self._cached_page = None
        self._column_data_lengths = []
        self._column_data_offsets = []
        self._current_row_in_file_index = 0
        self._current_row_on_page_index = 0
        self._current_row_in_file_index = 0

        self._path_or_buf, _, _ = get_filepath_or_buffer(path_or_buf)
        if isinstance(self._path_or_buf, compat.string_types):
            self._path_or_buf = open(self._path_or_buf, 'rb')

        self._get_properties()
        self._parse_metadata()

    def _get_properties(self):

        # Check magic number
        self._path_or_buf.seek(0)
        self._cached_page = self._path_or_buf.read(288)
        if self._cached_page[0:len(_magic)] != _magic:
            raise ValueError("magic number mismatch (not a SAS file?)")

        # Get alignment information
        align1, align2 = 0, 0
        buf = self._read_bytes(_align_1_offset, _align_1_length)
        if buf == _u64_byte_checker_value:
            align2 = _align_2_value
            self.U64 = True
            self._int_length = 8
            self._page_bit_offset = _page_bit_offset_x64
            self._subheader_pointer_length = _subheader_pointer_length_x64
        else:
            self.U64 = False
            self._page_bit_offset = _page_bit_offset_x86
            self._subheader_pointer_length = _subheader_pointer_length_x86
            self._int_length = 4
        buf = self._read_bytes(_align_2_offset, _align_2_length)
        if buf == _align_1_checker_value:
            align1 = _align_2_value
        total_align = align1 + align2

        # Get endianness information
        buf = self._read_bytes(_endianness_offset, _endianness_length)
        if buf == b'\x01':
            self.byte_order = "<"
        else:
            self.byte_order = ">"

        # Get encoding information
        buf = self._read_bytes(_encoding_offset, _encoding_length)[0]
        if buf in _encoding_names:
            self.file_encoding = _encoding_names[buf]
        else:
            self.file_encoding = "unknown (code=%s)" % str(buf)

        # Get platform information
        buf = self._read_bytes(_platform_offset, _platform_length)
        if buf == b'1':
            self.platform = "unix"
        elif buf == b'2':
            self.platform = "windows"
        else:
            self.platform = "unknown"

        buf = self._read_bytes(_dataset_offset, _dataset_length)
        self.name = buf.rstrip(b'\x00 ').decode()

        buf = self._read_bytes(_file_type_offset, _file_type_length)
        self.file_type = buf.rstrip(b'\x00 ').decode()

        # Timestamp is epoch 01/01/1960
        epoch = pd.datetime(1960, 1, 1)
        x = self._read_float(_date_created_offset + align1,
                             _date_created_length)
        self.date_created = epoch + pd.to_timedelta(x, unit='s')
        x = self._read_float(_date_modified_offset + align1,
                             _date_modified_length)
        self.date_modified = epoch + pd.to_timedelta(x, unit='s')

        self.header_length = self._read_int(_header_size_offset + align1,
                                            _header_size_length)

        # Read the rest of the header into cached_page.
        buf = self._path_or_buf.read(self.header_length - 288)
        self._cached_page += buf
        if len(self._cached_page) != self.header_length:
            raise ValueError("The SAS7BDAT file appears to be truncated.")

        self._page_length = self._read_int(_page_size_offset + align1,
                                           _page_size_length)
        self._page_count = self._read_int(_page_count_offset + align1,
                                          _page_count_length)

        buf = self._read_bytes(_sas_release_offset + total_align,
                               _sas_release_length)
        self.sas_release = buf.rstrip(b'\x00 ').decode()

        buf = self._read_bytes(_sas_server_type_offset + total_align,
                               _sas_server_type_length)
        self.server_type = buf.rstrip(b'\x00 ').decode()

        buf = self._read_bytes(_os_version_number_offset + total_align,
                               _os_version_number_length)
        self.os_version = buf.rstrip(b'\x00 ').decode()

        buf = self._read_bytes(
            _os_name_offset, _os_name_length).rstrip(b'\x00 ')
        if len(buf) > 0:
            self.os_name = buf.rstrip(b'\x00 ').decode()
        else:
            buf = self._path_or_buf.read(_os_maker_offset, _os_maker_length)
            self.os_name = buf.rstrip(b'\x00 ').decode()

    # Read a single float of the given width (4 or 8).
    def _read_float(self, offset, width):
        if width not in (4, 8):
            raise ValueError("invalid float width")
        buf = self._read_bytes(offset, width)
        fd = "f" if width == 4 else "d"
        return struct.unpack(self.byte_order + fd, buf)[0]

    # Read a single signed integer of the given width (1, 2, 4 or 8).
    def _read_int(self, offset, width):
        if width not in (1, 2, 4, 8):
            raise ValueError("invalid int width")
        buf = self._read_bytes(offset, width)
        it = {1: "b", 2: "h", 4: "l", 8: "q"}[width]
        iv = struct.unpack(self.byte_order + it, buf)[0]
        return iv

    def _read_bytes(self, offset, length):
        if self._cached_page is None:
            self._path_or_buf.seek(offset)
            buf = self._path_or_buf.read(length)
            if len(buf) < length:
                msg = "Unable to read {:d} bytes from file position {:d}."
                raise ValueError(msg.format(length, offset))
            return buf
        else:
            if offset + length > len(self._cached_page):
                raise ValueError("The cached page is too small.")
            return self._cached_page[offset:offset + length]

    def _parse_metadata(self):
        done = False
        while not done:
            self._cached_page = self._path_or_buf.read(self._page_length)
            if len(self._cached_page) <= 0:
                break
            if len(self._cached_page) != self._page_length:
                raise ValueError(
                    "Failed to read a meta data page from the SAS file.")
            done = self._process_page_meta()

    def _process_page_meta(self):
        self._read_page_header()
        pt = [_page_meta_type, _page_amd_type] + _page_mix_types
        if self._current_page_type in pt:
            self._process_page_metadata()
        return ((self._current_page_type in [256] + _page_mix_types) or
                (self._current_page_data_subheader_pointers is not None))

    def _read_page_header(self):
        bit_offset = self._page_bit_offset
        tx = _page_type_offset + bit_offset
        self._current_page_type = self._read_int(tx, _page_type_length)
        tx = _block_count_offset + bit_offset
        self._current_page_block_count = self._read_int(tx,
                                                        _block_count_length)
        tx = _subheader_count_offset + bit_offset
        self._current_page_subheaders_count = (
            self._read_int(tx, _subheader_count_length))

    def _process_page_metadata(self):
        bit_offset = self._page_bit_offset

        for i in range(self._current_page_subheaders_count):
            pointer = self._process_subheader_pointers(
                _subheader_pointers_offset + bit_offset, i)
            if pointer.length == 0:
                continue
            if pointer.compression == _truncated_subheader_id:
                continue
            subheader_signature = self._read_subheader_signature(
                pointer.offset)
            subheader_index = (
                self._get_subheader_index(subheader_signature,
                                          pointer.compression, pointer.ptype))
            self._process_subheader(subheader_index, pointer)

    def _get_subheader_index(self, signature, compression, ptype):
        index = _subheader_signature_to_index.get(signature)
        if index is None:
            f1 = ((compression == _compressed_subheader_id) or
                  (compression == 0))
            f2 = (ptype == _compressed_subheader_type)
            if (self.compression != "") and f1 and f2:
                index = _index.dataSubheaderIndex
            else:
                raise ValueError("Unknown subheader signature")
        return index

    def _process_subheader_pointers(self, offset, subheader_pointer_index):

        subheader_pointer_length = self._subheader_pointer_length
        total_offset = (offset +
                        subheader_pointer_length * subheader_pointer_index)

        subheader_offset = self._read_int(total_offset, self._int_length)
        total_offset += self._int_length

        subheader_length = self._read_int(total_offset, self._int_length)
        total_offset += self._int_length

        subheader_compression = self._read_int(total_offset, 1)
        total_offset += 1

        subheader_type = self._read_int(total_offset, 1)

        x = _subheader_pointer()
        x.offset = subheader_offset
        x.length = subheader_length
        x.compression = subheader_compression
        x.ptype = subheader_type

        return x

    def _read_subheader_signature(self, offset):
        subheader_signature = self._read_bytes(offset, self._int_length)
        return subheader_signature

    def _process_subheader(self, subheader_index, pointer):
        offset = pointer.offset
        length = pointer.length

        if subheader_index == _index.rowSizeIndex:
            processor = self._process_rowsize_subheader
        elif subheader_index == _index.columnSizeIndex:
            processor = self._process_columnsize_subheader
        elif subheader_index == _index.columnTextIndex:
            processor = self._process_columntext_subheader
        elif subheader_index == _index.columnNameIndex:
            processor = self._process_columnname_subheader
        elif subheader_index == _index.columnAttributesIndex:
            processor = self._process_columnattributes_subheader
        elif subheader_index == _index.formatAndLabelIndex:
            processor = self._process_format_subheader
        elif subheader_index == _index.columnListIndex:
            processor = self._process_columnlist_subheader
        elif subheader_index == _index.subheaderCountsIndex:
            processor = self._process_subheader_counts
        elif subheader_index == _index.dataSubheaderIndex:
            self._current_page_data_subheader_pointers.append(pointer)
            return
        else:
            raise ValueError("unknown subheader index")

        processor(offset, length)

    def _process_rowsize_subheader(self, offset, length):

        int_len = self._int_length
        lcs_offset = offset
        lcp_offset = offset
        if self.U64:
            lcs_offset += 682
            lcp_offset += 706
        else:
            lcs_offset += 354
            lcp_offset += 378

        self.row_length = self._read_int(
            offset + _row_length_offset_multiplier * int_len, int_len)
        self.row_count = self._read_int(
            offset + _row_count_offset_multiplier * int_len, int_len)
        self.col_count_p1 = self._read_int(
            offset + _col_count_p1_multiplier * int_len, int_len)
        self.col_count_p2 = self._read_int(
            offset + _col_count_p2_multiplier * int_len, int_len)
        mx = _row_count_on_mix_page_offset_multiplier * int_len
        self._mix_page_row_count = self._read_int(offset + mx, int_len)
        self._lcs = self._read_int(lcs_offset, 2)
        self._lcp = self._read_int(lcp_offset, 2)

    def _process_columnsize_subheader(self, offset, length):
        int_len = self._int_length
        offset += int_len
        self.column_count = self._read_int(offset, int_len)
        if (self.col_count_p1 + self.col_count_p2 !=
                self.column_count):
            print("Warning: column count mismatch (%d + %d != %d)\n",
                  self.col_count_p1, self.col_count_p2, self.column_count)

    # Unknown purpose
    def _process_subheader_counts(self, offset, length):
        pass

    def _process_columntext_subheader(self, offset, length):

        offset += self._int_length
        text_block_size = self._read_int(offset, _text_block_size_length)

        buf = self._read_bytes(offset, text_block_size)
        self.column_names_strings.append(
            buf[0:text_block_size].rstrip(b"\x00 ").decode())

        if len(self.column_names_strings) == 1:
            column_name = self.column_names_strings[0]
            compression_literal = ""
            for cl in _compression_literals:
                if cl in column_name:
                    compression_literal = cl
            self.compression = compression_literal
            offset -= self._int_length

            offset1 = offset + 16
            if self.U64:
                offset1 += 4

            buf = self._read_bytes(offset1, self._lcp)
            compression_literal = buf.rstrip(b"\x00")
            if compression_literal == "":
                self._lcs = 0
                offset1 = offset + 32
                if self.U64:
                    offset1 += 4
                buf = self._read_bytes(offset1, self._lcp)
                self.creator_proc = buf[0:self._lcp].decode()
            elif compression_literal == _rle_compression:
                offset1 = offset + 40
                if self.U64:
                    offset1 += 4
                buf = self._read_bytes(offset1, self._lcp)
                self.creator_proc = buf[0:self._lcp].decode()
            elif self._lcs > 0:
                self._lcp = 0
                offset1 = offset + 16
                if self.U64:
                    offset1 += 4
                buf = self._read_bytes(offset1, self._lcs)
                self.creator_proc = buf[0:self._lcp].decode()

    def _process_columnname_subheader(self, offset, length):
        int_len = self._int_length
        offset += int_len
        column_name_pointers_count = (length - 2 * int_len - 12) // 8
        for i in range(column_name_pointers_count):
            text_subheader = offset + _column_name_pointer_length * \
                (i + 1) + _column_name_text_subheader_offset
            col_name_offset = offset + _column_name_pointer_length * \
                (i + 1) + _column_name_offset_offset
            col_name_length = offset + _column_name_pointer_length * \
                (i + 1) + _column_name_length_offset

            idx = self._read_int(
                text_subheader, _column_name_text_subheader_length)
            col_offset = self._read_int(
                col_name_offset, _column_name_offset_length)
            col_len = self._read_int(
                col_name_length, _column_name_length_length)

            name_str = self.column_names_strings[idx]
            self.column_names.append(name_str[col_offset:col_offset + col_len])

    def _process_columnattributes_subheader(self, offset, length):
        int_len = self._int_length
        column_attributes_vectors_count = (
            length - 2 * int_len - 12) // (int_len + 8)
        self.column_types = np.empty(
            column_attributes_vectors_count, dtype=np.dtype('S1'))
        for i in range(column_attributes_vectors_count):
            col_data_offset = (offset + int_len +
                               _column_data_offset_offset + i * (int_len + 8))
            col_data_len = (offset + 2 * int_len +
                            _column_data_length_offset + i * (int_len + 8))
            col_types = (offset + 2 * int_len +
                         _column_type_offset + i * (int_len + 8))

            self._column_data_offsets.append(
                self._read_int(col_data_offset, int_len))

            x = self._read_int(col_data_len, _column_data_length_length)
            self._column_data_lengths.append(x)

            x = self._read_int(col_types, _column_type_length)
            if x == 1:
                self.column_types[i] = b'd'
            else:
                self.column_types[i] = b's'

    def _process_columnlist_subheader(self, offset, length):
        # unknown purpose
        pass

    def _process_format_subheader(self, offset, length):
        int_len = self._int_length
        text_subheader_format = offset + \
            _column_format_text_subheader_index_offset + 3 * int_len
        col_format_offset = offset + _column_format_offset_offset + 3 * int_len
        col_format_len = offset + _column_format_length_offset + 3 * int_len
        text_subheader_label = offset + \
            _column_label_text_subheader_index_offset + 3 * int_len
        col_label_offset = offset + _column_label_offset_offset + 3 * int_len
        col_label_len = offset + _column_label_length_offset + 3 * int_len

        x = self._read_int(text_subheader_format,
                           _column_format_text_subheader_index_length)
        format_idx = min(x, len(self.column_names_strings) - 1)

        format_start = self._read_int(
            col_format_offset, _column_format_offset_length)
        format_len = self._read_int(
            col_format_len, _column_format_length_length)

        label_idx = self._read_int(
            text_subheader_label, _column_label_text_subheader_index_length)
        label_idx = min(label_idx, len(self.column_names_strings) - 1)

        label_start = self._read_int(
            col_label_offset, _column_label_offset_length)
        label_len = self._read_int(col_label_len, _column_label_length_length)

        label_names = self.column_names_strings[label_idx]
        column_label = label_names[label_start: label_start + label_len]
        format_names = self.column_names_strings[format_idx]
        column_format = format_names[format_start: format_start + format_len]
        current_column_number = len(self.columns)

        col = _column()
        col.col_id = current_column_number
        col.name = self.column_names[current_column_number]
        col.label = column_label
        col.format = column_format
        col.ctype = self.column_types[current_column_number]
        col.length = self._column_data_lengths[current_column_number]

        self.column_formats.append(column_format)
        self.columns.append(col)

    def read(self, nrows=None):

        if (nrows is None) and (self.chunksize is not None):
            nrows = self.chunksize
        elif nrows is None:
            nrows = self.row_count

        if self._current_row_in_file_index >= self.row_count:
            return None

        nd = (self.column_types == b'd').sum()
        ns = (self.column_types == b's').sum()

        self._string_chunk = np.empty((ns, nrows), dtype=np.object)
        self._byte_chunk = np.empty((nd, 8 * nrows), dtype=np.uint8)

        self._current_row_in_chunk_index = 0
        for i in range(nrows):
            done = self._readline()
            if done:
                break

        rslt = self._chunk_to_dataframe()
        if self.index is not None:
            rslt = rslt.set_index(self.index)

        return rslt

    def _readline(self):

        bit_offset = self._page_bit_offset
        subheader_pointer_length = self._subheader_pointer_length

        # If there is no page, go to the end of the header and read a page.
        if self._cached_page is None:
            self._path_or_buf.seek(self.header_length)
            done = self._read_next_page()
            if done:
                return True

        # Loop until a data row is read
        while True:
            if self._current_page_type == _page_meta_type:
                flag = (self._current_row_on_page_index >=
                        len(self._current_page_data_subheader_pointers))
                if flag:
                    done = self._read_next_page()
                    if done:
                        return True
                    self._current_row_on_page_index = 0
                    continue
                current_subheader_pointer = (
                    self._current_page_data_subheader_pointers[
                        self._current_row_on_page_index])
                process_byte_array_with_data(self,
                                             current_subheader_pointer.offset,
                                             current_subheader_pointer.length,
                                             self._byte_chunk,
                                             self._string_chunk)
                return False
            elif self._current_page_type in _page_mix_types:
                align_correction = (bit_offset + _subheader_pointers_offset +
                                    self._current_page_subheaders_count *
                                    subheader_pointer_length)
                align_correction = align_correction % 8
                offset = bit_offset + align_correction
                offset += _subheader_pointers_offset
                offset += (self._current_page_subheaders_count *
                           subheader_pointer_length)
                offset += self._current_row_on_page_index * self.row_length
                process_byte_array_with_data(self, offset, self.row_length,
                                             self._byte_chunk,
                                             self._string_chunk)
                mn = min(self.row_count, self._mix_page_row_count)
                if self._current_row_on_page_index == mn:
                    done = self._read_next_page()
                    if done:
                        return True
                    self._current_row_on_page_index = 0
                return False
            elif self._current_page_type == _page_data_type:
                process_byte_array_with_data(self,
                                             bit_offset +
                                             _subheader_pointers_offset +
                                             self._current_row_on_page_index *
                                             self.row_length,
                                             self.row_length, self._byte_chunk,
                                             self._string_chunk)
                flag = (self._current_row_on_page_index ==
                        self._current_page_block_count)
                if flag:
                    done = self._read_next_page()
                    if done:
                        return True
                    self._current_row_on_page_index = 0
                return False
            else:
                raise ValueError("unknown page type: %s",
                                 self._current_page_type)

    def _read_next_page(self):
        self._current_page_data_subheader_pointers = []
        self._cached_page = self._path_or_buf.read(self._page_length)
        if len(self._cached_page) <= 0:
            return True
        elif len(self._cached_page) != self._page_length:
            msg = ("failed to read complete page from file "
                   "(read {:d} of {:d} bytes)")
            raise ValueError(msg.format(len(self._cached_page),
                                        self._page_length))

        self._read_page_header()
        if self._current_page_type == _page_meta_type:
            self._process_page_metadata()
        pt = [_page_meta_type, _page_data_type] + [_page_mix_types]
        if self._current_page_type not in pt:
            return self._read_next_page()

        return False

    def _decompress(self, row_length, page):
        page = np.frombuffer(page, dtype=np.uint8)
        if self.compression == _rle_compression:
            return _rle_decompress(row_length, page)
        elif self.compression == _rdc_compression:
            return _rdc_decompress(row_length, page)
        else:
            raise ValueError("unknown SAS compression method: %s" %
                             self.compression)

    def _chunk_to_dataframe(self):

        n = self._current_row_in_chunk_index
        m = self._current_row_in_file_index
        ix = range(m - n, m)
        rslt = pd.DataFrame(index=ix)

        js, jb = 0, 0
        for j in range(self.column_count):

            name = self.column_names[j]

            if self.column_types[j] == b'd':
                rslt[name] = self._byte_chunk[jb, :].view(
                    dtype=self.byte_order + 'd')
                rslt[name] = np.asarray(rslt[name], dtype=np.float64)
                if self.convert_dates and (self.column_formats[j] == "MMDDYY"):
                    epoch = pd.datetime(1960, 1, 1)
                    rslt[name] = epoch + pd.to_timedelta(rslt[name], unit='d')
                jb += 1
            elif self.column_types[j] == b's':
                rslt[name] = self._string_chunk[js, :]
                rslt[name] = rslt[name].apply(lambda x: x.rstrip(b'\x00 '))
                if self.encoding is not None:
                    rslt[name] = rslt[name].apply(
                        lambda x: x.decode(encoding=self.encoding))
                if self.blank_missing:
                    ii = rslt[name].str.len() == 0
                    rslt.loc[ii, name] = np.nan
                js += 1
            else:
                raise ValueError("unknown column type %s" %
                                 self.column_types[j])

        return rslt
