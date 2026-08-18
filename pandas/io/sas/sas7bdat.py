"""
Read SAS7BDAT files

Based on code written by Jared Hobbs:
  https://bitbucket.org/jaredhobbs/sas7bdat

See also:
  https://github.com/BioStatMatt/sas7bdat

Partial documentation of the file format:
  https://cran.r-project.org/package=sas7bdat/vignettes/sas7bdat.pdf

Reference for binary data compression:
  http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/1992/9210/ross/ross.htm
"""

from __future__ import annotations

import codecs
from datetime import datetime
import functools
import sys
from typing import TYPE_CHECKING
import warnings

import numpy as np

from pandas._config import using_string_dtype

from pandas._libs import lib
from pandas._libs.byteswap import (
    read_double_with_byteswap,
    read_float_with_byteswap,
    read_uint16_with_byteswap,
    read_uint32_with_byteswap,
    read_uint64_with_byteswap,
)
from pandas._libs.sas import (
    Parser,
    collect_page_subheaders,
)
from pandas._libs.tslibs.conversion import cast_from_unit_vectorized
from pandas.compat import HAS_PYARROW
from pandas.errors import (
    EmptyDataError,
    Pandas4Warning,
)
from pandas.util._exceptions import find_stack_level

import pandas as pd
from pandas import (
    DataFrame,
    Timestamp,
)
from pandas.core.arrays.string_ import StringDtype
from pandas.core.arrays.string_arrow import ArrowStringArray

from pandas.io.common import get_handle
import pandas.io.sas.sas_constants as const
from pandas.io.sas.sasreader import SASReader

if TYPE_CHECKING:
    from pandas._typing import (
        CompressionOptions,
        FilePath,
        ReadBuffer,
    )


_unix_origin = Timestamp("1970-01-01")
_sas_origin = Timestamp("1960-01-01")


@functools.cache
def _utf8_translation_table(
    encoding: str,
) -> tuple[np.ndarray, np.ndarray, bool] | None:
    """
    Build a byte -> utf-8 lookup table for a single-byte `encoding`.

    Returns ``(table, lengths, ascii_identity)``, where ``table`` is a flat
    ``(256 * width,)`` array whose ``b``-th row holds the utf-8 encoding of
    source byte ``b``, ``lengths[b]`` is that encoding's length (or
    ``const.undefined_byte`` if `encoding` does not define ``b``), and
    ``ascii_identity`` says whether every byte below 0x80 maps to itself, which
    lets the parser memcpy ascii runs instead of translating byte by byte.

    Returns None if `encoding` is not a simple single-byte encoding — the
    multi-byte CJK encodings SAS can declare (shift_jis, big5, cp936, ...) and
    utf-8 itself, none of which have a per-byte translation.

    Cached because the probing below costs far more than a chunk's parse; the
    tables are tiny and read-only.
    """
    # A single-byte encoding maps each of the 256 bytes to exactly one
    # character, so the whole range round-trips to 256 replacement-or-real
    # characters. Multi-byte encodings consume several bytes per character and
    # come up short.
    try:
        if len(bytes(range(256)).decode(encoding, errors="replace")) != 256:
            return None
    except (LookupError, UnicodeError):
        return None

    encoded: list[bytes | None] = []
    undefined: list[int] = []
    for value in range(256):
        try:
            char = bytes([value]).decode(encoding)
        except UnicodeDecodeError:
            encoded.append(None)
            undefined.append(value)
            continue
        if len(char) != 1:
            return None
        encoded.append(char.encode("utf-8"))

    # A byte that is genuinely *undefined* fails to decode in every context,
    # whereas a multi-byte *lead* byte decodes fine once a continuation byte
    # follows it. Without this probe utf-8 itself would look single-byte (its
    # lead bytes are undefined in isolation) and every non-ascii cell would
    # wrongly raise.
    for value in undefined:
        for follower in range(256):
            try:
                bytes([value, follower]).decode(encoding)
            except UnicodeDecodeError:
                continue
            return None

    # A byte's meaning must also not depend on the bytes that follow it. The
    # probes above miss escape-style codecs: raw_unicode_escape maps every byte
    # one to one in isolation, yet decodes b"\\u0041" as "A".
    probe = b"\\u0041"
    if probe.decode(encoding, "replace") != "".join(
        bytes([value]).decode(encoding, "replace") for value in probe
    ):
        return None

    width = max((len(utf8) for utf8 in encoded if utf8 is not None), default=1)
    table = np.zeros(256 * width, dtype=np.uint8)
    lengths = np.full(256, const.undefined_byte, dtype=np.uint8)
    for value, utf8 in enumerate(encoded):
        if utf8 is not None:
            lengths[value] = len(utf8)
            table[value * width : value * width + len(utf8)] = bytearray(utf8)
    ascii_identity = all(encoded[value] == bytes([value]) for value in range(0x80))
    # Every reader gets these same arrays, so no one may write to them.
    table.flags.writeable = False
    lengths.flags.writeable = False
    return table, lengths, ascii_identity


# SAS uses a modified Gregorian calendar where years divisible by 4000 are
# not leap years (unlike proleptic Gregorian). These are the SAS day counts
# for the first day affected by each 4000-year boundary.
# See https://communities.sas.com/t5/SAS-Programming/Leap-Years-divisible-by-4000/td-p/663467
_SAS_MARCH1_4000 = 745154  # SAS day count for March 1, 4000
_SAS_MARCH1_8000 = 2206123  # SAS day count for March 1, 8000


def _sas_to_gregorian_correction(values: np.ndarray, unit: str) -> np.ndarray:
    """
    Compute the additive correction (in `unit`) to convert SAS day/second counts
    to proleptic Gregorian day/second counts.

    SAS omits Feb 29 for years divisible by 4000 (unlike proleptic Gregorian);
    this adds back the missing days. `unit` must be "d" (days) or "s" (seconds).
    """
    scale = 86400 if unit == "s" else 1
    thresholds = np.array([_SAS_MARCH1_4000, _SAS_MARCH1_8000], dtype=np.int64) * scale
    correction = np.zeros(len(values), dtype=np.float64)
    valid = ~np.isnan(values)
    for threshold in thresholds:
        correction[valid] += (values[valid] >= threshold).astype(np.float64) * scale
    return correction


def _convert_datetimes(sas_datetimes: pd.Series, unit: str) -> pd.Series:
    """
    Convert to Timestamp if possible, otherwise to datetime.datetime.
    SAS float64 lacks precision for more than ms resolution so the fit
    to datetime.datetime is ok.

    Parameters
    ----------
    sas_datetimes : {Series, Sequence[float]}
       Dates or datetimes in SAS
    unit : {'d', 's'}
       "d" if the floats represent dates, "s" for datetimes

    Returns
    -------
    Series
       Series of datetime64 dtype or datetime.datetime.
    """
    td = (_sas_origin - _unix_origin).as_unit("s")
    if unit == "s":
        corrected = sas_datetimes._values + _sas_to_gregorian_correction(
            sas_datetimes._values, unit="s"
        )
        millis = cast_from_unit_vectorized(corrected, unit="s", out_unit="ms")
        dt64ms = millis.view("M8[ms]") + td
        return pd.Series(dt64ms, index=sas_datetimes.index, copy=False)
    else:
        corrected = sas_datetimes._values + _sas_to_gregorian_correction(
            sas_datetimes._values, unit="d"
        )
        vals = np.array(corrected, dtype="M8[D]") + td
        return pd.Series(vals, dtype="M8[s]", index=sas_datetimes.index, copy=False)


class _Column:
    col_id: int
    name: str | bytes
    label: str | bytes
    format: str | bytes
    ctype: bytes
    length: int

    def __init__(
        self,
        col_id: int,
        # These can be bytes when convert_header_text is False
        name: str | bytes,
        label: str | bytes,
        format: str | bytes,
        ctype: bytes,
        length: int,
    ) -> None:
        self.col_id = col_id
        self.name = name
        self.label = label
        self.format = format
        self.ctype = ctype
        self.length = length


# SAS7BDAT represents a SAS data file in SAS7BDAT format.
class SAS7BDATReader(SASReader):
    """
    Read SAS files in SAS7BDAT format.

    Parameters
    ----------
    path_or_buf : path name or buffer
        Name of SAS file or file-like object pointing to SAS file
        contents.
    index : column identifier, defaults to None
        Column to use as index.
    convert_dates : bool, defaults to True
        Attempt to convert dates to Pandas datetime values.  Note that
        some rarely used SAS date formats may be unsupported.
    blank_missing : bool, defaults to True
        Convert empty strings to missing values (SAS uses blanks to
        indicate missing character variables).
    chunksize : int, defaults to None
        Return SAS7BDATReader object for iterations, returns chunks
        with given number of lines.
    encoding : str, 'infer', defaults to None
        String encoding acc. to Python standard encodings,
        encoding='infer' decodes using the encoding recorded in the file
        header, falling back to latin-1 if the header does not name one
        pandas recognizes; encoding=None will leave the data in binary format.

        .. deprecated:: 3.1.0
            The default will change from ``None`` to ``'infer'``.
    convert_text : bool, defaults to True
        If False, text variables are left as raw bytes.
    convert_header_text : bool, defaults to True
        If False, header text, including column names, are left as raw
        bytes.
    """

    _int_length: int
    _cached_page: bytes | None
    encoding: str | None

    def __init__(
        self,
        path_or_buf: FilePath | ReadBuffer[bytes],
        index=None,
        convert_dates: bool = True,
        blank_missing: bool = True,
        chunksize: int | None = None,
        encoding: str | lib.NoDefault | None = lib.no_default,
        convert_text: bool = True,
        convert_header_text: bool = True,
        compression: CompressionOptions = "infer",
    ) -> None:
        self._encoding_specified = encoding is not lib.no_default
        self._non_ascii_header_text = False

        self.index = index
        self.convert_dates = convert_dates
        self.blank_missing = blank_missing
        self.chunksize = chunksize
        self.encoding = None if encoding is lib.no_default else encoding
        self.convert_text = convert_text
        self.convert_header_text = convert_header_text

        self.default_encoding = "latin-1"
        self.compression = b""
        self.column_names_raw: list[bytes] = []
        self.column_names: list[str | bytes] = []
        self.column_formats: list[str | bytes] = []
        self.columns: list[_Column] = []

        # (offset, length) of each data subheader on the current page, filled
        # in by collect_page_subheaders. Sized in _get_properties, once the
        # page length is known.
        self._data_subheader_offsets = np.empty(0, dtype=np.int64)
        self._data_subheader_lengths = np.empty(0, dtype=np.int64)
        self._data_subheader_count = 0
        self._cached_page = None
        self._column_data_lengths: list[int] = []
        self._column_data_offsets: list[int] = []
        self._column_types: list[bytes] = []

        self._current_row_in_file_index = 0
        self._current_row_on_page_index = 0
        self._current_row_in_file_index = 0

        self.handles = get_handle(
            path_or_buf, "rb", is_text=False, compression=compression
        )

        self._path_or_buf = self.handles.handle

        # Same order as const.SASIndex
        self._subheader_processors = [
            self._process_rowsize_subheader,
            self._process_columnsize_subheader,
            self._process_subheader_counts,
            self._process_columntext_subheader,
            self._process_columnname_subheader,
            self._process_columnattributes_subheader,
            self._process_format_subheader,
            self._process_columnlist_subheader,
            None,  # Data
        ]

        try:
            self._get_properties()
            self._parse_metadata()
            self._validate_column_data_ranges()
        except Exception:
            self.close()
            raise
        self._metadata_at_open = self._metadata_signature()

    def _metadata_signature(self) -> tuple:
        """
        Everything the parser and ``_chunk_to_dataframe`` read out of the file's
        metadata, in one comparable value.

        The column lists are only ever appended to, so their lengths stand in
        for their contents and this stays O(1) rather than O(ncols). The
        compression is in here for completeness rather than because a file can
        change it: only the first column-text subheader is read for it, and a
        file with none of those cannot be opened at all.
        """
        return (
            self.row_length,
            self.row_count,
            self._mix_page_row_count,
            self.column_count,
            self.compression,
            len(self._column_data_offsets),
            len(self.column_names),
        )

    def _validate_column_data_ranges(self) -> None:
        # The column offsets, lengths and types come straight out of the file's
        # column-attributes subheader, and every row is then read at those
        # positions out of a buffer only row_length bytes long. Validate them
        # once here rather than per row: a corrupt or hostile file would
        # otherwise read past the end of that buffer. These are Python ints, so
        # unlike the int64 the parser holds them in, the sum cannot overflow.
        for index, (offset, length, ctype) in enumerate(
            zip(
                self._column_data_offsets,
                self._column_data_lengths,
                self._column_types,
                strict=True,
            )
        ):
            # A numeric column is widened into 8 bytes of the parser's output
            # buffer, so a longer one would also write past the front of it.
            if ctype == b"d" and length > 8:
                raise ValueError(
                    f"Column {index} is numeric but declares {length} bytes of "
                    f"data, more than the 8 a SAS number can occupy; the file "
                    f"is corrupt"
                )
            if offset + length > self.row_length:
                raise ValueError(
                    f"Column {index} spans bytes {offset}-{offset + length} of a "
                    f"{self.row_length}-byte row; the file is corrupt"
                )

    def column_data_lengths(self) -> np.ndarray:
        """Return a numpy int64 array of the column data lengths"""
        return np.asarray(self._column_data_lengths, dtype=np.int64)

    def column_data_offsets(self) -> np.ndarray:
        """Return a numpy int64 array of the column offsets"""
        return np.asarray(self._column_data_offsets, dtype=np.int64)

    def column_types(self) -> np.ndarray:
        """
        Returns a numpy character array of the column types:
           s (string) or d (double)
        """
        return np.asarray(self._column_types, dtype=np.dtype("S1"))

    def close(self) -> None:
        self.handles.close()

    def _get_properties(self) -> None:
        # Check magic number
        self._path_or_buf.seek(0)
        self._cached_page = self._path_or_buf.read(288)
        if self._cached_page[0 : len(const.magic)] != const.magic:
            raise ValueError("magic number mismatch (not a SAS file?)")

        # Get alignment information
        buf = self._read_bytes(const.align_1_offset, const.align_1_length)
        if buf == const.u64_byte_checker_value:
            self.U64 = True
            self._int_length = 8
            self._page_bit_offset = const.page_bit_offset_x64
            self._subheader_pointer_length = const.subheader_pointer_length_x64
        else:
            self.U64 = False
            self._page_bit_offset = const.page_bit_offset_x86
            self._subheader_pointer_length = const.subheader_pointer_length_x86
            self._int_length = 4
        buf = self._read_bytes(const.align_2_offset, const.align_2_length)
        if buf == const.align_1_checker_value:
            align1 = const.align_2_value
        else:
            align1 = 0

        # Get endianness information
        buf = self._read_bytes(const.endianness_offset, const.endianness_length)
        if buf == b"\x01":
            self.byte_order = "<"
            self.need_byteswap = sys.byteorder == "big"
        else:
            self.byte_order = ">"
            self.need_byteswap = sys.byteorder == "little"

        # Get encoding information
        buf = self._read_bytes(const.encoding_offset, const.encoding_length)[0]
        if buf in const.encoding_names:
            self.inferred_encoding = const.encoding_names[buf]
        else:
            self.inferred_encoding = f"unknown (code={buf})"
        if self.encoding == "infer":
            # GH#66470 fall back rather than trying to decode with "infer".
            # default_encoding is latin-1, which cannot raise, so encoding=None
            # (opting out of decoding) keeps decoding header text infallibly.
            self.encoding = const.encoding_names.get(buf, self.default_encoding)

        # Timestamp is epoch 01/01/1960
        epoch = datetime(1960, 1, 1)
        x = self._read_float(
            const.date_created_offset + align1, const.date_created_length
        )
        self.date_created = epoch + pd.to_timedelta(x, unit="s")
        x = self._read_float(
            const.date_modified_offset + align1, const.date_modified_length
        )
        self.date_modified = epoch + pd.to_timedelta(x, unit="s")

        self.header_length = self._read_uint(
            const.header_size_offset + align1, const.header_size_length
        )

        # Read the rest of the header into cached_page.
        buf = self._path_or_buf.read(self.header_length - 288)
        self._cached_page += buf
        if len(self._cached_page) != self.header_length:
            raise ValueError("The SAS7BDAT file appears to be truncated.")

        self._page_length = self._read_uint(
            const.page_size_offset + align1, const.page_size_length
        )

        # A page cannot hold more subheader pointers than its own length divided
        # by the pointer size, nor more than its subheader-count field can express,
        # so those bound the data-subheader arrays. The second bound matters
        # because _page_length is an unvalidated header field.
        # Parser caches memoryviews of these, so never rebind them after that.
        max_subheaders = min(
            self._page_length // self._subheader_pointer_length + 1,
            1 << (8 * const.subheader_count_length),
        )
        self._data_subheader_offsets = np.empty(max_subheaders, dtype=np.int64)
        self._data_subheader_lengths = np.empty(max_subheaders, dtype=np.int64)

    def __next__(self) -> DataFrame:
        da = self.read(nrows=self.chunksize or 1)
        if da.empty:
            self.close()
            raise StopIteration
        return da

    # Read a single float of the given width (4 or 8).
    def _read_float(self, offset: int, width: int) -> float:
        assert self._cached_page is not None
        if width == 4:
            return read_float_with_byteswap(
                self._cached_page, offset, self.need_byteswap
            )
        elif width == 8:
            return read_double_with_byteswap(
                self._cached_page, offset, self.need_byteswap
            )
        else:
            self.close()
            raise ValueError("invalid float width")

    # Read a single unsigned integer of the given width (1, 2, 4 or 8).
    def _read_uint(self, offset: int, width: int) -> int:
        assert self._cached_page is not None
        if width == 1:
            return self._read_bytes(offset, 1)[0]
        elif width == 2:
            return read_uint16_with_byteswap(
                self._cached_page, offset, self.need_byteswap
            )
        elif width == 4:
            return read_uint32_with_byteswap(
                self._cached_page, offset, self.need_byteswap
            )
        elif width == 8:
            return read_uint64_with_byteswap(
                self._cached_page, offset, self.need_byteswap
            )
        else:
            self.close()
            raise ValueError("invalid int width")

    def _read_bytes(self, offset: int, length: int):
        assert self._cached_page is not None
        if offset + length > len(self._cached_page):
            self.close()
            raise ValueError("The cached page is too small.")
        return self._cached_page[offset : offset + length]

    def _parse_metadata(self) -> None:
        done = False
        while not done:
            self._cached_page = self._path_or_buf.read(self._page_length)
            if len(self._cached_page) <= 0:
                break
            if len(self._cached_page) != self._page_length:
                raise ValueError("Failed to read a meta data page from the SAS file.")
            done = self._process_page_meta()

        self._column_convert_types: list[str | None] = []
        for j in range(self.column_count):
            if self._column_types[j] == b"d" and self.convert_dates:
                fmt = self.column_formats[j]
                if fmt in const.sas_date_formats:
                    self._column_convert_types.append("d")
                elif fmt in const.sas_datetime_formats:
                    self._column_convert_types.append("s")
                else:
                    self._column_convert_types.append(None)
            else:
                self._column_convert_types.append(None)

        # The default also governs header text, which is decoded either way, so
        #  a numeric-only file with non-ascii column names changes too
        if not self._encoding_specified and (
            (self.convert_text and b"s" in self._column_types)
            or self._non_ascii_header_text
        ):
            warnings.warn(
                "The default value of 'encoding' in read_sas is deprecated. In a "
                "future version the default will change from None to 'infer', and "
                "text will be decoded using the encoding recorded in the file "
                "rather than returned as bytes. Pass encoding='infer' to adopt "
                "the future behavior, or encoding=None to keep the current one.",
                Pandas4Warning,
                stacklevel=find_stack_level(),
            )

    def _process_page_meta(self) -> bool:
        self._read_page_header()
        pt = [*const.page_meta_types, const.page_amd_type, const.page_mix_type]
        if self._current_page_type in pt:
            self._process_page_metadata()
        is_data_page = self._current_page_type == const.page_data_type
        is_mix_page = self._current_page_type == const.page_mix_type
        return bool(is_data_page or is_mix_page or self._data_subheader_count > 0)

    def _read_page_header(self) -> None:
        bit_offset = self._page_bit_offset
        tx = const.page_type_offset + bit_offset
        self._current_page_type = (
            self._read_uint(tx, const.page_type_length) & const.page_type_mask2
        )
        tx = const.block_count_offset + bit_offset
        self._current_page_block_count = self._read_uint(tx, const.block_count_length)
        tx = const.subheader_count_offset + bit_offset
        self._current_page_subheaders_count = self._read_uint(
            tx, const.subheader_count_length
        )

    def _process_page_metadata(self) -> None:
        # Walks the page's subheader pointers in Cython, collecting data
        # subheaders into self._data_subheader_offsets/_lengths and dispatching
        # the rest back to self._subheader_processors.
        try:
            self._data_subheader_count = collect_page_subheaders(self)
        except Exception:
            self.close()
            raise

    def _process_rowsize_subheader(self, offset: int, length: int) -> None:
        int_len = self._int_length
        lcs_offset = offset
        lcp_offset = offset
        if self.U64:
            lcs_offset += 682
            lcp_offset += 706
        else:
            lcs_offset += 354
            lcp_offset += 378

        self.row_length = self._read_uint(
            offset + const.row_length_offset_multiplier * int_len,
            int_len,
        )
        if self.row_length > self._page_length:
            raise ValueError(
                f"row_length ({self.row_length}) exceeds the page size "
                f"({self._page_length}); the file is corrupt"
            )
        self.row_count = self._read_uint(
            offset + const.row_count_offset_multiplier * int_len,
            int_len,
        )
        self.col_count_p1 = self._read_uint(
            offset + const.col_count_p1_multiplier * int_len, int_len
        )
        self.col_count_p2 = self._read_uint(
            offset + const.col_count_p2_multiplier * int_len, int_len
        )
        mx = const.row_count_on_mix_page_offset_multiplier * int_len
        self._mix_page_row_count = self._read_uint(offset + mx, int_len)
        self._lcs = self._read_uint(lcs_offset, 2)
        self._lcp = self._read_uint(lcp_offset, 2)

    def _process_columnsize_subheader(self, offset: int, length: int) -> None:
        int_len = self._int_length
        offset += int_len
        self.column_count = self._read_uint(offset, int_len)
        if self.col_count_p1 + self.col_count_p2 != self.column_count:
            print(
                f"Warning: column count mismatch ({self.col_count_p1} + "
                f"{self.col_count_p2} != {self.column_count})\n"
            )

    # Unknown purpose
    def _process_subheader_counts(self, offset: int, length: int) -> None:
        pass

    def _process_columntext_subheader(self, offset: int, length: int) -> None:
        offset += self._int_length
        text_block_size = self._read_uint(offset, const.text_block_size_length)

        buf = self._read_bytes(offset, text_block_size)
        cname_raw = buf[0:text_block_size].rstrip(b"\x00 ")
        self.column_names_raw.append(cname_raw)

        if len(self.column_names_raw) == 1:
            compression_literal = b""
            for cl in const.compression_literals:
                if cl in cname_raw:
                    compression_literal = cl
            self.compression = compression_literal
            offset -= self._int_length

            offset1 = offset + 16
            if self.U64:
                offset1 += 4

            buf = self._read_bytes(offset1, self._lcp)
            compression_literal = buf.rstrip(b"\x00")
            if compression_literal == b"":
                self._lcs = 0
                offset1 = offset + 32
                if self.U64:
                    offset1 += 4
                buf = self._read_bytes(offset1, self._lcp)
                self.creator_proc = buf[0 : self._lcp]
            elif compression_literal == const.rle_compression:
                offset1 = offset + 40
                if self.U64:
                    offset1 += 4
                buf = self._read_bytes(offset1, self._lcp)
                self.creator_proc = buf[0 : self._lcp]
            elif self._lcs > 0:
                self._lcp = 0
                offset1 = offset + 16
                if self.U64:
                    offset1 += 4
                buf = self._read_bytes(offset1, self._lcs)
                self.creator_proc = buf[0 : self._lcp]
            if hasattr(self, "creator_proc"):
                self.creator_proc = self._convert_header_text(self.creator_proc)  # pyright: ignore[reportArgumentType]

    def _process_columnname_subheader(self, offset: int, length: int) -> None:
        int_len = self._int_length
        offset += int_len
        column_name_pointers_count = (length - 2 * int_len - 12) // 8
        for i in range(column_name_pointers_count):
            text_subheader = (
                offset
                + const.column_name_pointer_length * (i + 1)
                + const.column_name_text_subheader_offset
            )
            col_name_offset = (
                offset
                + const.column_name_pointer_length * (i + 1)
                + const.column_name_offset_offset
            )
            col_name_length = (
                offset
                + const.column_name_pointer_length * (i + 1)
                + const.column_name_length_offset
            )

            idx = self._read_uint(
                text_subheader, const.column_name_text_subheader_length
            )
            col_offset = self._read_uint(
                col_name_offset, const.column_name_offset_length
            )
            col_len = self._read_uint(col_name_length, const.column_name_length_length)

            name_raw = self.column_names_raw[idx]
            cname = name_raw[col_offset : col_offset + col_len]
            self.column_names.append(self._convert_header_text(cname))

    def _process_columnattributes_subheader(self, offset: int, length: int) -> None:
        int_len = self._int_length
        column_attributes_vectors_count = (length - 2 * int_len - 12) // (int_len + 8)
        for i in range(column_attributes_vectors_count):
            col_data_offset = (
                offset + int_len + const.column_data_offset_offset + i * (int_len + 8)
            )
            col_data_len = (
                offset
                + 2 * int_len
                + const.column_data_length_offset
                + i * (int_len + 8)
            )
            col_types = (
                offset + 2 * int_len + const.column_type_offset + i * (int_len + 8)
            )

            x = self._read_uint(col_data_offset, int_len)
            self._column_data_offsets.append(x)

            x = self._read_uint(col_data_len, const.column_data_length_length)
            self._column_data_lengths.append(x)

            x = self._read_uint(col_types, const.column_type_length)
            self._column_types.append(b"d" if x == 1 else b"s")

    def _process_columnlist_subheader(self, offset: int, length: int) -> None:
        # unknown purpose
        pass

    def _process_format_subheader(self, offset: int, length: int) -> None:
        int_len = self._int_length
        text_subheader_format = (
            offset + const.column_format_text_subheader_index_offset + 3 * int_len
        )
        col_format_offset = offset + const.column_format_offset_offset + 3 * int_len
        col_format_len = offset + const.column_format_length_offset + 3 * int_len
        text_subheader_label = (
            offset + const.column_label_text_subheader_index_offset + 3 * int_len
        )
        col_label_offset = offset + const.column_label_offset_offset + 3 * int_len
        col_label_len = offset + const.column_label_length_offset + 3 * int_len

        x = self._read_uint(
            text_subheader_format, const.column_format_text_subheader_index_length
        )
        format_idx = min(x, len(self.column_names_raw) - 1)

        format_start = self._read_uint(
            col_format_offset, const.column_format_offset_length
        )
        format_len = self._read_uint(col_format_len, const.column_format_length_length)

        label_idx = self._read_uint(
            text_subheader_label, const.column_label_text_subheader_index_length
        )
        label_idx = min(label_idx, len(self.column_names_raw) - 1)

        label_start = self._read_uint(
            col_label_offset, const.column_label_offset_length
        )
        label_len = self._read_uint(col_label_len, const.column_label_length_length)

        label_names = self.column_names_raw[label_idx]
        column_label = self._convert_header_text(
            label_names[label_start : label_start + label_len]
        )
        format_names = self.column_names_raw[format_idx]
        column_format = self._convert_header_text(
            format_names[format_start : format_start + format_len]
        )
        current_column_number = len(self.columns)

        col = _Column(
            current_column_number,
            self.column_names[current_column_number],
            column_label,
            column_format,
            self._column_types[current_column_number],
            self._column_data_lengths[current_column_number],
        )

        self.column_formats.append(column_format)
        self.columns.append(col)

    def read(self, nrows: int | None = None) -> DataFrame:
        if (nrows is None) and (self.chunksize is not None):
            nrows = self.chunksize
        elif nrows is None:
            nrows = self.row_count

        if len(self._column_types) == 0:
            self.close()
            raise EmptyDataError("No columns to parse from file")

        if nrows > 0 and self._current_row_in_file_index >= self.row_count:
            return DataFrame()

        nrows = min(nrows, self.row_count - self._current_row_in_file_index)

        nd = self._column_types.count(b"d")
        ns = self._column_types.count(b"s")

        self._byte_chunk = np.zeros((nd, 8 * nrows), dtype=np.uint8)
        self._setup_string_buffers(ns, nrows)

        self._current_row_in_chunk_index = 0
        try:
            p = Parser(self)
            p.read(nrows)
            # A metadata page can follow the data pages, and its subheaders
            # rewrite the layout the rows just read were parsed with -- and that
            # the next chunk would be parsed with. The reader cannot reshape
            # itself mid-file, so a file that redefines its own layout is corrupt
            # however the redefinition reads: too wide a row and the columns come
            # from the wrong bytes, too few columns and one comes back all-NaN,
            # too many and the parser walks off the ends of the offset and length
            # arrays, and a rewritten row count pads the result with unread rows
            # or truncates it.
            if self._metadata_signature() != self._metadata_at_open:
                raise ValueError(
                    "The file changes the layout it declared partway through; "
                    "the file is corrupt"
                )
            if self._str_mode != const.string_mode_object:
                self._str_values = p.string_values()
            # Release the growable parse buffers before building the DataFrame,
            # so that their spare capacity is not held alongside the copies above.
            del p

            rslt = self._chunk_to_dataframe()
        except Exception:
            # However this chunk failed -- a page that does not hold what it
            # claims, a layout the file redefined, a cell that cannot be
            # represented -- the reader is no longer usable, and read_sas only
            # closes it for the caller when it reads the whole file itself.
            self.close()
            raise
        if self.index is not None:
            rslt = rslt.set_index(self.index)

        return rslt

    def _string_mode(
        self, ns: int
    ) -> tuple[int, tuple[np.ndarray, np.ndarray, bool] | None]:
        """
        Pick how the Cython parser should emit string cells for this chunk.

        Returns ``(mode, table)``; see the string modes in sas_constants.py. The
        pyarrow modes require that the result actually be a pyarrow-backed str
        column, and that the source encoding translate byte by byte to utf-8.
        """
        if ns == 0 or not self.convert_text or self.encoding is None:
            # Nothing to decode: string cells stay raw bytes.
            return const.string_mode_object, None
        if not (using_string_dtype() and HAS_PYARROW):
            return const.string_mode_object, None
        if StringDtype(na_value=np.nan).storage != "pyarrow":
            # mode.string_storage is pinned to "python"; building an Arrow
            # array here would only have to be converted back.
            return const.string_mode_object, None

        try:
            canonical = codecs.lookup(self.encoding).name
        except LookupError:
            # Leave the error to the per-cell decode, which words it better.
            return const.string_mode_object, None
        if canonical == "utf-8":
            return const.string_mode_utf8, None
        table = _utf8_translation_table(canonical)
        if table is None:
            # Multi-byte encoding (shift_jis, big5, ...): no per-byte mapping.
            return const.string_mode_object, None
        return const.string_mode_table, table

    def _setup_string_buffers(self, ns: int, nrows: int) -> None:
        mode, table = self._string_mode(ns)
        self._str_mode = mode
        # The parser reads every buffer below unconditionally, so hand it empty
        # arrays for whichever output this chunk does not use.
        empty_u1 = np.empty(0, dtype=np.uint8)
        self._string_chunk = np.empty(
            (ns if mode == const.string_mode_object else 0, nrows), dtype=object
        )
        self._str_table = empty_u1
        self._str_table_len = empty_u1
        self._str_table_width = 1
        self._str_ascii_identity = False
        self._str_encoding = self.encoding
        if mode == const.string_mode_object:
            self._str_offsets = np.empty((0, 0), dtype=np.int64)
            self._str_valid = np.empty((0, 0), dtype=np.uint8)
            return

        self._str_offsets = np.zeros((ns, nrows + 1), dtype=np.int64)
        # Cells the parser never reaches stay null, matching the object path's
        # np.empty(dtype=object) filling those cells with None.
        self._str_valid = np.zeros((ns, nrows), dtype=np.uint8)
        if table is not None:
            self._str_table, self._str_table_len, self._str_ascii_identity = table
            self._str_table_width = len(self._str_table) // 256

    def _read_page_data(self):
        """Read the next page from the file. Returns True if EOF."""
        self._data_subheader_count = 0
        self._cached_page = self._path_or_buf.read(self._page_length)
        if len(self._cached_page) <= 0:
            return True
        elif len(self._cached_page) != self._page_length:
            self.close()
            msg = (
                "failed to read complete page from file (read "
                f"{len(self._cached_page):d} of {self._page_length:d} bytes)"
            )
            raise ValueError(msg)
        return False

    def _read_next_page(self):
        done = self._read_page_data()
        if done:
            return True

        self._read_page_header()
        if self._current_page_type in const.page_meta_types:
            self._process_page_metadata()

        if self._current_page_type not in [
            *const.page_meta_types,
            const.page_data_type,
            const.page_mix_type,
        ]:
            return self._read_next_page()

        return False

    def _chunk_to_dataframe(self) -> DataFrame:
        n = self._current_row_in_chunk_index
        m = self._current_row_in_file_index
        ix = range(m - n, m)
        rslt = {}

        js, jb = 0, 0
        infer_string = using_string_dtype()
        for j in range(self.column_count):
            name = self.column_names[j]

            if self._column_types[j] == b"d":
                col_arr = self._byte_chunk[jb, :].view(dtype=self.byte_order + "d")
                rslt[name] = pd.Series(col_arr, dtype=np.float64, index=ix, copy=False)
                convert_type = self._column_convert_types[j]
                if convert_type is not None:
                    rslt[name] = _convert_datetimes(rslt[name], convert_type)
                jb += 1
            elif self._column_types[j] == b"s":
                if self._str_mode != const.string_mode_object:
                    rslt[name] = pd.Series(
                        self._string_column(js), index=ix, copy=False
                    )
                else:
                    rslt[name] = pd.Series(
                        self._string_chunk[js, :], index=ix, copy=False
                    )
                    if self.convert_text and (self.encoding is not None):
                        rslt[name] = self._decode_string(rslt[name].str)
                        if infer_string:
                            rslt[name] = rslt[name].astype("str")

                js += 1
            else:
                self.close()
                raise ValueError(f"unknown column type {self._column_types[j]!r}")

        df = DataFrame(rslt, columns=self.column_names, index=ix, copy=False)
        return df

    def _string_column(self, js: int) -> ArrowStringArray:
        """
        Wrap string column `js`'s parse buffers as a pyarrow-backed str column.

        The Cython parser has already written utf-8 into ``_str_values[js]``,
        so this only has to hand the buffers to pyarrow — there is no
        object-dtype array of bytes to decode.
        """
        import pyarrow as pa

        values = self._str_values[js]
        # Copied out of the shared 2-D offsets array so that keeping one column
        # of the result alive does not pin every column's offsets.
        offsets = self._str_offsets[js].copy()
        valid = self._str_valid[js]
        nrows = len(valid)

        parsed = self._current_row_in_chunk_index
        if parsed < nrows:
            # The parser stopped early (a file whose pages run out before the
            # row count in its header). Offsets for the rows it never reached
            # are still zero, which would make them non-monotonic; flatten the
            # tail so the array is well formed and the caller's length check
            # reports the truncation instead of pyarrow reporting bad offsets.
            offsets[parsed + 1 :] = offsets[parsed]

        if valid.all():
            validity = None
        else:
            validity = pa.py_buffer(np.packbits(valid, bitorder="little"))
        pa_arr = pa.Array.from_buffers(
            pa.large_string(),
            nrows,
            [validity, pa.py_buffer(offsets), pa.py_buffer(values)],
        )
        if self._str_mode == const.string_mode_utf8:
            # from_buffers does no validation, and a file declaring utf-8 can
            # still hold invalid bytes, so check here rather than deferring the
            # failure to first access. validate() checks each value separately,
            # which matters because a multi-byte character split across two
            # fixed-width SAS fields is invalid in both cells even though the
            # concatenated values buffer is well-formed.
            try:
                pa_arr.validate(full=True)
            except pa.lib.ArrowInvalid:
                # Re-raise as the UnicodeDecodeError the per-cell decode gave.
                for row in range(len(offsets) - 1):
                    bytes(values[offsets[row] : offsets[row + 1]]).decode("utf-8")
                raise
        return ArrowStringArray(
            pa.chunked_array([pa_arr]), dtype=StringDtype(na_value=np.nan)
        )

    def _decode_string(self, b):
        return b.decode(self.encoding or self.default_encoding)

    def _convert_header_text(self, b: bytes) -> str | bytes:
        if self.convert_header_text:
            if not b.isascii():
                # latin-1 and the file's own encoding disagree only here
                self._non_ascii_header_text = True
            return self._decode_string(b)
        else:
            return b
