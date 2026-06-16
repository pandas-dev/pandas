"""
Tests for the parallel read_csv implementation (C engine, large local files).

The parallel path is enabled automatically when:
  - reading from a local uncompressed file path
  - using the C engine (the default)
  - not in iterator / chunked mode
  - file size >= _PARALLEL_READ_MIN_BYTES

We test correctness (parallel == serial) rather than performance.
"""

from __future__ import annotations

import csv
import io
import os
from typing import TYPE_CHECKING

import numpy as np
import pytest

from pandas.errors import ParserWarning

if TYPE_CHECKING:
    from pathlib import Path

from pandas import (
    DataFrame,
    option_context,
    read_csv,
)
import pandas._testing as tm

from pandas.io.parsers.base_parser import ParserBase
from pandas.io.parsers.readers import (
    _can_parallelize_csv,
    _find_chunk_byte_offsets,
    _find_data_start_offset,
    _read_csv_parallel,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_large_csv(path: Path, n_rows: int = 10_000) -> None:
    """Write a CSV with enough rows to split into multiple parallel chunks."""
    rng = np.random.default_rng(42)
    df = DataFrame(
        {
            "int_col": rng.integers(0, 10_000, size=n_rows),
            "float_col": rng.random(n_rows),
            "str_col": rng.choice(["foo", "bar", "baz", "qux"], size=n_rows),
        }
    )
    df.to_csv(path, index=False)


# ---------------------------------------------------------------------------
# _can_parallelize_csv
# ---------------------------------------------------------------------------


class TestCanParallelizeCsv:
    """Unit tests for the eligibility predicate."""

    def _kwds(self, **overrides):
        base = {
            "engine": "c",
            "engine_specified": False,
            "iterator": False,
            "chunksize": None,
            "nrows": None,
            "skiprows": None,
            "header": 0,
            "skipfooter": 0,
            "index_col": None,
            "usecols": None,
            "compression": "infer",
        }
        base.update(overrides)
        return base

    def test_eligible_large_file(self, tmp_path, monkeypatch):
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "big.csv"
        _make_large_csv(path)
        # Lower the threshold so the test file qualifies regardless of its size.
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert _can_parallelize_csv(path, self._kwds())

    def test_accepts_default_engine(self, tmp_path, monkeypatch):
        """engine=None (the default) should be treated the same as engine='c'."""
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert _can_parallelize_csv(path, self._kwds(engine=None))

    def test_rejects_file_like_object(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        with open(path, encoding="utf-8") as f:
            assert not _can_parallelize_csv(f, self._kwds())

    def test_rejects_url(self):
        assert not _can_parallelize_csv("http://example.com/data.csv", self._kwds())
        assert not _can_parallelize_csv("s3://bucket/key.csv", self._kwds())

    def test_rejects_compressed(self, tmp_path):
        path = tmp_path / "data.csv.gz"
        path.write_bytes(b"fake")
        assert not _can_parallelize_csv(path, self._kwds(compression="gzip"))
        assert not _can_parallelize_csv(path, self._kwds(compression="infer"))

    def test_rejects_iterator_mode(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(iterator=True))

    def test_rejects_chunksize(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(chunksize=1000))

    def test_rejects_nrows(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(nrows=100))

    def test_rejects_non_c_engine(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(engine="python"))
        assert not _can_parallelize_csv(path, self._kwds(engine="pyarrow"))

    def test_rejects_custom_lineterminator(self, tmp_path, monkeypatch):
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        # The chunk splitter always scans for \n; a different lineterminator
        # would misalign chunk boundaries and produce silently wrong output.
        assert not _can_parallelize_csv(path, self._kwds(lineterminator="\r"))
        assert not _can_parallelize_csv(path, self._kwds(lineterminator="|"))
        # \n and None (default) are acceptable.
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert _can_parallelize_csv(path, self._kwds(lineterminator="\n"))
        assert _can_parallelize_csv(path, self._kwds(lineterminator=None))

    @pytest.mark.parametrize(
        "encoding", ["utf-16", "UTF-16", "utf_16_le", "utf-32", "UTF_32_BE"]
    )
    def test_rejects_utf16_32_encoding(self, tmp_path, monkeypatch, encoding):
        # UTF-16/32 encode \n as a multi-byte sequence (e.g. b"\x0a\x00" in
        # UTF-16LE).  Splitting on raw \n bytes would misalign chunks and
        # silently corrupt data.
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(encoding=encoding))

    @pytest.mark.parametrize("encoding", ["latin-1", "cp1252", "shift_jis"])
    def test_rejects_non_utf8_encoding(self, tmp_path, monkeypatch, encoding):
        # Chunk workers feed raw file bytes to the C tokenizer, which decodes
        # words as UTF-8; only UTF-8-compatible encodings are byte-safe (GH#64347).
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(encoding=encoding))
        assert _can_parallelize_csv(path, self._kwds(encoding="utf-8"))
        assert _can_parallelize_csv(path, self._kwds(encoding="UTF8"))

    def test_rejects_python_engine_seps(self, tmp_path, monkeypatch):
        # Multi-char/regex seps (other than r"\s+") and sep=None force the
        # python engine inside TextFileReader (GH#64347).
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(delimiter=";;"))
        assert not _can_parallelize_csv(path, self._kwds(delimiter=None))
        assert _can_parallelize_csv(path, self._kwds(delimiter=r"\s+"))

    def test_rejects_comment(self, tmp_path, monkeypatch):
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(comment="#"))

    def test_rejects_escapechar(self, tmp_path, monkeypatch):
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(escapechar="\\"))

    def test_rejects_dialect(self, tmp_path, monkeypatch):
        # dialect is merged into the kwds inside TextFileReader, i.e. after
        # this check runs; a dialect-specified escapechar would otherwise
        # bypass the escapechar check above (GH#64347).
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(dialect="excel"))

    def test_rejects_parse_dates(self, tmp_path, monkeypatch):
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(parse_dates=["a"]))
        assert not _can_parallelize_csv(path, self._kwds(parse_dates=True))
        assert _can_parallelize_csv(path, self._kwds(parse_dates=None))
        assert _can_parallelize_csv(path, self._kwds(parse_dates=False))

    def test_rejects_storage_options(self, tmp_path):
        # storage_options raises for local paths in the serial path; the
        # parallel path strips it and would mask that error (GH#64347).
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(
            path, self._kwds(storage_options={"anon": True})
        )

    def test_rejects_on_bad_lines_warn(self, tmp_path, monkeypatch):
        # "warn" includes line numbers in its warnings; chunk workers would
        # report chunk-relative (i.e. wrong) ones (GH#64347).
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        method = ParserBase.BadLineHandleMethod
        assert not _can_parallelize_csv(path, self._kwds(on_bad_lines=method.BLHM_WARN))
        assert _can_parallelize_csv(path, self._kwds(on_bad_lines=method.BLHM_ERROR))
        assert _can_parallelize_csv(path, self._kwds(on_bad_lines=method.BLHM_SKIP))

    def test_rejects_blank_line_in_preamble(self, tmp_path, monkeypatch):
        # A blank line before the header shifts where pandas locates the
        # header relative to the physical line count (GH#64347).
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("\na,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds())

    def test_rejects_callable_skiprows(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(skiprows=lambda x: x == 0))

    def test_rejects_list_skiprows(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(skiprows=[0, 2]))

    def test_rejects_multi_level_header(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(header=[0, 1]))

    def test_rejects_index_col(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(index_col=0))
        assert not _can_parallelize_csv(path, self._kwds(index_col="a"))

    def test_accepts_index_col_false(self, tmp_path):
        # index_col=False is allowed (suppresses implicit index)
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        # File is too small → False due to size, not index_col=False
        assert not _can_parallelize_csv(path, self._kwds(index_col=False))

    def test_rejects_usecols(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_text("a,b,c\n1,2,3\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds(usecols=["a", "b"]))

    def test_rejects_small_file(self, tmp_path):
        path = tmp_path / "small.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        assert not _can_parallelize_csv(path, self._kwds())


# ---------------------------------------------------------------------------
# _find_data_start_offset
# ---------------------------------------------------------------------------


class TestFindDataStartOffset:
    def test_no_skiprows_header_0(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_bytes(b"a,b,c\n1,2,3\n4,5,6\n")
        offset = _find_data_start_offset(str(path), header=0, skiprows=0)
        assert offset == len(b"a,b,c\n")

    def test_header_none(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_bytes(b"1,2,3\n4,5,6\n")
        offset = _find_data_start_offset(str(path), header=None, skiprows=0)
        assert offset == 0

    def test_skiprows_2_header_0(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_bytes(b"skip1\nskip2\na,b\n1,2\n")
        offset = _find_data_start_offset(str(path), header=0, skiprows=2)
        # 'skip1\n' + 'skip2\n' + 'a,b\n' = 6 + 6 + 4 = 16
        assert offset == len(b"skip1\nskip2\na,b\n")

    def test_header_2(self, tmp_path):
        # header=2 means lines 0,1,2 are consumed (lines 0-1 skipped, line 2 is header)
        path = tmp_path / "data.csv"
        path.write_bytes(b"r0\nr1\nr2_header\n1,2\n")
        offset = _find_data_start_offset(str(path), header=2, skiprows=0)
        assert offset == len(b"r0\nr1\nr2_header\n")

    def test_eof_within_preamble(self, tmp_path):
        # File ends before all preamble lines are consumed → no crash.
        path = tmp_path / "data.csv"
        path.write_bytes(b"only_one_line\n")
        offset = _find_data_start_offset(str(path), header=0, skiprows=5)
        assert offset == len(b"only_one_line\n")


# ---------------------------------------------------------------------------
# _find_chunk_byte_offsets
# ---------------------------------------------------------------------------


class TestFindChunkByteOffsets:
    def _write_csv(self, path: Path, n_rows: int = 10) -> None:
        lines = "a,b\n" + "".join(f"{i},{i}\n" for i in range(n_rows))
        path.write_text(lines, encoding="utf-8")

    def test_single_chunk(self, tmp_path):
        path = tmp_path / "data.csv"
        self._write_csv(path)
        data_start = len(b"a,b\n")
        offsets = _find_chunk_byte_offsets(str(path), 1, data_start)
        assert offsets[0] == data_start
        assert offsets[-1] == os.path.getsize(str(path))
        assert len(offsets) == 2

    def test_multiple_chunks_cover_full_file(self, tmp_path):
        path = tmp_path / "data.csv"
        self._write_csv(path, n_rows=100)
        data_start = len(b"a,b\n")
        file_size = os.path.getsize(str(path))
        offsets = _find_chunk_byte_offsets(str(path), 4, data_start)

        assert offsets[0] == data_start
        assert offsets[-1] == file_size
        # All offsets within valid range
        assert all(data_start <= off <= file_size for off in offsets)
        # Strictly increasing
        assert all(offsets[i] < offsets[i + 1] for i in range(len(offsets) - 1))

    def test_chunks_align_to_newlines(self, tmp_path):
        path = tmp_path / "data.csv"
        self._write_csv(path, n_rows=1000)
        data_start = len(b"a,b\n")
        raw_bytes = path.read_bytes()
        offsets = _find_chunk_byte_offsets(str(path), 4, data_start)

        for off in offsets[1:-1]:
            # Each interior split point must be at the start of a line.
            assert raw_bytes[off - 1 : off] == b"\n"

    def test_fewer_chunks_than_requested_when_file_tiny(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_bytes(b"a,b\n1,2\n")
        offsets = _find_chunk_byte_offsets(str(path), 8, data_start=4)
        # Can't produce 8 distinct split points in 4 bytes.
        assert len(offsets) >= 2


# ---------------------------------------------------------------------------
# _read_csv_parallel  (correctness: parallel == serial)
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestReadCsvParallel:
    """End-to-end correctness: parallel result must match serial result."""

    def _serial_read(self, path, **kwargs):
        return read_csv(path, engine="c", **kwargs)

    def _parallel_read(self, path, kwds, n_workers=4):
        """Call the internal helper directly so file-size guards don't apply."""
        result = _read_csv_parallel(str(path), kwds, n_workers)
        if result is None:
            pytest.skip("parallel read not applicable to this file")
        return result

    def _base_kwds(self, path, **overrides):
        """Return kwds as _read() would see them for a default read_csv call."""
        from pandas._libs import lib

        from pandas.io.parsers.readers import _refine_defaults_read

        kwds = {
            "header": "infer",
            "names": lib.no_default,
            "index_col": None,
            "usecols": None,
            "dtype": None,
            "converters": None,
            "true_values": None,
            "false_values": None,
            "skipinitialspace": False,
            "skiprows": None,
            "skipfooter": 0,
            "nrows": None,
            "na_values": None,
            "keep_default_na": True,
            "na_filter": True,
            "skip_blank_lines": True,
            "parse_dates": None,
            "date_format": None,
            "dayfirst": False,
            "cache_dates": True,
            "iterator": False,
            "chunksize": None,
            "compression": "infer",
            "thousands": None,
            "decimal": ".",
            "lineterminator": None,
            "quotechar": '"',
            "quoting": 0,
            "doublequote": True,
            "escapechar": None,
            "comment": None,
            "encoding": None,
            "encoding_errors": "strict",
            "dialect": None,
            "on_bad_lines": "error",
            "low_memory": True,
            "memory_map": False,
            "float_precision": None,
            "storage_options": None,
            "dtype_backend": lib.no_default,
        }
        kwds.update(overrides)
        # Apply the same defaulting logic as read_csv → _refine_defaults_read.
        kwds_defaults = _refine_defaults_read(
            dialect=None,
            delimiter=lib.no_default,
            engine=None,  # triggers engine="c"
            sep=lib.no_default,
            on_bad_lines=kwds.pop("on_bad_lines"),
            names=kwds.pop("names"),
            defaults={"delimiter": ","},
            dtype_backend=kwds.pop("dtype_backend"),
        )
        kwds.update(kwds_defaults)
        return kwds

    # ------------------------------------------------------------------
    # Basic dtypes
    # ------------------------------------------------------------------

    def test_ints_and_floats(self, tmp_path):
        path = tmp_path / "data.csv"
        n = 5_000
        rng = np.random.default_rng(0)
        df = DataFrame({"a": rng.integers(0, 100, n), "b": rng.random(n)})
        df.to_csv(path, index=False)
        kwds = self._base_kwds(path)
        result = self._parallel_read(path, kwds)
        expected = self._serial_read(path)
        tm.assert_frame_equal(result.reset_index(drop=True), expected)

    def test_strings(self, tmp_path):
        path = tmp_path / "data.csv"
        n = 5_000
        rng = np.random.default_rng(1)
        df = DataFrame(
            {
                "name": rng.choice(["alice", "bob", "carol"], n),
                "val": rng.integers(0, 10, n),
            }
        )
        df.to_csv(path, index=False)
        kwds = self._base_kwds(path)
        result = self._parallel_read(path, kwds)
        expected = self._serial_read(path)
        tm.assert_frame_equal(result.reset_index(drop=True), expected)

    def test_with_na_values(self, tmp_path):
        path = tmp_path / "data.csv"
        n = 5_000
        rng = np.random.default_rng(2)
        vals = rng.random(n)
        vals[rng.integers(0, n, 200)] = float("nan")
        df = DataFrame({"a": vals, "b": rng.integers(0, 5, n)})
        df.to_csv(path, index=False)
        kwds = self._base_kwds(path)
        result = self._parallel_read(path, kwds)
        expected = self._serial_read(path)
        tm.assert_frame_equal(result.reset_index(drop=True), expected)

    def test_header_none(self, tmp_path):
        path = tmp_path / "data.csv"
        n = 5_000
        rng = np.random.default_rng(3)
        df = DataFrame(rng.integers(0, 100, (n, 4)))
        df.to_csv(path, index=False, header=False)
        kwds = self._base_kwds(path, header=None)
        result = self._parallel_read(path, kwds)
        expected = self._serial_read(path, header=None)
        tm.assert_frame_equal(result.reset_index(drop=True), expected)

    def test_skiprows_int(self, tmp_path):
        path = tmp_path / "data.csv"
        n = 5_000
        lines = ["skip_me\n"] * 3 + ["a,b\n"] + [f"{i},{i}\n" for i in range(n)]
        path.write_text("".join(lines), encoding="utf-8")
        kwds = self._base_kwds(path, skiprows=3)
        result = self._parallel_read(path, kwds)
        expected = self._serial_read(path, skiprows=3)
        tm.assert_frame_equal(result.reset_index(drop=True), expected)

    def test_implicit_index_returns_none(self, tmp_path):
        # Data rows have one more field than the header (implicit index).
        # concat(ignore_index=True) would silently drop the index, so the
        # helper must bail out up front (GH#64347).
        path = tmp_path / "implicit.csv"
        lines = "a,b\n" + "".join(f"idx{i},{i},{i * 2}\n" for i in range(5_000))
        path.write_text(lines, encoding="utf-8")
        kwds = self._base_kwds(path)
        assert _read_csv_parallel(str(path), kwds, 4) is None

    def test_quoted_newline_in_header_returns_none(self, tmp_path):
        # A quoted embedded newline in the header makes the physical line
        # count disagree with the logical row count, so data_start would land
        # mid-header and chunk 0 would inject a garbage row (GH#64347).
        path = tmp_path / "hdr.csv"
        lines = '"col\n1",col2\n' + "".join(f"x{i},y{i}\n" for i in range(5_000))
        path.write_text(lines, encoding="utf-8")
        kwds = self._base_kwds(path)
        assert _read_csv_parallel(str(path), kwds, 4) is None

    def test_dtype_specified(self, tmp_path):
        path = tmp_path / "data.csv"
        n = 5_000
        rng = np.random.default_rng(5)
        df = DataFrame({"a": rng.integers(0, 100, n), "b": rng.integers(0, 100, n)})
        df.to_csv(path, index=False)
        kwds = self._base_kwds(path, dtype={"a": "float64", "b": "float64"})
        result = self._parallel_read(path, kwds)
        expected = self._serial_read(path, dtype={"a": "float64", "b": "float64"})
        tm.assert_frame_equal(result.reset_index(drop=True), expected)


# ---------------------------------------------------------------------------
# Embedded newlines in quoted fields
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "csv_bytes,expected_rows",
    [
        # Simple embedded newline
        pytest.param(
            b'a,b,c\n1,"hello\nworld",3\n4,normal,6\n',
            2,
            id="simple_embedded_newline",
        ),
        # Embedded newline where the text after looks like a complete data row
        # (the most dangerous case for silent corruption)
        pytest.param(
            b'a,b,c\n1,"hello\n2,foo,bar",3\n100,baz,200\n',
            2,
            id="adversarial_after_split_looks_like_full_row",
        ),
        # Multiple embedded newlines in the same file
        pytest.param(
            b'a,b\n1,"line1\nline2\nline3"\n2,normal\n3,"also\nembedded"\n',
            3,
            id="multiple_embedded_newlines",
        ),
    ],
)
def test_embedded_newline_falls_back_to_serial(
    tmp_path, monkeypatch, csv_bytes, expected_rows
):
    """
    CSV files with newlines inside quoted fields must produce the same result
    as a direct serial read.  The parallel path either parses correctly (when
    the split misses the embedded newline) or falls back to serial (when it
    hits it).  Silent data corruption must never occur.
    """
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)
    path = tmp_path / "embedded.csv"
    path.write_bytes(csv_bytes)

    serial = read_csv(io.BytesIO(csv_bytes))
    parallel = read_csv(path)  # may take parallel path or fall back

    assert len(parallel) == expected_rows
    tm.assert_frame_equal(
        parallel.reset_index(drop=True), serial.reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Integration: automatic activation through read_csv()
# ---------------------------------------------------------------------------


def test_read_csv_auto_parallel(tmp_path, monkeypatch):
    """
    read_csv() transparently uses the parallel path for large local files.
    Result must match the serial read obtained via engine='python'.
    """
    import pandas.io.parsers.readers as _readers

    path = tmp_path / "big.csv"
    _make_large_csv(path)
    # Lower the threshold so any file triggers the parallel path.
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)

    # Force the parallel path regardless of the platform default (it is off by
    # default on Windows) so this exercises parallel == serial on every platform.
    with option_context("mode.max_threads", 4):
        result = read_csv(path)  # auto-selects parallel path for C engine
    expected = read_csv(path, engine="python")
    tm.assert_frame_equal(result, expected)


def test_read_csv_parallel_vs_serial_large_file(tmp_path, monkeypatch):
    """
    For a file that exceeds the threshold, the parallel result equals the
    result from a direct serial C-engine read (with parallelism forced off).
    """
    import pandas.io.parsers.readers as _readers

    path = tmp_path / "big.csv"
    _make_large_csv(path)
    # Lower the threshold so any file triggers the parallel path.
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)

    serial = read_csv(path, engine="python")
    # Force the parallel path so this runs on every platform (off by default on
    # Windows).
    with option_context("mode.max_threads", 4):
        parallel = read_csv(path, engine="c")
    tm.assert_frame_equal(parallel, serial)


def test_parallel_default_off_on_windows(tmp_path, monkeypatch):
    """The parallel path is off by default on Windows but on elsewhere.

    Windows shows no speedup (and a slowdown at two threads) even with the file
    warm in the OS cache, so the default is serial there; users opt in via
    ``mode.max_threads`` (which is honoured on every platform).
    """
    import pandas.io.parsers.readers as _readers

    path = tmp_path / "big.csv"
    _make_large_csv(path)
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
    # Pin the non-Windows default so the test does not depend on the host's
    # actual core count.
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 4)

    # Stub the parallel reader so this exercises only the platform-gating
    # decision.  Calling the real one would start threads, which fails on
    # no-thread platforms such as Pyodide/WASM.
    calls = []

    def stub(*args, **kwargs):
        calls.append(args)
        return DataFrame()

    monkeypatch.setattr(_readers, "_read_csv_parallel", stub)

    monkeypatch.setattr(_readers.sys, "platform", "win32")
    read_csv(path)
    assert calls == []  # serial by default on Windows

    monkeypatch.setattr(_readers.sys, "platform", "linux")
    read_csv(path)
    assert len(calls) == 1  # parallel by default elsewhere


def test_parallel_default_thread_cap(tmp_path, monkeypatch):
    """The default worker count is capped at 4, regardless of core count."""
    import pandas.io.parsers.readers as _readers

    path = tmp_path / "big.csv"
    _make_large_csv(path)
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
    monkeypatch.setattr(_readers.sys, "platform", "linux")
    # More cores than the cap: the default should clamp down to 4.
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 16)

    workers = []

    def stub(_path, _kwds, n_workers):
        workers.append(n_workers)
        return DataFrame()

    monkeypatch.setattr(_readers, "_read_csv_parallel", stub)

    read_csv(path)
    assert workers == [4]

    # An explicit mode.max_threads still overrides the cap.
    workers.clear()
    with option_context("mode.max_threads", 8):
        read_csv(path)
    assert workers == [8]


# ---------------------------------------------------------------------------
# Regression tests: inputs the parallel path must hand back to serial
# (all results must be identical to a serial read of the same bytes)
# ---------------------------------------------------------------------------


def _read_forced_parallel(path, monkeypatch, **kwargs):
    """read_csv with the parallel path force-enabled on every platform."""
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)
    with option_context("mode.max_threads", 4):
        return read_csv(path, **kwargs)


def test_parallel_latin1_matches_serial(tmp_path, monkeypatch):
    # Non-UTF-8 encodings must take the serial path; the C tokenizer decodes
    # raw chunk bytes as UTF-8 (GH#64347).
    raw = b"col1,col2\n" + b"".join(
        f"se\xf1or{i},{i}\n".encode("latin-1") for i in range(500)
    )
    path = tmp_path / "latin.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, encoding="latin-1")
    expected = read_csv(io.BytesIO(raw), encoding="latin-1")
    tm.assert_frame_equal(result, expected)
    assert result.loc[0, "col1"] == "se\xf1or0"


def test_parallel_comment_before_header_matches_serial(tmp_path, monkeypatch):
    # A full-line comment above the header shifts the preamble byte offset;
    # the header row must not be ingested as data (GH#64347).
    raw = b"# a comment\ncol1,col2\n" + b"".join(
        f"{i},{i * 2}\n".encode() for i in range(500)
    )
    path = tmp_path / "comment.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, comment="#")
    expected = read_csv(io.BytesIO(raw), comment="#")
    tm.assert_frame_equal(result, expected)
    assert result["col1"].dtype == np.int64


def test_parallel_blank_line_before_header_matches_serial(tmp_path, monkeypatch):
    # Same as the comment case, but triggered by a plain blank line (GH#64347).
    raw = b"\ncol1,col2\n" + b"".join(f"{i},{i * 2}\n".encode() for i in range(500))
    path = tmp_path / "blank.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch)
    expected = read_csv(io.BytesIO(raw))
    tm.assert_frame_equal(result, expected)
    assert result["col1"].dtype == np.int64


def test_parallel_single_long_line(tmp_path, monkeypatch):
    # A data section with no interior newlines cannot be split; this must
    # fall back to serial instead of raising (GH#64347).
    raw = b"col1,col2\n" + b"a" * 5_000 + b",b"
    path = tmp_path / "longline.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch)
    expected = read_csv(io.BytesIO(raw))
    tm.assert_frame_equal(result, expected)


def test_parallel_mixed_dtype_column_matches_serial(tmp_path, monkeypatch):
    # A numeric column whose only non-numeric row sits in one chunk parses as
    # int64 in the clean chunks but object in the dirty one; serial gives
    # strings for the whole column.  The parallel path must detect the
    # disagreement and re-read serially (GH#64347).
    clean = b"".join(f"{i},{i}\n".encode() for i in range(500))
    raw = b"col1,col2\n" + clean + b"oops,999\n" + clean
    path = tmp_path / "mixed.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch)
    expected = read_csv(io.BytesIO(raw))
    tm.assert_frame_equal(result, expected)
    # serial semantics: the whole column stays as strings
    assert result.loc[0, "col1"] == "0"


def test_parallel_quoted_newline_in_header_matches_serial(tmp_path, monkeypatch):
    # A quoted embedded newline in the header shifts data_start mid-header;
    # without a guard, chunk 0 starts inside the header and silently injects
    # a garbage row (GH#64347).  All-string data so the per-column dtype
    # consistency check cannot catch it.
    raw = b'"col\n1",col2\n' + b"".join(f"x{i},y{i}\n".encode() for i in range(500))
    path = tmp_path / "hdr.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch)
    expected = read_csv(io.BytesIO(raw))
    tm.assert_frame_equal(result, expected)
    assert len(result) == 500


def test_parallel_storage_options_raises(tmp_path, monkeypatch):
    # storage_options with a local path raises ValueError in the serial path;
    # the parallel path strips it and must not mask that error (GH#64347).
    raw = b"a,b\n" + b"".join(f"{i},{i}\n".encode() for i in range(500))
    path = tmp_path / "so.csv"
    path.write_bytes(raw)

    msg = "storage_options passed with file object or non-fsspec file path"
    with pytest.raises(ValueError, match=msg):
        _read_forced_parallel(path, monkeypatch, storage_options={"anon": True})


def test_parallel_on_bad_lines_warn_line_numbers(tmp_path, monkeypatch):
    # on_bad_lines="warn" must report absolute line numbers; chunk workers
    # would report chunk-relative ones, so the serial path is used (GH#64347).
    raw = (
        b"a,b\n"
        + b"".join(f"{i},{i}\n".encode() for i in range(300))
        + b"1,2,3,4\n"
        + b"".join(f"{i},{i}\n".encode() for i in range(300, 600))
    )
    path = tmp_path / "bad.csv"
    path.write_bytes(raw)

    with tm.assert_produces_warning(ParserWarning, match="line 302"):
        result = _read_forced_parallel(path, monkeypatch, on_bad_lines="warn")
    with tm.assert_produces_warning(ParserWarning, match="line 302"):
        expected = read_csv(io.BytesIO(raw), on_bad_lines="warn")
    tm.assert_frame_equal(result, expected)


def test_parallel_implicit_index_matches_serial(tmp_path, monkeypatch):
    # Data rows with one more field than the header get an implicit index;
    # the parallel path must hand back to serial rather than drop it via
    # concat(ignore_index=True) (GH#64347).
    raw = b"col1,col2\n" + b"".join(
        f"idx{i},{i},{i * 2}\n".encode() for i in range(500)
    )
    path = tmp_path / "implicit.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch)
    expected = read_csv(io.BytesIO(raw))
    tm.assert_frame_equal(result, expected)
    assert list(result.index[:2]) == ["idx0", "idx1"]


def test_parallel_dialect_escapechar_matches_serial(tmp_path, monkeypatch):
    # dialect is merged into the kwds only inside TextFileReader, after the
    # eligibility check runs; a dialect-specified escapechar must still force
    # the serial path - an escaped newline split across chunks would
    # otherwise corrupt data (GH#64347).
    class EscDialect(csv.Dialect):
        delimiter = ","
        quotechar = '"'
        doublequote = True
        skipinitialspace = False
        lineterminator = "\r\n"
        quoting = csv.QUOTE_NONE
        escapechar = "\\"

    raw = b"col1,col2\n" + b"".join(
        f"text\\\nmore{i},{i}\n".encode() for i in range(500)
    )
    path = tmp_path / "dialect.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, dialect=EscDialect)
    expected = read_csv(io.BytesIO(raw), dialect=EscDialect)
    tm.assert_frame_equal(result, expected)
    assert result.loc[0, "col1"] == "text\nmore0"


def test_parallel_multichar_sep_matches_serial(tmp_path, monkeypatch):
    # Multi-char seps force the python engine inside TextFileReader; the
    # parallel path must step aside rather than hit the C-engine assert
    # (GH#64347).
    raw = b"col1;;col2\n" + b"".join(f"{i};;{i * 2}\n".encode() for i in range(500))
    path = tmp_path / "multisep.csv"
    path.write_bytes(raw)

    warn_msg = "Falling back to the 'python' engine"
    with tm.assert_produces_warning(ParserWarning, match=warn_msg):
        result = _read_forced_parallel(path, monkeypatch, sep=";;")
    with tm.assert_produces_warning(ParserWarning, match=warn_msg):
        expected = read_csv(io.BytesIO(raw), sep=";;")
    tm.assert_frame_equal(result, expected)
