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
import ctypes
import io
import os
from typing import TYPE_CHECKING

import numpy as np
import pytest

from pandas.compat._cpu import (
    _count_distinct_cores,
    _count_processor_core_records,
    _parse_cgroup_v2_quota,
    _parse_cpu_list,
    available_cpu_count,
    physical_core_count,
)
from pandas.errors import (
    ParserError,
    ParserWarning,
)

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
    _default_n_workers,
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

    def test_rejects_tar(self, tmp_path):
        # .tar infers as compression="tar" in get_handle; the parallel path
        # must not feed raw archive bytes to the tokenizer.
        path = tmp_path / "data.tar"
        path.write_bytes(b"fake")
        assert not _can_parallelize_csv(path, self._kwds())

    def test_rejects_dict_compression(self, tmp_path):
        path = tmp_path / "data.csv"
        path.write_bytes(b"fake")
        kwds = self._kwds(compression={"method": "gzip"})
        assert not _can_parallelize_csv(path, kwds)

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

    @pytest.mark.parametrize("encoding", ["latin-1", "cp1252", "shift_jis", "ascii"])
    def test_rejects_non_utf8_encoding(self, tmp_path, monkeypatch, encoding):
        # Chunk workers feed raw file bytes to the C tokenizer, which decodes
        # words as UTF-8; only UTF-8-compatible encodings are byte-safe.  "ascii"
        # is excluded too: the workers would decode a non-ASCII byte as UTF-8 and
        # succeed, masking the UnicodeDecodeError serial raises (GH#64347).
        import pandas.io.parsers.readers as _readers

        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(encoding=encoding))
        assert _can_parallelize_csv(path, self._kwds(encoding="utf-8"))
        assert _can_parallelize_csv(path, self._kwds(encoding="UTF8"))
        assert _can_parallelize_csv(path, self._kwds(encoding="utf-8-sig"))

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

    def test_dict_dtype_duplicate_names_falls_back(self, tmp_path):
        # A dict dtype keyed on a duplicated header name must return None so the
        # caller re-reads serially; the same dtype on unique names must still
        # parallelize (GH#64347).
        path = tmp_path / "data.csv"
        body = "".join(f"{i},{i + 1}\n" for i in range(5_000))

        path.write_text("a,a\n" + body, encoding="utf-8")
        kwds = self._base_kwds(path, dtype={"a": str})
        assert _read_csv_parallel(str(path), kwds, 4) is None

        path.write_text("a,b\n" + body, encoding="utf-8")
        kwds = self._base_kwds(path, dtype={"a": str})
        assert _read_csv_parallel(str(path), kwds, 4) is not None


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
    # actual core count / topology / cgroup limits.
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 4)
    monkeypatch.setattr(_readers, "physical_core_count", lambda: 4)
    monkeypatch.setattr(_readers, "available_cpu_count", lambda: None)

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


def test_parallel_default_uses_physical_cores(tmp_path, monkeypatch):
    """The default worker count follows the detected physical core count."""
    import pandas.io.parsers.readers as _readers

    path = tmp_path / "big.csv"
    _make_large_csv(path)
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
    monkeypatch.setattr(_readers.sys, "platform", "linux")
    # 6 physical cores out of 16 logical CPUs (SMT): the default follows the
    # physical core count, not the logical one.
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 16)
    monkeypatch.setattr(_readers, "physical_core_count", lambda: 6)
    monkeypatch.setattr(_readers, "available_cpu_count", lambda: None)

    workers = []

    def stub(_path, _kwds, n_workers):
        workers.append(n_workers)
        return DataFrame()

    monkeypatch.setattr(_readers, "_read_csv_parallel", stub)

    read_csv(path)
    assert workers == [6]

    # An explicit mode.max_threads still overrides the physical-core default.
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


def test_parallel_ascii_non_ascii_byte_matches_serial(tmp_path, monkeypatch):
    # encoding="ascii" with a non-ASCII (but valid UTF-8) byte deep in the
    # file: the serial path raises UnicodeDecodeError, but the workers decode
    # chunk bytes as UTF-8 and would silently succeed, so ascii must take the
    # serial path (GH#64347).
    raw = (
        b"a,b\n"
        + b"".join(f"{i},x{i}\n".encode() for i in range(500))
        + b"500,\xc3\xa9\n"  # U+00E9 (é): valid UTF-8, not ASCII
        + b"".join(f"{i},x{i}\n".encode() for i in range(501, 1000))
    )
    path = tmp_path / "ascii.csv"
    path.write_bytes(raw)

    msg = "'ascii' codec can't decode"
    with pytest.raises(UnicodeDecodeError, match=msg):
        _read_forced_parallel(path, monkeypatch, encoding="ascii")
    with pytest.raises(UnicodeDecodeError, match=msg):
        read_csv(io.BytesIO(raw), encoding="ascii")


def test_parallel_dup_names_dict_dtype_matches_serial(tmp_path, monkeypatch):
    # A dict dtype keyed on a duplicated header name is applied by serial to
    # every de-duplicated column ("a" and "a.1"), but the workers receive the
    # already de-duplicated names and would match only "a".  The parallel path
    # must recover the raw header names and hand back to serial (GH#64347).
    raw = b"a,a\n" + b"".join(f"{i},{i + 1}\n".encode() for i in range(500))
    path = tmp_path / "dup.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, dtype={"a": str})
    expected = read_csv(io.BytesIO(raw), dtype={"a": str})
    tm.assert_frame_equal(result, expected)
    # the dtype reaches BOTH duplicated columns, not just the first
    assert result.dtypes["a"] == result.dtypes["a.1"]


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


def _write_with_line_at_chunk_start(path, replacement: bytes) -> None:
    """Write a fixed-width CSV, then overwrite the line at the second chunk
    boundary with *replacement* (same byte length, so offsets stay valid)."""
    rows = [f"{i:06d},{i * 2:06d}" for i in range(4000)]
    # write bytes so line endings stay "\n" on every platform - the byte-offset
    # math below assumes single-byte terminators (Windows text mode adds "\r")
    path.write_bytes(("a,b\n" + "\n".join(rows) + "\n").encode("utf-8"))
    data_start = _find_data_start_offset(str(path), 0, 0)
    offsets = _find_chunk_byte_offsets(str(path), 4, data_start)
    boundary = offsets[1]
    raw = path.read_bytes()
    line_end = raw.index(b"\n", boundary)
    assert line_end - boundary == len(replacement)
    with open(path, "r+b") as fd:
        fd.seek(boundary)
        fd.write(replacement)


def test_parallel_ragged_line_at_chunk_start_raises(tmp_path, monkeypatch):
    # A line with extra fields sitting exactly at a chunk boundary: the chunk
    # worker's first line used to be exempt from the field-count check, so
    # parallel silently truncated the row while serial raised (GH#64347).
    path = tmp_path / "ragged.csv"
    _write_with_line_at_chunk_start(path, b"11,22,33,4444")

    with pytest.raises(ParserError, match="Expected 2 fields"):
        _read_forced_parallel(path, monkeypatch)


def test_parallel_ragged_line_at_chunk_start_skip_matches_serial(tmp_path, monkeypatch):
    # Same layout with on_bad_lines="skip": the chunk worker must skip the
    # bad line just like serial, not keep a truncated version of it.
    path = tmp_path / "ragged.csv"
    _write_with_line_at_chunk_start(path, b"11,22,33,4444")

    with option_context("mode.max_threads", 1):
        expected = read_csv(path, on_bad_lines="skip")
    result = _read_forced_parallel(path, monkeypatch, on_bad_lines="skip")
    tm.assert_frame_equal(result, expected)
    assert len(result) == 3999


def test_parallel_bom_bytes_mid_file_match_serial(tmp_path, monkeypatch):
    # Data lines beginning with the UTF-8 BOM byte sequence: each chunk
    # worker's fresh parser used to strip it from the chunk's first line,
    # while serial strips only at byte 0 of the file (GH#64347).
    lines = ["a,b"] + ["\ufeff" + f"x{i},{i}" for i in range(4000)]
    path = tmp_path / "bom.csv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with option_context("mode.max_threads", 1):
        expected = read_csv(path)
    result = _read_forced_parallel(path, monkeypatch)
    tm.assert_frame_equal(result, expected)
    assert (result["a"].str[0] == "\ufeff").all()


def test_parallel_bom_at_file_start_header_none(tmp_path, monkeypatch):
    # With header=None and no skiprows the first chunk starts at byte 0 and
    # must strip a leading BOM exactly like the serial parser does.
    body = "".join(f"{i},{i * 2}\n" for i in range(4000))
    path = tmp_path / "bom0.csv"
    path.write_bytes(b"\xef\xbb\xbf" + body.encode())

    with option_context("mode.max_threads", 1):
        expected = read_csv(path, header=None, names=["a", "b"])
    result = _read_forced_parallel(path, monkeypatch, header=None, names=["a", "b"])
    tm.assert_frame_equal(result, expected)
    # the stripped BOM means the first cell parses as an integer
    assert result["a"].iloc[0] == 0


# ---------------------------------------------------------------------------
# Performance-core detection and the default worker count
# ---------------------------------------------------------------------------


def test_physical_core_count_within_bounds():
    # Whatever the host topology, the detector returns a positive int no larger
    # than the logical CPU count, and never raises.
    physical_core_count.cache_clear()
    try:
        n_perf = physical_core_count()
    finally:
        physical_core_count.cache_clear()
    assert isinstance(n_perf, int)
    assert 1 <= n_perf <= (os.cpu_count() or 1)


@pytest.mark.parametrize(
    "platform_name, probe_name",
    [
        ("darwin", "_physical_cores_darwin"),
        ("linux", "_physical_cores_linux"),
        ("win32", "_physical_cores_windows"),
    ],
)
@pytest.mark.parametrize("failure", ["returns_none", "raises"])
def test_physical_core_count_falls_back(
    monkeypatch, platform_name, probe_name, failure
):
    # A probe that returns None or raises must fall back to os.cpu_count().
    from pandas.compat import _cpu

    def probe():
        if failure == "raises":
            raise OSError("probe failed")

    monkeypatch.setattr(_cpu.sys, "platform", platform_name)
    monkeypatch.setattr(_cpu, probe_name, probe)
    _cpu.physical_core_count.cache_clear()
    try:
        assert _cpu.physical_core_count() == (os.cpu_count() or 1)
    finally:
        _cpu.physical_core_count.cache_clear()


@pytest.mark.parametrize(
    "probe_result, expected",
    [
        (1, 1),  # a valid subset is used as-is
        (0, "total"),  # nonsensical count -> fall back
        (10**6, "total"),  # more cores than exist -> fall back
    ],
)
def test_physical_core_count_validates_probe(monkeypatch, probe_result, expected):
    from pandas.compat import _cpu

    total = os.cpu_count() or 1
    monkeypatch.setattr(_cpu.sys, "platform", "darwin")
    monkeypatch.setattr(_cpu, "_physical_cores_darwin", lambda: probe_result)
    _cpu.physical_core_count.cache_clear()
    try:
        result = _cpu.physical_core_count()
    finally:
        _cpu.physical_core_count.cache_clear()
    assert result == (total if expected == "total" else expected)


@pytest.mark.parametrize(
    "n_perf, available, expected",
    [
        (6, None, 6),  # unconstrained, below the cap -> detected count
        (1, None, 1),
        (24, None, 16),  # e.g. an M-series Ultra -> capped
        (32, None, 16),  # cap binds
        (8, 4, 4),  # cgroup/affinity tighter than the perf-core count
        (8, 12, 8),  # allocation looser than the perf-core count
        (32, 8, 8),  # allocation tighter than both cap and perf cores
    ],
)
def test_default_n_workers_combines_detection_allocation_cap(
    monkeypatch, n_perf, available, expected
):
    # Default = min(physical performance cores, available CPUs, _MAX_DEFAULT_WORKERS).
    import pandas.io.parsers.readers as _readers

    assert _readers._MAX_DEFAULT_WORKERS == 16
    monkeypatch.setattr(_readers.sys, "platform", "linux")
    monkeypatch.setattr(_readers, "physical_core_count", lambda: n_perf)
    monkeypatch.setattr(_readers, "available_cpu_count", lambda: available)
    with option_context("mode.max_threads", None):
        assert _default_n_workers() == expected


def test_default_n_workers_windows_is_serial(monkeypatch):
    import pandas.io.parsers.readers as _readers

    monkeypatch.setattr(_readers.sys, "platform", "win32")
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 12)
    monkeypatch.setattr(_readers, "physical_core_count", lambda: 6)
    with option_context("mode.max_threads", None):
        assert _default_n_workers() == 1


def test_default_n_workers_wasm_is_serial(monkeypatch):
    import pandas.io.parsers.readers as _readers

    monkeypatch.setattr(_readers.sys, "platform", "emscripten")
    # WASM stays serial even if a worker count is requested explicitly.
    with option_context("mode.max_threads", 8):
        assert _default_n_workers() == 1


@pytest.mark.parametrize("platform_name", ["linux", "win32"])
def test_default_n_workers_max_threads_wins(monkeypatch, platform_name):
    import pandas.io.parsers.readers as _readers

    monkeypatch.setattr(_readers.sys, "platform", platform_name)
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 12)
    monkeypatch.setattr(_readers, "physical_core_count", lambda: 6)
    with option_context("mode.max_threads", 3):
        assert _default_n_workers() == 3


@pytest.mark.parametrize(
    "spec, expected",
    [
        ("0-5", [0, 1, 2, 3, 4, 5]),
        ("0-3,8", [0, 1, 2, 3, 8]),
        ("0,2,4", [0, 2, 4]),
        ("3", [3]),
        ("", []),
    ],
)
def test_parse_cpu_list(spec, expected):
    assert _parse_cpu_list(spec) == expected


@pytest.mark.parametrize(
    "topology, expected",
    [
        # 4 SMT threads -> 2 physical cores (siblings collapse by (pkg, core)).
        ([(0, 0), (0, 0), (0, 1), (0, 1)], 2),
        # two packages, one core each
        ([(0, 0), (1, 0)], 2),
        # unreadable entries are ignored
        ([(0, 0), None, (0, 0)], 1),
        ([None, None], None),
        ([], None),
    ],
)
def test_count_distinct_cores(topology, expected):
    assert _count_distinct_cores(topology) == expected


@pytest.mark.parametrize(
    "text, expected",
    [
        ("max 100000", None),  # unlimited
        ("max", None),
        ("400000 100000", 4.0),
        ("150000 100000", 1.5),
        ("100000", 1.0),  # period defaults to 100000us
        ("0 100000", None),  # zero quota -> ignored
        ("garbage", None),
        ("", None),
    ],
)
def test_parse_cgroup_v2_quota(text, expected):
    assert _parse_cgroup_v2_quota(text) == expected


@pytest.mark.parametrize(
    "affinity, quota, expected",
    [
        ({0, 1, 2, 3, 4, 5, 6, 7}, None, 8),  # affinity only
        ({0, 1, 2, 3, 4, 5, 6, 7}, 4.0, 4),  # cgroup quota tighter
        ({0, 1, 2, 3}, 8.0, 4),  # affinity tighter
        ({0, 1}, 1.5, 1),  # fractional quota floors, min 1
    ],
)
def test_available_cpu_count(monkeypatch, affinity, quota, expected):
    from pandas.compat import _cpu

    monkeypatch.setattr(
        os, "sched_getaffinity", lambda pid: set(affinity), raising=False
    )
    monkeypatch.setattr(_cpu, "_cgroup_cpu_quota", lambda: quota)
    available_cpu_count.cache_clear()
    try:
        assert available_cpu_count() == expected
    finally:
        available_cpu_count.cache_clear()


def test_available_cpu_count_unconstrained(monkeypatch):
    # No affinity limit and no cgroup quota -> None (do not clamp).
    from pandas.compat import _cpu

    def raise_oserror(pid):
        raise OSError

    monkeypatch.setattr(os, "sched_getaffinity", raise_oserror, raising=False)
    monkeypatch.setattr(_cpu, "_cgroup_cpu_quota", lambda: None)
    available_cpu_count.cache_clear()
    try:
        assert available_cpu_count() is None
    finally:
        available_cpu_count.cache_clear()


def _fake_cpu_topology_reader(cores_per_cpu):
    """Build a fake ``_read_sysfs_int`` returning core/package ids from a map."""
    import re

    def read_int(path):
        match = re.search(r"/cpu(\d+)/topology/(\w+)", path)
        if match is None:
            return None
        cpu, field = int(match.group(1)), match.group(2)
        if field == "physical_package_id":
            return 0
        if field == "core_id":
            return cores_per_cpu.get(cpu)
        return None

    return read_int


def test_physical_cores_linux_hybrid_collapses_smt(monkeypatch):
    # Hybrid part: logical CPUs 0-7 are 4 physical P-cores with 2 SMT threads
    # each, and CPUs 8-11 are 4 efficiency cores.  SMT siblings collapse;
    # efficiency cores count.
    from pandas.compat import _cpu

    strs = {"/sys/devices/system/cpu/present": "0-11"}
    cores = {cpu: cpu // 2 for cpu in range(8)}
    cores.update({cpu: cpu - 4 for cpu in range(8, 12)})
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: strs.get(path))
    monkeypatch.setattr(_cpu, "_read_sysfs_int", _fake_cpu_topology_reader(cores))
    assert _cpu._physical_cores_linux() == 8


def test_physical_cores_linux_homogeneous_uses_present(monkeypatch):
    # All present CPUs, SMT collapsed (0-3 = 2 physical cores).
    from pandas.compat import _cpu

    strs = {"/sys/devices/system/cpu/present": "0-3"}
    cores = {cpu: cpu // 2 for cpu in range(4)}
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: strs.get(path))
    monkeypatch.setattr(_cpu, "_read_sysfs_int", _fake_cpu_topology_reader(cores))
    assert _cpu._physical_cores_linux() == 2


def test_cgroup_cpu_quota_v2(monkeypatch):
    from pandas.compat import _cpu

    monkeypatch.setattr(
        _cpu,
        "_read_sysfs_str",
        lambda path: "400000 100000" if path == "/sys/fs/cgroup/cpu.max" else None,
    )
    assert _cpu._cgroup_cpu_quota() == 4.0


def test_cgroup_cpu_quota_v1_fallback(monkeypatch):
    from pandas.compat import _cpu

    ints = {
        "/sys/fs/cgroup/cpu/cpu.cfs_quota_us": 300000,
        "/sys/fs/cgroup/cpu/cpu.cfs_period_us": 100000,
    }
    monkeypatch.setattr(_cpu, "_read_sysfs_str", lambda path: None)  # no cgroup v2
    monkeypatch.setattr(_cpu, "_read_sysfs_int", lambda path: ints.get(path))
    assert _cpu._cgroup_cpu_quota() == 3.0


def test_cgroup_cpu_quota_unlimited(monkeypatch):
    from pandas.compat import _cpu

    monkeypatch.setattr(
        _cpu,
        "_read_sysfs_str",
        lambda path: "max 100000" if path == "/sys/fs/cgroup/cpu.max" else None,
    )
    monkeypatch.setattr(_cpu, "_read_sysfs_int", lambda path: None)
    assert _cpu._cgroup_cpu_quota() is None


def _make_win_processor_buffer(records):
    """Build a synthetic ``GetLogicalProcessorInformationEx`` buffer.

    ``records`` is a list of ``(size, efficiency_class)`` pairs, one per
    physical-core relationship record.  ``Size`` sits at byte offset ``+4`` and
    ``EfficiencyClass`` at ``+9`` within each record.
    """
    total = sum(size for size, _ in records)
    buf = (ctypes.c_byte * total)()
    addr = ctypes.addressof(buf)
    offset = 0
    for size, eff in records:
        ctypes.c_uint32.from_address(addr + offset + 4).value = size
        ctypes.c_uint8.from_address(addr + offset + 9).value = eff
        offset += size
    return buf, total


def test_count_processor_core_records_hybrid():
    # 4 performance cores (class 1) + 4 efficiency cores (class 0): every
    # physical core counts, efficiency class does not matter.
    buf, length = _make_win_processor_buffer([(32, 1)] * 4 + [(32, 0)] * 4)
    assert _count_processor_core_records(buf, length) == 8


def test_count_processor_core_records_homogeneous():
    buf, length = _make_win_processor_buffer([(32, 0)] * 8)
    assert _count_processor_core_records(buf, length) == 8


def test_count_processor_core_records_variable_record_size():
    # Records advance by their own Size field, not a fixed stride.
    buf, length = _make_win_processor_buffer([(40, 1), (24, 1), (32, 0)])
    assert _count_processor_core_records(buf, length) == 3


def test_count_processor_core_records_empty():
    buf = (ctypes.c_byte * 0)()
    assert _count_processor_core_records(buf, 0) is None


def test_parallel_string_dtype_python_storage(tmp_path, monkeypatch):
    # Without the pyarrow string fast path, the gathered object column must
    # get the same string-dtype inference a serial read applies via the
    # DataFrame constructor (GH#66275).
    import pandas.io.parsers.readers as _readers

    path = tmp_path / "data.csv"
    rows = "\n".join(f"{i},s{i % 7}" for i in range(5000))
    path.write_text("a,b\n" + rows + "\n", encoding="utf-8")
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
    with option_context("mode.string_storage", "python"):
        with option_context("mode.max_threads", 1):
            serial = read_csv(path)
        with option_context("mode.max_threads", 4):
            parallel = read_csv(path)
    tm.assert_frame_equal(parallel, serial)
