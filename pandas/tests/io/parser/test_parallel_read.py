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
import mmap
import os
from typing import TYPE_CHECKING
import warnings

import numpy as np
import pytest

from pandas.compat import WASM
from pandas.errors import (
    EmptyDataError,
    ParserError,
    ParserWarning,
)

if TYPE_CHECKING:
    from pathlib import Path

import pandas as pd
import pandas._testing as tm

from pandas.io.parsers import readers as _readers
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
    df = pd.DataFrame(
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
        path = tmp_path / "big.csv"
        _make_large_csv(path)
        # Lower the threshold so the test file qualifies regardless of its size.
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert _can_parallelize_csv(path, self._kwds())

    def test_accepts_default_engine(self, tmp_path, monkeypatch):
        """engine=None (the default) should be treated the same as engine='c'."""
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
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(delimiter=";;"))
        assert not _can_parallelize_csv(path, self._kwds(delimiter=None))
        assert _can_parallelize_csv(path, self._kwds(delimiter=r"\s+"))

    def test_rejects_multibyte_sep_and_quotechar(self, tmp_path, monkeypatch):
        # A separator or quotechar wider than one byte forces the python
        # engine too, warning as it does so (GH#66259).
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(delimiter="§"))
        assert not _can_parallelize_csv(path, self._kwds(quotechar="»"))
        assert _can_parallelize_csv(path, self._kwds(quotechar="'"))

    def test_rejects_comment(self, tmp_path, monkeypatch):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(comment="#"))

    def test_rejects_escapechar(self, tmp_path, monkeypatch):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(escapechar="\\"))

    def test_rejects_dialect(self, tmp_path, monkeypatch):
        # dialect is merged into the kwds inside TextFileReader, i.e. after
        # this check runs; a dialect-specified escapechar would otherwise
        # bypass the escapechar check above (GH#64347).
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(dialect="excel"))

    def test_rejects_parse_dates(self, tmp_path, monkeypatch):
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(parse_dates=["a"]))
        assert not _can_parallelize_csv(path, self._kwds(parse_dates=True))
        assert _can_parallelize_csv(path, self._kwds(parse_dates=None))
        assert _can_parallelize_csv(path, self._kwds(parse_dates=False))

    @pytest.mark.parametrize("low_memory", [False, 0, np.False_])
    def test_rejects_low_memory_false(self, tmp_path, monkeypatch, low_memory):
        # A falsy low_memory guarantees whole-file type inference, which the
        # per-chunk parallel path cannot honour; low_memory=True (and the unset
        # default) document per-chunk divergence, so they stay eligible
        # (GH#64347).  low_memory is not coerced to bool, so falsy non-False
        # values must be rejected too (GH#66327).
        path = tmp_path / "data.csv"
        path.write_text("a,b\n1,2\n", encoding="utf-8")
        monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
        assert not _can_parallelize_csv(path, self._kwds(low_memory=low_memory))
        assert _can_parallelize_csv(path, self._kwds(low_memory=True))
        assert _can_parallelize_csv(path, self._kwds())

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

    def test_one_long_line_does_not_collapse_the_split(self, tmp_path):
        # A line long enough to swallow several split targets must only cost
        # those targets, not every boundary after it.
        path = tmp_path / "data.csv"
        lines = [b"a,b\n"] + [b"1,2\n"] * 20_000
        lines.insert(10, b'1,"' + b"x" * 40_000 + b'"\n')
        path.write_bytes(b"".join(lines))

        offsets = _find_chunk_byte_offsets(str(path), 32, data_start=4)
        # Skipping only the swallowed targets keeps most of the split; giving
        # up at the first one degenerates to 2 chunks.  The long line covers
        # roughly a third of the 32 targets.
        assert len(offsets) - 1 >= 20
        assert all(offsets[i] < offsets[i + 1] for i in range(len(offsets) - 1))


# ---------------------------------------------------------------------------
# _read_csv_parallel  (correctness: parallel == serial)
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestReadCsvParallel:
    """End-to-end correctness: parallel result must match serial result."""

    def _serial_read(self, path, **kwargs):
        return pd.read_csv(path, engine="c", **kwargs)

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
        df = pd.DataFrame({"a": rng.integers(0, 100, n), "b": rng.random(n)})
        df.to_csv(path, index=False)
        kwds = self._base_kwds(path)
        result = self._parallel_read(path, kwds)
        expected = self._serial_read(path)
        tm.assert_frame_equal(result.reset_index(drop=True), expected)

    def test_strings(self, tmp_path):
        path = tmp_path / "data.csv"
        n = 5_000
        rng = np.random.default_rng(1)
        df = pd.DataFrame(
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
        df = pd.DataFrame({"a": vals, "b": rng.integers(0, 5, n)})
        df.to_csv(path, index=False)
        kwds = self._base_kwds(path)
        result = self._parallel_read(path, kwds)
        expected = self._serial_read(path)
        tm.assert_frame_equal(result.reset_index(drop=True), expected)

    def test_header_none(self, tmp_path):
        path = tmp_path / "data.csv"
        n = 5_000
        rng = np.random.default_rng(3)
        df = pd.DataFrame(rng.integers(0, 100, (n, 4)))
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

    def test_long_first_line_keeps_fine_split(self, tmp_path, monkeypatch):
        # The chunk count comes from an estimated row count.  The estimate
        # samples several lines across the data section, so one huge first
        # data line must not collapse the split (or worse, force a serial
        # read via the few-rows bailout) (GH#66275).
        path = tmp_path / "data.csv"
        lines = ["a,b\n", f"0,s{'x' * 50_000}\n"]
        lines += [f"{i},s{i}\n" for i in range(30_000)]
        path.write_text("".join(lines), encoding="utf-8")

        n_targets = []

        def spy(filepath, n_chunks, data_start):
            n_targets.append(n_chunks)
            return _find_chunk_byte_offsets(filepath, n_chunks, data_start)

        monkeypatch.setattr("pandas.io.parsers.readers._find_chunk_byte_offsets", spy)
        kwds = self._base_kwds(path)
        result = _read_csv_parallel(str(path), kwds, 4)
        assert result is not None
        expected = self._serial_read(path)
        tm.assert_frame_equal(result.reset_index(drop=True), expected)
        # 30k rows over the 2000-row chunk floor supports the full 4*3
        # split; an estimate from the first line alone would give 2 chunks.
        assert n_targets[0] >= 10

    def test_few_estimated_rows_returns_none(self, tmp_path):
        # ~300 long rows: even a 2-way split would leave each chunk far
        # below _PARALLEL_MIN_CHUNK_ROWS, and the file is small enough that
        # the per-chunk per-column fixed costs would outweigh the
        # parallelism (GH#66275).
        path = tmp_path / "wide.csv"
        header = ",".join(f"c{i}" for i in range(500)) + "\n"
        row = ",".join(["123456789"] * 500) + "\n"
        path.write_text(header + row * 300, encoding="utf-8")
        kwds = self._base_kwds(path)
        assert _read_csv_parallel(str(path), kwds, 4) is None

    def test_blank_first_data_line_returns_none(self, tmp_path):
        # header=None with a blank first data line: name inference parses a
        # single blank line and raises EmptyDataError, which must be turned
        # into a serial fallback rather than propagate (GH#66259).
        path = tmp_path / "blank_first.csv"
        path.write_text(
            "\n" + "".join(f"{i},{i * 2}\n" for i in range(5_000)), encoding="utf-8"
        )
        kwds = self._base_kwds(path, header=None)
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
        df = pd.DataFrame({"a": rng.integers(0, 100, n), "b": rng.integers(0, 100, n)})
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

    serial = pd.read_csv(io.BytesIO(csv_bytes))
    parallel = pd.read_csv(path)  # may take parallel path or fall back

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
    path = tmp_path / "big.csv"
    _make_large_csv(path)
    # Lower the threshold so any file triggers the parallel path.
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)

    # Force the parallel path regardless of the host's CPU allocation (a
    # single usable CPU defaults to serial) so this exercises parallel ==
    # serial everywhere.
    with pd.option_context("mode.max_threads", 4):
        result = pd.read_csv(path)  # auto-selects parallel path for C engine
    expected = pd.read_csv(path, engine="python")
    tm.assert_frame_equal(result, expected)


def test_read_csv_parallel_vs_serial_large_file(tmp_path, monkeypatch):
    """
    For a file that exceeds the threshold, the parallel result equals the
    result from a direct serial C-engine read (with parallelism forced off).
    """
    path = tmp_path / "big.csv"
    _make_large_csv(path)
    # Lower the threshold so any file triggers the parallel path.
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)

    serial = pd.read_csv(path, engine="python")
    # Force the parallel path so this runs even on a single-CPU allocation.
    with pd.option_context("mode.max_threads", 4):
        parallel = pd.read_csv(path, engine="c")
    tm.assert_frame_equal(parallel, serial)


@pytest.mark.parametrize("platform_name", ["linux", "darwin", "win32"])
def test_parallel_on_by_default(tmp_path, monkeypatch, platform_name):
    """The parallel path is on by default on every threaded platform."""
    path = tmp_path / "big.csv"
    _make_large_csv(path)
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
    # Pin the default so the test does not depend on the host's actual core
    # count, or on how many CPUs a container/affinity mask leaves the runner --
    # one usable CPU would make the default serial everywhere.
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 4)
    monkeypatch.setattr(_readers, "available_cpu_count", lambda: None)

    # Stub the parallel reader so this exercises only the platform-gating
    # decision.  Calling the real one would start threads, which fails on
    # no-thread platforms such as Pyodide/WASM.
    calls = []

    def stub(*args, **kwargs):
        calls.append(args)
        return pd.DataFrame()

    monkeypatch.setattr(_readers, "_read_csv_parallel", stub)
    monkeypatch.setattr(_readers.sys, "platform", platform_name)

    pd.read_csv(path)
    assert len(calls) == 1


def test_parallel_default_off_on_wasm(tmp_path, monkeypatch):
    """Emscripten cannot spawn threads, so it stays serial."""
    path = tmp_path / "big.csv"
    _make_large_csv(path)
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 4)
    monkeypatch.setattr(_readers, "available_cpu_count", lambda: None)

    calls = []

    def stub(*args, **kwargs):
        calls.append(args)
        return pd.DataFrame()

    monkeypatch.setattr(_readers, "_read_csv_parallel", stub)
    monkeypatch.setattr(_readers.sys, "platform", "emscripten")

    pd.read_csv(path)
    assert calls == []


def test_parallel_default_thread_cap(tmp_path, monkeypatch):
    """The default worker count is capped at 4, regardless of core count."""
    path = tmp_path / "big.csv"
    _make_large_csv(path)
    monkeypatch.setattr(_readers, "_PARALLEL_READ_MIN_BYTES", 1)
    monkeypatch.setattr(_readers.sys, "platform", "linux")
    # More cores than the cap: the default should clamp down to 4.
    monkeypatch.setattr(_readers.os, "cpu_count", lambda: 16)
    # Otherwise a CI runner with fewer than 4 usable CPUs clamps below the cap
    # and this test measures the runner, not the cap.
    monkeypatch.setattr(_readers, "available_cpu_count", lambda: None)

    workers = []

    def stub(_path, _kwds, n_workers):
        workers.append(n_workers)
        return pd.DataFrame()

    monkeypatch.setattr(_readers, "_read_csv_parallel", stub)

    pd.read_csv(path)
    assert workers == [4]

    # An explicit mode.max_threads still overrides the cap.
    workers.clear()
    with pd.option_context("mode.max_threads", 8):
        pd.read_csv(path)
    assert workers == [8]


# ---------------------------------------------------------------------------
# _default_n_workers
# ---------------------------------------------------------------------------


class TestDefaultNWorkers:
    """Unit tests for the default worker count."""

    @pytest.mark.parametrize("platform_name", ["linux", "darwin", "win32"])
    @pytest.mark.parametrize(
        "cpu_count, available, expected",
        [
            (2, None, 2),  # unconstrained, below the cap -> logical CPU count
            (16, None, 4),  # cap binds
            (16, 1, 1),  # single-CPU container
            (16, 2, 2),  # cgroup/affinity tighter than the cap
            (16, 8, 4),  # allocation looser than the cap -> cap still binds
            (2, 8, 2),  # allocation looser than the machine
        ],
    )
    def test_combines_cap_and_allocation(
        self, monkeypatch, cpu_count, available, expected, platform_name
    ):
        # Default = min(logical CPUs, available CPUs, _MAX_DEFAULT_WORKERS),
        # on every threaded platform.
        assert _readers._MAX_DEFAULT_WORKERS == 4
        monkeypatch.setattr(_readers.sys, "platform", platform_name)
        monkeypatch.setattr(_readers.os, "cpu_count", lambda: cpu_count)
        monkeypatch.setattr(_readers, "available_cpu_count", lambda: available)
        with pd.option_context("mode.max_threads", None):
            assert _default_n_workers() == expected

    def test_wasm_is_serial(self, monkeypatch):
        monkeypatch.setattr(_readers.sys, "platform", "emscripten")
        # WASM stays serial even if a worker count is requested explicitly.
        with pd.option_context("mode.max_threads", 8):
            assert _default_n_workers() == 1

    def test_max_threads_exceeds_cap(self, monkeypatch):
        # _MAX_DEFAULT_WORKERS and the availability clamp bound the *default*
        # only.  The mode.max_threads docs promise an explicit setting still
        # wins.
        monkeypatch.setattr(_readers.sys, "platform", "linux")
        monkeypatch.setattr(_readers.os, "cpu_count", lambda: 2)
        monkeypatch.setattr(_readers, "available_cpu_count", lambda: 1)
        with pd.option_context("mode.max_threads", 32):
            assert _default_n_workers() == 32

    @pytest.mark.parametrize("platform_name", ["linux", "win32"])
    def test_max_threads_wins(self, monkeypatch, platform_name):
        monkeypatch.setattr(_readers.sys, "platform", platform_name)
        with pd.option_context("mode.max_threads", 3):
            assert _default_n_workers() == 3


# ---------------------------------------------------------------------------
# Regression tests: inputs the parallel path must hand back to serial
# (all results must be identical to a serial read of the same bytes)
# ---------------------------------------------------------------------------


def _read_forced_parallel(path, monkeypatch, **kwargs):
    """read_csv with the parallel path force-enabled on every platform."""
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)
    with pd.option_context("mode.max_threads", 4):
        return pd.read_csv(path, **kwargs)


def test_parallel_latin1_matches_serial(tmp_path, monkeypatch):
    # Non-UTF-8 encodings must take the serial path; the C tokenizer decodes
    # raw chunk bytes as UTF-8 (GH#64347).
    raw = b"col1,col2\n" + b"".join(
        f"se\xf1or{i},{i}\n".encode("latin-1") for i in range(500)
    )
    path = tmp_path / "latin.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, encoding="latin-1")
    expected = pd.read_csv(io.BytesIO(raw), encoding="latin-1")
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
        pd.read_csv(io.BytesIO(raw), encoding="ascii")


def test_parallel_dup_names_dict_dtype_matches_serial(tmp_path, monkeypatch):
    # A dict dtype keyed on a duplicated header name is applied by serial to
    # every de-duplicated column ("a" and "a.1"), but the workers receive the
    # already de-duplicated names and would match only "a".  The parallel path
    # must recover the raw header names and hand back to serial (GH#64347).
    raw = b"a,a\n" + b"".join(f"{i},{i + 1}\n".encode() for i in range(500))
    path = tmp_path / "dup.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, dtype={"a": str})
    expected = pd.read_csv(io.BytesIO(raw), dtype={"a": str})
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
    expected = pd.read_csv(io.BytesIO(raw), comment="#")
    tm.assert_frame_equal(result, expected)
    assert result["col1"].dtype == np.int64


def test_parallel_blank_line_before_header_matches_serial(tmp_path, monkeypatch):
    # Same as the comment case, but triggered by a plain blank line (GH#64347).
    raw = b"\ncol1,col2\n" + b"".join(f"{i},{i * 2}\n".encode() for i in range(500))
    path = tmp_path / "blank.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch)
    expected = pd.read_csv(io.BytesIO(raw))
    tm.assert_frame_equal(result, expected)
    assert result["col1"].dtype == np.int64


@pytest.mark.parametrize(
    "preamble,skiprows",
    [
        (b"\n", 0),
        (b"   \n", 0),
        (b"junk1\njunk2\n\n", 2),
    ],
)
def test_parallel_blank_first_data_line_matches_serial(
    tmp_path, monkeypatch, preamble, skiprows
):
    # header=None puts no header line in the preamble, so a blank first data
    # line ends up as the single line name inference parses, and that raises
    # EmptyDataError where the serial path just skips it (GH#66259).  A
    # whitespace-only line must not tokenize as a single field either: that
    # would infer one column and silently read the rest as an implicit index.
    raw = preamble + b"".join(f"{i},{i * 2}\n".encode() for i in range(500))
    path = tmp_path / "blank_first.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, header=None, skiprows=skiprows)
    expected = pd.read_csv(io.BytesIO(raw), header=None, skiprows=skiprows)
    tm.assert_frame_equal(result, expected)
    assert result.shape == (500, 2)


@pytest.mark.skipif(WASM, reason="WASM stays serial, so the spy sees no call")
def test_parallel_all_blank_lines_raises(tmp_path, monkeypatch):
    # The fallback for a blank first data line must not mask a file that has
    # no columns at all: serial raises, so the parallel path must too
    # (GH#66259).
    path = tmp_path / "blank_only.csv"
    path.write_bytes(b"\n" * 500)

    returned = []

    def spy(filepath, kwds, n_workers):
        result = _read_csv_parallel(filepath, kwds, n_workers)
        returned.append(result)
        return result

    monkeypatch.setattr("pandas.io.parsers.readers._read_csv_parallel", spy)

    msg = "No columns to parse from file"
    with pytest.raises(EmptyDataError, match=msg):
        _read_forced_parallel(path, monkeypatch, header=None)
    # The serial read raises the same error, so without this the test would
    # pass even if the parallel path propagated it directly.
    assert returned == [None], "the parallel path did not fall back"


def test_parallel_single_long_line(tmp_path, monkeypatch):
    # A data section with no interior newlines cannot be split; this must
    # fall back to serial instead of raising (GH#64347).
    raw = b"col1,col2\n" + b"a" * 5_000 + b",b"
    path = tmp_path / "longline.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch)
    expected = pd.read_csv(io.BytesIO(raw))
    tm.assert_frame_equal(result, expected)


@pytest.mark.skipif(WASM, reason="WASM cannot spawn threads, so no split happens")
def test_parallel_long_line_keeps_the_split(tmp_path, monkeypatch):
    # A line covering the middle of the file swallows the interior split
    # targets.  The boundaries after it must survive: the results are correct
    # either way, so a collapse is invisible except as lost parallelism, which
    # is what this asserts (GH#66392).
    path = tmp_path / "data.csv"
    blob = "1," + "9" * 39_995 + "\n"  # 40 KB, covering 20-60% of the data
    path.write_text(
        "a,b\n" + "1,2\n" * 5_000 + blob + "1,2\n" * 10_000, encoding="utf-8"
    )

    chunk_counts = []

    def spy(filepath, n_chunks, data_start):
        offsets = _find_chunk_byte_offsets(filepath, n_chunks, data_start)
        chunk_counts.append(len(offsets) - 1)
        return offsets

    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)
    monkeypatch.setattr("pandas.io.parsers.readers._find_chunk_byte_offsets", spy)

    with pd.option_context("mode.max_threads", 1):
        expected = pd.read_csv(path)
    with pd.option_context("mode.max_threads", 4):
        result = pd.read_csv(path)
    tm.assert_frame_equal(result, expected)
    # Two of the interior targets land inside the long line; giving up at the
    # first leaves 2 chunks no matter how many workers were asked for.
    assert chunk_counts, "the parallel path did not run"
    assert chunk_counts[0] >= 3


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
    expected = pd.read_csv(io.BytesIO(raw))
    tm.assert_frame_equal(result, expected)
    # serial semantics: the whole column stays as strings
    assert result.loc[0, "col1"] == "0"


def test_parallel_deferred_strings_pyarrow_backend(tmp_path, monkeypatch):
    # dtype_backend="pyarrow" builds pa.string() (int32 offset) pending columns
    # rather than the large_string ones the default string dtype uses; the
    # gather must reconcile them to ArrowDtype and match a serial read
    # (GH#66277).
    pytest.importorskip("pyarrow")
    path = tmp_path / "data.csv"
    rows = "\n".join(f"s{i % 13},{i}" for i in range(5000))
    path.write_text("a,b\n" + rows + "\n", encoding="utf-8")
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)

    with pd.option_context("mode.max_threads", 1):
        serial = pd.read_csv(path, dtype_backend="pyarrow")
    with pd.option_context("mode.max_threads", 4):
        parallel = pd.read_csv(path, dtype_backend="pyarrow")
    tm.assert_frame_equal(parallel, serial)


def test_parallel_deferred_strings_token_width_tiers(tmp_path, monkeypatch):
    # The deferred string path fills each chunk's data buffer with the same
    # fixed-width token copy the serial path uses, so widths straddling that
    # copy's 16- and 32-byte tiers have to survive a chunked read too, where
    # each worker sizes and fills its own chunk's buffer (GH#66756).
    pytest.importorskip("pyarrow")
    widths = [1, 2, 15, 16, 17, 31, 32, 33, 64]
    values = [
        "".join(chr(ord("a") + pos % 26) for pos in range(width)) for width in widths
    ]
    path = tmp_path / "tiers.csv"
    rows = "".join(f"{value},{value.upper()}\n" for value in values * 400)
    path.write_text("a,b\n" + rows, encoding="utf-8")
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)

    # Both options are pinned rather than inherited: the fast path needs
    # infer_string on *and* pyarrow storage, and under either default flipped
    # (PANDAS_FUTURE_INFER_STRING=0, mode.string_storage="python") the columns
    # would come back from the object path instead; the dtype check below
    # fails loudly if that happens.
    with pd.option_context(
        "future.infer_string",
        True,
        "mode.string_storage",
        "pyarrow",
        "mode.max_threads",
        1,
    ):
        serial = pd.read_csv(path)
    with pd.option_context(
        "future.infer_string",
        True,
        "mode.string_storage",
        "pyarrow",
        "mode.max_threads",
        4,
    ):
        parallel = pd.read_csv(path)
    tm.assert_frame_equal(parallel, serial)
    assert serial["a"].dtype == pd.StringDtype("pyarrow", na_value=np.nan)
    assert serial["a"].tolist() == values * 400


def test_parallel_deferred_strings_mixed_chunk_dtypes(tmp_path, monkeypatch):
    # One chunk of a column parses as strings while the rest parse as int64,
    # so the gather sees pending string columns alongside numeric ones and must
    # fall back to a serial re-read rather than mixing the two (GH#66277).
    chunk = "".join(f"{i},{i}\n" for i in range(2000))
    raw = "col1,col2\n" + chunk + "text,7\n" + chunk
    path = tmp_path / "mixed_str.csv"
    path.write_text(raw, encoding="utf-8")
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)

    with pd.option_context("mode.max_threads", 1):
        serial = pd.read_csv(path)
    with pd.option_context("mode.max_threads", 4):
        parallel = pd.read_csv(path)
    tm.assert_frame_equal(parallel, serial)


@pytest.mark.parametrize("low_memory", [False, 0, np.False_])
def test_parallel_low_memory_false_matches_serial(tmp_path, monkeypatch, low_memory):
    # A falsy low_memory promises whole-file type inference.  A chunk holding
    # only NA tokens and a uint64-range int hits the C parser's na_filter=0
    # uint64-conflict path (gh-14983) and keeps "nan" as a literal string, while
    # the serial whole-file pass masks it to NaN.  Both chunk results are object
    # dtype, so the per-column dtype-consistency check cannot catch it - the
    # parallel path must step aside (GH#64347).  low_memory is not coerced to
    # bool and the serial path branches on truthiness, so falsy non-False values
    # must step aside too (GH#66327).
    raw = b"col\n" + b"hello\n" * 500 + (b"nan\n" + b"9223372036854775808\n") * 3000
    path = tmp_path / "uint64.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, low_memory=low_memory)
    expected = pd.read_csv(io.BytesIO(raw), low_memory=low_memory)
    tm.assert_frame_equal(result, expected)
    # serial semantics: every "nan" token is masked to NaN, none survive as text
    assert (result["col"] == "nan").sum() == 0
    assert result["col"].isna().sum() == 3000


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"na_values": ["MISSING"]},
        {"na_values": {"col": ["MISSING"]}},
        {"dtype_backend": "numpy_nullable"},
    ],
)
def test_parallel_uint64_na_conflict_matches_serial(tmp_path, monkeypatch, kwargs):
    # A chunk holding only NA tokens and uint64-range ints converts to no
    # numeric dtype and is emitted with its NA tokens left as literal strings
    # (gh-14983).  A chunk that also sees ordinary text takes the object branch
    # and masks them instead, and both results are str dtype, so the gather's
    # dtype reconciliation cannot tell them apart - the parallel path must step
    # aside rather than return values a serial read never produces (GH#66259).
    token = "MISSING" if kwargs.get("na_values") else "nan"
    body = f"{token}\n9223372036854775808\n" * 3000
    path = tmp_path / "uint64_default_lm.csv"
    path.write_text("col\n" + "text\n" * 500 + body, encoding="utf-8")

    result = _read_forced_parallel(path, monkeypatch, **kwargs)
    with pd.option_context("mode.max_threads", 1):
        expected = pd.read_csv(path, **kwargs)
    tm.assert_frame_equal(result, expected)


@pytest.mark.skipif(WASM, reason="WASM stays serial, so the spy sees no call")
def test_parallel_uint64_na_conflict_na_filter_false(tmp_path, monkeypatch):
    # With na_filter=False the literal NA tokens are what a serial read gives
    # too, so such a file keeps using the parallel path.  GH#66259
    body = "nan\n9223372036854775808\n" * 3000
    path = tmp_path / "uint64_no_filter.csv"
    path.write_text("col\n" + "text\n" * 500 + body, encoding="utf-8")

    fell_back = []
    orig = _readers._read_csv_parallel

    def spy(*args, **kwargs):
        result = orig(*args, **kwargs)
        fell_back.append(result is None)
        return result

    monkeypatch.setattr(_readers, "_read_csv_parallel", spy)

    result = _read_forced_parallel(path, monkeypatch, na_filter=False)
    assert fell_back == [False]
    with pd.option_context("mode.max_threads", 1):
        expected = pd.read_csv(path, na_filter=False)
    tm.assert_frame_equal(result, expected)


def test_parallel_quoted_newline_in_header_matches_serial(tmp_path, monkeypatch):
    # A quoted embedded newline in the header shifts data_start mid-header;
    # without a guard, chunk 0 starts inside the header and silently injects
    # a garbage row (GH#64347).  All-string data so the per-column dtype
    # consistency check cannot catch it.
    raw = b'"col\n1",col2\n' + b"".join(f"x{i},y{i}\n".encode() for i in range(500))
    path = tmp_path / "hdr.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch)
    expected = pd.read_csv(io.BytesIO(raw))
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
        expected = pd.read_csv(io.BytesIO(raw), on_bad_lines="warn")
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
    expected = pd.read_csv(io.BytesIO(raw))
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
    expected = pd.read_csv(io.BytesIO(raw), dialect=EscDialect)
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
        expected = pd.read_csv(io.BytesIO(raw), sep=";;")
    tm.assert_frame_equal(result, expected)


def _write_with_line_at_chunk_start(path, replacement: bytes, monkeypatch) -> int:
    """Write a fixed-width CSV, then overwrite the line at the second chunk
    boundary with *replacement* (same byte length, so offsets stay valid).

    Returns the boundary, for :func:`_assert_chunk_starts_at` to check against
    the split the read actually plans - hard-coding a count here that the read
    does not use would put the line in the middle of a chunk instead.
    """
    rows = [f"{i:06d},{i * 2:06d}" for i in range(4000)]
    # write bytes so line endings stay "\n" on every platform - the byte-offset
    # math below assumes single-byte terminators (Windows text mode adds "\r")
    path.write_bytes(("a,b\n" + "\n".join(rows) + "\n").encode("utf-8"))
    # 4000 rows is only two chunks' worth at the row floor, and the split has to
    # be finer than that for a boundary to land mid-file.  _read_forced_parallel
    # uses 4 workers, so this yields the 4 * 3 oversubscription.
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_MIN_CHUNK_ROWS", 1)
    data_start = _find_data_start_offset(str(path), 0, 0)
    offsets = _find_chunk_byte_offsets(str(path), 12, data_start)
    boundary = offsets[1]
    raw = path.read_bytes()
    line_end = raw.index(b"\n", boundary)
    assert line_end - boundary == len(replacement)
    with open(path, "r+b") as fd:
        fd.seek(boundary)
        fd.write(replacement)
    return boundary


def _spy_on_chunk_offsets(monkeypatch) -> list:
    """Record the byte offsets each planned split actually uses."""
    seen: list = []
    real = _find_chunk_byte_offsets

    def spy(filepath, n_chunks, data_start):
        offsets = real(filepath, n_chunks, data_start)
        seen.append(offsets)
        return offsets

    monkeypatch.setattr("pandas.io.parsers.readers._find_chunk_byte_offsets", spy)
    return seen


def _assert_chunk_starts_at(seen: list, boundary: int) -> None:
    # Without this the tests below still pass on a coarser split - a ragged line
    # raises (and a skipped one is skipped) wherever it sits - while no longer
    # placing it at the chunk start that is the point of the fixture.
    assert seen, "the parallel path did not run"
    assert boundary in seen[0], f"{boundary} is not a chunk start in {seen[0]}"


@pytest.mark.skipif(WASM, reason="WASM stays serial, so the spy sees no call")
def test_parallel_ragged_line_at_chunk_start_raises(tmp_path, monkeypatch):
    # A line with extra fields sitting exactly at a chunk boundary: the chunk
    # worker's first line used to be exempt from the field-count check, so
    # parallel silently truncated the row while serial raised (GH#64347).
    path = tmp_path / "ragged.csv"
    boundary = _write_with_line_at_chunk_start(path, b"11,22,33,4444", monkeypatch)
    seen = _spy_on_chunk_offsets(monkeypatch)

    with pytest.raises(ParserError, match="Expected 2 fields"):
        _read_forced_parallel(path, monkeypatch)
    _assert_chunk_starts_at(seen, boundary)


@pytest.mark.skipif(WASM, reason="WASM stays serial, so the spy sees no call")
def test_parallel_ragged_line_at_chunk_start_skip_matches_serial(tmp_path, monkeypatch):
    # Same layout with on_bad_lines="skip": the chunk worker must skip the
    # bad line just like serial, not keep a truncated version of it.
    path = tmp_path / "ragged.csv"
    boundary = _write_with_line_at_chunk_start(path, b"11,22,33,4444", monkeypatch)
    seen = _spy_on_chunk_offsets(monkeypatch)

    with pd.option_context("mode.max_threads", 1):
        expected = pd.read_csv(path, on_bad_lines="skip")
    result = _read_forced_parallel(path, monkeypatch, on_bad_lines="skip")
    tm.assert_frame_equal(result, expected)
    assert len(result) == 3999
    _assert_chunk_starts_at(seen, boundary)


def test_parallel_bom_bytes_mid_file_match_serial(tmp_path, monkeypatch):
    # Data lines beginning with the UTF-8 BOM byte sequence: each chunk
    # worker's fresh parser used to strip it from the chunk's first line,
    # while serial strips only at byte 0 of the file (GH#64347).
    lines = ["a,b"] + ["\ufeff" + f"x{i},{i}" for i in range(4000)]
    path = tmp_path / "bom.csv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pd.option_context("mode.max_threads", 1):
        expected = pd.read_csv(path)
    result = _read_forced_parallel(path, monkeypatch)
    tm.assert_frame_equal(result, expected)
    assert (result["a"].str[0] == "\ufeff").all()


def test_parallel_bom_at_file_start_header_none(tmp_path, monkeypatch):
    # With header=None and no skiprows the first chunk starts at byte 0 and
    # must strip a leading BOM exactly like the serial parser does.
    body = "".join(f"{i},{i * 2}\n" for i in range(4000))
    path = tmp_path / "bom0.csv"
    path.write_bytes(b"\xef\xbb\xbf" + body.encode())

    with pd.option_context("mode.max_threads", 1):
        expected = pd.read_csv(path, header=None, names=["a", "b"])
    result = _read_forced_parallel(path, monkeypatch, header=None, names=["a", "b"])
    tm.assert_frame_equal(result, expected)
    # the stripped BOM means the first cell parses as an integer
    assert result["a"].iloc[0] == 0


def test_parallel_string_dtype_python_storage(tmp_path, monkeypatch):
    # Without the pyarrow string fast path, the gathered object column must
    # get the same string-dtype inference a serial read applies via the
    # DataFrame constructor (GH#66275).
    path = tmp_path / "data.csv"
    rows = "\n".join(f"{i},s{i % 7}" for i in range(5000))
    path.write_text("a,b\n" + rows + "\n", encoding="utf-8")
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)

    with pd.option_context("mode.string_storage", "python"):
        with pd.option_context("mode.max_threads", 1):
            serial = pd.read_csv(path)
        with pd.option_context("mode.max_threads", 4):
            parallel = pd.read_csv(path)
    tm.assert_frame_equal(parallel, serial)


def test_parallel_blank_line_run(tmp_path, monkeypatch):
    # A run of blank lines spanning entire chunks parses those chunks to zero
    # rows; the reused worker parsers must not leak StopIteration out of
    # read_csv (GH#66275).
    path = tmp_path / "data.csv"
    rows = "\n".join(f"{i},{i * 2}" for i in range(100))
    path.write_text(
        "a,b\n" + rows + "\n" + "\n" * 20000 + rows + "\n", encoding="utf-8"
    )
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)

    with pd.option_context("mode.max_threads", 1):
        expected = pd.read_csv(path)
    with pd.option_context("mode.max_threads", 4):
        result = pd.read_csv(path)
    tm.assert_frame_equal(result, expected)


def test_parallel_bad_line_run_skip(tmp_path, monkeypatch):
    # Same for chunks whose lines are all skipped via on_bad_lines="skip"
    # (GH#66275).
    path = tmp_path / "data.csv"
    rows = "\n".join(f"{i},{i * 2}" for i in range(100))
    bad = "\n".join("x,y,z,w,v" for _ in range(5000))
    path.write_text("a,b\n" + rows + "\n" + bad + "\n" + rows + "\n", encoding="utf-8")
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)
    # The good rows sit at both ends of the file, so at the two chunks the row
    # floor allows here every chunk has something to parse and the all-skipped
    # chunk this is about never happens.
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_MIN_CHUNK_ROWS", 1)

    with pd.option_context("mode.max_threads", 1):
        expected = pd.read_csv(path, on_bad_lines="skip")
    with pd.option_context("mode.max_threads", 4):
        result = pd.read_csv(path, on_bad_lines="skip")
    tm.assert_frame_equal(result, expected)


def test_parallel_embedded_nul_boolean_column(tmp_path, monkeypatch):
    # GH#66524: a NUL-bearing field must not match a true/false value on its
    # pre-NUL prefix in the parallel path either, which feeds the tokenizer
    # through load_buffer() rather than parser_buffer_bytes().
    path = tmp_path / "nul_bool.csv"
    path.write_bytes(b"a\n" + b"True\n" * 200 + b'"True\x00xyz"\n' + b"False\n" * 200)
    monkeypatch.setattr("pandas.io.parsers.readers._PARALLEL_READ_MIN_BYTES", 1)

    with pd.option_context("mode.max_threads", 4):
        result = pd.read_csv(path)
    expected = pd.read_csv(path, engine="python")
    tm.assert_frame_equal(result, expected)
    assert result["a"][200] == "True\x00xyz"


def _warnings_from(func, *args, **kwargs):
    """Run *func*, returning every warning it raises (no de-duplication)."""
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        result = func(*args, **kwargs)
    return result, recorded


def _track_parallel(monkeypatch) -> list[str]:
    """Record how each parallel attempt ends: used, declined, or raised."""
    outcomes: list[str] = []
    original = _readers._read_csv_parallel

    def tracked(*args, **kwargs):
        try:
            result = original(*args, **kwargs)
        except Exception:
            outcomes.append("raised")
            raise
        outcomes.append("used" if result is not None else "declined")
        return result

    monkeypatch.setattr(_readers, "_read_csv_parallel", tracked)
    return outcomes


def _converter_dtype_warning(name: str) -> str:
    return (
        f"Both a converter and dtype were specified for column {name} - "
        "only the converter will be used."
    )


# one column, then two: de-duplicating the collected warnings must not
# collapse the distinct ones
@pytest.mark.parametrize("names", [["col1"], ["col1", "col2"]])
@pytest.mark.skipif(WASM, reason="WASM stays serial, so no chunk repeats the warning")
def test_parallel_converter_dtype_warns_once(tmp_path, monkeypatch, names):
    # The converter+dtype ParserWarning was raised once per chunk, from a pool
    # thread whose stacklevel walk lands in threading internals (GH#66259).
    raw = b"col1,col2\n" + b"".join(f"{i},{i * 2}\n".encode() for i in range(500))
    path = tmp_path / "converter.csv"
    path.write_bytes(raw)
    kwargs = {
        "converters": dict.fromkeys(names, int),
        "dtype": dict.fromkeys(names, "int64"),
    }
    outcomes = _track_parallel(monkeypatch)

    result, recorded = _warnings_from(
        _read_forced_parallel, path, monkeypatch, **kwargs
    )

    expected, _ = _warnings_from(pd.read_csv, io.BytesIO(raw), **kwargs)
    assert outcomes == ["used"]
    tm.assert_frame_equal(result, expected)
    assert [str(warning.message) for warning in recorded] == [
        _converter_dtype_warning(name) for name in names
    ]
    assert recorded[0].category is ParserWarning
    # the caller's frame, not a worker thread's
    assert recorded[0].filename == __file__


@pytest.mark.skipif(WASM, reason="WASM stays serial, so no worker raises")
def test_parallel_worker_exception_still_warns(tmp_path, monkeypatch):
    # An exception the caller does not answer with a serial read must not carry
    # the collected warnings off with it - this read is their only chance to be
    # raised (GH#66259).
    raw = b"col1,col2\n" + b"".join(f"{i},{i * 2}\n".encode() for i in range(1000))
    path = tmp_path / "boom.csv"
    path.write_bytes(raw)

    def boom(value):
        if value == "1400":
            raise RuntimeError("boom")
        return value

    kwargs = {"converters": {"col1": int, "col2": boom}, "dtype": {"col1": "int64"}}
    outcomes = _track_parallel(monkeypatch)
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        with pytest.raises(RuntimeError, match="boom"):
            _read_forced_parallel(path, monkeypatch, **kwargs)

    assert outcomes == ["raised"]
    assert [str(warning.message) for warning in recorded] == [
        _converter_dtype_warning("col1")
    ]


@pytest.mark.skipif(WASM, reason="WASM stays serial, so no worker raises")
def test_parallel_worker_exception_unmaps_the_file(tmp_path, monkeypatch):
    # A worker exception must not leave the file mapped until its traceback is
    # collected: the workers hold memoryview slices of the mapping, so they are
    # closed before it is (GH#66259).
    raw = b"col1,col2\n" + b"".join(f"{i},{i * 2}\n".encode() for i in range(1000))
    path = tmp_path / "boom.csv"
    path.write_bytes(raw)

    mapped: list[mmap.mmap] = []

    class TrackingMmap(mmap.mmap):
        def __init__(self, *args, **kwargs) -> None:
            mapped.append(self)

    monkeypatch.setattr(mmap, "mmap", TrackingMmap)

    def boom(value):
        if value == "1400":
            raise RuntimeError("boom")
        return value

    with pytest.raises(RuntimeError, match="boom"):
        _read_forced_parallel(path, monkeypatch, converters={"col2": boom})

    assert len(mapped) == 1
    assert mapped[0].closed


@pytest.mark.skipif(WASM, reason="WASM stays serial, so there is no attempt to decline")
def test_parallel_fallback_does_not_repeat_c_layer_warning(tmp_path, monkeypatch):
    # A parallel attempt that hands back to serial must not leave its own copy
    # of a warning behind: the name-inference read hits the same converter+dtype
    # warning the replacing serial read raises (GH#66259).
    raw = (
        b"col1,col2\n"
        + b"".join(f"{i},{i * 2}\n".encode() for i in range(950))
        # col2 is str in the last chunk only, so the per-chunk dtypes disagree
        # and the parallel path declines after inferring the column names
        + b"".join(f"{i},x{i}\n".encode() for i in range(50))
    )
    path = tmp_path / "fallback.csv"
    path.write_bytes(raw)
    kwargs = {"converters": {"col1": int}, "dtype": {"col1": "int64"}}
    outcomes = _track_parallel(monkeypatch)

    result, recorded = _warnings_from(
        _read_forced_parallel, path, monkeypatch, **kwargs
    )

    expected, _ = _warnings_from(pd.read_csv, io.BytesIO(raw), **kwargs)
    assert outcomes == ["declined"]
    tm.assert_frame_equal(result, expected)
    assert [str(warning.message) for warning in recorded] == [
        _converter_dtype_warning("col1")
    ]


@pytest.mark.skipif(
    WASM, reason="WASM stays serial, so there is no attempt to fall back"
)
def test_parallel_fallback_does_not_repeat_python_layer_warning(tmp_path, monkeypatch):
    # Same for a warning raised by the Python layer rather than the C parser:
    # index_col=False with more fields than header names warns in the
    # name-inference read, and again in the serial read that a worker's
    # ParserError falls back to (GH#66259).
    raw = b"col1,col2\n" + b"".join(
        f"{i},{i * 2},{i * 3}\n".encode() for i in range(500)
    )
    path = tmp_path / "mismatch.csv"
    path.write_bytes(raw)
    outcomes = _track_parallel(monkeypatch)

    result, recorded = _warnings_from(
        _read_forced_parallel, path, monkeypatch, index_col=False
    )

    expected, _ = _warnings_from(pd.read_csv, io.BytesIO(raw), index_col=False)
    assert outcomes == ["raised"]
    tm.assert_frame_equal(result, expected)
    assert [str(warning.message) for warning in recorded] == [
        "Length of header or names does not match length of data. This leads "
        "to a loss of data with index_col=False."
    ]


@pytest.mark.parametrize("kwargs", [{"quotechar": "»"}, {"sep": "§"}])
def test_multibyte_sep_or_quotechar_warns_once(tmp_path, monkeypatch, kwargs):
    # Both force TextFileReader onto the python engine, which warns while
    # falling back.  The parallel path has to decline them up front: its
    # name-inference read would raise that warning before it can notice the
    # engine it got, and the serial read then raises it again (GH#66259).
    raw = b"col1,col2\n" + b"".join(f"{i},{i * 2}\n".encode() for i in range(500))
    path = tmp_path / "multibyte.csv"
    path.write_bytes(raw)
    outcomes = _track_parallel(monkeypatch)

    result, recorded = _warnings_from(
        _read_forced_parallel, path, monkeypatch, **kwargs
    )

    expected, _ = _warnings_from(pd.read_csv, io.BytesIO(raw), **kwargs)
    assert outcomes == []
    tm.assert_frame_equal(result, expected)
    assert len(recorded) == 1
    assert recorded[0].category is ParserWarning
    assert "Falling back to the 'python' engine" in str(recorded[0].message)


# 2**63 - 1 gathers through pyarrow's checked int-to-double cast; 2**64 does not
# even fit an integer chunk, so it fails earlier, in a worker's own conversion
@pytest.mark.parametrize(
    ("big", "outcome"),
    [("9223372036854775807", "declined"), ("18446744073709551616", "raised")],
)
@pytest.mark.skipif(WASM, reason="WASM stays serial, so the chunks never disagree")
def test_parallel_pyarrow_huge_int_chunk_matches_serial(
    tmp_path, monkeypatch, big, outcome
):
    # A chunk holding only huge integers converts as an integer column, while
    # the whole-file serial read sees the trailing float too and infers double
    # (or leaves the column a string).  The parallel path used to propagate the
    # resulting ArrowInvalid / OverflowError instead of re-reading serially
    # (GH#66259).
    pytest.importorskip("pyarrow")
    raw = b"col1\n" + f"{big}\n".encode() * 500 + b"1.5\n"
    path = tmp_path / "huge_int.csv"
    path.write_bytes(raw)
    outcomes = _track_parallel(monkeypatch)

    result = _read_forced_parallel(path, monkeypatch, dtype_backend="pyarrow")

    expected = pd.read_csv(io.BytesIO(raw), dtype_backend="pyarrow")
    assert outcomes == [outcome]
    tm.assert_frame_equal(result, expected)


@pytest.mark.skipif(WASM, reason="WASM stays serial, so nothing infers per chunk")
def test_parallel_pyarrow_huge_int_column_stays_parallel(tmp_path, monkeypatch):
    # The guard above must not cost the parallel path columns it handles fine:
    # with no float row every chunk agrees on int64 and nothing falls back
    # (GH#66259).
    pytest.importorskip("pyarrow")
    raw = b"col1\n" + b"9223372036854775807\n" * 500
    path = tmp_path / "huge_int_only.csv"
    path.write_bytes(raw)
    outcomes = _track_parallel(monkeypatch)

    result = _read_forced_parallel(path, monkeypatch, dtype_backend="pyarrow")

    expected = pd.read_csv(io.BytesIO(raw), dtype_backend="pyarrow")
    assert outcomes == ["used"]
    tm.assert_frame_equal(result, expected)


@pytest.mark.skipif(WASM, reason="WASM stays serial, so no sample line is parsed")
def test_parallel_huge_int_first_line_matches_serial(tmp_path, monkeypatch):
    # The name-inference read parses one data line only to learn the column
    # names, so a value that converts on its own but not for the whole column
    # must not raise there (GH#66259).
    pytest.importorskip("pyarrow")
    raw = b"col1\n18446744073709551616\n" + b"".join(
        f"{i}.5\n".encode() for i in range(500)
    )
    path = tmp_path / "huge_first.csv"
    path.write_bytes(raw)

    result = _read_forced_parallel(path, monkeypatch, dtype_backend="pyarrow")

    tm.assert_frame_equal(result, pd.read_csv(io.BytesIO(raw), dtype_backend="pyarrow"))
