"""
Module for formatting output data into CSV files.
"""

from __future__ import annotations

import codecs
from collections.abc import (
    Hashable,
    Iterable,
    Iterator,
    Sequence,
)
import csv as csvlib
import io
import os
from typing import (
    IO,
    TYPE_CHECKING,
    Any,
    AnyStr,
    cast,
)
import warnings

import numpy as np

from pandas._libs import writers as libwriters
from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import cache_readonly
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
    DatetimeTZDtype,
)
from pandas.core.dtypes.generic import (
    ABCDatetimeIndex,
    ABCIndex,
    ABCMultiIndex,
    ABCPeriodIndex,
)
from pandas.core.dtypes.missing import notna

from pandas.core.indexes.api import Index

from pandas.io.common import get_handle

if TYPE_CHECKING:
    from pandas._typing import (
        CompressionOptions,
        DtypeObj,
        FilePath,
        FloatFormatType,
        IndexLabel,
        SequenceNotStr,
        StorageOptions,
        WriteBuffer,
        npt,
    )

    from pandas import DataFrame

    from pandas.io.formats.format import DataFrameFormatter


_DEFAULT_CHUNKSIZE_CELLS = 100_000

# Objects that expect str (not bytes) to be written to them, even though
# they may wrap/proxy an underlying binary handle (e.g. codecs.StreamWriter
# wraps a binary file but itself only accepts str). The pyarrow CSV writer
# only ever writes bytes, so none of these are usable as its destination.
_TEXT_LIKE_CLASSES = (io.TextIOBase, codecs.StreamWriter, codecs.StreamReaderWriter)


class CSVFormatter:
    cols: npt.NDArray[np.object_]

    def __init__(
        self,
        formatter: DataFrameFormatter,
        path_or_buf: FilePath | WriteBuffer[str] | WriteBuffer[bytes] | None = None,
        engine: str = "auto",
        sep: str = ",",
        cols: Sequence[Hashable] | None = None,
        index_label: IndexLabel | None = None,
        mode: str | None = None,
        encoding: str | None = None,
        errors: str = "strict",
        compression: CompressionOptions = "infer",
        quoting: int | None = None,
        lineterminator: str | None = "\n",
        chunksize: int | None = None,
        quotechar: str | None = '"',
        date_format: str | None = None,
        doublequote: bool = True,
        escapechar: str | None = None,
        storage_options: StorageOptions | None = None,
    ) -> None:
        self.fmt = formatter

        self.obj = self.fmt.frame

        self.filepath_or_buffer = path_or_buf
        self.encoding = encoding
        self.compression: CompressionOptions = compression
        self.mode = mode
        self.storage_options = storage_options

        self.sep = sep
        self.index_label = self._initialize_index_label(index_label)
        self.errors = errors
        self.quoting = quoting or csvlib.QUOTE_MINIMAL
        self.doublequote = doublequote
        self.escapechar = escapechar
        self.quotechar = self._initialize_quotechar(quotechar)
        # Remember whether the caller actually requested a lineterminator,
        # separately from the resolved value below: pyarrow always writes
        # "\n" regardless of platform, which is fine to use silently when
        # nothing specific was requested, but not when the caller
        # explicitly asked for something else (see
        # _pyarrow_option_incompatibility).
        self._lineterminator_requested = lineterminator
        self.lineterminator = lineterminator or os.linesep
        self.date_format = date_format
        self.cols = self._initialize_columns(cols)
        self.chunksize = self._initialize_chunksize(chunksize)

        # only set once engine resolution succeeds with engine="pyarrow"
        self._arrow_table: Any | None = None
        self.engine = self._resolve_engine(engine)

    @property
    def na_rep(self) -> str:
        return self.fmt.na_rep

    @property
    def float_format(self) -> FloatFormatType | None:
        return self.fmt.float_format

    @property
    def decimal(self) -> str:
        return self.fmt.decimal

    @property
    def header(self) -> bool | SequenceNotStr[str]:
        return self.fmt.header

    @property
    def index(self) -> bool:
        return self.fmt.index

    def _initialize_index_label(self, index_label: IndexLabel | None) -> IndexLabel:
        if index_label is not False:
            if index_label is None:
                return self._get_index_label_from_obj()
            elif not isinstance(index_label, (list, tuple, np.ndarray, ABCIndex)):
                # given a string for a DF with Index
                return [index_label]
        return index_label

    def _get_index_label_from_obj(self) -> Sequence[Hashable]:
        if isinstance(self.obj.index, ABCMultiIndex):
            return self._get_index_label_multiindex()
        else:
            return self._get_index_label_flat()

    def _get_index_label_multiindex(self) -> Sequence[Hashable]:
        return [name or "" for name in self.obj.index.names]

    def _get_index_label_flat(self) -> Sequence[Hashable]:
        index_label = self.obj.index.name
        return [""] if index_label is None else [index_label]

    def _initialize_quotechar(self, quotechar: str | None) -> str | None:
        if self.quoting != csvlib.QUOTE_NONE or self.escapechar is not None:
            # prevents crash in _csv
            return quotechar
        return None

    @property
    def has_mi_columns(self) -> bool:
        return bool(isinstance(self.obj.columns, ABCMultiIndex))

    def _initialize_columns(
        self, cols: Iterable[Hashable] | None
    ) -> npt.NDArray[np.object_]:
        # validate mi options
        if self.has_mi_columns:
            if cols is not None:
                msg = "cannot specify cols with a MultiIndex on the columns"
                raise TypeError(msg)

        if cols is not None:
            if isinstance(cols, ABCIndex):
                cols = cols._get_values_for_csv(**self._number_format)
            else:
                cols = list(cols)
            self.obj = self.obj.loc[:, cols]

        # update columns to include possible multiplicity of dupes
        # and make sure cols is just a list of labels
        new_cols = self.obj.columns
        return new_cols._get_values_for_csv(**self._number_format)

    def _initialize_chunksize(self, chunksize: int | None) -> int:
        if chunksize is None:
            return (_DEFAULT_CHUNKSIZE_CELLS // (len(self.cols) or 1)) or 1
        return int(chunksize)

    @property
    def _number_format(self) -> dict[str, Any]:
        """Dictionary used for storing number formatting settings."""
        return {
            "na_rep": self.na_rep,
            "float_format": self.float_format,
            "date_format": self.date_format,
            "quoting": self.quoting,
            "decimal": self.decimal,
        }

    @cache_readonly
    def data_index(self) -> Index:
        data_index = self.obj.index
        if (
            isinstance(data_index, (ABCDatetimeIndex, ABCPeriodIndex))
            and self.date_format is not None
        ):
            data_index = Index(
                [x.strftime(self.date_format) if notna(x) else "" for x in data_index]
            )
        elif isinstance(data_index, ABCMultiIndex):
            data_index = data_index.remove_unused_levels()
        return data_index

    @property
    def nlevels(self) -> int:
        if self.index:
            return getattr(self.data_index, "nlevels", 1)
        else:
            return 0

    @property
    def _has_aliases(self) -> bool:
        return isinstance(self.header, (tuple, list, np.ndarray, ABCIndex))

    @property
    def _need_to_save_header(self) -> bool:
        return bool(self._has_aliases or self.header)

    @property
    def write_cols(self) -> SequenceNotStr[Hashable]:
        if self._has_aliases:
            assert not isinstance(self.header, bool)
            if len(self.header) != len(self.cols):
                raise ValueError(
                    f"Writing {len(self.cols)} cols but got {len(self.header)} aliases"
                )
            return self.header
        else:
            # self.cols is an ndarray derived from Index._get_values_for_csv,
            #  so its entries are strings, i.e. hashable
            return cast("SequenceNotStr[Hashable]", self.cols)

    @property
    def encoded_labels(self) -> list[Hashable]:
        encoded_labels: list[Hashable] = []

        if self.index and self.index_label:
            assert isinstance(self.index_label, Sequence)
            encoded_labels = list(self.index_label)

        if not self.has_mi_columns or self._has_aliases:
            encoded_labels += list(self.write_cols)

        return encoded_labels

    def _pyarrow_option_incompatibility(self) -> str | None:
        """
        Return a reason the pyarrow engine cannot honor the requested
        options, or None if the options are compatible. These are options
        for which the pyarrow CSV writer has no equivalent, so silently
        using pyarrow would produce output that quietly ignores the option
        rather than applying it.
        """
        if isinstance(self.filepath_or_buffer, _TEXT_LIKE_CLASSES):
            return "The pyarrow engine can only write to a binary buffer or file path."
        if self.mode is not None and "b" not in self.mode:
            return "The pyarrow engine can only write in binary mode."

        if self.quotechar is not None and self.quotechar != '"':
            return 'The pyarrow engine only supports " as a quotechar.'

        # pyarrow's CSV writer always writes "\n" as its line terminator; it
        # has no option to configure this, and does not follow os.linesep
        # (e.g. it still writes "\n" on Windows). When the caller didn't
        # request a specific lineterminator, silently writing "\n" instead
        # of the python engine's platform-native default is a cosmetic-only
        # difference we accept (like always-quoted headers); but an
        # explicit request for anything other than "\n" cannot be honored.
        if (
            self._lineterminator_requested is not None
            and self._lineterminator_requested != "\n"
        ):
            return "The lineterminator option is not supported with the pyarrow engine."

        unsupported_options = [
            # each pair is (option value, default, option name)
            (self.decimal, ".", "decimal"),
            (self.float_format, None, "float_format"),
            (self.na_rep, "", "na_rep"),
            (self.date_format, None, "date_format"),
            (self.errors, "strict", "errors"),
            (self.escapechar, None, "escapechar"),
            (self.doublequote, True, "doublequote"),
        ]
        for opt_val, default, option in unsupported_options:
            if opt_val != default:
                return f"The {option} option is not supported with the pyarrow engine."

        # pyarrow's CSV writer always writes utf-8, so an explicit request
        # for utf-8 is compatible even though it differs from the sentinel
        # (None) used to detect a non-default `encoding`.
        if self.encoding is not None and self.encoding.replace("-", "").replace(
            "_", ""
        ).lower() not in ("utf8", "u8"):
            return "The encoding option is not supported with the pyarrow engine."

        if self.has_mi_columns:
            return "The pyarrow engine does not support a MultiIndex on the columns."

        return None

    def _iter_written_dtypes(self) -> Iterator[DtypeObj]:
        """Dtypes of every column that will actually be written, incl. the index."""
        if self.index:
            index = self.obj.index
            if isinstance(index, ABCMultiIndex):
                yield from index.dtypes
            else:
                yield index.dtype
        yield from self.obj.dtypes

    @staticmethod
    def _is_tz_aware_dtype(dtype: DtypeObj) -> bool:
        if isinstance(dtype, DatetimeTZDtype):
            return True
        if isinstance(dtype, ArrowDtype):
            pa = import_optional_dependency("pyarrow")
            pa_dtype = dtype.pyarrow_dtype
            return pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None
        return False

    def _pyarrow_renders_differently(self) -> bool:
        """
        Whether any column/index level being written has a dtype that
        the pyarrow engine renders in a materially different (though
        still value-preserving) textual form than the python engine:
        bool -> "true"/"false" instead of "True"/"False", timedelta64 ->
        raw integer nanoseconds instead of a human-readable duration, and
        timezone-aware datetime64 -> a different offset/precision format.
        """
        for dtype in self._iter_written_dtypes():
            if dtype.kind in ("b", "m") or self._is_tz_aware_dtype(dtype):
                return True
        return False

    @staticmethod
    def _has_whole_number_float_column(obj: DataFrame) -> bool:
        """
        Whether any float-dtype column in `obj` (which already has any
        index moved into a column by _build_arrow_obj) holds only
        whole-number values (ignoring missing values). The python engine
        always writes a trailing ".0" for such values (e.g. "1.0"), which
        is what lets read_csv infer a float dtype back on round-trip; the
        pyarrow engine omits it (writing "1"), which can silently change
        the dtype read_csv infers to an integer type instead.
        """
        for _, col in obj.items():
            if col.dtype.kind != "f":
                continue
            non_null = col[col.notna()]
            if len(non_null) and (non_null % 1 == 0).all():
                return True
        return False

    def _build_arrow_obj(self) -> DataFrame:
        """
        Return self.obj with the index moved into columns (mirroring how
        the python engine writes the index), ready for pa.Table.from_pandas.
        """
        obj = self.obj
        if self._has_aliases:
            # header=[...] aliases only rename the CSV output, not obj
            # itself, but pyarrow's writer has no separate notion of a
            # header distinct from the table's own schema, so give it
            # column names that already match the requested aliases.
            obj = obj.rename(
                columns=dict(zip(self.obj.columns, self.write_cols, strict=True))
            )

        if not self.index:
            return obj
        # Use unique placeholder names to reset_index with, since the
        # index levels may be unnamed (and thus collide with each other,
        # or with "") once converted to their final header labels.
        nlevels = self.obj.index.nlevels
        placeholders = [f"__pandas_index_level_{i}__" for i in range(nlevels)]
        obj = obj.reset_index(names=placeholders)
        target_names = [
            name if name is not None else "" for name in self.obj.index.names
        ]
        return obj.rename(columns=dict(zip(placeholders, target_names, strict=True)))

    def _resolve_engine(self, requested_engine: str) -> str:
        """
        Resolve the user-requested engine ("python", "pyarrow", or "auto")
        to the concrete engine ("python" or "pyarrow") that will actually
        be used, validating/warning/falling back as appropriate.
        """
        if requested_engine not in ("python", "pyarrow", "auto"):
            raise ValueError(
                "engine must be one of 'python', 'pyarrow', or 'auto', "
                f"got {requested_engine!r}"
            )
        if requested_engine == "python":
            return "python"

        explicit = requested_engine == "pyarrow"

        option_incompatibility = self._pyarrow_option_incompatibility()
        if option_incompatibility is not None:
            if explicit:
                raise ValueError(option_incompatibility)
            return "python"

        pa = import_optional_dependency(
            "pyarrow", errors="raise" if explicit else "ignore"
        )
        if pa is None:
            return "python"

        arrow_obj = self._build_arrow_obj()

        if self._pyarrow_renders_differently() or self._has_whole_number_float_column(
            arrow_obj
        ):
            if explicit:
                warnings.warn(
                    "The pyarrow engine renders bool, timedelta64, and "
                    "timezone-aware datetime64 columns differently than "
                    "the python engine (e.g. 'true'/'false' instead of "
                    "'True'/'False'), and omits the trailing '.0' on a "
                    "whole-number float column (e.g. '1' instead of "
                    "'1.0'), which can change the dtype read_csv infers "
                    "on round-trip.",
                    UserWarning,
                    stacklevel=find_stack_level(),
                )
            else:
                return "python"

        if arrow_obj.isna().all(axis=1).any():
            reason = (
                "The pyarrow engine cannot represent a row that is "
                "entirely missing values without writing a blank output "
                "line, which would be silently dropped when the file is "
                "read back with the default skip_blank_lines=True."
            )
            if explicit:
                raise ValueError(reason)
            return "python"

        try:
            # preserve_index=False: the index (if any) is already moved
            # into a column by _build_arrow_obj above; letting pyarrow's
            # own preserve_index heuristic run too would either duplicate
            # it or -- when self.index is False -- add it back as an
            # unwanted extra column.
            table = pa.Table.from_pandas(arrow_obj, preserve_index=False)
        except (pa.lib.ArrowException, TypeError, ValueError):
            if explicit:
                raise
            return "python"

        # Some dtypes (e.g. Period, Interval extension types, or Python
        # list/dict objects that pyarrow infers as a list/struct type)
        # convert to a pyarrow Table successfully, but pyarrow's CSV writer
        # cannot serialize them and only raises once actually writing.
        # Detect that up front with a cheap trial write of a small slice
        # (never more than a handful of rows), rather than letting it raise
        # mid-write into the real destination -- and rather than trying to
        # maintain a list of every pyarrow type the CSV writer rejects.
        pa_csv = import_optional_dependency("pyarrow.csv")

        # NotImplementedError signals an unsupported *quoting* value (e.g.
        # QUOTE_NONNUMERIC) -- an engine capability gap, handled like any
        # other pyarrow-incompatibility below. A TypeError here (quotechar
        # not set while quoting is enabled) is a genuine usage error that
        # the python engine would reject too, so it is left to propagate
        # as-is rather than being treated as an auto-fallback trigger.
        try:
            write_options = self._build_pyarrow_write_options(pa_csv)
        except NotImplementedError:
            reason = f"The pyarrow engine does not support quoting={self.quoting!r}."
            if explicit:
                raise ValueError(reason) from None
            return "python"

        try:
            pa_csv.write_csv(table.slice(0, 1), io.BytesIO(), write_options)
        except (pa.lib.ArrowException, TypeError, ValueError, NotImplementedError):
            reason = (
                "The pyarrow engine cannot write one or more columns in "
                "this data (e.g. Period/Interval dtype, or an object "
                "column of list/dict values)."
            )
            if explicit:
                raise ValueError(reason) from None
            return "python"

        self._arrow_table = table
        return "pyarrow"

    def save(self) -> None:
        """
        Create the writer & save.
        """
        if self.mode is None:
            self.mode = "w" + ("b" if self.engine == "pyarrow" else "")

        if self.filepath_or_buffer is None:
            self.filepath_or_buffer = (
                io.BytesIO() if self.engine == "pyarrow" else io.StringIO()
            )

        if self.engine == "pyarrow" and (
            "b" not in self.mode
            or isinstance(self.filepath_or_buffer, _TEXT_LIKE_CLASSES)
        ):
            raise ValueError("The pyarrow engine can only open files in binary mode.")

        # apply compression and byte/text conversion
        with get_handle(
            self.filepath_or_buffer,
            self.mode,
            encoding=self.encoding,
            errors=self.errors,
            compression=self.compression,
            storage_options=self.storage_options,
            # pyarrow engine exclusively writes bytes
            is_text=self.engine == "python",
        ) as handles:
            # Note: self.encoding is irrelevant here

            # This is a mypy bug?
            # error: Cannot infer type argument 1 of "_save" of "CSVFormatter"  [misc]
            self._save(handles.handle)  # type: ignore[misc]

    def _save_pyarrow(self, handle: IO[AnyStr]) -> None:
        pa_csv = import_optional_dependency("pyarrow.csv")

        # built and validated during engine resolution in __init__
        table = self._arrow_table
        write_options = self._build_pyarrow_write_options(pa_csv)
        pa_csv.write_csv(table, handle, write_options)

    def _build_pyarrow_write_options(self, pa_csv: Any) -> Any:
        # Map quoting arg to pyarrow equivalents. QUOTE_NONE is checked
        # before the quotechar-is-None case below, since quotechar is
        # legitimately None whenever quoting=QUOTE_NONE and no escapechar
        # was given (see _initialize_quotechar).
        if self.quoting == csvlib.QUOTE_MINIMAL:
            pa_quoting = "needed"
        elif self.quoting == csvlib.QUOTE_NONE:
            pa_quoting = "none"
        elif self.quotechar is None:
            raise TypeError("quotechar must be set if quoting enabled")
        elif self.quoting == csvlib.QUOTE_ALL:
            # TODO: Is this a 1-1 mapping?
            # This doesn't quote nulls, check if Python does this
            pa_quoting = "all_valid"
        else:
            raise NotImplementedError(
                f"Quoting option {self.quoting} is not supported with engine='pyarrow'"
            )

        kwargs: dict[str, Any] = {
            "include_header": self._need_to_save_header,
            "batch_size": self.chunksize,
        }
        kwargs["delimiter"] = self.sep
        kwargs["quoting_style"] = pa_quoting

        return pa_csv.WriteOptions(**kwargs)

    def _save(self, handle: IO[AnyStr]) -> None:
        if self.engine == "pyarrow":
            self._save_pyarrow(handle)
        else:
            self.writer = csvlib.writer(
                # error: Argument of type "IO[AnyStr@_save]" cannot be assigned
                # to parameter "csvfile" of type "SupportsWrite[str]"
                # in function "writer"
                # error: Argument "quoting" to "writer" has incompatible type "int";
                # expected "Literal[0, 1, 2, 3]"
                handle,  # type: ignore[arg-type]
                lineterminator=self.lineterminator,
                delimiter=self.sep,
                quoting=self.quoting,  # type: ignore[arg-type]
                doublequote=self.doublequote,
                escapechar=self.escapechar,
                quotechar=self.quotechar,
            )
            if self._need_to_save_header:
                self._save_header()
            self._save_body()

    def _save_header(self) -> None:
        if not self.has_mi_columns or self._has_aliases:
            self.writer.writerow(self.encoded_labels)
        else:
            for row in self._generate_multiindex_header_rows():
                self.writer.writerow(row)

    def _generate_multiindex_header_rows(self) -> Iterator[list[Hashable]]:
        columns = self.obj.columns
        for i in range(columns.nlevels):
            # we need at least 1 index column to write our col names
            col_line = []
            if self.index:
                # name is the first column
                col_line.append(columns.names[i])

                if isinstance(self.index_label, list) and len(self.index_label) > 1:
                    col_line.extend([""] * (len(self.index_label) - 1))

            col_line.extend(columns._get_level_values(i))
            yield col_line

        # Write out the index line if it's not empty.
        # Otherwise, we will print out an extraneous
        # blank line between the mi and the data rows.
        if self.encoded_labels and set(self.encoded_labels) != {""}:
            yield self.encoded_labels + [""] * len(columns)

    def _save_body(self) -> None:
        nrows = len(self.data_index)
        chunks = (nrows // self.chunksize) + 1
        for i in range(chunks):
            start_i = i * self.chunksize
            end_i = min(start_i + self.chunksize, nrows)
            if start_i >= end_i:
                break
            self._save_chunk(start_i, end_i)

    def _save_chunk(self, start_i: int, end_i: int) -> None:
        # create the data for a chunk
        slicer = slice(start_i, end_i)
        df = self.obj.iloc[slicer]

        res = df._get_values_for_csv(**self._number_format)
        data = list(res._iter_column_arrays())

        ix = (
            self.data_index[slicer]._get_values_for_csv(**self._number_format)
            if self.nlevels != 0
            else np.empty(end_i - start_i)
        )
        libwriters.write_csv_rows(
            data,
            ix,
            self.nlevels,
            self.cols,
            self.writer,
        )
