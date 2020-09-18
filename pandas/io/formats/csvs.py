"""
Module for formatting output data into CSV files.
"""

import csv as csvlib
from io import StringIO, TextIOWrapper
import os
from typing import Any, Dict, Hashable, Iterator, List, Optional, Sequence, Union

import numpy as np

from pandas._libs import writers as libwriters
from pandas._typing import (
    CompressionOptions,
    FilePathOrBuffer,
    IndexLabel,
    Label,
    StorageOptions,
)

from pandas.core.dtypes.generic import (
    ABCDatetimeIndex,
    ABCIndexClass,
    ABCMultiIndex,
    ABCPeriodIndex,
)
from pandas.core.dtypes.missing import notna

from pandas.core.indexes.api import Index

from pandas.io.common import get_filepath_or_buffer, get_handle


class CSVFormatter:
    def __init__(
        self,
        obj,
        path_or_buf: Optional[FilePathOrBuffer[str]] = None,
        sep: str = ",",
        na_rep: str = "",
        float_format: Optional[str] = None,
        cols: Optional[Sequence[Label]] = None,
        header: Union[bool, Sequence[Hashable]] = True,
        index: bool = True,
        index_label: Optional[IndexLabel] = None,
        mode: str = "w",
        encoding: Optional[str] = None,
        errors: str = "strict",
        compression: CompressionOptions = "infer",
        quoting: Optional[int] = None,
        line_terminator="\n",
        chunksize: Optional[int] = None,
        quotechar: Optional[str] = '"',
        date_format: Optional[str] = None,
        doublequote: bool = True,
        escapechar: Optional[str] = None,
        decimal=".",
        storage_options: StorageOptions = None,
    ):
        self.obj = obj

        self.encoding = encoding or "utf-8"

        if path_or_buf is None:
            path_or_buf = StringIO()

        ioargs = get_filepath_or_buffer(
            path_or_buf,
            encoding=self.encoding,
            compression=compression,
            mode=mode,
            storage_options=storage_options,
        )

        self.compression = ioargs.compression.pop("method")
        self.compression_args = ioargs.compression
        self.path_or_buf = ioargs.filepath_or_buffer
        self.should_close = ioargs.should_close
        self.mode = ioargs.mode

        self.sep = sep
        self.na_rep = na_rep
        self.float_format = float_format
        self.decimal = decimal
        self.header = header
        self.index = index
        self.index_label = index_label
        self.errors = errors
        self.quoting = quoting or csvlib.QUOTE_MINIMAL
        self.quotechar = quotechar
        self.doublequote = doublequote
        self.escapechar = escapechar
        self.line_terminator = line_terminator or os.linesep
        self.date_format = date_format
        self.cols = cols  # type: ignore[assignment]
        self.chunksize = chunksize  # type: ignore[assignment]

    @property
    def index_label(self) -> IndexLabel:
        return self._index_label

    @index_label.setter
    def index_label(self, index_label: Optional[IndexLabel]) -> None:
        if index_label is not False:
            if index_label is None:
                index_label = self._get_index_label_from_obj()
            elif not isinstance(index_label, (list, tuple, np.ndarray, ABCIndexClass)):
                # given a string for a DF with Index
                index_label = [index_label]
        self._index_label = index_label

    def _get_index_label_from_obj(self) -> List[str]:
        if isinstance(self.obj.index, ABCMultiIndex):
            return self._get_index_label_multiindex()
        else:
            return self._get_index_label_flat()

    def _get_index_label_multiindex(self) -> List[str]:
        return [name or "" for name in self.obj.index.names]

    def _get_index_label_flat(self) -> List[str]:
        index_label = self.obj.index.name
        return [""] if index_label is None else [index_label]

    @property
    def quotechar(self) -> Optional[str]:
        if self.quoting != csvlib.QUOTE_NONE:
            # prevents crash in _csv
            return self._quotechar
        return None

    @quotechar.setter
    def quotechar(self, quotechar: Optional[str]) -> None:
        self._quotechar = quotechar

    @property
    def has_mi_columns(self) -> bool:
        return bool(isinstance(self.obj.columns, ABCMultiIndex))

    @property
    def cols(self) -> Sequence[Label]:
        return self._cols

    @cols.setter
    def cols(self, cols: Optional[Sequence[Label]]) -> None:
        self._cols = self._refine_cols(cols)

    def _refine_cols(self, cols: Optional[Sequence[Label]]) -> Sequence[Label]:
        # validate mi options
        if self.has_mi_columns:
            if cols is not None:
                msg = "cannot specify cols with a MultiIndex on the columns"
                raise TypeError(msg)

        if cols is not None:
            if isinstance(cols, ABCIndexClass):
                cols = cols.to_native_types(**self._number_format)
            else:
                cols = list(cols)
            self.obj = self.obj.loc[:, cols]

        # update columns to include possible multiplicity of dupes
        # and make sure sure cols is just a list of labels
        cols = self.obj.columns
        if isinstance(cols, ABCIndexClass):
            return cols.to_native_types(**self._number_format)
        else:
            assert isinstance(cols, Sequence)
            return list(cols)

    @property
    def _number_format(self) -> Dict[str, Any]:
        """Dictionary used for storing number formatting settings."""
        return dict(
            na_rep=self.na_rep,
            float_format=self.float_format,
            date_format=self.date_format,
            quoting=self.quoting,
            decimal=self.decimal,
        )

    @property
    def chunksize(self) -> int:
        return self._chunksize

    @chunksize.setter
    def chunksize(self, chunksize: Optional[int]) -> None:
        if chunksize is None:
            chunksize = (100000 // (len(self.cols) or 1)) or 1
        assert chunksize is not None
        self._chunksize = int(chunksize)

    @property
    def data_index(self) -> Index:
        data_index = self.obj.index
        if (
            isinstance(data_index, (ABCDatetimeIndex, ABCPeriodIndex))
            and self.date_format is not None
        ):
            data_index = Index(
                [x.strftime(self.date_format) if notna(x) else "" for x in data_index]
            )
        return data_index

    @property
    def nlevels(self) -> int:
        if self.index:
            return getattr(self.data_index, "nlevels", 1)
        else:
            return 0

    @property
    def _has_aliases(self) -> bool:
        return isinstance(self.header, (tuple, list, np.ndarray, ABCIndexClass))

    @property
    def _need_to_save_header(self) -> bool:
        return bool(self._has_aliases or self.header)

    @property
    def write_cols(self) -> Sequence[Label]:
        if self._has_aliases:
            assert not isinstance(self.header, bool)
            if len(self.header) != len(self.cols):
                raise ValueError(
                    f"Writing {len(self.cols)} cols but got {len(self.header)} aliases"
                )
            else:
                return self.header
        else:
            return self.cols

    @property
    def encoded_labels(self) -> List[Label]:
        encoded_labels: List[Label] = []

        if self.index and self.index_label:
            assert isinstance(self.index_label, Sequence)
            encoded_labels = list(self.index_label)

        if not self.has_mi_columns or self._has_aliases:
            encoded_labels += list(self.write_cols)

        return encoded_labels

    def save(self) -> None:
        """
        Create the writer & save.
        """
        # get a handle or wrap an existing handle to take care of 1) compression and
        # 2) text -> byte conversion
        f, handles = get_handle(
            self.path_or_buf,
            self.mode,
            encoding=self.encoding,
            errors=self.errors,
            compression=dict(self.compression_args, method=self.compression),
        )

        try:
            # Note: self.encoding is irrelevant here
            self.writer = csvlib.writer(
                f,
                lineterminator=self.line_terminator,
                delimiter=self.sep,
                quoting=self.quoting,
                doublequote=self.doublequote,
                escapechar=self.escapechar,
                quotechar=self.quotechar,
            )

            self._save()

        finally:
            if self.should_close:
                f.close()
            elif (
                isinstance(f, TextIOWrapper)
                and not f.closed
                and f != self.path_or_buf
                and hasattr(self.path_or_buf, "write")
            ):
                # get_handle uses TextIOWrapper for non-binary handles. TextIOWrapper
                # closes the wrapped handle if it is not detached.
                f.flush()  # make sure everything is written
                f.detach()  # makes f unusable
                del f
            elif f != self.path_or_buf:
                f.close()
            for _fh in handles:
                _fh.close()

    def _save(self) -> None:
        if self._need_to_save_header:
            self._save_header()
        self._save_body()

    def _save_header(self) -> None:
        if not self.has_mi_columns or self._has_aliases:
            self.writer.writerow(self.encoded_labels)
        else:
            for row in self._generate_multiindex_header_rows():
                self.writer.writerow(row)

    def _generate_multiindex_header_rows(self) -> Iterator[List[Label]]:
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
        chunks = int(nrows / self.chunksize) + 1
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

        res = df._mgr.to_native_types(**self._number_format)
        data = [res.iget_values(i) for i in range(len(res.items))]

        ix = self.data_index.to_native_types(slicer=slicer, **self._number_format)
        libwriters.write_csv_rows(data, ix, self.nlevels, self.cols, self.writer)
