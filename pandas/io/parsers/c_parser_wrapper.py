import pandas._libs.parsers as parsers
from pandas._typing import FilePathOrBuffer

from pandas.core.indexes.api import ensure_index_from_sequences

from pandas.io.parsers.base_parser import (
    ParserBase,
    is_index_col,
)


class CParserWrapper(ParserBase):
    def __init__(self, src: FilePathOrBuffer, **kwds):
        self.kwds = kwds
        kwds = kwds.copy()

        ParserBase.__init__(self, kwds)

        # #2442
        kwds["allow_leading_cols"] = self.index_col is not False

        # GH20529, validate usecol arg before TextReader
        kwds["usecols"] = self.usecols

        # open handles
        self._open_handles(src, kwds)
        assert self.handles is not None
        for key in ("storage_options", "encoding", "memory_map", "compression"):
            kwds.pop(key, None)
        if self.handles.is_mmap and hasattr(self.handles.handle, "mmap"):
            # error: Item "IO[Any]" of "Union[IO[Any], RawIOBase, BufferedIOBase,
            # TextIOBase, TextIOWrapper, mmap]" has no attribute "mmap"

            # error: Item "RawIOBase" of "Union[IO[Any], RawIOBase, BufferedIOBase,
            # TextIOBase, TextIOWrapper, mmap]" has no attribute "mmap"

            # error: Item "BufferedIOBase" of "Union[IO[Any], RawIOBase, BufferedIOBase,
            # TextIOBase, TextIOWrapper, mmap]" has no attribute "mmap"

            # error: Item "TextIOBase" of "Union[IO[Any], RawIOBase, BufferedIOBase,
            # TextIOBase, TextIOWrapper, mmap]" has no attribute "mmap"

            # error: Item "TextIOWrapper" of "Union[IO[Any], RawIOBase, BufferedIOBase,
            # TextIOBase, TextIOWrapper, mmap]" has no attribute "mmap"

            # error: Item "mmap" of "Union[IO[Any], RawIOBase, BufferedIOBase,
            # TextIOBase, TextIOWrapper, mmap]" has no attribute "mmap"
            self.handles.handle = self.handles.handle.mmap  # type: ignore[union-attr]

        try:
            self._reader = parsers.TextReader(self.handles.handle, **kwds)
        except Exception:
            self.handles.close()
            raise
        self.unnamed_cols = self._reader.unnamed_cols

        passed_names = self.names is None

        if self._reader.header is None:
            self.names = None
        else:
            if len(self._reader.header) > 1:
                # we have a multi index in the columns
                (
                    self.names,
                    self.index_names,
                    self.col_names,
                    passed_names,
                ) = self._extract_multi_indexer_columns(
                    self._reader.header, self.index_names, self.col_names, passed_names
                )
            else:
                self.names = list(self._reader.header[0])

        if self.names is None:
            if self.prefix:
                self.names = [
                    f"{self.prefix}{i}" for i in range(self._reader.table_width)
                ]
            else:
                self.names = list(range(self._reader.table_width))

        # gh-9755
        #
        # need to set orig_names here first
        # so that proper indexing can be done
        # with _set_noconvert_columns
        #
        # once names has been filtered, we will
        # then set orig_names again to names
        self.orig_names = self.names[:]

        if self.usecols:
            usecols = self._evaluate_usecols(self.usecols, self.orig_names)

            # GH 14671
            # assert for mypy, orig_names is List or None, None would error in issubset
            assert self.orig_names is not None
            if self.usecols_dtype == "string" and not set(usecols).issubset(
                self.orig_names
            ):
                self._validate_usecols_names(usecols, self.orig_names)

            if len(self.names) > len(usecols):
                self.names = [
                    n
                    for i, n in enumerate(self.names)
                    if (i in usecols or n in usecols)
                ]

            if len(self.names) < len(usecols):
                self._validate_usecols_names(usecols, self.names)

        self._validate_parse_dates_presence(self.names)
        self._set_noconvert_columns()

        self.orig_names = self.names

        if not self._has_complex_date_col:
            if self._reader.leading_cols == 0 and is_index_col(self.index_col):

                self._name_processed = True
                (index_names, self.names, self.index_col) = self._clean_index_names(
                    self.names, self.index_col, self.unnamed_cols
                )

                if self.index_names is None:
                    self.index_names = index_names

            if self._reader.header is None and not passed_names:
                assert self.index_names is not None
                self.index_names = [None] * len(self.index_names)

        self._implicit_index = self._reader.leading_cols > 0

    def close(self) -> None:
        super().close()

        # close additional handles opened by C parser
        try:
            self._reader.close()
        except ValueError:
            pass

    def _set_noconvert_columns(self):
        """
        Set the columns that should not undergo dtype conversions.

        Currently, any column that is involved with date parsing will not
        undergo such conversions.
        """
        assert self.orig_names is not None
        col_indices = [self.orig_names.index(x) for x in self.names]
        noconvert_columns = self._set_noconvert_dtype_columns(col_indices, self.names)
        for col in noconvert_columns:
            self._reader.set_noconvert(col)

    def set_error_bad_lines(self, status):
        self._reader.set_error_bad_lines(int(status))

    def read(self, nrows=None):
        try:
            data = self._reader.read(nrows)
        except StopIteration:
            if self._first_chunk:
                self._first_chunk = False
                names = self._maybe_dedup_names(self.orig_names)
                index, columns, col_dict = self._get_empty_meta(
                    names,
                    self.index_col,
                    self.index_names,
                    dtype=self.kwds.get("dtype"),
                )
                columns = self._maybe_make_multi_index_columns(columns, self.col_names)

                if self.usecols is not None:
                    columns = self._filter_usecols(columns)

                col_dict = {k: v for k, v in col_dict.items() if k in columns}

                return index, columns, col_dict

            else:
                self.close()
                raise

        # Done with first read, next time raise StopIteration
        self._first_chunk = False

        names = self.names

        if self._reader.leading_cols:
            if self._has_complex_date_col:
                raise NotImplementedError("file structure not yet supported")

            # implicit index, no index names
            arrays = []

            for i in range(self._reader.leading_cols):
                if self.index_col is None:
                    values = data.pop(i)
                else:
                    values = data.pop(self.index_col[i])

                values = self._maybe_parse_dates(values, i, try_parse_dates=True)
                arrays.append(values)

            index = ensure_index_from_sequences(arrays)

            if self.usecols is not None:
                names = self._filter_usecols(names)

            names = self._maybe_dedup_names(names)

            # rename dict keys
            data = sorted(data.items())
            data = {k: v for k, (i, v) in zip(names, data)}

            names, data = self._do_date_conversions(names, data)

        else:
            # rename dict keys
            data = sorted(data.items())

            # ugh, mutation

            # assert for mypy, orig_names is List or None, None would error in list(...)
            assert self.orig_names is not None
            names = list(self.orig_names)
            names = self._maybe_dedup_names(names)

            if self.usecols is not None:
                names = self._filter_usecols(names)

            # columns as list
            alldata = [x[1] for x in data]

            data = {k: v for k, (i, v) in zip(names, data)}

            names, data = self._do_date_conversions(names, data)
            index, names = self._make_index(data, alldata, names)

        # maybe create a mi on the columns
        names = self._maybe_make_multi_index_columns(names, self.col_names)

        return index, names, data

    def _filter_usecols(self, names):
        # hackish
        usecols = self._evaluate_usecols(self.usecols, names)
        if usecols is not None and len(names) != len(usecols):
            names = [
                name for i, name in enumerate(names) if i in usecols or name in usecols
            ]
        return names

    def _get_index_names(self):
        names = list(self._reader.header[0])
        idx_names = None

        if self._reader.leading_cols == 0 and self.index_col is not None:
            (idx_names, names, self.index_col) = self._clean_index_names(
                names, self.index_col, self.unnamed_cols
            )

        return names, idx_names

    def _maybe_parse_dates(self, values, index, try_parse_dates=True):
        if try_parse_dates and self._should_parse_dates(index):
            values = self._date_conv(values)
        return values
