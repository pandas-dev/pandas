# -*- coding: utf-8 -*-
"""
Module for formatting output data into CSV files.
"""

from __future__ import print_function

import warnings

import csv as csvlib
from zipfile import ZipFile

import numpy as np

from pandas._libs import writers as libwriters

from pandas import compat
from pandas.compat import StringIO, range, zip

from pandas.core.dtypes.missing import notna
from pandas.core.dtypes.generic import (
    ABCMultiIndex, ABCPeriodIndex, ABCDatetimeIndex, ABCIndexClass)

from pandas.io.common import (
    _expand_user,
    _get_handle,
    _infer_compression,
    _stringify_path,
    UnicodeWriter,
)


class CSVFormatter(object):

    def __init__(self, obj, path_or_buf=None, sep=",", na_rep='',
                 float_format=None, cols=None, header=True, index=True,
                 index_label=None, mode='w', nanRep=None, encoding=None,
                 compression='infer', quoting=None, line_terminator='\n',
                 chunksize=None, tupleize_cols=False, quotechar='"',
                 date_format=None, doublequote=True, escapechar=None,
                 decimal='.'):

        self.obj = obj

        if path_or_buf is None:
            path_or_buf = StringIO()

        self.path_or_buf = _expand_user(_stringify_path(path_or_buf))
        self.sep = sep
        self.na_rep = na_rep
        self.float_format = float_format
        self.decimal = decimal

        self.header = header
        self.index = index
        self.index_label = index_label
        self.mode = mode
        if encoding is None:
            encoding = 'ascii' if compat.PY2 else 'utf-8'
        self.encoding = encoding
        self.compression = _infer_compression(self.path_or_buf, compression)

        if quoting is None:
            quoting = csvlib.QUOTE_MINIMAL
        self.quoting = quoting

        if quoting == csvlib.QUOTE_NONE:
            # prevents crash in _csv
            quotechar = None
        self.quotechar = quotechar

        self.doublequote = doublequote
        self.escapechar = escapechar

        self.line_terminator = line_terminator

        self.date_format = date_format

        self.tupleize_cols = tupleize_cols
        self.has_mi_columns = (isinstance(obj.columns, ABCMultiIndex) and
                               not self.tupleize_cols)

        # validate mi options
        if self.has_mi_columns:
            if cols is not None:
                raise TypeError("cannot specify cols with a MultiIndex on the "
                                "columns")

        if cols is not None:
            if isinstance(cols, ABCIndexClass):
                cols = cols.to_native_types(na_rep=na_rep,
                                            float_format=float_format,
                                            date_format=date_format,
                                            quoting=self.quoting)
            else:
                cols = list(cols)
            self.obj = self.obj.loc[:, cols]

        # update columns to include possible multiplicity of dupes
        # and make sure sure cols is just a list of labels
        cols = self.obj.columns
        if isinstance(cols, ABCIndexClass):
            cols = cols.to_native_types(na_rep=na_rep,
                                        float_format=float_format,
                                        date_format=date_format,
                                        quoting=self.quoting)
        else:
            cols = list(cols)

        # save it
        self.cols = cols

        # preallocate data 2d list
        self.blocks = self.obj._data.blocks
        ncols = sum(b.shape[0] for b in self.blocks)
        self.data = [None] * ncols

        if chunksize is None:
            chunksize = (100000 // (len(self.cols) or 1)) or 1
        self.chunksize = int(chunksize)

        self.data_index = obj.index
        if (isinstance(self.data_index, (ABCDatetimeIndex, ABCPeriodIndex)) and
                date_format is not None):
            from pandas import Index
            self.data_index = Index([x.strftime(date_format) if notna(x) else
                                     '' for x in self.data_index])

        self.nlevels = getattr(self.data_index, 'nlevels', 1)
        if not index:
            self.nlevels = 0

    def save(self):
        """
        Create the writer & save
        """
        # GH21227 internal compression is not used when file-like passed.
        if self.compression and hasattr(self.path_or_buf, 'write'):
            msg = ("compression has no effect when passing file-like "
                   "object as input.")
            warnings.warn(msg, RuntimeWarning, stacklevel=2)

        # when zip compression is called.
        is_zip = isinstance(self.path_or_buf, ZipFile) or (
            not hasattr(self.path_or_buf, 'write')
            and self.compression == 'zip')

        if is_zip:
            # zipfile doesn't support writing string to archive. uses string
            # buffer to receive csv writing and dump into zip compression
            # file handle. GH21241, GH21118
            f = StringIO()
            close = False
        elif hasattr(self.path_or_buf, 'write'):
            f = self.path_or_buf
            close = False
        else:
            f, handles = _get_handle(self.path_or_buf, self.mode,
                                     encoding=self.encoding,
                                     compression=self.compression)
            close = True

        try:
            writer_kwargs = dict(lineterminator=self.line_terminator,
                                 delimiter=self.sep, quoting=self.quoting,
                                 doublequote=self.doublequote,
                                 escapechar=self.escapechar,
                                 quotechar=self.quotechar)
            if self.encoding == 'ascii':
                self.writer = csvlib.writer(f, **writer_kwargs)
            else:
                writer_kwargs['encoding'] = self.encoding
                self.writer = UnicodeWriter(f, **writer_kwargs)

            self._save()

        finally:
            if is_zip:
                # GH17778 handles zip compression separately.
                buf = f.getvalue()
                if hasattr(self.path_or_buf, 'write'):
                    self.path_or_buf.write(buf)
                else:
                    f, handles = _get_handle(self.path_or_buf, self.mode,
                                             encoding=self.encoding,
                                             compression=self.compression)
                    f.write(buf)
                    close = True
            if close:
                f.close()
                for _fh in handles:
                    _fh.close()

    def _save_header(self):

        writer = self.writer
        obj = self.obj
        index_label = self.index_label
        cols = self.cols
        has_mi_columns = self.has_mi_columns
        header = self.header
        encoded_labels = []

        has_aliases = isinstance(header, (tuple, list, np.ndarray,
                                          ABCIndexClass))
        if not (has_aliases or self.header):
            return
        if has_aliases:
            if len(header) != len(cols):
                raise ValueError(('Writing {ncols} cols but got {nalias} '
                                 'aliases'.format(ncols=len(cols),
                                                  nalias=len(header))))
            else:
                write_cols = header
        else:
            write_cols = cols

        if self.index:
            # should write something for index label
            if index_label is not False:
                if index_label is None:
                    if isinstance(obj.index, ABCMultiIndex):
                        index_label = []
                        for i, name in enumerate(obj.index.names):
                            if name is None:
                                name = ''
                            index_label.append(name)
                    else:
                        index_label = obj.index.name
                        if index_label is None:
                            index_label = ['']
                        else:
                            index_label = [index_label]
                elif not isinstance(index_label,
                                    (list, tuple, np.ndarray, ABCIndexClass)):
                    # given a string for a DF with Index
                    index_label = [index_label]

                encoded_labels = list(index_label)
            else:
                encoded_labels = []

        if not has_mi_columns or has_aliases:
            encoded_labels += list(write_cols)
            writer.writerow(encoded_labels)
        else:
            # write out the mi
            columns = obj.columns

            # write out the names for each level, then ALL of the values for
            # each level
            for i in range(columns.nlevels):

                # we need at least 1 index column to write our col names
                col_line = []
                if self.index:

                    # name is the first column
                    col_line.append(columns.names[i])

                    if isinstance(index_label, list) and len(index_label) > 1:
                        col_line.extend([''] * (len(index_label) - 1))

                col_line.extend(columns._get_level_values(i))

                writer.writerow(col_line)

            # Write out the index line if it's not empty.
            # Otherwise, we will print out an extraneous
            # blank line between the mi and the data rows.
            if encoded_labels and set(encoded_labels) != set(['']):
                encoded_labels.extend([''] * len(columns))
                writer.writerow(encoded_labels)

    def _save(self):

        self._save_header()

        nrows = len(self.data_index)

        # write in chunksize bites
        chunksize = self.chunksize
        chunks = int(nrows / chunksize) + 1

        for i in range(chunks):
            start_i = i * chunksize
            end_i = min((i + 1) * chunksize, nrows)
            if start_i >= end_i:
                break

            self._save_chunk(start_i, end_i)

    def _save_chunk(self, start_i, end_i):

        data_index = self.data_index

        # create the data for a chunk
        slicer = slice(start_i, end_i)
        for i in range(len(self.blocks)):
            b = self.blocks[i]
            d = b.to_native_types(slicer=slicer, na_rep=self.na_rep,
                                  float_format=self.float_format,
                                  decimal=self.decimal,
                                  date_format=self.date_format,
                                  quoting=self.quoting)

            for col_loc, col in zip(b.mgr_locs, d):
                # self.data is a preallocated list
                self.data[col_loc] = col

        ix = data_index.to_native_types(slicer=slicer, na_rep=self.na_rep,
                                        float_format=self.float_format,
                                        decimal=self.decimal,
                                        date_format=self.date_format,
                                        quoting=self.quoting)

        libwriters.write_csv_rows(self.data, ix, self.nlevels,
                                  self.cols, self.writer)
