# -*- coding: utf-8 -*-
from __future__ import print_function
from distutils.version import LooseVersion
# pylint: disable=W0141

import sys

from pandas.core.base import PandasObject
from pandas.core.common import adjoin, notnull
from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas import compat
from pandas.compat import(StringIO, lzip, range, map, zip, reduce, u,
                          OrderedDict)
from pandas.util.terminal import get_terminal_size
from pandas.core.config import get_option, set_option
import pandas.core.common as com
import pandas.lib as lib
from pandas.tslib import iNaT, Timestamp, Timedelta, format_array_from_datetime
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex
import pandas as pd
import numpy as np

import itertools
import csv

common_docstring = """
    Parameters
    ----------
    buf : StringIO-like, optional
        buffer to write to
    columns : sequence, optional
        the subset of columns to write; default None writes all columns
    col_space : int, optional
        the minimum width of each column
    header : bool, optional
        whether to print column labels, default True
    index : bool, optional
        whether to print index (row) labels, default True
    na_rep : string, optional
        string representation of NAN to use, default 'NaN'
    formatters : list or dict of one-parameter functions, optional
        formatter functions to apply to columns' elements by position or name,
        default None. The result of each function must be a unicode string.
        List must be of length equal to the number of columns.
    float_format : one-parameter function, optional
        formatter function to apply to columns' elements if they are floats,
        default None. The result of this function must be a unicode string.
    sparsify : bool, optional
        Set to False for a DataFrame with a hierarchical index to print every
        multiindex key at each row, default True
    index_names : bool, optional
        Prints the names of the indexes, default True"""

justify_docstring  = """
    justify : {'left', 'right'}, default None
        Left or right-justify the column labels. If None uses the option from
        the print configuration (controlled by set_option), 'right' out
        of the box."""

return_docstring = """

    Returns
    -------
    formatted : string (or unicode, depending on data and options)"""

docstring_to_string = common_docstring + justify_docstring + return_docstring

class CategoricalFormatter(object):

    def __init__(self, categorical, buf=None, length=True,
                 na_rep='NaN', footer=True):
        self.categorical = categorical
        self.buf = buf if buf is not None else StringIO(u(""))
        self.na_rep = na_rep
        self.length = length
        self.footer = footer

    def _get_footer(self):
        footer = ''

        if self.length:
            if footer:
                footer += ', '
            footer += "Length: %d" % len(self.categorical)

        level_info = self.categorical._repr_categories_info()

        # Levels are added in a newline
        if footer:
            footer += '\n'
        footer += level_info

        return compat.text_type(footer)

    def _get_formatted_values(self):
        return format_array(self.categorical.get_values(), None,
                            float_format=None,
                            na_rep=self.na_rep)

    def to_string(self):
        categorical = self.categorical

        if len(categorical) == 0:
            if self.footer:
                return self._get_footer()
            else:
                return u('')

        fmt_values = self._get_formatted_values()

        result = ['%s' % i for i in fmt_values]
        result = [i.strip() for i in result]
        result = u(', ').join(result)
        result = [u('[')+result+u(']')]
        if self.footer:
            footer = self._get_footer()
            if footer:
                result.append(footer)

        return compat.text_type(u('\n').join(result))


class SeriesFormatter(object):

    def __init__(self, series, buf=None, length=True, header=True,
                 na_rep='NaN', name=False, float_format=None, dtype=True,
                 max_rows=None):
        self.series = series
        self.buf = buf if buf is not None else StringIO()
        self.name = name
        self.na_rep = na_rep
        self.header = header
        self.length = length
        self.max_rows = max_rows

        if float_format is None:
            float_format = get_option("display.float_format")
        self.float_format = float_format
        self.dtype = dtype
        self.adj = _get_adjustment()

        self._chk_truncate()

    def _chk_truncate(self):
        from pandas.tools.merge import concat
        max_rows = self.max_rows
        truncate_v = max_rows and (len(self.series) > max_rows)
        series = self.series
        if truncate_v:
            if max_rows == 1:
                row_num = max_rows
                series = series.iloc[:max_rows]
            else:
                row_num = max_rows // 2
                series = concat((series.iloc[:row_num], series.iloc[-row_num:]))
            self.tr_row_num = row_num
        self.tr_series = series
        self.truncate_v = truncate_v

    def _get_footer(self):
        name = self.series.name
        footer = u('')

        if getattr(self.series.index, 'freq', None) is not None:
            footer += 'Freq: %s' % self.series.index.freqstr

        if self.name is not False and name is not None:
            if footer:
                footer += ', '

            series_name = com.pprint_thing(name,
                                           escape_chars=('\t', '\r', '\n'))
            footer += ("Name: %s" %
                       series_name) if name is not None else ""

        if self.length:
            if footer:
                footer += ', '
            footer += 'Length: %d' % len(self.series)

        if self.dtype is not False and self.dtype is not None:
            name = getattr(self.tr_series.dtype, 'name', None)
            if name:
                if footer:
                    footer += ', '
                footer += 'dtype: %s' % com.pprint_thing(name)

        # level infos are added to the end and in a new line, like it is done for Categoricals
        if com.is_categorical_dtype(self.tr_series.dtype):
            level_info = self.tr_series._values._repr_categories_info()
            if footer:
                footer += "\n"
            footer += level_info

        return compat.text_type(footer)

    def _get_formatted_index(self):
        index = self.tr_series.index
        is_multi = isinstance(index, MultiIndex)

        if is_multi:
            have_header = any(name for name in index.names)
            fmt_index = index.format(names=True)
        else:
            have_header = index.name is not None
            fmt_index = index.format(name=True)
        return fmt_index, have_header

    def _get_formatted_values(self):
        return format_array(self.tr_series._values, None,
                            float_format=self.float_format,
                            na_rep=self.na_rep)

    def to_string(self):
        series = self.tr_series
        footer = self._get_footer()

        if len(series) == 0:
            return 'Series([], ' + footer + ')'

        fmt_index, have_header = self._get_formatted_index()
        fmt_values = self._get_formatted_values()

        maxlen = max(self.adj.len(x) for x in fmt_index)  # max index len
        pad_space = min(maxlen, 60)

        if self.truncate_v:
            n_header_rows = 0
            row_num = self.tr_row_num
            width = self.adj.len(fmt_values[row_num-1])
            if width > 3:
                dot_str = '...'
            else:
                dot_str = '..'
            # Series uses mode=center because it has single value columns
            # DataFrame uses mode=left
            dot_str = self.adj.justify([dot_str], width, mode='center')[0]
            fmt_values.insert(row_num + n_header_rows, dot_str)
            fmt_index.insert(row_num + 1, '')

        result = self.adj.adjoin(3, *[fmt_index[1:], fmt_values])

        if self.header and have_header:
            result = fmt_index[0] + '\n' + result

        if footer:
            result += '\n' + footer

        return compat.text_type(u('').join(result))


class TextAdjustment(object):

    def __init__(self):
        self.encoding = get_option("display.encoding")

    def len(self, text):
        return compat.strlen(text, encoding=self.encoding)

    def justify(self, texts, max_len, mode='right'):
        return com._justify(texts, max_len, mode=mode)

    def adjoin(self, space, *lists, **kwargs):
        return com.adjoin(space, *lists, strlen=self.len,
                          justfunc=self.justify, **kwargs)


class EastAsianTextAdjustment(TextAdjustment):

    def __init__(self):
        super(EastAsianTextAdjustment, self).__init__()
        if get_option("display.unicode.ambiguous_as_wide"):
            self.ambiguous_width = 2
        else:
            self.ambiguous_width = 1

    def len(self, text):
        return compat.east_asian_len(text, encoding=self.encoding,
                                     ambiguous_width=self.ambiguous_width)

    def justify(self, texts, max_len, mode='right'):
        # re-calculate padding space per str considering East Asian Width
        def _get_pad(t):
            return max_len - self.len(t) + len(t)

        if mode == 'left':
            return [x.ljust(_get_pad(x)) for x in texts]
        elif mode == 'center':
            return [x.center(_get_pad(x)) for x in texts]
        else:
            return [x.rjust(_get_pad(x)) for x in texts]


def _get_adjustment():
    use_east_asian_width = get_option("display.unicode.east_asian_width")
    if use_east_asian_width:
        return EastAsianTextAdjustment()
    else:
        return TextAdjustment()


class TableFormatter(object):
    is_truncated = False
    show_dimensions = None

    @property
    def should_show_dimensions(self):
        return self.show_dimensions is True or (self.show_dimensions == 'truncate' and
                                                self.is_truncated)

    def _get_formatter(self, i):
        if isinstance(self.formatters, (list, tuple)):
            if com.is_integer(i):
                return self.formatters[i]
            else:
                return None
        else:
            if com.is_integer(i) and i not in self.columns:
                i = self.columns[i]
            return self.formatters.get(i, None)


class DataFrameFormatter(TableFormatter):

    """
    Render a DataFrame

    self.to_string() : console-friendly tabular output
    self.to_html()   : html table
    self.to_latex()   : LaTeX tabular environment table

    """

    __doc__ = __doc__ if __doc__ else ''
    __doc__ += common_docstring + justify_docstring + return_docstring

    def __init__(self, frame, buf=None, columns=None, col_space=None,
                 header=True, index=True, na_rep='NaN', formatters=None,
                 justify=None, float_format=None, sparsify=None,
                 index_names=True, line_width=None, max_rows=None,
                 max_cols=None, show_dimensions=False, **kwds):
        self.frame = frame
        self.buf = buf if buf is not None else StringIO()
        self.show_index_names = index_names

        if sparsify is None:
            sparsify = get_option("display.multi_sparse")

        self.sparsify = sparsify

        self.float_format = float_format
        self.formatters = formatters if formatters is not None else {}
        self.na_rep = na_rep
        self.col_space = col_space
        self.header = header
        self.index = index
        self.line_width = line_width
        self.max_rows = max_rows
        self.max_cols = max_cols
        self.max_rows_displayed = min(max_rows or len(self.frame),
                                      len(self.frame))
        self.show_dimensions = show_dimensions

        if justify is None:
            self.justify = get_option("display.colheader_justify")
        else:
            self.justify = justify

        self.kwds = kwds

        if columns is not None:
            self.columns = _ensure_index(columns)
            self.frame = self.frame[self.columns]
        else:
            self.columns = frame.columns

        self._chk_truncate()
        self.adj = _get_adjustment()

    def _chk_truncate(self):
        '''
        Checks whether the frame should be truncated. If so, slices
        the frame up.
        '''
        from pandas.tools.merge import concat

        # Column of which first element is used to determine width of a dot col
        self.tr_size_col = -1

        # Cut the data to the information actually printed
        max_cols = self.max_cols
        max_rows = self.max_rows

        if max_cols == 0 or max_rows == 0:  # assume we are in the terminal (why else = 0)
            (w, h) = get_terminal_size()
            self.w = w
            self.h = h
            if self.max_rows == 0:
                dot_row = 1
                prompt_row = 1
                if self.show_dimensions:
                    show_dimension_rows = 3
                n_add_rows = self.header + dot_row + show_dimension_rows + prompt_row
                max_rows_adj = self.h - n_add_rows  # rows available to fill with actual data
                self.max_rows_adj = max_rows_adj

            # Format only rows and columns that could potentially fit the screen
            if max_cols == 0 and len(self.frame.columns) > w:
                max_cols = w
            if max_rows == 0 and len(self.frame) > h:
                max_rows = h

        if not hasattr(self, 'max_rows_adj'):
            self.max_rows_adj = max_rows
        if not hasattr(self, 'max_cols_adj'):
            self.max_cols_adj = max_cols

        max_cols_adj = self.max_cols_adj
        max_rows_adj = self.max_rows_adj

        truncate_h = max_cols_adj and (len(self.columns) > max_cols_adj)
        truncate_v = max_rows_adj and (len(self.frame) > max_rows_adj)

        frame = self.frame
        if truncate_h:
            if max_cols_adj == 0:
                col_num = len(frame.columns)
            elif max_cols_adj == 1:
                frame = frame.iloc[:, :max_cols]
                col_num = max_cols
            else:
                col_num = (max_cols_adj // 2)
                frame = concat((frame.iloc[:, :col_num], frame.iloc[:, -col_num:]), axis=1)
            self.tr_col_num = col_num
        if truncate_v:
            if max_rows_adj == 0:
                row_num = len(frame)
            if max_rows_adj == 1:
                row_num = max_rows
                frame = frame.iloc[:max_rows, :]
            else:
                row_num = max_rows_adj // 2
                frame = concat((frame.iloc[:row_num, :], frame.iloc[-row_num:, :]))
            self.tr_row_num = row_num

        self.tr_frame = frame
        self.truncate_h = truncate_h
        self.truncate_v = truncate_v
        self.is_truncated = self.truncate_h or self.truncate_v

    def _to_str_columns(self):
        """
        Render a DataFrame to a list of columns (as lists of strings).
        """
        frame = self.tr_frame

        # may include levels names also

        str_index = self._get_formatted_index(frame)
        str_columns = self._get_formatted_column_labels(frame)

        if self.header:
            stringified = []
            for i, c in enumerate(frame):
                cheader = str_columns[i]
                max_colwidth = max(self.col_space or 0,
                                   *(self.adj.len(x) for x in cheader))
                fmt_values = self._format_col(i)
                fmt_values = _make_fixed_width(fmt_values, self.justify,
                                               minimum=max_colwidth,
                                               adj=self.adj)

                max_len = max(np.max([self.adj.len(x) for x in fmt_values]),
                              max_colwidth)
                cheader = self.adj.justify(cheader, max_len, mode=self.justify)
                stringified.append(cheader + fmt_values)
        else:
            stringified = []
            for i, c in enumerate(frame):
                fmt_values = self._format_col(i)
                fmt_values = _make_fixed_width(fmt_values, self.justify,
                                               minimum=(self.col_space or 0),
                                               adj=self.adj)

                stringified.append(fmt_values)

        strcols = stringified
        if self.index:
            strcols.insert(0, str_index)

        # Add ... to signal truncated
        truncate_h = self.truncate_h
        truncate_v = self.truncate_v

        if truncate_h:
            col_num = self.tr_col_num
            col_width = self.adj.len(strcols[self.tr_size_col][0])  # infer from column header
            strcols.insert(self.tr_col_num + 1, ['...'.center(col_width)] * (len(str_index)))
        if truncate_v:
            n_header_rows = len(str_index) - len(frame)
            row_num = self.tr_row_num
            for ix, col in enumerate(strcols):
                cwidth = self.adj.len(strcols[ix][row_num])  # infer from above row
                is_dot_col = False
                if truncate_h:
                    is_dot_col = ix == col_num + 1
                if cwidth > 3 or is_dot_col:
                    my_str = '...'
                else:
                    my_str = '..'

                if ix == 0:
                    dot_mode = 'left'
                elif is_dot_col:
                    cwidth = self.adj.len(strcols[self.tr_size_col][0])
                    dot_mode = 'center'
                else:
                    dot_mode = 'right'
                dot_str = self.adj.justify([my_str], cwidth, mode=dot_mode)[0]
                strcols[ix].insert(row_num + n_header_rows, dot_str)
        return strcols

    def to_string(self):
        """
        Render a DataFrame to a console-friendly tabular output.
        """
        from pandas import Series

        frame = self.frame

        if len(frame.columns) == 0 or len(frame.index) == 0:
            info_line = (u('Empty %s\nColumns: %s\nIndex: %s')
                         % (type(self.frame).__name__,
                            com.pprint_thing(frame.columns),
                            com.pprint_thing(frame.index)))
            text = info_line
        else:
            strcols = self._to_str_columns()
            if self.line_width is None:  # no need to wrap around just print the whole frame
                text = self.adj.adjoin(1, *strcols)
            elif not isinstance(self.max_cols, int) or self.max_cols > 0:  # need to wrap around
                text = self._join_multiline(*strcols)
            else:  # max_cols == 0. Try to fit frame to terminal
                text = self.adj.adjoin(1, *strcols).split('\n')
                row_lens = Series(text).apply(len)
                max_len_col_ix = np.argmax(row_lens)
                max_len = row_lens[max_len_col_ix]
                headers = [ele[0] for ele in strcols]
                # Size of last col determines dot col size. See `self._to_str_columns
                size_tr_col = len(headers[self.tr_size_col])
                max_len += size_tr_col  # Need to make space for largest row plus truncate dot col
                dif = max_len - self.w
                adj_dif = dif
                col_lens = Series([Series(ele).apply(len).max() for ele in strcols])
                n_cols = len(col_lens)
                counter = 0
                while adj_dif > 0 and n_cols > 1:
                    counter += 1
                    mid = int(round(n_cols / 2.))
                    mid_ix = col_lens.index[mid]
                    col_len = col_lens[mid_ix]
                    adj_dif -= (col_len + 1)  # adjoin adds one
                    col_lens = col_lens.drop(mid_ix)
                    n_cols = len(col_lens)
                max_cols_adj = n_cols - self.index  # subtract index column
                self.max_cols_adj = max_cols_adj

                # Call again _chk_truncate to cut frame appropriately
                # and then generate string representation
                self._chk_truncate()
                strcols = self._to_str_columns()
                text = self.adj.adjoin(1, *strcols)

        self.buf.writelines(text)

        if self.should_show_dimensions:
            self.buf.write("\n\n[%d rows x %d columns]"
                           % (len(frame), len(frame.columns)))

    def _join_multiline(self, *strcols):
        lwidth = self.line_width
        adjoin_width = 1
        strcols = list(strcols)
        if self.index:
            idx = strcols.pop(0)
            lwidth -= np.array([self.adj.len(x) for x in idx]).max() + adjoin_width

        col_widths = [np.array([self.adj.len(x) for x in col]).max()
                      if len(col) > 0 else 0
                      for col in strcols]
        col_bins = _binify(col_widths, lwidth)
        nbins = len(col_bins)

        if self.truncate_v:
            nrows = self.max_rows_adj + 1
        else:
            nrows = len(self.frame)

        str_lst = []
        st = 0
        for i, ed in enumerate(col_bins):
            row = strcols[st:ed]
            row.insert(0, idx)
            if nbins > 1:
                if ed <= len(strcols) and i < nbins - 1:
                    row.append([' \\'] + ['  '] * (nrows - 1))
                else:
                    row.append([' '] * nrows)
            str_lst.append(self.adj.adjoin(adjoin_width, *row))
            st = ed
        return '\n\n'.join(str_lst)

    def to_latex(self, column_format=None, longtable=False):
        """
        Render a DataFrame to a LaTeX tabular/longtable environment output.
        """
        self.escape = self.kwds.get('escape', True)

        def get_col_type(dtype):
            if issubclass(dtype.type, np.number):
                return 'r'
            else:
                return 'l'

        frame = self.frame

        if len(frame.columns) == 0 or len(frame.index) == 0:
            info_line = (u('Empty %s\nColumns: %s\nIndex: %s')
                         % (type(self.frame).__name__,
                            frame.columns, frame.index))
            strcols = [[info_line]]
        else:
            strcols = self._to_str_columns()

        if self.index and isinstance(self.frame.index, MultiIndex):
            clevels = self.frame.columns.nlevels
            strcols.pop(0)
            name = any(self.frame.columns.names)
            for i, lev in enumerate(self.frame.index.levels):
                lev2 = lev.format(name=name)
                blank = ' ' * len(lev2[0])
                lev3 = [blank] * clevels
                for level_idx, group in itertools.groupby(
                        self.frame.index.labels[i]):
                    count = len(list(group))
                    lev3.extend([lev2[level_idx]] + [blank] * (count - 1))
                strcols.insert(i, lev3)

        if column_format is None:
            dtypes = self.frame.dtypes._values
            column_format = ''.join(map(get_col_type, dtypes))
            if self.index:
                index_format = 'l' * self.frame.index.nlevels
                column_format = index_format + column_format
        elif not isinstance(column_format,
                            compat.string_types):  # pragma: no cover
            raise AssertionError('column_format must be str or unicode, not %s'
                                 % type(column_format))

        def write(buf, frame, column_format, strcols, longtable=False):
            if not longtable:
                buf.write('\\begin{tabular}{%s}\n' % column_format)
                buf.write('\\toprule\n')
            else:
                buf.write('\\begin{longtable}{%s}\n' % column_format)
                buf.write('\\toprule\n')

            nlevels = frame.columns.nlevels
            for i, row in enumerate(zip(*strcols)):
                if i == nlevels:
                    buf.write('\\midrule\n')  # End of header
                    if longtable:
                        buf.write('\\endhead\n')
                        buf.write('\\midrule\n')
                        buf.write('\\multicolumn{3}{r}{{Continued on next '
                                  'page}} \\\\\n')
                        buf.write('\midrule\n')
                        buf.write('\endfoot\n\n')
                        buf.write('\\bottomrule\n')
                        buf.write('\\endlastfoot\n')
                if self.escape:
                    crow = [(x.replace('\\', '\\textbackslash')  # escape backslashes first
                             .replace('_', '\\_')
                             .replace('%', '\\%')
                             .replace('$', '\\$')
                             .replace('#', '\\#')
                             .replace('{', '\\{')
                             .replace('}', '\\}')
                             .replace('~', '\\textasciitilde')
                             .replace('^', '\\textasciicircum')
                             .replace('&', '\\&') if x else '{}') for x in row]
                else:
                    crow = [x if x else '{}' for x in row]
                buf.write(' & '.join(crow))
                buf.write(' \\\\\n')

            if not longtable:
                buf.write('\\bottomrule\n')
                buf.write('\\end{tabular}\n')
            else:
                buf.write('\\end{longtable}\n')

        if hasattr(self.buf, 'write'):
            write(self.buf, frame, column_format, strcols, longtable)
        elif isinstance(self.buf, compat.string_types):
            with open(self.buf, 'w') as f:
                write(f, frame, column_format, strcols, longtable)
        else:
            raise TypeError('buf is not a file name and it has no write '
                            'method')

    def _format_col(self, i):
        frame = self.tr_frame
        formatter = self._get_formatter(i)
        return format_array(
            frame.iloc[:, i]._values,
            formatter, float_format=self.float_format, na_rep=self.na_rep,
            space=self.col_space
        )

    def to_html(self, classes=None, notebook=False):
        """
        Render a DataFrame to a html table.

        Parameters
        ----------
        notebook : {True, False}, optional, default False
            Whether the generated HTML is for IPython Notebook.

        """
        html_renderer = HTMLFormatter(self, classes=classes,
                                      max_rows=self.max_rows,
                                      max_cols=self.max_cols,
                                      notebook=notebook)
        if hasattr(self.buf, 'write'):
            html_renderer.write_result(self.buf)
        elif isinstance(self.buf, compat.string_types):
            with open(self.buf, 'w') as f:
                html_renderer.write_result(f)
        else:
            raise TypeError('buf is not a file name and it has no write '
                            ' method')

    def _get_formatted_column_labels(self, frame):
        from pandas.core.index import _sparsify

        def is_numeric_dtype(dtype):
            return issubclass(dtype.type, np.number)

        columns = frame.columns

        if isinstance(columns, MultiIndex):
            fmt_columns = columns.format(sparsify=False, adjoin=False)
            fmt_columns = lzip(*fmt_columns)
            dtypes = self.frame.dtypes._values

            # if we have a Float level, they don't use leading space at all
            restrict_formatting = any([l.is_floating for l in columns.levels])
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))

            def space_format(x, y):
                if y not in self.formatters and need_leadsp[x] and not restrict_formatting:
                    return ' ' + y
                return y

            str_columns = list(zip(*[[space_format(x, y) for y in x] for x in fmt_columns]))
            if self.sparsify:
                str_columns = _sparsify(str_columns)

            str_columns = [list(x) for x in zip(*str_columns)]
        else:
            fmt_columns = columns.format()
            dtypes = self.frame.dtypes
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))
            str_columns = [[' ' + x
                            if not self._get_formatter(i) and need_leadsp[x]
                            else x]
                           for i, (col, x) in
                           enumerate(zip(columns, fmt_columns))]

        if self.show_index_names and self.has_index_names:
            for x in str_columns:
                x.append('')

        # self.str_columns = str_columns
        return str_columns

    @property
    def has_index_names(self):
        return _has_names(self.frame.index)

    @property
    def has_column_names(self):
        return _has_names(self.frame.columns)

    def _get_formatted_index(self, frame):
        # Note: this is only used by to_string() and to_latex(), not by to_html().
        index = frame.index
        columns = frame.columns

        show_index_names = self.show_index_names and self.has_index_names
        show_col_names = (self.show_index_names and self.has_column_names)

        fmt = self._get_formatter('__index__')

        if isinstance(index, MultiIndex):
            fmt_index = index.format(sparsify=self.sparsify, adjoin=False,
                                     names=show_index_names,
                                     formatter=fmt)
        else:
            fmt_index = [index.format(name=show_index_names, formatter=fmt)]
        fmt_index = [tuple(_make_fixed_width(list(x), justify='left',
                                             minimum=(self.col_space or 0),
                                             adj=self.adj))
                     for x in fmt_index]

        adjoined = self.adj.adjoin(1, *fmt_index).split('\n')

        # empty space for columns
        if show_col_names:
            col_header = ['%s' % x for x in self._get_column_name_list()]
        else:
            col_header = [''] * columns.nlevels

        if self.header:
            return col_header + adjoined
        else:
            return adjoined

    def _get_column_name_list(self):
        names = []
        columns = self.frame.columns
        if isinstance(columns, MultiIndex):
            names.extend('' if name is None else name
                         for name in columns.names)
        else:
            names.append('' if columns.name is None else columns.name)
        return names


class HTMLFormatter(TableFormatter):

    indent_delta = 2

    def __init__(self, formatter, classes=None, max_rows=None, max_cols=None,
                 notebook=False):
        self.fmt = formatter
        self.classes = classes

        self.frame = self.fmt.frame
        self.columns = self.fmt.tr_frame.columns
        self.elements = []
        self.bold_rows = self.fmt.kwds.get('bold_rows', False)
        self.escape = self.fmt.kwds.get('escape', True)

        self.max_rows = max_rows or len(self.fmt.frame)
        self.max_cols = max_cols or len(self.fmt.columns)
        self.show_dimensions = self.fmt.show_dimensions
        self.is_truncated = (self.max_rows < len(self.fmt.frame) or
                             self.max_cols < len(self.fmt.columns))
        self.notebook = notebook

    def write(self, s, indent=0):
        rs = com.pprint_thing(s)
        self.elements.append(' ' * indent + rs)

    def write_th(self, s, indent=0, tags=None):
        if (self.fmt.col_space is not None
                and self.fmt.col_space > 0):
            tags = (tags or "")
            tags += 'style="min-width: %s;"' % self.fmt.col_space

        return self._write_cell(s, kind='th', indent=indent, tags=tags)

    def write_td(self, s, indent=0, tags=None):
        return self._write_cell(s, kind='td', indent=indent, tags=tags)

    def _write_cell(self, s, kind='td', indent=0, tags=None):
        if tags is not None:
            start_tag = '<%s %s>' % (kind, tags)
        else:
            start_tag = '<%s>' % kind

        if self.escape:
            # escape & first to prevent double escaping of &
            esc = OrderedDict(
                [('&', r'&amp;'), ('<', r'&lt;'), ('>', r'&gt;')]
            )
        else:
            esc = {}
        rs = com.pprint_thing(s, escape_chars=esc).strip()
        self.write(
            '%s%s</%s>' % (start_tag, rs, kind), indent)

    def write_tr(self, line, indent=0, indent_delta=4, header=False,
                 align=None, tags=None, nindex_levels=0):
        if tags is None:
            tags = {}

        if align is None:
            self.write('<tr>', indent)
        else:
            self.write('<tr style="text-align: %s;">' % align, indent)
        indent += indent_delta

        for i, s in enumerate(line):
            val_tag = tags.get(i, None)
            if header or (self.bold_rows and i < nindex_levels):
                self.write_th(s, indent, tags=val_tag)
            else:
                self.write_td(s, indent, tags=val_tag)

        indent -= indent_delta
        self.write('</tr>', indent)

    def write_result(self, buf):
        indent = 0
        frame = self.frame

        _classes = ['dataframe']  # Default class.
        if self.classes is not None:
            if isinstance(self.classes, str):
                self.classes = self.classes.split()
            if not isinstance(self.classes, (list, tuple)):
                raise AssertionError(('classes must be list or tuple, '
                                      'not %s') % type(self.classes))
            _classes.extend(self.classes)

        if self.notebook:
            div_style = ''
            try:
                import IPython
                if IPython.__version__ < LooseVersion('3.0.0'):
                    div_style = ' style="max-width:1500px;overflow:auto;"'
            except ImportError:
                pass

            self.write('<div{0}>'.format(div_style))

        self.write('<table border="1" class="%s">' % ' '.join(_classes),
                   indent)

        indent += self.indent_delta
        indent = self._write_header(indent)
        indent = self._write_body(indent)

        self.write('</table>', indent)
        if self.should_show_dimensions:
            by = chr(215) if compat.PY3 else unichr(215)  # ×
            self.write(u('<p>%d rows %s %d columns</p>') %
                       (len(frame), by, len(frame.columns)))

        if self.notebook:
            self.write('</div>')

        _put_lines(buf, self.elements)

    def _write_header(self, indent):
        truncate_h = self.fmt.truncate_h
        row_levels = self.frame.index.nlevels
        if not self.fmt.header:
            # write nothing
            return indent

        def _column_header():
            if self.fmt.index:
                row = [''] * (self.frame.index.nlevels - 1)
            else:
                row = []

            if isinstance(self.columns, MultiIndex):
                if self.fmt.has_column_names and self.fmt.index:
                    row.append(single_column_table(self.columns.names))
                else:
                    row.append('')
                style = "text-align: %s;" % self.fmt.justify
                row.extend([single_column_table(c, self.fmt.justify, style) for
                            c in self.columns])
            else:
                if self.fmt.index:
                    row.append(self.columns.name or '')
                row.extend(self.columns)
            return row

        self.write('<thead>', indent)
        row = []

        indent += self.indent_delta

        if isinstance(self.columns, MultiIndex):
            template = 'colspan="%d" halign="left"'

            if self.fmt.sparsify:
                # GH3547
                sentinel = com.sentinel_factory()
            else:
                sentinel = None
            levels = self.columns.format(sparsify=sentinel,
                                         adjoin=False, names=False)
            level_lengths = _get_level_lengths(levels, sentinel)
            inner_lvl = len(level_lengths) - 1
            for lnum, (records, values) in enumerate(zip(level_lengths,
                                                         levels)):
                if truncate_h:
                    # modify the header lines
                    ins_col = self.fmt.tr_col_num
                    if self.fmt.sparsify:
                        recs_new = {}
                        # Increment tags after ... col.
                        for tag, span in list(records.items()):
                            if tag >= ins_col:
                                recs_new[tag + 1] = span
                            elif tag + span > ins_col:
                                recs_new[tag] = span + 1
                                if lnum == inner_lvl:
                                    values = values[:ins_col] + (u('...'),) + \
                                        values[ins_col:]
                                else:  # sparse col headers do not receive a ...
                                    values = (values[:ins_col] + (values[ins_col - 1],) +
                                              values[ins_col:])
                            else:
                                recs_new[tag] = span
                            # if ins_col lies between tags, all col headers get ...
                            if tag + span == ins_col:
                                recs_new[ins_col] = 1
                                values = values[:ins_col] + (u('...'),) + \
                                    values[ins_col:]
                        records = recs_new
                        inner_lvl = len(level_lengths) - 1
                        if lnum == inner_lvl:
                            records[ins_col] = 1
                    else:
                        recs_new = {}
                        for tag, span in list(records.items()):
                            if tag >= ins_col:
                                recs_new[tag + 1] = span
                            else:
                                recs_new[tag] = span
                        recs_new[ins_col] = 1
                        records = recs_new
                        values = values[:ins_col] + [u('...')] + values[ins_col:]

                name = self.columns.names[lnum]
                row = [''] * (row_levels - 1) + ['' if name is None
                                                 else com.pprint_thing(name)]

                if row == [""] and self.fmt.index is False:
                    row = []

                tags = {}
                j = len(row)
                for i, v in enumerate(values):
                    if i in records:
                        if records[i] > 1:
                            tags[j] = template % records[i]
                    else:
                        continue
                    j += 1
                    row.append(v)
                self.write_tr(row, indent, self.indent_delta, tags=tags,
                              header=True)
        else:
            col_row = _column_header()
            align = self.fmt.justify

            if truncate_h:
                ins_col = row_levels + self.fmt.tr_col_num
                col_row.insert(ins_col, '...')

            self.write_tr(col_row, indent, self.indent_delta, header=True,
                          align=align)

        if self.fmt.has_index_names and self.fmt.index:
            row = [
                x if x is not None else '' for x in self.frame.index.names
            ] + [''] * min(len(self.columns), self.max_cols)
            if truncate_h:
                ins_col = row_levels + self.fmt.tr_col_num
                row.insert(ins_col, '')
            self.write_tr(row, indent, self.indent_delta, header=True)

        indent -= self.indent_delta
        self.write('</thead>', indent)

        return indent

    def _write_body(self, indent):
        self.write('<tbody>', indent)
        indent += self.indent_delta

        fmt_values = {}
        for i in range(min(len(self.columns), self.max_cols)):
            fmt_values[i] = self.fmt._format_col(i)

        # write values
        if self.fmt.index:
            if isinstance(self.frame.index, MultiIndex):
                self._write_hierarchical_rows(fmt_values, indent)
            else:
                self._write_regular_rows(fmt_values, indent)
        else:
            for i in range(len(self.frame)):
                row = [fmt_values[j][i] for j in range(len(self.columns))]
                self.write_tr(row, indent, self.indent_delta, tags=None)

        indent -= self.indent_delta
        self.write('</tbody>', indent)
        indent -= self.indent_delta

        return indent

    def _write_regular_rows(self, fmt_values, indent):
        truncate_h = self.fmt.truncate_h
        truncate_v = self.fmt.truncate_v

        ncols = len(self.fmt.tr_frame.columns)
        nrows = len(self.fmt.tr_frame)
        fmt = self.fmt._get_formatter('__index__')
        if fmt is not None:
            index_values = self.fmt.tr_frame.index.map(fmt)
        else:
            index_values = self.fmt.tr_frame.index.format()

        row = []
        for i in range(nrows):

            if truncate_v and i == (self.fmt.tr_row_num):
                str_sep_row = ['...' for ele in row]
                self.write_tr(str_sep_row, indent, self.indent_delta, tags=None,
                              nindex_levels=1)

            row = []
            row.append(index_values[i])
            row.extend(fmt_values[j][i] for j in range(ncols))

            if truncate_h:
                dot_col_ix = self.fmt.tr_col_num + 1
                row.insert(dot_col_ix, '...')
            self.write_tr(row, indent, self.indent_delta, tags=None,
                          nindex_levels=1)

    def _write_hierarchical_rows(self, fmt_values, indent):
        template = 'rowspan="%d" valign="top"'

        truncate_h = self.fmt.truncate_h
        truncate_v = self.fmt.truncate_v
        frame = self.fmt.tr_frame
        ncols = len(frame.columns)
        nrows = len(frame)
        row_levels = self.frame.index.nlevels

        idx_values = frame.index.format(sparsify=False, adjoin=False, names=False)
        idx_values = lzip(*idx_values)

        if self.fmt.sparsify:
            # GH3547
            sentinel = com.sentinel_factory()
            levels = frame.index.format(sparsify=sentinel, adjoin=False, names=False)

            level_lengths = _get_level_lengths(levels, sentinel)
            inner_lvl = len(level_lengths) - 1
            if truncate_v:
                # Insert ... row and adjust idx_values and
                # level_lengths to take this into account.
                ins_row = self.fmt.tr_row_num
                for lnum, records in enumerate(level_lengths):
                    rec_new = {}
                    for tag, span in list(records.items()):
                        if tag >= ins_row:
                            rec_new[tag + 1] = span
                        elif tag + span > ins_row:
                            rec_new[tag] = span + 1
                            dot_row = list(idx_values[ins_row - 1])
                            dot_row[-1] = u('...')
                            idx_values.insert(ins_row, tuple(dot_row))
                        else:
                            rec_new[tag] = span
                        # If ins_row lies between tags, all cols idx cols receive ...
                        if tag + span == ins_row:
                            rec_new[ins_row] = 1
                            if lnum == 0:
                                idx_values.insert(ins_row, tuple([u('...')]*len(level_lengths)))
                    level_lengths[lnum] = rec_new

                level_lengths[inner_lvl][ins_row] = 1
                for ix_col in range(len(fmt_values)):
                    fmt_values[ix_col].insert(ins_row, '...')
                nrows += 1

            for i in range(nrows):
                row = []
                tags = {}

                sparse_offset = 0
                j = 0
                for records, v in zip(level_lengths, idx_values[i]):
                    if i in records:
                        if records[i] > 1:
                            tags[j] = template % records[i]
                    else:
                        sparse_offset += 1
                        continue

                    j += 1
                    row.append(v)

                row.extend(fmt_values[j][i] for j in range(ncols))
                if truncate_h:
                    row.insert(row_levels - sparse_offset + self.fmt.tr_col_num, '...')
                self.write_tr(row, indent, self.indent_delta, tags=tags,
                              nindex_levels=len(levels) - sparse_offset)
        else:
            for i in range(len(frame)):
                idx_values = list(zip(*frame.index.format(sparsify=False,
                                                          adjoin=False,
                                                          names=False)))
                row = []
                row.extend(idx_values[i])
                row.extend(fmt_values[j][i] for j in range(ncols))
                if truncate_h:
                    row.insert(row_levels + self.fmt.tr_col_num, '...')
                self.write_tr(row, indent, self.indent_delta, tags=None,
                              nindex_levels=frame.index.nlevels)


def _get_level_lengths(levels, sentinel=''):
    from itertools import groupby

    def _make_grouper():
        record = {'count': 0}

        def grouper(x):
            if x != sentinel:
                record['count'] += 1
            return record['count']
        return grouper

    result = []
    for lev in levels:
        i = 0
        f = _make_grouper()
        recs = {}
        for key, gpr in groupby(lev, f):
            values = list(gpr)
            recs[i] = len(values)
            i += len(values)

        result.append(recs)

    return result


class CSVFormatter(object):

    def __init__(self, obj, path_or_buf=None, sep=",", na_rep='', float_format=None,
                 cols=None, header=True, index=True, index_label=None,
                 mode='w', nanRep=None, encoding=None, quoting=None,
                 line_terminator='\n', chunksize=None, engine=None,
                 tupleize_cols=False, quotechar='"', date_format=None,
                 doublequote=True, escapechar=None, decimal='.'):

        self.engine = engine  # remove for 0.13
        self.obj = obj

        if path_or_buf is None:
            path_or_buf = StringIO()

        self.path_or_buf = path_or_buf
        self.sep = sep
        self.na_rep = na_rep
        self.float_format = float_format
        self.decimal = decimal

        self.header = header
        self.index = index
        self.index_label = index_label
        self.mode = mode
        self.encoding = encoding

        if quoting is None:
            quoting = csv.QUOTE_MINIMAL
        self.quoting = quoting

        if quoting == csv.QUOTE_NONE:
            # prevents crash in _csv
            quotechar = None
        self.quotechar = quotechar

        self.doublequote = doublequote
        self.escapechar = escapechar

        self.line_terminator = line_terminator

        self.date_format = date_format

        # GH3457
        if not self.obj.columns.is_unique and engine == 'python':
            raise NotImplementedError("columns.is_unique == False not "
                                      "supported with engine='python'")

        self.tupleize_cols = tupleize_cols
        self.has_mi_columns = isinstance(obj.columns, MultiIndex
                                         ) and not self.tupleize_cols

        # validate mi options
        if self.has_mi_columns:
            if cols is not None:
                raise TypeError("cannot specify cols with a MultiIndex on the "
                                "columns")

        if cols is not None:
            if isinstance(cols, Index):
                cols = cols.to_native_types(na_rep=na_rep,
                                            float_format=float_format,
                                            date_format=date_format,
                                            quoting=self.quoting)
            else:
                cols = np.asarray(list(cols))
            self.obj = self.obj.loc[:, cols]

        # update columns to include possible multiplicity of dupes
        # and make sure sure cols is just a list of labels
        cols = self.obj.columns
        if isinstance(cols, Index):
            cols = cols.to_native_types(na_rep=na_rep,
                                        float_format=float_format,
                                        date_format=date_format,
                                        quoting=self.quoting)
        else:
            cols = np.asarray(list(cols))

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
        if isinstance(obj.index, PeriodIndex):
            self.data_index = obj.index.to_timestamp()

        if (isinstance(self.data_index, DatetimeIndex) and
                date_format is not None):
            self.data_index = Index([x.strftime(date_format)
                                     if notnull(x) else ''
                                     for x in self.data_index])

        self.nlevels = getattr(self.data_index, 'nlevels', 1)
        if not index:
            self.nlevels = 0

    # original python implem. of df.to_csv
    # invoked by df.to_csv(engine=python)
    def _helper_csv(self, writer, na_rep=None, cols=None,
                    header=True, index=True,
                    index_label=None, float_format=None, date_format=None):
        if cols is None:
            cols = self.columns

        has_aliases = isinstance(header, (tuple, list, np.ndarray, Index))
        if has_aliases or header:
            if index:
                # should write something for index label
                if index_label is not False:
                    if index_label is None:
                        if isinstance(self.obj.index, MultiIndex):
                            index_label = []
                            for i, name in enumerate(self.obj.index.names):
                                if name is None:
                                    name = ''
                                index_label.append(name)
                        else:
                            index_label = self.obj.index.name
                            if index_label is None:
                                index_label = ['']
                            else:
                                index_label = [index_label]
                    elif not isinstance(index_label,
                                        (list, tuple, np.ndarray, Index)):
                        # given a string for a DF with Index
                        index_label = [index_label]

                    encoded_labels = list(index_label)
                else:
                    encoded_labels = []

                if has_aliases:
                    if len(header) != len(cols):
                        raise ValueError(('Writing %d cols but got %d aliases'
                                          % (len(cols), len(header))))
                    else:
                        write_cols = header
                else:
                    write_cols = cols
                encoded_cols = list(write_cols)

                writer.writerow(encoded_labels + encoded_cols)
            else:
                encoded_cols = list(cols)
                writer.writerow(encoded_cols)

        if date_format is None:
            date_formatter = lambda x: Timestamp(x)._repr_base
        else:
            def strftime_with_nulls(x):
                x = Timestamp(x)
                if notnull(x):
                    return x.strftime(date_format)

            date_formatter = lambda x: strftime_with_nulls(x)

        data_index = self.obj.index

        if isinstance(self.obj.index, PeriodIndex):
            data_index = self.obj.index.to_timestamp()

        if isinstance(data_index, DatetimeIndex) and date_format is not None:
            data_index = Index([date_formatter(x) for x in data_index])

        values = self.obj.copy()
        values.index = data_index
        values.columns = values.columns.to_native_types(
            na_rep=na_rep,
            float_format=float_format,
            date_format=date_format,
            quoting=self.quoting)
        values = values[cols]

        series = {}
        for k, v in compat.iteritems(values._series):
            series[k] = v._values

        nlevels = getattr(data_index, 'nlevels', 1)
        for j, idx in enumerate(data_index):
            row_fields = []
            if index:
                if nlevels == 1:
                    row_fields = [idx]
                else:  # handle MultiIndex
                    row_fields = list(idx)
            for i, col in enumerate(cols):
                val = series[col][j]
                if lib.checknull(val):
                    val = na_rep

                if float_format is not None and com.is_float(val):
                    val = float_format % val
                elif isinstance(val, (np.datetime64, Timestamp)):
                    val = date_formatter(val)

                row_fields.append(val)

            writer.writerow(row_fields)

    def save(self):
        # create the writer & save
        if hasattr(self.path_or_buf, 'write'):
            f = self.path_or_buf
            close = False
        else:
            f = com._get_handle(self.path_or_buf, self.mode,
                                encoding=self.encoding)
            close = True

        try:
            writer_kwargs = dict(lineterminator=self.line_terminator,
                                 delimiter=self.sep, quoting=self.quoting,
                                 doublequote=self.doublequote,
                                 escapechar=self.escapechar,
                                 quotechar=self.quotechar)
            if self.encoding is not None:
                writer_kwargs['encoding'] = self.encoding
                self.writer = com.UnicodeWriter(f, **writer_kwargs)
            else:
                self.writer = csv.writer(f, **writer_kwargs)

            if self.engine == 'python':
                # to be removed in 0.13
                self._helper_csv(self.writer, na_rep=self.na_rep,
                                 float_format=self.float_format,
                                 cols=self.cols, header=self.header,
                                 index=self.index,
                                 index_label=self.index_label,
                                 date_format=self.date_format)

            else:
                self._save()

        finally:
            if close:
                f.close()

    def _save_header(self):

        writer = self.writer
        obj = self.obj
        index_label = self.index_label
        cols = self.cols
        has_mi_columns = self.has_mi_columns
        header = self.header
        encoded_labels = []

        has_aliases = isinstance(header, (tuple, list, np.ndarray, Index))
        if not (has_aliases or self.header):
            return
        if has_aliases:
            if len(header) != len(cols):
                raise ValueError(('Writing %d cols but got %d aliases'
                                  % (len(cols), len(header))))
            else:
                write_cols = header
        else:
            write_cols = cols

        if self.index:
            # should write something for index label
            if index_label is not False:
                if index_label is None:
                    if isinstance(obj.index, MultiIndex):
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
                elif not isinstance(index_label, (list, tuple, np.ndarray, Index)):
                    # given a string for a DF with Index
                    index_label = [index_label]

                encoded_labels = list(index_label)
            else:
                encoded_labels = []

        if not has_mi_columns:
            encoded_labels += list(write_cols)

        # write out the mi
        if has_mi_columns:
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

                col_line.extend(columns.get_level_values(i))

                writer.writerow(col_line)

            # add blanks for the columns, so that we
            # have consistent seps
            encoded_labels.extend([''] * len(columns))

        # write out the index label line
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
            d = b.to_native_types(slicer=slicer,
                                  na_rep=self.na_rep,
                                  float_format=self.float_format,
                                  decimal=self.decimal,
                                  date_format=self.date_format,
                                  quoting=self.quoting)

            for col_loc, col in zip(b.mgr_locs, d):
                # self.data is a preallocated list
                self.data[col_loc] = col

        ix = data_index.to_native_types(slicer=slicer,
                                        na_rep=self.na_rep,
                                        float_format=self.float_format,
                                        date_format=self.date_format,
                                        quoting=self.quoting)

        lib.write_csv_rows(self.data, ix, self.nlevels, self.cols, self.writer)

# from collections import namedtuple
# ExcelCell = namedtuple("ExcelCell",
#                        'row, col, val, style, mergestart, mergeend')


class ExcelCell(object):
    __fields__ = ('row', 'col', 'val', 'style', 'mergestart', 'mergeend')
    __slots__ = __fields__

    def __init__(self, row, col, val,
                 style=None, mergestart=None, mergeend=None):
        self.row = row
        self.col = col
        self.val = val
        self.style = style
        self.mergestart = mergestart
        self.mergeend = mergeend


header_style = {"font": {"bold": True},
                "borders": {"top": "thin",
                            "right": "thin",
                            "bottom": "thin",
                            "left": "thin"},
                "alignment": {"horizontal": "center", "vertical": "top"}}


class ExcelFormatter(object):

    """
    Class for formatting a DataFrame to a list of ExcelCells,

    Parameters
    ----------
    df : dataframe
    na_rep: na representation
    float_format : string, default None
            Format string for floating point numbers
    cols : sequence, optional
        Columns to write
    header : boolean or list of string, default True
        Write out column names. If a list of string is given it is
        assumed to be aliases for the column names
    index : boolean, default True
        output row names (index)
    index_label : string or sequence, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.
    merge_cells : boolean, default False
            Format MultiIndex and Hierarchical Rows as merged cells.
    inf_rep : string, default `'inf'`
        representation for np.inf values (which aren't representable in Excel)
        A `'-'` sign will be added in front of -inf.
    """

    def __init__(self, df, na_rep='', float_format=None, cols=None,
                 header=True, index=True, index_label=None, merge_cells=False,
                 inf_rep='inf'):
        self.df = df
        self.rowcounter = 0
        self.na_rep = na_rep
        self.columns = cols
        if cols is None:
            self.columns = df.columns
        self.float_format = float_format
        self.index = index
        self.index_label = index_label
        self.header = header
        self.merge_cells = merge_cells
        self.inf_rep = inf_rep

    def _format_value(self, val):
        if lib.checknull(val):
            val = self.na_rep
        elif com.is_float(val):
            if np.isposinf(val):
                val = self.inf_rep
            elif np.isneginf(val):
                val = '-%s' % self.inf_rep
            elif self.float_format is not None:
                val = float(self.float_format % val)
        return val

    def _format_header_mi(self):

        if self.columns.nlevels > 1:
            if not self.index:
                raise NotImplementedError("Writing to Excel with MultiIndex"
                                          " columns and no index ('index'=False) "
                                          "is not yet implemented.")

        has_aliases = isinstance(self.header, (tuple, list, np.ndarray, Index))
        if not(has_aliases or self.header):
            return

        columns = self.columns
        level_strs = columns.format(sparsify=True, adjoin=False, names=False)
        level_lengths = _get_level_lengths(level_strs)
        coloffset = 0
        lnum = 0

        if self.index and isinstance(self.df.index, MultiIndex):
            coloffset = len(self.df.index[0]) - 1

        if self.merge_cells:
            # Format multi-index as a merged cells.
            for lnum in range(len(level_lengths)):
                name = columns.names[lnum]
                yield ExcelCell(lnum, coloffset, name, header_style)

            for lnum, (spans, levels, labels) in enumerate(zip(level_lengths,
                                                               columns.levels,
                                                               columns.labels)
                                                           ):
                values = levels.take(labels)
                for i in spans:
                    if spans[i] > 1:
                        yield ExcelCell(lnum,
                                        coloffset + i + 1,
                                        values[i],
                                        header_style,
                                        lnum,
                                        coloffset + i + spans[i])
                    else:
                        yield ExcelCell(lnum,
                                        coloffset + i + 1,
                                        values[i],
                                        header_style)
        else:
            # Format in legacy format with dots to indicate levels.
            for i, values in enumerate(zip(*level_strs)):
                v = ".".join(map(com.pprint_thing, values))
                yield ExcelCell(lnum, coloffset + i + 1, v, header_style)

        self.rowcounter = lnum

    def _format_header_regular(self):
        has_aliases = isinstance(self.header, (tuple, list, np.ndarray, Index))
        if has_aliases or self.header:
            coloffset = 0

            if self.index:
                coloffset = 1
                if isinstance(self.df.index, MultiIndex):
                    coloffset = len(self.df.index[0])

            colnames = self.columns
            if has_aliases:
                if len(self.header) != len(self.columns):
                    raise ValueError(('Writing %d cols but got %d aliases'
                                      % (len(self.columns), len(self.header))))
                else:
                    colnames = self.header

            for colindex, colname in enumerate(colnames):
                yield ExcelCell(self.rowcounter, colindex + coloffset, colname,
                                header_style)

    def _format_header(self):
        if isinstance(self.columns, MultiIndex):
            gen = self._format_header_mi()
        else:
            gen = self._format_header_regular()

        gen2 = ()
        if self.df.index.names:
            row = [x if x is not None else ''
                   for x in self.df.index.names] + [''] * len(self.columns)
            if reduce(lambda x, y: x and y, map(lambda x: x != '', row)):
                gen2 = (ExcelCell(self.rowcounter, colindex, val, header_style)
                        for colindex, val in enumerate(row))
                self.rowcounter += 1
        return itertools.chain(gen, gen2)

    def _format_body(self):

        if isinstance(self.df.index, MultiIndex):
            return self._format_hierarchical_rows()
        else:
            return self._format_regular_rows()

    def _format_regular_rows(self):
        has_aliases = isinstance(self.header, (tuple, list, np.ndarray, Index))
        if has_aliases or self.header:
            self.rowcounter += 1

        coloffset = 0
        # output index and index_label?
        if self.index:
            # chek aliases
            # if list only take first as this is not a MultiIndex
            if self.index_label and isinstance(self.index_label,
                                               (list, tuple, np.ndarray, Index)):
                index_label = self.index_label[0]
            # if string good to go
            elif self.index_label and isinstance(self.index_label, str):
                index_label = self.index_label
            else:
                index_label = self.df.index.names[0]

            if isinstance(self.columns, MultiIndex):
                self.rowcounter += 1

            if index_label and self.header is not False:
                yield ExcelCell(self.rowcounter - 1,
                                0,
                                index_label,
                                header_style)

            # write index_values
            index_values = self.df.index
            if isinstance(self.df.index, PeriodIndex):
                index_values = self.df.index.to_timestamp()

            coloffset = 1
            for idx, idxval in enumerate(index_values):
                yield ExcelCell(self.rowcounter + idx, 0, idxval, header_style)

        # Get a frame that will account for any duplicates in the column names.
        col_mapped_frame = self.df.loc[:, self.columns]

        # Write the body of the frame data series by series.
        for colidx in range(len(self.columns)):
            series = col_mapped_frame.iloc[:, colidx]
            for i, val in enumerate(series):
                yield ExcelCell(self.rowcounter + i, colidx + coloffset, val)

    def _format_hierarchical_rows(self):
        has_aliases = isinstance(self.header, (tuple, list, np.ndarray, Index))
        if has_aliases or self.header:
            self.rowcounter += 1

        gcolidx = 0

        if self.index:
            index_labels = self.df.index.names
            # check for aliases
            if self.index_label and isinstance(self.index_label,
                                               (list, tuple, np.ndarray, Index)):
                index_labels = self.index_label

            # MultiIndex columns require an extra row
            # with index names (blank if None) for
            # unambigous round-trip
            if isinstance(self.columns, MultiIndex):
                self.rowcounter += 1

            # if index labels are not empty go ahead and dump
            if (any(x is not None for x in index_labels)
                    and self.header is not False):

                for cidx, name in enumerate(index_labels):
                    yield ExcelCell(self.rowcounter - 1,
                                    cidx,
                                    name,
                                    header_style)

            if self.merge_cells:
                # Format hierarchical rows as merged cells.
                level_strs = self.df.index.format(sparsify=True, adjoin=False,
                                                  names=False)
                level_lengths = _get_level_lengths(level_strs)

                for spans, levels, labels in zip(level_lengths,
                                                 self.df.index.levels,
                                                 self.df.index.labels):
                    values = levels.take(labels)
                    for i in spans:
                        if spans[i] > 1:
                            yield ExcelCell(self.rowcounter + i,
                                            gcolidx,
                                            values[i],
                                            header_style,
                                            self.rowcounter + i + spans[i] - 1,
                                            gcolidx)
                        else:
                            yield ExcelCell(self.rowcounter + i,
                                            gcolidx,
                                            values[i],
                                            header_style)
                    gcolidx += 1

            else:
                # Format hierarchical rows with non-merged values.
                for indexcolvals in zip(*self.df.index):
                    for idx, indexcolval in enumerate(indexcolvals):
                        yield ExcelCell(self.rowcounter + idx,
                                        gcolidx,
                                        indexcolval,
                                        header_style)
                    gcolidx += 1

        # Get a frame that will account for any duplicates in the column names.
        col_mapped_frame = self.df.loc[:, self.columns]

        # Write the body of the frame data series by series.
        for colidx in range(len(self.columns)):
            series = col_mapped_frame.iloc[:, colidx]
            for i, val in enumerate(series):
                yield ExcelCell(self.rowcounter + i, gcolidx + colidx, val)

    def get_formatted_cells(self):
        for cell in itertools.chain(self._format_header(),
                                    self._format_body()):
            cell.val = self._format_value(cell.val)
            yield cell

# ----------------------------------------------------------------------
# Array formatters


def format_array(values, formatter, float_format=None, na_rep='NaN',
                 digits=None, space=None, justify='right'):

    if com.is_categorical_dtype(values):
        fmt_klass = CategoricalArrayFormatter
    elif com.is_float_dtype(values.dtype):
        fmt_klass = FloatArrayFormatter
    elif com.is_period_arraylike(values):
        fmt_klass = PeriodArrayFormatter
    elif com.is_integer_dtype(values.dtype):
        fmt_klass = IntArrayFormatter
    elif com.is_datetimetz(values):
        fmt_klass = Datetime64TZFormatter
    elif com.is_datetime64_dtype(values.dtype):
        fmt_klass = Datetime64Formatter
    elif com.is_timedelta64_dtype(values.dtype):
        fmt_klass = Timedelta64Formatter
    else:
        fmt_klass = GenericArrayFormatter

    if space is None:
        space = get_option("display.column_space")

    if float_format is None:
        float_format = get_option("display.float_format")

    if digits is None:
        digits = get_option("display.precision")

    fmt_obj = fmt_klass(values, digits=digits, na_rep=na_rep,
                        float_format=float_format,
                        formatter=formatter, space=space,
                        justify=justify)

    return fmt_obj.get_result()


class GenericArrayFormatter(object):

    def __init__(self, values, digits=7, formatter=None, na_rep='NaN',
                 space=12, float_format=None, justify='right'):
        self.values = values
        self.digits = digits
        self.na_rep = na_rep
        self.space = space
        self.formatter = formatter
        self.float_format = float_format
        self.justify = justify

    def get_result(self):
        fmt_values = self._format_strings()
        return _make_fixed_width(fmt_values, self.justify)

    def _format_strings(self):
        if self.float_format is None:
            float_format = get_option("display.float_format")
            if float_format is None:
                fmt_str = '%% .%dg' % get_option("display.precision")
                float_format = lambda x: fmt_str % x
        else:
            float_format = self.float_format

        formatter = self.formatter if self.formatter is not None else \
            (lambda x: com.pprint_thing(x, escape_chars=('\t', '\r', '\n')))

        def _format(x):
            if self.na_rep is not None and lib.checknull(x):
                if x is None:
                    return 'None'
                elif x is pd.NaT:
                    return 'NaT'
                return self.na_rep
            elif isinstance(x, PandasObject):
                return '%s' % x
            else:
                # object dtype
                return '%s' % formatter(x)

        vals = self.values
        if isinstance(vals, Index):
            vals = vals._values

        is_float = lib.map_infer(vals, com.is_float) & notnull(vals)
        leading_space = is_float.any()

        fmt_values = []
        for i, v in enumerate(vals):
            if not is_float[i] and leading_space:
                fmt_values.append(' %s' % _format(v))
            elif is_float[i]:
                fmt_values.append(float_format(v))
            else:
                fmt_values.append(' %s' % _format(v))

        return fmt_values


class FloatArrayFormatter(GenericArrayFormatter):

    """

    """

    def __init__(self, *args, **kwargs):
        GenericArrayFormatter.__init__(self, *args, **kwargs)

        if self.float_format is not None and self.formatter is None:
            self.formatter = self.float_format

    def _format_with(self, fmt_str):
        def _val(x, threshold):
            if notnull(x):
                if (threshold is None or
                        abs(x) > get_option("display.chop_threshold")):
                    return fmt_str % x
                else:
                    if fmt_str.endswith("e"):  # engineering format
                        return "0"
                    else:
                        return fmt_str % 0
            else:

                return self.na_rep

        threshold = get_option("display.chop_threshold")
        fmt_values = [_val(x, threshold) for x in self.values]
        return _trim_zeros(fmt_values, self.na_rep)

    def _format_strings(self):
        if self.formatter is not None:
            fmt_values = [self.formatter(x) for x in self.values]
        else:
            fmt_str = '%% .%df' % self.digits
            fmt_values = self._format_with(fmt_str)

            if len(fmt_values) > 0:
                maxlen = max(len(x) for x in fmt_values)
            else:
                maxlen = 0

            too_long = maxlen > self.digits + 6

            abs_vals = np.abs(self.values)

            # this is pretty arbitrary for now
            has_large_values = (abs_vals > 1e8).any()
            has_small_values = ((abs_vals < 10 ** (-self.digits)) &
                                (abs_vals > 0)).any()

            if too_long and has_large_values:
                fmt_str = '%% .%de' % self.digits
                fmt_values = self._format_with(fmt_str)
            elif has_small_values:
                fmt_str = '%% .%de' % self.digits
                fmt_values = self._format_with(fmt_str)

        return fmt_values


class IntArrayFormatter(GenericArrayFormatter):

    def _format_strings(self):
        formatter = self.formatter or (lambda x: '% d' % x)
        fmt_values = [formatter(x) for x in self.values]
        return fmt_values


class Datetime64Formatter(GenericArrayFormatter):
    def __init__(self, values, nat_rep='NaT', date_format=None, **kwargs):
        super(Datetime64Formatter, self).__init__(values, **kwargs)
        self.nat_rep = nat_rep
        self.date_format = date_format

    def _format_strings(self):
        """ we by definition have DO NOT have a TZ """

        values = self.values
        if not isinstance(values, DatetimeIndex):
            values = DatetimeIndex(values)

        fmt_values = format_array_from_datetime(values.asi8.ravel(),
                                                format=_get_format_datetime64_from_values(values, self.date_format),
                                                na_rep=self.nat_rep).reshape(values.shape)
        return fmt_values.tolist()


class PeriodArrayFormatter(IntArrayFormatter):

    def _format_strings(self):
        values = PeriodIndex(self.values).to_native_types()
        formatter = self.formatter or (lambda x: '%s' % x)
        fmt_values = [formatter(x) for x in values]
        return fmt_values


class CategoricalArrayFormatter(GenericArrayFormatter):

    def __init__(self, values, *args, **kwargs):
        GenericArrayFormatter.__init__(self, values, *args, **kwargs)

    def _format_strings(self):
        fmt_values = format_array(self.values.get_values(), self.formatter,
                                  float_format=self.float_format,
                                  na_rep=self.na_rep, digits=self.digits,
                                  space=self.space, justify=self.justify)
        return fmt_values


def _is_dates_only(values):
    # return a boolean if we are only dates (and don't have a timezone)
    values = DatetimeIndex(values)
    if values.tz is not None:
        return False

    values_int = values.asi8
    consider_values = values_int != iNaT
    one_day_nanos = (86400 * 1e9)
    even_days = np.logical_and(consider_values, values_int % one_day_nanos != 0).sum() == 0
    if even_days:
        return True
    return False

def _format_datetime64(x, tz=None, nat_rep='NaT'):
    if x is None or lib.checknull(x):
        return nat_rep

    if tz is not None or not isinstance(x, Timestamp):
        x = Timestamp(x, tz=tz)

    return str(x)


def _format_datetime64_dateonly(x, nat_rep='NaT', date_format=None):
    if x is None or lib.checknull(x):
        return nat_rep

    if not isinstance(x, Timestamp):
        x = Timestamp(x)

    if date_format:
        return x.strftime(date_format)
    else:
        return x._date_repr

def _get_format_datetime64(is_dates_only, nat_rep='NaT', date_format=None):

    if is_dates_only:
        return lambda x, tz=None: _format_datetime64_dateonly(x,
                                                              nat_rep=nat_rep,
                                                              date_format=date_format)
    else:
        return lambda x, tz=None: _format_datetime64(x, tz=tz, nat_rep=nat_rep)


def _get_format_datetime64_from_values(values, date_format):
    """ given values and a date_format, return a string format """
    is_dates_only = _is_dates_only(values)
    if is_dates_only:
        return date_format or "%Y-%m-%d"
    return date_format


class Datetime64TZFormatter(Datetime64Formatter):

    def _format_strings(self):
        """ we by definition have a TZ """

        values = self.values.asobject
        is_dates_only = _is_dates_only(values)
        formatter = (self.formatter or _get_format_datetime64(is_dates_only, date_format=self.date_format))
        fmt_values = [ formatter(x) for x in values ]

        return fmt_values

class Timedelta64Formatter(GenericArrayFormatter):

    def __init__(self, values, nat_rep='NaT', box=False, **kwargs):
        super(Timedelta64Formatter, self).__init__(values, **kwargs)
        self.nat_rep = nat_rep
        self.box = box

    def _format_strings(self):
        formatter = self.formatter or _get_format_timedelta64(self.values, nat_rep=self.nat_rep,
                                                              box=self.box)
        fmt_values = np.array([formatter(x) for x in self.values])
        return fmt_values


def _get_format_timedelta64(values, nat_rep='NaT', box=False):
    """
    Return a formatter function for a range of timedeltas.
    These will all have the same format argument

    If box, then show the return in quotes
    """

    values_int = values.astype(np.int64)

    consider_values = values_int != iNaT

    one_day_nanos = (86400 * 1e9)
    even_days = np.logical_and(consider_values, values_int % one_day_nanos != 0).sum() == 0
    all_sub_day = np.logical_and(consider_values, np.abs(values_int) >= one_day_nanos).sum() == 0

    if even_days:
        format = 'even_day'
    elif all_sub_day:
        format = 'sub_day'
    else:
        format = 'long'

    def _formatter(x):
        if x is None or lib.checknull(x):
            return nat_rep

        if not isinstance(x, Timedelta):
            x = Timedelta(x)
        result = x._repr_base(format=format)
        if box:
            result = "'{0}'".format(result)
        return result

    return _formatter


def _make_fixed_width(strings, justify='right', minimum=None,
                      adj=None):

    if len(strings) == 0 or justify == 'all':
        return strings

    if adj is None:
        adj = _get_adjustment()

    max_len = np.max([adj.len(x) for x in strings])

    if minimum is not None:
        max_len = max(minimum, max_len)

    conf_max = get_option("display.max_colwidth")
    if conf_max is not None and max_len > conf_max:
        max_len = conf_max

    def just(x):
        if conf_max is not None:
            if (conf_max > 3) & (adj.len(x) > max_len):
                x = x[:max_len - 3] + '...'
        return x

    strings = [just(x) for x in strings]
    result = adj.justify(strings, max_len, mode=justify)
    return result


def _trim_zeros(str_floats, na_rep='NaN'):
    """
    Trims zeros and decimal points.
    """
    trimmed = str_floats

    def _cond(values):
        non_na = [x for x in values if x != na_rep]
        return (len(non_na) > 0 and all([x.endswith('0') for x in non_na]) and
                not(any([('e' in x) or ('E' in x) for x in non_na])))

    while _cond(trimmed):
        trimmed = [x[:-1] if x != na_rep else x for x in trimmed]

    # trim decimal points
    return [x[:-1] if x.endswith('.') and x != na_rep else x for x in trimmed]


def single_column_table(column, align=None, style=None):
    table = '<table'
    if align is not None:
        table += (' align="%s"' % align)
    if style is not None:
        table += (' style="%s"' % style)
    table += '><tbody>'
    for i in column:
        table += ('<tr><td>%s</td></tr>' % str(i))
    table += '</tbody></table>'
    return table


def single_row_table(row):  # pragma: no cover
    table = '<table><tbody><tr>'
    for i in row:
        table += ('<td>%s</td>' % str(i))
    table += '</tr></tbody></table>'
    return table


def _has_names(index):
    if isinstance(index, MultiIndex):
        return any([x is not None for x in index.names])
    else:
        return index.name is not None


# ------------------------------------------------------------------------------
# Global formatting options

_initial_defencoding = None


def detect_console_encoding():
    """
    Try to find the most capable encoding supported by the console.
    slighly modified from the way IPython handles the same issue.
    """
    import locale
    global _initial_defencoding

    encoding = None
    try:
        encoding = sys.stdout.encoding or sys.stdin.encoding
    except AttributeError:
        pass

    # try again for something better
    if not encoding or 'ascii' in encoding.lower():
        try:
            encoding = locale.getpreferredencoding()
        except Exception:
            pass

    # when all else fails. this will usually be "ascii"
    if not encoding or 'ascii' in encoding.lower():
        encoding = sys.getdefaultencoding()

    # GH3360, save the reported defencoding at import time
    # MPL backends may change it. Make available for debugging.
    if not _initial_defencoding:
        _initial_defencoding = sys.getdefaultencoding()

    return encoding


def get_console_size():
    """Return console size as tuple = (width, height).

    Returns (None,None) in non-interactive session.
    """
    display_width = get_option('display.width')
    # deprecated.
    display_height = get_option('display.height', silent=True)

    # Consider
    # interactive shell terminal, can detect term size
    # interactive non-shell terminal (ipnb/ipqtconsole), cannot detect term
    # size non-interactive script, should disregard term size

    # in addition
    # width,height have default values, but setting to 'None' signals
    # should use Auto-Detection, But only in interactive shell-terminal.
    # Simple. yeah.

    if com.in_interactive_session():
        if com.in_ipython_frontend():
            # sane defaults for interactive non-shell terminal
            # match default for width,height in config_init
            from pandas.core.config import get_default_val
            terminal_width = get_default_val('display.width')
            terminal_height = get_default_val('display.height')
        else:
            # pure terminal
            terminal_width, terminal_height = get_terminal_size()
    else:
        terminal_width, terminal_height = None, None

    # Note if the User sets width/Height to None (auto-detection)
    # and we're in a script (non-inter), this will return (None,None)
    # caller needs to deal.
    return (display_width or terminal_width, display_height or terminal_height)


class EngFormatter(object):

    """
    Formats float values according to engineering format.

    Based on matplotlib.ticker.EngFormatter
    """

    # The SI engineering prefixes
    ENG_PREFIXES = {
        -24: "y",
        -21: "z",
        -18: "a",
        -15: "f",
        -12: "p",
        -9: "n",
        -6: "u",
        -3: "m",
        0: "",
        3: "k",
        6: "M",
        9: "G",
        12: "T",
        15: "P",
        18: "E",
        21: "Z",
        24: "Y"
    }

    def __init__(self, accuracy=None, use_eng_prefix=False):
        self.accuracy = accuracy
        self.use_eng_prefix = use_eng_prefix

    def __call__(self, num):
        """ Formats a number in engineering notation, appending a letter
        representing the power of 1000 of the original number. Some examples:

        >>> format_eng(0)       # for self.accuracy = 0
        ' 0'

        >>> format_eng(1000000) # for self.accuracy = 1,
                                #     self.use_eng_prefix = True
        ' 1.0M'

        >>> format_eng("-1e-6") # for self.accuracy = 2
                                #     self.use_eng_prefix = False
        '-1.00E-06'

        @param num: the value to represent
        @type num: either a numeric value or a string that can be converted to
                   a numeric value (as per decimal.Decimal constructor)

        @return: engineering formatted string
        """
        import decimal
        import math
        dnum = decimal.Decimal(str(num))

        sign = 1

        if dnum < 0:  # pragma: no cover
            sign = -1
            dnum = -dnum

        if dnum != 0:
            pow10 = decimal.Decimal(int(math.floor(dnum.log10() / 3) * 3))
        else:
            pow10 = decimal.Decimal(0)

        pow10 = pow10.min(max(self.ENG_PREFIXES.keys()))
        pow10 = pow10.max(min(self.ENG_PREFIXES.keys()))
        int_pow10 = int(pow10)

        if self.use_eng_prefix:
            prefix = self.ENG_PREFIXES[int_pow10]
        else:
            if int_pow10 < 0:
                prefix = 'E-%02d' % (-int_pow10)
            else:
                prefix = 'E+%02d' % int_pow10

        mant = sign * dnum / (10 ** pow10)

        if self.accuracy is None:  # pragma: no cover
            format_str = u("% g%s")
        else:
            format_str = (u("%% .%if%%s") % self.accuracy)

        formatted = format_str % (mant, prefix)

        return formatted  # .strip()


def set_eng_float_format(accuracy=3, use_eng_prefix=False):
    """
    Alter default behavior on how float is formatted in DataFrame.
    Format float in engineering format. By accuracy, we mean the number of
    decimal digits after the floating point.

    See also EngFormatter.
    """

    set_option("display.float_format", EngFormatter(accuracy, use_eng_prefix))
    set_option("display.column_space", max(12, accuracy + 9))


def _put_lines(buf, lines):
    if any(isinstance(x, compat.text_type) for x in lines):
        lines = [compat.text_type(x) for x in lines]
    buf.write('\n'.join(lines))


def _binify(cols, line_width):
    adjoin_width = 1
    bins = []
    curr_width = 0
    i_last_column = len(cols) - 1
    for i, w in enumerate(cols):
        w_adjoined = w + adjoin_width
        curr_width += w_adjoined
        if i_last_column == i:
            wrap = curr_width + 1 > line_width and i > 0
        else:
            wrap = curr_width + 2 > line_width and i > 0
        if wrap:
            bins.append(i)
            curr_width = w_adjoined

    bins.append(len(cols))
    return bins

if __name__ == '__main__':
    arr = np.array([746.03, 0.00, 5620.00, 1592.36])
    # arr = np.array([11111111.1, 1.55])
    # arr = [314200.0034, 1.4125678]
    arr = np.array([327763.3119, 345040.9076, 364460.9915, 398226.8688,
                    383800.5172, 433442.9262, 539415.0568, 568590.4108,
                    599502.4276, 620921.8593, 620898.5294, 552427.1093,
                    555221.2193, 519639.7059, 388175.7, 379199.5854,
                    614898.25, 504833.3333, 560600., 941214.2857,
                    1134250., 1219550., 855736.85, 1042615.4286,
                    722621.3043, 698167.1818, 803750.])
    fmt = FloatArrayFormatter(arr, digits=7)
    print(fmt.get_result())
