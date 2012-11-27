# pylint: disable=W0141

from itertools import izip
import sys

try:
    from StringIO import StringIO
except:
    from io import StringIO

from pandas.core.common import adjoin, isnull, notnull
from pandas.core.index import MultiIndex, _ensure_index
from pandas.util import py3compat

import pandas.core.common as com
import pandas.lib as lib

import numpy as np

import itertools


docstring_to_string = """
     Parameters
     ----------
     frame : DataFrame
         object to render
    buf : StringIO-like, optional
        buffer to write to
    columns : sequence, optional
        the subset of columns to write; default None writes all columns
    col_space : int, optional
        the width of each columns
    header : bool, optional
        whether to print column labels, default True
    index : bool, optional
        whether to print index (row) labels, default True
    na_rep : string, optional
        string representation of NAN to use, default 'NaN'
    formatters : list or dict of one-parameter functions, optional
        formatter functions to apply to columns' elements by position or name,
        default None, if the result is a string , it must be a unicode string.
    float_format : one-parameter function, optional
        formatter function to apply to columns' elements if they are floats
        default None
    sparsify : bool, optional
        Set to False for a DataFrame with a hierarchical index to print every
        multiindex key at each row, default True
    justify : {'left', 'right'}, default None
        Left or right-justify the column labels. If None uses the option from
        the print configuration (controlled by set_printoptions), 'right' out
        of the box.
    index_names : bool, optional
        Prints the names of the indexes, default True
    force_unicode : bool, default False
        Always return a unicode result

    Returns
    -------
    formatted : string (or unicode, depending on data and options)"""


class SeriesFormatter(object):

    def __init__(self, series, buf=None, header=True, length=True,
                 na_rep='NaN', name=False, float_format=None):
        self.series = series
        self.buf = buf if buf is not None else StringIO(u"")
        self.name = name
        self.na_rep = na_rep
        self.length = length
        self.header = header

        if float_format is None:
            float_format = print_config.float_format
        self.float_format = float_format

    def _get_footer(self):
        footer = u''

        if self.name:
            if getattr(self.series.index, 'freq', None):
                footer += 'Freq: %s' % self.series.index.freqstr

            if footer and self.series.name is not None:
                footer += ', '

            series_name = com.pprint_thing(self.series.name)
            footer += ("Name: %s" % series_name) if self.series.name is not None else ""

        if self.length:
            if footer:
                footer += ', '
            footer += 'Length: %d' % len(self.series)

        return unicode(footer)

    def _get_formatted_index(self):
        index = self.series.index
        is_multi = isinstance(index, MultiIndex)
        if is_multi:
            have_header = any(name for name in index.names)
            fmt_index = index.format(names=True)
        else:
            have_header = index.name is not None
            fmt_index = index.format(name=True)
        return fmt_index, have_header

    def _get_formatted_values(self):
        return format_array(self.series.values, None,
                            float_format=self.float_format,
                            na_rep=self.na_rep)

    def to_string(self):
        series = self.series

        if len(series) == 0:
            return u''

        fmt_index, have_header = self._get_formatted_index()
        fmt_values = self._get_formatted_values()

        maxlen = max(len(x) for x in fmt_index)
        pad_space = min(maxlen, 60)

        result = ['%s   %s'] * len(fmt_values)
        for i, (k, v) in enumerate(izip(fmt_index[1:], fmt_values)):
            try:
                idx = k.ljust(pad_space + _encode_diff(k))
            except UnicodeEncodeError:
                idx = k.ljust(pad_space)
            result[i] = result[i] % (idx, v)

        if self.header and have_header:
            result.insert(0, fmt_index[0])

        footer = self._get_footer()
        if footer:
            result.append(footer)

        return unicode(u'\n'.join(result))

if py3compat.PY3:  # pragma: no cover
    _encode_diff = lambda x: 0

    _strlen = len
else:
    def _encode_diff(x):
        return len(x) - len(x.decode(print_config.encoding))

    def _strlen(x):
        try:
            return len(x.decode(print_config.encoding))
        except UnicodeError:
            return len(x)


class DataFrameFormatter(object):
    """
    Render a DataFrame

    self.to_string() : console-friendly tabular output
    self.to_html()   : html table
    self.to_latex()   : LaTeX tabular environment table

    """

    __doc__ = __doc__ if __doc__ else ''
    __doc__ += docstring_to_string

    def __init__(self, frame, buf=None, columns=None, col_space=None,
                 header=True, index=True, na_rep='NaN', formatters=None,
                 justify=None, float_format=None, sparsify=None,
                 index_names=True, **kwds):
        self.frame = frame
        self.buf = buf if buf is not None else StringIO()
        self.show_index_names = index_names

        if sparsify is None:
            sparsify = print_config.multi_sparse

        self.sparsify = sparsify

        self.float_format = float_format
        self.formatters = formatters if formatters is not None else {}
        self.na_rep = na_rep
        self.col_space = col_space
        self.header = header
        self.index = index

        if justify is None:
            self.justify = print_config.colheader_justify
        else:
            self.justify = justify

        self.kwds = kwds

        if columns is not None:
            self.columns = _ensure_index(columns)
            self.frame = self.frame[self.columns]
        else:
            self.columns = frame.columns

    def _to_str_columns(self, force_unicode=None):
        """
        Render a DataFrame to a list of columns (as lists of strings).
        """
        import warnings
        if force_unicode is not None:  # pragma: no cover
            warnings.warn("force_unicode is deprecated, it will have no effect",
                          FutureWarning)

        # may include levels names also
        str_index = self._get_formatted_index()
        str_columns = self._get_formatted_column_labels()

        stringified = []

        for i, c in enumerate(self.columns):
            if self.header:
                fmt_values = self._format_col(i)
                cheader = str_columns[i]

                max_colwidth = max(_strlen(x) for x in cheader)

                fmt_values = _make_fixed_width(fmt_values, self.justify,
                                               minimum=max_colwidth)

                max_len = max(max(_strlen(x) for x in fmt_values),
                              max_colwidth)
                if self.justify == 'left':
                    cheader = [x.ljust(max_len) for x in cheader]
                else:
                    cheader = [x.rjust(max_len) for x in cheader]

                stringified.append(cheader + fmt_values)
            else:
                stringified = [_make_fixed_width(self._format_col(i),
                                                 self.justify)
                               for i, c in enumerate(self.columns)]

        strcols = stringified
        if self.index:
            strcols.insert(0, str_index)

        return strcols

    def to_string(self, force_unicode=None):
        """
        Render a DataFrame to a console-friendly tabular output.
        """
        import warnings
        if force_unicode is not None:  # pragma: no cover
            warnings.warn("force_unicode is deprecated, it will have no effect",
                          FutureWarning)

        frame = self.frame

        if len(frame.columns) == 0 or len(frame.index) == 0:
            info_line = (u'Empty %s\nColumns: %s\nIndex: %s'
                         % (type(self.frame).__name__,
                            com.pprint_thing(frame.columns),
                            com.pprint_thing(frame.index)))
            text = info_line
        else:
            strcols = self._to_str_columns()
            text = adjoin(1, *strcols)

        self.buf.writelines(text)

    def to_latex(self, force_unicode=None, column_format=None):
        """
        Render a DataFrame to a LaTeX tabular environment output.
        """
        import warnings
        if force_unicode is not None:  # pragma: no cover
            warnings.warn("force_unicode is deprecated, it will have no effect",
                          FutureWarning)

        frame = self.frame

        if len(frame.columns) == 0 or len(frame.index) == 0:
            info_line = (u'Empty %s\nColumns: %s\nIndex: %s'
                         % (type(self.frame).__name__,
                            frame.columns, frame.index))
            strcols = [[info_line]]
        else:
            strcols = self._to_str_columns()

        if column_format is None:
            column_format = '|l|%s|' % '|'.join('c' for _ in strcols)
        elif not isinstance(column_format, basestring):
            raise AssertionError(('column_format must be str or unicode, not %s'
                                  % type(column_format)))

        self.buf.write('\\begin{tabular}{%s}\n' % column_format)
        self.buf.write('\\hline\n')

        nlevels = frame.index.nlevels
        for i, row in enumerate(izip(*strcols)):
            if i == nlevels:
                self.buf.write('\\hline\n')  # End of header
            crow = [(x.replace('_', '\\_')
                     .replace('%', '\\%')
                     .replace('&', '\\&') if x else '{}') for x in row]
            self.buf.write(' & '.join(crow))
            self.buf.write(' \\\\\n')

        self.buf.write('\\hline\n')
        self.buf.write('\\end{tabular}\n')

    def _format_col(self, i):
        col = self.columns[i]
        formatter = self.formatters.get(col)
        return format_array(self.frame.icol(i).values, formatter,
                            float_format=self.float_format,
                            na_rep=self.na_rep,
                            space=self.col_space)

    def to_html(self, classes=None):
        """
        Render a DataFrame to a html table.
        """
        html_renderer = HTMLFormatter(self, classes=classes)
        html_renderer.write_result(self.buf)

    def _get_formatted_column_labels(self):
        from pandas.core.index import _sparsify

        def is_numeric_dtype(dtype):
            return issubclass(dtype.type, np.number)

        if isinstance(self.columns, MultiIndex):
            fmt_columns = self.columns.format(sparsify=False, adjoin=False)
            fmt_columns = zip(*fmt_columns)
            dtypes = self.frame.dtypes.values
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))
            str_columns = zip(*[[' ' + y
                                if y not in self.formatters and need_leadsp[x]
                                else y for y in x]
                               for x in fmt_columns])
            if self.sparsify:
                str_columns = _sparsify(str_columns)

            str_columns = [list(x) for x in zip(*str_columns)]
        else:
            fmt_columns = self.columns.format()
            dtypes = self.frame.dtypes
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))
            str_columns = [[' ' + x
                            if col not in self.formatters and need_leadsp[x]
                            else x]
                           for col, x in zip(self.columns, fmt_columns)]

        if self.show_index_names and self.has_index_names:
            for x in str_columns:
                x.append('')

        return str_columns

    @property
    def has_index_names(self):
        return _has_names(self.frame.index)

    @property
    def has_column_names(self):
        return _has_names(self.frame.columns)

    def _get_formatted_index(self):
        # Note: this is only used by to_string(), not by to_html().
        index = self.frame.index
        columns = self.frame.columns

        show_index_names = self.show_index_names and self.has_index_names
        show_col_names = (self.show_index_names and self.has_column_names)

        if isinstance(index, MultiIndex):
            fmt_index = index.format(sparsify=self.sparsify, adjoin=False,
                                     names=show_index_names)
        else:
            fmt_index = [index.format(name=show_index_names)]

        adjoined = adjoin(1, *fmt_index).split('\n')

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


class HTMLFormatter(object):

    indent_delta = 2

    def __init__(self, formatter, classes=None):
        self.fmt = formatter
        self.classes = classes

        self.frame = self.fmt.frame
        self.columns = formatter.columns
        self.elements = []

        _bold_row = self.fmt.kwds.get('bold_rows', False)
        _temp = '<strong>%s</strong>'

        def _maybe_bold_row(x):
            if _bold_row:
                return ([_temp % y for y in x] if isinstance(x, tuple)
                        else _temp % x)
            else:
                return x
        self._maybe_bold_row = _maybe_bold_row

    def write(self, s, indent=0):
        self.elements.append(' ' * indent + com.pprint_thing(s))

    def write_th(self, s, indent=0, tags=None):
        return self._write_cell(s, kind='th', indent=indent, tags=tags)

    def write_td(self, s, indent=0, tags=None):
        return self._write_cell(s, kind='td', indent=indent, tags=tags)

    def _write_cell(self, s, kind='td', indent=0, tags=None):
        if tags is not None:
            start_tag = '<%s %s>' % (kind, tags)
        else:
            start_tag = '<%s>' % kind
        self.write('%s%s</%s>' % (start_tag, com.pprint_thing(s), kind), indent)

    def write_tr(self, line, indent=0, indent_delta=4, header=False,
                 align=None, tags=None):
        if tags is None:
            tags = {}

        if align is None:
            self.write('<tr>', indent)
        else:
            self.write('<tr style="text-align: %s;">' % align, indent)
        indent += indent_delta

        for i, s in enumerate(line):
            val_tag = tags.get(i, None)
            if header:
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

        self.write('<table border="1" class="%s">' % ' '.join(_classes),
                   indent)

        if len(frame.columns) == 0 or len(frame.index) == 0:
            self.write('<tbody>', indent + self.indent_delta)
            self.write_tr([repr(frame.index),
                           'Empty %s' % type(frame).__name__],
                          indent + (2 * self.indent_delta),
                          self.indent_delta)
            self.write('</tbody>', indent + self.indent_delta)
        else:
            indent += self.indent_delta
            indent = self._write_header(indent)
            indent = self._write_body(indent)

        self.write('</table>', indent)

        _put_lines(buf, self.elements)

    def _write_header(self, indent):
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

            levels = self.columns.format(sparsify=True, adjoin=False,
                                         names=False)
            level_lengths = _get_level_lengths(levels)

            row_levels = self.frame.index.nlevels

            for lnum, (records, values) in enumerate(zip(level_lengths, levels)):
                name = self.columns.names[lnum]
                row = [''] * (row_levels - 1) + ['' if name is None
                                                 else str(name)]

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

            self.write_tr(col_row, indent, self.indent_delta, header=True,
                    align=align)

        if self.fmt.has_index_names:
            row = [x if x is not None else ''
                   for x in self.frame.index.names] + [''] * len(self.columns)
            self.write_tr(row, indent, self.indent_delta, header=True)

        indent -= self.indent_delta
        self.write('</thead>', indent)

        return indent

    def _write_body(self, indent):
        self.write('<tbody>', indent)
        indent += self.indent_delta

        fmt_values = {}
        for i in range(len(self.columns)):
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
        ncols = len(self.columns)

        if '__index__' in self.fmt.formatters:
            f = self.fmt.formatters['__index__']
            index_values = self.frame.index.values.map(f)
        else:
            index_values = self.frame.index.format()

        for i in range(len(self.frame)):
            row = []
            row.append(self._maybe_bold_row(index_values[i]))
            row.extend(fmt_values[j][i] for j in range(ncols))
            self.write_tr(row, indent, self.indent_delta, tags=None)

    def _write_hierarchical_rows(self, fmt_values, indent):
        template = 'rowspan="%d" valign="top"'

        frame = self.frame
        ncols = len(self.columns)

        idx_values = frame.index.format(sparsify=False, adjoin=False,
                                        names=False)
        idx_values = zip(*idx_values)

        if self.fmt.sparsify:
            levels = frame.index.format(sparsify=True, adjoin=False,
                                        names=False)
            level_lengths = _get_level_lengths(levels)

            for i in range(len(frame)):
                row = []
                tags = {}

                j = 0
                for records, v in zip(level_lengths, idx_values[i]):
                    if i in records:
                        if records[i] > 1:
                            tags[j] = template % records[i]
                    else:
                        continue
                    j += 1
                    row.append(self._maybe_bold_row(v))

                row.extend(fmt_values[j][i] for j in range(ncols))
                self.write_tr(row, indent, self.indent_delta, tags=tags)
        else:
            for i in range(len(frame)):
                idx_values = zip(*frame.index.format(sparsify=False,
                                                     adjoin=False,
                                                     names=False))
                row = []
                row.extend(self._maybe_bold_row(x) for x in idx_values[i])
                row.extend(fmt_values[j][i] for j in range(ncols))
                self.write_tr(row, indent, self.indent_delta, tags=None)


def _get_level_lengths(levels):
    from itertools import groupby

    def _make_grouper():
        record = {'count': 0}

        def grouper(x):
            if x != '':
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


#from collections import namedtuple
# ExcelCell = namedtuple("ExcelCell",
#                        'row, col, val, style, mergestart, mergeend')

class ExcelCell:
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
              "alignment": {"horizontal": "center"}}


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
    index : boolean, default True
        output row names (index)
    """

    def __init__(self,
                 df,
                 na_rep='',
                 float_format=None,
                 cols=None
                 ):
        self.df = df
        self.rowcounter = 0
        self.na_rep = na_rep
        self.columns = cols
        if cols is None:
            self.columns = df.columns
        self.float_format = float_format

    def _format_value(self, val):
        if lib.checknull(val):
            val = self.na_rep
        if self.float_format is not None and com.is_float(val):
            val = self.float_format % val
        return val

    def _format_header_mi(self):
        levels = self.columns.format(sparsify=True, adjoin=False,
                                   names=False)
        level_lenghts = _get_level_lengths(levels)
        for lnum, (records, values) in enumerate(zip(level_lenghts,
                                                     levels)):
            name = self.columns.names[lnum]
            yield ExcelCell(lnum, 0, name, header_style)
            for i in records:
                if records[i] > 1:
                    yield ExcelCell(lnum, i + 1, values[i],
                            header_style, lnum, i + records[i])
                else:
                    yield ExcelCell(lnum, i + 1, values[i], header_style)

            self.rowcounter = lnum

    def _format_header_regular(self):
        for colindex, colname in enumerate(self.columns):
            yield ExcelCell(self.rowcounter, colindex, colname, header_style)

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

        if isinstance(self.df.columns, MultiIndex):
            return self._format_hierarchical_rows()
        else:
            return self._format_regular_rows()

    def _format_regular_rows(self):
        self.rowcounter += 1
        for colidx, colname in enumerate(self.columns):
            series = self.df[colname]
            for i, val in enumerate(series):
                yield ExcelCell(self.rowcounter + i, colidx, val)

    def _format_hierarchical_rows(self):
        self.rowcounter += 1
        for idx, idxval in enumerate(self.df.index):
            yield ExcelCell(self.rowcounter + idx, 0, idxval, header_style)

        for colidx, colname in enumerate(self.columns):
            series = self.df[colname]
            for i, val in enumerate(series):
                yield ExcelCell(self.rowcounter + i, colidx + 1, val)

    def get_formatted_cells(self):
        for cell in itertools.chain(self._format_header(),
                                    self._format_body()):
            cell.val = self._format_value(cell.val)
            yield cell

#----------------------------------------------------------------------
# Array formatters


def format_array(values, formatter, float_format=None, na_rep='NaN',
                 digits=None, space=None, justify='right'):
    if com.is_float_dtype(values.dtype):
        fmt_klass = FloatArrayFormatter
    elif com.is_integer_dtype(values.dtype):
        fmt_klass = IntArrayFormatter
    elif com.is_datetime64_dtype(values.dtype):
        fmt_klass = Datetime64Formatter
    else:
        fmt_klass = GenericArrayFormatter

    if space is None:
        space = print_config.column_space

    if float_format is None:
        float_format = print_config.float_format

    if digits is None:
        digits = print_config.precision

    fmt_obj = fmt_klass(values, digits, na_rep=na_rep,
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
            float_format = print_config.float_format
            if float_format is None:
                fmt_str = '%% .%dg' % print_config.precision
                float_format = lambda x: fmt_str % x
        else:
            float_format = self.float_format

        formatter = com.pprint_thing if self.formatter is None else self.formatter

        def _format(x):
            if self.na_rep is not None and lib.checknull(x):
                if x is None:
                    return 'None'
                return self.na_rep
            else:
                # object dtype
                return '%s' % formatter(x)

        vals = self.values

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
        fmt_values = [fmt_str % x if notnull(x) else self.na_rep
                      for x in self.values]
        return _trim_zeros(fmt_values, self.na_rep)

    def get_result(self):
        if self.formatter is not None:
            fmt_values = [self.formatter(x) for x in self.values]
        else:
            fmt_str = '%% .%df' % (self.digits - 1)
            fmt_values = self._format_with(fmt_str)

            if len(fmt_values) > 0:
                maxlen = max(len(x) for x in fmt_values)
            else:
                maxlen = 0

            too_long = maxlen > self.digits + 5

            abs_vals = np.abs(self.values)

            # this is pretty arbitrary for now
            has_large_values = (abs_vals > 1e8).any()
            has_small_values = ((abs_vals < 10 ** (-self.digits)) &
                                (abs_vals > 0)).any()

            if too_long and has_large_values:
                fmt_str = '%% .%de' % (self.digits - 1)
                fmt_values = self._format_with(fmt_str)
            elif has_small_values:
                fmt_str = '%% .%de' % (self.digits - 1)
                fmt_values = self._format_with(fmt_str)

        return _make_fixed_width(fmt_values, self.justify)


class IntArrayFormatter(GenericArrayFormatter):

    def get_result(self):
        if self.formatter:
            formatter = self.formatter
        else:
            formatter = lambda x: '% d' % x

        fmt_values = [formatter(x) for x in self.values]

        return _make_fixed_width(fmt_values, self.justify)


class Datetime64Formatter(GenericArrayFormatter):

    def get_result(self):
        if self.formatter:
            formatter = self.formatter
        else:
            formatter = _format_datetime64

        fmt_values = [formatter(x) for x in self.values]
        return _make_fixed_width(fmt_values, self.justify)


def _format_datetime64(x, tz=None):
    if isnull(x):
        return 'NaT'

    stamp = lib.Timestamp(x, tz=tz)
    return stamp._repr_base


def _make_fixed_width(strings, justify='right', minimum=None):
    if len(strings) == 0:
        return strings

    max_len = max(_strlen(x) for x in strings)

    if minimum is not None:
        max_len = max(minimum, max_len)

    conf_max = print_config.max_colwidth
    if conf_max is not None and max_len > conf_max:
        max_len = conf_max

    if justify == 'left':
        justfunc = lambda self, x: self.ljust(x)
    else:
        justfunc = lambda self, x: self.rjust(x)

    def just(x):
        try:
            eff_len = max_len + _encode_diff(x)
        except UnicodeError:
            eff_len = max_len

        if conf_max is not None:
            if (conf_max > 3) & (_strlen(x) > max_len):
                x = x[:eff_len - 3] + '...'

        return justfunc(x, eff_len)

    return [just(x) for x in strings]


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


#------------------------------------------------------------------------------
# Global formatting options


def set_printoptions(precision=None, column_space=None, max_rows=None,
                     max_columns=None, colheader_justify=None,
                     max_colwidth=None, notebook_repr_html=None,
                     date_dayfirst=None, date_yearfirst=None,
                     multi_sparse=None, encoding=None):
    """
    Alter default behavior of DataFrame.toString

    precision : int
        Floating point output precision (number of significant digits). This is
        only a suggestion
    column_space : int
        Default space for DataFrame columns, defaults to 12
    max_rows : int
    max_columns : int
        max_rows and max_columns are used in __repr__() methods to decide if
        to_string() or info() is used to render an object to a string.
        Either one, or both can be set to 0 (experimental). Pandas will figure
        out how big the terminal is and will not display more rows or/and
        columns that can fit on it.
    colheader_justify
    notebook_repr_html : boolean
        When True (default), IPython notebook will use html representation for
        pandas objects (if it is available).
    date_dayfirst : boolean
        When True, prints and parses dates with the day first, eg 20/01/2005
    date_yearfirst : boolean
        When True, prints and parses dates with the year first, eg 2005/01/20
    multi_sparse : boolean
        Default True, "sparsify" MultiIndex display (don't display repeated
        elements in outer levels within groups)
    """
    if precision is not None:
        print_config.precision = precision
    if column_space is not None:
        print_config.column_space = column_space
    if max_rows is not None:
        print_config.max_rows = max_rows
    if max_colwidth is not None:
        print_config.max_colwidth = max_colwidth
    if max_columns is not None:
        print_config.max_columns = max_columns
    if colheader_justify is not None:
        print_config.colheader_justify = colheader_justify
    if notebook_repr_html is not None:
        print_config.notebook_repr_html = notebook_repr_html
    if date_dayfirst is not None:
        print_config.date_dayfirst = date_dayfirst
    if date_yearfirst is not None:
        print_config.date_yearfirst = date_yearfirst
    if multi_sparse is not None:
        print_config.multi_sparse = multi_sparse
    if encoding is not None:
        print_config.encoding = encoding


def reset_printoptions():
    print_config.reset()


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
            format_str = u"% g%s"
        else:
            format_str = (u"%% .%if%%s" % self.accuracy)

        formatted = format_str % (mant, prefix)

        return formatted  #.strip()


def set_eng_float_format(precision=None, accuracy=3, use_eng_prefix=False):
    """
    Alter default behavior on how float is formatted in DataFrame.
    Format float in engineering format. By accuracy, we mean the number of
    decimal digits after the floating point.

    See also EngFormatter.
    """
    if precision is not None:  # pragma: no cover
        import warnings
        warnings.warn("'precision' parameter in set_eng_float_format is "
                      "being renamed to 'accuracy'", FutureWarning)
        accuracy = precision

    print_config.float_format = EngFormatter(accuracy, use_eng_prefix)
    print_config.column_space = max(12, accuracy + 9)


class _GlobalPrintConfig(object):
    """
    Holds the console formatting settings for DataFrame and friends
    """

    def __init__(self):
        self.precision = self.digits = 7
        self.float_format = None
        self.column_space = 12
        self.max_rows = 200
        self.max_colwidth = 50
        self.max_columns = 0
        self.colheader_justify = 'right'
        self.notebook_repr_html = True
        self.date_dayfirst = False
        self.date_yearfirst = False
        self.pprint_nest_depth = 3
        self.multi_sparse = True
        self.encoding = self.detect_encoding()

    def detect_encoding(self):
        """
        Try to find the most capable encoding supported by the console.
        slighly modified from the way IPython handles the same issue.
        """
        import locale

        encoding = None
        try:
            encoding = sys.stdin.encoding
        except AttributeError:
            pass

        if not encoding or encoding == 'ascii':  # try again for better
            try:
                encoding = locale.getpreferredencoding()
            except Exception:
                pass

        if not encoding:  # when all else fails. this will usually be "ascii"
            encoding = sys.getdefaultencoding()

        return encoding

    def reset(self):
        self.__init__()

print_config = _GlobalPrintConfig()


def _put_lines(buf, lines):
    if any(isinstance(x, unicode) for x in lines):
        lines = [unicode(x) for x in lines]
    buf.write('\n'.join(lines))


if __name__ == '__main__':
    arr = np.array([746.03, 0.00, 5620.00, 1592.36])
    # arr = np.array([11111111.1, 1.55])
    # arr = [314200.0034, 1.4125678]
    arr = np.array([  327763.3119,   345040.9076,   364460.9915,   398226.8688,
                      383800.5172,   433442.9262,   539415.0568,   568590.4108,
                      599502.4276,   620921.8593,   620898.5294,   552427.1093,
                      555221.2193,   519639.7059,   388175.7   ,   379199.5854,
                      614898.25  ,   504833.3333,   560600.    ,   941214.2857,
                      1134250.    ,  1219550.    ,   855736.85 ,  1042615.4286,
                      722621.3043,   698167.1818,   803750.    ])
    fmt = FloatArrayFormatter(arr, digits=7)
    print fmt.get_result()
