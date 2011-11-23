from StringIO import StringIO
from pandas.core.common import adjoin, _pfixed
from pandas.core.index import MultiIndex, _ensure_index


class DataFrameFormatter(object):
    """
    Render a DataFrame

    self.to_string() : console-friendly tabular output
    self.to_html() : html table
    """
    def __init__(self, frame, buf=None, columns=None, col_space=None,
                 na_rep='NaN', formatters=None, float_format=None,
                 sparsify=True, index_names=True):

        self.frame = frame
        self.buf = buf if buf is not None else StringIO()
        self.show_index_names = index_names
        self.sparsify = sparsify
        self.float_format = float_format
        self.formatters = formatters
        self.na_rep = na_rep
        self.col_space = col_space

        if columns is not None:
            self.columns = _ensure_index(columns)
        else:
            self.columns = frame.columns

    def to_string(self):
        """
        Render a DataFrame to a console-friendly tabular output.
        """
        frame = self.frame
        format_col = self._get_column_formatter()

        to_write = []

        if len(frame.columns) == 0 or len(frame.index) == 0:
            info_line = 'Empty %s\nColumns: %s\nIndex: %s'
            to_write.append(info_line % (type(self.frame).__name__,
                                         repr(frame.columns),
                                         repr(frame.index)))
        else:
            # may include levels names also
            str_index = self._get_formatted_index()
            str_columns = self._get_formatted_column_labels()

            stringified = [str_columns[i] + format_col(c)
                           for i, c in enumerate(self.columns)]

            to_write.append(adjoin(1, str_index, *stringified))

        for s in to_write:
            if isinstance(s, unicode):
                to_write = [unicode(s) for s in to_write]
                break

        self.buf.writelines(to_write)

    def to_html(self):
        """
        Render a DataFrame to a html table.
        """
        def write(buf, s, indent=0):
            buf.write(unicode((' ' * indent) + str(s) + '\n'))

        def write_th(buf, s, indent=0):
            write(buf, '<th>%s</th>' % str(s), indent)

        def write_td(buf, s, indent=0):
            write(buf, '<td>%s</td>' % str(s), indent)

        def write_tr(buf, l, indent=0, indent_delta=4, header=False):
            write(buf, '<tr>', indent)
            indent += indent_delta
            if header:
                for s in l:
                    write_th(buf, s, indent)
            else:
                for s in l:
                    write_td(buf, s, indent)
            indent -= indent_delta
            write(buf, '</tr>', indent)

        def single_column_table(column):
            table = '<table><tbody>'
            for i in column:
                table += ('<tr><td>%s</td></tr>' % str(i))
            table += '</tbody></table>'
            return table

        def single_row_table(row):
            table = '<table><tbody><tr>'
            for i in row:
                table += ('<td>%s</td>' % str(i))
            table += '</tr></tbody></table>'
            return table

        indent = 0
        indent_delta = 2
        frame = self.frame
        buf = self.buf
        format_col = self._get_column_formatter()

        write(buf, '<table border="1">', indent)

        def _column_header():
            row = [''] * (frame.index.nlevels - 1)

            if isinstance(frame.columns, MultiIndex):
                if self.has_column_names:
                    row.append(single_column_table(frame.columns.names))
                row.extend([single_column_table(c) for c in frame.columns])
            else:
                row.append(frame.columns.name or '')
                row.extend(frame.columns)
            return row

        if len(frame.columns) == 0 or len(frame.index) == 0:
            write(buf, '<tbody>', indent  + indent_delta)
            write_tr(buf,
                     [repr(frame.index),
                      'Empty %s' % type(self.frame).__name__],
                     indent + (2 * indent_delta),
                     indent_delta)
            write(buf, '</tbody>', indent  + indent_delta)
        else:
            indent += indent_delta
            write(buf, '<thead>', indent)
            row = []

            # header row
            col_row = _column_header()
            indent += indent_delta
            write_tr(buf, col_row, indent, indent_delta, header=True)
            if self.has_index_names:
                row = frame.index.names + [''] * len(frame.columns)
                write_tr(buf, row, indent, indent_delta, header=True)
            write(buf, '</thead>', indent)
            write(buf, '<tbody>', indent)

            # write values
            for i in range(len(frame)):
                row = []
                try:
                    row.extend(frame.index[i])
                except TypeError:
                    row.append(frame.index[i])
                for column in frame.columns:
                    row.append(format_col(column, i))
                write_tr(buf, row, indent, indent_delta)
            indent -= indent_delta
            write(buf, '</body>', indent)
            indent -= indent_delta

        write(buf, '</table>', indent)

    def _get_column_formatter(self):
        from pandas.core.common import _format

        col_space = self.col_space

        if col_space is None:
            def _myformat(v):
                return _format(v, na_rep=self.na_rep,
                               float_format=self.float_format)
        else:
            def _myformat(v):
                return _pfixed(v, col_space, na_rep=self.na_rep,
                               float_format=self.float_format)

        formatters = {} if self.formatters is None else self.formatters

        def _format_col(col, i=None):
            formatter = formatters.get(col, _myformat)
            if i == None:
                return [formatter(x) for x in self.frame[col]]
            else:
                return formatter(self.frame[col][i])

        return _format_col

    def _get_formatted_column_labels(self):
        from pandas.core.index import _sparsify

        if isinstance(self.columns, MultiIndex):
            fmt_columns = self.columns.format(sparsify=False, adjoin=False)
            str_columns = zip(*[[' %s' % y for y in x]
                                for x in zip(*fmt_columns)])
            if self.sparsify:
                str_columns = _sparsify(str_columns)

            str_columns = [list(x) for x in zip(*str_columns)]
        else:
            str_columns = [[' %s' % x] for x in self.columns.format()]

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
        index = self.frame.index
        columns = self.frame.columns

        show_index_names = self.show_index_names and self.has_index_names
        show_col_names = self.show_index_names and self.has_column_names

        if isinstance(index, MultiIndex):
            fmt_index = index.format(sparsify=self.sparsify, adjoin=False,
                                     names=show_index_names)
        else:
            fmt_index = [index.format(name=show_index_names)]

        adjoined = adjoin(1, *fmt_index).split('\n')

        # empty space for columns
        if show_col_names:
            col_header = ['  %s' % x for x in self._get_column_name_list()]
        else:
            col_header = [''] * columns.nlevels

        return col_header + adjoined

    def _get_column_name_list(self):
        names = []
        columns = self.frame.columns
        if isinstance(columns, MultiIndex):
            names.extend('' if name is None else name
                         for name in columns.names)
        else:
            names.append('' if columns.name is None else columns.name)
        return names

def _has_names(index):
    if isinstance(index, MultiIndex):
        return any([x is not None for x in index.names])
    else:
        return index.name is not None
