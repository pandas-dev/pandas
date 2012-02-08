from itertools import izip

from StringIO import StringIO
from pandas.core.common import adjoin, isnull, notnull, _stringify
from pandas.core.index import MultiIndex, _ensure_index
from pandas.util import py3compat

import pandas.core.common as com
import pandas._tseries as lib

import numpy as np

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
        default None
    float_format : one-parameter function, optional
        formatter function to apply to columns' elements if they are floats
        default None
    sparsify : bool, optional
        Set to False for a DataFrame with a hierarchical index to print every
        multiindex key at each row, default True
    justify : {'left', 'right'}, default None
        Left or right-justify the column labels. If None uses the option from
        the configuration in pandas.core.common, 'left' out of the box
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
        self.buf = buf if buf is not None else StringIO()
        self.name = name
        self.na_rep = na_rep
        self.length = length
        self.header = header

        if float_format is None:
            float_format = print_config.float_format
        self.float_format = float_format

    def _get_footer(self):
        footer = ''
        if self.name:
            footer += ("Name: %s" % str(self.series.name)
                       if self.series.name else '')

        if self.length:
            if footer:
                footer += ', '
            footer += 'Length: %d' % len(self.series)
        return footer

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
            return ''

        fmt_index, have_header = self._get_formatted_index()
        fmt_values = self._get_formatted_values()

        maxlen = max(len(x) for x in fmt_index)
        pad_space = min(maxlen, 60)
        result = ['%s   %s' % (k.ljust(pad_space), v)
                  for (k, v) in izip(fmt_index[1:], fmt_values)]

        if self.header and have_header:
            result.insert(0, fmt_index[0])

        footer = self._get_footer()
        if footer:
            result.append(footer)

        return '\n'.join(result)


class DataFrameFormatter(object):
    """
    Render a DataFrame

    self.to_string() : console-friendly tabular output
    self.to_html()   : html table

    """

    __doc__ += docstring_to_string

    def __init__(self, frame, buf=None, columns=None, col_space=None,
                 header=True, index=True, na_rep='NaN', formatters=None,
                 justify=None, float_format=None, sparsify=True,
                 index_names=True, **kwds):
        self.frame = frame
        self.buf = buf if buf is not None else StringIO()
        self.show_index_names = index_names
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
        else:
            self.columns = frame.columns

    def to_string(self, force_unicode=False):
        """
        Render a DataFrame to a console-friendly tabular output.
        """
        frame = self.frame

        to_write = []

        if len(frame.columns) == 0 or len(frame.index) == 0:
            info_line = (u'Empty %s\nColumns: %s\nIndex: %s'
                         % (type(self.frame).__name__,
                            frame.columns, frame.index))
            to_write.append(info_line)
        else:
            # may include levels names also
            str_index = self._get_formatted_index()
            str_columns = self._get_formatted_column_labels()

            stringified = []

            for i, c in enumerate(self.columns):
                if self.header:
                    fmt_values = self._format_col(c)
                    cheader = str_columns[i]
                    max_len = max(max(len(x) for x in fmt_values),
                                  max(len(x) for x in cheader))
                    if self.justify == 'left':
                        cheader = [x.ljust(max_len) for x in cheader]
                    else:
                        cheader = [x.rjust(max_len) for x in cheader]
                    fmt_values = cheader + fmt_values
                    stringified.append(_make_fixed_width(fmt_values,
                                                         self.justify))
                else:
                    stringified = [_make_fixed_width(self._format_col(c),
                                                     self.justify)
                                   for c in self.columns]

            if self.index:
                to_write.append(adjoin(1, str_index, *stringified))
            else:
                to_write.append(adjoin(1, *stringified))

        if not py3compat.PY3:
            if force_unicode:
                to_write = [unicode(s) for s in to_write]
            else:
                # generally everything is plain strings, which has ascii
                # encoding.  problem is when there is a char with value over 127
                # - everything then gets converted to unicode.
                try:
                    for s in to_write:
                        str(s)
                except UnicodeError:
                    to_write = [unicode(s) for s in to_write]

        self.buf.writelines(to_write)

    def _format_col(self, col):
        formatter = self.formatters.get(col)
        return format_array(self.frame[col].values, formatter,
                            float_format=self.float_format,
                            na_rep=self.na_rep,
                            space=self.col_space)

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

        indent = 0
        indent_delta = 2
        frame = self.frame
        buf = self.buf

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

            # header row
            if self.header:
                write(buf, '<thead>', indent)
                row = []

                col_row = _column_header()
                indent += indent_delta
                write_tr(buf, col_row, indent, indent_delta, header=True)
                if self.has_index_names:
                    row = frame.index.names + [''] * len(frame.columns)
                    write_tr(buf, row, indent, indent_delta, header=True)

                write(buf, '</thead>', indent)

            write(buf, '<tbody>', indent)

            _bold_row = self.kwds.get('bold_rows', False)
            def _maybe_bold_row(x):
                temp = '<strong>%s</strong>'
                if _bold_row:
                    return ([temp % y for y in x] if isinstance(x, tuple)
                            else temp % x)
                else:
                    return x

            fmt_values = {}
            for col in frame.columns:
                fmt_values[col] = self._format_col(col)

            # write values
            for i in range(len(frame)):
                row = []
                if isinstance(frame.index, MultiIndex):
                    row.extend(_maybe_bold_row(frame.index[i]))
                else:
                    row.append(_maybe_bold_row(frame.index[i]))
                for col in frame.columns:
                    row.append(fmt_values[col][i])
                write_tr(buf, row, indent, indent_delta)
            indent -= indent_delta
            write(buf, '</tbody>', indent)
            indent -= indent_delta

        write(buf, '</table>', indent)

    def _get_formatted_column_labels(self):
        from pandas.core.index import _sparsify

        def is_numeric_dtype(dtype):
            return issubclass(dtype.type, np.number)

        if isinstance(self.columns, MultiIndex):
            fmt_columns = self.columns.format(sparsify=False, adjoin=False)
            fmt_columns = zip(*fmt_columns)
            dtypes = self.frame.dtypes.values
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))
            str_columns = zip(*[[u' %s' % y
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
            str_columns = [[u' %s' % x
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

#----------------------------------------------------------------------
# Array formatters


def format_array(values, formatter, float_format=None, na_rep='NaN',
                 digits=None, space=None, justify='right'):
    if com.is_float_dtype(values.dtype):
        fmt_klass = FloatArrayFormatter
    elif com.is_integer_dtype(values.dtype):
        fmt_klass = IntArrayFormatter
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
        if self.float_format is None:
            float_format = print_config.float_format
            if float_format is None:
                fmt_str = '%% .%dg' % print_config.precision
                float_format = lambda x: fmt_str % x
        else:
            float_format = self.float_format

        formatter = _stringify if self.formatter is None else self.formatter

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

        return _make_fixed_width(fmt_values, self.justify)

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
                maxlen =0

            too_long = maxlen > self.digits + 5

            # this is pretty arbitrary for now
            has_large_values = (np.abs(self.values) > 1e8).any()

            if too_long and has_large_values:
                fmt_str = '%% .%de' % (self.digits - 1)
                fmt_values = self._format_with(fmt_str)

        return _make_fixed_width(fmt_values, self.justify)


class IntArrayFormatter(GenericArrayFormatter):


    def get_result(self):
        fmt_values = ['% d' % x for x in self.values]
        return _make_fixed_width(fmt_values, self.justify)


def _make_fixed_width(strings, justify='right'):
    if len(strings) == 0:
        return strings

    max_len = max(len(x) for x in strings)
    if justify == 'left':
        justfunc = lambda self, x: self.ljust(x)
    else:
        justfunc = lambda self, x: self.rjust(x)

    def just(x):
        return justfunc(x[:max_len], max_len)

    return [just(x) for x in strings]

def _trim_zeros(str_floats, na_rep='NaN'):
    """
    Trims zeros and decimal points
    """
    # TODO: what if exponential?
    trimmed = str_floats

    def _cond(values):
        non_na = [x for x in values if x != na_rep]
        return len(non_na) > 0 and all([x.endswith('0') for x in non_na])

    while _cond(trimmed):
        trimmed = [x[:-1] if x != na_rep else x for x in trimmed]

    # trim decimal points
    return [x[:-1] if x.endswith('.') and x != na_rep else x for x in trimmed]


def single_column_table(column):
    table = '<table><tbody>'
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



#-------------------------------------------------------------------------------
# Global formatting options

def set_printoptions(precision=None, column_space=None, max_rows=None,
                     max_columns=None, colheader_justify='right'):
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
    """
    if precision is not None:
        print_config.precision = precision
    if column_space is not None:
        print_config.column_space = column_space
    if max_rows is not None:
        print_config.max_rows = max_rows
    if max_columns is not None:
        print_config.max_columns = max_columns
    if colheader_justify is not None:
        print_config.colheader_justify = colheader_justify

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
            pow10 = decimal.Decimal(int(math.floor(dnum.log10()/3)*3))
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

        mant = sign*dnum/(10**pow10)

        if self.accuracy is None:  # pragma: no cover
            format_str = u"% g%s"
        else:
            format_str = (u"%% .%if%%s" % self.accuracy )

        formatted = format_str % (mant, prefix)

        return formatted #.strip()

def set_eng_float_format(precision=None, accuracy=3, use_eng_prefix=False):
    """
    Alter default behavior on how float is formatted in DataFrame.
    Format float in engineering format. By accuracy, we mean the number of
    decimal digits after the floating point.

    See also EngFormatter.
    """
    if precision is not None: # pragma: no cover
        import warnings
        warnings.warn("'precision' parameter in set_eng_float_format is "
                      "being renamed to 'accuracy'" , FutureWarning)
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
        self.max_columns = 0
        self.colheader_justify = 'right'

    def reset(self):
        self.__init__()

print_config = _GlobalPrintConfig()


if __name__ == '__main__':
    arr = np.array([746.03, 0.00, 5620.00, 1592.36])
    # arr = np.array([11111111.1, 1.55])
    # arr = [314200.0034, 1.4125678]
    arr = np.array([  327763.3119,   345040.9076,   364460.9915,   398226.8688,
                      383800.5172,   433442.9262,   539415.0568,   568590.4108,
                      599502.4276,   620921.8593,   620898.5294,   552427.1093,
                      555221.2193,   519639.7059,   388175.7   ,   379199.5854,
                      614898.25  ,   504833.3333,   560600.    ,   941214.2857,
                      1134250.    ,  1219550.    ,   855736.85  ,  1042615.4286,
                      722621.3043,   698167.1818,   803750.    ])
    fmt = FloatArrayFormatter(arr, digits=7)
    print fmt.get_result()
