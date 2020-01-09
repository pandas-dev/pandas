"""
Internal module for formatting output data in csv, html,
and latex files. This module also applies to display formatting.
"""

from contextlib import contextmanager
from datetime import tzinfo
import decimal
from functools import partial
from io import StringIO
import math
import re
from shutil import get_terminal_size
from typing import (
    IO,
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
    cast,
)
from unicodedata import east_asian_width

import numpy as np

from pandas._config.config import get_option, set_option

from pandas._libs import lib
from pandas._libs.missing import NA
from pandas._libs.tslib import format_array_from_datetime
from pandas._libs.tslibs import NaT, Timedelta, Timestamp, iNaT
from pandas._libs.tslibs.nattype import NaTType
from pandas._typing import FilePathOrBuffer
from pandas.errors import AbstractMethodError

from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_complex_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_extension_array_dtype,
    is_float,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_list_like,
    is_numeric_dtype,
    is_scalar,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.generic import (
    ABCIndexClass,
    ABCMultiIndex,
    ABCSeries,
    ABCSparseArray,
)
from pandas.core.dtypes.missing import isna, notna

from pandas.core.arrays.datetimes import DatetimeArray
from pandas.core.arrays.timedeltas import TimedeltaArray
from pandas.core.base import PandasObject
import pandas.core.common as com
from pandas.core.indexes.api import Index, ensure_index
from pandas.core.indexes.datetimes import DatetimeIndex
from pandas.core.indexes.timedeltas import TimedeltaIndex

from pandas.io.common import stringify_path
from pandas.io.formats.printing import adjoin, justify, pprint_thing

if TYPE_CHECKING:
    from pandas import Series, DataFrame, Categorical

formatters_type = Union[
    List[Callable], Tuple[Callable, ...], Mapping[Union[str, int], Callable]
]
float_format_type = Union[str, Callable, "EngFormatter"]

common_docstring = """
        Parameters
        ----------
        buf : str, Path or StringIO-like, optional, default None
            Buffer to write to. If None, the output is returned as a string.
        columns : sequence, optional, default None
            The subset of columns to write. Writes all columns by default.
        col_space : %(col_space_type)s, optional
            %(col_space)s.
        header : %(header_type)s, optional
            %(header)s.
        index : bool, optional, default True
            Whether to print index (row) labels.
        na_rep : str, optional, default 'NaN'
            String representation of NAN to use.
        formatters : list, tuple or dict of one-param. functions, optional
            Formatter functions to apply to columns' elements by position or
            name.
            The result of each function must be a unicode string.
            List/tuple must be of length equal to the number of columns.
        float_format : one-parameter function, optional, default None
            Formatter function to apply to columns' elements if they are
            floats. The result of this function must be a unicode string.
        sparsify : bool, optional, default True
            Set to False for a DataFrame with a hierarchical index to print
            every multiindex key at each row.
        index_names : bool, optional, default True
            Prints the names of the indexes.
        justify : str, default None
            How to justify the column labels. If None uses the option from
            the print configuration (controlled by set_option), 'right' out
            of the box. Valid values are

            * left
            * right
            * center
            * justify
            * justify-all
            * start
            * end
            * inherit
            * match-parent
            * initial
            * unset.
        max_rows : int, optional
            Maximum number of rows to display in the console.
        min_rows : int, optional
            The number of rows to display in the console in a truncated repr
            (when number of rows is above `max_rows`).
        max_cols : int, optional
            Maximum number of columns to display in the console.
        show_dimensions : bool, default False
            Display DataFrame dimensions (number of rows by number of columns).
        decimal : str, default '.'
            Character recognized as decimal separator, e.g. ',' in Europe.
    """

_VALID_JUSTIFY_PARAMETERS = (
    "left",
    "right",
    "center",
    "justify",
    "justify-all",
    "start",
    "end",
    "inherit",
    "match-parent",
    "initial",
    "unset",
)

return_docstring = """
        Returns
        -------
        str or None
            If buf is None, returns the result as a string. Otherwise returns
            None.
    """


class CategoricalFormatter:
    def __init__(
        self,
        categorical: "Categorical",
        buf: Optional[IO[str]] = None,
        length: bool = True,
        na_rep: str = "NaN",
        footer: bool = True,
    ):
        self.categorical = categorical
        self.buf = buf if buf is not None else StringIO("")
        self.na_rep = na_rep
        self.length = length
        self.footer = footer

    def _get_footer(self) -> str:
        footer = ""

        if self.length:
            if footer:
                footer += ", "
            footer += "Length: {length}".format(length=len(self.categorical))

        level_info = self.categorical._repr_categories_info()

        # Levels are added in a newline
        if footer:
            footer += "\n"
        footer += level_info

        return str(footer)

    def _get_formatted_values(self) -> List[str]:
        return format_array(
            self.categorical._internal_get_values(),
            None,
            float_format=None,
            na_rep=self.na_rep,
        )

    def to_string(self) -> str:
        categorical = self.categorical

        if len(categorical) == 0:
            if self.footer:
                return self._get_footer()
            else:
                return ""

        fmt_values = self._get_formatted_values()

        fmt_values = ["{i}".format(i=i) for i in fmt_values]
        fmt_values = [i.strip() for i in fmt_values]
        values = ", ".join(fmt_values)
        result = ["[" + values + "]"]
        if self.footer:
            footer = self._get_footer()
            if footer:
                result.append(footer)

        return str("\n".join(result))


class SeriesFormatter:
    def __init__(
        self,
        series: "Series",
        buf: Optional[IO[str]] = None,
        length: Union[bool, str] = True,
        header: bool = True,
        index: bool = True,
        na_rep: str = "NaN",
        name: bool = False,
        float_format: Optional[str] = None,
        dtype: bool = True,
        max_rows: Optional[int] = None,
        min_rows: Optional[int] = None,
    ):
        self.series = series
        self.buf = buf if buf is not None else StringIO()
        self.name = name
        self.na_rep = na_rep
        self.header = header
        self.length = length
        self.index = index
        self.max_rows = max_rows
        self.min_rows = min_rows

        if float_format is None:
            float_format = get_option("display.float_format")
        self.float_format = float_format
        self.dtype = dtype
        self.adj = _get_adjustment()

        self._chk_truncate()

    def _chk_truncate(self) -> None:
        from pandas.core.reshape.concat import concat

        self.tr_row_num: Optional[int]

        min_rows = self.min_rows
        max_rows = self.max_rows
        # truncation determined by max_rows, actual truncated number of rows
        # used below by min_rows
        truncate_v = max_rows and (len(self.series) > max_rows)
        series = self.series
        if truncate_v:
            max_rows = cast(int, max_rows)
            if min_rows:
                # if min_rows is set (not None or 0), set max_rows to minimum
                # of both
                max_rows = min(min_rows, max_rows)
            if max_rows == 1:
                row_num = max_rows
                series = series.iloc[:max_rows]
            else:
                row_num = max_rows // 2
                series = series._ensure_type(
                    concat((series.iloc[:row_num], series.iloc[-row_num:]))
                )
            self.tr_row_num = row_num
        else:
            self.tr_row_num = None
        self.tr_series = series
        self.truncate_v = truncate_v

    def _get_footer(self) -> str:
        name = self.series.name
        footer = ""

        if getattr(self.series.index, "freq", None) is not None:
            footer += "Freq: {freq}".format(freq=self.series.index.freqstr)

        if self.name is not False and name is not None:
            if footer:
                footer += ", "

            series_name = pprint_thing(name, escape_chars=("\t", "\r", "\n"))
            footer += (
                ("Name: {sname}".format(sname=series_name)) if name is not None else ""
            )

        if self.length is True or (self.length == "truncate" and self.truncate_v):
            if footer:
                footer += ", "
            footer += "Length: {length}".format(length=len(self.series))

        if self.dtype is not False and self.dtype is not None:
            name = getattr(self.tr_series.dtype, "name", None)
            if name:
                if footer:
                    footer += ", "
                footer += "dtype: {typ}".format(typ=pprint_thing(name))

        # level infos are added to the end and in a new line, like it is done
        # for Categoricals
        if is_categorical_dtype(self.tr_series.dtype):
            level_info = self.tr_series._values._repr_categories_info()
            if footer:
                footer += "\n"
            footer += level_info

        return str(footer)

    def _get_formatted_index(self) -> Tuple[List[str], bool]:
        index = self.tr_series.index
        is_multi = isinstance(index, ABCMultiIndex)

        if is_multi:
            have_header = any(name for name in index.names)
            fmt_index = index.format(names=True)
        else:
            have_header = index.name is not None
            fmt_index = index.format(name=True)
        return fmt_index, have_header

    def _get_formatted_values(self) -> List[str]:
        return format_array(
            self.tr_series._values,
            None,
            float_format=self.float_format,
            na_rep=self.na_rep,
        )

    def to_string(self) -> str:
        series = self.tr_series
        footer = self._get_footer()

        if len(series) == 0:
            return "{name}([], {footer})".format(
                name=type(self.series).__name__, footer=footer
            )

        fmt_index, have_header = self._get_formatted_index()
        fmt_values = self._get_formatted_values()

        if self.truncate_v:
            n_header_rows = 0
            row_num = self.tr_row_num
            row_num = cast(int, row_num)
            width = self.adj.len(fmt_values[row_num - 1])
            if width > 3:
                dot_str = "..."
            else:
                dot_str = ".."
            # Series uses mode=center because it has single value columns
            # DataFrame uses mode=left
            dot_str = self.adj.justify([dot_str], width, mode="center")[0]
            fmt_values.insert(row_num + n_header_rows, dot_str)
            fmt_index.insert(row_num + 1, "")

        if self.index:
            result = self.adj.adjoin(3, *[fmt_index[1:], fmt_values])
        else:
            result = self.adj.adjoin(3, fmt_values)

        if self.header and have_header:
            result = fmt_index[0] + "\n" + result

        if footer:
            result += "\n" + footer

        return str("".join(result))


class TextAdjustment:
    def __init__(self):
        self.encoding = get_option("display.encoding")

    def len(self, text: str) -> int:
        return len(text)

    def justify(self, texts: Any, max_len: int, mode: str = "right") -> List[str]:
        return justify(texts, max_len, mode=mode)

    def adjoin(self, space: int, *lists, **kwargs) -> str:
        return adjoin(space, *lists, strlen=self.len, justfunc=self.justify, **kwargs)


class EastAsianTextAdjustment(TextAdjustment):
    def __init__(self):
        super().__init__()
        if get_option("display.unicode.ambiguous_as_wide"):
            self.ambiguous_width = 2
        else:
            self.ambiguous_width = 1

        # Definition of East Asian Width
        # http://unicode.org/reports/tr11/
        # Ambiguous width can be changed by option
        self._EAW_MAP = {"Na": 1, "N": 1, "W": 2, "F": 2, "H": 1}

    def len(self, text: str) -> int:
        """
        Calculate display width considering unicode East Asian Width
        """
        if not isinstance(text, str):
            return len(text)

        return sum(
            self._EAW_MAP.get(east_asian_width(c), self.ambiguous_width) for c in text
        )

    def justify(
        self, texts: Iterable[str], max_len: int, mode: str = "right"
    ) -> List[str]:
        # re-calculate padding space per str considering East Asian Width
        def _get_pad(t):
            return max_len - self.len(t) + len(t)

        if mode == "left":
            return [x.ljust(_get_pad(x)) for x in texts]
        elif mode == "center":
            return [x.center(_get_pad(x)) for x in texts]
        else:
            return [x.rjust(_get_pad(x)) for x in texts]


def _get_adjustment() -> TextAdjustment:
    use_east_asian_width = get_option("display.unicode.east_asian_width")
    if use_east_asian_width:
        return EastAsianTextAdjustment()
    else:
        return TextAdjustment()


class TableFormatter:

    show_dimensions: Union[bool, str]
    is_truncated: bool
    formatters: formatters_type
    columns: Index

    @property
    def should_show_dimensions(self) -> bool:
        return self.show_dimensions is True or (
            self.show_dimensions == "truncate" and self.is_truncated
        )

    def _get_formatter(self, i: Union[str, int]) -> Optional[Callable]:
        if isinstance(self.formatters, (list, tuple)):
            if is_integer(i):
                i = cast(int, i)
                return self.formatters[i]
            else:
                return None
        else:
            if is_integer(i) and i not in self.columns:
                i = self.columns[i]
            return self.formatters.get(i, None)

    @contextmanager
    def get_buffer(
        self, buf: Optional[FilePathOrBuffer[str]], encoding: Optional[str] = None
    ):
        """
        Context manager to open, yield and close buffer for filenames or Path-like
        objects, otherwise yield buf unchanged.
        """
        if buf is not None:
            buf = stringify_path(buf)
        else:
            buf = StringIO()

        if encoding is None:
            encoding = "utf-8"
        elif not isinstance(buf, str):
            raise ValueError("buf is not a file name and encoding is specified.")

        if hasattr(buf, "write"):
            yield buf
        elif isinstance(buf, str):
            with open(buf, "w", encoding=encoding, newline="") as f:
                # GH#30034 open instead of codecs.open prevents a file leak
                #  if we have an invalid encoding argument.
                # newline="" is needed to roundtrip correctly on
                #  windows test_to_latex_filename
                yield f
        else:
            raise TypeError("buf is not a file name and it has no write method")

    def write_result(self, buf: IO[str]) -> None:
        """
        Write the result of serialization to buf.
        """
        raise AbstractMethodError(self)

    def get_result(
        self,
        buf: Optional[FilePathOrBuffer[str]] = None,
        encoding: Optional[str] = None,
    ) -> Optional[str]:
        """
        Perform serialization. Write to buf or return as string if buf is None.
        """
        with self.get_buffer(buf, encoding=encoding) as f:
            self.write_result(buf=f)
            if buf is None:
                return f.getvalue()
            return None


class DataFrameFormatter(TableFormatter):
    """
    Render a DataFrame

    self.to_string() : console-friendly tabular output
    self.to_html()   : html table
    self.to_latex()   : LaTeX tabular environment table

    """

    __doc__ = __doc__ if __doc__ else ""
    __doc__ += common_docstring + return_docstring

    def __init__(
        self,
        frame: "DataFrame",
        columns: Optional[Sequence[str]] = None,
        col_space: Optional[Union[str, int]] = None,
        header: Union[bool, Sequence[str]] = True,
        index: bool = True,
        na_rep: str = "NaN",
        formatters: Optional[formatters_type] = None,
        justify: Optional[str] = None,
        float_format: Optional[float_format_type] = None,
        sparsify: Optional[bool] = None,
        index_names: bool = True,
        line_width: Optional[int] = None,
        max_rows: Optional[int] = None,
        min_rows: Optional[int] = None,
        max_cols: Optional[int] = None,
        show_dimensions: Union[bool, str] = False,
        decimal: str = ".",
        table_id: Optional[str] = None,
        render_links: bool = False,
        bold_rows: bool = False,
        escape: bool = True,
    ):
        self.frame = frame
        self.show_index_names = index_names

        if sparsify is None:
            sparsify = get_option("display.multi_sparse")

        self.sparsify = sparsify

        self.float_format = float_format
        if formatters is None:
            self.formatters = {}
        elif len(frame.columns) == len(formatters) or isinstance(formatters, dict):
            self.formatters = formatters
        else:
            raise ValueError(
                (
                    "Formatters length({flen}) should match "
                    "DataFrame number of columns({dlen})"
                ).format(flen=len(formatters), dlen=len(frame.columns))
            )
        self.na_rep = na_rep
        self.decimal = decimal
        self.col_space = col_space
        self.header = header
        self.index = index
        self.line_width = line_width
        self.max_rows = max_rows
        self.min_rows = min_rows
        self.max_cols = max_cols
        self.max_rows_displayed = min(max_rows or len(self.frame), len(self.frame))
        self.show_dimensions = show_dimensions
        self.table_id = table_id
        self.render_links = render_links

        if justify is None:
            self.justify = get_option("display.colheader_justify")
        else:
            self.justify = justify

        self.bold_rows = bold_rows
        self.escape = escape

        if columns is not None:
            self.columns = ensure_index(columns)
            self.frame = self.frame[self.columns]
        else:
            self.columns = frame.columns

        self._chk_truncate()
        self.adj = _get_adjustment()

    def _chk_truncate(self) -> None:
        """
        Checks whether the frame should be truncated. If so, slices
        the frame up.
        """
        from pandas.core.reshape.concat import concat

        # Cut the data to the information actually printed
        max_cols = self.max_cols
        max_rows = self.max_rows
        self.max_rows_adj: Optional[int]
        max_rows_adj: Optional[int]

        if max_cols == 0 or max_rows == 0:  # assume we are in the terminal
            (w, h) = get_terminal_size()
            self.w = w
            self.h = h
            if self.max_rows == 0:
                dot_row = 1
                prompt_row = 1
                if self.show_dimensions:
                    show_dimension_rows = 3
                # assume we only get here if self.header is boolean.
                # i.e. not to_latex() where self.header may be List[str]
                self.header = cast(bool, self.header)
                n_add_rows = self.header + dot_row + show_dimension_rows + prompt_row
                # rows available to fill with actual data
                max_rows_adj = self.h - n_add_rows
                self.max_rows_adj = max_rows_adj

            # Format only rows and columns that could potentially fit the
            # screen
            if max_cols == 0 and len(self.frame.columns) > w:
                max_cols = w
            if max_rows == 0 and len(self.frame) > h:
                max_rows = h

        if not hasattr(self, "max_rows_adj"):
            if max_rows:
                if (len(self.frame) > max_rows) and self.min_rows:
                    # if truncated, set max_rows showed to min_rows
                    max_rows = min(self.min_rows, max_rows)
            self.max_rows_adj = max_rows
        if not hasattr(self, "max_cols_adj"):
            self.max_cols_adj = max_cols

        max_cols_adj = self.max_cols_adj
        max_rows_adj = self.max_rows_adj

        truncate_h = max_cols_adj and (len(self.columns) > max_cols_adj)
        truncate_v = max_rows_adj and (len(self.frame) > max_rows_adj)

        frame = self.frame
        if truncate_h:
            # cast here since if truncate_h is True, max_cols_adj is not None
            max_cols_adj = cast(int, max_cols_adj)
            if max_cols_adj == 0:
                col_num = len(frame.columns)
            elif max_cols_adj == 1:
                max_cols = cast(int, max_cols)
                frame = frame.iloc[:, :max_cols]
                col_num = max_cols
            else:
                col_num = max_cols_adj // 2
                frame = concat(
                    (frame.iloc[:, :col_num], frame.iloc[:, -col_num:]), axis=1
                )
                # truncate formatter
                if isinstance(self.formatters, (list, tuple)):
                    truncate_fmt = self.formatters
                    self.formatters = [
                        *truncate_fmt[:col_num],
                        *truncate_fmt[-col_num:],
                    ]
            self.tr_col_num = col_num
        if truncate_v:
            # cast here since if truncate_v is True, max_rows_adj is not None
            max_rows_adj = cast(int, max_rows_adj)
            if max_rows_adj == 1:
                row_num = max_rows
                frame = frame.iloc[:max_rows, :]
            else:
                row_num = max_rows_adj // 2
                frame = concat((frame.iloc[:row_num, :], frame.iloc[-row_num:, :]))
            self.tr_row_num = row_num
        else:
            self.tr_row_num = None

        self.tr_frame = frame
        self.truncate_h = truncate_h
        self.truncate_v = truncate_v
        self.is_truncated = bool(self.truncate_h or self.truncate_v)

    def _to_str_columns(self) -> List[List[str]]:
        """
        Render a DataFrame to a list of columns (as lists of strings).
        """
        # this method is not used by to_html where self.col_space
        # could be a string so safe to cast
        self.col_space = cast(int, self.col_space)

        frame = self.tr_frame
        # may include levels names also

        str_index = self._get_formatted_index(frame)

        if not is_list_like(self.header) and not self.header:
            stringified = []
            for i, c in enumerate(frame):
                fmt_values = self._format_col(i)
                fmt_values = _make_fixed_width(
                    fmt_values,
                    self.justify,
                    minimum=(self.col_space or 0),
                    adj=self.adj,
                )
                stringified.append(fmt_values)
        else:
            if is_list_like(self.header):
                # cast here since can't be bool if is_list_like
                self.header = cast(List[str], self.header)
                if len(self.header) != len(self.columns):
                    raise ValueError(
                        (
                            "Writing {ncols} cols but got {nalias} "
                            "aliases".format(
                                ncols=len(self.columns), nalias=len(self.header)
                            )
                        )
                    )
                str_columns = [[label] for label in self.header]
            else:
                str_columns = self._get_formatted_column_labels(frame)

            if self.show_row_idx_names:
                for x in str_columns:
                    x.append("")

            stringified = []
            for i, c in enumerate(frame):
                cheader = str_columns[i]
                header_colwidth = max(
                    self.col_space or 0, *(self.adj.len(x) for x in cheader)
                )
                fmt_values = self._format_col(i)
                fmt_values = _make_fixed_width(
                    fmt_values, self.justify, minimum=header_colwidth, adj=self.adj
                )

                max_len = max(max(self.adj.len(x) for x in fmt_values), header_colwidth)
                cheader = self.adj.justify(cheader, max_len, mode=self.justify)
                stringified.append(cheader + fmt_values)

        strcols = stringified
        if self.index:
            strcols.insert(0, str_index)

        # Add ... to signal truncated
        truncate_h = self.truncate_h
        truncate_v = self.truncate_v

        if truncate_h:
            col_num = self.tr_col_num
            strcols.insert(self.tr_col_num + 1, [" ..."] * (len(str_index)))
        if truncate_v:
            n_header_rows = len(str_index) - len(frame)
            row_num = self.tr_row_num
            # cast here since if truncate_v is True, self.tr_row_num is not None
            row_num = cast(int, row_num)
            for ix, col in enumerate(strcols):
                # infer from above row
                cwidth = self.adj.len(strcols[ix][row_num])
                is_dot_col = False
                if truncate_h:
                    is_dot_col = ix == col_num + 1
                if cwidth > 3 or is_dot_col:
                    my_str = "..."
                else:
                    my_str = ".."

                if ix == 0:
                    dot_mode = "left"
                elif is_dot_col:
                    cwidth = 4
                    dot_mode = "right"
                else:
                    dot_mode = "right"
                dot_str = self.adj.justify([my_str], cwidth, mode=dot_mode)[0]
                strcols[ix].insert(row_num + n_header_rows, dot_str)
        return strcols

    def write_result(self, buf: IO[str]) -> None:
        """
        Render a DataFrame to a console-friendly tabular output.
        """
        from pandas import Series

        frame = self.frame

        if len(frame.columns) == 0 or len(frame.index) == 0:
            info_line = "Empty {name}\nColumns: {col}\nIndex: {idx}".format(
                name=type(self.frame).__name__,
                col=pprint_thing(frame.columns),
                idx=pprint_thing(frame.index),
            )
            text = info_line
        else:

            strcols = self._to_str_columns()
            if self.line_width is None:  # no need to wrap around just print
                # the whole frame
                text = self.adj.adjoin(1, *strcols)
            elif (
                not isinstance(self.max_cols, int) or self.max_cols > 0
            ):  # need to wrap around
                text = self._join_multiline(*strcols)
            else:  # max_cols == 0. Try to fit frame to terminal
                lines = self.adj.adjoin(1, *strcols).split("\n")
                max_len = Series(lines).str.len().max()
                # plus truncate dot col
                dif = max_len - self.w
                # '+ 1' to avoid too wide repr (GH PR #17023)
                adj_dif = dif + 1
                col_lens = Series([Series(ele).apply(len).max() for ele in strcols])
                n_cols = len(col_lens)
                counter = 0
                while adj_dif > 0 and n_cols > 1:
                    counter += 1
                    mid = int(round(n_cols / 2.0))
                    mid_ix = col_lens.index[mid]
                    col_len = col_lens[mid_ix]
                    # adjoin adds one
                    adj_dif -= col_len + 1
                    col_lens = col_lens.drop(mid_ix)
                    n_cols = len(col_lens)
                # subtract index column
                max_cols_adj = n_cols - self.index
                # GH-21180. Ensure that we print at least two.
                max_cols_adj = max(max_cols_adj, 2)
                self.max_cols_adj = max_cols_adj

                # Call again _chk_truncate to cut frame appropriately
                # and then generate string representation
                self._chk_truncate()
                strcols = self._to_str_columns()
                text = self.adj.adjoin(1, *strcols)
        buf.writelines(text)

        if self.should_show_dimensions:
            buf.write(
                "\n\n[{nrows} rows x {ncols} columns]".format(
                    nrows=len(frame), ncols=len(frame.columns)
                )
            )

    def _join_multiline(self, *args) -> str:
        lwidth = self.line_width
        adjoin_width = 1
        strcols = list(args)
        if self.index:
            idx = strcols.pop(0)
            lwidth -= np.array([self.adj.len(x) for x in idx]).max() + adjoin_width

        col_widths = [
            np.array([self.adj.len(x) for x in col]).max() if len(col) > 0 else 0
            for col in strcols
        ]

        assert lwidth is not None
        col_bins = _binify(col_widths, lwidth)
        nbins = len(col_bins)

        if self.truncate_v:
            # cast here since if truncate_v is True, max_rows_adj is not None
            self.max_rows_adj = cast(int, self.max_rows_adj)
            nrows = self.max_rows_adj + 1
        else:
            nrows = len(self.frame)

        str_lst = []
        st = 0
        for i, ed in enumerate(col_bins):
            row = strcols[st:ed]
            if self.index:
                row.insert(0, idx)
            if nbins > 1:
                if ed <= len(strcols) and i < nbins - 1:
                    row.append([" \\"] + ["  "] * (nrows - 1))
                else:
                    row.append([" "] * nrows)
            str_lst.append(self.adj.adjoin(adjoin_width, *row))
            st = ed
        return "\n\n".join(str_lst)

    def to_string(
        self,
        buf: Optional[FilePathOrBuffer[str]] = None,
        encoding: Optional[str] = None,
    ) -> Optional[str]:
        return self.get_result(buf=buf, encoding=encoding)

    def to_latex(
        self,
        buf: Optional[FilePathOrBuffer[str]] = None,
        column_format: Optional[str] = None,
        longtable: bool = False,
        encoding: Optional[str] = None,
        multicolumn: bool = False,
        multicolumn_format: Optional[str] = None,
        multirow: bool = False,
        caption: Optional[str] = None,
        label: Optional[str] = None,
    ) -> Optional[str]:
        """
        Render a DataFrame to a LaTeX tabular/longtable environment output.
        """

        from pandas.io.formats.latex import LatexFormatter

        return LatexFormatter(
            self,
            column_format=column_format,
            longtable=longtable,
            multicolumn=multicolumn,
            multicolumn_format=multicolumn_format,
            multirow=multirow,
            caption=caption,
            label=label,
        ).get_result(buf=buf, encoding=encoding)

    def _format_col(self, i: int) -> List[str]:
        frame = self.tr_frame
        formatter = self._get_formatter(i)
        return format_array(
            frame.iloc[:, i]._values,
            formatter,
            float_format=self.float_format,
            na_rep=self.na_rep,
            space=self.col_space,
            decimal=self.decimal,
        )

    def to_html(
        self,
        buf: Optional[FilePathOrBuffer[str]] = None,
        encoding: Optional[str] = None,
        classes: Optional[Union[str, List, Tuple]] = None,
        notebook: bool = False,
        border: Optional[int] = None,
    ) -> Optional[str]:
        """
        Render a DataFrame to a html table.

        Parameters
        ----------
        classes : str or list-like
            classes to include in the `class` attribute of the opening
            ``<table>`` tag, in addition to the default "dataframe".
        notebook : {True, False}, optional, default False
            Whether the generated HTML is for IPython Notebook.
        border : int
            A ``border=border`` attribute is included in the opening
            ``<table>`` tag. Default ``pd.options.display.html.border``.
         """
        from pandas.io.formats.html import HTMLFormatter, NotebookFormatter

        Klass = NotebookFormatter if notebook else HTMLFormatter
        return Klass(self, classes=classes, border=border).get_result(
            buf=buf, encoding=encoding
        )

    def _get_formatted_column_labels(self, frame: "DataFrame") -> List[List[str]]:
        from pandas.core.indexes.multi import _sparsify

        columns = frame.columns

        if isinstance(columns, ABCMultiIndex):
            fmt_columns = columns.format(sparsify=False, adjoin=False)
            fmt_columns = list(zip(*fmt_columns))
            dtypes = self.frame.dtypes._values

            # if we have a Float level, they don't use leading space at all
            restrict_formatting = any(l.is_floating for l in columns.levels)
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))

            def space_format(x, y):
                if (
                    y not in self.formatters
                    and need_leadsp[x]
                    and not restrict_formatting
                ):
                    return " " + y
                return y

            str_columns = list(
                zip(*[[space_format(x, y) for y in x] for x in fmt_columns])
            )
            if self.sparsify and len(str_columns):
                str_columns = _sparsify(str_columns)

            str_columns = [list(x) for x in zip(*str_columns)]
        else:
            fmt_columns = columns.format()
            dtypes = self.frame.dtypes
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))
            str_columns = [
                [" " + x if not self._get_formatter(i) and need_leadsp[x] else x]
                for i, (col, x) in enumerate(zip(columns, fmt_columns))
            ]
        # self.str_columns = str_columns
        return str_columns

    @property
    def has_index_names(self) -> bool:
        return _has_names(self.frame.index)

    @property
    def has_column_names(self) -> bool:
        return _has_names(self.frame.columns)

    @property
    def show_row_idx_names(self) -> bool:
        return all((self.has_index_names, self.index, self.show_index_names))

    @property
    def show_col_idx_names(self) -> bool:
        return all((self.has_column_names, self.show_index_names, self.header))

    def _get_formatted_index(self, frame: "DataFrame") -> List[str]:
        # Note: this is only used by to_string() and to_latex(), not by
        # to_html(). so safe to cast col_space here.
        self.col_space = cast(int, self.col_space)
        index = frame.index
        columns = frame.columns
        fmt = self._get_formatter("__index__")

        if isinstance(index, ABCMultiIndex):
            fmt_index = index.format(
                sparsify=self.sparsify,
                adjoin=False,
                names=self.show_row_idx_names,
                formatter=fmt,
            )
        else:
            fmt_index = [index.format(name=self.show_row_idx_names, formatter=fmt)]

        fmt_index = [
            tuple(
                _make_fixed_width(
                    list(x), justify="left", minimum=(self.col_space or 0), adj=self.adj
                )
            )
            for x in fmt_index
        ]

        adjoined = self.adj.adjoin(1, *fmt_index).split("\n")

        # empty space for columns
        if self.show_col_idx_names:
            col_header = ["{x}".format(x=x) for x in self._get_column_name_list()]
        else:
            col_header = [""] * columns.nlevels

        if self.header:
            return col_header + adjoined
        else:
            return adjoined

    def _get_column_name_list(self) -> List[str]:
        names: List[str] = []
        columns = self.frame.columns
        if isinstance(columns, ABCMultiIndex):
            names.extend("" if name is None else name for name in columns.names)
        else:
            names.append("" if columns.name is None else columns.name)
        return names


# ----------------------------------------------------------------------
# Array formatters


def format_array(
    values: Any,
    formatter: Optional[Callable],
    float_format: Optional[float_format_type] = None,
    na_rep: str = "NaN",
    digits: Optional[int] = None,
    space: Optional[Union[str, int]] = None,
    justify: str = "right",
    decimal: str = ".",
    leading_space: Optional[bool] = None,
) -> List[str]:
    """
    Format an array for printing.

    Parameters
    ----------
    values
    formatter
    float_format
    na_rep
    digits
    space
    justify
    decimal
    leading_space : bool, optional
        Whether the array should be formatted with a leading space.
        When an array as a column of a Series or DataFrame, we do want
        the leading space to pad between columns.

        When formatting an Index subclass
        (e.g. IntervalIndex._format_native_types), we don't want the
        leading space since it should be left-aligned.

    Returns
    -------
    List[str]
    """

    fmt_klass: Type[GenericArrayFormatter]
    if is_datetime64_dtype(values.dtype):
        fmt_klass = Datetime64Formatter
    elif is_datetime64tz_dtype(values):
        fmt_klass = Datetime64TZFormatter
    elif is_timedelta64_dtype(values.dtype):
        fmt_klass = Timedelta64Formatter
    elif is_extension_array_dtype(values.dtype):
        fmt_klass = ExtensionArrayFormatter
    elif is_float_dtype(values.dtype) or is_complex_dtype(values.dtype):
        fmt_klass = FloatArrayFormatter
    elif is_integer_dtype(values.dtype):
        fmt_klass = IntArrayFormatter
    else:
        fmt_klass = GenericArrayFormatter

    if space is None:
        space = get_option("display.column_space")

    if float_format is None:
        float_format = get_option("display.float_format")

    if digits is None:
        digits = get_option("display.precision")

    fmt_obj = fmt_klass(
        values,
        digits=digits,
        na_rep=na_rep,
        float_format=float_format,
        formatter=formatter,
        space=space,
        justify=justify,
        decimal=decimal,
        leading_space=leading_space,
    )

    return fmt_obj.get_result()


class GenericArrayFormatter:
    def __init__(
        self,
        values: Any,
        digits: int = 7,
        formatter: Optional[Callable] = None,
        na_rep: str = "NaN",
        space: Union[str, int] = 12,
        float_format: Optional[float_format_type] = None,
        justify: str = "right",
        decimal: str = ".",
        quoting: Optional[int] = None,
        fixed_width: bool = True,
        leading_space: Optional[bool] = None,
    ):
        self.values = values
        self.digits = digits
        self.na_rep = na_rep
        self.space = space
        self.formatter = formatter
        self.float_format = float_format
        self.justify = justify
        self.decimal = decimal
        self.quoting = quoting
        self.fixed_width = fixed_width
        self.leading_space = leading_space

    def get_result(self) -> List[str]:
        fmt_values = self._format_strings()
        return _make_fixed_width(fmt_values, self.justify)

    def _format_strings(self) -> List[str]:
        if self.float_format is None:
            float_format = get_option("display.float_format")
            if float_format is None:
                fmt_str = "{{x: .{prec:d}g}}".format(
                    prec=get_option("display.precision")
                )
                float_format = lambda x: fmt_str.format(x=x)
        else:
            float_format = self.float_format

        formatter = (
            self.formatter
            if self.formatter is not None
            else (lambda x: pprint_thing(x, escape_chars=("\t", "\r", "\n")))
        )

        def _format(x):
            if self.na_rep is not None and is_scalar(x) and isna(x):
                try:
                    # try block for np.isnat specifically
                    # determine na_rep if x is None or NaT-like
                    if x is None:
                        return "None"
                    elif x is NA:
                        return str(NA)
                    elif x is NaT or np.isnat(x):
                        return "NaT"
                except (TypeError, ValueError):
                    # np.isnat only handles datetime or timedelta objects
                    pass
                return self.na_rep
            elif isinstance(x, PandasObject):
                return "{x}".format(x=x)
            else:
                # object dtype
                return "{x}".format(x=formatter(x))

        vals = self.values
        if isinstance(vals, Index):
            vals = vals._values
        elif isinstance(vals, ABCSparseArray):
            vals = vals.values

        is_float_type = lib.map_infer(vals, is_float) & notna(vals)
        leading_space = self.leading_space
        if leading_space is None:
            leading_space = is_float_type.any()

        fmt_values = []
        for i, v in enumerate(vals):
            if not is_float_type[i] and leading_space:
                fmt_values.append(" {v}".format(v=_format(v)))
            elif is_float_type[i]:
                fmt_values.append(float_format(v))
            else:
                if leading_space is False:
                    # False specifically, so that the default is
                    # to include a space if we get here.
                    tpl = "{v}"
                else:
                    tpl = " {v}"
                fmt_values.append(tpl.format(v=_format(v)))

        return fmt_values


class FloatArrayFormatter(GenericArrayFormatter):
    """

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # float_format is expected to be a string
        # formatter should be used to pass a function
        if self.float_format is not None and self.formatter is None:
            # GH21625, GH22270
            self.fixed_width = False
            if callable(self.float_format):
                self.formatter = self.float_format
                self.float_format = None

    def _value_formatter(
        self,
        float_format: Optional[float_format_type] = None,
        threshold: Optional[Union[float, int]] = None,
    ) -> Callable:
        """Returns a function to be applied on each value to format it
        """

        # the float_format parameter supersedes self.float_format
        if float_format is None:
            float_format = self.float_format

        # we are going to compose different functions, to first convert to
        # a string, then replace the decimal symbol, and finally chop according
        # to the threshold

        # when there is no float_format, we use str instead of '%g'
        # because str(0.0) = '0.0' while '%g' % 0.0 = '0'
        if float_format:

            def base_formatter(v):
                return float_format(value=v) if notna(v) else self.na_rep

        else:

            def base_formatter(v):
                return str(v) if notna(v) else self.na_rep

        if self.decimal != ".":

            def decimal_formatter(v):
                return base_formatter(v).replace(".", self.decimal, 1)

        else:
            decimal_formatter = base_formatter

        if threshold is None:
            return decimal_formatter

        def formatter(value):
            if notna(value):
                if abs(value) > threshold:
                    return decimal_formatter(value)
                else:
                    return decimal_formatter(0.0)
            else:
                return self.na_rep

        return formatter

    def get_result_as_array(self) -> np.ndarray:
        """
        Returns the float values converted into strings using
        the parameters given at initialisation, as a numpy array
        """

        if self.formatter is not None:
            return np.array([self.formatter(x) for x in self.values])

        if self.fixed_width:
            threshold = get_option("display.chop_threshold")
        else:
            threshold = None

        # if we have a fixed_width, we'll need to try different float_format
        def format_values_with(float_format):
            formatter = self._value_formatter(float_format, threshold)

            # default formatter leaves a space to the left when formatting
            # floats, must be consistent for left-justifying NaNs (GH #25061)
            if self.justify == "left":
                na_rep = " " + self.na_rep
            else:
                na_rep = self.na_rep

            # separate the wheat from the chaff
            values = self.values
            is_complex = is_complex_dtype(values)
            mask = isna(values)
            if hasattr(values, "to_dense"):  # sparse numpy ndarray
                values = values.to_dense()
            values = np.array(values, dtype="object")
            values[mask] = na_rep
            imask = (~mask).ravel()
            values.flat[imask] = np.array(
                [formatter(val) for val in values.ravel()[imask]]
            )

            if self.fixed_width:
                if is_complex:
                    result = _trim_zeros_complex(values, na_rep)
                else:
                    result = _trim_zeros_float(values, na_rep)
                return np.asarray(result, dtype="object")

            return values

        # There is a special default string when we are fixed-width
        # The default is otherwise to use str instead of a formatting string
        float_format: Optional[float_format_type]
        if self.float_format is None:
            if self.fixed_width:
                float_format = partial(
                    "{value: .{digits:d}f}".format, digits=self.digits
                )
            else:
                float_format = self.float_format
        else:
            float_format = lambda value: self.float_format % value

        formatted_values = format_values_with(float_format)

        if not self.fixed_width:
            return formatted_values

        # we need do convert to engineering format if some values are too small
        # and would appear as 0, or if some values are too big and take too
        # much space

        if len(formatted_values) > 0:
            maxlen = max(len(x) for x in formatted_values)
            too_long = maxlen > self.digits + 6
        else:
            too_long = False

        with np.errstate(invalid="ignore"):
            abs_vals = np.abs(self.values)
            # this is pretty arbitrary for now
            # large values: more that 8 characters including decimal symbol
            # and first digit, hence > 1e6
            has_large_values = (abs_vals > 1e6).any()
            has_small_values = (
                (abs_vals < 10 ** (-self.digits)) & (abs_vals > 0)
            ).any()

        if has_small_values or (too_long and has_large_values):
            float_format = partial("{value: .{digits:d}e}".format, digits=self.digits)
            formatted_values = format_values_with(float_format)

        return formatted_values

    def _format_strings(self) -> List[str]:
        # shortcut
        if self.formatter is not None:
            return [self.formatter(x) for x in self.values]

        return list(self.get_result_as_array())


class IntArrayFormatter(GenericArrayFormatter):
    def _format_strings(self) -> List[str]:
        formatter = self.formatter or (lambda x: "{x: d}".format(x=x))
        fmt_values = [formatter(x) for x in self.values]
        return fmt_values


class Datetime64Formatter(GenericArrayFormatter):
    def __init__(
        self,
        values: Union[np.ndarray, "Series", DatetimeIndex, DatetimeArray],
        nat_rep: str = "NaT",
        date_format: None = None,
        **kwargs,
    ):
        super().__init__(values, **kwargs)
        self.nat_rep = nat_rep
        self.date_format = date_format

    def _format_strings(self) -> List[str]:
        """ we by definition have DO NOT have a TZ """

        values = self.values

        if not isinstance(values, DatetimeIndex):
            values = DatetimeIndex(values)

        if self.formatter is not None and callable(self.formatter):
            return [self.formatter(x) for x in values]

        fmt_values = format_array_from_datetime(
            values.asi8.ravel(),
            format=_get_format_datetime64_from_values(values, self.date_format),
            na_rep=self.nat_rep,
        ).reshape(values.shape)
        return fmt_values.tolist()


class ExtensionArrayFormatter(GenericArrayFormatter):
    def _format_strings(self) -> List[str]:
        values = self.values
        if isinstance(values, (ABCIndexClass, ABCSeries)):
            values = values._values

        formatter = values._formatter(boxed=True)

        if is_categorical_dtype(values.dtype):
            # Categorical is special for now, so that we can preserve tzinfo
            array = values._internal_get_values()
        else:
            array = np.asarray(values)

        fmt_values = format_array(
            array,
            formatter,
            float_format=self.float_format,
            na_rep=self.na_rep,
            digits=self.digits,
            space=self.space,
            justify=self.justify,
            leading_space=self.leading_space,
        )
        return fmt_values


def format_percentiles(
    percentiles: Union[
        np.ndarray, List[Union[int, float]], List[float], List[Union[str, float]]
    ]
) -> List[str]:
    """
    Outputs rounded and formatted percentiles.

    Parameters
    ----------
    percentiles : list-like, containing floats from interval [0,1]

    Returns
    -------
    formatted : list of strings

    Notes
    -----
    Rounding precision is chosen so that: (1) if any two elements of
    ``percentiles`` differ, they remain different after rounding
    (2) no entry is *rounded* to 0% or 100%.
    Any non-integer is always rounded to at least 1 decimal place.

    Examples
    --------
    Keeps all entries different after rounding:

    >>> format_percentiles([0.01999, 0.02001, 0.5, 0.666666, 0.9999])
    ['1.999%', '2.001%', '50%', '66.667%', '99.99%']

    No element is rounded to 0% or 100% (unless already equal to it).
    Duplicates are allowed:

    >>> format_percentiles([0, 0.5, 0.02001, 0.5, 0.666666, 0.9999])
    ['0%', '50%', '2.0%', '50%', '66.67%', '99.99%']
    """

    percentiles = np.asarray(percentiles)

    # It checks for np.NaN as well
    with np.errstate(invalid="ignore"):
        if (
            not is_numeric_dtype(percentiles)
            or not np.all(percentiles >= 0)
            or not np.all(percentiles <= 1)
        ):
            raise ValueError("percentiles should all be in the interval [0,1]")

    percentiles = 100 * percentiles
    int_idx = np.isclose(percentiles.astype(int), percentiles)

    if np.all(int_idx):
        out = percentiles.astype(int).astype(str)
        return [i + "%" for i in out]

    unique_pcts = np.unique(percentiles)
    to_begin = unique_pcts[0] if unique_pcts[0] > 0 else None
    to_end = 100 - unique_pcts[-1] if unique_pcts[-1] < 100 else None

    # Least precision that keeps percentiles unique after rounding
    prec = -np.floor(
        np.log10(np.min(np.ediff1d(unique_pcts, to_begin=to_begin, to_end=to_end)))
    ).astype(int)
    prec = max(1, prec)
    out = np.empty_like(percentiles, dtype=object)
    out[int_idx] = percentiles[int_idx].astype(int).astype(str)
    out[~int_idx] = percentiles[~int_idx].round(prec).astype(str)
    return [i + "%" for i in out]


def _is_dates_only(
    values: Union[np.ndarray, DatetimeArray, Index, DatetimeIndex]
) -> bool:
    # return a boolean if we are only dates (and don't have a timezone)
    assert values.ndim == 1

    values = DatetimeIndex(values)
    if values.tz is not None:
        return False

    values_int = values.asi8
    consider_values = values_int != iNaT
    one_day_nanos = 86400 * 1e9
    even_days = (
        np.logical_and(consider_values, values_int % int(one_day_nanos) != 0).sum() == 0
    )
    if even_days:
        return True
    return False


def _format_datetime64(
    x: Union[NaTType, Timestamp], tz: Optional[tzinfo] = None, nat_rep: str = "NaT"
) -> str:
    if x is None or (is_scalar(x) and isna(x)):
        return nat_rep

    if tz is not None or not isinstance(x, Timestamp):
        if getattr(x, "tzinfo", None) is not None:
            x = Timestamp(x).tz_convert(tz)
        else:
            x = Timestamp(x).tz_localize(tz)

    return str(x)


def _format_datetime64_dateonly(
    x: Union[NaTType, Timestamp], nat_rep: str = "NaT", date_format: None = None
) -> str:
    if x is None or (is_scalar(x) and isna(x)):
        return nat_rep

    if not isinstance(x, Timestamp):
        x = Timestamp(x)

    if date_format:
        return x.strftime(date_format)
    else:
        return x._date_repr


def _get_format_datetime64(
    is_dates_only: bool, nat_rep: str = "NaT", date_format: None = None
) -> Callable:

    if is_dates_only:
        return lambda x, tz=None: _format_datetime64_dateonly(
            x, nat_rep=nat_rep, date_format=date_format
        )
    else:
        return lambda x, tz=None: _format_datetime64(x, tz=tz, nat_rep=nat_rep)


def _get_format_datetime64_from_values(
    values: Union[np.ndarray, DatetimeArray, DatetimeIndex], date_format: Optional[str]
) -> Optional[str]:
    """ given values and a date_format, return a string format """

    if isinstance(values, np.ndarray) and values.ndim > 1:
        # We don't actually care about the order of values, and DatetimeIndex
        #  only accepts 1D values
        values = values.ravel()

    is_dates_only = _is_dates_only(values)
    if is_dates_only:
        return date_format or "%Y-%m-%d"
    return date_format


class Datetime64TZFormatter(Datetime64Formatter):
    def _format_strings(self) -> List[str]:
        """ we by definition have a TZ """

        values = self.values.astype(object)
        is_dates_only = _is_dates_only(values)
        formatter = self.formatter or _get_format_datetime64(
            is_dates_only, date_format=self.date_format
        )
        fmt_values = [formatter(x) for x in values]

        return fmt_values


class Timedelta64Formatter(GenericArrayFormatter):
    def __init__(
        self,
        values: Union[np.ndarray, TimedeltaIndex],
        nat_rep: str = "NaT",
        box: bool = False,
        **kwargs,
    ):
        super().__init__(values, **kwargs)
        self.nat_rep = nat_rep
        self.box = box

    def _format_strings(self) -> List[str]:
        formatter = self.formatter or _get_format_timedelta64(
            self.values, nat_rep=self.nat_rep, box=self.box
        )
        return [formatter(x) for x in self.values]


def _get_format_timedelta64(
    values: Union[np.ndarray, TimedeltaIndex, TimedeltaArray],
    nat_rep: str = "NaT",
    box: bool = False,
) -> Callable:
    """
    Return a formatter function for a range of timedeltas.
    These will all have the same format argument

    If box, then show the return in quotes
    """

    values_int = values.astype(np.int64)

    consider_values = values_int != iNaT

    one_day_nanos = 86400 * 1e9
    even_days = (
        np.logical_and(consider_values, values_int % one_day_nanos != 0).sum() == 0
    )
    all_sub_day = (
        np.logical_and(consider_values, np.abs(values_int) >= one_day_nanos).sum() == 0
    )

    if even_days:
        format = None
    elif all_sub_day:
        format = "sub_day"
    else:
        format = "long"

    def _formatter(x):
        if x is None or (is_scalar(x) and isna(x)):
            return nat_rep

        if not isinstance(x, Timedelta):
            x = Timedelta(x)
        result = x._repr_base(format=format)
        if box:
            result = "'{res}'".format(res=result)
        return result

    return _formatter


def _make_fixed_width(
    strings: List[str],
    justify: str = "right",
    minimum: Optional[int] = None,
    adj: Optional[TextAdjustment] = None,
) -> List[str]:

    if len(strings) == 0 or justify == "all":
        return strings

    if adj is None:
        adj = _get_adjustment()

    max_len = max(adj.len(x) for x in strings)

    if minimum is not None:
        max_len = max(minimum, max_len)

    conf_max = get_option("display.max_colwidth")
    if conf_max is not None and max_len > conf_max:
        max_len = conf_max

    def just(x):
        if conf_max is not None:
            if (conf_max > 3) & (adj.len(x) > max_len):
                x = x[: max_len - 3] + "..."
        return x

    strings = [just(x) for x in strings]
    result = adj.justify(strings, max_len, mode=justify)
    return result


def _trim_zeros_complex(str_complexes: np.ndarray, na_rep: str = "NaN") -> List[str]:
    """
    Separates the real and imaginary parts from the complex number, and
    executes the _trim_zeros_float method on each of those.
    """
    return [
        "".join(_trim_zeros_float(re.split(r"([j+-])", x), na_rep))
        for x in str_complexes
    ]


def _trim_zeros_float(
    str_floats: Union[np.ndarray, List[str]], na_rep: str = "NaN"
) -> List[str]:
    """
    Trims zeros, leaving just one before the decimal points if need be.
    """
    trimmed = str_floats

    def _is_number(x):
        return x != na_rep and not x.endswith("inf")

    def _cond(values):
        finite = [x for x in values if _is_number(x)]
        return (
            len(finite) > 0
            and all(x.endswith("0") for x in finite)
            and not (any(("e" in x) or ("E" in x) for x in finite))
        )

    while _cond(trimmed):
        trimmed = [x[:-1] if _is_number(x) else x for x in trimmed]

    # leave one 0 after the decimal points if need be.
    return [x + "0" if x.endswith(".") and _is_number(x) else x for x in trimmed]


def _has_names(index: Index) -> bool:
    if isinstance(index, ABCMultiIndex):
        return com.any_not_none(*index.names)
    else:
        return index.name is not None


class EngFormatter:
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
        24: "Y",
    }

    def __init__(self, accuracy: Optional[int] = None, use_eng_prefix: bool = False):
        self.accuracy = accuracy
        self.use_eng_prefix = use_eng_prefix

    def __call__(self, num: Union[int, float]) -> str:
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
        dnum = decimal.Decimal(str(num))

        if decimal.Decimal.is_nan(dnum):
            return "NaN"

        if decimal.Decimal.is_infinite(dnum):
            return "inf"

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
                prefix = "E-{pow10:02d}".format(pow10=-int_pow10)
            else:
                prefix = "E+{pow10:02d}".format(pow10=int_pow10)

        mant = sign * dnum / (10 ** pow10)

        if self.accuracy is None:  # pragma: no cover
            format_str = "{mant: g}{prefix}"
        else:
            format_str = "{{mant: .{acc:d}f}}{{prefix}}".format(acc=self.accuracy)

        formatted = format_str.format(mant=mant, prefix=prefix)

        return formatted


def set_eng_float_format(accuracy: int = 3, use_eng_prefix: bool = False) -> None:
    """
    Alter default behavior on how float is formatted in DataFrame.
    Format float in engineering format. By accuracy, we mean the number of
    decimal digits after the floating point.

    See also EngFormatter.
    """

    set_option("display.float_format", EngFormatter(accuracy, use_eng_prefix))
    set_option("display.column_space", max(12, accuracy + 9))


def _binify(cols: List[int], line_width: int) -> List[int]:
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


def get_level_lengths(
    levels: Any, sentinel: Union[bool, object, str] = ""
) -> List[Dict[int, int]]:
    """For each index in each level the function returns lengths of indexes.

    Parameters
    ----------
    levels : list of lists
        List of values on for level.
    sentinel : string, optional
        Value which states that no new index starts on there.

    Returns
    -------
    Returns list of maps. For each level returns map of indexes (key is index
    in row and value is length of index).
    """
    if len(levels) == 0:
        return []

    control = [True] * len(levels[0])

    result = []
    for level in levels:
        last_index = 0

        lengths = {}
        for i, key in enumerate(level):
            if control[i] and key == sentinel:
                pass
            else:
                control[i] = False
                lengths[last_index] = i - last_index
                last_index = i

        lengths[last_index] = len(level) - last_index

        result.append(lengths)

    return result


def buffer_put_lines(buf: IO[str], lines: List[str]) -> None:
    """
    Appends lines to a buffer.

    Parameters
    ----------
    buf
        The buffer to write to
    lines
        The lines to append.
    """
    if any(isinstance(x, str) for x in lines):
        lines = [str(x) for x in lines]
    buf.write("\n".join(lines))
