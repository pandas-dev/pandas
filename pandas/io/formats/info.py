from abc import ABC, abstractmethod
import sys
from typing import IO, List, Optional, Union

from pandas._config import get_option

from pandas._typing import Dtype, FrameOrSeries

from pandas.core.indexes.api import Index
from pandas.core.series import Series

from pandas.io.formats import format as fmt
from pandas.io.formats.printing import pprint_thing


def _put_str(s: Union[str, Dtype], space: int) -> str:
    """
    Make string of specified length, padding to the right if necessary.

    Parameters
    ----------
    s : Union[str, Dtype]
        String to be formatted.
    space : int
        Length to force string to be of.

    Returns
    -------
    str
        String coerced to given length.

    Examples
    --------
    >>> pd.io.formats.info._put_str("panda", 6)
    'panda '
    >>> pd.io.formats.info._put_str("panda", 4)
    'pand'
    """
    return str(s)[:space].ljust(space)


def _sizeof_fmt(num: Union[int, float], size_qualifier: str) -> str:
    """
    Return size in human readable format.

    Parameters
    ----------
    num : int
        Size in bytes.
    size_qualifier : str
        Either empty, or '+' (if lower bound).

    Returns
    -------
    str
        Size in human readable format.

    Examples
    --------
    >>> _sizeof_fmt(23028, '')
    '22.5 KB'

    >>> _sizeof_fmt(23028, '+')
    '22.5+ KB'
    """
    for x in ["bytes", "KB", "MB", "GB", "TB"]:
        if num < 1024.0:
            return f"{num:3.1f}{size_qualifier} {x}"
        num /= 1024.0
    return f"{num:3.1f}{size_qualifier} PB"


class BaseInfo(ABC):
    """Base class for DataFrameInfo and SeriesInfo.

    Parameters
    ----------
    data : FrameOrSeries
        Either dataframe or series.
    memory_usage : bool or str, optional
        If "deep", introspect the data deeply by interrogating object dtypes
        for system-level memory consumption, and include it in the returned
        values.
    """

    def __init__(
        self,
        data: FrameOrSeries,
        memory_usage: Optional[Union[bool, str]] = None,
    ):
        self.data = data
        self.memory_usage = self._initialize_memory_usage(memory_usage)

    def _initialize_memory_usage(
        self,
        memory_usage: Optional[Union[bool, str]] = None,
    ) -> Union[bool, str]:
        if memory_usage is None:
            memory_usage = get_option("display.memory_usage")
        return memory_usage

    @property
    @abstractmethod
    def ids(self) -> Index:
        pass

    @property
    @abstractmethod
    def dtypes(self) -> Series:
        """Dtypes.

        Returns
        -------
        dtypes : Series
            Dtype of each of the DataFrame's columns.
        """
        return self.data.dtypes

    @property
    def mem_usage(self) -> int:
        """Memory usage in bytes.

        Returns
        -------
        mem_usage : int
            Object's total memory usage in bytes.
        """
        if self.memory_usage == "deep":
            deep = True
        else:
            deep = False
        return self.data.memory_usage(index=True, deep=deep).sum()

    @property
    def size_qualifier(self) -> str:
        size_qualifier = ""
        if self.memory_usage:
            if self.memory_usage != "deep":
                # size_qualifier is just a best effort; not guaranteed to catch
                # all cases (e.g., it misses categorical data even with object
                # categories)
                if (
                    "object" in self.counts
                    or self.data.index._is_memory_usage_qualified()
                ):
                    size_qualifier = "+"
        return size_qualifier


class DataFrameInfo(BaseInfo):
    """Class storing dataframe-specific info."""

    @property
    def ids(self) -> Index:
        """Column names.

        Returns
        -------
        ids : Index
            DataFrame's column names.
        """
        return self.data.columns

    @property
    def dtypes(self) -> Series:
        """Dtypes.

        Returns
        -------
        dtypes : Series
            Dtype of each of the DataFrame's columns.
        """
        return self.data.dtypes

    @property
    def counts(self):
        # groupby dtype.name to collect e.g. Categorical columns
        return self.dtypes.value_counts().groupby(lambda x: x.name).sum()

    def to_buffer(self, *, buf, max_cols, verbose, null_counts) -> None:
        """
        Print a concise summary of a %(klass)s.

        This method prints information about a %(klass)s including
        the index dtype%(type_sub)s, non-null values and memory usage.

        Parameters
        ----------
        data : %(klass)s
            %(klass)s to print information about.
        verbose : bool, optional
            Whether to print the full summary. By default, the setting in
            ``pandas.options.display.max_info_columns`` is followed.
        buf : writable buffer, defaults to sys.stdout
            Where to send the output. By default, the output is printed to
            sys.stdout. Pass a writable buffer if you need to further process
            the output.
        %(max_cols_sub)s
        memory_usage : bool, str, optional
            Specifies whether total memory usage of the %(klass)s
            elements (including the index) should be displayed. By default,
            this follows the ``pandas.options.display.memory_usage`` setting.

            True always show memory usage. False never shows memory usage.
            A value of 'deep' is equivalent to "True with deep introspection".
            Memory usage is shown in human-readable units (base-2
            representation). Without deep introspection a memory estimation is
            made based in column dtype and number of rows assuming values
            consume the same memory amount for corresponding dtypes. With deep
            memory introspection, a real memory usage calculation is performed
            at the cost of computational resources.
        null_counts : bool, optional
            Whether to show the non-null counts. By default, this is shown
            only if the %(klass)s is smaller than
            ``pandas.options.display.max_info_rows`` and
            ``pandas.options.display.max_info_columns``. A value of True always
            shows the counts, and False never shows the counts.

        Returns
        -------
        None
            This method prints a summary of a %(klass)s and returns None.

        See Also
        --------
        %(see_also_sub)s

        Examples
        --------
        %(examples_sub)s
        """
        printer = InfoPrinter(
            info=self,
            max_cols=max_cols,
            verbose=verbose,
            null_counts=null_counts,
        )
        printer.to_buffer(buf)


class InfoPrinter:
    """Class for printing dataframe or series info."""

    def __init__(
        self,
        info: DataFrameInfo,
        max_cols: Optional[int] = None,
        verbose: Optional[bool] = None,
        null_counts: Optional[bool] = None,
    ):
        self.info = info
        self.data = info.data
        self.max_cols = max_cols
        self.verbose = verbose
        self.null_counts = null_counts

    @property
    def max_cols(self):
        return self._max_cols

    @max_cols.setter
    def max_cols(self, max_cols):
        # hack
        if max_cols is None:
            max_cols = get_option("display.max_info_columns", self.col_count + 1)
        self._max_cols = max_cols

    @property
    def max_rows(self):
        return get_option("display.max_info_rows", len(self.data) + 1)

    @property
    def exceeds_info_cols(self):
        return self.col_count > self.max_cols

    @property
    def show_counts(self) -> bool:
        if self.null_counts is None:
            return bool(
                (self.col_count <= self.max_cols) and (len(self.data) < self.max_rows)
            )
        else:
            return self.null_counts

    @property
    def col_count(self):
        return len(self.info.ids)

    def to_buffer(self, buf: Optional[IO[str]] = None) -> None:
        klass = self._select_table_builder()
        table_builder = klass(info=self.info, printer=self)
        lines = table_builder.get_lines()
        if buf is None:  # pragma: no cover
            buf = sys.stdout
        fmt.buffer_put_lines(buf, lines)

    def _select_table_builder(self):
        if self.verbose:
            return self._select_verbose_table_builder()
        elif self.verbose is False:  # specifically set to False, not necessarily None
            return DataFrameTableBuilderNonVerbose
        else:
            if self.exceeds_info_cols:
                return DataFrameTableBuilderNonVerbose
            else:
                return self._select_verbose_table_builder()

    def _select_verbose_table_builder(self):
        if self.show_counts:
            return DataFrameTableBuilderVerboseWithCounts
        else:
            return DataFrameTableBuilderVerboseNoCounts


class TableBuilderAbstract(ABC):
    """Abstract builder for info table."""

    _lines: List[str]

    def __init__(self, *, info, printer):
        self.info = info
        self.printer = printer

    def get_lines(self):
        self._lines = []
        if self.col_count == 0:
            self._fill_empty_info()
        else:
            self._fill_non_empty_info()
        return self._lines

    @abstractmethod
    def _fill_empty_info(self):
        pass

    @abstractmethod
    def _fill_non_empty_info(self):
        pass


class DataFrameTableBuilder(TableBuilderAbstract):
    """Abstract builder for dataframe info table."""

    def _fill_empty_info(self):
        self.add_object_type_line()
        self.add_index_range_line()
        self._lines.append(f"Empty {type(self.data).__name__}")

    def _fill_non_empty_info(self):
        self.add_object_type_line()
        self.add_index_range_line()
        self.add_columns_summary_line()
        self.add_header_line()
        self.add_separator_line()
        self.add_body_lines()
        self.add_dtypes_line()
        if self.memory_usage:
            self.add_memory_usage_line()

    @property
    def data(self):
        return self.info.data

    @property
    def counts(self):
        return self.info.counts

    @property
    def memory_usage(self):
        return self.info.memory_usage

    @property
    def ids(self):
        return self.info.ids

    @property
    def dtypes(self):
        return self.info.dtypes

    @property
    def show_counts(self):
        return self.printer.show_counts

    @property
    def col_count(self):
        return self.printer.col_count

    def add_object_type_line(self):
        self._lines.append(str(type(self.data)))

    def add_index_range_line(self):
        self._lines.append(self.data.index._summary())

    @abstractmethod
    def add_columns_summary_line(self):
        pass

    @abstractmethod
    def add_header_line(self):
        pass

    @abstractmethod
    def add_separator_line(self):
        pass

    @abstractmethod
    def add_body_lines(self):
        pass

    def add_dtypes_line(self):
        collected_dtypes = [
            f"{key}({val:d})" for key, val in sorted(self.counts.items())
        ]
        self._lines.append(f"dtypes: {', '.join(collected_dtypes)}")

    def add_memory_usage_line(self):
        self._lines.append(
            "memory usage: "
            f"{_sizeof_fmt(self.info.mem_usage, self.info.size_qualifier)}\n"
        )


class DataFrameTableBuilderNonVerbose(DataFrameTableBuilder):
    """Info table builder for non-verbose output."""

    def add_columns_summary_line(self):
        self._lines.append(self.ids._summary(name="Columns"))

    def add_header_line(self):
        pass

    def add_separator_line(self):
        pass

    def add_body_lines(self):
        pass


class DataFrameTableBuilderVerbose(DataFrameTableBuilder):
    """Info table builder for verbose output."""

    COL_SPACE = 2
    SPACING = " " * COL_SPACE
    HEADERS: List[str]

    def __init__(self, *, info, printer):
        super().__init__(info=info, printer=printer)
        self.strcols: List[List[str]] = self._get_strcols()

    @abstractmethod
    def _get_strcols(self) -> List[List[str]]:
        pass

    def add_columns_summary_line(self):
        self._lines.append(f"Data columns (total {self.col_count} columns):")

    @property
    def header_column_widths(self):
        return [len(col) for col in self.HEADERS]

    @property
    def body_column_widths(self):
        return [max(len(x) for x in col) for col in self.strcols]

    @property
    def gross_column_widths(self):
        return [
            max(header_colwidth, body_colwidth)
            for header_colwidth, body_colwidth in zip(
                self.header_column_widths, self.body_column_widths
            )
        ]

    def add_header_line(self):
        header_line = self.SPACING.join(
            [
                _put_str(header, col_width)
                for header, col_width in zip(self.HEADERS, self.gross_column_widths)
            ]
        )
        self._lines.append(header_line)

    def add_separator_line(self):
        separator_line = self.SPACING.join(
            [
                _put_str("-" * header_colwidth, gross_colwidth)
                for header_colwidth, gross_colwidth in zip(
                    self.header_column_widths, self.gross_column_widths
                )
            ]
        )
        self._lines.append(separator_line)

    def add_body_lines(self):
        strrows = list(zip(*self.strcols))
        for row in strrows:
            body_line = self.SPACING.join(
                [
                    _put_str(col, gross_colwidth)
                    for col, gross_colwidth in zip(row, self.gross_column_widths)
                ]
            )
            self._lines.append(body_line)


class DataFrameTableBuilderVerboseNoCounts(DataFrameTableBuilderVerbose):
    """Verbose info table builder without non-null counts column."""

    HEADERS = [
        " # ",
        "Column",
        "Dtype",
    ]

    def _get_strcols(self) -> List[List[str]]:
        line_numbers = [f" {i}" for i, _ in enumerate(self.ids)]
        columns = [pprint_thing(col) for col in self.ids]
        dtypes = [pprint_thing(dtype) for dtype in self.dtypes]
        return [line_numbers, columns, dtypes]


class DataFrameTableBuilderVerboseWithCounts(DataFrameTableBuilderVerbose):
    """Verbose info table builder with non-null counts column."""

    HEADERS = [
        " # ",
        "Column",
        "Non-Null Count",
        "Dtype",
    ]

    @property
    def count_non_null(self):
        return "{count} non-null"

    def _get_strcols(self) -> List[List[str]]:
        line_numbers = [f" {i}" for i, _ in enumerate(self.ids)]
        columns = [pprint_thing(col) for col in self.ids]
        non_null_counts = [
            self.count_non_null.format(count=count) for count in self.data.count()
        ]
        dtypes = [pprint_thing(dtype) for dtype in self.dtypes]
        return [line_numbers, columns, non_null_counts, dtypes]
