from abc import ABCMeta, abstractmethod
import sys
from typing import IO, TYPE_CHECKING, List, Optional, Tuple, Union

from pandas._config import get_option

from pandas._typing import Dtype, FrameOrSeries

from pandas.core.indexes.api import Index

from pandas.io.formats import format as fmt
from pandas.io.formats.printing import pprint_thing

if TYPE_CHECKING:
    from pandas.core.series import Series  # noqa: F401


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


class BaseInfo(metaclass=ABCMeta):
    def __init__(
        self,
        data: FrameOrSeries,
        verbose: Optional[bool] = None,
        buf: Optional[IO[str]] = None,
        max_cols: Optional[int] = None,
        memory_usage: Optional[Union[bool, str]] = None,
        null_counts: Optional[bool] = None,
    ):
        if buf is None:  # pragma: no cover
            buf = sys.stdout
        if memory_usage is None:
            memory_usage = get_option("display.memory_usage")

        self.data = data
        self.verbose = verbose
        self.buf = buf
        self.max_cols = max_cols
        self.memory_usage = memory_usage
        self.null_counts = null_counts

    @abstractmethod
    def _get_mem_usage(self, deep: bool) -> int:
        """
        Get memory usage in bytes.

        Parameters
        ----------
        deep : bool
            If True, introspect the data deeply by interrogating object dtypes
            for system-level memory consumption, and include it in the returned
            values.

        Returns
        -------
        mem_usage : int
            Object's total memory usage in bytes.
        """
        pass

    @abstractmethod
    def _get_ids_and_dtypes(self) -> Tuple["Index", "Series"]:
        """
        Get column names and dtypes.

        Returns
        -------
        ids : Index
            DataFrame's column names.
        dtypes : Series
            Dtype of each of the DataFrame's columns.
        """
        pass

    @abstractmethod
    def _verbose_repr(
        self, lines: List[str], ids: "Index", dtypes: "Series", show_counts: bool
    ) -> None:
        """
        Append name, non-null count (optional), and dtype for each column to `lines`.

        Parameters
        ----------
        lines : List[str]
            Lines that will contain `info` representation.
        ids : Index
            The DataFrame's column names.
        dtypes : Series
            The DataFrame's columns' dtypes.
        show_counts : bool
            If True, count of non-NA cells for each column will be appended to `lines`.
        """
        pass

    @abstractmethod
    def _non_verbose_repr(self, lines: List[str], ids: "Index") -> None:
        """
        Append short summary of columns' names to `lines`.

        Parameters
        ----------
        lines : List[str]
            Lines that will contain `info` representation.
        ids : Index
            The DataFrame's column names.
        """
        pass

    def info(self) -> None:
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
        lines = []

        lines.append(str(type(self.data)))
        lines.append(self.data.index._summary())

        ids, dtypes = self._get_ids_and_dtypes()
        col_count = len(ids)

        if col_count == 0:
            lines.append(f"Empty {type(self.data).__name__}")
            fmt.buffer_put_lines(self.buf, lines)
            return

        # hack
        max_cols = self.max_cols
        if max_cols is None:
            max_cols = get_option("display.max_info_columns", col_count + 1)

        max_rows = get_option("display.max_info_rows", len(self.data) + 1)

        if self.null_counts is None:
            show_counts = (col_count <= max_cols) and (len(self.data) < max_rows)
        else:
            show_counts = self.null_counts
        exceeds_info_cols = col_count > max_cols

        if self.verbose:
            self._verbose_repr(lines, ids, dtypes, show_counts)
        elif self.verbose is False:  # specifically set to False, not necessarily None
            self._non_verbose_repr(lines, ids)
        else:
            if exceeds_info_cols:
                self._non_verbose_repr(lines, ids)
            else:
                self._verbose_repr(lines, ids, dtypes, show_counts)

        # groupby dtype.name to collect e.g. Categorical columns
        counts = dtypes.value_counts().groupby(lambda x: x.name).sum()
        collected_dtypes = [f"{k[0]}({k[1]:d})" for k in sorted(counts.items())]
        lines.append(f"dtypes: {', '.join(collected_dtypes)}")

        if self.memory_usage:
            # append memory usage of df to display
            size_qualifier = ""
            if self.memory_usage == "deep":
                deep = True
            else:
                # size_qualifier is just a best effort; not guaranteed to catch
                # all cases (e.g., it misses categorical data even with object
                # categories)
                deep = False
                if "object" in counts or self.data.index._is_memory_usage_qualified():
                    size_qualifier = "+"
            mem_usage = self._get_mem_usage(deep=deep)
            lines.append(f"memory usage: {_sizeof_fmt(mem_usage, size_qualifier)}\n")
        fmt.buffer_put_lines(self.buf, lines)


class DataFrameInfo(BaseInfo):
    def _get_mem_usage(self, deep: bool) -> int:
        return self.data.memory_usage(index=True, deep=deep).sum()

    def _get_ids_and_dtypes(self) -> Tuple["Index", "Series"]:
        return self.data.columns, self.data.dtypes

    def _verbose_repr(
        self, lines: List[str], ids: "Index", dtypes: "Series", show_counts: bool
    ) -> None:
        col_count = len(ids)
        lines.append(f"Data columns (total {col_count} columns):")

        id_head = " # "
        column_head = "Column"
        col_space = 2

        max_col = max(len(pprint_thing(k)) for k in ids)
        len_column = len(pprint_thing(column_head))
        space = max(max_col, len_column) + col_space

        max_id = len(pprint_thing(col_count))
        len_id = len(pprint_thing(id_head))
        space_num = max(max_id, len_id) + col_space

        header = _put_str(id_head, space_num) + _put_str(column_head, space)
        if show_counts:
            counts = self.data.count()
            if col_count != len(counts):  # pragma: no cover
                raise AssertionError(
                    f"Columns must equal counts ({col_count} != {len(counts)})"
                )
            count_header = "Non-Null Count"
            len_count = len(count_header)
            non_null = " non-null"
            max_count = max(len(pprint_thing(k)) for k in counts) + len(non_null)
            space_count = max(len_count, max_count) + col_space
            count_temp = "{count}" + non_null
        else:
            count_header = ""
            space_count = len(count_header)
            len_count = space_count
            count_temp = "{count}"

        dtype_header = "Dtype"
        len_dtype = len(dtype_header)
        max_dtypes = max(len(pprint_thing(k)) for k in dtypes)
        space_dtype = max(len_dtype, max_dtypes)
        header += _put_str(count_header, space_count) + _put_str(
            dtype_header, space_dtype
        )

        lines.append(header)
        lines.append(
            _put_str("-" * len_id, space_num)
            + _put_str("-" * len_column, space)
            + _put_str("-" * len_count, space_count)
            + _put_str("-" * len_dtype, space_dtype)
        )

        for i, col in enumerate(ids):
            dtype = dtypes[i]
            col = pprint_thing(col)

            line_no = _put_str(f" {i}", space_num)
            count = ""
            if show_counts:
                count = counts[i]

            lines.append(
                line_no
                + _put_str(col, space)
                + _put_str(count_temp.format(count=count), space_count)
                + _put_str(dtype, space_dtype)
            )

    def _non_verbose_repr(self, lines: List[str], ids: "Index") -> None:
        lines.append(ids._summary(name="Columns"))
