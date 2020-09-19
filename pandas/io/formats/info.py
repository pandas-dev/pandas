from abc import ABCMeta, abstractmethod
import sys
from typing import IO, TYPE_CHECKING, List, NamedTuple, Optional, Tuple, Union, cast

from pandas._config import get_option

from pandas._typing import Dtype, FrameOrSeries

from pandas.core.indexes.api import Index

from pandas.io.formats import format as fmt
from pandas.io.formats.printing import pprint_thing

if TYPE_CHECKING:
    from pandas.core.series import Series  # noqa: F401


class CountConfigs(NamedTuple):
    """
    Configs with which to display counts.

    Attributes
    ----------
    counts : Series
        Non-null count of Series (or of each column of DataFrame).
    count_header : str
        Header that will be printed out above non-null counts in output.
    space_count : int
        Number of spaces that count_header should occupy
        (including space before `dtypes` column).
    len_count : int
        Length of count header.
    count_temp : str
        String that can be formatted to include non-null count.
    """

    counts: "Series"
    count_header: str
    space_count: int
    len_count: int
    count_temp: str


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


def _get_count_configs(
    counts: "Series", col_space: int, show_counts: bool, col_count: Optional[int] = None
) -> CountConfigs:
    """
    Get configs for displaying counts, depending on the value of `show_counts`.

    Parameters
    ----------
    counts : Series
        Non-null count of Series (or of each column of DataFrame).
    col_space : int
        How many space to leave between non-null count and dtype columns.
    show_counts : bool
        Whether to display non-null counts.
    col_count : int, optional
        Number of columns in DataFrame.

    Returns
    -------
    CountConfigs
    """
    if show_counts:
        if col_count is not None and col_count != len(counts):  # pragma: no cover
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
    return CountConfigs(counts, count_header, space_count, len_count, count_temp)


def _display_counts_and_dtypes(
    lines: List[str],
    ids: "Index",
    dtypes: "Series",
    show_counts: bool,
    count_configs: CountConfigs,
    space_dtype: int,
    space: int = 0,
    space_num: int = 0,
) -> None:
    """
    Append count and dtype of Series (or of each column of Frame) to `lines`.

    Parameters
    ----------
    lines : List[str]
        At this stage, this contains the main header and the info table headers.
    ids : Index
        Series name (or names of DataFrame columns).
    dtypes : Series
        Series dtype (or dtypes of DataFrame columns).
    show_counts : bool
        Whether to show non-null counts.
    count_configs: CountConfigs
        Configs with which to display counts.
    space_dtype : int
        Number of spaces that `dtypes` column should occupy.
    space : int = 0
        Number of spaces that `Column` header should occupy
        (including space before `non-null count` column).
    space_num : int = 0
        Number of spaces that ` # ` header should occupy (including space
        before `Column` column), only applicable for `DataFrame.info`.
    """
    for i, col in enumerate(ids):
        dtype = dtypes[i]
        col = pprint_thing(col)

        line_no = _put_str(f" {i}", space_num)
        count = ""
        if show_counts:
            count = count_configs.counts[i]

        lines.append(
            line_no
            + _put_str(col, space)
            + _put_str(
                count_configs.count_temp.format(count=count), count_configs.space_count
            )
            + _put_str(dtype, space_dtype)
        )


def _get_header_and_spaces(
    dtypes: "Series", space_count: int, count_header: str, header: str = ""
) -> Tuple[int, str, int]:
    """
    Append extra columns (count and type) to header, if applicable.

    Parameters
    ----------
    dtypes : Series
        Series dtype (or dtypes of DataFrame columns).
    space_count : int
        Number of spaces that count_header should occupy
        (including space before `dtypes` column).
    count_header : str
        Header that will be printed out above non-null counts in output.
    header : str
        Current header.

    Returns
    -------
    space_dtype : int
        Number of spaces that `dtypes` column should occupy.
    header : str
        Header with extra columns (count and type) appended.
    len_dtype : int
        Length of dtype header.
    """
    dtype_header = "Dtype"
    len_dtype = len(dtype_header)
    max_dtypes = max(len(pprint_thing(k)) for k in dtypes)
    space_dtype = max(len_dtype, max_dtypes)
    header += _put_str(count_header, space_count) + _put_str(dtype_header, space_dtype)
    return space_dtype, header, len_dtype


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
        counts = self.data.count()
        count_configs = _get_count_configs(counts, col_space, show_counts, col_count)

        space_dtype, header, len_dtype = _get_header_and_spaces(
            dtypes, count_configs.space_count, count_configs.count_header, header
        )

        lines.append(header)
        lines.append(
            _put_str("-" * len_id, space_num)
            + _put_str("-" * len_column, space)
            + _put_str("-" * count_configs.len_count, count_configs.space_count)
            + _put_str("-" * len_dtype, space_dtype,)
        )

        _display_counts_and_dtypes(
            lines,
            ids,
            dtypes,
            show_counts,
            count_configs,
            space_dtype,
            space,
            space_num,
        )

    def _non_verbose_repr(self, lines: List[str], ids: "Index") -> None:
        lines.append(ids._summary(name="Columns"))


class SeriesInfo(BaseInfo):
    def _get_mem_usage(self, deep: bool) -> int:
        return self.data.memory_usage(index=True, deep=deep)

    def _get_ids_and_dtypes(self) -> Tuple["Index", "Series"]:
        ids = Index([self.data.name])
        dtypes = cast("Series", self.data._constructor(self.data.dtypes))
        return ids, dtypes

    def _verbose_repr(
        self, lines: List[str], ids: "Index", dtypes: "Series", show_counts: bool
    ) -> None:
        lines.append(f"Series name: {self.data.name}")

        id_space = 2

        counts = cast("Series", self.data._constructor(self.data.count()))
        count_configs = _get_count_configs(counts, id_space, show_counts)

        space_dtype, header, len_dtype = _get_header_and_spaces(
            dtypes, count_configs.space_count, count_configs.count_header
        )

        lines.append(header)
        lines.append(
            _put_str("-" * count_configs.len_count, count_configs.space_count)
            + _put_str("-" * len_dtype, space_dtype)
        )

        _display_counts_and_dtypes(
            lines, ids, dtypes, show_counts, count_configs, space_dtype,
        )

    def _non_verbose_repr(self, lines: List[str], ids: "Index") -> None:
        pass
