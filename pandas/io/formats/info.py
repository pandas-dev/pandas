import sys
from typing import IO, TYPE_CHECKING, Optional, Tuple, Union

from pandas._config import get_option

from pandas._typing import Dtype, FrameOrSeries

from pandas.io.formats import format as fmt
from pandas.io.formats.printing import pprint_thing

if TYPE_CHECKING:
    from pandas.core.indexes.api import Index  # noqa: F401
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


def _get_ids_and_dtypes(data: FrameOrSeries) -> Tuple["Index", "Series"]:
    """
    Get DataFrame's columns and dtypes.

    Parameters
    ----------
    data : DataFrame
        Object that `info` was called on.

    Returns
    -------
    ids : Index
        DataFrame's columns.
    dtypes : Series
        Dtype of each of the DataFrame's columns.
    """
    ids = data.columns
    dtypes = data.dtypes
    return ids, dtypes


def info(
    data: FrameOrSeries,
    verbose: Optional[bool] = None,
    buf: Optional[IO[str]] = None,
    max_cols: Optional[int] = None,
    memory_usage: Optional[Union[bool, str]] = None,
    null_counts: Optional[bool] = None,
) -> None:
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
    if buf is None:  # pragma: no cover
        buf = sys.stdout

    lines = []

    lines.append(str(type(data)))
    lines.append(data.index._summary())

    ids, dtypes = _get_ids_and_dtypes(data)
    col_count = len(ids)

    if col_count == 0:
        lines.append(f"Empty {type(data).__name__}")
        fmt.buffer_put_lines(buf, lines)
        return

    # hack
    if max_cols is None:
        max_cols = get_option("display.max_info_columns", col_count + 1)

    max_rows = get_option("display.max_info_rows", len(data) + 1)

    if null_counts is None:
        show_counts = (col_count <= max_cols) and (len(data) < max_rows)
    else:
        show_counts = null_counts
    exceeds_info_cols = col_count > max_cols

    def _verbose_repr():
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
            counts = data.count()
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

    def _non_verbose_repr():
        lines.append(ids._summary(name="Columns"))

    def _sizeof_fmt(num, size_qualifier):
        # returns size in human readable format
        for x in ["bytes", "KB", "MB", "GB", "TB"]:
            if num < 1024.0:
                return f"{num:3.1f}{size_qualifier} {x}"
            num /= 1024.0
        return f"{num:3.1f}{size_qualifier} PB"

    if verbose:
        _verbose_repr()
    elif verbose is False:  # specifically set to False, not necessarily None
        _non_verbose_repr()
    else:
        if exceeds_info_cols:
            _non_verbose_repr()
        else:
            _verbose_repr()

    # groupby dtype.name to collect e.g. Categorical columns
    counts = dtypes.value_counts().groupby(lambda x: x.name).sum()
    collected_dtypes = [f"{k[0]}({k[1]:d})" for k in sorted(counts.items())]
    lines.append(f"dtypes: {', '.join(collected_dtypes)}")

    if memory_usage is None:
        memory_usage = get_option("display.memory_usage")
    if memory_usage:
        # append memory usage of df to display
        size_qualifier = ""
        if memory_usage == "deep":
            deep = True
        else:
            # size_qualifier is just a best effort; not guaranteed to catch
            # all cases (e.g., it misses categorical data even with object
            # categories)
            deep = False
            if "object" in counts or data.index._is_memory_usage_qualified():
                size_qualifier = "+"
        mem_usage = data.memory_usage(index=True, deep=deep).sum()
        lines.append(f"memory usage: {_sizeof_fmt(mem_usage, size_qualifier)}\n")
    fmt.buffer_put_lines(buf, lines)
