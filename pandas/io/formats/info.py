import sys

from pandas._config import get_option

from pandas.io.formats import format as fmt
from pandas.io.formats.printing import pprint_thing


def _put_str(s, space):
    return str(s)[:space].ljust(space)


def info(
    data, verbose=None, buf=None, max_cols=None, memory_usage=None, null_counts=None
) -> None:
    """
    Print a concise summary of a DataFrame.

    This method prints information about a DataFrame including
    the index dtype and column dtypes, non-null values and memory usage.

    Parameters
    ----------
    data : DataFrame
        DataFrame to print information about.
    verbose : bool, optional
        Whether to print the full summary. By default, the setting in
        ``pandas.options.display.max_info_columns`` is followed.
    buf : writable buffer, defaults to sys.stdout
        Where to send the output. By default, the output is printed to
        sys.stdout. Pass a writable buffer if you need to further process
        the output.
    max_cols : int, optional
        When to switch from the verbose to the truncated output. If the
        DataFrame has more than `max_cols` columns, the truncated output
        is used. By default, the setting in
        ``pandas.options.display.max_info_columns`` is used.
    memory_usage : bool, str, optional
        Specifies whether total memory usage of the DataFrame
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
        only if the frame is smaller than
        ``pandas.options.display.max_info_rows`` and
        ``pandas.options.display.max_info_columns``. A value of True always
        shows the counts, and False never shows the counts.

    Returns
    -------
    None
        This method prints a summary of a DataFrame and returns None.

    See Also
    --------
    DataFrame.describe: Generate descriptive statistics of DataFrame
        columns.
    DataFrame.memory_usage: Memory usage of DataFrame columns.

    Examples
    --------
    >>> int_values = [1, 2, 3, 4, 5]
    >>> text_values = ['alpha', 'beta', 'gamma', 'delta', 'epsilon']
    >>> float_values = [0.0, 0.25, 0.5, 0.75, 1.0]
    >>> df = pd.DataFrame({"int_col": int_values, "text_col": text_values,
    ...                   "float_col": float_values})
    >>> df
        int_col text_col  float_col
    0        1    alpha       0.00
    1        2     beta       0.25
    2        3    gamma       0.50
    3        4    delta       0.75
    4        5  epsilon       1.00

    Prints information of all columns:

    >>> df.info(verbose=True)
    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 5 entries, 0 to 4
    Data columns (total 3 columns):
        #   Column     Non-Null Count  Dtype
    ---  ------     --------------  -----
        0   int_col    5 non-null      int64
        1   text_col   5 non-null      object
        2   float_col  5 non-null      float64
    dtypes: float64(1), int64(1), object(1)
    memory usage: 248.0+ bytes

    Prints a summary of columns count and its dtypes but not per column
    information:

    >>> df.info(verbose=False)
    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 5 entries, 0 to 4
    Columns: 3 entries, int_col to float_col
    dtypes: float64(1), int64(1), object(1)
    memory usage: 248.0+ bytes

    Pipe output of DataFrame.info to buffer instead of sys.stdout, get
    buffer content and writes to a text file:

    >>> import io
    >>> buffer = io.StringIO()
    >>> df.info(buf=buffer)
    >>> s = buffer.getvalue()
    >>> with open("df_info.txt", "w",
    ...           encoding="utf-8") as f:  # doctest: +SKIP
    ...     f.write(s)
    260

    The `memory_usage` parameter allows deep introspection mode, specially
    useful for big DataFrames and fine-tune memory optimization:

    >>> random_strings_array = np.random.choice(['a', 'b', 'c'], 10 ** 6)
    >>> df = pd.DataFrame({
    ...     'column_1': np.random.choice(['a', 'b', 'c'], 10 ** 6),
    ...     'column_2': np.random.choice(['a', 'b', 'c'], 10 ** 6),
    ...     'column_3': np.random.choice(['a', 'b', 'c'], 10 ** 6)
    ... })
    >>> df.info()
    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 1000000 entries, 0 to 999999
    Data columns (total 3 columns):
        #   Column    Non-Null Count    Dtype
    ---  ------    --------------    -----
        0   column_1  1000000 non-null  object
        1   column_2  1000000 non-null  object
        2   column_3  1000000 non-null  object
    dtypes: object(3)
    memory usage: 22.9+ MB

    >>> df.info(memory_usage='deep')
    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 1000000 entries, 0 to 999999
    Data columns (total 3 columns):
        #   Column    Non-Null Count    Dtype
    ---  ------    --------------    -----
        0   column_1  1000000 non-null  object
        1   column_2  1000000 non-null  object
        2   column_3  1000000 non-null  object
    dtypes: object(3)
    memory usage: 188.8 MB
    """
    if buf is None:  # pragma: no cover
        buf = sys.stdout

    lines = []

    lines.append(str(type(data)))
    lines.append(data.index._summary())

    if len(data.columns) == 0:
        lines.append(f"Empty {type(data).__name__}")
        fmt.buffer_put_lines(buf, lines)
        return

    cols = data.columns
    col_count = len(data.columns)

    # hack
    if max_cols is None:
        max_cols = get_option("display.max_info_columns", len(data.columns) + 1)

    max_rows = get_option("display.max_info_rows", len(data) + 1)

    if null_counts is None:
        show_counts = (col_count <= max_cols) and (len(data) < max_rows)
    else:
        show_counts = null_counts
    exceeds_info_cols = col_count > max_cols

    def _verbose_repr():
        lines.append(f"Data columns (total {len(data.columns)} columns):")

        id_head = " # "
        column_head = "Column"
        col_space = 2

        max_col = max(len(pprint_thing(k)) for k in cols)
        len_column = len(pprint_thing(column_head))
        space = max(max_col, len_column) + col_space

        max_id = len(pprint_thing(col_count))
        len_id = len(pprint_thing(id_head))
        space_num = max(max_id, len_id) + col_space

        header = _put_str(id_head, space_num) + _put_str(column_head, space)
        if show_counts:
            counts = data.count()
            if len(cols) != len(counts):  # pragma: no cover
                raise AssertionError(
                    f"Columns must equal counts ({len(cols)} != {len(counts)})"
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
        max_dtypes = max(len(pprint_thing(k)) for k in data.dtypes)
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

        for i, col in enumerate(data.columns):
            dtype = data.dtypes.iloc[i]
            col = pprint_thing(col)

            line_no = _put_str(f" {i}", space_num)
            count = ""
            if show_counts:
                count = counts.iloc[i]

            lines.append(
                line_no
                + _put_str(col, space)
                + _put_str(count_temp.format(count=count), space_count)
                + _put_str(dtype, space_dtype)
            )

    def _non_verbose_repr():
        lines.append(data.columns._summary(name="Columns"))

    def _sizeof_fmt(num, size_qualifier):
        # returns size in human readable format
        for x in ["bytes", "KB", "MB", "GB", "TB"]:
            if num < 1024.0:
                return f"{num:3.1f}{size_qualifier} {x}"
            num /= 1024.0
        return f"{num:3.1f}{size_qualifier} PB"

    if verbose:
        _verbose_repr()
    elif verbose is False:  # specifically set to False, not nesc None
        _non_verbose_repr()
    else:
        if exceeds_info_cols:
            _non_verbose_repr()
        else:
            _verbose_repr()

    counts = data._data.get_dtype_counts()
    dtypes = [f"{k[0]}({k[1]:d})" for k in sorted(counts.items())]
    lines.append(f"dtypes: {', '.join(dtypes)}")

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
