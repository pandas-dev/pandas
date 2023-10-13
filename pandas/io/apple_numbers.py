"""Apple Numbers format"""
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

from pandas._libs import lib
from pandas.compat._optional import import_optional_dependency

if TYPE_CHECKING:
    from collections.abc import (
        Hashable,
        Iterable,
        Sequence,
    )

    from pandas._typing import (
        DtypeBackend,
        FilePath,
        IndexLabel,
        ReadBuffer,
        WriteBuffer,
    )

    from pandas.core.api import DataFrame


def read_apple_numbers(
    io: FilePath | ReadBuffer[str],
    *,
    tables: dict | None = None,
    header: int | Sequence[int] | None = None,
    index_col: int | Sequence[int] | None = None,
    skiprows: int | Sequence[int] | slice | None = None,
    parse_dates: bool = False,
    thousands: str | None = ",",
    decimal: str = ".",
    converters: dict | None = None,
    na_values: Iterable[object] | None = None,
    keep_default_na: bool = True,
    dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
) -> DataFrame | list[DataFrame]:
    r"""
    Read Numbers Document tables ``DataFrame`` objects.

    Supports an option to read a single table from one sheet or a multiple tables
    returned as a list of ``DataFrame`` objects

    Parameters
    ----------
    io : path object, or file-like object
        Path object (implementing ``os.PathLike[str]``), or file-like
        object implementing a string ``read()`` function.

    sheets : dict, or None, default None
        The sheets and tables to load. Strings are used for sheet and table names.
        Integers are used in zero-indexed sheet and table positions.  Lists of
        strings/integers are used to request multiple sheets or tables. Specify
        None to get all tables from all sheets.

        Example cases:

        * ``{0: 0}``: the first table of the first sheet
        * ``{0: [0, "Table 4"]}``: the first table and "Table 4" from the
        first sheet
        * ``{"Sheet 2: [1, 2]}``: the second and third table from the sheet
        named "Sheet 2"
        * ``{0: 0, "Sheet 2": ["Table 1", "Table 2"]}``: the first table of
        the first sheet, and tables named "Table 1" and "Table 2" from the
        sheet named "Sheet 2"
        * None: All tables from all sheets.

    header : int or list-like, optional
        The row (or list of rows for a :class:`~pandas.MultiIndex`) to use to
        make the columns headers.

    index_col : int or list-like, optional
        The column (or list of columns) to use to create the index.

    skiprows : int, list-like or slice, optional
        Number of rows to skip after parsing the column integer. 0-based. If a
        sequence of integers or a slice is given, will skip the rows indexed by
        that sequence.  Note that a single element sequence means 'skip the nth
        row' whereas an integer means 'skip n rows'.

    parse_dates : bool, optional
        See :func:`~read_csv` for more details.

    thousands : str, optional
        Separator to use to parse thousands. Defaults to ``','``.

    decimal : str, default '.'
        Character to recognize as decimal point (e.g. use ',' for European
        data).

    converters : dict, default None
        Dict of functions for converting values in certain columns. Keys can
        either be integers or column labels, values are functions that take one
        input argument, the cell (not column) content, and return the
        transformed content.

    na_values : iterable, default None
        Custom NA values.

    keep_default_na : bool, default True
        If na_values are specified and keep_default_na is False the default NaN
        values are overridden, otherwise they're appended to.

    dtype_backend : {{'numpy_nullable', 'pyarrow'}}, default 'numpy_nullable'
        Back-end data type applied to the resultant :class:`DataFrame`
        (still experimental). Behaviour is as follows:

        * ``"numpy_nullable"``: returns nullable-dtype-backed :class:`DataFrame`
          (default).
        * ``"pyarrow"``: returns pyarrow-backed nullable :class:`ArrowDtype`
          DataFrame.
    """
    import_optional_dependency("numbers_parser")

    raise NotImplementedError("Apple Numbers support is in development")


def to_apple_numbers(
    self,
    target: FilePath | WriteBuffer[str] | object,
    sheet_name: str = "Sheet 1",
    table_name: str = "Table 1",
    na_rep: str = "",
    columns: Sequence[Hashable] | None = None,
    header: Sequence[Hashable] | bool = True,
    index: bool = True,
    index_label: IndexLabel | None = None,
    startrow: int = 0,
    startcol: int = 0,
    merge_cells: bool = True,
    engine_kwargs: dict[str, Any] | None = None,
) -> None:
    """
    Write a DataFrame to an Apple Numbers document.

    To write a single object to a Numbers document, it is only necessary to
    specify the target filename. To write mltiple tables, you must first
    create a `numbers_parser.Document` with a target file name, and then
    specify a sheet and table in the file to write to.

    Multiple tables may be written to by specifying unique `sheet_name`
    and `table_name`.

    Parameters
    ----------
    target : path-like, file-like, or Document object
        File path or existing Document.
    sheet_name : str, default 'Sheet 1'
        Name of sheet which will contain DataFrame.
    table_name : str, default 'Table 1'
        Name of table which will contain DataFrame.
    na_rep : str, default ''
        Missing data representation.
    columns : sequence or list of str, optional
        Columns to write.
    header : bool or list of str, default True
        Write out the column names. If a list of string is given it is
        assumed to be aliases for the column names.
    index : bool, default True
        Write row names (index).
    index_label : str or sequence, optional
        Column label for index column(s) if desired. If not specified, and
        `header` and `index` are True, then the index names are used. A
        sequence should be given if the DataFrame uses MultiIndex.
    startrow : int, default 0
        Upper left cell row to dump data frame.
    startcol : int, default 0
        Upper left cell column to dump data frame.
    merge_cells : bool, default True
        Write MultiIndex and Hierarchical Rows as merged cells.

    engine_kwargs : dict, optional
        Arbitrary keyword arguments passed to Numbers writer

    Documents can be resaved, but once a `DataFrame` is written to the document,
    that table cannot be overwritten with another `DataFrame`. Other objects
    can be written to new sheets and tables.
    """
    import_optional_dependency("numbers_parser")

    raise NotImplementedError("Apple Numbers support is in development")
