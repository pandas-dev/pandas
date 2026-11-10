"""
Shared building blocks for the pandas SQL backends: result wrapping,
date parsing, and the ``SQLTable`` / ``PandasSQL`` base classes.
"""

from __future__ import annotations

from abc import (
    ABC,
    abstractmethod,
)
from datetime import (
    date,
    datetime,
)
from functools import partial
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Self,
    cast,
)
import warnings

import numpy as np

from pandas._config import using_string_dtype

from pandas._libs import lib
from pandas.compat._optional import (
    import_optional_dependency,
)
from pandas.errors import (
    Pandas4Warning,
)
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    is_dict_like,
    is_list_like,
    is_object_dtype,
    is_string_dtype,
)
from pandas.core.dtypes.dtypes import DatetimeTZDtype
from pandas.core.dtypes.missing import isna

from pandas.core.api import (
    DataFrame,
    Series,
)
from pandas.core.arrays import ArrowExtensionArray
from pandas.core.arrays.string_ import StringDtype
from pandas.core.base import PandasObject
import pandas.core.common as com
from pandas.core.internals.construction import convert_object_array
from pandas.core.tools.datetimes import (
    stringify_numeric_column,
    to_datetime,
)

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Generator,
        Hashable,
        Iterator,
        Mapping,
    )
    from contextlib import (
        ExitStack,
    )

    from sqlalchemy.sql.expression import (
        Select,
        TextClause,
    )

    from pandas._typing import (
        DtypeArg,
        DtypeBackend,
    )

    from pandas import Index


# -----------------------------------------------------------------------------
# -- Helper functions


def _process_parse_dates_argument(parse_dates):
    """Process parse_dates argument for read_sql functions"""
    # handle non-list entries for parse_dates gracefully
    if parse_dates is True or parse_dates is None or parse_dates is False:
        parse_dates = []

    elif not hasattr(parse_dates, "__iter__"):
        parse_dates = [parse_dates]
    return parse_dates


def _warn_and_stringify_numeric_column(col):
    # GH#55663 read_sql historically cast integer/float columns to strings so a
    # ``format`` could be applied to them. That is deprecated in favor of the
    # user casting explicitly, which lets us drop this shim in a future version.
    warnings.warn(
        "Passing a format for integer or float columns to read_sql via "
        "parse_dates is deprecated. Cast the column to string and apply "
        "to_datetime after reading, e.g. "
        "to_datetime(df[col].astype('string'), format=...).",
        Pandas4Warning,
        stacklevel=find_stack_level(),
    )
    return stringify_numeric_column(col)


def _handle_date_column(
    col, utc: bool = False, format: str | dict[str, Any] | None = None
):
    if isinstance(format, dict):
        # GH35185 Allow custom error values in parse_dates argument of
        # read_sql like functions.
        # Format can take on custom to_datetime argument values such as
        # {"errors": "coerce"} or {"dayfirst": True}
        if format.get("format") is not None and col.dtype.kind in "iuf":
            col = _warn_and_stringify_numeric_column(col)
        return to_datetime(col, **format)
    else:
        # Allow passing of formatting string for integers
        # GH17855
        if format is None and (
            issubclass(col.dtype.type, np.floating)
            or issubclass(col.dtype.type, np.integer)
        ):
            format = "s"
        if format in ["D", "d", "h", "m", "s", "ms", "us", "ns"]:
            return to_datetime(col, errors="coerce", unit=format, utc=utc)
        elif isinstance(col.dtype, DatetimeTZDtype):
            # coerce to UTC timezone
            # GH11216
            return to_datetime(col, utc=True)
        else:
            if format is not None and col.dtype.kind in "iuf":
                col = _warn_and_stringify_numeric_column(col)
            return to_datetime(col, errors="coerce", format=format, utc=utc)


def parse_date_columns(data_frame: DataFrame, parse_dates) -> DataFrame:
    """
    Force non-datetime columns to be read as such.
    Supports both string formatted and integer timestamp columns.
    """
    parse_dates = _process_parse_dates_argument(parse_dates)

    # we want to coerce datetime64_tz dtypes for now to UTC
    # we could in theory do a 'nice' conversion from a FixedOffset tz
    # GH11216
    for i, (col_name, df_col) in enumerate(data_frame.items()):
        if isinstance(df_col.dtype, DatetimeTZDtype) or col_name in parse_dates:
            try:
                fmt = parse_dates[col_name]
            except (KeyError, TypeError):
                fmt = None
            data_frame.isetitem(i, _handle_date_column(df_col, format=fmt))

    return data_frame


def _convert_arrays_to_dataframe(
    data,
    columns,
    coerce_float: bool = True,
    dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
) -> DataFrame:
    if data and is_dict_like(data[0]) and not isinstance(data[0], (list, tuple)):
        # DBAPI cursors returning dict rows, e.g. pymysql's DictCursor. Use
        # values() rather than lookup by column name, since such cursors may
        # rename duplicate columns (GH#53028). Skip list-based rows like
        # psycopg2's DictRow, whose values() drops duplicate columns.
        data = [tuple(row.values()) for row in data]
    content = lib.to_object_array_tuples(data)
    idx_len = content.shape[0]
    arrays = convert_object_array(
        list(content.T),
        dtype=None,
        coerce_float=coerce_float,
        dtype_backend=dtype_backend,
    )
    if dtype_backend == "pyarrow":
        pa = import_optional_dependency("pyarrow")

        result_arrays = []
        for arr in arrays:
            pa_array = pa.array(arr, from_pandas=True)
            if arr.dtype == "string":
                # TODO: Arrow still infers strings arrays as regular strings instead
                # of large_string, which is what we preserver everywhere else for
                # dtype_backend="pyarrow". We may want to reconsider this
                pa_array = pa_array.cast(pa.string())
            result_arrays.append(ArrowExtensionArray(pa_array))
        arrays = result_arrays  # type: ignore[assignment]
    if arrays:
        return DataFrame._from_arrays(
            arrays, columns=columns, index=range(idx_len), verify_integrity=False
        )
    else:
        return DataFrame(columns=columns)


def wrap_result(
    data,
    columns,
    index_col=None,
    coerce_float: bool = True,
    parse_dates=None,
    dtype: DtypeArg | None = None,
    dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
) -> DataFrame:
    """Wrap result set of a SQLAlchemy query in a DataFrame."""
    frame = _convert_arrays_to_dataframe(data, columns, coerce_float, dtype_backend)

    if dtype:
        frame = frame.astype(dtype)

    frame = parse_date_columns(frame, parse_dates)

    if index_col is not None:
        frame = frame.set_index(index_col)

    return frame


def get_unicode_name(name: object) -> str:
    try:
        uname = str(name).encode("utf-8", "strict").decode("utf-8")
    except UnicodeError as err:
        raise ValueError(f"Cannot convert identifier to UTF-8: '{name}'") from err
    return uname


class SQLTable(PandasObject):
    """
    For mapping Pandas tables to SQL tables.
    Uses fact that table is reflected by SQLAlchemy to
    do better type conversions.
    Also holds various flags needed to avoid having to
    pass them between functions all the time.
    """

    # TODO: support for multiIndex

    def __init__(
        self,
        name: str,
        pandas_sql_engine,
        frame=None,
        index: bool | str | list[str] | None = True,
        if_exists: Literal["fail", "replace", "append", "delete_rows"] = "fail",
        prefix: str = "pandas",
        index_label=None,
        schema=None,
        keys=None,
        dtype: DtypeArg | None = None,
    ) -> None:
        self.name = name
        self.pd_sql = pandas_sql_engine
        self.prefix = prefix
        self.frame = frame
        self.index = self._index_name(index, index_label)
        self.schema = schema
        self.if_exists = if_exists
        self.keys = keys
        self.dtype = dtype

        if frame is not None:
            # We want to initialize based on a dataframe
            self.table = self._create_table_setup()
        else:
            # no data provided, read-only mode
            self.table = self.pd_sql.get_table(self.name, self.schema)

        if self.table is None:
            raise ValueError(f"Could not init table '{name}'")

        if not len(self.name):
            raise ValueError("Empty table name specified")

    def exists(self):
        return self.pd_sql.has_table(self.name, self.schema)

    def sql_schema(self) -> str:
        from sqlalchemy.schema import CreateTable

        return str(CreateTable(self.table).compile(self.pd_sql.con))

    def _execute_create(self) -> None:
        # Inserting table into database, add to MetaData object
        self.table = self.table.to_metadata(self.pd_sql.meta)
        with self.pd_sql.run_transaction():
            self.table.create(bind=self.pd_sql.con)

    def create(self) -> bool:
        """
        Create the table if needed; return whether a table was (re)created.
        """
        if self.exists():
            if self.if_exists == "fail":
                raise ValueError(f"Table '{self.name}' already exists.")
            elif self.if_exists == "replace":
                self.pd_sql.drop_table(self.name, self.schema)
                self._execute_create()
                return True
            elif self.if_exists == "append":
                pass
            elif self.if_exists == "delete_rows":
                self.pd_sql.delete_rows(self.name, self.schema)
            else:
                raise ValueError(f"'{self.if_exists}' is not valid for if_exists")
            return False
        else:
            self._execute_create()
            return True

    def _execute_insert(self, conn, keys: list[str], data_iter) -> int:
        """
        Execute SQL statement inserting data

        Parameters
        ----------
        conn : sqlalchemy.engine.Engine or sqlalchemy.engine.Connection
        keys : list of str
           Column names
        data_iter : generator of list
           Each item contains a list of values to be inserted
        """
        data = [dict(zip(keys, row, strict=True)) for row in data_iter]
        result = self.pd_sql.execute(self.table.insert(), data)
        return result.rowcount

    def _execute_insert_multi(self, conn, keys: list[str], data_iter) -> int:
        """
        Alternative to _execute_insert for DBs support multi-value INSERT.

        Note: multi-value insert is usually faster for analytics DBs
        and tables containing a few columns
        but performance degrades quickly with increase of columns.

        """

        from sqlalchemy import insert

        data = [dict(zip(keys, row, strict=True)) for row in data_iter]
        stmt = insert(self.table).values(data)
        result = self.pd_sql.execute(stmt)
        return result.rowcount

    def insert_data(self) -> tuple[list[str], list[np.ndarray]]:
        if self.index is not None:
            temp = self.frame.copy(deep=False)
            temp.index.names = self.index
            try:
                temp = temp.reset_index()
            except ValueError as err:
                raise ValueError(f"duplicate name in index/columns: {err}") from err
        else:
            temp = self.frame

        column_names = list(map(str, temp.columns))
        ncols = len(column_names)
        # this just pre-allocates the list: None's will be replaced with ndarrays
        # error: List item 0 has incompatible type "None"; expected "ndarray"
        data_list: list[np.ndarray] = [None] * ncols  # type: ignore[list-item]

        for i, (_, ser) in enumerate(temp.items()):
            if ser.dtype.kind == "M":
                if isinstance(ser._values, ArrowExtensionArray):
                    import pyarrow as pa

                    if pa.types.is_date(ser.dtype.pyarrow_dtype):
                        # GH#53854 to_pydatetime not supported for pyarrow date dtypes
                        d = ser._values.to_numpy(dtype=object)
                    else:
                        d = ser.dt.to_pydatetime()._values
                else:
                    d = ser._values.to_pydatetime()
            elif ser.dtype.kind == "m":
                vals = ser._values
                if isinstance(vals, ArrowExtensionArray):
                    vals = vals.to_numpy(dtype=np.dtype("m8[ns]"))
                else:
                    vals = vals.as_unit("ns")
                # store as integers, see GH#6921, GH#7076
                d = vals.view("i8").astype(object)
            else:
                d = ser._values.astype(object)

            assert isinstance(d, np.ndarray), type(d)

            if ser._can_hold_na:
                # Note: this will miss timedeltas since they are converted to int
                mask = isna(d)
                d[mask] = None

            data_list[i] = d

        return column_names, data_list

    def insert(
        self,
        chunksize: int | None = None,
        method: Literal["multi"] | Callable[..., Any] | None = None,
    ) -> int | None:
        # set insert method
        if method is None:
            exec_insert = self._execute_insert
        elif method == "multi":
            exec_insert = self._execute_insert_multi
        elif callable(method):
            exec_insert = partial(method, self)
        else:
            raise ValueError(f"Invalid parameter `method`: {method}")

        keys, data_list = self.insert_data()

        nrows = len(self.frame)

        if nrows == 0:
            return 0

        if chunksize is None:
            chunksize = nrows
        elif chunksize == 0:
            raise ValueError("chunksize argument should be non-zero")

        chunks = (nrows // chunksize) + 1
        total_inserted = None
        with self.pd_sql.run_transaction() as conn:
            for i in range(chunks):
                start_i = i * chunksize
                end_i = min((i + 1) * chunksize, nrows)
                if start_i >= end_i:
                    break

                chunk_iter = zip(
                    *(arr[start_i:end_i] for arr in data_list), strict=True
                )
                num_inserted = exec_insert(conn, keys, chunk_iter)
                # GH 46891
                if num_inserted is not None:
                    if total_inserted is None:
                        total_inserted = num_inserted
                    else:
                        total_inserted += num_inserted
        return total_inserted

    def _query_iterator(
        self,
        result,
        exit_stack: ExitStack,
        chunksize: int | None,
        columns,
        coerce_float: bool = True,
        parse_dates=None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> Generator[DataFrame]:
        """Return generator through chunked result set."""
        has_read_data = False
        with exit_stack:
            while True:
                data = result.fetchmany(chunksize)
                if not data:
                    if not has_read_data:
                        yield DataFrame.from_records(
                            [], columns=columns, coerce_float=coerce_float
                        )
                    break

                has_read_data = True
                self.frame = _convert_arrays_to_dataframe(
                    data, columns, coerce_float, dtype_backend
                )

                self._harmonize_columns(
                    parse_dates=parse_dates, dtype_backend=dtype_backend
                )

                if self.index is not None:
                    self.frame = self.frame.set_index(self.index)

                yield self.frame

    def read(
        self,
        exit_stack: ExitStack,
        coerce_float: bool = True,
        parse_dates=None,
        columns=None,
        chunksize: int | None = None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> DataFrame | Iterator[DataFrame]:
        from sqlalchemy import select

        if columns is not None and len(columns) > 0:
            cols = [self.table.c[n] for n in columns]
            if self.index is not None:
                for idx in self.index[::-1]:
                    cols.insert(0, self.table.c[idx])
            sql_select = select(*cols)
        else:
            sql_select = select(self.table)
        result = self.pd_sql.execute(sql_select)
        # SQLAlchemy keys are str subclasses that need it to unpickle (GH#52816)
        column_names = [str(col) for col in result.keys()]

        if chunksize is not None:
            return self._query_iterator(
                result,
                exit_stack,
                chunksize,
                column_names,
                coerce_float=coerce_float,
                parse_dates=parse_dates,
                dtype_backend=dtype_backend,
            )
        else:
            data = result.fetchall()
            self.frame = _convert_arrays_to_dataframe(
                data, column_names, coerce_float, dtype_backend
            )

            self._harmonize_columns(
                parse_dates=parse_dates, dtype_backend=dtype_backend
            )

            if self.index is not None:
                self.frame = self.frame.set_index(self.index)

            return self.frame

    def _index_name(self, index, index_label):
        # for writing: index=True to include index in sql table
        if index is True:
            nlevels = self.frame.index.nlevels
            # if index_label is specified, set this as index name(s)
            if index_label is not None:
                if not isinstance(index_label, list):
                    index_label = [index_label]
                if len(index_label) != nlevels:
                    raise ValueError(
                        "Length of 'index_label' should match number of "
                        f"levels, which is {nlevels}"
                    )
                return index_label
            # return the used column labels for the index columns
            if (
                nlevels == 1
                and "index" not in self.frame.columns
                and self.frame.index.name is None
            ):
                return ["index"]
            else:
                return com.fill_missing_names(self.frame.index.names)

        # for reading: index=(list of) string to specify column to set as index
        elif isinstance(index, str):
            return [index]
        elif isinstance(index, list):
            return index
        else:
            return None

    def _get_column_names_and_types(self, dtype_mapper):
        column_names_and_types = []
        if self.index is not None:
            for i, idx_label in enumerate(self.index):
                idx_type = dtype_mapper(self.frame.index._get_level_values(i))
                column_names_and_types.append((str(idx_label), idx_type, True))

        column_names_and_types += [
            (str(self.frame.columns[i]), dtype_mapper(self.frame.iloc[:, i]), False)
            for i in range(len(self.frame.columns))
        ]

        return column_names_and_types

    def _create_table_setup(self):
        from sqlalchemy import (
            Column,
            PrimaryKeyConstraint,
            Table,
        )
        from sqlalchemy.schema import MetaData

        column_names_and_types = self._get_column_names_and_types(self._sqlalchemy_type)

        columns: list[Any] = [
            Column(name, typ, index=is_index)
            for name, typ, is_index in column_names_and_types
        ]

        if self.keys is not None:
            if not is_list_like(self.keys):
                keys = [self.keys]
            else:
                keys = self.keys
            pkc = PrimaryKeyConstraint(*keys, name=self.name + "_pk")
            columns.append(pkc)

        schema = self.schema or self.pd_sql.meta.schema

        # At this point, attach to new metadata, only attach to self.meta
        # once table is created.
        meta = MetaData()
        return Table(self.name, meta, *columns, schema=schema)

    def _harmonize_columns(
        self,
        parse_dates=None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> None:
        """
        Make the DataFrame's column types align with the SQL table
        column types.
        Need to work around limited NA value support. Floats are always
        fine, ints must always be floats if there are Null values.
        Booleans are hard because converting bool column with None replaces
        all Nones with false. Therefore only convert bool if there are no
        NA values.
        Datetimes should already be converted to np.datetime64 if supported,
        but here we also force conversion if required.
        """
        parse_dates = _process_parse_dates_argument(parse_dates)

        assert self.frame is not None  # caller always sets frame before this
        for sql_col in self.table.columns:
            col_name = sql_col.name
            try:
                df_col = self.frame[col_name]

                # Handle date parsing upfront; don't try to convert columns
                # twice
                if col_name in parse_dates:
                    try:
                        fmt = parse_dates[col_name]
                    except TypeError:
                        fmt = None
                    self.frame[col_name] = _handle_date_column(df_col, format=fmt)
                    continue

                # the type the dataframe column should have
                col_type = self._get_dtype(sql_col.type)

                if (
                    col_type is datetime
                    or col_type is date
                    or col_type is DatetimeTZDtype
                ):
                    # Convert tz-aware Datetime SQL columns to UTC
                    utc = col_type is DatetimeTZDtype
                    self.frame[col_name] = _handle_date_column(df_col, utc=utc)
                elif dtype_backend == "numpy" and col_type is float:
                    # floats support NA, can always convert!
                    self.frame[col_name] = df_col.astype(col_type)
                elif (
                    using_string_dtype()
                    and is_string_dtype(col_type)
                    and is_object_dtype(self.frame[col_name])
                ):
                    self.frame[col_name] = df_col.astype(col_type)
                elif dtype_backend == "numpy" and len(df_col) == df_col.count():
                    # No NA values, can convert ints and bools
                    if col_type is np.dtype("int64") or col_type is bool:
                        self.frame[col_name] = df_col.astype(col_type)
            except KeyError:
                pass  # this column not in results

    def _sqlalchemy_type(self, col: Index | Series):
        dtype: DtypeArg = self.dtype or {}
        if is_dict_like(dtype):
            dtype = cast("Mapping[Hashable, Any]", dtype)
            if col.name in dtype:
                return dtype[col.name]

        # Infer type of column, while ignoring missing values.
        # Needed for inserting typed data containing NULLs, GH 8778.
        col_type = lib.infer_dtype(col, skipna=True)

        from sqlalchemy.dialects import mysql
        from sqlalchemy.types import (
            TIMESTAMP,
            BigInteger,
            Boolean,
            Date,
            DateTime,
            Float,
            Integer,
            SmallInteger,
            Text,
            Time,
        )

        # GH#48435 MySQL defaults to whole seconds, dropping fractional ones
        mysql_dialects = ("mysql", "mariadb")

        if col_type in ("datetime64", "datetime"):
            # GH 9086: TIMESTAMP is the suggested type if the column contains
            # timezone information
            try:
                # error: Item "Index" of "Union[Index, Series]" has no attribute "dt"
                is_tz_aware = col.dt.tz is not None  # type: ignore[union-attr]
            except AttributeError:
                # The column is actually a DatetimeIndex
                # GH 26761 or an Index with date-like data e.g. 9999-01-01
                is_tz_aware = getattr(col, "tz", None) is not None
            if is_tz_aware:
                return TIMESTAMP(timezone=True).with_variant(
                    mysql.TIMESTAMP(timezone=True, fsp=6), *mysql_dialects
                )
            return DateTime().with_variant(mysql.DATETIME(fsp=6), *mysql_dialects)
        if col_type == "timedelta64":
            warnings.warn(
                "the 'timedelta' type is not supported, and will be "
                "written as integer values (ns frequency) to the database.",
                UserWarning,
                stacklevel=find_stack_level(),
            )
            return BigInteger
        elif col_type == "floating":
            if col.dtype == "float32":
                return Float(precision=23)
            else:
                return Float(precision=53)
        elif col_type == "integer":
            # GH35076 Map pandas integer to optimal SQLAlchemy integer type
            if col.dtype.name.lower() in ("int8", "uint8", "int16"):
                return SmallInteger
            elif col.dtype.name.lower() in ("uint16", "int32"):
                return Integer
            elif col.dtype.name.lower() == "uint64":
                raise ValueError("Unsigned 64 bit integer datatype is not supported")
            else:
                return BigInteger
        elif col_type == "boolean":
            return Boolean
        elif col_type == "date":
            return Date
        elif col_type == "time":
            return Time().with_variant(mysql.TIME(fsp=6), *mysql_dialects)
        elif col_type == "complex":
            raise ValueError("Complex datatypes not supported")

        return Text

    def _get_dtype(self, sqltype):
        from sqlalchemy.types import (
            TIMESTAMP,
            Boolean,
            Date,
            DateTime,
            Float,
            Integer,
            String,
        )

        if isinstance(sqltype, Float):
            return float
        elif isinstance(sqltype, Integer):
            # TODO: Refine integer size.
            return np.dtype("int64")
        elif isinstance(sqltype, TIMESTAMP):
            # we have a timezone capable type
            if not sqltype.timezone:
                return datetime
            return DatetimeTZDtype
        elif isinstance(sqltype, DateTime):
            # Caution: np.datetime64 is also a subclass of np.number.
            return datetime
        elif isinstance(sqltype, Date):
            return date
        elif isinstance(sqltype, Boolean):
            return bool
        elif isinstance(sqltype, String):
            if using_string_dtype():
                return StringDtype(na_value=np.nan)

        return object


class PandasSQL(PandasObject, ABC):
    """
    Subclasses Should define read_query and to_sql.
    """

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *args) -> None:
        pass

    def read_table(
        self,
        table_name: str,
        index_col: str | list[str] | None = None,
        coerce_float: bool = True,
        parse_dates=None,
        columns=None,
        schema: str | None = None,
        chunksize: int | None = None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> DataFrame | Iterator[DataFrame]:
        raise NotImplementedError(
            "read_sql_table not supported for a DBAPI connection. Use "
            "read_sql_query, or pass a SQLAlchemy connectable."
        )

    @abstractmethod
    def read_query(
        self,
        sql: str,
        index_col: str | list[str] | None = None,
        coerce_float: bool = True,
        parse_dates=None,
        params=None,
        chunksize: int | None = None,
        dtype: DtypeArg | None = None,
        dtype_backend: DtypeBackend | Literal["numpy"] = "numpy",
    ) -> DataFrame | Iterator[DataFrame]:
        pass

    @abstractmethod
    def to_sql(
        self,
        frame,
        name: str,
        if_exists: Literal["fail", "replace", "append", "delete_rows"] = "fail",
        index: bool = True,
        index_label=None,
        schema=None,
        chunksize: int | None = None,
        dtype: DtypeArg | None = None,
        method: Literal["multi"] | Callable[..., Any] | None = None,
        engine: str = "auto",
        **engine_kwargs,
    ) -> int | None:
        pass

    @abstractmethod
    def execute(self, sql: str | Select[Any] | TextClause, params=None):
        pass

    @abstractmethod
    def has_table(self, name: str, schema: str | None = None) -> bool:
        pass

    @abstractmethod
    def _create_sql_schema(
        self,
        frame: DataFrame,
        table_name: str,
        keys: list[str] | None = None,
        dtype: DtypeArg | None = None,
        schema: str | None = None,
    ) -> str:
        pass
