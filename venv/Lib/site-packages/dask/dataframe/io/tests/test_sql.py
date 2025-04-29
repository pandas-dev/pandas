from __future__ import annotations

import io
import sys
from contextlib import contextmanager

import pytest

from dask.dataframe._compat import PANDAS_GE_300
from dask.dataframe.io.sql import read_sql, read_sql_query, read_sql_table
from dask.dataframe.utils import assert_eq, get_string_dtype
from dask.utils import tmpfile

pd = pytest.importorskip("pandas")
dd = pytest.importorskip("dask.dataframe")
pytest.importorskip("sqlalchemy")
pytest.importorskip("sqlite3")
np = pytest.importorskip("numpy")


data = """
name,number,age,negish
Alice,0,33,-5
Bob,1,40,-3
Chris,2,22,3
Dora,3,16,5
Edith,4,53,0
Francis,5,30,0
Garreth,6,20,0
"""

df = pd.read_csv(io.StringIO(data), index_col="number")


@pytest.fixture
def db():
    with tmpfile() as f:
        uri = "sqlite:///%s" % f
        df.to_sql("test", uri, index=True, if_exists="replace")
        yield uri


def test_empty(db):
    from sqlalchemy import Column, Integer, MetaData, Table, create_engine

    with tmpfile() as f:
        uri = "sqlite:///%s" % f
        metadata = MetaData()
        engine = create_engine(uri)
        table = Table(
            "empty_table",
            metadata,
            Column("id", Integer, primary_key=True),
            Column("col2", Integer),
        )
        metadata.create_all(engine)

        dask_df = read_sql_table(table.name, uri, index_col="id", npartitions=1)
        assert dask_df.index.name == "id"
        # The dtype of the empty result might no longer be as expected
        # assert dask_df.col2.dtype == np.dtype("int64")
        pd_dataframe = dask_df.compute()
        assert pd_dataframe.empty is True


@pytest.mark.filterwarnings(
    "ignore:The default dtype for empty Series will be 'object' instead of 'float64'"
)
@pytest.mark.parametrize("use_head", [True, False])
def test_single_column(db, use_head):
    from sqlalchemy import Column, Integer, MetaData, Table, create_engine

    with tmpfile() as f:
        uri = "sqlite:///%s" % f
        metadata = MetaData()
        engine = create_engine(uri)
        table = Table(
            "single_column",
            metadata,
            Column("id", Integer, primary_key=True),
        )
        metadata.create_all(engine)
        test_data = pd.DataFrame({"id": list(range(50))}).set_index("id")
        test_data.to_sql(table.name, uri, index=True, if_exists="replace")

        if use_head:
            dask_df = read_sql_table(table.name, uri, index_col="id", npartitions=2)
        else:
            dask_df = read_sql_table(
                table.name,
                uri,
                head_rows=0,
                npartitions=2,
                meta=test_data.iloc[:0],
                index_col="id",
            )
        assert dask_df.index.name == "id"
        assert dask_df.npartitions == 2
        pd_dataframe = dask_df.compute()
        assert_eq(test_data, pd_dataframe)


def test_passing_engine_as_uri_raises_helpful_error(db):
    # https://github.com/dask/dask/issues/6473
    from sqlalchemy import create_engine

    df = pd.DataFrame([{"i": i, "s": str(i) * 2} for i in range(4)])
    ddf = dd.from_pandas(df, npartitions=2)

    with tmpfile() as f:
        db = "sqlite:///%s" % f
        engine = create_engine(db)
        with pytest.raises(ValueError, match="Expected URI to be a string"):
            ddf.to_sql("test", engine, if_exists="replace")


@pytest.mark.skip(
    reason="Requires a postgres server. Sqlite does not support multiple schemas."
)
def test_empty_other_schema():
    from sqlalchemy import DDL, Column, Integer, MetaData, Table, create_engine, event

    # Database configurations.
    pg_host = "localhost"
    pg_port = "5432"
    pg_user = "user"
    pg_pass = "pass"
    pg_db = "db"
    db_url = f"postgresql://{pg_user}:{pg_pass}@{pg_host}:{pg_port}/{pg_db}"

    # Create an empty table in a different schema.
    table_name = "empty_table"
    schema_name = "other_schema"
    engine = create_engine(db_url)
    metadata = MetaData()
    table = Table(
        table_name,
        metadata,
        Column("id", Integer, primary_key=True),
        Column("col2", Integer),
        schema=schema_name,
    )
    # Create the schema and the table.
    event.listen(
        metadata, "before_create", DDL("CREATE SCHEMA IF NOT EXISTS %s" % schema_name)
    )
    metadata.create_all(engine)

    # Read the empty table from the other schema.
    dask_df = read_sql_table(
        table.name, db_url, index_col="id", schema=table.schema, npartitions=1
    )

    # Validate that the retrieved table is empty.
    assert dask_df.index.name == "id"
    assert dask_df.col2.dtype == np.dtype("int64")
    pd_dataframe = dask_df.compute()
    assert pd_dataframe.empty is True

    # Drop the schema and the table.
    engine.execute("DROP SCHEMA IF EXISTS %s CASCADE" % schema_name)


def test_needs_rational(db):
    import datetime

    now = datetime.datetime.now()
    d = datetime.timedelta(seconds=1)
    df = pd.DataFrame(
        {
            "a": list("ghjkl"),
            "b": [now + i * d for i in range(5)],
            "c": [True, True, False, True, True],
        }
    )
    df = pd.concat(
        [
            df,
            pd.DataFrame(
                [
                    {"a": "x", "b": now + d * 1000, "c": None},
                    {"a": None, "b": now + d * 1001, "c": None},
                ]
            ),
        ]
    )
    string_dtype = get_string_dtype()
    with tmpfile() as f:
        uri = "sqlite:///%s" % f
        df.to_sql("test", uri, index=False, if_exists="replace")

        # one partition contains NULL
        data = read_sql_table("test", uri, npartitions=2, index_col="b")
        df2 = df.set_index("b")
        assert_eq(data, df2.astype({"c": bool}))  # bools are coerced

        # one partition contains NULL, but big enough head
        data = read_sql_table("test", uri, npartitions=2, index_col="b", head_rows=12)
        df2 = df.set_index("b")
        assert_eq(data, df2)

        # empty partitions
        data = read_sql_table("test", uri, npartitions=20, index_col="b")
        part = data.get_partition(12).compute()
        assert part.dtypes.tolist() == [string_dtype, bool]
        assert part.empty
        df2 = df.set_index("b")
        assert_eq(data, df2.astype({"c": bool}))

        # explicit meta
        data = read_sql_table("test", uri, npartitions=2, index_col="b", meta=df2[:0])
        part = data.get_partition(1).compute()
        assert part.dtypes.tolist() == [string_dtype, string_dtype]
        df2 = df.set_index("b")
        assert_eq(data, df2, check_dtype=False)


def test_simple(db):
    # single chunk
    data = read_sql_table("test", db, npartitions=2, index_col="number").compute()
    assert (data.name == df.name).all()
    assert data.index.name == "number"
    assert_eq(data, df)


def test_npartitions(db):
    data = read_sql_table(
        "test", db, columns=list(df.columns), npartitions=2, index_col="number"
    )
    assert len(data.divisions) == 3
    assert (data.name.compute() == df.name).all()
    data = read_sql_table(
        "test", db, columns=["name"], npartitions=6, index_col="number"
    )
    assert_eq(data, df[["name"]])
    data = read_sql_table(
        "test",
        db,
        columns=list(df.columns),
        bytes_per_chunk="2 GiB",
        index_col="number",
    )
    assert data.npartitions == 1
    assert (data.name.compute() == df.name).all()

    data_1 = read_sql_table(
        "test",
        db,
        columns=list(df.columns),
        bytes_per_chunk=2**30,
        index_col="number",
        head_rows=1,
    )
    assert data_1.npartitions == 1
    assert (data_1.name.compute() == df.name).all()

    data = read_sql_table(
        "test",
        db,
        columns=list(df.columns),
        bytes_per_chunk=250,
        index_col="number",
        head_rows=1,
    )
    assert (
        (data.memory_usage_per_partition(deep=True, index=True) < 400).compute().all()
    )
    assert (data.name.compute() == df.name).all()


def test_divisions(db):
    data = read_sql_table(
        "test", db, columns=["name"], divisions=[0, 2, 4], index_col="number"
    )
    assert data.divisions == (0, 2, 4)
    assert data.index.max().compute() == 4
    assert_eq(data, df[["name"]][df.index <= 4])


@pytest.mark.xfail(PANDAS_GE_300, reason="memory doesn't match")
def test_division_or_partition(db):
    with pytest.raises(TypeError, match="either 'divisions' or 'npartitions'"):
        read_sql_table(
            "test",
            db,
            columns=["name"],
            index_col="number",
            divisions=[0, 2, 4],
            npartitions=3,
        )

    out = read_sql_table("test", db, index_col="number", bytes_per_chunk=100)
    m = out.memory_usage_per_partition(deep=True, index=True).compute()
    assert (50 < m).all() and (m < 200).all()
    assert_eq(out, df)


def test_meta(db):
    data = read_sql_table(
        "test", db, index_col="number", meta=dd.from_pandas(df, npartitions=1)
    ).compute()
    assert (data.name == df.name).all()
    assert data.index.name == "number"
    assert_eq(data, df)


def test_meta_no_head_rows(db):
    data = read_sql_table(
        "test",
        db,
        index_col="number",
        meta=dd.from_pandas(df, npartitions=1),
        npartitions=2,
        head_rows=0,
    )
    assert len(data.divisions) == 3
    data = data.compute()
    assert (data.name == df.name).all()
    assert data.index.name == "number"
    assert_eq(data, df)

    data = read_sql_table(
        "test",
        db,
        index_col="number",
        meta=dd.from_pandas(df, npartitions=1),
        divisions=[0, 3, 6],
        head_rows=0,
    )
    assert len(data.divisions) == 3
    data = data.compute()
    assert (data.name == df.name).all()
    assert data.index.name == "number"
    assert_eq(data, df)


def test_no_meta_no_head_rows(db):
    with pytest.raises(ValueError):
        read_sql_table("test", db, index_col="number", head_rows=0, npartitions=1)


def test_limits(db):
    data = read_sql_table("test", db, npartitions=2, index_col="number", limits=[1, 4])
    assert data.index.min().compute() == 1
    assert data.index.max().compute() == 4


def test_datetimes():
    import datetime

    now = datetime.datetime.now()
    d = datetime.timedelta(seconds=1)
    df = pd.DataFrame(
        {"a": list("ghjkl"), "b": [now + i * d for i in range(2, -3, -1)]}
    )
    with tmpfile() as f:
        uri = "sqlite:///%s" % f
        df.to_sql("test", uri, index=False, if_exists="replace")
        data = read_sql_table("test", uri, npartitions=2, index_col="b")
        assert data.index.dtype.kind == "M"
        assert data.divisions[0] == df.b.min()
        df2 = df.set_index("b")
        assert_eq(
            data.map_partitions(lambda x: x.sort_index()),
            df2.sort_index(),
            check_dtype=False,
        )


def test_extra_connection_engine_keywords(caplog, db):
    data = read_sql_table(
        "test", db, npartitions=2, index_col="number", engine_kwargs={"echo": False}
    ).compute()
    # no captured message from the stdout with the echo=False parameter (this is the default)
    out = "\n".join(r.message for r in caplog.records)
    assert out == ""
    assert_eq(data, df)
    # with the echo=True sqlalchemy parameter, you should get all SQL queries in the stdout
    data = read_sql_table(
        "test", db, npartitions=2, index_col="number", engine_kwargs={"echo": True}
    ).compute()
    out = "\n".join(r.message for r in caplog.records)
    assert "WHERE" in out
    assert "FROM" in out
    assert "SELECT" in out
    assert "AND" in out
    assert ">= ?" in out
    assert "< ?" in out
    assert "<= ?" in out
    assert_eq(data, df)


def test_query(db):
    import sqlalchemy as sa
    from sqlalchemy import sql

    s1 = sql.select(sql.column("number"), sql.column("name")).select_from(
        sql.table("test")
    )
    out = read_sql_query(s1, db, npartitions=2, index_col="number")
    assert_eq(out, df[["name"]])

    s2 = (
        sql.select(
            sa.cast(sql.column("number"), sa.types.BigInteger).label("number"),
            sql.column("name"),
        )
        .where(sql.column("number") >= 5)
        .select_from(sql.table("test"))
    )

    out = read_sql_query(s2, db, npartitions=2, index_col="number")
    assert_eq(out, df.loc[5:, ["name"]])


def test_query_index_from_query(db):
    from sqlalchemy import sql

    number = sql.column("number")
    name = sql.column("name")
    s1 = sql.select(number, name, sql.func.length(name).label("lenname")).select_from(
        sql.table("test")
    )
    out = read_sql_query(s1, db, npartitions=2, index_col="lenname")

    lenname_df = df.copy()
    lenname_df["lenname"] = lenname_df["name"].str.len()
    lenname_df = lenname_df.reset_index().set_index("lenname")
    assert_eq(out, lenname_df.loc[:, ["number", "name"]])


def test_query_with_meta(db):
    from sqlalchemy import sql

    data = {
        "name": pd.Series([], name="name", dtype="str"),
        "age": pd.Series([], name="age", dtype="int"),
    }
    index = pd.Index([], name="number", dtype="int")
    meta = pd.DataFrame(data, index=index)

    s1 = sql.select(
        sql.column("number"), sql.column("name"), sql.column("age")
    ).select_from(sql.table("test"))
    out = read_sql_query(s1, db, npartitions=2, index_col="number", meta=meta)
    # Don't check dtype for windows https://github.com/dask/dask/issues/8620
    assert_eq(out, df[["name", "age"]], check_dtype=sys.platform != "win32")


def test_no_character_index_without_divisions(db):
    # attempt to read the sql table with a character index and no divisions
    with pytest.raises(TypeError):
        read_sql_table("test", db, npartitions=2, index_col="name", divisions=None)


def test_read_sql(db):
    from sqlalchemy import sql

    s = sql.select(sql.column("number"), sql.column("name")).select_from(
        sql.table("test")
    )
    out = read_sql(s, db, npartitions=2, index_col="number")
    assert_eq(out, df[["name"]])

    data = read_sql_table("test", db, npartitions=2, index_col="number").compute()
    assert (data.name == df.name).all()
    assert data.index.name == "number"
    assert_eq(data, df)


@contextmanager
def tmp_db_uri():
    with tmpfile() as f:
        yield "sqlite:///%s" % f


@pytest.mark.parametrize("npartitions", (1, 2))
@pytest.mark.parametrize("parallel", (False, True))
def test_to_sql(npartitions, parallel):
    df_by_age = df.set_index("age")
    df_appended = pd.concat(
        [
            df,
            df,
        ]
    )

    ddf = dd.from_pandas(df, npartitions)
    ddf_by_age = ddf.set_index("age")

    # Simple round trip test: use existing "number" index_col
    with tmp_db_uri() as uri:
        ddf.to_sql("test", uri, parallel=parallel)
        result = read_sql_table("test", uri, "number")
        assert_eq(df, result)

    # Test writing no index, and reading back in with one of the other columns as index (`read_sql_table` requires
    # an index_col)
    with tmp_db_uri() as uri:
        ddf.to_sql("test", uri, parallel=parallel, index=False)

        result = read_sql_table("test", uri, "negish")
        assert_eq(df.set_index("negish"), result)

        result = read_sql_table("test", uri, "age")
        assert_eq(df_by_age, result)

    # Index by "age" instead
    with tmp_db_uri() as uri:
        ddf_by_age.to_sql("test", uri, parallel=parallel)
        result = read_sql_table("test", uri, "age")
        assert_eq(df_by_age, result)

    # Index column can't have "object" dtype if no partitions are provided
    with tmp_db_uri() as uri:
        ddf.set_index("name").to_sql("test", uri)
        with pytest.raises(
            TypeError,
            match='Provided index column is of type "object".  If divisions is not provided the index column type must be numeric or datetime.',  # noqa: E501
        ):
            read_sql_table("test", uri, "name")

    # Test various "if_exists" values
    with tmp_db_uri() as uri:
        ddf.to_sql("test", uri)

        # Writing a table that already exists fails
        with pytest.raises(ValueError, match="Table 'test' already exists"):
            ddf.to_sql("test", uri)

        ddf.to_sql("test", uri, parallel=parallel, if_exists="append")
        result = read_sql_table("test", uri, "number")

        assert_eq(df_appended, result)

        ddf_by_age.to_sql("test", uri, parallel=parallel, if_exists="replace")
        result = read_sql_table("test", uri, "age")
        assert_eq(df_by_age, result)

    # Verify number of partitions returned, when compute=False
    with tmp_db_uri() as uri:
        result = ddf.to_sql("test", uri, parallel=parallel, compute=False)

        # the first result is from the "meta" insert
        actual = len(result.compute())

        assert actual == npartitions


def test_to_sql_kwargs():
    ddf = dd.from_pandas(df, 2)
    with tmp_db_uri() as uri:
        ddf.to_sql("test", uri, method="multi")
        with pytest.raises(
            TypeError, match="to_sql\\(\\) got an unexpected keyword argument 'unknown'"
        ):
            ddf.to_sql("test", uri, unknown=None)


def test_to_sql_engine_kwargs(caplog):
    ddf = dd.from_pandas(df, 2)
    with tmp_db_uri() as uri:
        ddf.to_sql("test", uri, engine_kwargs={"echo": False})
        logs = "\n".join(r.message for r in caplog.records)
        assert logs == ""
        assert_eq(df, read_sql_table("test", uri, "number"))

    with tmp_db_uri() as uri:
        ddf.to_sql("test", uri, engine_kwargs={"echo": True})
        logs = "\n".join(r.message for r in caplog.records)
        assert "CREATE" in logs
        assert "INSERT" in logs

        assert_eq(df, read_sql_table("test", uri, "number"))
