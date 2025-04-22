from __future__ import annotations

import contextlib
from contextlib import closing
import csv
from datetime import (
    date,
    datetime,
    time,
    timedelta,
)
from io import StringIO
from pathlib import Path
import sqlite3
from typing import TYPE_CHECKING
import uuid

import numpy as np
import pytest

from pandas._config import using_string_dtype

from pandas._libs import lib
from pandas.compat import pa_version_under14p1
from pandas.compat._optional import import_optional_dependency
import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    Timestamp,
    concat,
    date_range,
    isna,
    to_datetime,
    to_timedelta,
)
import pandas._testing as tm
from pandas.util.version import Version

from pandas.io import sql
from pandas.io.sql import (
    SQLAlchemyEngine,
    SQLDatabase,
    SQLiteDatabase,
    get_engine,
    pandasSQL_builder,
    read_sql_query,
    read_sql_table,
)

if TYPE_CHECKING:
    import sqlalchemy


pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
    ),
    pytest.mark.single_cpu,
]


@pytest.fixture
def sql_strings():
    return {
        "read_parameters": {
            "sqlite": "SELECT * FROM iris WHERE Name=? AND SepalLength=?",
            "mysql": "SELECT * FROM iris WHERE `Name`=%s AND `SepalLength`=%s",
            "postgresql": 'SELECT * FROM iris WHERE "Name"=%s AND "SepalLength"=%s',
        },
        "read_named_parameters": {
            "sqlite": """
                SELECT * FROM iris WHERE Name=:name AND SepalLength=:length
                """,
            "mysql": """
                SELECT * FROM iris WHERE
                `Name`=%(name)s AND `SepalLength`=%(length)s
                """,
            "postgresql": """
                SELECT * FROM iris WHERE
                "Name"=%(name)s AND "SepalLength"=%(length)s
                """,
        },
        "read_no_parameters_with_percent": {
            "sqlite": "SELECT * FROM iris WHERE Name LIKE '%'",
            "mysql": "SELECT * FROM iris WHERE `Name` LIKE '%'",
            "postgresql": "SELECT * FROM iris WHERE \"Name\" LIKE '%'",
        },
    }


def iris_table_metadata():
    import sqlalchemy
    from sqlalchemy import (
        Column,
        Double,
        Float,
        MetaData,
        String,
        Table,
    )

    dtype = Double if Version(sqlalchemy.__version__) >= Version("2.0.0") else Float
    metadata = MetaData()
    iris = Table(
        "iris",
        metadata,
        Column("SepalLength", dtype),
        Column("SepalWidth", dtype),
        Column("PetalLength", dtype),
        Column("PetalWidth", dtype),
        Column("Name", String(200)),
    )
    return iris


def create_and_load_iris_sqlite3(conn, iris_file: Path):
    stmt = """CREATE TABLE iris (
            "SepalLength" REAL,
            "SepalWidth" REAL,
            "PetalLength" REAL,
            "PetalWidth" REAL,
            "Name" TEXT
        )"""

    cur = conn.cursor()
    cur.execute(stmt)
    with iris_file.open(newline=None, encoding="utf-8") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        stmt = "INSERT INTO iris VALUES(?, ?, ?, ?, ?)"
        # ADBC requires explicit types - no implicit str -> float conversion
        records = []
        records = [
            (
                float(row[0]),
                float(row[1]),
                float(row[2]),
                float(row[3]),
                row[4],
            )
            for row in reader
        ]

        cur.executemany(stmt, records)
    cur.close()

    conn.commit()


def create_and_load_iris_postgresql(conn, iris_file: Path):
    stmt = """CREATE TABLE iris (
            "SepalLength" DOUBLE PRECISION,
            "SepalWidth" DOUBLE PRECISION,
            "PetalLength" DOUBLE PRECISION,
            "PetalWidth" DOUBLE PRECISION,
            "Name" TEXT
        )"""
    with conn.cursor() as cur:
        cur.execute(stmt)
        with iris_file.open(newline=None, encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            next(reader)
            stmt = "INSERT INTO iris VALUES($1, $2, $3, $4, $5)"
            # ADBC requires explicit types - no implicit str -> float conversion
            records = [
                (
                    float(row[0]),
                    float(row[1]),
                    float(row[2]),
                    float(row[3]),
                    row[4],
                )
                for row in reader
            ]

            cur.executemany(stmt, records)

    conn.commit()


def create_and_load_iris(conn, iris_file: Path):
    from sqlalchemy import insert

    iris = iris_table_metadata()

    with iris_file.open(newline=None, encoding="utf-8") as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)
        params = [dict(zip(header, row)) for row in reader]
        stmt = insert(iris).values(params)
        with conn.begin() as con:
            iris.drop(con, checkfirst=True)
            iris.create(bind=con)
            con.execute(stmt)


def create_and_load_iris_view(conn):
    stmt = "CREATE VIEW iris_view AS SELECT * FROM iris"
    if isinstance(conn, sqlite3.Connection):
        cur = conn.cursor()
        cur.execute(stmt)
    else:
        adbc = import_optional_dependency("adbc_driver_manager.dbapi", errors="ignore")
        if adbc and isinstance(conn, adbc.Connection):
            with conn.cursor() as cur:
                cur.execute(stmt)
            conn.commit()
        else:
            from sqlalchemy import text

            stmt = text(stmt)
            with conn.begin() as con:
                con.execute(stmt)


def types_table_metadata(dialect: str):
    from sqlalchemy import (
        TEXT,
        Boolean,
        Column,
        DateTime,
        Float,
        Integer,
        MetaData,
        Table,
    )

    date_type = TEXT if dialect == "sqlite" else DateTime
    bool_type = Integer if dialect == "sqlite" else Boolean
    metadata = MetaData()
    types = Table(
        "types",
        metadata,
        Column("TextCol", TEXT),
        # error: Cannot infer type argument 1 of "Column"
        Column("DateCol", date_type),  # type: ignore[misc]
        Column("IntDateCol", Integer),
        Column("IntDateOnlyCol", Integer),
        Column("FloatCol", Float),
        Column("IntCol", Integer),
        # error: Cannot infer type argument 1 of "Column"
        Column("BoolCol", bool_type),  # type: ignore[misc]
        Column("IntColWithNull", Integer),
        # error: Cannot infer type argument 1 of "Column"
        Column("BoolColWithNull", bool_type),  # type: ignore[misc]
    )
    return types


def create_and_load_types_sqlite3(conn, types_data: list[dict]):
    stmt = """CREATE TABLE types (
                    "TextCol" TEXT,
                    "DateCol" TEXT,
                    "IntDateCol" INTEGER,
                    "IntDateOnlyCol" INTEGER,
                    "FloatCol" REAL,
                    "IntCol" INTEGER,
                    "BoolCol" INTEGER,
                    "IntColWithNull" INTEGER,
                    "BoolColWithNull" INTEGER
                )"""

    ins_stmt = """
                INSERT INTO types
                VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)
                """

    if isinstance(conn, sqlite3.Connection):
        cur = conn.cursor()
        cur.execute(stmt)
        cur.executemany(ins_stmt, types_data)
    else:
        with conn.cursor() as cur:
            cur.execute(stmt)
            cur.executemany(ins_stmt, types_data)

        conn.commit()


def create_and_load_types_postgresql(conn, types_data: list[dict]):
    with conn.cursor() as cur:
        stmt = """CREATE TABLE types (
                        "TextCol" TEXT,
                        "DateCol" TIMESTAMP,
                        "IntDateCol" INTEGER,
                        "IntDateOnlyCol" INTEGER,
                        "FloatCol" DOUBLE PRECISION,
                        "IntCol" INTEGER,
                        "BoolCol" BOOLEAN,
                        "IntColWithNull" INTEGER,
                        "BoolColWithNull" BOOLEAN
                    )"""
        cur.execute(stmt)

        stmt = """
                INSERT INTO types
                VALUES($1, $2::timestamp, $3, $4, $5, $6, $7, $8, $9)
                """

        cur.executemany(stmt, types_data)

    conn.commit()


def create_and_load_types(conn, types_data: list[dict], dialect: str):
    from sqlalchemy import insert
    from sqlalchemy.engine import Engine

    types = types_table_metadata(dialect)

    stmt = insert(types).values(types_data)
    if isinstance(conn, Engine):
        with conn.connect() as conn:
            with conn.begin():
                types.drop(conn, checkfirst=True)
                types.create(bind=conn)
                conn.execute(stmt)
    else:
        with conn.begin():
            types.drop(conn, checkfirst=True)
            types.create(bind=conn)
            conn.execute(stmt)


def create_and_load_postgres_datetz(conn):
    from sqlalchemy import (
        Column,
        DateTime,
        MetaData,
        Table,
        insert,
    )
    from sqlalchemy.engine import Engine

    metadata = MetaData()
    datetz = Table("datetz", metadata, Column("DateColWithTz", DateTime(timezone=True)))
    datetz_data = [
        {
            "DateColWithTz": "2000-01-01 00:00:00-08:00",
        },
        {
            "DateColWithTz": "2000-06-01 00:00:00-07:00",
        },
    ]
    stmt = insert(datetz).values(datetz_data)
    if isinstance(conn, Engine):
        with conn.connect() as conn:
            with conn.begin():
                datetz.drop(conn, checkfirst=True)
                datetz.create(bind=conn)
                conn.execute(stmt)
    else:
        with conn.begin():
            datetz.drop(conn, checkfirst=True)
            datetz.create(bind=conn)
            conn.execute(stmt)

    # "2000-01-01 00:00:00-08:00" should convert to
    # "2000-01-01 08:00:00"
    # "2000-06-01 00:00:00-07:00" should convert to
    # "2000-06-01 07:00:00"
    # GH 6415
    expected_data = [
        Timestamp("2000-01-01 08:00:00", tz="UTC"),
        Timestamp("2000-06-01 07:00:00", tz="UTC"),
    ]
    return Series(expected_data, name="DateColWithTz").astype("M8[us, UTC]")


def check_iris_frame(frame: DataFrame):
    pytype = frame.dtypes.iloc[0].type
    row = frame.iloc[0]
    assert issubclass(pytype, np.floating)
    tm.assert_series_equal(
        row, Series([5.1, 3.5, 1.4, 0.2, "Iris-setosa"], index=frame.columns, name=0)
    )
    assert frame.shape in ((150, 5), (8, 5))


def count_rows(conn, table_name: str):
    stmt = f"SELECT count(*) AS count_1 FROM {table_name}"
    adbc = import_optional_dependency("adbc_driver_manager.dbapi", errors="ignore")
    if isinstance(conn, sqlite3.Connection):
        cur = conn.cursor()
        return cur.execute(stmt).fetchone()[0]
    elif adbc and isinstance(conn, adbc.Connection):
        with conn.cursor() as cur:
            cur.execute(stmt)
            return cur.fetchone()[0]
    else:
        from sqlalchemy import create_engine
        from sqlalchemy.engine import Engine

        if isinstance(conn, str):
            try:
                engine = create_engine(conn)
                with engine.connect() as conn:
                    return conn.exec_driver_sql(stmt).scalar_one()
            finally:
                engine.dispose()
        elif isinstance(conn, Engine):
            with conn.connect() as conn:
                return conn.exec_driver_sql(stmt).scalar_one()
        else:
            return conn.exec_driver_sql(stmt).scalar_one()


@pytest.fixture
def iris_path(datapath):
    iris_path = datapath("io", "data", "csv", "iris.csv")
    return Path(iris_path)


@pytest.fixture
def types_data():
    return [
        {
            "TextCol": "first",
            "DateCol": "2000-01-03 00:00:00",
            "IntDateCol": 535852800,
            "IntDateOnlyCol": 20101010,
            "FloatCol": 10.10,
            "IntCol": 1,
            "BoolCol": False,
            "IntColWithNull": 1,
            "BoolColWithNull": False,
        },
        {
            "TextCol": "first",
            "DateCol": "2000-01-04 00:00:00",
            "IntDateCol": 1356998400,
            "IntDateOnlyCol": 20101212,
            "FloatCol": 10.10,
            "IntCol": 1,
            "BoolCol": False,
            "IntColWithNull": None,
            "BoolColWithNull": None,
        },
    ]


@pytest.fixture
def types_data_frame(types_data):
    dtypes = {
        "TextCol": "str",
        "DateCol": "str",
        "IntDateCol": "int64",
        "IntDateOnlyCol": "int64",
        "FloatCol": "float",
        "IntCol": "int64",
        "BoolCol": "int64",
        "IntColWithNull": "float",
        "BoolColWithNull": "float",
    }
    df = DataFrame(types_data)
    return df[dtypes.keys()].astype(dtypes)


@pytest.fixture
def test_frame1():
    columns = ["index", "A", "B", "C", "D"]
    data = [
        (
            "2000-01-03 00:00:00",
            0.980268513777,
            3.68573087906,
            -0.364216805298,
            -1.15973806169,
        ),
        (
            "2000-01-04 00:00:00",
            1.04791624281,
            -0.0412318367011,
            -0.16181208307,
            0.212549316967,
        ),
        (
            "2000-01-05 00:00:00",
            0.498580885705,
            0.731167677815,
            -0.537677223318,
            1.34627041952,
        ),
        (
            "2000-01-06 00:00:00",
            1.12020151869,
            1.56762092543,
            0.00364077397681,
            0.67525259227,
        ),
    ]
    return DataFrame(data, columns=columns)


@pytest.fixture
def test_frame3():
    columns = ["index", "A", "B"]
    data = [
        ("2000-01-03 00:00:00", 2**31 - 1, -1.987670),
        ("2000-01-04 00:00:00", -29, -0.0412318367011),
        ("2000-01-05 00:00:00", 20000, 0.731167677815),
        ("2000-01-06 00:00:00", -290867, 1.56762092543),
    ]
    return DataFrame(data, columns=columns)


def get_all_views(conn):
    if isinstance(conn, sqlite3.Connection):
        c = conn.execute("SELECT name FROM sqlite_master WHERE type='view'")
        return [view[0] for view in c.fetchall()]
    else:
        adbc = import_optional_dependency("adbc_driver_manager.dbapi", errors="ignore")
        if adbc and isinstance(conn, adbc.Connection):
            results = []
            info = conn.adbc_get_objects().read_all().to_pylist()
            for catalog in info:
                catalog["catalog_name"]
                for schema in catalog["catalog_db_schemas"]:
                    schema["db_schema_name"]
                    for table in schema["db_schema_tables"]:
                        if table["table_type"] == "view":
                            view_name = table["table_name"]
                            results.append(view_name)

            return results
        else:
            from sqlalchemy import inspect

            return inspect(conn).get_view_names()


def get_all_tables(conn):
    if isinstance(conn, sqlite3.Connection):
        c = conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
        return [table[0] for table in c.fetchall()]
    else:
        adbc = import_optional_dependency("adbc_driver_manager.dbapi", errors="ignore")

        if adbc and isinstance(conn, adbc.Connection):
            results = []
            info = conn.adbc_get_objects().read_all().to_pylist()
            for catalog in info:
                for schema in catalog["catalog_db_schemas"]:
                    for table in schema["db_schema_tables"]:
                        if table["table_type"] == "table":
                            table_name = table["table_name"]
                            results.append(table_name)

            return results
        else:
            from sqlalchemy import inspect

            return inspect(conn).get_table_names()


def drop_table(
    table_name: str,
    conn: sqlite3.Connection | sqlalchemy.engine.Engine | sqlalchemy.engine.Connection,
):
    if isinstance(conn, sqlite3.Connection):
        conn.execute(f"DROP TABLE IF EXISTS {sql._get_valid_sqlite_name(table_name)}")
        conn.commit()

    else:
        adbc = import_optional_dependency("adbc_driver_manager.dbapi", errors="ignore")
        if adbc and isinstance(conn, adbc.Connection):
            with conn.cursor() as cur:
                cur.execute(f'DROP TABLE IF EXISTS "{table_name}"')
        else:
            with conn.begin() as con:
                with sql.SQLDatabase(con) as db:
                    db.drop_table(table_name)


def drop_view(
    view_name: str,
    conn: sqlite3.Connection | sqlalchemy.engine.Engine | sqlalchemy.engine.Connection,
):
    import sqlalchemy

    if isinstance(conn, sqlite3.Connection):
        conn.execute(f"DROP VIEW IF EXISTS {sql._get_valid_sqlite_name(view_name)}")
        conn.commit()
    else:
        adbc = import_optional_dependency("adbc_driver_manager.dbapi", errors="ignore")
        if adbc and isinstance(conn, adbc.Connection):
            with conn.cursor() as cur:
                cur.execute(f'DROP VIEW IF EXISTS "{view_name}"')
        else:
            quoted_view = conn.engine.dialect.identifier_preparer.quote_identifier(
                view_name
            )
            stmt = sqlalchemy.text(f"DROP VIEW IF EXISTS {quoted_view}")
            with conn.begin() as con:
                con.execute(stmt)  # type: ignore[union-attr]


@pytest.fixture
def mysql_pymysql_engine():
    pytest.skip("Skipping MySQL tests")
    # sqlalchemy = pytest.importorskip("sqlalchemy")
    # pymysql = pytest.importorskip("pymysql")
    # engine = sqlalchemy.create_engine(
    #     "mysql+pymysql://root@localhost:3306/pandas",
    #     connect_args={"client_flag": pymysql.constants.CLIENT.MULTI_STATEMENTS},
    #     poolclass=sqlalchemy.pool.NullPool,
    # )
    # yield engine
    # for view in get_all_views(engine):
    #     drop_view(view, engine)
    # for tbl in get_all_tables(engine):
    #     drop_table(tbl, engine)
    # engine.dispose()


@pytest.fixture
def mysql_pymysql_engine_iris(mysql_pymysql_engine, iris_path):
    pytest.skip("Skipping MySQL tests")
    # create_and_load_iris(mysql_pymysql_engine, iris_path)
    # create_and_load_iris_view(mysql_pymysql_engine)
    # return mysql_pymysql_engine


@pytest.fixture
def mysql_pymysql_engine_types(mysql_pymysql_engine, types_data):
    pytest.skip("Skipping MySQL tests")
    # create_and_load_types(mysql_pymysql_engine, types_data, "mysql")
    # return mysql_pymysql_engine


@pytest.fixture
def mysql_pymysql_conn(mysql_pymysql_engine):
    pytest.skip("Skipping MySQL tests")
    # with mysql_pymysql_engine.connect() as conn:
    #     yield conn


@pytest.fixture
def mysql_pymysql_conn_iris(mysql_pymysql_engine_iris):
    pytest.skip("Skipping MySQL tests")
    # with mysql_pymysql_engine_iris.connect() as conn:
    #     yield conn


@pytest.fixture
def mysql_pymysql_conn_types(mysql_pymysql_engine_types):
    pytest.skip("Skipping MySQL tests")
    # with mysql_pymysql_engine_types.connect() as conn:
    #     yield conn


@pytest.fixture
def postgresql_psycopg2_engine():
    pytest.skip("Skipping MySQL tests")
    # sqlalchemy = pytest.importorskip("sqlalchemy")
    # pytest.importorskip("psycopg2")
    # engine = sqlalchemy.create_engine(
    #     "postgresql+psycopg2://postgres:postgres@localhost:5432/pandas",
    #     poolclass=sqlalchemy.pool.NullPool,
    # )
    # yield engine
    # for view in get_all_views(engine):
    #     drop_view(view, engine)
    # for tbl in get_all_tables(engine):
    #     drop_table(tbl, engine)
    # engine.dispose()


@pytest.fixture
def postgresql_psycopg2_engine_iris(postgresql_psycopg2_engine, iris_path):
    pytest.skip("Skipping MySQL tests")
    # create_and_load_iris(postgresql_psycopg2_engine, iris_path)
    # create_and_load_iris_view(postgresql_psycopg2_engine)
    # return postgresql_psycopg2_engine


@pytest.fixture
def postgresql_psycopg2_engine_types(postgresql_psycopg2_engine, types_data):
    pytest.skip("Skipping MySQL tests")
    # create_and_load_types(postgresql_psycopg2_engine, types_data, "postgres")
    # return postgresql_psycopg2_engine


@pytest.fixture
def postgresql_psycopg2_conn(postgresql_psycopg2_engine):
    pytest.skip("Skipping MySQL tests")
    # with postgresql_psycopg2_engine.connect() as conn:
    #     yield conn


@pytest.fixture
def postgresql_adbc_conn():
    pytest.skip("Skipping MySQL tests")
    # pytest.importorskip("pyarrow")
    # pytest.importorskip("adbc_driver_postgresql")
    # from adbc_driver_postgresql import dbapi

    # uri = "postgresql://postgres:postgres@localhost:5432/pandas"
    # with dbapi.connect(uri) as conn:
    #     yield conn
    #     for view in get_all_views(conn):
    #         drop_view(view, conn)
    #     for tbl in get_all_tables(conn):
    #         drop_table(tbl, conn)
    #     conn.commit()


@pytest.fixture
def postgresql_adbc_iris(postgresql_adbc_conn, iris_path):
    pytest.skip("Skipping MySQL tests")
    # import adbc_driver_manager as mgr

    # conn = postgresql_adbc_conn

    # try:
    #     conn.adbc_get_table_schema("iris")
    # except mgr.ProgrammingError:
    #     conn.rollback()
    #     create_and_load_iris_postgresql(conn, iris_path)
    # try:
    #     conn.adbc_get_table_schema("iris_view")
    # except mgr.ProgrammingError:  # note arrow-adbc issue 1022
    #     conn.rollback()
    #     create_and_load_iris_view(conn)
    # return conn


@pytest.fixture
def postgresql_adbc_types(postgresql_adbc_conn, types_data):
    pytest.skip("Skipping MySQL tests")
    # import adbc_driver_manager as mgr

    # conn = postgresql_adbc_conn

    # try:
    #     conn.adbc_get_table_schema("types")
    # except mgr.ProgrammingError:
    #     conn.rollback()
    #     new_data = [tuple(entry.values()) for entry in types_data]

    #     create_and_load_types_postgresql(conn, new_data)

    # return conn


@pytest.fixture
def postgresql_psycopg2_conn_iris(postgresql_psycopg2_engine_iris):
    pytest.skip("Skipping MySQL tests")
    # with postgresql_psycopg2_engine_iris.connect() as conn:
    #     yield conn


@pytest.fixture
def postgresql_psycopg2_conn_types(postgresql_psycopg2_engine_types):
    pytest.skip("Skipping MySQL tests")
    # with postgresql_psycopg2_engine_types.connect() as conn:
    #     yield conn


@pytest.fixture
def sqlite_str():
    pytest.importorskip("sqlalchemy")
    with tm.ensure_clean() as name:
        yield f"sqlite:///{name}"


@pytest.fixture
def sqlite_engine(sqlite_str):
    sqlalchemy = pytest.importorskip("sqlalchemy")
    engine = sqlalchemy.create_engine(sqlite_str, poolclass=sqlalchemy.pool.NullPool)
    yield engine
    for view in get_all_views(engine):
        drop_view(view, engine)
    for tbl in get_all_tables(engine):
        drop_table(tbl, engine)
    engine.dispose()


@pytest.fixture
def sqlite_conn(sqlite_engine):
    with sqlite_engine.connect() as conn:
        yield conn


@pytest.fixture
def sqlite_str_iris(sqlite_str, iris_path):
    sqlalchemy = pytest.importorskip("sqlalchemy")
    engine = sqlalchemy.create_engine(sqlite_str)
    create_and_load_iris(engine, iris_path)
    create_and_load_iris_view(engine)
    engine.dispose()
    return sqlite_str


@pytest.fixture
def sqlite_engine_iris(sqlite_engine, iris_path):
    create_and_load_iris(sqlite_engine, iris_path)
    create_and_load_iris_view(sqlite_engine)
    return sqlite_engine


@pytest.fixture
def sqlite_conn_iris(sqlite_engine_iris):
    with sqlite_engine_iris.connect() as conn:
        yield conn


@pytest.fixture
def sqlite_str_types(sqlite_str, types_data):
    sqlalchemy = pytest.importorskip("sqlalchemy")
    engine = sqlalchemy.create_engine(sqlite_str)
    create_and_load_types(engine, types_data, "sqlite")
    engine.dispose()
    return sqlite_str


@pytest.fixture
def sqlite_engine_types(sqlite_engine, types_data):
    create_and_load_types(sqlite_engine, types_data, "sqlite")
    return sqlite_engine


@pytest.fixture
def sqlite_conn_types(sqlite_engine_types):
    with sqlite_engine_types.connect() as conn:
        yield conn


@pytest.fixture
def sqlite_adbc_conn():
    pytest.skip("Skipping MySQL tests")
    # pytest.importorskip("pyarrow")
    # pytest.importorskip("adbc_driver_sqlite")
    # from adbc_driver_sqlite import dbapi

    # with tm.ensure_clean() as name:
    #     uri = f"file:{name}"
    #     with dbapi.connect(uri) as conn:
    #         yield conn
    #         for view in get_all_views(conn):
    #             drop_view(view, conn)
    #         for tbl in get_all_tables(conn):
    #             drop_table(tbl, conn)
    #         conn.commit()


@pytest.fixture
def sqlite_adbc_iris(sqlite_adbc_conn, iris_path):
    pytest.skip("Skipping MySQL tests")
    # import adbc_driver_manager as mgr

    # conn = sqlite_adbc_conn
    # try:
    #     conn.adbc_get_table_schema("iris")
    # except mgr.ProgrammingError:
    #     conn.rollback()
    #     create_and_load_iris_sqlite3(conn, iris_path)
    # try:
    #     conn.adbc_get_table_schema("iris_view")
    # except mgr.ProgrammingError:
    #     conn.rollback()
    #     create_and_load_iris_view(conn)
    # return conn


@pytest.fixture
def sqlite_adbc_types(sqlite_adbc_conn, types_data):
    pytest.skip("Skipping MySQL tests")
    # import adbc_driver_manager as mgr

    # conn = sqlite_adbc_conn
    # try:
    #     conn.adbc_get_table_schema("types")
    # except mgr.ProgrammingError:
    #     conn.rollback()
    #     new_data = []
    #     for entry in types_data:
    #         entry["BoolCol"] = int(entry["BoolCol"])
    #         if entry["BoolColWithNull"] is not None:
    #             entry["BoolColWithNull"] = int(entry["BoolColWithNull"])
    #         new_data.append(tuple(entry.values()))

    #     create_and_load_types_sqlite3(conn, new_data)
    #     conn.commit()

    # return conn


@pytest.fixture
def sqlite_buildin():
    with contextlib.closing(sqlite3.connect(":memory:")) as closing_conn:
        with closing_conn as conn:
            yield conn


@pytest.fixture
def sqlite_buildin_iris(sqlite_buildin, iris_path):
    create_and_load_iris_sqlite3(sqlite_buildin, iris_path)
    create_and_load_iris_view(sqlite_buildin)
    return sqlite_buildin


@pytest.fixture
def sqlite_buildin_types(sqlite_buildin, types_data):
    types_data = [tuple(entry.values()) for entry in types_data]
    create_and_load_types_sqlite3(sqlite_buildin, types_data)
    return sqlite_buildin


mysql_connectable = [
    pytest.param("mysql_pymysql_engine", marks=pytest.mark.db),
    pytest.param("mysql_pymysql_conn", marks=pytest.mark.db),
]

mysql_connectable_iris = [
    pytest.param("mysql_pymysql_engine_iris", marks=pytest.mark.db),
    pytest.param("mysql_pymysql_conn_iris", marks=pytest.mark.db),
]

mysql_connectable_types = [
    pytest.param("mysql_pymysql_engine_types", marks=pytest.mark.db),
    pytest.param("mysql_pymysql_conn_types", marks=pytest.mark.db),
]

postgresql_connectable = [
    pytest.param("postgresql_psycopg2_engine", marks=pytest.mark.db),
    pytest.param("postgresql_psycopg2_conn", marks=pytest.mark.db),
]

postgresql_connectable_iris = [
    pytest.param("postgresql_psycopg2_engine_iris", marks=pytest.mark.db),
    pytest.param("postgresql_psycopg2_conn_iris", marks=pytest.mark.db),
]

postgresql_connectable_types = [
    pytest.param("postgresql_psycopg2_engine_types", marks=pytest.mark.db),
    pytest.param("postgresql_psycopg2_conn_types", marks=pytest.mark.db),
]

sqlite_connectable = [
    "sqlite_engine",
    "sqlite_conn",
    "sqlite_str",
]

sqlite_connectable_iris = [
    "sqlite_engine_iris",
    "sqlite_conn_iris",
    "sqlite_str_iris",
]

sqlite_connectable_types = [
    "sqlite_engine_types",
    "sqlite_conn_types",
    "sqlite_str_types",
]

sqlalchemy_connectable = mysql_connectable + postgresql_connectable + sqlite_connectable

sqlalchemy_connectable_iris = (
    mysql_connectable_iris + postgresql_connectable_iris + sqlite_connectable_iris
)

sqlalchemy_connectable_types = (
    mysql_connectable_types + postgresql_connectable_types + sqlite_connectable_types
)

adbc_connectable = [
    "sqlite_adbc_conn",
    pytest.param("postgresql_adbc_conn", marks=pytest.mark.db),
]

adbc_connectable_iris = [
    pytest.param("postgresql_adbc_iris", marks=pytest.mark.db),
    "sqlite_adbc_iris",
]

adbc_connectable_types = [
    pytest.param("postgresql_adbc_types", marks=pytest.mark.db),
    "sqlite_adbc_types",
]


all_connectable = sqlalchemy_connectable + ["sqlite_buildin"] + adbc_connectable

all_connectable_iris = (
    sqlalchemy_connectable_iris + ["sqlite_buildin_iris"] + adbc_connectable_iris
)

all_connectable_types = (
    sqlalchemy_connectable_types + ["sqlite_buildin_types"] + adbc_connectable_types
)


@pytest.mark.parametrize("conn", all_connectable)
def test_dataframe_to_sql(conn, test_frame1, request):
    # GH 51086 if conn is sqlite_engine
    conn = request.getfixturevalue(conn)
    test_frame1.to_sql(name="test", con=conn, if_exists="append", index=False)


@pytest.mark.parametrize("conn", all_connectable)
def test_dataframe_to_sql_empty(conn, test_frame1, request):
    if conn == "postgresql_adbc_conn" and not using_string_dtype():
        request.node.add_marker(
            pytest.mark.xfail(
                reason="postgres ADBC driver < 1.2 cannot insert index with null type",
            )
        )

    # GH 51086 if conn is sqlite_engine
    conn = request.getfixturevalue(conn)
    empty_df = test_frame1.iloc[:0]
    empty_df.to_sql(name="test", con=conn, if_exists="append", index=False)


@pytest.mark.parametrize("conn", all_connectable)
def test_dataframe_to_sql_arrow_dtypes(conn, request):
    # GH 52046
    pytest.importorskip("pyarrow")
    df = DataFrame(
        {
            "int": pd.array([1], dtype="int8[pyarrow]"),
            "datetime": pd.array(
                [datetime(2023, 1, 1)], dtype="timestamp[ns][pyarrow]"
            ),
            "date": pd.array([date(2023, 1, 1)], dtype="date32[day][pyarrow]"),
            "timedelta": pd.array([timedelta(1)], dtype="duration[ns][pyarrow]"),
            "string": pd.array(["a"], dtype="string[pyarrow]"),
        }
    )

    if "adbc" in conn:
        if conn == "sqlite_adbc_conn":
            df = df.drop(columns=["timedelta"])
        if pa_version_under14p1:
            exp_warning = DeprecationWarning
            msg = "is_sparse is deprecated"
        else:
            exp_warning = None
            msg = ""
    else:
        exp_warning = UserWarning
        msg = "the 'timedelta'"

    conn = request.getfixturevalue(conn)
    with tm.assert_produces_warning(exp_warning, match=msg, check_stacklevel=False):
        df.to_sql(name="test_arrow", con=conn, if_exists="replace", index=False)


@pytest.mark.parametrize("conn", all_connectable)
def test_dataframe_to_sql_arrow_dtypes_missing(conn, request, nulls_fixture):
    # GH 52046
    pytest.importorskip("pyarrow")
    df = DataFrame(
        {
            "datetime": pd.array(
                [datetime(2023, 1, 1), nulls_fixture], dtype="timestamp[ns][pyarrow]"
            ),
        }
    )
    conn = request.getfixturevalue(conn)
    df.to_sql(name="test_arrow", con=conn, if_exists="replace", index=False)


@pytest.mark.parametrize("conn", all_connectable)
@pytest.mark.parametrize("method", [None, "multi"])
def test_to_sql(conn, method, test_frame1, request):
    if method == "multi" and "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'method' not implemented for ADBC drivers", strict=True
            )
        )

    conn = request.getfixturevalue(conn)
    with pandasSQL_builder(conn, need_transaction=True) as pandasSQL:
        pandasSQL.to_sql(test_frame1, "test_frame", method=method)
        assert pandasSQL.has_table("test_frame")
    assert count_rows(conn, "test_frame") == len(test_frame1)


@pytest.mark.parametrize("conn", all_connectable)
@pytest.mark.parametrize(
    "mode, num_row_coef", [("replace", 1), ("append", 2), ("delete_rows", 1)]
)
def test_to_sql_exist(conn, mode, num_row_coef, test_frame1, request):
    conn = request.getfixturevalue(conn)
    with pandasSQL_builder(conn, need_transaction=True) as pandasSQL:
        pandasSQL.to_sql(test_frame1, "test_frame", if_exists="fail")
        pandasSQL.to_sql(test_frame1, "test_frame", if_exists=mode)
        assert pandasSQL.has_table("test_frame")
    assert count_rows(conn, "test_frame") == num_row_coef * len(test_frame1)


@pytest.mark.parametrize("conn", all_connectable)
def test_to_sql_exist_fail(conn, test_frame1, request):
    conn = request.getfixturevalue(conn)
    with pandasSQL_builder(conn, need_transaction=True) as pandasSQL:
        pandasSQL.to_sql(test_frame1, "test_frame", if_exists="fail")
        assert pandasSQL.has_table("test_frame")

        msg = "Table 'test_frame' already exists"
        with pytest.raises(ValueError, match=msg):
            pandasSQL.to_sql(test_frame1, "test_frame", if_exists="fail")


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_read_iris_query(conn, request):
    conn = request.getfixturevalue(conn)
    iris_frame = read_sql_query("SELECT * FROM iris", conn)
    check_iris_frame(iris_frame)
    iris_frame = pd.read_sql("SELECT * FROM iris", conn)
    check_iris_frame(iris_frame)
    iris_frame = pd.read_sql("SELECT * FROM iris where 0=1", conn)
    assert iris_frame.shape == (0, 5)
    assert "SepalWidth" in iris_frame.columns


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_read_iris_query_chunksize(conn, request):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'chunksize' not implemented for ADBC drivers",
                strict=True,
            )
        )
    conn = request.getfixturevalue(conn)
    iris_frame = concat(read_sql_query("SELECT * FROM iris", conn, chunksize=7))
    check_iris_frame(iris_frame)
    iris_frame = concat(pd.read_sql("SELECT * FROM iris", conn, chunksize=7))
    check_iris_frame(iris_frame)
    iris_frame = concat(pd.read_sql("SELECT * FROM iris where 0=1", conn, chunksize=7))
    assert iris_frame.shape == (0, 5)
    assert "SepalWidth" in iris_frame.columns


@pytest.mark.parametrize("conn", sqlalchemy_connectable_iris)
def test_read_iris_query_expression_with_parameter(conn, request):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'chunksize' not implemented for ADBC drivers",
                strict=True,
            )
        )
    conn = request.getfixturevalue(conn)
    from sqlalchemy import (
        MetaData,
        Table,
        create_engine,
        select,
    )

    metadata = MetaData()
    autoload_con = create_engine(conn) if isinstance(conn, str) else conn
    iris = Table("iris", metadata, autoload_with=autoload_con)
    iris_frame = read_sql_query(
        select(iris), conn, params={"name": "Iris-setosa", "length": 5.1}
    )
    check_iris_frame(iris_frame)
    if isinstance(conn, str):
        autoload_con.dispose()


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_read_iris_query_string_with_parameter(conn, request, sql_strings):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'chunksize' not implemented for ADBC drivers",
                strict=True,
            )
        )

    for db, query in sql_strings["read_parameters"].items():
        if db in conn:
            break
    else:
        raise KeyError(f"No part of {conn} found in sql_strings['read_parameters']")
    conn = request.getfixturevalue(conn)
    iris_frame = read_sql_query(query, conn, params=("Iris-setosa", 5.1))
    check_iris_frame(iris_frame)


@pytest.mark.parametrize("conn", sqlalchemy_connectable_iris)
def test_read_iris_table(conn, request):
    # GH 51015 if conn = sqlite_iris_str
    conn = request.getfixturevalue(conn)
    iris_frame = read_sql_table("iris", conn)
    check_iris_frame(iris_frame)
    iris_frame = pd.read_sql("iris", conn)
    check_iris_frame(iris_frame)


@pytest.mark.parametrize("conn", sqlalchemy_connectable_iris)
def test_read_iris_table_chunksize(conn, request):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(reason="chunksize argument NotImplemented with ADBC")
        )
    conn = request.getfixturevalue(conn)
    iris_frame = concat(read_sql_table("iris", conn, chunksize=7))
    check_iris_frame(iris_frame)
    iris_frame = concat(pd.read_sql("iris", conn, chunksize=7))
    check_iris_frame(iris_frame)


@pytest.mark.parametrize("conn", sqlalchemy_connectable)
def test_to_sql_callable(conn, test_frame1, request):
    conn = request.getfixturevalue(conn)

    check = []  # used to double check function below is really being used

    def sample(pd_table, conn, keys, data_iter):
        check.append(1)
        data = [dict(zip(keys, row)) for row in data_iter]
        conn.execute(pd_table.table.insert(), data)

    with pandasSQL_builder(conn, need_transaction=True) as pandasSQL:
        pandasSQL.to_sql(test_frame1, "test_frame", method=sample)
        assert pandasSQL.has_table("test_frame")
    assert check == [1]
    assert count_rows(conn, "test_frame") == len(test_frame1)


@pytest.mark.parametrize("conn", all_connectable_types)
def test_default_type_conversion(conn, request):
    conn_name = conn
    if conn_name == "sqlite_buildin_types":
        request.applymarker(
            pytest.mark.xfail(
                reason="sqlite_buildin connection does not implement read_sql_table"
            )
        )

    conn = request.getfixturevalue(conn)
    df = sql.read_sql_table("types", conn)

    assert issubclass(df.FloatCol.dtype.type, np.floating)
    assert issubclass(df.IntCol.dtype.type, np.integer)

    # MySQL/sqlite has no real BOOL type
    if "postgresql" in conn_name:
        assert issubclass(df.BoolCol.dtype.type, np.bool_)
    else:
        assert issubclass(df.BoolCol.dtype.type, np.integer)

    # Int column with NA values stays as float
    assert issubclass(df.IntColWithNull.dtype.type, np.floating)

    # Bool column with NA = int column with NA values => becomes float
    if "postgresql" in conn_name:
        assert issubclass(df.BoolColWithNull.dtype.type, object)
    else:
        assert issubclass(df.BoolColWithNull.dtype.type, np.floating)


@pytest.mark.parametrize("conn", mysql_connectable)
def test_read_procedure(conn, request):
    conn = request.getfixturevalue(conn)

    # GH 7324
    # Although it is more an api test, it is added to the
    # mysql tests as sqlite does not have stored procedures
    from sqlalchemy import text
    from sqlalchemy.engine import Engine

    df = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3]})
    df.to_sql(name="test_frame", con=conn, index=False)

    proc = """DROP PROCEDURE IF EXISTS get_testdb;

    CREATE PROCEDURE get_testdb ()

    BEGIN
        SELECT * FROM test_frame;
    END"""
    proc = text(proc)
    if isinstance(conn, Engine):
        with conn.connect() as engine_conn:
            with engine_conn.begin():
                engine_conn.execute(proc)
    else:
        with conn.begin():
            conn.execute(proc)

    res1 = sql.read_sql_query("CALL get_testdb();", conn)
    tm.assert_frame_equal(df, res1)

    # test delegation to read_sql_query
    res2 = sql.read_sql("CALL get_testdb();", conn)
    tm.assert_frame_equal(df, res2)


@pytest.mark.parametrize("conn", postgresql_connectable)
@pytest.mark.parametrize("expected_count", [2, "Success!"])
def test_copy_from_callable_insertion_method(conn, expected_count, request):
    # GH 8953
    # Example in io.rst found under _io.sql.method
    # not available in sqlite, mysql
    def psql_insert_copy(table, conn, keys, data_iter):
        # gets a DBAPI connection that can provide a cursor
        dbapi_conn = conn.connection
        with dbapi_conn.cursor() as cur:
            s_buf = StringIO()
            writer = csv.writer(s_buf)
            writer.writerows(data_iter)
            s_buf.seek(0)

            columns = ", ".join([f'"{k}"' for k in keys])
            if table.schema:
                table_name = f"{table.schema}.{table.name}"
            else:
                table_name = table.name

            sql_query = f"COPY {table_name} ({columns}) FROM STDIN WITH CSV"
            cur.copy_expert(sql=sql_query, file=s_buf)
        return expected_count

    conn = request.getfixturevalue(conn)
    expected = DataFrame({"col1": [1, 2], "col2": [0.1, 0.2], "col3": ["a", "n"]})
    result_count = expected.to_sql(
        name="test_frame", con=conn, index=False, method=psql_insert_copy
    )
    # GH 46891
    if expected_count is None:
        assert result_count is None
    else:
        assert result_count == expected_count
    result = sql.read_sql_table("test_frame", conn)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("conn", postgresql_connectable)
def test_insertion_method_on_conflict_do_nothing(conn, request):
    # GH 15988: Example in to_sql docstring
    conn = request.getfixturevalue(conn)

    from sqlalchemy.dialects.postgresql import insert
    from sqlalchemy.engine import Engine
    from sqlalchemy.sql import text

    def insert_on_conflict(table, conn, keys, data_iter):
        data = [dict(zip(keys, row)) for row in data_iter]
        stmt = (
            insert(table.table)
            .values(data)
            .on_conflict_do_nothing(index_elements=["a"])
        )
        result = conn.execute(stmt)
        return result.rowcount

    create_sql = text(
        """
    CREATE TABLE test_insert_conflict (
        a  integer PRIMARY KEY,
        b  numeric,
        c  text
    );
    """
    )
    if isinstance(conn, Engine):
        with conn.connect() as con:
            with con.begin():
                con.execute(create_sql)
    else:
        with conn.begin():
            conn.execute(create_sql)

    expected = DataFrame([[1, 2.1, "a"]], columns=list("abc"))
    expected.to_sql(
        name="test_insert_conflict", con=conn, if_exists="append", index=False
    )

    df_insert = DataFrame([[1, 3.2, "b"]], columns=list("abc"))
    inserted = df_insert.to_sql(
        name="test_insert_conflict",
        con=conn,
        index=False,
        if_exists="append",
        method=insert_on_conflict,
    )
    result = sql.read_sql_table("test_insert_conflict", conn)
    tm.assert_frame_equal(result, expected)
    assert inserted == 0

    # Cleanup
    with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
        pandasSQL.drop_table("test_insert_conflict")


@pytest.mark.parametrize("conn", all_connectable)
def test_to_sql_on_public_schema(conn, request):
    if "sqlite" in conn or "mysql" in conn:
        request.applymarker(
            pytest.mark.xfail(
                reason="test for public schema only specific to postgresql"
            )
        )

    conn = request.getfixturevalue(conn)

    test_data = DataFrame([[1, 2.1, "a"], [2, 3.1, "b"]], columns=list("abc"))
    test_data.to_sql(
        name="test_public_schema",
        con=conn,
        if_exists="append",
        index=False,
        schema="public",
    )

    df_out = sql.read_sql_table("test_public_schema", conn, schema="public")
    tm.assert_frame_equal(test_data, df_out)


@pytest.mark.parametrize("conn", mysql_connectable)
def test_insertion_method_on_conflict_update(conn, request):
    # GH 14553: Example in to_sql docstring
    conn = request.getfixturevalue(conn)

    from sqlalchemy.dialects.mysql import insert
    from sqlalchemy.engine import Engine
    from sqlalchemy.sql import text

    def insert_on_conflict(table, conn, keys, data_iter):
        data = [dict(zip(keys, row)) for row in data_iter]
        stmt = insert(table.table).values(data)
        stmt = stmt.on_duplicate_key_update(b=stmt.inserted.b, c=stmt.inserted.c)
        result = conn.execute(stmt)
        return result.rowcount

    create_sql = text(
        """
    CREATE TABLE test_insert_conflict (
        a INT PRIMARY KEY,
        b FLOAT,
        c VARCHAR(10)
    );
    """
    )
    if isinstance(conn, Engine):
        with conn.connect() as con:
            with con.begin():
                con.execute(create_sql)
    else:
        with conn.begin():
            conn.execute(create_sql)

    df = DataFrame([[1, 2.1, "a"]], columns=list("abc"))
    df.to_sql(name="test_insert_conflict", con=conn, if_exists="append", index=False)

    expected = DataFrame([[1, 3.2, "b"]], columns=list("abc"))
    inserted = expected.to_sql(
        name="test_insert_conflict",
        con=conn,
        index=False,
        if_exists="append",
        method=insert_on_conflict,
    )
    result = sql.read_sql_table("test_insert_conflict", conn)
    tm.assert_frame_equal(result, expected)
    assert inserted == 2

    # Cleanup
    with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
        pandasSQL.drop_table("test_insert_conflict")


@pytest.mark.parametrize("conn", postgresql_connectable)
def test_read_view_postgres(conn, request):
    # GH 52969
    conn = request.getfixturevalue(conn)

    from sqlalchemy.engine import Engine
    from sqlalchemy.sql import text

    table_name = f"group_{uuid.uuid4().hex}"
    view_name = f"group_view_{uuid.uuid4().hex}"

    sql_stmt = text(
        f"""
    CREATE TABLE {table_name} (
        group_id INTEGER,
        name TEXT
    );
    INSERT INTO {table_name} VALUES
        (1, 'name');
    CREATE VIEW {view_name}
    AS
    SELECT * FROM {table_name};
    """
    )
    if isinstance(conn, Engine):
        with conn.connect() as con:
            with con.begin():
                con.execute(sql_stmt)
    else:
        with conn.begin():
            conn.execute(sql_stmt)
    result = read_sql_table(view_name, conn)
    expected = DataFrame({"group_id": [1], "name": "name"})
    tm.assert_frame_equal(result, expected)


def test_read_view_sqlite(sqlite_buildin):
    # GH 52969
    create_table = """
CREATE TABLE groups (
   group_id INTEGER,
   name TEXT
);
"""
    insert_into = """
INSERT INTO groups VALUES
    (1, 'name');
"""
    create_view = """
CREATE VIEW group_view
AS
SELECT * FROM groups;
"""
    sqlite_buildin.execute(create_table)
    sqlite_buildin.execute(insert_into)
    sqlite_buildin.execute(create_view)
    result = pd.read_sql("SELECT * FROM group_view", sqlite_buildin)
    expected = DataFrame({"group_id": [1], "name": "name"})
    tm.assert_frame_equal(result, expected)


def flavor(conn_name):
    if "postgresql" in conn_name:
        return "postgresql"
    elif "sqlite" in conn_name:
        return "sqlite"
    elif "mysql" in conn_name:
        return "mysql"

    raise ValueError(f"unsupported connection: {conn_name}")


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_read_sql_iris_parameter(conn, request, sql_strings):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'params' not implemented for ADBC drivers",
                strict=True,
            )
        )
    conn_name = conn
    conn = request.getfixturevalue(conn)
    query = sql_strings["read_parameters"][flavor(conn_name)]
    params = ("Iris-setosa", 5.1)
    with pandasSQL_builder(conn) as pandasSQL:
        with pandasSQL.run_transaction():
            iris_frame = pandasSQL.read_query(query, params=params)
    check_iris_frame(iris_frame)


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_read_sql_iris_named_parameter(conn, request, sql_strings):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'params' not implemented for ADBC drivers",
                strict=True,
            )
        )

    conn_name = conn
    conn = request.getfixturevalue(conn)
    query = sql_strings["read_named_parameters"][flavor(conn_name)]
    params = {"name": "Iris-setosa", "length": 5.1}
    with pandasSQL_builder(conn) as pandasSQL:
        with pandasSQL.run_transaction():
            iris_frame = pandasSQL.read_query(query, params=params)
    check_iris_frame(iris_frame)


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_read_sql_iris_no_parameter_with_percent(conn, request, sql_strings):
    if "mysql" in conn or ("postgresql" in conn and "adbc" not in conn):
        request.applymarker(pytest.mark.xfail(reason="broken test"))

    conn_name = conn
    conn = request.getfixturevalue(conn)

    query = sql_strings["read_no_parameters_with_percent"][flavor(conn_name)]
    with pandasSQL_builder(conn) as pandasSQL:
        with pandasSQL.run_transaction():
            iris_frame = pandasSQL.read_query(query, params=None)
    check_iris_frame(iris_frame)


# -----------------------------------------------------------------------------
# -- Testing the public API


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_api_read_sql_view(conn, request):
    conn = request.getfixturevalue(conn)
    iris_frame = sql.read_sql_query("SELECT * FROM iris_view", conn)
    check_iris_frame(iris_frame)


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_api_read_sql_with_chunksize_no_result(conn, request):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(reason="chunksize argument NotImplemented with ADBC")
        )
    conn = request.getfixturevalue(conn)
    query = 'SELECT * FROM iris_view WHERE "SepalLength" < 0.0'
    with_batch = sql.read_sql_query(query, conn, chunksize=5)
    without_batch = sql.read_sql_query(query, conn)
    tm.assert_frame_equal(concat(with_batch), without_batch)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_to_sql(conn, request, test_frame1):
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_frame1", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_frame1")

    sql.to_sql(test_frame1, "test_frame1", conn)
    assert sql.has_table("test_frame1", conn)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_to_sql_fail(conn, request, test_frame1):
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_frame2", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_frame2")

    sql.to_sql(test_frame1, "test_frame2", conn, if_exists="fail")
    assert sql.has_table("test_frame2", conn)

    msg = "Table 'test_frame2' already exists"
    with pytest.raises(ValueError, match=msg):
        sql.to_sql(test_frame1, "test_frame2", conn, if_exists="fail")


@pytest.mark.parametrize("conn", all_connectable)
def test_api_to_sql_replace(conn, request, test_frame1):
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_frame3", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_frame3")

    sql.to_sql(test_frame1, "test_frame3", conn, if_exists="fail")
    # Add to table again
    sql.to_sql(test_frame1, "test_frame3", conn, if_exists="replace")
    assert sql.has_table("test_frame3", conn)

    num_entries = len(test_frame1)
    num_rows = count_rows(conn, "test_frame3")

    assert num_rows == num_entries


@pytest.mark.parametrize("conn", all_connectable)
def test_api_to_sql_append(conn, request, test_frame1):
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_frame4", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_frame4")

    assert sql.to_sql(test_frame1, "test_frame4", conn, if_exists="fail") == 4

    # Add to table again
    assert sql.to_sql(test_frame1, "test_frame4", conn, if_exists="append") == 4
    assert sql.has_table("test_frame4", conn)

    num_entries = 2 * len(test_frame1)
    num_rows = count_rows(conn, "test_frame4")

    assert num_rows == num_entries


@pytest.mark.parametrize("conn", all_connectable)
def test_api_to_sql_type_mapping(conn, request, test_frame3):
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_frame5", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_frame5")

    sql.to_sql(test_frame3, "test_frame5", conn, index=False)
    result = sql.read_sql("SELECT * FROM test_frame5", conn)

    tm.assert_frame_equal(test_frame3, result)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_to_sql_series(conn, request):
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_series", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_series")

    s = Series(np.arange(5, dtype="int64"), name="series")
    sql.to_sql(s, "test_series", conn, index=False)
    s2 = sql.read_sql_query("SELECT * FROM test_series", conn)
    tm.assert_frame_equal(s.to_frame(), s2)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_roundtrip(conn, request, test_frame1):
    conn_name = conn
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_frame_roundtrip", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_frame_roundtrip")

    sql.to_sql(test_frame1, "test_frame_roundtrip", con=conn)
    result = sql.read_sql_query("SELECT * FROM test_frame_roundtrip", con=conn)

    # HACK!
    if "adbc" in conn_name:
        result = result.drop(columns="__index_level_0__")
    else:
        result = result.drop(columns="level_0")
    tm.assert_frame_equal(result, test_frame1)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_roundtrip_chunksize(conn, request, test_frame1):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(reason="chunksize argument NotImplemented with ADBC")
        )
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_frame_roundtrip", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_frame_roundtrip")

    sql.to_sql(
        test_frame1,
        "test_frame_roundtrip",
        con=conn,
        index=False,
        chunksize=2,
    )
    result = sql.read_sql_query("SELECT * FROM test_frame_roundtrip", con=conn)
    tm.assert_frame_equal(result, test_frame1)


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_api_execute_sql(conn, request):
    # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
    conn = request.getfixturevalue(conn)
    with sql.pandasSQL_builder(conn) as pandas_sql:
        iris_results = pandas_sql.execute("SELECT * FROM iris")
        row = iris_results.fetchone()
        iris_results.close()
    assert list(row) == [5.1, 3.5, 1.4, 0.2, "Iris-setosa"]


@pytest.mark.parametrize("conn", all_connectable_types)
def test_api_date_parsing(conn, request):
    conn_name = conn
    conn = request.getfixturevalue(conn)
    # Test date parsing in read_sql
    # No Parsing
    df = sql.read_sql_query("SELECT * FROM types", conn)
    if not ("mysql" in conn_name or "postgres" in conn_name):
        assert not issubclass(df.DateCol.dtype.type, np.datetime64)

    df = sql.read_sql_query("SELECT * FROM types", conn, parse_dates=["DateCol"])
    assert issubclass(df.DateCol.dtype.type, np.datetime64)
    assert df.DateCol.tolist() == [
        Timestamp(2000, 1, 3, 0, 0, 0),
        Timestamp(2000, 1, 4, 0, 0, 0),
    ]

    df = sql.read_sql_query(
        "SELECT * FROM types",
        conn,
        parse_dates={"DateCol": "%Y-%m-%d %H:%M:%S"},
    )
    assert issubclass(df.DateCol.dtype.type, np.datetime64)
    assert df.DateCol.tolist() == [
        Timestamp(2000, 1, 3, 0, 0, 0),
        Timestamp(2000, 1, 4, 0, 0, 0),
    ]

    df = sql.read_sql_query("SELECT * FROM types", conn, parse_dates=["IntDateCol"])
    assert issubclass(df.IntDateCol.dtype.type, np.datetime64)
    assert df.IntDateCol.tolist() == [
        Timestamp(1986, 12, 25, 0, 0, 0),
        Timestamp(2013, 1, 1, 0, 0, 0),
    ]

    df = sql.read_sql_query(
        "SELECT * FROM types", conn, parse_dates={"IntDateCol": "s"}
    )
    assert issubclass(df.IntDateCol.dtype.type, np.datetime64)
    assert df.IntDateCol.tolist() == [
        Timestamp(1986, 12, 25, 0, 0, 0),
        Timestamp(2013, 1, 1, 0, 0, 0),
    ]

    df = sql.read_sql_query(
        "SELECT * FROM types",
        conn,
        parse_dates={"IntDateOnlyCol": "%Y%m%d"},
    )
    assert issubclass(df.IntDateOnlyCol.dtype.type, np.datetime64)
    assert df.IntDateOnlyCol.tolist() == [
        Timestamp("2010-10-10"),
        Timestamp("2010-12-12"),
    ]


@pytest.mark.parametrize("conn", all_connectable_types)
@pytest.mark.parametrize("error", ["raise", "coerce"])
@pytest.mark.parametrize(
    "read_sql, text, mode",
    [
        (sql.read_sql, "SELECT * FROM types", ("sqlalchemy", "fallback")),
        (sql.read_sql, "types", ("sqlalchemy")),
        (
            sql.read_sql_query,
            "SELECT * FROM types",
            ("sqlalchemy", "fallback"),
        ),
        (sql.read_sql_table, "types", ("sqlalchemy")),
    ],
)
def test_api_custom_dateparsing_error(
    conn, request, read_sql, text, mode, error, types_data_frame
):
    conn_name = conn
    conn = request.getfixturevalue(conn)
    if text == "types" and conn_name == "sqlite_buildin_types":
        request.applymarker(
            pytest.mark.xfail(reason="failing combination of arguments")
        )

    expected = types_data_frame.astype({"DateCol": "datetime64[s]"})

    result = read_sql(
        text,
        con=conn,
        parse_dates={
            "DateCol": {"errors": error},
        },
    )
    if "postgres" in conn_name:
        # TODO: clean up types_data_frame fixture
        result["BoolCol"] = result["BoolCol"].astype(int)
        result["BoolColWithNull"] = result["BoolColWithNull"].astype(float)

    if conn_name == "postgresql_adbc_types":
        expected = expected.astype(
            {
                "IntDateCol": "int32",
                "IntDateOnlyCol": "int32",
                "IntCol": "int32",
            }
        )

    if conn_name == "postgresql_adbc_types" and pa_version_under14p1:
        expected["DateCol"] = expected["DateCol"].astype("datetime64[ns]")
    elif "postgres" in conn_name or "mysql" in conn_name:
        expected["DateCol"] = expected["DateCol"].astype("datetime64[us]")
    else:
        expected["DateCol"] = expected["DateCol"].astype("datetime64[s]")
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("conn", all_connectable_types)
def test_api_date_and_index(conn, request):
    # Test case where same column appears in parse_date and index_col
    conn = request.getfixturevalue(conn)
    df = sql.read_sql_query(
        "SELECT * FROM types",
        conn,
        index_col="DateCol",
        parse_dates=["DateCol", "IntDateCol"],
    )

    assert issubclass(df.index.dtype.type, np.datetime64)
    assert issubclass(df.IntDateCol.dtype.type, np.datetime64)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_timedelta(conn, request):
    # see #6921
    conn_name = conn
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_timedelta", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_timedelta")

    df = to_timedelta(Series(["00:00:01", "00:00:03"], name="foo")).to_frame()

    if conn_name == "sqlite_adbc_conn":
        request.node.add_marker(
            pytest.mark.xfail(
                reason="sqlite ADBC driver doesn't implement timedelta",
            )
        )

    if "adbc" in conn_name:
        if pa_version_under14p1:
            exp_warning = DeprecationWarning
        else:
            exp_warning = None
    else:
        exp_warning = UserWarning

    with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
        result_count = df.to_sql(name="test_timedelta", con=conn)
    assert result_count == 2
    result = sql.read_sql_query("SELECT * FROM test_timedelta", conn)

    if conn_name == "postgresql_adbc_conn":
        # TODO: Postgres stores an INTERVAL, which ADBC reads as a Month-Day-Nano
        # Interval; the default pandas type mapper maps this to a DateOffset
        # but maybe we should try and restore the timedelta here?
        expected = Series(
            [
                pd.DateOffset(months=0, days=0, microseconds=1000000, nanoseconds=0),
                pd.DateOffset(months=0, days=0, microseconds=3000000, nanoseconds=0),
            ],
            name="foo",
        )
    else:
        expected = df["foo"].astype("int64")
    tm.assert_series_equal(result["foo"], expected)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_complex_raises(conn, request):
    conn_name = conn
    conn = request.getfixturevalue(conn)
    df = DataFrame({"a": [1 + 1j, 2j]})

    if "adbc" in conn_name:
        msg = "datatypes not supported"
    else:
        msg = "Complex datatypes not supported"
    with pytest.raises(ValueError, match=msg):
        assert df.to_sql("test_complex", con=conn) is None


@pytest.mark.parametrize("conn", all_connectable)
@pytest.mark.parametrize(
    "index_name,index_label,expected",
    [
        # no index name, defaults to 'index'
        (None, None, "index"),
        # specifying index_label
        (None, "other_label", "other_label"),
        # using the index name
        ("index_name", None, "index_name"),
        # has index name, but specifying index_label
        ("index_name", "other_label", "other_label"),
        # index name is integer
        (0, None, "0"),
        # index name is None but index label is integer
        (None, 0, "0"),
    ],
)
def test_api_to_sql_index_label(conn, request, index_name, index_label, expected):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(reason="index_label argument NotImplemented with ADBC")
        )
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_index_label", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_index_label")

    temp_frame = DataFrame({"col1": range(4)})
    temp_frame.index.name = index_name
    query = "SELECT * FROM test_index_label"
    sql.to_sql(temp_frame, "test_index_label", conn, index_label=index_label)
    frame = sql.read_sql_query(query, conn)
    assert frame.columns[0] == expected


@pytest.mark.parametrize("conn", all_connectable)
def test_api_to_sql_index_label_multiindex(conn, request):
    conn_name = conn
    if "mysql" in conn_name:
        request.applymarker(
            pytest.mark.xfail(
                reason="MySQL can fail using TEXT without length as key", strict=False
            )
        )
    elif "adbc" in conn_name:
        request.node.add_marker(
            pytest.mark.xfail(reason="index_label argument NotImplemented with ADBC")
        )

    conn = request.getfixturevalue(conn)
    if sql.has_table("test_index_label", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_index_label")

    expected_row_count = 4
    temp_frame = DataFrame(
        {"col1": range(4)},
        index=MultiIndex.from_product([("A0", "A1"), ("B0", "B1")]),
    )

    # no index name, defaults to 'level_0' and 'level_1'
    result = sql.to_sql(temp_frame, "test_index_label", conn)
    assert result == expected_row_count
    frame = sql.read_sql_query("SELECT * FROM test_index_label", conn)
    assert frame.columns[0] == "level_0"
    assert frame.columns[1] == "level_1"

    # specifying index_label
    result = sql.to_sql(
        temp_frame,
        "test_index_label",
        conn,
        if_exists="replace",
        index_label=["A", "B"],
    )
    assert result == expected_row_count
    frame = sql.read_sql_query("SELECT * FROM test_index_label", conn)
    assert frame.columns[:2].tolist() == ["A", "B"]

    # using the index name
    temp_frame.index.names = ["A", "B"]
    result = sql.to_sql(temp_frame, "test_index_label", conn, if_exists="replace")
    assert result == expected_row_count
    frame = sql.read_sql_query("SELECT * FROM test_index_label", conn)
    assert frame.columns[:2].tolist() == ["A", "B"]

    # has index name, but specifying index_label
    result = sql.to_sql(
        temp_frame,
        "test_index_label",
        conn,
        if_exists="replace",
        index_label=["C", "D"],
    )
    assert result == expected_row_count
    frame = sql.read_sql_query("SELECT * FROM test_index_label", conn)
    assert frame.columns[:2].tolist() == ["C", "D"]

    msg = "Length of 'index_label' should match number of levels, which is 2"
    with pytest.raises(ValueError, match=msg):
        sql.to_sql(
            temp_frame,
            "test_index_label",
            conn,
            if_exists="replace",
            index_label="C",
        )


@pytest.mark.parametrize("conn", all_connectable)
def test_api_multiindex_roundtrip(conn, request):
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_multiindex_roundtrip", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_multiindex_roundtrip")

    df = DataFrame.from_records(
        [(1, 2.1, "line1"), (2, 1.5, "line2")],
        columns=["A", "B", "C"],
        index=["A", "B"],
    )

    df.to_sql(name="test_multiindex_roundtrip", con=conn)
    result = sql.read_sql_query(
        "SELECT * FROM test_multiindex_roundtrip", conn, index_col=["A", "B"]
    )
    tm.assert_frame_equal(df, result, check_index_type=True)


@pytest.mark.parametrize("conn", all_connectable)
@pytest.mark.parametrize(
    "dtype",
    [
        None,
        int,
        float,
        {"A": int, "B": float},
    ],
)
def test_api_dtype_argument(conn, request, dtype):
    # GH10285 Add dtype argument to read_sql_query
    conn_name = conn
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_dtype_argument", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_dtype_argument")

    df = DataFrame([[1.2, 3.4], [5.6, 7.8]], columns=["A", "B"])
    assert df.to_sql(name="test_dtype_argument", con=conn) == 2

    expected = df.astype(dtype)

    if "postgres" in conn_name:
        query = 'SELECT "A", "B" FROM test_dtype_argument'
    else:
        query = "SELECT A, B FROM test_dtype_argument"
    result = sql.read_sql_query(query, con=conn, dtype=dtype)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_integer_col_names(conn, request):
    conn = request.getfixturevalue(conn)
    df = DataFrame([[1, 2], [3, 4]], columns=[0, 1])
    sql.to_sql(df, "test_frame_integer_col_names", conn, if_exists="replace")


@pytest.mark.parametrize("conn", all_connectable)
def test_api_get_schema(conn, request, test_frame1):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'get_schema' not implemented for ADBC drivers",
                strict=True,
            )
        )
    conn = request.getfixturevalue(conn)
    create_sql = sql.get_schema(test_frame1, "test", con=conn)
    assert "CREATE" in create_sql


@pytest.mark.parametrize("conn", all_connectable)
def test_api_get_schema_with_schema(conn, request, test_frame1):
    # GH28486
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'get_schema' not implemented for ADBC drivers",
                strict=True,
            )
        )
    conn = request.getfixturevalue(conn)
    create_sql = sql.get_schema(test_frame1, "test", con=conn, schema="pypi")
    assert "CREATE TABLE pypi." in create_sql


@pytest.mark.parametrize("conn", all_connectable)
def test_api_get_schema_dtypes(conn, request):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'get_schema' not implemented for ADBC drivers",
                strict=True,
            )
        )
    conn_name = conn
    conn = request.getfixturevalue(conn)
    float_frame = DataFrame({"a": [1.1, 1.2], "b": [2.1, 2.2]})

    if conn_name == "sqlite_buildin":
        dtype = "INTEGER"
    else:
        from sqlalchemy import Integer

        dtype = Integer
    create_sql = sql.get_schema(float_frame, "test", con=conn, dtype={"b": dtype})
    assert "CREATE" in create_sql
    assert "INTEGER" in create_sql


@pytest.mark.parametrize("conn", all_connectable)
def test_api_get_schema_keys(conn, request, test_frame1):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="'get_schema' not implemented for ADBC drivers",
                strict=True,
            )
        )
    conn_name = conn
    conn = request.getfixturevalue(conn)
    frame = DataFrame({"Col1": [1.1, 1.2], "Col2": [2.1, 2.2]})
    create_sql = sql.get_schema(frame, "test", con=conn, keys="Col1")

    if "mysql" in conn_name:
        constraint_sentence = "CONSTRAINT test_pk PRIMARY KEY (`Col1`)"
    else:
        constraint_sentence = 'CONSTRAINT test_pk PRIMARY KEY ("Col1")'
    assert constraint_sentence in create_sql

    # multiple columns as key (GH10385)
    create_sql = sql.get_schema(test_frame1, "test", con=conn, keys=["A", "B"])
    if "mysql" in conn_name:
        constraint_sentence = "CONSTRAINT test_pk PRIMARY KEY (`A`, `B`)"
    else:
        constraint_sentence = 'CONSTRAINT test_pk PRIMARY KEY ("A", "B")'
    assert constraint_sentence in create_sql


@pytest.mark.parametrize("conn", all_connectable)
def test_api_chunksize_read(conn, request):
    if "adbc" in conn:
        request.node.add_marker(
            pytest.mark.xfail(reason="chunksize argument NotImplemented with ADBC")
        )
    conn_name = conn
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_chunksize", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_chunksize")

    df = DataFrame(
        np.random.default_rng(2).standard_normal((22, 5)), columns=list("abcde")
    )
    df.to_sql(name="test_chunksize", con=conn, index=False)

    # reading the query in one time
    res1 = sql.read_sql_query("select * from test_chunksize", conn)

    # reading the query in chunks with read_sql_query
    res2 = DataFrame()
    i = 0
    sizes = [5, 5, 5, 5, 2]

    for chunk in sql.read_sql_query("select * from test_chunksize", conn, chunksize=5):
        res2 = concat([res2, chunk], ignore_index=True)
        assert len(chunk) == sizes[i]
        i += 1

    tm.assert_frame_equal(res1, res2)

    # reading the query in chunks with read_sql_query
    if conn_name == "sqlite_buildin":
        with pytest.raises(NotImplementedError, match=""):
            sql.read_sql_table("test_chunksize", conn, chunksize=5)
    else:
        res3 = DataFrame()
        i = 0
        sizes = [5, 5, 5, 5, 2]

        for chunk in sql.read_sql_table("test_chunksize", conn, chunksize=5):
            res3 = concat([res3, chunk], ignore_index=True)
            assert len(chunk) == sizes[i]
            i += 1

        tm.assert_frame_equal(res1, res3)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_categorical(conn, request):
    if conn == "postgresql_adbc_conn":
        adbc = import_optional_dependency("adbc_driver_postgresql", errors="ignore")
        if adbc is not None and Version(adbc.__version__) < Version("0.9.0"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="categorical dtype not implemented for ADBC postgres driver",
                    strict=True,
                )
            )
    # GH8624
    # test that categorical gets written correctly as dense column
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_categorical", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_categorical")

    df = DataFrame(
        {
            "person_id": [1, 2, 3],
            "person_name": ["John P. Doe", "Jane Dove", "John P. Doe"],
        }
    )
    df2 = df.copy()
    df2["person_name"] = df2["person_name"].astype("category")

    df2.to_sql(name="test_categorical", con=conn, index=False)
    res = sql.read_sql_query("SELECT * FROM test_categorical", conn)

    tm.assert_frame_equal(res, df)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_unicode_column_name(conn, request):
    # GH 11431
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_unicode", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_unicode")

    df = DataFrame([[1, 2], [3, 4]], columns=["\xe9", "b"])
    df.to_sql(name="test_unicode", con=conn, index=False)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_escaped_table_name(conn, request):
    # GH 13206
    conn_name = conn
    conn = request.getfixturevalue(conn)
    if sql.has_table("d1187b08-4943-4c8d-a7f6", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("d1187b08-4943-4c8d-a7f6")

    df = DataFrame({"A": [0, 1, 2], "B": [0.2, np.nan, 5.6]})
    df.to_sql(name="d1187b08-4943-4c8d-a7f6", con=conn, index=False)

    if "postgres" in conn_name:
        query = 'SELECT * FROM "d1187b08-4943-4c8d-a7f6"'
    else:
        query = "SELECT * FROM `d1187b08-4943-4c8d-a7f6`"
    res = sql.read_sql_query(query, conn)

    tm.assert_frame_equal(res, df)


@pytest.mark.parametrize("conn", all_connectable)
def test_api_read_sql_duplicate_columns(conn, request):
    # GH#53117
    if "adbc" in conn:
        pa = pytest.importorskip("pyarrow")
        if not (
            Version(pa.__version__) >= Version("16.0")
            and conn in ["sqlite_adbc_conn", "postgresql_adbc_conn"]
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="pyarrow->pandas throws ValueError", strict=True
                )
            )
    conn = request.getfixturevalue(conn)
    if sql.has_table("test_table", conn):
        with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
            pandasSQL.drop_table("test_table")

    df = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3], "c": 1})
    df.to_sql(name="test_table", con=conn, index=False)

    result = pd.read_sql("SELECT a, b, a +1 as a, c FROM test_table", conn)
    expected = DataFrame(
        [[1, 0.1, 2, 1], [2, 0.2, 3, 1], [3, 0.3, 4, 1]],
        columns=["a", "b", "a", "c"],
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("conn", all_connectable)
def test_read_table_columns(conn, request, test_frame1):
    # test columns argument in read_table
    conn_name = conn
    if conn_name == "sqlite_buildin":
        request.applymarker(pytest.mark.xfail(reason="Not Implemented"))

    conn = request.getfixturevalue(conn)
    sql.to_sql(test_frame1, "test_frame", conn)

    cols = ["A", "B"]

    result = sql.read_sql_table("test_frame", conn, columns=cols)
    assert result.columns.tolist() == cols


@pytest.mark.parametrize("conn", all_connectable)
def test_read_table_index_col(conn, request, test_frame1):
    # test columns argument in read_table
    conn_name = conn
    if conn_name == "sqlite_buildin":
        request.applymarker(pytest.mark.xfail(reason="Not Implemented"))

    conn = request.getfixturevalue(conn)
    sql.to_sql(test_frame1, "test_frame", conn)

    result = sql.read_sql_table("test_frame", conn, index_col="index")
    assert result.index.names == ["index"]

    result = sql.read_sql_table("test_frame", conn, index_col=["A", "B"])
    assert result.index.names == ["A", "B"]

    result = sql.read_sql_table(
        "test_frame", conn, index_col=["A", "B"], columns=["C", "D"]
    )
    assert result.index.names == ["A", "B"]
    assert result.columns.tolist() == ["C", "D"]


@pytest.mark.parametrize("conn", all_connectable_iris)
def test_read_sql_delegate(conn, request):
    if conn == "sqlite_buildin_iris":
        request.applymarker(
            pytest.mark.xfail(
                reason="sqlite_buildin connection does not implement read_sql_table"
            )
        )

    conn = request.getfixturevalue(conn)
    iris_frame1 = sql.read_sql_query("SELECT * FROM iris", conn)
    iris_frame2 = sql.read_sql("SELECT * FROM iris", conn)
    tm.assert_frame_equal(iris_frame1, iris_frame2)

    iris_frame1 = sql.read_sql_table("iris", conn)
    iris_frame2 = sql.read_sql("iris", conn)
    tm.assert_frame_equal(iris_frame1, iris_frame2)


def test_not_reflect_all_tables(sqlite_conn):
    conn = sqlite_conn
    from sqlalchemy import text
    from sqlalchemy.engine import Engine

    # create invalid table
    query_list = [
        text("CREATE TABLE invalid (x INTEGER, y UNKNOWN);"),
        text("CREATE TABLE other_table (x INTEGER, y INTEGER);"),
    ]

    for query in query_list:
        if isinstance(conn, Engine):
            with conn.connect() as conn:
                with conn.begin():
                    conn.execute(query)
        else:
            with conn.begin():
                conn.execute(query)

    with tm.assert_produces_warning(None):
        sql.read_sql_table("other_table", conn)
        sql.read_sql_query("SELECT * FROM other_table", conn)


@pytest.mark.parametrize("conn", all_connectable)
def test_warning_case_insensitive_table_name(conn, request, test_frame1):
    conn_name = conn
    if conn_name == "sqlite_buildin" or "adbc" in conn_name:
        request.applymarker(pytest.mark.xfail(reason="Does not raise warning"))

    conn = request.getfixturevalue(conn)
    # see gh-7815
    with tm.assert_produces_warning(
        UserWarning,
        match=(
            r"The provided table name 'TABLE1' is not found exactly as such in "
            r"the database after writing the table, possibly due to case "
            r"sensitivity issues. Consider using lower case table names."
        ),
    ):
        with sql.SQLDatabase(conn) as db:
            db.check_case_sensitive("TABLE1", "")

    # Test that the warning is certainly NOT triggered in a normal case.
    with tm.assert_produces_warning(None):
        test_frame1.to_sql(name="CaseSensitive", con=conn)


@pytest.mark.parametrize("conn", sqlalchemy_connectable)
def test_sqlalchemy_type_mapping(conn, request):
    conn = request.getfixturevalue(conn)
    from sqlalchemy import TIMESTAMP

    # Test Timestamp objects (no datetime64 because of timezone) (GH9085)
    df = DataFrame(
        {"time": to_datetime(["2014-12-12 01:54", "2014-12-11 02:54"], utc=True)}
    )
    with sql.SQLDatabase(conn) as db:
        table = sql.SQLTable("test_type", db, frame=df)
        # GH 9086: TIMESTAMP is the suggested type for datetimes with timezones
        assert isinstance(table.table.c["time"].type, TIMESTAMP)


@pytest.mark.parametrize("conn", sqlalchemy_connectable)
@pytest.mark.parametrize(
    "integer, expected",
    [
        ("int8", "SMALLINT"),
        ("Int8", "SMALLINT"),
        ("uint8", "SMALLINT"),
        ("UInt8", "SMALLINT"),
        ("int16", "SMALLINT"),
        ("Int16", "SMALLINT"),
        ("uint16", "INTEGER"),
        ("UInt16", "INTEGER"),
        ("int32", "INTEGER"),
        ("Int32", "INTEGER"),
        ("uint32", "BIGINT"),
        ("UInt32", "BIGINT"),
        ("int64", "BIGINT"),
        ("Int64", "BIGINT"),
        (int, "BIGINT" if np.dtype(int).name == "int64" else "INTEGER"),
    ],
)
def test_sqlalchemy_integer_mapping(conn, request, integer, expected):
    # GH35076 Map pandas integer to optimal SQLAlchemy integer type
    conn = request.getfixturevalue(conn)
    df = DataFrame([0, 1], columns=["a"], dtype=integer)
    with sql.SQLDatabase(conn) as db:
        table = sql.SQLTable("test_type", db, frame=df)

        result = str(table.table.c.a.type)
    assert result == expected


@pytest.mark.parametrize("conn", sqlalchemy_connectable)
@pytest.mark.parametrize("integer", ["uint64", "UInt64"])
def test_sqlalchemy_integer_overload_mapping(conn, request, integer):
    conn = request.getfixturevalue(conn)
    # GH35076 Map pandas integer to optimal SQLAlchemy integer type
    df = DataFrame([0, 1], columns=["a"], dtype=integer)
    with sql.SQLDatabase(conn) as db:
        with pytest.raises(
            ValueError, match="Unsigned 64 bit integer datatype is not supported"
        ):
            sql.SQLTable("test_type", db, frame=df)


@pytest.mark.parametrize("conn", all_connectable)
def test_database_uri_string(conn, request, test_frame1):
    pytest.importorskip("sqlalchemy")
    conn = request.getfixturevalue(conn)
    # Test read_sql and .to_sql method with a database URI (GH10654)
    # db_uri = 'sqlite:///:memory:' # raises
    # sqlalchemy.exc.OperationalError: (sqlite3.OperationalError) near
    # "iris": syntax error [SQL: 'iris']
    with tm.ensure_clean() as name:
        db_uri = "sqlite:///" + name
        table = "iris"
        test_frame1.to_sql(name=table, con=db_uri, if_exists="replace", index=False)
        test_frame2 = sql.read_sql(table, db_uri)
        test_frame3 = sql.read_sql_table(table, db_uri)
        query = "SELECT * FROM iris"
        test_frame4 = sql.read_sql_query(query, db_uri)
    tm.assert_frame_equal(test_frame1, test_frame2)
    tm.assert_frame_equal(test_frame1, test_frame3)
    tm.assert_frame_equal(test_frame1, test_frame4)


@td.skip_if_installed("pg8000")
@pytest.mark.parametrize("conn", all_connectable)
def test_pg8000_sqlalchemy_passthrough_error(conn, request):
    pytest.importorskip("sqlalchemy")
    conn = request.getfixturevalue(conn)
    # using driver that will not be installed on CI to trigger error
    # in sqlalchemy.create_engine -> test passing of this error to user
    db_uri = "postgresql+pg8000://user:pass@host/dbname"
    with pytest.raises(ImportError, match="pg8000"):
        sql.read_sql("select * from table", db_uri)


@pytest.mark.parametrize("conn", sqlalchemy_connectable_iris)
def test_query_by_text_obj(conn, request):
    # WIP : GH10846
    conn_name = conn
    conn = request.getfixturevalue(conn)
    from sqlalchemy import text

    if "postgres" in conn_name:
        name_text = text('select * from iris where "Name"=:name')
    else:
        name_text = text("select * from iris where name=:name")
    iris_df = sql.read_sql(name_text, conn, params={"name": "Iris-versicolor"})
    all_names = set(iris_df["Name"])
    assert all_names == {"Iris-versicolor"}


@pytest.mark.parametrize("conn", sqlalchemy_connectable_iris)
def test_query_by_select_obj(conn, request):
    conn = request.getfixturevalue(conn)
    # WIP : GH10846
    from sqlalchemy import (
        bindparam,
        select,
    )

    iris = iris_table_metadata()
    name_select = select(iris).where(iris.c.Name == bindparam("name"))
    iris_df = sql.read_sql(name_select, conn, params={"name": "Iris-setosa"})
    all_names = set(iris_df["Name"])
    assert all_names == {"Iris-setosa"}


@pytest.mark.parametrize("conn", all_connectable)
def test_column_with_percentage(conn, request):
    # GH 37157
    conn_name = conn
    if conn_name == "sqlite_buildin":
        request.applymarker(pytest.mark.xfail(reason="Not Implemented"))

    conn = request.getfixturevalue(conn)
    df = DataFrame({"A": [0, 1, 2], "%_variation": [3, 4, 5]})
    df.to_sql(name="test_column_percentage", con=conn, index=False)

    res = sql.read_sql_table("test_column_percentage", conn)

    tm.assert_frame_equal(res, df)


def test_sql_open_close(test_frame3):
    # Test if the IO in the database still work if the connection closed
    # between the writing and reading (as in many real situations).

    with tm.ensure_clean() as name:
        with closing(sqlite3.connect(name)) as conn:
            assert sql.to_sql(test_frame3, "test_frame3_legacy", conn, index=False) == 4

        with closing(sqlite3.connect(name)) as conn:
            result = sql.read_sql_query("SELECT * FROM test_frame3_legacy;", conn)

    tm.assert_frame_equal(test_frame3, result)


@td.skip_if_installed("sqlalchemy")
def test_con_string_import_error():
    conn = "mysql://root@localhost/pandas"
    msg = "Using URI string without sqlalchemy installed"
    with pytest.raises(ImportError, match=msg):
        sql.read_sql("SELECT * FROM iris", conn)


@td.skip_if_installed("sqlalchemy")
def test_con_unknown_dbapi2_class_does_not_error_without_sql_alchemy_installed():
    class MockSqliteConnection:
        def __init__(self, *args, **kwargs) -> None:
            self.conn = sqlite3.Connection(*args, **kwargs)

        def __getattr__(self, name):
            return getattr(self.conn, name)

        def close(self):
            self.conn.close()

    with contextlib.closing(MockSqliteConnection(":memory:")) as conn:
        with tm.assert_produces_warning(UserWarning, match="only supports SQLAlchemy"):
            sql.read_sql("SELECT 1", conn)


def test_sqlite_read_sql_delegate(sqlite_buildin_iris):
    conn = sqlite_buildin_iris
    iris_frame1 = sql.read_sql_query("SELECT * FROM iris", conn)
    iris_frame2 = sql.read_sql("SELECT * FROM iris", conn)
    tm.assert_frame_equal(iris_frame1, iris_frame2)

    msg = "Execution failed on sql 'iris': near \"iris\": syntax error"
    with pytest.raises(sql.DatabaseError, match=msg):
        sql.read_sql("iris", conn)


def test_get_schema2(test_frame1):
    # without providing a connection object (available for backwards comp)
    create_sql = sql.get_schema(test_frame1, "test")
    assert "CREATE" in create_sql


def test_sqlite_type_mapping(sqlite_buildin):
    # Test Timestamp objects (no datetime64 because of timezone) (GH9085)
    conn = sqlite_buildin
    df = DataFrame(
        {"time": to_datetime(["2014-12-12 01:54", "2014-12-11 02:54"], utc=True)}
    )
    db = sql.SQLiteDatabase(conn)
    table = sql.SQLiteTable("test_type", db, frame=df)
    schema = table.sql_schema()
    for col in schema.split("\n"):
        if col.split()[0].strip('"') == "time":
            assert col.split()[1] == "TIMESTAMP"


# -----------------------------------------------------------------------------
# -- Database flavor specific tests


@pytest.mark.parametrize("conn", sqlalchemy_connectable)
def test_create_table(conn, request):
    if conn == "sqlite_str":
        pytest.skip("sqlite_str has no inspection system")

    conn = request.getfixturevalue(conn)

    from sqlalchemy import inspect

    temp_frame = DataFrame({"one": [1.0, 2.0, 3.0, 4.0], "two": [4.0, 3.0, 2.0, 1.0]})
    with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
        assert pandasSQL.to_sql(temp_frame, "temp_frame") == 4

    insp = inspect(conn)
    assert insp.has_table("temp_frame")

    # Cleanup
    with sql.SQLDatabase(conn, need_transaction=True) as pandasSQL:
        pandasSQL.drop_table("temp_frame")


@pytest.mark.parametrize("conn", sqlalchemy_connectable)
def test_drop_table(conn, request):
    if conn == "sqlite_str":
        pytest.skip("sqlite_str has no inspection system")

    conn = request.getfixturevalue(conn)

    from sqlalchemy import inspect

    temp_frame = DataFrame({"one": [1.0, 2.0, 3.0, 4.0], "two": [4.0, 3.0, 2.0, 1.0]})
    with sql.SQLDatabase(conn) as pandasSQL:
        with pandasSQL.run_transaction():
            assert pandasSQL.to_sql(temp_frame, "temp_frame") == 4

        insp = inspect(conn)
        assert insp.has_table("temp_frame")

        with pandasSQL.run_transaction():
            pandasSQL.drop_table("temp_frame")
        try:
            insp.clear_cache()  # needed with SQLAlchemy 2.0, unavailable prior
        except AttributeError:
            pass
        assert not insp.has_table("temp_frame")


@pytest.mark.parametrize("conn_name", all_connectable)
def test_delete_rows_success(conn_name, test_frame1, request):
    table_name = "temp_delete_rows_frame"
    conn = request.getfixturevalue(conn_name)

    with pandasSQL_builder(conn) as pandasSQL:
        with pandasSQL.run_transaction():
            assert pandasSQL.to_sql(test_frame1, table_name) == test_frame1.shape[0]

        with pandasSQL.run_transaction():
            assert pandasSQL.delete_rows(table_name) is None

        assert count_rows(conn, table_name) == 0
        assert pandasSQL.has_table(table_name)


@pytest.mark.parametrize("conn_name", all_connectable)
def test_delete_rows_is_atomic(conn_name, request):
    sqlalchemy = pytest.importorskip("sqlalchemy")

    table_name = "temp_delete_rows_atomic_frame"
    table_stmt = f"CREATE TABLE {table_name} (a INTEGER, b INTEGER UNIQUE NOT NULL)"

    if conn_name != "sqlite_buildin" and "adbc" not in conn_name:
        table_stmt = sqlalchemy.text(table_stmt)

    # setting dtype is mandatory for adbc related tests
    original_df = DataFrame({"a": [1, 2], "b": [3, 4]}, dtype="int32")
    replacing_df = DataFrame({"a": [5, 6, 7], "b": [8, 8, 8]}, dtype="int32")

    conn = request.getfixturevalue(conn_name)
    pandasSQL = pandasSQL_builder(conn)

    with pandasSQL.run_transaction() as cur:
        cur.execute(table_stmt)

    with pandasSQL.run_transaction():
        pandasSQL.to_sql(original_df, table_name, if_exists="append", index=False)

    # inserting duplicated values in a UNIQUE constraint column
    with pytest.raises(pd.errors.DatabaseError):
        with pandasSQL.run_transaction():
            pandasSQL.to_sql(
                replacing_df, table_name, if_exists="delete_rows", index=False
            )

    # failed "delete_rows" is rolled back preserving original data
    with pandasSQL.run_transaction():
        result_df = pandasSQL.read_query(f"SELECT * FROM {table_name}", dtype="int32")
        tm.assert_frame_equal(result_df, original_df)


@pytest.mark.parametrize("conn", all_connectable)
def test_roundtrip(conn, request, test_frame1):
    if conn == "sqlite_str":
        pytest.skip("sqlite_str has no inspection system")

    conn_name = conn
    conn = request.getfixturevalue(conn)
    pandasSQL = pandasSQL_builder(conn)
    with pandasSQL.run_transaction():
        assert pandasSQL.to_sql(test_frame1, "test_frame_roundtrip") == 4
        result = pandasSQL.read_query("SELECT * FROM test_frame_roundtrip")

    if "adbc" in conn_name:
        result = result.rename(columns={"__index_level_0__": "level_0"})
    result.set_index("level_0", inplace=True)
    # result.index.astype(int)

    result.index.name = None

    tm.assert_frame_equal(result, test_frame1)
