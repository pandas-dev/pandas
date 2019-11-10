### REPRODUCIBLE SQLTable Creation:table
import sqlalchemy
from sqlalchemy import Table, create_engine, select
from sqlalchemy.engine.base import Connection
from sqlalchemy.sql import tuple_

import pandas as pd
from vortexa_utils.database import ProdFactory

from pandas.io.sql import SQLDatabase, SQLTable


def get_pkey(table: Table):
    return [pkey.name for pkey in table.primary_key.columns.values()]


def get_pkey_values(table: Table, conn: Connection):
    pkeys = get_pkey(table)
    statement = select([table.c[name] for name in pkeys])
    return [row for row in conn.execute(statement)]
    # for row in conn.execute(statement):
    #   yield row


def pkey_generator(table, engine):
    pkeys = get_pkey(table)
    statement = select([table.c[name] for name in pkeys])
    with engine.connect() as conn:
        for row in conn.execute(statement):
            yield row


# Leaves connection open
def pkey_results_proxy(table, engine):
    pkeys = get_pkey(table)
    statement = select([table.c[name] for name in pkeys])
    with engine.connect() as conn:
        result = conn.execute(statement)
    return result


def pkey_generator2(table, engine):
    pkeys = get_pkey(table)
    statement = select([table.c[name] for name in pkeys])
    with engine.connect() as conn:
        result = conn.execute(statement)
    try:
        for row in result:
            yield result.fetchone()
    finally:
        result.close()


# replace table with self
def get_pkey_values(table: SQLTable):
    pkeys = [pkey.name for pkey in table.table.primary_key.columns.values()]
    statement = select([table.table.c[name] for name in pkeys])
    table.pd_sql.execute(statement)


def generate_mask(df, dictionary):
    return [df[key] == value for key, value in dictionary.items()]


def generate_mask_of_masks(list_of_masks):
    return pd.concat([mask for mask in list_of_masks], axis=1).all(1)


engine = sqlalchemy.create_engine("enter string here")
meta = MetaData(engine)
table_name = "charterers"  # or wtv
meta.reflect(only=[table_name], views=True)
db = SQLDatabase(engine, meta=meta)
table = SQLTable(table_name, db, index=None, schema=None)


engine_v = ProdFactory().engine()
engine = create_engine("sqlite:///:memory:")
table_name = "charterers"
df = pd.read_sql_table(table_name, engine_v)
df_test = df.head().copy()
df_test["name"] = df_test["name"].apply(lambda x: x + "_TEST")
engine.execute(
    "create table charterers(id text primary key, name text, energy integer)"
)
def create_test_df(df):
    df2 = df.head().copy()
    df2['name'] = df2['name'].apply(lambda x: x + '_NEW')
    return df2

def read_table(table):
    with engine.connect() as conn:
        result = conn.execute(f'select * from {table}')
        return result.fetchall()

def clear_table(table):
    with engine.connect() as conn:
        conn.execute(f'delete from {table}')

def top_up_table(table):
    df.to_sql(table, con=engine, if_exists='append', index=False)
    return read_table()

def reset_table(table):
    clear_table(table)
    top_up_table(table)

df.to_sql(table_name, index=False, if_exists="append", con=engine)

db = SQLDatabase(engine, schema=None, meta=None)
new_data = SQLTable(table_name, db, frame=df_test, index=False)


def delete_matching_keys(sql_table, key_columns, value_iter):
    delete_expression = sql_table.table.delete().where(
        tuple_(*(table.table.c[col] for col in key_columns)).in_(list(zip(value_iter)))
    )
    with sql_table.pd_sql.run_transaction() as conn:
        conn.execute(delete_expression)
