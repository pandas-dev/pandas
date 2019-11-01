from sqlalchemy import Table, select
from sqlalchemy.engine.base import Connection


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
    for row in result:
        yield result.fetchone()

# replace table with self
def get_pkey_values(table: SQLTable):
    pkeys = [pkey.name for pkey in table.table.primary_key.columns.values()]
    statement = select([table.table.c[name] for name in pkeys])
    table.pd_sql.execute(statement)


### REPRODUCIBLE SQLTable Creation:
import sqlalchemy

engine = sqlalchemy.create_engine('enter string here')
meta = MetaData(engine)
table_name = 'charterers' # or wtv
meta.reflect(only=[table_name], views=True)
db = SQLDatabase(engine, meta=meta)
table = SQLTable(table_name, db, index=None, schema=None)