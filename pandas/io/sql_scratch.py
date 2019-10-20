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
