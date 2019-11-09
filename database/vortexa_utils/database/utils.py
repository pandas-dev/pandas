from typing import Iterable, List

import sqlalchemy
from pandas.io.sql import SQLTable
from sqlalchemy.engine import Connectable


def upsert(
    table: SQLTable, conn: Connectable, keys: List[str], data_iter: Iterable
):
    """Upsert method to be used with `pandas.DataFrame.to_sql`.

    In pandas > 0.24.0 you can specify a method to control the insertion clause
    used by `pandas.DataFrame.to_sql`.

    Parameters
    ----------
    table : pandas.io.sql.SQLTable
        Description of parameter `table`.
    conn : sqlalchemy.engine.Connectable
        Description of parameter `conn`.
    keys : List[str]
        Description of parameter `keys`.
    data_iter : Iterable
        Description of parameter `data_iter`.

    Returns
    -------
    type
        Description of returned object.

    Examples
    -------
    Examples should be written in doctest format, and
    should illustrate how to use the function/class.
    >>>

    """
    cols = ", ".join(f'"{k}"' for k in keys)
    if table.schema:
        tname = "{}.{}".format(table.schema, table.name)
    else:
        tname = table.name

    # placeholder = ", ".join(["?"] * len(keys))
    placeholder = ", ".join([f":{k}" for k in keys])
    datas = ({k: d for k, d in zip(keys, data)} for data in data_iter)
    if conn.engine.driver.endswith("sqlite"):
        # sqlite
        sql = f"INSERT or IGNORE INTO {tname} ({cols}) VALUES ({placeholder})"
    else:
        # postgresql
        sql = sqlalchemy.text(
            f"""
                INSERT INTO {tname}
                ({cols})
                VALUES ({placeholder})
                ON CONFLICT DO NOTHING
        """
        )

    conn.execute(sql, *datas)
