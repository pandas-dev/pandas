from pandas.io.json._json import (
    read_json,
    to_json,
    ujson_dumps,
    ujson_loads,
)
from pandas.io.json._table_schema import build_table_schema

ujson_dumps.__module__ = "pandas.io.json"
ujson_loads.__module__ = "pandas.io.json"

__all__ = [
    "build_table_schema",
    "read_json",
    "to_json",
    "ujson_dumps",
    "ujson_loads",
]
