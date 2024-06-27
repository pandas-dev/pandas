from functools import partial
from .parquet.ttypes import *


def __getattr__(name):
    # for compatability with coe that calls, e.g., parquet_thrift.RowGroup(...)
    from fastparquet.cencoding import ThriftObject
    if name[0].isupper():
        return partial(ThriftObject.from_fields, thrift_name=name)
    raise AttributeError(name)
