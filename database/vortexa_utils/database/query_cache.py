import hashlib
import logging
import os
import time
from base64 import urlsafe_b64encode
from dataclasses import dataclass
from functools import wraps
from typing import Union

import pandas as pd
from sqlalchemy.engine import Connection, Engine

from pyarrow.lib import ArrowIOError

logger = logging.getLogger(__name__)


@dataclass
class QueryCache(object):
    result_extension: str = ".parquet.snappy"
    cache_dir: str = os.path.join(
        "/tmp", "python_utils", "database", "query_cache", "df_cache"
    )
    ttl: int = 3600

    def __post_init__(self):
        os.makedirs(self.cache_dir, exist_ok=True)

    def path(self, url):
        return os.path.join(
            self.cache_dir,
            url.drivername,
            f"{url.host}:{url.port}",
            url.database,
        )

    def filename(self, query):
        query_digest = urlsafe_b64encode(
            hashlib.blake2s(str(query).encode(), digest_size=8).digest()
        )
        return query_digest.decode("ascii") + self.result_extension

    @wraps(pd.read_sql)
    def read_sql(
        self,
        query: str,
        con: Union[Engine, Connection],
        ttl: int = None,
        invalidate_cache: bool = False,
        *args,
        **kwargs,
    ) -> pd.DataFrame:

        # formulate a path
        path = self.path(con.engine.url)
        filename = self.filename(query)
        filepath = os.path.join(path, filename)
        os.makedirs(path, exist_ok=True)

        # check if the cache exists and is valid
        ttl = self.ttl if ttl is None else ttl

        if (
            os.path.isfile(filepath)
            and time.time() - os.path.getmtime(filepath) < ttl
        ):
            try:
                logger.debug("reading from cache %s", filepath)
                df = pd.read_parquet(filepath)
            except ArrowIOError as e:
                logger.error("Invalid Cache file, error: %s", e)
            else:
                return df
        logger.debug("reading from database")
        df = pd.read_sql(query, con, *args, **kwargs)
        df.to_parquet(filepath)
        return df
