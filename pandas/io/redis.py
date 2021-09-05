import logging
import pickle
import uuid

from pandas.compat._optional import import_optional_dependency

import pandas as pd


def _try_import():
    # since pandas is a dependency of redis
    # we need to import on first use
    msg = (
        "The redis package is required to load data from Redis. "
        "See the docs: https://pypi.org/project/redis/"
    )
    redis = import_optional_dependency("redis", extra=msg)
    return redis


redis = _try_import()


class Redis_IO:
    """Core module for caching Pandas DataFrames to Redis."""

    def __init__(self, pandas_obj: pd.DataFrame):
        self._validate(pandas_obj)
        self.obj = pandas_obj
        self.alias = None
        self.logger = logging.getLogger("Redis_IO")
        self.to_cache = pandas_obj

    @staticmethod
    def _validate(obj: pd.DataFrame):
        """Simple method to validate that the method is being attached to a DataFame and not any other type.

        :param obj: object being chained via the namespace.
        :type obj: pd.DataFrame
        :raises TypeError: Not a DataFrame
        """
        if not isinstance(obj, pd.DataFrame):
            raise TypeError(
                f"Object to cache should be of type DataFrame, but is {type(obj)}"
            )

    def _validate_conn(self, redis_conn: redis.client.Redis):
        """Helper method to verify Redis connectivity

        :param redis_conn:
        :type redis_conn: redis.client.Redis
        :raises TypeError: Not a Redis connection object
        :raises ConnectionError: Can't connect to Redis, wrapper to the specific connection error from the Redis library
        """
        if not isinstance(redis_conn, redis.client.Redis):
            raise TypeError(
                f"Redis connection object should be of type redis.Client.Redis but is {type(redis_conn)}"
            )
        try:
            redis_conn.ping()
            self.logger.info("Connected to Redis")
        except redis.exceptions.ConnectionError as r_con_error:
            raise ConnectionError(
                f"Cannot connect to Redis instance: f{str(r_con_error)}"
            )

    def to_redis(
        self,
        redis_conn: redis.client.Redis,
        alias: str = None,
        if_exists: str = "Overwrite",
    ):
        """Pandas IO method to cache a Dataframe to Redis. This uses pickle (previously PyArrow, but the serialization was deprecated) to serialize the Dataframe and stores it as a string value in Redis.

        *Note* Redis has a 512 MB limit on values, so If the serialized Dataframe is over 512MB, a ValueError will be raised.

        If you are using if_exists="Append", the following logic is applied:
            - If the existing DataFrame has columns that don't exist in the new DataFrame, we just append which will add values to the columns that exist in the new DataFrame.
            - If the new DataFrame has columns that don't exist in the existing DataFrame, we throw an error. This is due to a desire to maintain integrity of what has already been cached.

        This operation is registered with Pandas under the redis namespace once installed. See example usage below:

        Parameters
        ----------
        :param redis_conn: Redis connection object created by the python-redis library
        :type redis_conn: redis.client.Redis
        :param alias: String used to reference the cached dataframe in redis, defaults to a randomly generated UUID. This can be accessed after caching via df.redis.alias
        :type alias: str, optional
        :param if_exists: How to handle an alias that already exists, options are "Overwrite", "Append", "Duplicate", "Quit". Note that duplicate will generate a UUID alias for the new cached DataFrame, defaults to "Overwrite"
        :type if_exists: str, optional

        .. highlight::

            import pandas as pd
            import redis

            df = pd.DataFrame(data=[1.2.3], columns=['A'. 'B', 'C'])

            redis_conn = redis.StrictRedis(host="your host", port=6379, db=0)

            df.redis.to_redis(redis_conn, alias="test")

            alias = df.redis.alias

        .. highlight::
        """

        if if_exists.upper() not in ["OVERWRITE", "APPEND", "DUPLICATE", "QUIT"]:
            raise ValueError(
                f"Keyword argument if_exists must be one of 'Overwrite', 'Append', 'Quit' or  'Duplicate', but recieved {if_exists}"
            )

        alias = alias if alias != None else str(uuid.uuid4())

        self._validate_conn(redis_conn)

        if redis_conn.exists(alias):
            logging.info(f"{alias} exists")
            if if_exists.upper() == "APPEND":
                existing_df = pd.read_redis(redis_conn, alias)
                extra_cols = [
                    col for col in self.obj.columns if col not in existing_df.columns
                ]
                if extra_cols:
                    raise IndexError(
                        f"Cannot append new DataFrame to one currently in cache. The new Dataframe contains the following extra columns{str(extra_cols)}"
                    )
                appended_df = existing_df.append(self.obj, ignore_index=True)
                self.to_cache = appended_df
            elif if_exists.upper() == "DUPLICATE":
                alias = str(uuid.uuid4())
            elif if_exists.upper() == "QUIT":
                return

        df_serialzed = pickle.dumps(self.obj)

        res = redis_conn.set(alias, df_serialzed)
        if res == True:
            self.alias = alias
            logging.info("Dataframe cached")
            logging.info(f"Key {alias}")


def read_redis(redis_conn: redis.client.Redis, key: str) -> pd.DataFrame:
    """Pandas IO method to retrieve a cached Dataframe.

    Parameters
    ----------
    :param redis_conn: Redis connection object
    :type redis_conn: redis.client.Redis
    :param key: String used to refer to the cached DataFrame you wish to retrieve.
    :type key: str
    :raises KeyError: These are not the droids you're looking for.
    :return: DataFrame retrieved from cache.
    :rtype: pd.DataFrame

    .. highlight::
            import pandas as pd
            import redis

            redis_conn = redis.StrictRedis(host="your host", port=6379, db=0)

            df = pd.read_redis(redis_conn, alias="test")

    .. highlight::

    """

    redis_pandas = Redis_IO(pd.DataFrame())
    redis_pandas._validate_conn(redis_conn)

    if not redis_conn.exists(key):
        raise KeyError(f"{key} does not exist in the specifed Redis instance")

    serialized_df = redis_conn.get(key)
    dataframe = pickle.loads(serialized_df)
    return dataframe
