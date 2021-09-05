import random
import sys 

import pandas
import pytest

import redis
import fakeredis

def generate_df(n_rows: int, n_cols: int) -> pandas.DataFrame:
    """Method for generating a random Dataframe of n number of rows by n number of columns.

    :param n_rows: Number of rows to populate.
    :type n_rows: int
    :param n_cols: Number of columns to populate
    :type n_cols: int
    :return: Randomly generated DataFrame
    :rtype: pd.DataFrame
    """

    data = []
    cols = [str(x) for x in range(1, n_cols)]
    for i in range(1, n_rows):
        col_vals = [random.randint(1000000, 9999999) for y in range(1, n_cols)]
        data.append(col_vals)

    df = pandas.DataFrame(data, columns=cols)
    return df

@pytest.fixture()
def rp():
    df = generate_df(50, 8)
    rp = pandas.io.redis.Redis_IO(df)
    return rp


def test_object_validation_success():
    """Test that the dataframe validation works on a dataframe
    """
    df = generate_df(50, 8)
    assert pandas.io.redis.Redis_IO(df)


def test_object_validation_fails():
    """Test that the dataframe validation fails on objects other than dataframe
    """
    with pytest.raises(TypeError):
        pandas.io.redis.Redis_IO(5)


def test_redis_conn_validation_success(rp):
    """Test that the Redis connection validation works successfully
    """
    assert rp._validate_conn(fakeredis.FakeStrictRedis()) is None


def test_redis_conn_validation_connection_failure(rp):
    """Test that the Redis connection validation successfully errors when it can't connect to the provided Redis instance 
    """
    conn = redis.StrictRedis(host="fake_host", port=6379, db=0)

    with pytest.raises(ConnectionError):
        rp._validate_conn(conn)


def test_to_redis_invalid_exists_param(rp):
    """Test that the to_redis method enforces allowed values for the if_exists param
    """

    with pytest.raises(ValueError):
        rp.to_redis(fakeredis.FakeStrictRedis(), if_exists="Do Nothing")

def test_to_redis_if_exists_quit_success(rp):
    """Test that the to_redis method successfully returns when if_exists="Append"
    """
    conn = fakeredis.FakeStrictRedis()
    rp.to_redis(conn, alias='quit')

    df2 = generate_df(50, 8)
    rp2 = pandas.io.redis.Redis_IO(df2)
    assert rp2.to_redis(conn, alias='quit', if_exists="Quit") is None



def test_to_redis_if_exists_append_success(rp):
    """Test that the to_redis method successfully appends a given dataframe to the existing DataFrame when if_exists="Append"
    """
    conn = fakeredis.FakeStrictRedis()
    rp.to_redis(conn, alias='append')

    df2 = generate_df(50, 8)
    rp2 = pandas.io.redis.Redis_IO(df2)
    assert rp2.to_redis(conn, alias='append', if_exists="Append") is None
    assert rp2.alias == "append"


def test_to_redis_if_exists_append_less_columns_success(rp):
    """Test that the to_redis method successfully appends a given dataframe to the existing DataFrame when if_exists="Append" and the new DataFrame has a subset of columns of the cached DataFrame. 
    """
    conn = fakeredis.FakeStrictRedis()
    rp.to_redis(conn, alias='append')

    df2 = generate_df(50, 6) 
    rp2 = pandas.io.redis.Redis_IO(df2)
    assert rp2.to_redis(conn, alias='append', if_exists="Append") is None

def test_to_redis_if_exists_append_more_columns_error(rp):
    """Test that the to_redis method fails to append a given dataframe to the existing DataFrame when if_exists="Append" and the new DataFrame has columns that are not present in the cached DataFrame. 
    """
    conn = fakeredis.FakeStrictRedis()
    rp.to_redis(conn, alias='append')

    df2 = generate_df(50, 10) 
    rp2 = pandas.io.redis.Redis_IO(df2)
    with pytest.raises(IndexError):
        rp2.to_redis(conn, alias='append', if_exists="Append")

def test_to_redis_if_exists_duplicate_success(rp):
    """Test that the to_redis method successfully duplicates a given dataframe when if_exists="Duplicate". 
    """
    conn = fakeredis.FakeStrictRedis()
    rp.to_redis(conn, alias='duplicate')

    df2 = generate_df(50, 8)
    rp2 = pandas.io.redis.Redis_IO(df2)
    assert rp2.to_redis(conn, alias='duplicate', if_exists="Duplicate") is None

def test_to_redis_success(rp):
    """Test that the to_redis method successfully caches a given DataFrame. 
    """
    assert rp.to_redis(fakeredis.FakeStrictRedis(), alias='SuccessfulTest') is None
    assert rp.alias == 'SuccessfulTest'


def test_read_redis_doesnt_exist_error():
    """Test that the read_redis method throws an error when the requested key does not exist
    """
    with pytest.raises(KeyError):
        pandas.io.redis.read_redis(fakeredis.FakeStrictRedis(), "DoesntExist")
    

def test_read_redis_success(rp):
    """Test that the read_redis method successfully retrieves a cached DataFrame.
    """
    conn = fakeredis.FakeStrictRedis()
    rp.to_redis(conn, alias='readRedis')
    
    assert pandas.io.redis.read_redis(redis_conn=conn, key="readRedis") is not None

