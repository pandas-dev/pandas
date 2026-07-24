from __future__ import annotations

import os
import random

import pytest

import dask.bag as db

fastavro = pytest.importorskip("fastavro")

expected = [
    {
        "name": random.choice(["fred", "wilma", "barney", "betty"]),
        "number": random.randint(0, 100),
    }
    for _ in range(1000)
]
schema = {
    "doc": "Descr",
    "name": "Random",
    "namespace": "test",
    "type": "record",
    "fields": [{"name": "name", "type": "string"}, {"name": "number", "type": "int"}],
}


def test_onefile_oneblock(tmpdir):
    tmpdir = str(tmpdir)
    fn = os.path.join(tmpdir, "one.avro")
    with open(fn, "wb") as f:
        fastavro.writer(f, records=expected, schema=schema)
    b = db.read_avro(fn, blocksize=None)
    assert b.npartitions == 1
    assert b.compute() == expected


def test_twofile_oneblock(tmpdir):
    tmpdir = str(tmpdir)
    fn1 = os.path.join(tmpdir, "one.avro")
    fn2 = os.path.join(tmpdir, "two.avro")
    with open(fn1, "wb") as f:
        fastavro.writer(f, records=expected[:500], schema=schema)
    with open(fn2, "wb") as f:
        fastavro.writer(f, records=expected[500:], schema=schema)
    b = db.read_avro(os.path.join(tmpdir, "*.avro"), blocksize=None)
    assert b.npartitions == 2
    assert b.compute() == expected


def test_twofile_multiblock(tmpdir):
    tmpdir = str(tmpdir)
    fn1 = os.path.join(tmpdir, "one.avro")
    fn2 = os.path.join(tmpdir, "two.avro")
    with open(fn1, "wb") as f:
        fastavro.writer(f, records=expected[:500], schema=schema, sync_interval=100)
    with open(fn2, "wb") as f:
        fastavro.writer(f, records=expected[500:], schema=schema, sync_interval=100)
    b = db.read_avro(os.path.join(tmpdir, "*.avro"), blocksize=None)
    assert b.npartitions == 2
    assert b.compute() == expected

    b = db.read_avro(os.path.join(tmpdir, "*.avro"), blocksize=1000)
    assert b.npartitions > 2
    assert b.compute() == expected


def test_roundtrip_simple(tmpdir):
    from dask.delayed import Delayed

    tmpdir = str(tmpdir)
    fn = os.path.join(tmpdir, "out*.avro")
    b = db.from_sequence([{"a": i} for i in [1, 2, 3, 4, 5]], npartitions=2)
    schema = {
        "name": "Test",
        "type": "record",
        "fields": [{"name": "a", "type": "int"}],
    }
    out = b.to_avro(fn, schema, compute=False)
    assert isinstance(out[0], Delayed)
    out = b.to_avro(fn, schema)
    assert len(out) == 2
    b2 = db.read_avro(fn)
    assert b.compute() == b2.compute()


@pytest.mark.parametrize("codec", ["null", "deflate", "snappy"])
def test_roundtrip(tmpdir, codec):
    tmpdir = str(tmpdir)
    if codec == "snappy":
        pytest.importorskip("snappy")
    fn = os.path.join(tmpdir, "out*.avro")
    b = db.from_sequence(expected, npartitions=3)
    b.to_avro(fn, schema=schema, codec=codec)
    b2 = db.read_avro(fn)
    assert b.compute() == b2.compute()


def test_invalid_schema(tmpdir):
    tmpdir = str(tmpdir)
    b = db.from_sequence(expected, npartitions=3)
    fn = os.path.join(tmpdir, "out*.avro")
    with pytest.raises(AssertionError):
        b.to_avro(fn, schema=[])
    with pytest.raises(AssertionError):
        b.to_avro(fn, schema={})
    with pytest.raises(AssertionError):
        b.to_avro(fn, schema={"doc": "unknown"})
    with pytest.raises(AssertionError):
        b.to_avro(fn, schema={"name": "test"})
    with pytest.raises(AssertionError):
        b.to_avro(fn, schema={"name": "test", "type": "wrong"})
    with pytest.raises(AssertionError):
        b.to_avro(fn, schema={"name": "test", "type": "record"})
    with pytest.raises(AssertionError):
        b.to_avro(fn, schema={"name": "test", "type": "record"})
    with pytest.raises(AssertionError):
        b.to_avro(
            fn, schema={"name": "test", "type": "record", "fields": [{"name": "a"}]}
        )
