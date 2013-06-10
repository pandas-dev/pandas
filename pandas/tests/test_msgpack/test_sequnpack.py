#!/usr/bin/env python
# coding: utf-8

import six
from pandas.msgpack import Unpacker, BufferFull
from pandas.msgpack import OutOfData
from pytest import raises


def test_partialdata():
    unpacker = Unpacker()
    unpacker.feed(b'\xa5')
    with raises(StopIteration): next(iter(unpacker))
    unpacker.feed(b'h')
    with raises(StopIteration): next(iter(unpacker))
    unpacker.feed(b'a')
    with raises(StopIteration): next(iter(unpacker))
    unpacker.feed(b'l')
    with raises(StopIteration): next(iter(unpacker))
    unpacker.feed(b'l')
    with raises(StopIteration): next(iter(unpacker))
    unpacker.feed(b'o')
    assert next(iter(unpacker)) == b'hallo'

def test_foobar():
    unpacker = Unpacker(read_size=3, use_list=1)
    unpacker.feed(b'foobar')
    assert unpacker.unpack() == ord(b'f')
    assert unpacker.unpack() == ord(b'o')
    assert unpacker.unpack() == ord(b'o')
    assert unpacker.unpack() == ord(b'b')
    assert unpacker.unpack() == ord(b'a')
    assert unpacker.unpack() == ord(b'r')
    with raises(OutOfData):
        unpacker.unpack()

    unpacker.feed(b'foo')
    unpacker.feed(b'bar')

    k = 0
    for o, e in zip(unpacker, 'foobarbaz'):
        assert o == ord(e)
        k += 1
    assert k == len(b'foobar')

def test_foobar_skip():
    unpacker = Unpacker(read_size=3, use_list=1)
    unpacker.feed(b'foobar')
    assert unpacker.unpack() == ord(b'f')
    unpacker.skip()
    assert unpacker.unpack() == ord(b'o')
    unpacker.skip()
    assert unpacker.unpack() == ord(b'a')
    unpacker.skip()
    with raises(OutOfData):
        unpacker.unpack()

def test_maxbuffersize():
    with raises(ValueError):
        Unpacker(read_size=5, max_buffer_size=3)
    unpacker = Unpacker(read_size=3, max_buffer_size=3, use_list=1)
    unpacker.feed(b'fo')
    with raises(BufferFull):
        unpacker.feed(b'ob')
    unpacker.feed(b'o')
    assert ord('f') == next(unpacker)
    unpacker.feed(b'b')
    assert ord('o') == next(unpacker)
    assert ord('o') == next(unpacker)
    assert ord('b') == next(unpacker)


def test_readbytes():
    unpacker = Unpacker(read_size=3)
    unpacker.feed(b'foobar')
    assert unpacker.unpack() == ord(b'f')
    assert unpacker.read_bytes(3) == b'oob'
    assert unpacker.unpack() == ord(b'a')
    assert unpacker.unpack() == ord(b'r')

    # Test buffer refill
    unpacker = Unpacker(six.BytesIO(b'foobar'), read_size=3)
    assert unpacker.unpack() == ord(b'f')
    assert unpacker.read_bytes(3) == b'oob'
    assert unpacker.unpack() == ord(b'a')
    assert unpacker.unpack() == ord(b'r')
