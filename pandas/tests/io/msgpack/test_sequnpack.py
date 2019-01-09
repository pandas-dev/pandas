# coding: utf-8

import pytest

from pandas import compat

from pandas.io.msgpack import BufferFull, OutOfData, Unpacker


class TestPack(object):

    def test_partial_data(self):
        unpacker = Unpacker()
        msg = "No more data to unpack"

        for data in [b"\xa5", b"h", b"a", b"l", b"l"]:
            unpacker.feed(data)
            with pytest.raises(StopIteration, match=msg):
                next(iter(unpacker))

        unpacker.feed(b"o")
        assert next(iter(unpacker)) == b"hallo"

    def test_foobar(self):
        unpacker = Unpacker(read_size=3, use_list=1)
        unpacker.feed(b'foobar')
        assert unpacker.unpack() == ord(b'f')
        assert unpacker.unpack() == ord(b'o')
        assert unpacker.unpack() == ord(b'o')
        assert unpacker.unpack() == ord(b'b')
        assert unpacker.unpack() == ord(b'a')
        assert unpacker.unpack() == ord(b'r')
        msg = "No more data to unpack"
        with pytest.raises(OutOfData, match=msg):
            unpacker.unpack()

        unpacker.feed(b'foo')
        unpacker.feed(b'bar')

        k = 0
        for o, e in zip(unpacker, 'foobarbaz'):
            assert o == ord(e)
            k += 1
        assert k == len(b'foobar')

    def test_foobar_skip(self):
        unpacker = Unpacker(read_size=3, use_list=1)
        unpacker.feed(b'foobar')
        assert unpacker.unpack() == ord(b'f')
        unpacker.skip()
        assert unpacker.unpack() == ord(b'o')
        unpacker.skip()
        assert unpacker.unpack() == ord(b'a')
        unpacker.skip()
        msg = "No more data to unpack"
        with pytest.raises(OutOfData, match=msg):
            unpacker.unpack()

    def test_maxbuffersize_read_size_exceeds_max_buffer_size(self):
        msg = "read_size should be less or equal to max_buffer_size"
        with pytest.raises(ValueError, match=msg):
            Unpacker(read_size=5, max_buffer_size=3)

    def test_maxbuffersize_bufferfull(self):
        unpacker = Unpacker(read_size=3, max_buffer_size=3, use_list=1)
        unpacker.feed(b'foo')
        with pytest.raises(BufferFull, match=r'^$'):
            unpacker.feed(b'b')

    def test_maxbuffersize(self):
        unpacker = Unpacker(read_size=3, max_buffer_size=3, use_list=1)
        unpacker.feed(b'foo')
        assert ord('f') == next(unpacker)
        unpacker.feed(b'b')
        assert ord('o') == next(unpacker)
        assert ord('o') == next(unpacker)
        assert ord('b') == next(unpacker)

    def test_readbytes(self):
        unpacker = Unpacker(read_size=3)
        unpacker.feed(b'foobar')
        assert unpacker.unpack() == ord(b'f')
        assert unpacker.read_bytes(3) == b'oob'
        assert unpacker.unpack() == ord(b'a')
        assert unpacker.unpack() == ord(b'r')

        # Test buffer refill
        unpacker = Unpacker(compat.BytesIO(b'foobar'), read_size=3)
        assert unpacker.unpack() == ord(b'f')
        assert unpacker.read_bytes(3) == b'oob'
        assert unpacker.unpack() == ord(b'a')
        assert unpacker.unpack() == ord(b'r')

    def test_issue124(self):
        unpacker = Unpacker()
        unpacker.feed(b'\xa1?\xa1!')
        assert tuple(unpacker) == (b'?', b'!')
        assert tuple(unpacker) == ()
        unpacker.feed(b"\xa1?\xa1")
        assert tuple(unpacker) == (b'?', )
        assert tuple(unpacker) == ()
        unpacker.feed(b"!")
        assert tuple(unpacker) == (b'!', )
        assert tuple(unpacker) == ()
