# coding: utf-8
from __future__ import (
    absolute_import, division, print_function, unicode_literals)

import pytest

from pandas.io.msgpack import ExtType, Packer, Unpacker, packb, unpackb


class TestLimits(object):

    def test_integer(self):
        x = -(2 ** 63)
        assert unpackb(packb(x)) == x
        msg = (r"((long |Python )?(int )?too (big|large) to convert"
               r"( to C (unsigned )?long))?")
        with pytest.raises((OverflowError, ValueError), match=msg):
            packb(x - 1)
        x = 2 ** 64 - 1
        assert unpackb(packb(x)) == x
        with pytest.raises((OverflowError, ValueError), match=msg):
            packb(x + 1)

    def test_array_header(self):
        packer = Packer()
        packer.pack_array_header(2 ** 32 - 1)
        with pytest.raises((OverflowError, ValueError)):
            packer.pack_array_header(2 ** 32)

    def test_map_header(self):
        packer = Packer()
        packer.pack_map_header(2 ** 32 - 1)
        with pytest.raises((OverflowError, ValueError)):
            packer.pack_array_header(2 ** 32)

    def test_max_str_len(self):
        d = 'x' * 3
        packed = packb(d)

        unpacker = Unpacker(max_str_len=3, encoding='utf-8')
        unpacker.feed(packed)
        assert unpacker.unpack() == d

        unpacker = Unpacker(max_str_len=2, encoding='utf-8')
        unpacker.feed(packed)

        msg = "3 exceeds max_str_len"
        with pytest.raises(ValueError, match=msg):
            unpacker.unpack()

    def test_max_bin_len(self):
        d = b'x' * 3
        packed = packb(d, use_bin_type=True)

        unpacker = Unpacker(max_bin_len=3)
        unpacker.feed(packed)
        assert unpacker.unpack() == d

        unpacker = Unpacker(max_bin_len=2)
        unpacker.feed(packed)

        msg = "3 exceeds max_bin_len"
        with pytest.raises(ValueError, match=msg):
            unpacker.unpack()

    def test_max_array_len(self):
        d = [1, 2, 3]
        packed = packb(d)

        unpacker = Unpacker(max_array_len=3)
        unpacker.feed(packed)
        assert unpacker.unpack() == d

        unpacker = Unpacker(max_array_len=2)
        unpacker.feed(packed)

        msg = "3 exceeds max_array_len"
        with pytest.raises(ValueError, match=msg):
            unpacker.unpack()

    def test_max_map_len(self):
        d = {1: 2, 3: 4, 5: 6}
        packed = packb(d)

        unpacker = Unpacker(max_map_len=3)
        unpacker.feed(packed)
        assert unpacker.unpack() == d

        unpacker = Unpacker(max_map_len=2)
        unpacker.feed(packed)

        msg = "3 exceeds max_map_len"
        with pytest.raises(ValueError, match=msg):
            unpacker.unpack()

    def test_max_ext_len(self):
        d = ExtType(42, b"abc")
        packed = packb(d)

        unpacker = Unpacker(max_ext_len=3)
        unpacker.feed(packed)
        assert unpacker.unpack() == d

        unpacker = Unpacker(max_ext_len=2)
        unpacker.feed(packed)

        msg = "4 exceeds max_ext_len"
        with pytest.raises(ValueError, match=msg):
            unpacker.unpack()
