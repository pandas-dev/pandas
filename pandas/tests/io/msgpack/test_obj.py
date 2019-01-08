# coding: utf-8

import pytest

from pandas.io.msgpack import packb, unpackb


class DecodeError(Exception):
    pass


class TestObj(object):

    def _arr_to_str(self, arr):
        return ''.join(str(c) for c in arr)

    def bad_complex_decoder(self, o):
        raise DecodeError("Ooops!")

    def _decode_complex(self, obj):
        if b'__complex__' in obj:
            return complex(obj[b'real'], obj[b'imag'])
        return obj

    def _encode_complex(self, obj):
        if isinstance(obj, complex):
            return {b'__complex__': True, b'real': 1, b'imag': 2}
        return obj

    def test_encode_hook(self):
        packed = packb([3, 1 + 2j], default=self._encode_complex)
        unpacked = unpackb(packed, use_list=1)
        assert unpacked[1] == {b'__complex__': True, b'real': 1, b'imag': 2}

    def test_decode_hook(self):
        packed = packb([3, {b'__complex__': True, b'real': 1, b'imag': 2}])
        unpacked = unpackb(packed, object_hook=self._decode_complex,
                           use_list=1)
        assert unpacked[1] == 1 + 2j

    def test_decode_pairs_hook(self):
        packed = packb([3, {1: 2, 3: 4}])
        prod_sum = 1 * 2 + 3 * 4
        unpacked = unpackb(
            packed, object_pairs_hook=lambda l: sum(k * v for k, v in l),
            use_list=1)
        assert unpacked[1] == prod_sum

    def test_only_one_obj_hook(self):
        msg = "object_pairs_hook and object_hook are mutually exclusive"
        with pytest.raises(TypeError, match=msg):
            unpackb(b'', object_hook=lambda x: x,
                    object_pairs_hook=lambda x: x)

    def test_bad_hook(self):
        msg = r"can't serialize \(1\+2j\)"
        with pytest.raises(TypeError, match=msg):
            packed = packb([3, 1 + 2j], default=lambda o: o)
            unpacked = unpackb(packed, use_list=1)  # noqa

    def test_array_hook(self):
        packed = packb([1, 2, 3])
        unpacked = unpackb(packed, list_hook=self._arr_to_str, use_list=1)
        assert unpacked == '123'

    def test_an_exception_in_objecthook1(self):
        with pytest.raises(DecodeError, match='Ooops!'):
            packed = packb({1: {'__complex__': True, 'real': 1, 'imag': 2}})
            unpackb(packed, object_hook=self.bad_complex_decoder)

    def test_an_exception_in_objecthook2(self):
        with pytest.raises(DecodeError, match='Ooops!'):
            packed = packb({1: [{'__complex__': True, 'real': 1, 'imag': 2}]})
            unpackb(packed, list_hook=self.bad_complex_decoder, use_list=1)
