#!/usr/bin/env python
# coding: utf-8

from pytest import raises
from pandas.msgpack import packb, unpackb

def _decode_complex(obj):
    if b'__complex__' in obj:
        return complex(obj[b'real'], obj[b'imag'])
    return obj

def _encode_complex(obj):
    if isinstance(obj, complex):
        return {b'__complex__': True, b'real': 1, b'imag': 2}
    return obj

def test_encode_hook():
    packed = packb([3, 1+2j], default=_encode_complex)
    unpacked = unpackb(packed, use_list=1)
    assert unpacked[1] == {b'__complex__': True, b'real': 1, b'imag': 2}

def test_decode_hook():
    packed = packb([3, {b'__complex__': True, b'real': 1, b'imag': 2}])
    unpacked = unpackb(packed, object_hook=_decode_complex, use_list=1)
    assert unpacked[1] == 1+2j

def test_decode_pairs_hook():
    packed = packb([3, {1: 2, 3: 4}])
    prod_sum = 1 * 2 + 3 * 4
    unpacked = unpackb(packed, object_pairs_hook=lambda l: sum(k * v for k, v in l), use_list=1)
    assert unpacked[1] == prod_sum

def test_only_one_obj_hook():
    with raises(ValueError):
        unpackb(b'', object_hook=lambda x: x, object_pairs_hook=lambda x: x)

def test_bad_hook():
    with raises(ValueError):
        packed = packb([3, 1+2j], default=lambda o: o)
        unpacked = unpackb(packed, use_list=1)

def _arr_to_str(arr):
    return ''.join(str(c) for c in arr)

def test_array_hook():
    packed = packb([1,2,3])
    unpacked = unpackb(packed, list_hook=_arr_to_str, use_list=1)
    assert unpacked == '123'


class DecodeError(Exception):
    pass

def bad_complex_decoder(o):
    raise DecodeError("Ooops!")


def test_an_exception_in_objecthook1():
    with raises(DecodeError):
        packed = packb({1: {'__complex__': True, 'real': 1, 'imag': 2}})
        unpackb(packed, object_hook=bad_complex_decoder)


def test_an_exception_in_objecthook2():
    with raises(DecodeError):
        packed = packb({1: [{'__complex__': True, 'real': 1, 'imag': 2}]})
        unpackb(packed, list_hook=bad_complex_decoder, use_list=1)
