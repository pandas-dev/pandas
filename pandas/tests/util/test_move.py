# -*- coding: utf-8 -*-
import pytest

from pandas.util._move import BadMove, move_into_mutable_buffer, stolenbuf


def test_cannot_create_instance_of_stolen_buffer():
    # Stolen buffers need to be created through the smart constructor
    # "move_into_mutable_buffer," which has a bunch of checks in it.

    msg = "cannot create 'pandas.util._move.stolenbuf' instances"
    with pytest.raises(TypeError, match=msg):
        stolenbuf()


def test_more_than_one_ref():
    # Test case for when we try to use "move_into_mutable_buffer"
    # when the object being moved has other references.

    b = b"testing"

    with pytest.raises(BadMove, match="testing") as e:
        def handle_success(type_, value, tb):
            assert value.args[0] is b
            return type(e).handle_success(e, type_, value, tb)  # super

        e.handle_success = handle_success
        move_into_mutable_buffer(b)


def test_exactly_one_ref():
    # Test case for when the object being moved has exactly one reference.

    b = b"testing"

    # We need to pass an expression on the stack to ensure that there are
    # not extra references hanging around. We cannot rewrite this test as
    #   buf = b[:-3]
    #   as_stolen_buf = move_into_mutable_buffer(buf)
    # because then we would have more than one reference to buf.
    as_stolen_buf = move_into_mutable_buffer(b[:-3])

    # Materialize as byte-array to show that it is mutable.
    assert bytearray(as_stolen_buf) == b"test"
