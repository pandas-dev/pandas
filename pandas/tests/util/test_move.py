# -*- coding: utf-8 -*-
import sys
from uuid import uuid4

import pytest

from pandas.compat import PY3, intern
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


@pytest.mark.skipif(PY3, reason="bytes objects cannot be interned in PY3")
def test_interned():
    salt = uuid4().hex

    def make_string():
        # We need to actually create a new string so that it has refcount
        # one. We use a uuid so that we know the string could not already
        # be in the intern table.
        return "".join(("testing: ", salt))

    # This should work, the string has one reference on the stack.
    move_into_mutable_buffer(make_string())
    refcount = [None]  # nonlocal

    def ref_capture(ob):
        # Subtract two because those are the references owned by this frame:
        #   1. The local variables of this stack frame.
        #   2. The python data stack of this stack frame.
        refcount[0] = sys.getrefcount(ob) - 2
        return ob

    with pytest.raises(BadMove, match="testing"):
        # If we intern the string, it will still have one reference. Now,
        # it is in the intern table, so if other people intern the same
        # string while the mutable buffer holds the first string they will
        # be the same instance.
        move_into_mutable_buffer(ref_capture(intern(make_string())))  # noqa

    assert refcount[0] == 1
