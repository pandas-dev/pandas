# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import Bool, Tuple, List

from .utils import setup

from ..widget import Widget

from ..._version import __control_protocol_version__

# A widget with simple traits
class SimpleWidget(Widget):
    a = Bool().tag(sync=True)
    b = Tuple(Bool(), Bool(), Bool(), default_value=(False, False, False)).tag(
        sync=True
    )
    c = List(Bool()).tag(sync=True)


def test_empty_send_state():
    w = SimpleWidget()
    w.send_state([])
    assert w.comm.messages == []


def test_empty_hold_sync():
    w = SimpleWidget()
    with w.hold_sync():
        pass
    assert w.comm.messages == []
