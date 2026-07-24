# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from ipywidgets import Widget
import ipywidgets.widgets.widget

# The new comm package is not available in our Python 3.7 CI (older ipykernel version)
try:
    import comm
    NEW_COMM_PACKAGE = True
except ImportError:
    NEW_COMM_PACKAGE = False

import ipykernel.comm
import pytest

class DummyComm():
    comm_id = 'a-b-c-d'
    kernel = 'Truthy'

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.messages = []

    def open(self, *args, **kwargs):
        pass

    def on_msg(self, *args, **kwargs):
        pass

    def send(self, *args, **kwargs):
        self.messages.append((args, kwargs))

    def close(self, *args, **kwargs):
        pass


def dummy_create_comm(**kwargs):
    return DummyComm()


def dummy_get_comm_manager(**kwargs):
    return {}


_widget_attrs = {}
undefined = object()

if NEW_COMM_PACKAGE:
    orig_comm = ipykernel.comm.comm.BaseComm
else:
    orig_comm = ipykernel.comm.Comm
orig_create_comm = None
orig_get_comm_manager = None

if NEW_COMM_PACKAGE:
    orig_create_comm = comm.create_comm
    orig_get_comm_manager = comm.get_comm_manager

def setup_test_comm():
    if NEW_COMM_PACKAGE:
        comm.create_comm = dummy_create_comm
        comm.get_comm_manager = dummy_get_comm_manager
        ipykernel.comm.comm.BaseComm = DummyComm
    else:
        ipykernel.comm.Comm = DummyComm
    Widget.comm.klass = DummyComm
    ipywidgets.widgets.widget.Comm = DummyComm
    _widget_attrs['_repr_mimebundle_'] = Widget._repr_mimebundle_
    def raise_not_implemented(*args, **kwargs):
        raise NotImplementedError()
    Widget._repr_mimebundle_ = raise_not_implemented

def teardown_test_comm():
    if NEW_COMM_PACKAGE:
        comm.create_comm = orig_create_comm
        comm.get_comm_manager = orig_get_comm_manager
        ipykernel.comm.comm.BaseComm = orig_comm
    else:
        ipykernel.comm.Comm = orig_comm
    Widget.comm.klass = orig_comm
    ipywidgets.widgets.widget.Comm = orig_comm
    for attr, value in _widget_attrs.items():
        if value is undefined:
            delattr(Widget, attr)
        else:
            setattr(Widget, attr, value)
    _widget_attrs.clear()

@pytest.fixture(autouse=True)
def setup():
    setup_test_comm()
    yield
    teardown_test_comm()

def call_method(method, *args, **kwargs):
    method(*args, **kwargs)
