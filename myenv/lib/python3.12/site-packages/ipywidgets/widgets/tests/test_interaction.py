# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Test interact and interactive."""

from unittest.mock import patch

import os
from collections import OrderedDict
import pytest

import ipywidgets as widgets

from traitlets import TraitError, Float
from ipywidgets import (interact, interact_manual, interactive,
                        interaction, Output, Widget)

#-----------------------------------------------------------------------------
# Utility stuff
#-----------------------------------------------------------------------------

from .utils import setup, teardown

def f(**kwargs):
    pass

displayed = []
@pytest.fixture()
def clear_display():
    global displayed
    displayed = []

def record_display(*args):
    displayed.extend(args)

#-----------------------------------------------------------------------------
# Actual tests
#-----------------------------------------------------------------------------

def check_widget(w, **d):
    """Check a single widget against a dict"""
    for attr, expected in d.items():
        if attr == 'cls':
            assert w.__class__ is expected
        else:
            value = getattr(w, attr)
            assert value == expected, "{}.{} = {!r} != {!r}".format(w.__class__.__name__, attr, value, expected)

            # For numeric values, the types should match too
            if isinstance(value, (int, float)):
                tv = type(value)
                te = type(expected)
                assert tv is te, "type({}.{}) = {!r} != {!r}".format(w.__class__.__name__, attr, tv, te)

def check_widget_children(container, **to_check):
    """Check that widgets are created as expected"""
    # build a widget dictionary, so it matches
    widgets = {}
    for w in container.children:
        if not isinstance(w, Output):
            widgets[w.description] = w

    for key, d in to_check.items():
        assert key in widgets
        check_widget(widgets[key], **d)


def test_single_value_string():
    a = 'hello'
    c = interactive(f, a=a)
    w = c.children[0]
    check_widget(w,
        cls=widgets.Text,
        description='a',
        value=a,
    )

def test_single_value_bool():
    for a in (True, False):
        c = interactive(f, a=a)
        w = c.children[0]
        check_widget(w,
            cls=widgets.Checkbox,
            description='a',
            value=a,
        )

def test_single_value_float():
    for a in (2.25, 1.0, -3.5, 0.0):
        if not a:
            expected_min = 0.0
            expected_max = 1.0
        elif a > 0:
            expected_min = -a
            expected_max = 3*a
        else:
            expected_min = 3*a
            expected_max = -a
        c = interactive(f, a=a)
        w = c.children[0]
        check_widget(w,
            cls=widgets.FloatSlider,
            description='a',
            value=a,
            min=expected_min,
            max=expected_max,
            step=0.1,
            readout=True,
        )

def test_single_value_int():
    for a in (1, 5, -3, 0):
        if not a:
            expected_min = 0
            expected_max = 1
        elif a > 0:
            expected_min = -a
            expected_max = 3*a
        else:
            expected_min = 3*a
            expected_max = -a
        c = interactive(f, a=a)
        assert len(c.children) == 2
        w = c.children[0]
        check_widget(w,
            cls=widgets.IntSlider,
            description='a',
            value=a,
            min=expected_min,
            max=expected_max,
            step=1,
            readout=True,
        )

def test_list_str():
    values = ['hello', 'there', 'guy']
    first = values[0]
    c = interactive(f, lis=values)
    assert len(c.children) == 2
    d = dict(
        cls=widgets.Dropdown,
        value=first,
        options=tuple(values),
        _options_labels=tuple(values),
        _options_values=tuple(values),
    )
    check_widget_children(c, lis=d)

def test_list_int():
    values = [3, 1, 2]
    first = values[0]
    c = interactive(f, lis=values)
    assert len(c.children) == 2
    d = dict(
        cls=widgets.Dropdown,
        value=first,
        options=tuple(values),
        _options_labels=tuple(str(v) for v in values),
        _options_values=tuple(values),
    )
    check_widget_children(c, lis=d)

def test_list_tuple():
    values = [(3, 300), (1, 100), (2, 200)]
    first = values[0][1]
    c = interactive(f, lis=values)
    assert len(c.children) == 2
    d = dict(
        cls=widgets.Dropdown,
        value=first,
        options=tuple(values),
        _options_labels=("3", "1", "2"),
        _options_values=(300, 100, 200),
    )
    check_widget_children(c, lis=d)

def test_list_tuple_invalid():
    for bad in [
        (),
    ]:
        with pytest.raises(ValueError):
            print(bad) # because there is no custom message in assert_raises
            c = interactive(f, tup=bad)

def test_dict():
    for d in [
        dict(a=5),
        dict(a=5, b='b', c=dict),
    ]:
        c = interactive(f, d=d)
        w = c.children[0]
        check = dict(
            cls=widgets.Dropdown,
            description='d',
            value=next(iter(d.values())),
            options=d,
            _options_labels=tuple(d.keys()),
            _options_values=tuple(d.values()),
        )
        check_widget(w, **check)

def test_ordereddict():
    from collections import OrderedDict
    items = [(3, 300), (1, 100), (2, 200)]
    first = items[0][1]
    values = OrderedDict(items)
    c = interactive(f, lis=values)
    assert len(c.children) == 2
    d = dict(
        cls=widgets.Dropdown,
        value=first,
        options=values,
        _options_labels=("3", "1", "2"),
        _options_values=(300, 100, 200),
    )
    check_widget_children(c, lis=d)

def test_iterable():
    def yield_values():
        yield 3
        yield 1
        yield 2
    first = next(yield_values())
    c = interactive(f, lis=yield_values())
    assert len(c.children) == 2
    d = dict(
        cls=widgets.Dropdown,
        value=first,
        options=(3, 1, 2),
        _options_labels=("3", "1", "2"),
        _options_values=(3, 1, 2),
    )
    check_widget_children(c, lis=d)

def test_iterable_tuple():
    values = [(3, 300), (1, 100), (2, 200)]
    first = values[0][1]
    c = interactive(f, lis=iter(values))
    assert len(c.children) == 2
    d = dict(
        cls=widgets.Dropdown,
        value=first,
        options=tuple(values),
        _options_labels=("3", "1", "2"),
        _options_values=(300, 100, 200),
    )
    check_widget_children(c, lis=d)

def test_mapping():
    from collections.abc import Mapping
    from collections import OrderedDict
    class TestMapping(Mapping):
        def __init__(self, values):
            self.values = values
        def __getitem__(self):
            raise NotImplementedError
        def __len__(self):
            raise NotImplementedError
        def __iter__(self):
            raise NotImplementedError
        def items(self):
            return self.values

    items = [(3, 300), (1, 100), (2, 200)]
    first = items[0][1]
    values = TestMapping(items)
    c = interactive(f, lis=values)
    assert len(c.children) == 2
    d = dict(
        cls=widgets.Dropdown,
        value=first,
        options=tuple(items),
        _options_labels=("3", "1", "2"),
        _options_values=(300, 100, 200),
    )
    check_widget_children(c, lis=d)

def test_decorator_kwarg(clear_display):
    with patch.object(interaction, 'display', record_display):
        @interact(a=5)
        def foo(a):
            pass
    assert len(displayed) == 1
    w = displayed[0].children[0]
    check_widget(w,
        cls=widgets.IntSlider,
        value=5,
    )

def test_interact_instancemethod(clear_display):
    class Foo:
        def show(self, x):
            print(x)

    f = Foo()

    with patch.object(interaction, 'display', record_display):
        g = interact(f.show, x=(1,10))
    assert len(displayed) == 1
    w = displayed[0].children[0]
    check_widget(w,
        cls=widgets.IntSlider,
        value=5,
    )

def test_decorator_no_call(clear_display):
    with patch.object(interaction, 'display', record_display):
        @interact
        def foo(a='default'):
            pass
    assert len(displayed) == 1
    w = displayed[0].children[0]
    check_widget(w,
        cls=widgets.Text,
        value='default',
    )

def test_call_interact(clear_display):
    def foo(a='default'):
        pass
    with patch.object(interaction, 'display', record_display):
        ifoo = interact(foo)
    assert len(displayed) == 1
    w = displayed[0].children[0]
    check_widget(w,
        cls=widgets.Text,
        value='default',
    )

def test_call_interact_on_trait_changed_none_return(clear_display):
    def foo(a='default'):
        pass
    with patch.object(interaction, 'display', record_display):
        ifoo = interact(foo)
    assert len(displayed) == 1
    w = displayed[0].children[0]
    check_widget(w,
        cls=widgets.Text,
        value='default',
    )
    with patch.object(interaction, 'display', record_display):
        w.value = 'called'
    assert len(displayed) == 1

def test_call_interact_kwargs(clear_display):
    def foo(a='default'):
        pass
    with patch.object(interaction, 'display', record_display):
        ifoo = interact(foo, a=10)
    assert len(displayed) == 1
    w = displayed[0].children[0]
    check_widget(w,
        cls=widgets.IntSlider,
        value=10,
    )

def test_call_decorated_on_trait_change(clear_display):
    """test calling @interact decorated functions"""
    d = {}
    with patch.object(interaction, 'display', record_display):
        @interact
        def foo(a='default'):
            d['a'] = a
            return a
    assert len(displayed) == 2 # display the result and the interact
    w = displayed[1].children[0]
    check_widget(w,
        cls=widgets.Text,
        value='default',
    )
    # test calling the function directly
    a = foo('hello')
    assert a == 'hello'
    assert d['a'] == 'hello'

    # test that setting trait values calls the function
    with patch.object(interaction, 'display', record_display):
        w.value = 'called'
    assert d['a'] == 'called'
    assert len(displayed) == 3
    assert w.value == displayed[-1]

def test_call_decorated_kwargs_on_trait_change(clear_display):
    """test calling @interact(foo=bar) decorated functions"""
    d = {}
    with patch.object(interaction, 'display', record_display):
        @interact(a='kwarg')
        def foo(a='default'):
            d['a'] = a
            return a
    assert len(displayed) == 2 # display the result and the interact
    w = displayed[1].children[0]
    check_widget(w,
        cls=widgets.Text,
        value='kwarg',
    )
    # test calling the function directly
    a = foo('hello')
    assert a == 'hello'
    assert d['a'] == 'hello'

    # test that setting trait values calls the function
    with patch.object(interaction, 'display', record_display):
        w.value = 'called'
    assert d['a'] == 'called'
    assert len(displayed) == 3
    assert w.value == displayed[-1]



def test_fixed():
    c = interactive(f, a=widgets.fixed(5), b='text')
    assert len(c.children) == 2
    w = c.children[0]
    check_widget(w,
        cls=widgets.Text,
        value='text',
        description='b',
    )

def test_default_description():
    c = interactive(f, b='text')
    w = c.children[0]
    check_widget(w,
        cls=widgets.Text,
        value='text',
        description='b',
    )

def test_custom_description():
    d = {}
    def record_kwargs(**kwargs):
        d.clear()
        d.update(kwargs)

    c = interactive(record_kwargs, b=widgets.Text(value='text', description='foo'))
    w = c.children[0]
    check_widget(w,
        cls=widgets.Text,
        value='text',
        description='foo',
    )
    w.value = 'different text'
    assert d == {'b': 'different text'}

def test_raises_on_non_value_widget():
    """ Test that passing in a non-value widget raises an error """

    class BadWidget(Widget):
        """ A widget that contains a `value` traitlet """
        value = Float()

    with pytest.raises(TypeError, match=".* not a ValueWidget.*"):
        interactive(f, b=BadWidget())

def test_interact_manual_button():
    c = interact.options(manual=True).widget(f)
    w = c.children[0]
    check_widget(w, cls=widgets.Button)

def test_interact_manual_nocall():
    callcount = 0
    def calltest(testarg):
        callcount += 1
    c = interact.options(manual=True)(calltest, testarg=5).widget
    c.children[0].value = 10
    assert callcount == 0

def test_interact_call():
    w = interact.widget(f)
    w.update()

    w = interact_manual.widget(f)
    w.update()

def test_interact_options():
    def f(x):
        return x
    w = interact.options(manual=False).options(manual=True)(f, x=21).widget
    assert w.manual == True

    w = interact_manual.options(manual=False).options()(x=21).widget(f)
    assert w.manual == False

    w = interact(x=21)().options(manual=True)(f).widget
    assert w.manual == True

def test_interact_options_bad():
    with pytest.raises(ValueError):
        interact.options(bad="foo")

def test_int_range_logic():
    irsw = widgets.IntRangeSlider
    w = irsw(value=(2, 4), min=0, max=6)
    check_widget(w, cls=irsw, value=(2, 4), min=0, max=6)
    w.upper = 3
    w.max = 3
    check_widget(w, cls=irsw, value=(2, 3), min=0, max=3)

    w.min = 0
    w.max = 6
    w.lower = 2
    w.upper = 4
    check_widget(w, cls=irsw, value=(2, 4), min=0, max=6)
    w.value = (0, 1) #lower non-overlapping range
    check_widget(w, cls=irsw, value=(0, 1), min=0, max=6)
    w.value = (5, 6) #upper non-overlapping range
    check_widget(w, cls=irsw, value=(5, 6), min=0, max=6)
    w.lower = 2
    check_widget(w, cls=irsw, value=(2, 6), min=0, max=6)

    with pytest.raises(TraitError):
        w.min = 7
    with pytest.raises(TraitError):
        w.max = -1

    w = irsw(min=2, max=3, value=(2, 3))
    check_widget(w, min=2, max=3, value=(2, 3))
    w = irsw(min=100, max=200, value=(125, 175))
    check_widget(w, value=(125, 175))

    with pytest.raises(TraitError):
        irsw(min=2, max=1)


def test_float_range_logic():
    frsw = widgets.FloatRangeSlider
    w = frsw(value=(.2, .4), min=0., max=.6)
    check_widget(w, cls=frsw, value=(.2, .4), min=0., max=.6)

    w.min = 0.
    w.max = .6
    w.lower = .2
    w.upper = .4
    check_widget(w, cls=frsw, value=(.2, .4), min=0., max=.6)
    w.value = (0., .1) #lower non-overlapping range
    check_widget(w, cls=frsw, value=(0., .1), min=0., max=.6)
    w.value = (.5, .6) #upper non-overlapping range
    check_widget(w, cls=frsw, value=(.5, .6), min=0., max=.6)
    w.lower = .2
    check_widget(w, cls=frsw, value=(.2, .6), min=0., max=.6)

    with pytest.raises(TraitError):
        w.min = .7
    with pytest.raises(TraitError):
        w.max = -.1

    w = frsw(min=2, max=3, value=(2.2, 2.5))
    check_widget(w, min=2., max=3.)

    with pytest.raises(TraitError):
        frsw(min=.2, max=.1)


def test_multiple_selection():
    smw = widgets.SelectMultiple

    # degenerate multiple select
    w = smw()
    check_widget(w, value=tuple())

    # don't accept random other value when no options
    with pytest.raises(TraitError):
        w.value = (2,)
    check_widget(w, value=tuple())

    # basic multiple select
    w = smw(options=[(1, 1)], value=[1])
    check_widget(w, cls=smw, value=(1,), options=((1, 1),))

    # don't accept random other value
    with pytest.raises(TraitError):
        w.value = w.value + (2,)
    check_widget(w, value=(1,))

    # change options, which resets value
    w.options = w.options + ((2, 2),)
    check_widget(w, options=((1, 1), (2,2)), value=())

    # change value
    w.value = (1,2)
    check_widget(w, value=(1, 2))

    # dict style
    w.options = {1: 1}
    check_widget(w, options={1:1})

    # updating
    w.options = (1,)
    with pytest.raises(TraitError):
        w.value = (2,)
    check_widget(w, options=(1,) )

def test_interact_noinspect():
    a = 'hello'
    c = interactive(dict, a=a)
    w = c.children[0]
    check_widget(w,
        cls=widgets.Text,
        description='a',
        value=a,
    )


def test_get_interact_value():
    from ipywidgets.widgets import ValueWidget
    from traitlets import Unicode
    class TheAnswer(ValueWidget):
        _model_name = Unicode('TheAnswer')
        description = Unicode()
        def get_interact_value(self):
            return 42
    w = TheAnswer()
    c = interactive(lambda v: v, v=w)
    c.update()
    assert c.result == 42

def test_state_schema():
    from ipywidgets.widgets import IntSlider, Widget
    import json
    import jsonschema
    s = IntSlider()
    state = Widget.get_manager_state(drop_defaults=True)
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../', 'state.schema.json')) as f:
        schema = json.load(f)
    jsonschema.validate(state, schema)
