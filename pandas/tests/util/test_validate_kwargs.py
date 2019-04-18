# -*- coding: utf-8 -*-
from collections import OrderedDict

import pytest

from pandas.util._validators import validate_bool_kwarg, validate_kwargs

_fname = "func"


def test_bad_kwarg():
    good_arg = "f"
    bad_arg = good_arg + "o"

    compat_args = OrderedDict()
    compat_args[good_arg] = "foo"
    compat_args[bad_arg + "o"] = "bar"
    kwargs = {good_arg: "foo", bad_arg: "bar"}

    msg = (r"{fname}\(\) got an unexpected "
           r"keyword argument '{arg}'".format(fname=_fname, arg=bad_arg))

    with pytest.raises(TypeError, match=msg):
        validate_kwargs(_fname, kwargs, compat_args)


@pytest.mark.parametrize("i", range(1, 3))
def test_not_all_none(i):
    bad_arg = "foo"
    msg = (r"the '{arg}' parameter is not supported "
           r"in the pandas implementation of {func}\(\)".
           format(arg=bad_arg, func=_fname))

    compat_args = OrderedDict()
    compat_args["foo"] = 1
    compat_args["bar"] = "s"
    compat_args["baz"] = None

    kwarg_keys = ("foo", "bar", "baz")
    kwarg_vals = (2, "s", None)

    kwargs = dict(zip(kwarg_keys[:i], kwarg_vals[:i]))

    with pytest.raises(ValueError, match=msg):
        validate_kwargs(_fname, kwargs, compat_args)


def test_validation():
    # No exceptions should be raised.
    compat_args = OrderedDict()
    compat_args["f"] = None
    compat_args["b"] = 1
    compat_args["ba"] = "s"

    kwargs = dict(f=None, b=1)
    validate_kwargs(_fname, kwargs, compat_args)


@pytest.mark.parametrize("name", ["inplace", "copy"])
@pytest.mark.parametrize("value", [1, "True", [1, 2, 3], 5.0])
def test_validate_bool_kwarg_fail(name, value):
    msg = ("For argument \"%s\" expected type bool, received type %s" %
           (name, type(value).__name__))

    with pytest.raises(ValueError, match=msg):
        validate_bool_kwarg(value, name)


@pytest.mark.parametrize("name", ["inplace", "copy"])
@pytest.mark.parametrize("value", [True, False, None])
def test_validate_bool_kwarg(name, value):
    assert validate_bool_kwarg(value, name) == value
