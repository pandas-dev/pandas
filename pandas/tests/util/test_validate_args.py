# -*- coding: utf-8 -*-
from collections import OrderedDict

import pytest

from pandas.util._validators import validate_args

_fname = "func"


def test_bad_min_fname_arg_count():
    msg = "'max_fname_arg_count' must be non-negative"

    with pytest.raises(ValueError, match=msg):
        validate_args(_fname, (None,), -1, "foo")


def test_bad_arg_length_max_value_single():
    args = (None, None)
    compat_args = ("foo",)

    min_fname_arg_count = 0
    max_length = len(compat_args) + min_fname_arg_count
    actual_length = len(args) + min_fname_arg_count
    msg = (r"{fname}\(\) takes at most {max_length} "
           r"argument \({actual_length} given\)"
           .format(fname=_fname, max_length=max_length,
                   actual_length=actual_length))

    with pytest.raises(TypeError, match=msg):
        validate_args(_fname, args, min_fname_arg_count, compat_args)


def test_bad_arg_length_max_value_multiple():
    args = (None, None)
    compat_args = dict(foo=None)

    min_fname_arg_count = 2
    max_length = len(compat_args) + min_fname_arg_count
    actual_length = len(args) + min_fname_arg_count
    msg = (r"{fname}\(\) takes at most {max_length} "
           r"arguments \({actual_length} given\)"
           .format(fname=_fname, max_length=max_length,
                   actual_length=actual_length))

    with pytest.raises(TypeError, match=msg):
        validate_args(_fname, args, min_fname_arg_count, compat_args)


@pytest.mark.parametrize("i", range(1, 3))
def test_not_all_defaults(i):
    bad_arg = "foo"
    msg = ("the '{arg}' parameter is not supported "
           r"in the pandas implementation of {func}\(\)".
           format(arg=bad_arg, func=_fname))

    compat_args = OrderedDict()
    compat_args["foo"] = 2
    compat_args["bar"] = -1
    compat_args["baz"] = 3

    arg_vals = (1, -1, 3)

    with pytest.raises(ValueError, match=msg):
        validate_args(_fname, arg_vals[:i], 2, compat_args)


def test_validation():
    # No exceptions should be raised.
    validate_args(_fname, (None,), 2, dict(out=None))

    compat_args = OrderedDict()
    compat_args["axis"] = 1
    compat_args["out"] = None

    validate_args(_fname, (1, None), 2, compat_args)
