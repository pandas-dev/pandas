"""
Module that contains many useful utilities
for validating data or function arguments
"""


def validate_args(args, min_length=0, max_length=None, msg=""):
    """
    Checks whether the length of the `*args` argument passed into a function
    has at least `min_length` arguments. If `max_length` is an integer, checks
    whether `*args` has at most `max_length` arguments inclusive. Raises a
    ValueError if any of the aforementioned conditions are False.

    Parameters
    ----------
    args: tuple
        The `*args` parameter passed into a function

    min_length: int, optional
        The minimum number of arguments that should be contained in the `args`.
        tuple. This number must be non-negative. The default is '0'.

    max_length: int, optional
        If not `None`, the maximum number of arguments that should be contained
        in the `args` parameter. This number must be at least as large as the
        provided `min_length` value. The default is None.

    msg: str, optional
        Error message to display when a custom check of args fails. For
        example, pandas does not support a non-None argument for `out`
        when rounding a `Series` or `DataFrame` object. `msg` in this
        case can be "Inplace rounding is not supported".

    Raises
    ------
    ValueError if `args` fails to have a length that is at least `min_length`
    and at most `max_length` inclusive (provided `max_length` is not None)

    """
    length = len(args)

    if min_length < 0:
        raise ValueError("'min_length' must be non-negative")

    if max_length is None:
        if length < min_length:
            raise ValueError(("expected at least {min_length} arguments "
                              "but got {length} arguments instead".
                              format(min_length=min_length, length=length)))

    if min_length > max_length:
        raise ValueError("'min_length' > 'max_length'")

    if (length < min_length) or (length > max_length):
        raise ValueError(("expected between {min_length} and {max_length} "
                          "arguments inclusive but got {length} arguments "
                          "instead".format(min_length=min_length,
                                           length=length,
                                           max_length=max_length)))

    # See gh-12600; this is to allow compatibility with NumPy,
    # which passes in an 'out' parameter as a positional argument
    if args:
        args = list(filter(lambda elt: elt is not None, args))

        if args:
            raise ValueError(msg)


def validate_kwargs(fname, kwargs, *compat_args):
    """
    Checks whether parameters passed to the **kwargs argument in a
    function 'fname' are valid parameters as specified in *compat_args

    Parameters
    ----------
    fname: str
        The name of the function being passed the `**kwargs` parameter

    kwargs: dict
        The `**kwargs` parameter passed into `fname`

    compat_args: *args
        A tuple of keys that `kwargs` is allowed to have

    Raises
    ------
    ValueError if `kwargs` contains keys not in `compat_args`

    """
    list(map(kwargs.__delitem__, filter(
        kwargs.__contains__, compat_args)))
    if kwargs:
        bad_arg = list(kwargs)[0]  # first 'key' element
        raise TypeError(("{fname}() got an unexpected "
                         "keyword argument '{arg}'".
                         format(fname=fname, arg=bad_arg)))
