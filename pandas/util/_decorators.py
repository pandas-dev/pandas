from functools import wraps
import inspect
from textwrap import dedent
from typing import Any, Callable, Dict, List, Optional, Tuple, Type, Union
import warnings

from pandas._libs.properties import cache_readonly  # noqa


def deprecate(
    name: str,
    alternative: Callable,
    version: str,
    alt_name: Optional[str] = None,
    klass: Optional[Type[Warning]] = None,
    stacklevel: int = 2,
    msg: Optional[str] = None,
) -> Callable:
    """
    Return a new function that emits a deprecation warning on use.

    To use this method for a deprecated function, another function
    `alternative` with the same signature must exist. The deprecated
    function will emit a deprecation warning, and in the docstring
    it will contain the deprecation directive with the provided version
    so it can be detected for future removal.

    Parameters
    ----------
    name : str
        Name of function to deprecate.
    alternative : func
        Function to use instead.
    version : str
        Version of pandas in which the method has been deprecated.
    alt_name : str, optional
        Name to use in preference of alternative.__name__.
    klass : Warning, default FutureWarning
    stacklevel : int, default 2
    msg : str
        The message to display in the warning.
        Default is '{name} is deprecated. Use {alt_name} instead.'
    """

    alt_name = alt_name or alternative.__name__
    klass = klass or FutureWarning
    warning_msg = msg or "{} is deprecated, use {} instead".format(name, alt_name)

    @wraps(alternative)
    def wrapper(*args, **kwargs):
        warnings.warn(warning_msg, klass, stacklevel=stacklevel)
        return alternative(*args, **kwargs)

    # adding deprecated directive to the docstring
    msg = msg or "Use `{alt_name}` instead.".format(alt_name=alt_name)
    doc_error_msg = (
        "deprecate needs a correctly formatted docstring in "
        "the target function (should have a one liner short "
        "summary, and opening quotes should be in their own "
        "line). Found:\n{}".format(alternative.__doc__)
    )

    # when python is running in optimized mode (i.e. `-OO`), docstrings are
    # removed, so we check that a docstring with correct formatting is used
    # but we allow empty docstrings
    if alternative.__doc__:
        if alternative.__doc__.count("\n") < 3:
            raise AssertionError(doc_error_msg)
        empty1, summary, empty2, doc = alternative.__doc__.split("\n", 3)
        if empty1 or empty2 and not summary:
            raise AssertionError(doc_error_msg)
        wrapper.__doc__ = dedent(
            """
        {summary}

        .. deprecated:: {depr_version}
            {depr_msg}

        {rest_of_docstring}"""
        ).format(
            summary=summary.strip(),
            depr_version=version,
            depr_msg=msg,
            rest_of_docstring=dedent(doc),
        )

    return wrapper


def deprecate_kwarg(
    old_arg_name: str,
    new_arg_name: Optional[str],
    mapping: Optional[Union[Dict, Callable[[Any], Any]]] = None,
    stacklevel: int = 2,
) -> Callable:
    """
    Decorator to deprecate a keyword argument of a function.

    Parameters
    ----------
    old_arg_name : str
        Name of argument in function to deprecate
    new_arg_name : str or None
        Name of preferred argument in function. Use None to raise warning that
        ``old_arg_name`` keyword is deprecated.
    mapping : dict or callable
        If mapping is present, use it to translate old arguments to
        new arguments. A callable must do its own value checking;
        values not found in a dict will be forwarded unchanged.

    Examples
    --------
    The following deprecates 'cols', using 'columns' instead

    >>> @deprecate_kwarg(old_arg_name='cols', new_arg_name='columns')
    ... def f(columns=''):
    ...     print(columns)
    ...
    >>> f(columns='should work ok')
    should work ok

    >>> f(cols='should raise warning')
    FutureWarning: cols is deprecated, use columns instead
      warnings.warn(msg, FutureWarning)
    should raise warning

    >>> f(cols='should error', columns="can\'t pass do both")
    TypeError: Can only specify 'cols' or 'columns', not both

    >>> @deprecate_kwarg('old', 'new', {'yes': True, 'no': False})
    ... def f(new=False):
    ...     print('yes!' if new else 'no!')
    ...
    >>> f(old='yes')
    FutureWarning: old='yes' is deprecated, use new=True instead
      warnings.warn(msg, FutureWarning)
    yes!

    To raise a warning that a keyword will be removed entirely in the future

    >>> @deprecate_kwarg(old_arg_name='cols', new_arg_name=None)
    ... def f(cols='', another_param=''):
    ...     print(cols)
    ...
    >>> f(cols='should raise warning')
    FutureWarning: the 'cols' keyword is deprecated and will be removed in a
    future version please takes steps to stop use of 'cols'
    should raise warning
    >>> f(another_param='should not raise warning')
    should not raise warning

    >>> f(cols='should raise warning', another_param='')
    FutureWarning: the 'cols' keyword is deprecated and will be removed in a
    future version please takes steps to stop use of 'cols'
    should raise warning
    """

    if mapping is not None and not hasattr(mapping, "get") and not callable(mapping):
        raise TypeError(
            "mapping from old to new argument values " "must be dict or callable!"
        )

    def _deprecate_kwarg(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            old_arg_value = kwargs.pop(old_arg_name, None)

            if new_arg_name is None and old_arg_value is not None:
                msg = (
                    "the '{old_name}' keyword is deprecated and will be "
                    "removed in a future version. "
                    "Please take steps to stop the use of '{old_name}'"
                ).format(old_name=old_arg_name)
                warnings.warn(msg, FutureWarning, stacklevel=stacklevel)
                kwargs[old_arg_name] = old_arg_value
                return func(*args, **kwargs)

            if old_arg_value is not None:
                if mapping is not None:
                    if hasattr(mapping, "get"):
                        new_arg_value = mapping.get(old_arg_value, old_arg_value)
                    else:
                        new_arg_value = mapping(old_arg_value)
                    msg = (
                        "the {old_name}={old_val!r} keyword is deprecated, "
                        "use {new_name}={new_val!r} instead"
                    ).format(
                        old_name=old_arg_name,
                        old_val=old_arg_value,
                        new_name=new_arg_name,
                        new_val=new_arg_value,
                    )
                else:
                    new_arg_value = old_arg_value
                    msg = (
                        "the '{old_name}' keyword is deprecated, "
                        "use '{new_name}' instead"
                    ).format(old_name=old_arg_name, new_name=new_arg_name)

                warnings.warn(msg, FutureWarning, stacklevel=stacklevel)
                if kwargs.get(new_arg_name, None) is not None:
                    msg = (
                        "Can only specify '{old_name}' or '{new_name}', " "not both"
                    ).format(old_name=old_arg_name, new_name=new_arg_name)
                    raise TypeError(msg)
                else:
                    kwargs[new_arg_name] = new_arg_value
            return func(*args, **kwargs)

        return wrapper

    return _deprecate_kwarg


def rewrite_axis_style_signature(
    name: str, extra_params: List[Tuple[str, Any]]
) -> Callable:
    def decorate(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
        params = [
            inspect.Parameter("self", kind),
            inspect.Parameter(name, kind, default=None),
            inspect.Parameter("index", kind, default=None),
            inspect.Parameter("columns", kind, default=None),
            inspect.Parameter("axis", kind, default=None),
        ]

        for pname, default in extra_params:
            params.append(inspect.Parameter(pname, kind, default=default))

        sig = inspect.Signature(params)

        func.__signature__ = sig
        return wrapper

    return decorate


# Substitution and Appender are derived from matplotlib.docstring (1.1.0)
# module http://matplotlib.org/users/license.html


class Substitution:
    """
    A decorator to take a function's docstring and perform string
    substitution on it.

    This decorator should be robust even if func.__doc__ is None
    (for example, if -OO was passed to the interpreter)

    Usage: construct a docstring.Substitution with a sequence or
    dictionary suitable for performing substitution; then
    decorate a suitable function with the constructed object. e.g.

    sub_author_name = Substitution(author='Jason')

    @sub_author_name
    def some_function(x):
        "%(author)s wrote this function"

    # note that some_function.__doc__ is now "Jason wrote this function"

    One can also use positional arguments.

    sub_first_last_names = Substitution('Edgar Allen', 'Poe')

    @sub_first_last_names
    def some_function(x):
        "%s %s wrote the Raven"
    """

    def __init__(self, *args, **kwargs):
        if args and kwargs:
            raise AssertionError("Only positional or keyword args are allowed")

        self.params = args or kwargs

    def __call__(self, func: Callable) -> Callable:
        func.__doc__ = func.__doc__ and func.__doc__ % self.params
        return func

    def update(self, *args, **kwargs) -> None:
        """
        Update self.params with supplied args.

        If called, we assume self.params is a dict.
        """

        self.params.update(*args, **kwargs)


class Appender:
    """
    A function decorator that will append an addendum to the docstring
    of the target function.

    This decorator should be robust even if func.__doc__ is None
    (for example, if -OO was passed to the interpreter).

    Usage: construct a docstring.Appender with a string to be joined to
    the original docstring. An optional 'join' parameter may be supplied
    which will be used to join the docstring and addendum. e.g.

    add_copyright = Appender("Copyright (c) 2009", join='\n')

    @add_copyright
    def my_dog(has='fleas'):
        "This docstring will have a copyright below"
        pass
    """

    def __init__(self, addendum: Optional[str], join: str = "", indents: int = 0):
        if indents > 0:
            self.addendum = indent(addendum, indents=indents)  # type: Optional[str]
        else:
            self.addendum = addendum
        self.join = join

    def __call__(self, func: Callable) -> Callable:
        func.__doc__ = func.__doc__ if func.__doc__ else ""
        self.addendum = self.addendum if self.addendum else ""
        docitems = [func.__doc__, self.addendum]
        func.__doc__ = dedent(self.join.join(docitems))
        return func


def indent(text: Optional[str], indents: int = 1) -> str:
    if not text or not isinstance(text, str):
        return ""
    jointext = "".join(["\n"] + ["    "] * indents)
    return jointext.join(text.split("\n"))
