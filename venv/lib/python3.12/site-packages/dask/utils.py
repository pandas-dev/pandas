from __future__ import annotations

import codecs
import functools
import gc
import inspect
import os
import re
import shutil
import sys
import tempfile
import types
import uuid
import warnings
from collections.abc import Callable, Hashable, Iterable, Iterator, Mapping, Set
from contextlib import ContextDecorator, contextmanager, nullcontext, suppress
from datetime import datetime, timedelta
from errno import ENOENT
from functools import wraps
from importlib import import_module
from numbers import Integral, Number
from operator import add
from threading import Lock
from typing import Any, ClassVar, Literal, TypeVar, cast, overload
from weakref import WeakValueDictionary

import tlz as toolz

from dask import config
from dask.typing import no_default

K = TypeVar("K")
V = TypeVar("V")
T = TypeVar("T")

# used in decorators to preserve the signature of the function it decorates
# see https://mypy.readthedocs.io/en/stable/generics.html#declaring-decorators
FuncType = Callable[..., Any]
F = TypeVar("F", bound=FuncType)

system_encoding = sys.getdefaultencoding()
if system_encoding == "ascii":
    system_encoding = "utf-8"


def apply(func, args, kwargs=None):
    """Apply a function given its positional and keyword arguments.

    Equivalent to ``func(*args, **kwargs)``
    Most Dask users will never need to use the ``apply`` function.
    It is typically only used by people who need to inject
    keyword argument values into a low level Dask task graph.

    Parameters
    ----------
    func : callable
        The function you want to apply.
    args : tuple
        A tuple containing all the positional arguments needed for ``func``
        (eg: ``(arg_1, arg_2, arg_3)``)
    kwargs : dict, optional
        A dictionary mapping the keyword arguments
        (eg: ``{"kwarg_1": value, "kwarg_2": value}``

    Examples
    --------
    >>> from dask.utils import apply
    >>> def add(number, second_number=5):
    ...     return number + second_number
    ...
    >>> apply(add, (10,), {"second_number": 2})  # equivalent to add(*args, **kwargs)
    12

    >>> task = apply(add, (10,), {"second_number": 2})
    >>> dsk = {'task-name': task}  # adds the task to a low level Dask task graph
    """
    if kwargs:
        return func(*args, **kwargs)
    else:
        return func(*args)


def _deprecated(
    *,
    version: str | None = None,
    after_version: str | None = None,
    message: str | None = None,
    use_instead: str | None = None,
    category: type[Warning] = FutureWarning,
):
    """Decorator to mark a function as deprecated

    Parameters
    ----------
    version : str, optional
        Version of Dask in which the function was deprecated. If specified, the version
        will be included in the default warning message. This should no longer be used
        after the introduction of automated versioning system.
    after_version : str, optional
        Version of Dask after which the function was deprecated. If specified, the
        version will be included in the default warning message.
    message : str, optional
        Custom warning message to raise.
    use_instead : str, optional
        Name of function to use in place of the deprecated function.
        If specified, this will be included in the default warning
        message.
    category : type[Warning], optional
        Type of warning to raise. Defaults to ``FutureWarning``.

    Examples
    --------

    >>> from dask.utils import _deprecated
    >>> @_deprecated(after_version="X.Y.Z", use_instead="bar")
    ... def foo():
    ...     return "baz"
    """

    def decorator(func):
        if message is None:
            msg = f"{func.__name__} "
            if after_version is not None:
                msg += f"was deprecated after version {after_version} "
            elif version is not None:
                msg += f"was deprecated in version {version} "
            else:
                msg += "is deprecated "
            msg += "and will be removed in a future release."

            if use_instead is not None:
                msg += f" Please use {use_instead} instead."
        else:
            msg = message

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            warnings.warn(msg, category=category, stacklevel=2)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def _deprecated_kwarg(
    old_arg_name: str,
    new_arg_name: str | None = None,
    mapping: Mapping[Any, Any] | Callable[[Any], Any] | None = None,
    stacklevel: int = 2,
    comment: str | None = None,
) -> Callable[[F], F]:
    """
    Decorator to deprecate a keyword argument of a function.

    Parameters
    ----------
    old_arg_name : str
        Name of argument in function to deprecate
    new_arg_name : str, optional
        Name of preferred argument in function. Omit to warn that
        ``old_arg_name`` keyword is deprecated.
    mapping : dict or callable, optional
        If mapping is present, use it to translate old arguments to
        new arguments. A callable must do its own value checking;
        values not found in a dict will be forwarded unchanged.
    comment :  str, optional
        Additional message to deprecation message. Useful to pass
        on suggestions with the deprecation warning.

    Examples
    --------
    The following deprecates 'cols', using 'columns' instead

    >>> @_deprecated_kwarg(old_arg_name='cols', new_arg_name='columns')
    ... def f(columns=''):
    ...     print(columns)
    ...
    >>> f(columns='should work ok')
    should work ok

    >>> f(cols='should raise warning')  # doctest: +SKIP
    FutureWarning: cols is deprecated, use columns instead
      warnings.warn(msg, FutureWarning)
    should raise warning

    >>> f(cols='should error', columns="can\'t pass do both")  # doctest: +SKIP
    TypeError: Can only specify 'cols' or 'columns', not both

    >>> @_deprecated_kwarg('old', 'new', {'yes': True, 'no': False})
    ... def f(new=False):
    ...     print('yes!' if new else 'no!')
    ...
    >>> f(old='yes')  # doctest: +SKIP
    FutureWarning: old='yes' is deprecated, use new=True instead
      warnings.warn(msg, FutureWarning)
    yes!

    To raise a warning that a keyword will be removed entirely in the future

    >>> @_deprecated_kwarg(old_arg_name='cols', new_arg_name=None)
    ... def f(cols='', another_param=''):
    ...     print(cols)
    ...
    >>> f(cols='should raise warning')  # doctest: +SKIP
    FutureWarning: the 'cols' keyword is deprecated and will be removed in a
    future version please takes steps to stop use of 'cols'
    should raise warning
    >>> f(another_param='should not raise warning')  # doctest: +SKIP
    should not raise warning

    >>> f(cols='should raise warning', another_param='')  # doctest: +SKIP
    FutureWarning: the 'cols' keyword is deprecated and will be removed in a
    future version please takes steps to stop use of 'cols'
    should raise warning
    """
    if mapping is not None and not hasattr(mapping, "get") and not callable(mapping):
        raise TypeError(
            "mapping from old to new argument values must be dict or callable!"
        )

    comment_ = f"\n{comment}" or ""

    def _deprecated_kwarg(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Callable[..., Any]:
            old_arg_value = kwargs.pop(old_arg_name, no_default)

            if old_arg_value is not no_default:
                if new_arg_name is None:
                    msg = (
                        f"the {repr(old_arg_name)} keyword is deprecated and "
                        "will be removed in a future version. Please take "
                        f"steps to stop the use of {repr(old_arg_name)}"
                    ) + comment_
                    warnings.warn(msg, FutureWarning, stacklevel=stacklevel)
                    kwargs[old_arg_name] = old_arg_value
                    return func(*args, **kwargs)

                elif mapping is not None:
                    if callable(mapping):
                        new_arg_value = mapping(old_arg_value)
                    else:
                        new_arg_value = mapping.get(old_arg_value, old_arg_value)
                    msg = (
                        f"the {old_arg_name}={repr(old_arg_value)} keyword is "
                        "deprecated, use "
                        f"{new_arg_name}={repr(new_arg_value)} instead."
                    )
                else:
                    new_arg_value = old_arg_value
                    msg = (
                        f"the {repr(old_arg_name)} keyword is deprecated, "
                        f"use {repr(new_arg_name)} instead."
                    )

                warnings.warn(msg + comment_, FutureWarning, stacklevel=stacklevel)
                if kwargs.get(new_arg_name) is not None:
                    msg = (
                        f"Can only specify {repr(old_arg_name)} "
                        f"or {repr(new_arg_name)}, not both."
                    )
                    raise TypeError(msg)
                kwargs[new_arg_name] = new_arg_value
            return func(*args, **kwargs)

        return cast(F, wrapper)

    return _deprecated_kwarg


def deepmap(func, *seqs):
    """Apply function inside nested lists

    >>> inc = lambda x: x + 1
    >>> deepmap(inc, [[1, 2], [3, 4]])
    [[2, 3], [4, 5]]

    >>> add = lambda x, y: x + y
    >>> deepmap(add, [[1, 2], [3, 4]], [[10, 20], [30, 40]])
    [[11, 22], [33, 44]]
    """
    if isinstance(seqs[0], (list, Iterator)):
        return [deepmap(func, *items) for items in zip(*seqs)]
    else:
        return func(*seqs)


@_deprecated()
def homogeneous_deepmap(func, seq):
    if not seq:
        return seq
    n = 0
    tmp = seq
    while isinstance(tmp, list):
        n += 1
        tmp = tmp[0]

    return ndeepmap(n, func, seq)


def ndeepmap(n, func, seq):
    """Call a function on every element within a nested container

    >>> def inc(x):
    ...     return x + 1
    >>> L = [[1, 2], [3, 4, 5]]
    >>> ndeepmap(2, inc, L)
    [[2, 3], [4, 5, 6]]
    """
    if n == 1:
        return [func(item) for item in seq]
    elif n > 1:
        return [ndeepmap(n - 1, func, item) for item in seq]
    elif isinstance(seq, list):
        return func(seq[0])
    else:
        return func(seq)


def import_required(mod_name, error_msg):
    """Attempt to import a required dependency.

    Raises a RuntimeError if the requested module is not available.
    """
    try:
        return import_module(mod_name)
    except ImportError as e:
        raise RuntimeError(error_msg) from e


@contextmanager
def tmpfile(extension="", dir=None):
    """
    Function to create and return a unique temporary file with the given extension, if provided.

    Parameters
    ----------
    extension : str
        The extension of the temporary file to be created
    dir : str
        If ``dir`` is not None, the file will be created in that directory; otherwise,
        Python's default temporary directory is used.

    Returns
    -------
    out : str
        Path to the temporary file

    See Also
    --------
    NamedTemporaryFile : Built-in alternative for creating temporary files
    tmp_path : pytest fixture for creating a temporary directory unique to the test invocation

    Notes
    -----
    This context manager is particularly useful on Windows for opening temporary files multiple times.
    """
    extension = extension.lstrip(".")
    if extension:
        extension = "." + extension
    handle, filename = tempfile.mkstemp(extension, dir=dir)
    os.close(handle)
    os.remove(filename)

    try:
        yield filename
    finally:
        if os.path.exists(filename):
            with suppress(OSError):  # sometimes we can't remove a generated temp file
                if os.path.isdir(filename):
                    shutil.rmtree(filename)
                else:
                    os.remove(filename)


@contextmanager
def tmpdir(dir=None):
    """
    Function to create and return a unique temporary directory.

    Parameters
    ----------
    dir : str
        If ``dir`` is not None, the directory will be created in that directory; otherwise,
        Python's default temporary directory is used.

    Returns
    -------
    out : str
        Path to the temporary directory

    Notes
    -----
    This context manager is particularly useful on Windows for opening temporary directories multiple times.
    """
    dirname = tempfile.mkdtemp(dir=dir)

    try:
        yield dirname
    finally:
        if os.path.exists(dirname):
            if os.path.isdir(dirname):
                with suppress(OSError):
                    shutil.rmtree(dirname)
            else:
                with suppress(OSError):
                    os.remove(dirname)


@contextmanager
def filetext(text, extension="", open=open, mode="w"):
    with tmpfile(extension=extension) as filename:
        f = open(filename, mode=mode)
        try:
            f.write(text)
        finally:
            try:
                f.close()
            except AttributeError:
                pass

        yield filename


@contextmanager
def changed_cwd(new_cwd):
    old_cwd = os.getcwd()
    os.chdir(new_cwd)
    try:
        yield
    finally:
        os.chdir(old_cwd)


@contextmanager
def tmp_cwd(dir=None):
    with tmpdir(dir) as dirname:
        with changed_cwd(dirname):
            yield dirname


class IndexCallable:
    """Provide getitem syntax for functions

    >>> def inc(x):
    ...     return x + 1

    >>> I = IndexCallable(inc)
    >>> I[3]
    4
    """

    __slots__ = ("fn",)

    def __init__(self, fn):
        self.fn = fn

    def __getitem__(self, key):
        return self.fn(key)


@contextmanager
def filetexts(d, open=open, mode="t", use_tmpdir=True):
    """Dumps a number of textfiles to disk

    Parameters
    ----------
    d : dict
        a mapping from filename to text like {'a.csv': '1,1\n2,2'}

    Since this is meant for use in tests, this context manager will
    automatically switch to a temporary current directory, to avoid
    race conditions when running tests in parallel.
    """
    with tmp_cwd() if use_tmpdir else nullcontext():
        for filename, text in d.items():
            try:
                os.makedirs(os.path.dirname(filename))
            except OSError:
                pass
            f = open(filename, "w" + mode)
            try:
                f.write(text)
            finally:
                try:
                    f.close()
                except AttributeError:
                    pass

        yield list(d)

        for filename in d:
            if os.path.exists(filename):
                with suppress(OSError):
                    os.remove(filename)


def concrete(seq):
    """Make nested iterators concrete lists

    >>> data = [[1, 2], [3, 4]]
    >>> seq = iter(map(iter, data))
    >>> concrete(seq)
    [[1, 2], [3, 4]]
    """
    if isinstance(seq, Iterator):
        seq = list(seq)
    if isinstance(seq, (tuple, list)):
        seq = list(map(concrete, seq))
    return seq


def pseudorandom(n: int, p, random_state=None):
    """Pseudorandom array of integer indexes

    >>> pseudorandom(5, [0.5, 0.5], random_state=123)
    array([1, 0, 0, 1, 1], dtype=int8)

    >>> pseudorandom(10, [0.5, 0.2, 0.2, 0.1], random_state=5)
    array([0, 2, 0, 3, 0, 1, 2, 1, 0, 0], dtype=int8)
    """
    import numpy as np

    p = list(p)
    cp = np.cumsum([0] + p)
    assert np.allclose(1, cp[-1])
    assert len(p) < 256

    if not isinstance(random_state, np.random.RandomState):
        random_state = np.random.RandomState(random_state)

    x = random_state.random_sample(n)
    out = np.empty(n, dtype="i1")

    for i, (low, high) in enumerate(zip(cp[:-1], cp[1:])):
        out[(x >= low) & (x < high)] = i
    return out


def random_state_data(n: int, random_state=None) -> list:
    """Return a list of arrays that can initialize
    ``np.random.RandomState``.

    Parameters
    ----------
    n : int
        Number of arrays to return.
    random_state : int or np.random.RandomState, optional
        If an int, is used to seed a new ``RandomState``.
    """
    import numpy as np

    if not all(
        hasattr(random_state, attr) for attr in ["normal", "beta", "bytes", "uniform"]
    ):
        random_state = np.random.RandomState(random_state)

    random_data = random_state.bytes(624 * n * 4)  # `n * 624` 32-bit integers
    l = list(np.frombuffer(random_data, dtype="<u4").reshape((n, -1)))
    assert len(l) == n
    return l


def is_integer(i) -> bool:
    """
    >>> is_integer(6)
    True
    >>> is_integer(42.0)
    True
    >>> is_integer('abc')
    False
    """
    return isinstance(i, Integral) or (isinstance(i, float) and i.is_integer())


ONE_ARITY_BUILTINS = {
    abs,
    all,
    any,
    ascii,
    bool,
    bytearray,
    bytes,
    callable,
    chr,
    classmethod,
    complex,
    dict,
    dir,
    enumerate,
    eval,
    float,
    format,
    frozenset,
    hash,
    hex,
    id,
    int,
    iter,
    len,
    list,
    max,
    min,
    next,
    oct,
    open,
    ord,
    range,
    repr,
    reversed,
    round,
    set,
    slice,
    sorted,
    staticmethod,
    str,
    sum,
    tuple,
    type,
    vars,
    zip,
    memoryview,
}
MULTI_ARITY_BUILTINS = {
    compile,
    delattr,
    divmod,
    filter,
    getattr,
    hasattr,
    isinstance,
    issubclass,
    map,
    pow,
    setattr,
}


def getargspec(func):
    """Version of inspect.getargspec that works with partial and warps."""
    if isinstance(func, functools.partial):
        return getargspec(func.func)

    func = getattr(func, "__wrapped__", func)
    if isinstance(func, type):
        return inspect.getfullargspec(func.__init__)
    else:
        return inspect.getfullargspec(func)


def takes_multiple_arguments(func, varargs=True):
    """Does this function take multiple arguments?

    >>> def f(x, y): pass
    >>> takes_multiple_arguments(f)
    True

    >>> def f(x): pass
    >>> takes_multiple_arguments(f)
    False

    >>> def f(x, y=None): pass
    >>> takes_multiple_arguments(f)
    False

    >>> def f(*args): pass
    >>> takes_multiple_arguments(f)
    True

    >>> class Thing:
    ...     def __init__(self, a): pass
    >>> takes_multiple_arguments(Thing)
    False

    """
    if func in ONE_ARITY_BUILTINS:
        return False
    elif func in MULTI_ARITY_BUILTINS:
        return True

    try:
        spec = getargspec(func)
    except Exception:
        return False

    try:
        is_constructor = spec.args[0] == "self" and isinstance(func, type)
    except Exception:
        is_constructor = False

    if varargs and spec.varargs:
        return True

    ndefaults = 0 if spec.defaults is None else len(spec.defaults)
    return len(spec.args) - ndefaults - is_constructor > 1


def get_named_args(func) -> list[str]:
    """Get all non ``*args/**kwargs`` arguments for a function"""
    s = inspect.signature(func)
    return [
        n
        for n, p in s.parameters.items()
        if p.kind in [p.POSITIONAL_OR_KEYWORD, p.POSITIONAL_ONLY, p.KEYWORD_ONLY]
    ]


class Dispatch:
    """Simple single dispatch."""

    def __init__(self, name=None):
        self._lookup = {}
        self._lazy = {}
        if name:
            self.__name__ = name

    def register(self, type, func=None):
        """Register dispatch of `func` on arguments of type `type`"""

        def wrapper(func):
            if isinstance(type, tuple):
                for t in type:
                    self.register(t, func)
            else:
                self._lookup[type] = func
            return func

        return wrapper(func) if func is not None else wrapper

    def register_lazy(self, toplevel, func=None):
        """
        Register a registration function which will be called if the
        *toplevel* module (e.g. 'pandas') is ever loaded.
        """

        def wrapper(func):
            self._lazy[toplevel] = func
            return func

        return wrapper(func) if func is not None else wrapper

    def dispatch(self, cls):
        """Return the function implementation for the given ``cls``"""
        lk = self._lookup
        if cls in lk:
            return lk[cls]
        for cls2 in cls.__mro__:
            # Is a lazy registration function present?
            try:
                toplevel, _, _ = cls2.__module__.partition(".")
            except Exception:
                continue
            try:
                register = self._lazy[toplevel]
            except KeyError:
                pass
            else:
                register()
                self._lazy.pop(toplevel, None)
                meth = self.dispatch(cls)  # recurse
                lk[cls] = meth
                lk[cls2] = meth
                return meth
            try:
                impl = lk[cls2]
            except KeyError:
                pass
            else:
                if cls is not cls2:
                    # Cache lookup
                    lk[cls] = impl
                return impl
        raise TypeError(f"No dispatch for {cls}")

    def __call__(self, arg, *args, **kwargs):
        """
        Call the corresponding method based on type of argument.
        """
        meth = self.dispatch(type(arg))
        return meth(arg, *args, **kwargs)

    @property
    def __doc__(self):
        try:
            func = self.dispatch(object)
            return func.__doc__
        except TypeError:
            return "Single Dispatch for %s" % self.__name__


def ensure_not_exists(filename) -> None:
    """
    Ensure that a file does not exist.
    """
    try:
        os.unlink(filename)
    except OSError as e:
        if e.errno != ENOENT:
            raise


def _skip_doctest(line):
    # NumPy docstring contains cursor and comment only example
    stripped = line.strip()
    if stripped == ">>>" or stripped.startswith(">>> #"):
        return line
    elif ">>>" in stripped and "+SKIP" not in stripped:
        if "# doctest:" in line:
            return line + ", +SKIP"
        else:
            return line + "  # doctest: +SKIP"
    else:
        return line


def skip_doctest(doc):
    if doc is None:
        return ""
    return "\n".join([_skip_doctest(line) for line in doc.split("\n")])


def extra_titles(doc):
    lines = doc.split("\n")
    titles = {
        i: lines[i].strip()
        for i in range(len(lines) - 1)
        if lines[i + 1].strip() and all(c == "-" for c in lines[i + 1].strip())
    }

    seen = set()
    for i, title in sorted(titles.items()):
        if title in seen:
            new_title = "Extra " + title
            lines[i] = lines[i].replace(title, new_title)
            lines[i + 1] = lines[i + 1].replace("-" * len(title), "-" * len(new_title))
        else:
            seen.add(title)

    return "\n".join(lines)


def ignore_warning(doc, cls, name, extra="", skipblocks=0, inconsistencies=None):
    """Expand docstring by adding disclaimer and extra text"""
    import inspect

    if inspect.isclass(cls):
        l1 = "This docstring was copied from {}.{}.{}.\n\n".format(
            cls.__module__,
            cls.__name__,
            name,
        )
    else:
        l1 = f"This docstring was copied from {cls.__name__}.{name}.\n\n"
    l2 = "Some inconsistencies with the Dask version may exist."

    i = doc.find("\n\n")
    if i != -1:
        # Insert our warning
        head = doc[: i + 2]
        tail = doc[i + 2 :]
        while skipblocks > 0:
            i = tail.find("\n\n")
            head = tail[: i + 2]
            tail = tail[i + 2 :]
            skipblocks -= 1
        # Indentation of next line
        indent = re.match(r"\s*", tail).group(0)
        # Insert the warning, indented, with a blank line before and after
        if extra:
            more = [indent, extra.rstrip("\n") + "\n\n"]
        else:
            more = []
        if inconsistencies is not None:
            l3 = f"Known inconsistencies: \n {inconsistencies}"
            bits = [head, indent, l1, l2, "\n\n", l3, "\n\n"] + more + [tail]
        else:
            bits = [head, indent, l1, indent, l2, "\n\n"] + more + [tail]
        doc = "".join(bits)

    return doc


def unsupported_arguments(doc, args):
    """Mark unsupported arguments with a disclaimer"""
    lines = doc.split("\n")
    for arg in args:
        subset = [
            (i, line)
            for i, line in enumerate(lines)
            if re.match(r"^\s*" + arg + " ?:", line)
        ]
        if len(subset) == 1:
            [(i, line)] = subset
            lines[i] = line + "  (Not supported in Dask)"
    return "\n".join(lines)


def _derived_from(
    cls, method, ua_args=None, extra="", skipblocks=0, inconsistencies=None
):
    """Helper function for derived_from to ease testing"""
    ua_args = ua_args or []

    # do not use wraps here, as it hides keyword arguments displayed
    # in the doc
    original_method = getattr(cls, method.__name__)

    doc = getattr(original_method, "__doc__", None)

    if isinstance(original_method, property):
        # some things like SeriesGroupBy.unique are generated.
        original_method = original_method.fget
        if not doc:
            doc = getattr(original_method, "__doc__", None)

    if isinstance(original_method, functools.cached_property):
        original_method = original_method.func
        if not doc:
            doc = getattr(original_method, "__doc__", None)

    if doc is None:
        doc = ""

    # pandas DataFrame/Series sometimes override methods without setting __doc__
    if not doc and cls.__name__ in {"DataFrame", "Series"}:
        for obj in cls.mro():
            obj_method = getattr(obj, method.__name__, None)
            if obj_method is not None and obj_method.__doc__:
                doc = obj_method.__doc__
                break

    # Insert disclaimer that this is a copied docstring
    if doc:
        doc = ignore_warning(
            doc,
            cls,
            method.__name__,
            extra=extra,
            skipblocks=skipblocks,
            inconsistencies=inconsistencies,
        )
    elif extra:
        doc += extra.rstrip("\n") + "\n\n"

    # Mark unsupported arguments
    try:
        method_args = get_named_args(method)
        original_args = get_named_args(original_method)
        not_supported = [m for m in original_args if m not in method_args]
    except ValueError:
        not_supported = []
    if len(ua_args) > 0:
        not_supported.extend(ua_args)
    if len(not_supported) > 0:
        doc = unsupported_arguments(doc, not_supported)

    doc = skip_doctest(doc)
    doc = extra_titles(doc)

    return doc


def derived_from(
    original_klass, version=None, ua_args=None, skipblocks=0, inconsistencies=None
):
    """Decorator to attach original class's docstring to the wrapped method.

    The output structure will be: top line of docstring, disclaimer about this
    being auto-derived, any extra text associated with the method being patched,
    the body of the docstring and finally, the list of keywords that exist in
    the original method but not in the dask version.

    Parameters
    ----------
    original_klass: type
        Original class which the method is derived from
    version : str
        Original package version which supports the wrapped method
    ua_args : list
        List of keywords which Dask doesn't support. Keywords existing in
        original but not in Dask will automatically be added.
    skipblocks : int
        How many text blocks (paragraphs) to skip from the start of the
        docstring. Useful for cases where the target has extra front-matter.
    inconsistencies: list
        List of known inconsistencies with method whose docstrings are being
        copied.
    """
    ua_args = ua_args or []

    def wrapper(method):
        try:
            extra = getattr(method, "__doc__", None) or ""
            method.__doc__ = _derived_from(
                original_klass,
                method,
                ua_args=ua_args,
                extra=extra,
                skipblocks=skipblocks,
                inconsistencies=inconsistencies,
            )
            return method

        except AttributeError:
            module_name = original_klass.__module__.split(".")[0]

            @functools.wraps(method)
            def wrapped(*args, **kwargs):
                msg = f"Base package doesn't support '{method.__name__}'."
                if version is not None:
                    msg2 = " Use {0} {1} or later to use this method."
                    msg += msg2.format(module_name, version)
                raise NotImplementedError(msg)

            return wrapped

    return wrapper


def funcname(func) -> str:
    """Get the name of a function."""
    # functools.partial
    if isinstance(func, functools.partial):
        return funcname(func.func)
    # methodcaller
    if isinstance(func, methodcaller):
        return func.method[:50]

    module_name = getattr(func, "__module__", None) or ""
    type_name = getattr(type(func), "__name__", None) or ""

    # toolz.curry
    if "toolz" in module_name and "curry" == type_name:
        return func.func_name[:50]
    # multipledispatch objects
    if "multipledispatch" in module_name and "Dispatcher" == type_name:
        return func.name[:50]
    # numpy.vectorize objects
    if "numpy" in module_name and "vectorize" == type_name:
        return ("vectorize_" + funcname(func.pyfunc))[:50]

    # All other callables
    try:
        name = func.__name__
        if name == "<lambda>":
            return "lambda"
        return name[:50]
    except AttributeError:
        return str(func)[:50]


def typename(typ: Any, short: bool = False) -> str:
    """
    Return the name of a type

    Examples
    --------
    >>> typename(int)
    'int'

    >>> from dask.core import literal
    >>> typename(literal)
    'dask.core.literal'
    >>> typename(literal, short=True)
    'dask.literal'
    """
    if not isinstance(typ, type):
        return typename(type(typ))
    try:
        if not typ.__module__ or typ.__module__ == "builtins":
            return typ.__name__
        else:
            if short:
                module, *_ = typ.__module__.split(".")
            else:
                module = typ.__module__
            return module + "." + typ.__name__
    except AttributeError:
        return str(typ)


def ensure_bytes(s) -> bytes:
    """Attempt to turn `s` into bytes.

    Parameters
    ----------
    s : Any
        The object to be converted. Will correctly handled
        * str
        * bytes
        * objects implementing the buffer protocol (memoryview, ndarray, etc.)

    Returns
    -------
    b : bytes

    Raises
    ------
    TypeError
        When `s` cannot be converted

    Examples
    --------
    >>> ensure_bytes('123')
    b'123'
    >>> ensure_bytes(b'123')
    b'123'
    >>> ensure_bytes(bytearray(b'123'))
    b'123'
    """
    if isinstance(s, bytes):
        return s
    elif hasattr(s, "encode"):
        return s.encode()
    else:
        try:
            return bytes(s)
        except Exception as e:
            raise TypeError(
                f"Object {s} is neither a bytes object nor can be encoded to bytes"
            ) from e


def ensure_unicode(s) -> str:
    """Turn string or bytes to string

    >>> ensure_unicode('123')
    '123'
    >>> ensure_unicode(b'123')
    '123'
    """
    if isinstance(s, str):
        return s
    elif hasattr(s, "decode"):
        return s.decode()
    else:
        try:
            return codecs.decode(s)
        except Exception as e:
            raise TypeError(
                f"Object {s} is neither a str object nor can be decoded to str"
            ) from e


def digit(n, k, base):
    """

    >>> digit(1234, 0, 10)
    4
    >>> digit(1234, 1, 10)
    3
    >>> digit(1234, 2, 10)
    2
    >>> digit(1234, 3, 10)
    1
    """
    return n // base**k % base


def insert(tup, loc, val):
    """

    >>> insert(('a', 'b', 'c'), 0, 'x')
    ('x', 'b', 'c')
    """
    L = list(tup)
    L[loc] = val
    return tuple(L)


def memory_repr(num):
    for x in ["bytes", "KB", "MB", "GB", "TB"]:
        if num < 1024.0:
            return f"{num:3.1f} {x}"
        num /= 1024.0


def asciitable(columns, rows):
    """Formats an ascii table for given columns and rows.

    Parameters
    ----------
    columns : list
        The column names
    rows : list of tuples
        The rows in the table. Each tuple must be the same length as
        ``columns``.
    """
    rows = [tuple(str(i) for i in r) for r in rows]
    columns = tuple(str(i) for i in columns)
    widths = tuple(max(*map(len, x), len(c)) for x, c in zip(zip(*rows), columns))
    row_template = ("|" + (" %%-%ds |" * len(columns))) % widths
    header = row_template % tuple(columns)
    bar = "+%s+" % "+".join("-" * (w + 2) for w in widths)
    data = "\n".join(row_template % r for r in rows)
    return "\n".join([bar, header, bar, data, bar])


def put_lines(buf, lines):
    if any(not isinstance(x, str) for x in lines):
        lines = [str(x) for x in lines]
    buf.write("\n".join(lines))


_method_cache: dict[str, methodcaller] = {}


class methodcaller:
    """
    Return a callable object that calls the given method on its operand.

    Unlike the builtin `operator.methodcaller`, instances of this class are
    cached and arguments are passed at call time instead of build time.
    """

    __slots__ = ("method",)
    method: str

    @property
    def func(self) -> str:
        # For `funcname` to work
        return self.method

    def __new__(cls, method: str):
        try:
            return _method_cache[method]
        except KeyError:
            self = object.__new__(cls)
            self.method = method
            _method_cache[method] = self
            return self

    def __call__(self, __obj, *args, **kwargs):
        return getattr(__obj, self.method)(*args, **kwargs)

    def __reduce__(self):
        return (methodcaller, (self.method,))

    def __str__(self):
        return f"<{self.__class__.__name__}: {self.method}>"

    __repr__ = __str__


class itemgetter:
    """Variant of operator.itemgetter that supports equality tests"""

    __slots__ = ("index",)

    def __init__(self, index):
        self.index = index

    def __call__(self, x):
        return x[self.index]

    def __reduce__(self):
        return (itemgetter, (self.index,))

    def __eq__(self, other):
        return type(self) is type(other) and self.index == other.index


class MethodCache:
    """Attribute access on this object returns a methodcaller for that
    attribute.

    Examples
    --------
    >>> a = [1, 3, 3]
    >>> M.count(a, 3) == a.count(3)
    True
    """

    def __getattr__(self, item):
        return methodcaller(item)

    def __dir__(self):
        return list(_method_cache)


M = MethodCache()


class SerializableLock:
    """A Serializable per-process Lock

    This wraps a normal ``threading.Lock`` object and satisfies the same
    interface.  However, this lock can also be serialized and sent to different
    processes.  It will not block concurrent operations between processes (for
    this you should look at ``multiprocessing.Lock`` or ``locket.lock_file``
    but will consistently deserialize into the same lock.

    So if we make a lock in one process::

        lock = SerializableLock()

    And then send it over to another process multiple times::

        bytes = pickle.dumps(lock)
        a = pickle.loads(bytes)
        b = pickle.loads(bytes)

    Then the deserialized objects will operate as though they were the same
    lock, and collide as appropriate.

    This is useful for consistently protecting resources on a per-process
    level.

    The creation of locks is itself not threadsafe.
    """

    _locks: ClassVar[WeakValueDictionary[Hashable, Lock]] = WeakValueDictionary()
    token: Hashable
    lock: Lock

    def __init__(self, token: Hashable | None = None):
        self.token = token or str(uuid.uuid4())
        if self.token in SerializableLock._locks:
            self.lock = SerializableLock._locks[self.token]
        else:
            self.lock = Lock()
            SerializableLock._locks[self.token] = self.lock

    def acquire(self, *args, **kwargs):
        return self.lock.acquire(*args, **kwargs)

    def release(self, *args, **kwargs):
        return self.lock.release(*args, **kwargs)

    def __enter__(self):
        self.lock.__enter__()

    def __exit__(self, *args):
        self.lock.__exit__(*args)

    def locked(self):
        return self.lock.locked()

    def __getstate__(self):
        return self.token

    def __setstate__(self, token):
        self.__init__(token)

    def __str__(self):
        return f"<{self.__class__.__name__}: {self.token}>"

    __repr__ = __str__


def get_scheduler_lock(collection=None, scheduler=None):
    """Get an instance of the appropriate lock for a certain situation based on
    scheduler used."""
    from dask import multiprocessing
    from dask.base import get_scheduler

    actual_get = get_scheduler(collections=[collection], scheduler=scheduler)

    if actual_get == multiprocessing.get:
        return multiprocessing.get_context().Manager().Lock()
    else:
        # if this is a distributed client, we need to lock on
        # the level between processes, SerializableLock won't work
        try:
            import distributed.lock
            from distributed.worker import get_client

            client = get_client()
        except (ImportError, ValueError):
            pass
        else:
            if actual_get == client.get:
                return distributed.lock.Lock()

    return SerializableLock()


def ensure_dict(d: Mapping[K, V], *, copy: bool = False) -> dict[K, V]:
    """Convert a generic Mapping into a dict.
    Optimize use case of :class:`~dask.highlevelgraph.HighLevelGraph`.

    Parameters
    ----------
    d : Mapping
    copy : bool
        If True, guarantee that the return value is always a shallow copy of d;
        otherwise it may be the input itself.
    """
    if type(d) is dict:
        return d.copy() if copy else d
    try:
        layers = d.layers  # type: ignore
    except AttributeError:
        return dict(d)

    result = {}
    for layer in toolz.unique(layers.values(), key=id):
        result.update(layer)
    return result


def ensure_set(s: Set[T], *, copy: bool = False) -> set[T]:
    """Convert a generic Set into a set.

    Parameters
    ----------
    s : Set
    copy : bool
        If True, guarantee that the return value is always a shallow copy of s;
        otherwise it may be the input itself.
    """
    if type(s) is set:
        return s.copy() if copy else s
    return set(s)


class OperatorMethodMixin:
    """A mixin for dynamically implementing operators"""

    __slots__ = ()

    @classmethod
    def _bind_operator(cls, op):
        """bind operator to this class"""
        name = op.__name__

        if name.endswith("_"):
            # for and_ and or_
            name = name[:-1]
        elif name == "inv":
            name = "invert"

        meth = f"__{name}__"

        if name in ("abs", "invert", "neg", "pos"):
            setattr(cls, meth, cls._get_unary_operator(op))
        else:
            setattr(cls, meth, cls._get_binary_operator(op))

            if name in ("eq", "gt", "ge", "lt", "le", "ne", "getitem"):
                return

            rmeth = f"__r{name}__"
            setattr(cls, rmeth, cls._get_binary_operator(op, inv=True))

    @classmethod
    def _get_unary_operator(cls, op):
        """Must return a method used by unary operator"""
        raise NotImplementedError

    @classmethod
    def _get_binary_operator(cls, op, inv=False):
        """Must return a method used by binary operator"""
        raise NotImplementedError


def partial_by_order(*args, **kwargs):
    """

    >>> from operator import add
    >>> partial_by_order(5, function=add, other=[(1, 10)])
    15
    """
    function = kwargs.pop("function")
    other = kwargs.pop("other")
    args2 = list(args)
    for i, arg in other:
        args2.insert(i, arg)
    return function(*args2, **kwargs)


def is_arraylike(x) -> bool:
    """Is this object a numpy array or something similar?

    This function tests specifically for an object that already has
    array attributes (e.g. np.ndarray, dask.array.Array, cupy.ndarray,
    sparse.COO), **NOT** for something that can be coerced into an
    array object (e.g. Python lists and tuples). It is meant for dask
    developers and developers of downstream libraries.

    Note that this function does not correspond with NumPy's
    definition of array_like, which includes any object that can be
    coerced into an array (see definition in the NumPy glossary):
    https://numpy.org/doc/stable/glossary.html

    Examples
    --------
    >>> import numpy as np
    >>> is_arraylike(np.ones(5))
    True
    >>> is_arraylike(np.ones(()))
    True
    >>> is_arraylike(5)
    False
    >>> is_arraylike('cat')
    False
    """
    from dask.base import is_dask_collection

    is_duck_array = hasattr(x, "__array_function__") or hasattr(x, "__array_ufunc__")

    return bool(
        hasattr(x, "shape")
        and isinstance(x.shape, tuple)
        and hasattr(x, "dtype")
        and not any(is_dask_collection(n) for n in x.shape)
        # We special case scipy.sparse and cupyx.scipy.sparse arrays as having partial
        # support for them is useful in scenarios where we mostly call `map_partitions`
        # or `map_blocks` with scikit-learn functions on dask arrays and dask dataframes.
        # https://github.com/dask/dask/pull/3738
        and (is_duck_array or "scipy.sparse" in typename(type(x)))
    )


def is_dataframe_like(df) -> bool:
    """Looks like a Pandas DataFrame"""
    if (df.__class__.__module__, df.__class__.__name__) == (
        "pandas.core.frame",
        "DataFrame",
    ):
        # fast exec for most likely input
        return True
    typ = df.__class__
    return (
        all(hasattr(typ, name) for name in ("groupby", "head", "merge", "mean"))
        and all(hasattr(df, name) for name in ("dtypes", "columns"))
        and not any(hasattr(typ, name) for name in ("name", "dtype"))
    )


def is_series_like(s) -> bool:
    """Looks like a Pandas Series"""
    typ = s.__class__
    return (
        all(hasattr(typ, name) for name in ("groupby", "head", "mean"))
        and all(hasattr(s, name) for name in ("dtype", "name"))
        and "index" not in typ.__name__.lower()
    )


def is_index_like(s) -> bool:
    """Looks like a Pandas Index"""
    typ = s.__class__
    return (
        all(hasattr(s, name) for name in ("name", "dtype"))
        and "index" in typ.__name__.lower()
    )


def is_cupy_type(x) -> bool:
    # TODO: avoid explicit reference to CuPy
    return "cupy" in str(type(x))


def natural_sort_key(s: str) -> list[str | int]:
    """
    Sorting `key` function for performing a natural sort on a collection of
    strings

    See https://en.wikipedia.org/wiki/Natural_sort_order

    Parameters
    ----------
    s : str
        A string that is an element of the collection being sorted

    Returns
    -------
    tuple[str or int]
        Tuple of the parts of the input string where each part is either a
        string or an integer

    Examples
    --------
    >>> a = ['f0', 'f1', 'f2', 'f8', 'f9', 'f10', 'f11', 'f19', 'f20', 'f21']
    >>> sorted(a)
    ['f0', 'f1', 'f10', 'f11', 'f19', 'f2', 'f20', 'f21', 'f8', 'f9']
    >>> sorted(a, key=natural_sort_key)
    ['f0', 'f1', 'f2', 'f8', 'f9', 'f10', 'f11', 'f19', 'f20', 'f21']
    """
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", s)]


def parse_bytes(s: float | str) -> int:
    """Parse byte string to numbers

    >>> from dask.utils import parse_bytes
    >>> parse_bytes('100')
    100
    >>> parse_bytes('100 MB')
    100000000
    >>> parse_bytes('100M')
    100000000
    >>> parse_bytes('5kB')
    5000
    >>> parse_bytes('5.4 kB')
    5400
    >>> parse_bytes('1kiB')
    1024
    >>> parse_bytes('1e6')
    1000000
    >>> parse_bytes('1e6 kB')
    1000000000
    >>> parse_bytes('MB')
    1000000
    >>> parse_bytes(123)
    123
    >>> parse_bytes('5 foos')
    Traceback (most recent call last):
        ...
    ValueError: Could not interpret 'foos' as a byte unit
    """
    if isinstance(s, (int, float)):
        return int(s)
    s = s.replace(" ", "")
    if not any(char.isdigit() for char in s):
        s = "1" + s

    for i in range(len(s) - 1, -1, -1):
        if not s[i].isalpha():
            break
    index = i + 1

    prefix = s[:index]
    suffix = s[index:]

    try:
        n = float(prefix)
    except ValueError as e:
        raise ValueError("Could not interpret '%s' as a number" % prefix) from e

    try:
        multiplier = byte_sizes[suffix.lower()]
    except KeyError as e:
        raise ValueError("Could not interpret '%s' as a byte unit" % suffix) from e

    result = n * multiplier
    return int(result)


byte_sizes = {
    "kB": 10**3,
    "MB": 10**6,
    "GB": 10**9,
    "TB": 10**12,
    "PB": 10**15,
    "KiB": 2**10,
    "MiB": 2**20,
    "GiB": 2**30,
    "TiB": 2**40,
    "PiB": 2**50,
    "B": 1,
    "": 1,
}
byte_sizes = {k.lower(): v for k, v in byte_sizes.items()}
byte_sizes.update({k[0]: v for k, v in byte_sizes.items() if k and "i" not in k})
byte_sizes.update({k[:-1]: v for k, v in byte_sizes.items() if k and "i" in k})


def format_time(n: float) -> str:
    """format integers as time

    >>> from dask.utils import format_time
    >>> format_time(1)
    '1.00 s'
    >>> format_time(0.001234)
    '1.23 ms'
    >>> format_time(0.00012345)
    '123.45 us'
    >>> format_time(123.456)
    '123.46 s'
    >>> format_time(1234.567)
    '20m 34s'
    >>> format_time(12345.67)
    '3hr 25m'
    >>> format_time(123456.78)
    '34hr 17m'
    >>> format_time(1234567.89)
    '14d 6hr'
    """
    if n > 24 * 60 * 60 * 2:
        d = int(n / 3600 / 24)
        h = int((n - d * 3600 * 24) / 3600)
        return f"{d}d {h}hr"
    if n > 60 * 60 * 2:
        h = int(n / 3600)
        m = int((n - h * 3600) / 60)
        return f"{h}hr {m}m"
    if n > 60 * 10:
        m = int(n / 60)
        s = int(n - m * 60)
        return f"{m}m {s}s"
    if n >= 1:
        return "%.2f s" % n
    if n >= 1e-3:
        return "%.2f ms" % (n * 1e3)
    return "%.2f us" % (n * 1e6)


def format_time_ago(n: datetime) -> str:
    """Calculate a '3 hours ago' type string from a Python datetime.

    Examples
    --------
    >>> from datetime import datetime, timedelta

    >>> now = datetime.now()
    >>> format_time_ago(now)
    'Just now'

    >>> past = datetime.now() - timedelta(minutes=1)
    >>> format_time_ago(past)
    '1 minute ago'

    >>> past = datetime.now() - timedelta(minutes=2)
    >>> format_time_ago(past)
    '2 minutes ago'

    >>> past = datetime.now() - timedelta(hours=1)
    >>> format_time_ago(past)
    '1 hour ago'

    >>> past = datetime.now() - timedelta(hours=6)
    >>> format_time_ago(past)
    '6 hours ago'

    >>> past = datetime.now() - timedelta(days=1)
    >>> format_time_ago(past)
    '1 day ago'

    >>> past = datetime.now() - timedelta(days=5)
    >>> format_time_ago(past)
    '5 days ago'

    >>> past = datetime.now() - timedelta(days=8)
    >>> format_time_ago(past)
    '1 week ago'

    >>> past = datetime.now() - timedelta(days=16)
    >>> format_time_ago(past)
    '2 weeks ago'

    >>> past = datetime.now() - timedelta(days=190)
    >>> format_time_ago(past)
    '6 months ago'

    >>> past = datetime.now() - timedelta(days=800)
    >>> format_time_ago(past)
    '2 years ago'

    """
    units = {
        "years": lambda diff: diff.days / 365,
        "months": lambda diff: diff.days / 30.436875,  # Average days per month
        "weeks": lambda diff: diff.days / 7,
        "days": lambda diff: diff.days,
        "hours": lambda diff: diff.seconds / 3600,
        "minutes": lambda diff: diff.seconds % 3600 / 60,
    }
    diff = datetime.now() - n
    for unit, func in units.items():
        dur = int(func(diff))
        if dur > 0:
            if dur == 1:  # De-pluralize
                unit = unit[:-1]
            return f"{dur} {unit} ago"
    return "Just now"


def format_bytes(n: int) -> str:
    """Format bytes as text

    >>> from dask.utils import format_bytes
    >>> format_bytes(1)
    '1 B'
    >>> format_bytes(1234)
    '1.21 kiB'
    >>> format_bytes(12345678)
    '11.77 MiB'
    >>> format_bytes(1234567890)
    '1.15 GiB'
    >>> format_bytes(1234567890000)
    '1.12 TiB'
    >>> format_bytes(1234567890000000)
    '1.10 PiB'

    For all values < 2**60, the output is always <= 10 characters.
    """
    for prefix, k in (
        ("Pi", 2**50),
        ("Ti", 2**40),
        ("Gi", 2**30),
        ("Mi", 2**20),
        ("ki", 2**10),
    ):
        if n >= k * 0.9:
            return f"{n / k:.2f} {prefix}B"
    return f"{n} B"


timedelta_sizes = {
    "s": 1,
    "ms": 1e-3,
    "us": 1e-6,
    "ns": 1e-9,
    "m": 60,
    "h": 3600,
    "d": 3600 * 24,
    "w": 7 * 3600 * 24,
}

tds2 = {
    "second": 1,
    "minute": 60,
    "hour": 60 * 60,
    "day": 60 * 60 * 24,
    "week": 7 * 60 * 60 * 24,
    "millisecond": 1e-3,
    "microsecond": 1e-6,
    "nanosecond": 1e-9,
}
tds2.update({k + "s": v for k, v in tds2.items()})
timedelta_sizes.update(tds2)
timedelta_sizes.update({k.upper(): v for k, v in timedelta_sizes.items()})


@overload
def parse_timedelta(s: None, default: str | Literal[False] = "seconds") -> None: ...


@overload
def parse_timedelta(
    s: str | float | timedelta, default: str | Literal[False] = "seconds"
) -> float: ...


def parse_timedelta(s, default="seconds"):
    """Parse timedelta string to number of seconds

    Parameters
    ----------
    s : str, float, timedelta, or None
    default: str or False, optional
        Unit of measure if s  does not specify one. Defaults to seconds.
        Set to False to require s to explicitly specify its own unit.

    Examples
    --------
    >>> from datetime import timedelta
    >>> from dask.utils import parse_timedelta
    >>> parse_timedelta('3s')
    3
    >>> parse_timedelta('3.5 seconds')
    3.5
    >>> parse_timedelta('300ms')
    0.3
    >>> parse_timedelta(timedelta(seconds=3))  # also supports timedeltas
    3
    """
    if s is None:
        return None
    if isinstance(s, timedelta):
        s = s.total_seconds()
        return int(s) if int(s) == s else s
    if isinstance(s, Number):
        s = str(s)
    s = s.replace(" ", "")
    if not s[0].isdigit():
        s = "1" + s

    for i in range(len(s) - 1, -1, -1):
        if not s[i].isalpha():
            break
    index = i + 1

    prefix = s[:index]
    suffix = s[index:] or default
    if suffix is False:
        raise ValueError(f"Missing time unit: {s}")
    if not isinstance(suffix, str):
        raise TypeError(f"default must be str or False, got {default!r}")

    n = float(prefix)

    try:
        multiplier = timedelta_sizes[suffix.lower()]
    except KeyError:
        valid_units = ", ".join(timedelta_sizes.keys())
        raise KeyError(
            f"Invalid time unit: {suffix}. Valid units are: {valid_units}"
        ) from None

    result = n * multiplier
    if int(result) == result:
        result = int(result)
    return result


def has_keyword(func, keyword):
    try:
        return keyword in inspect.signature(func).parameters
    except Exception:
        return False


def ndimlist(seq):
    if not isinstance(seq, (list, tuple)):
        return 0
    elif not seq:
        return 1
    else:
        return 1 + ndimlist(seq[0])


def iter_chunks(sizes, max_size):
    """Split sizes into chunks of total max_size each

    Parameters
    ----------
    sizes : iterable of numbers
        The sizes to be chunked
    max_size : number
        Maximum total size per chunk.
        It must be greater or equal than each size in sizes
    """
    chunk, chunk_sum = [], 0
    iter_sizes = iter(sizes)
    size = next(iter_sizes, None)
    while size is not None:
        assert size <= max_size
        if chunk_sum + size <= max_size:
            chunk.append(size)
            chunk_sum += size
            size = next(iter_sizes, None)
        else:
            assert chunk
            yield chunk
            chunk, chunk_sum = [], 0
    if chunk:
        yield chunk


hex_pattern = re.compile("[a-f]+")


@functools.lru_cache(100000)
def key_split(s):
    """
    >>> key_split('x')
    'x'
    >>> key_split('x-1')
    'x'
    >>> key_split('x-1-2-3')
    'x'
    >>> key_split(('x-2', 1))
    'x'
    >>> key_split("('x-2', 1)")
    'x'
    >>> key_split("('x', 1)")
    'x'
    >>> key_split('hello-world-1')
    'hello-world'
    >>> key_split(b'hello-world-1')
    'hello-world'
    >>> key_split('ae05086432ca935f6eba409a8ecd4896')
    'data'
    >>> key_split('<module.submodule.myclass object at 0xdaf372')
    'myclass'
    >>> key_split(None)
    'Other'
    >>> key_split('x-abcdefab')  # ignores hex
    'x'
    >>> key_split('_(x)')  # strips unpleasant characters
    'x'
    """
    # If we convert the key, recurse to utilize LRU cache better
    if type(s) is bytes:
        return key_split(s.decode())
    if type(s) is tuple:
        return key_split(s[0])
    try:
        words = s.split("-")
        if not words[0][0].isalpha():
            result = words[0].split(",")[0].strip("_'()\"")
        else:
            result = words[0]
        for word in words[1:]:
            if word.isalpha() and not (
                len(word) == 8 and hex_pattern.match(word) is not None
            ):
                result += "-" + word
            else:
                break
        if len(result) == 32 and re.match(r"[a-f0-9]{32}", result):
            return "data"
        else:
            if result[0] == "<":
                result = result.strip("<>").split()[0].split(".")[-1]
            return sys.intern(result)
    except Exception:
        return "Other"


def stringify(obj, exclusive: Iterable | None = None):
    """Convert an object to a string

    If ``exclusive`` is specified, search through `obj` and convert
    values that are in ``exclusive``.

    Note that when searching through dictionaries, only values are
    converted, not the keys.

    Parameters
    ----------
    obj : Any
        Object (or values within) to convert to string
    exclusive: Iterable, optional
        Set of values to search for when converting values to strings

    Returns
    -------
    result : type(obj)
        Stringified copy of ``obj`` or ``obj`` itself if it is already a
        string or bytes.

    Examples
    --------
    >>> stringify(b'x')
    b'x'
    >>> stringify('x')
    'x'
    >>> stringify({('a',0):('a',0), ('a',1): ('a',1)})
    "{('a', 0): ('a', 0), ('a', 1): ('a', 1)}"
    >>> stringify({('a',0):('a',0), ('a',1): ('a',1)}, exclusive={('a',0)})
    {('a', 0): "('a', 0)", ('a', 1): ('a', 1)}
    """

    typ = type(obj)
    if typ is str or typ is bytes:
        return obj
    elif exclusive is None:
        return str(obj)

    if typ is list:
        return [stringify(v, exclusive) for v in obj]
    if typ is dict:
        return {k: stringify(v, exclusive) for k, v in obj.items()}
    try:
        if obj in exclusive:
            return stringify(obj)
    except TypeError:  # `obj` not hashable
        pass
    if typ is tuple:  # If the tuple itself isn't a key, check its elements
        return tuple(stringify(v, exclusive) for v in obj)
    return obj


class cached_property(functools.cached_property):
    """Read only version of functools.cached_property."""

    def __set__(self, instance, val):
        """Raise an error when attempting to set a cached property."""
        raise AttributeError("Can't set attribute")


class _HashIdWrapper:
    """Hash and compare a wrapped object by identity instead of value"""

    def __init__(self, wrapped):
        self.wrapped = wrapped

    def __eq__(self, other):
        if not isinstance(other, _HashIdWrapper):
            return NotImplemented
        return self.wrapped is other.wrapped

    def __ne__(self, other):
        if not isinstance(other, _HashIdWrapper):
            return NotImplemented
        return self.wrapped is not other.wrapped

    def __hash__(self):
        return id(self.wrapped)


@functools.lru_cache
def _cumsum(seq, initial_zero):
    if isinstance(seq, _HashIdWrapper):
        seq = seq.wrapped
    if initial_zero:
        return tuple(toolz.accumulate(add, seq, 0))
    else:
        return tuple(toolz.accumulate(add, seq))


@functools.lru_cache
def _max(seq):
    if isinstance(seq, _HashIdWrapper):
        seq = seq.wrapped
    return max(seq)


def cached_max(seq):
    """Compute max with caching.

    Caching is by the identity of `seq` rather than the value. It is thus
    important that `seq` is a tuple of immutable objects, and this function
    is intended for use where `seq` is a value that will persist (generally
    block sizes).

    Parameters
    ----------
    seq : tuple
        Values to reduce

    Returns
    -------
    tuple
    """
    assert isinstance(seq, tuple)
    # Look up by identity first, to avoid a linear-time __hash__
    # if we've seen this tuple object before.
    result = _max(_HashIdWrapper(seq))
    return result


def cached_cumsum(seq, initial_zero=False):
    """Compute :meth:`toolz.accumulate` with caching.

    Caching is by the identify of `seq` rather than the value. It is thus
    important that `seq` is a tuple of immutable objects, and this function
    is intended for use where `seq` is a value that will persist (generally
    block sizes).

    Parameters
    ----------
    seq : tuple
        Values to cumulatively sum.
    initial_zero : bool, optional
        If true, the return value is prefixed with a zero.

    Returns
    -------
    tuple
    """
    if isinstance(seq, tuple):
        # Look up by identity first, to avoid a linear-time __hash__
        # if we've seen this tuple object before.
        result = _cumsum(_HashIdWrapper(seq), initial_zero)
    else:
        # Construct a temporary tuple, and look up by value.
        result = _cumsum(tuple(seq), initial_zero)
    return result


def show_versions() -> None:
    """Provide version information for bug reports."""

    from json import dumps
    from platform import uname
    from sys import stdout, version_info

    from dask._compatibility import importlib_metadata

    try:
        from distributed import __version__ as distributed_version
    except ImportError:
        distributed_version = None

    from dask import __version__ as dask_version

    deps = [
        "numpy",
        "pandas",
        "cloudpickle",
        "fsspec",
        "bokeh",
        "pyarrow",
        "zarr",
    ]

    result: dict[str, str | None] = {
        # note: only major, minor, micro are extracted
        "Python": ".".join([str(i) for i in version_info[:3]]),
        "Platform": uname().system,
        "dask": dask_version,
        "distributed": distributed_version,
    }

    for modname in deps:
        try:
            result[modname] = importlib_metadata.version(modname)
        except importlib_metadata.PackageNotFoundError:
            result[modname] = None

    stdout.writelines(dumps(result, indent=2))

    return


def maybe_pluralize(count, noun, plural_form=None):
    """Pluralize a count-noun string pattern when necessary"""
    if count == 1:
        return f"{count} {noun}"
    else:
        return f"{count} {plural_form or noun + 's'}"


def is_namedtuple_instance(obj: Any) -> bool:
    """Returns True if obj is an instance of a namedtuple.

    Note: This function checks for the existence of the methods and
    attributes that make up the namedtuple API, so it will return True
    IFF obj's type implements that API.
    """
    return (
        isinstance(obj, tuple)
        and hasattr(obj, "_make")
        and hasattr(obj, "_asdict")
        and hasattr(obj, "_replace")
        and hasattr(obj, "_fields")
        and hasattr(obj, "_field_defaults")
    )


def get_default_shuffle_method() -> str:
    if d := config.get("dataframe.shuffle.method", None):
        return d
    try:
        from distributed import default_client

        default_client()
    except (ImportError, ValueError):
        return "disk"

    try:
        from distributed.shuffle import check_minimal_arrow_version

        check_minimal_arrow_version()
    except ModuleNotFoundError:
        return "tasks"
    return "p2p"


def get_meta_library(like):
    if hasattr(like, "_meta"):
        like = like._meta

    return import_module(typename(like).partition(".")[0])


class shorten_traceback:
    """Context manager that removes irrelevant stack elements from traceback.

    * omits frames from modules that match `admin.traceback.shorten`
    * always keeps the first and last frame.
    """

    __slots__ = ()

    def __enter__(self) -> None:
        pass

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: types.TracebackType | None,
    ) -> None:
        if exc_val and exc_tb:
            exc_val.__traceback__ = self.shorten(exc_tb)

    @staticmethod
    def shorten(exc_tb: types.TracebackType) -> types.TracebackType:
        paths = config.get("admin.traceback.shorten")
        if not paths:
            return exc_tb

        exp = re.compile(".*(" + "|".join(paths) + ")")
        curr: types.TracebackType | None = exc_tb
        prev: types.TracebackType | None = None

        while curr:
            if prev is None:
                prev = curr  # first frame
            elif not curr.tb_next:
                # always keep last frame
                prev.tb_next = curr
                prev = prev.tb_next
            elif not exp.match(curr.tb_frame.f_code.co_filename):
                # keep if module is not listed in config
                prev.tb_next = curr
                prev = curr
            curr = curr.tb_next

        # Uncomment to remove the first frame, which is something you don't want to keep
        # if it matches the regexes. Requires Python >=3.11.
        # if exc_tb.tb_next and exp.match(exc_tb.tb_frame.f_code.co_filename):
        #     return exc_tb.tb_next

        return exc_tb


def unzip(ls, nout):
    """Unzip a list of lists into ``nout`` outputs."""
    out = list(zip(*ls))
    if not out:
        out = [()] * nout
    return out


class disable_gc(ContextDecorator):
    """Context manager to disable garbage collection."""

    def __init__(self, collect=False):
        self.collect = collect
        self._gc_enabled = gc.isenabled()

    def __enter__(self):
        gc.disable()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._gc_enabled:
            gc.enable()
        return False
