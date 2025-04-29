from __future__ import annotations

import operator
import types
import uuid
import warnings
from collections.abc import Sequence
from dataclasses import fields, is_dataclass, replace
from functools import partial

from tlz import concat, curry, merge, unique

from dask import config
from dask.base import (
    DaskMethodsMixin,
    dont_optimize,
    is_dask_collection,
    named_schedulers,
    replace_name_in_key,
)
from dask.base import tokenize as _tokenize
from dask.context import globalmethod
from dask.core import flatten, quote
from dask.highlevelgraph import HighLevelGraph
from dask.optimization import fuse
from dask.typing import Graph, NestedKeys
from dask.utils import (
    OperatorMethodMixin,
    apply,
    ensure_dict,
    funcname,
    is_namedtuple_instance,
    methodcaller,
    unzip,
)

__all__ = ["Delayed", "delayed"]


DEFAULT_GET = named_schedulers.get("threads", named_schedulers["sync"])


def finalize(collection):
    assert is_dask_collection(collection)

    name = "finalize-" + tokenize(collection)
    keys = collection.__dask_keys__()
    finalize, args = collection.__dask_postcompute__()
    layer = {name: (finalize, keys) + args}
    graph = HighLevelGraph.from_collections(name, layer, dependencies=[collection])
    return Delayed(name, graph)


def unpack_collections(expr):
    """Normalize a python object and merge all sub-graphs.

    - Replace ``Delayed`` with their keys
    - Convert literals to things the schedulers can handle
    - Extract dask graphs from all enclosed values

    Parameters
    ----------
    expr : object
        The object to be normalized. This function knows how to handle
        dask collections, as well as most builtin python types.

    Returns
    -------
    task : normalized task to be run
    collections : a tuple of collections

    Examples
    --------
    >>> import dask
    >>> a = delayed(1, 'a')
    >>> b = delayed(2, 'b')
    >>> task, collections = unpack_collections([a, b, 3])
    >>> task
    ['a', 'b', 3]
    >>> collections
    (Delayed('a'), Delayed('b'))

    >>> task, collections = unpack_collections({a: 1, b: 2})
    >>> task
    (<class 'dict'>, [['a', 1], ['b', 2]])
    >>> collections
    (Delayed('a'), Delayed('b'))
    """
    # FIXME: See blockwise._unpack_collections for a TaskSpec version
    if isinstance(expr, Delayed):
        return expr._key, (expr,)

    if is_dask_collection(expr):
        if hasattr(expr, "optimize"):
            # Optimize dask-expr collections
            expr = expr.optimize()

        finalized = finalize(expr)
        return finalized._key, (finalized,)

    if type(expr) is type(iter(list())):
        expr = list(expr)
    elif type(expr) is type(iter(tuple())):
        expr = tuple(expr)
    elif type(expr) is type(iter(set())):
        expr = set(expr)

    typ = type(expr)

    if typ in (list, tuple, set):
        args, collections = unzip((unpack_collections(e) for e in expr), 2)
        args = list(args)
        collections = tuple(unique(concat(collections), key=id))
        # Ensure output type matches input type
        if typ is not list:
            args = (typ, args)
        return args, collections

    if typ is dict:
        args, collections = unpack_collections([[k, v] for k, v in expr.items()])
        return (dict, args), collections

    if typ is slice:
        args, collections = unpack_collections([expr.start, expr.stop, expr.step])
        return (slice, *args), collections

    if is_dataclass(expr):
        args, collections = unpack_collections(
            [
                [f.name, getattr(expr, f.name)]
                for f in fields(expr)
                if hasattr(expr, f.name)  # if init=False, field might not exist
            ]
        )
        if not collections:
            return expr, ()
        try:
            _fields = {
                f.name: getattr(expr, f.name)
                for f in fields(expr)
                if hasattr(expr, f.name)
            }
            replace(expr, **_fields)
        except (TypeError, ValueError) as e:
            if isinstance(e, ValueError) or "is declared with init=False" in str(e):
                raise ValueError(
                    f"Failed to unpack {typ} instance. "
                    "Note that using fields with `init=False` are not supported."
                ) from e
            else:
                raise TypeError(
                    f"Failed to unpack {typ} instance. "
                    "Note that using a custom __init__ is not supported."
                ) from e
        return (apply, typ, (), (dict, args)), collections

    if is_namedtuple_instance(expr):
        args, collections = unpack_collections([v for v in expr])
        return (typ, *args), collections

    return expr, ()


def to_task_dask(expr):
    """Normalize a python object and merge all sub-graphs.

    - Replace ``Delayed`` with their keys
    - Convert literals to things the schedulers can handle
    - Extract dask graphs from all enclosed values

    Parameters
    ----------
    expr : object
        The object to be normalized. This function knows how to handle
        ``Delayed``s, as well as most builtin python types.

    Returns
    -------
    task : normalized task to be run
    dask : a merged dask graph that forms the dag for this task

    Examples
    --------
    >>> import dask
    >>> a = delayed(1, 'a')
    >>> b = delayed(2, 'b')
    >>> task, dask = to_task_dask([a, b, 3])  # doctest: +SKIP
    >>> task  # doctest: +SKIP
    ['a', 'b', 3]
    >>> dict(dask)  # doctest: +SKIP
    {'a': 1, 'b': 2}

    >>> task, dasks = to_task_dask({a: 1, b: 2})  # doctest: +SKIP
    >>> task  # doctest: +SKIP
    (dict, [['a', 1], ['b', 2]])
    >>> dict(dask)  # doctest: +SKIP
    {'a': 1, 'b': 2}
    """
    warnings.warn(
        "The dask.delayed.to_dask_dask function has been "
        "Deprecated in favor of unpack_collections",
        stacklevel=2,
    )

    if isinstance(expr, Delayed):
        return expr.key, expr.dask

    if is_dask_collection(expr):
        name = "finalize-" + tokenize(expr, pure=True)
        keys = expr.__dask_keys__()
        opt = getattr(expr, "__dask_optimize__", dont_optimize)
        finalize, args = expr.__dask_postcompute__()
        dsk = {name: (finalize, keys) + args}
        dsk.update(opt(expr.__dask_graph__(), keys))
        return name, dsk

    if type(expr) is type(iter(list())):
        expr = list(expr)
    elif type(expr) is type(iter(tuple())):
        expr = tuple(expr)
    elif type(expr) is type(iter(set())):
        expr = set(expr)
    typ = type(expr)

    if typ in (list, tuple, set):
        args, dasks = unzip((to_task_dask(e) for e in expr), 2)
        args = list(args)
        dsk = merge(dasks)
        # Ensure output type matches input type
        return (args, dsk) if typ is list else ((typ, args), dsk)

    if typ is dict:
        args, dsk = to_task_dask([[k, v] for k, v in expr.items()])
        return (dict, args), dsk

    if is_dataclass(expr):
        args, dsk = to_task_dask(
            [
                [f.name, getattr(expr, f.name)]
                for f in fields(expr)
                if hasattr(expr, f.name)  # if init=False, field might not exist
            ]
        )

        return (apply, typ, (), (dict, args)), dsk

    if is_namedtuple_instance(expr):
        args, dsk = to_task_dask([v for v in expr])
        return (typ, *args), dsk

    if typ is slice:
        args, dsk = to_task_dask([expr.start, expr.stop, expr.step])
        return (slice,) + tuple(args), dsk

    return expr, {}


def tokenize(*args, pure=None, **kwargs):
    """Mapping function from task -> consistent name.

    Parameters
    ----------
    args : object
        Python objects that summarize the task.
    pure : boolean, optional
        If True, a consistent hash function is tried on the input. If this
        fails, then a unique identifier is used. If False (default), then a
        unique identifier is always used.
    """
    if pure is None:
        pure = config.get("delayed_pure", False)

    if pure:
        return _tokenize(*args, **kwargs)
    else:
        return str(uuid.uuid4())


@curry
def delayed(obj, name=None, pure=None, nout=None, traverse=True):
    """Wraps a function or object to produce a ``Delayed``.

    ``Delayed`` objects act as proxies for the object they wrap, but all
    operations on them are done lazily by building up a dask graph internally.

    Parameters
    ----------
    obj : object
        The function or object to wrap
    name : Dask key, optional
        The key to use in the underlying graph for the wrapped object. Defaults
        to hashing content. Note that this only affects the name of the object
        wrapped by this call to delayed, and *not* the output of delayed
        function calls - for that use ``dask_key_name=`` as described below.

        .. note::

           Because this ``name`` is used as the key in task graphs, you should
           ensure that it uniquely identifies ``obj``. If you'd like to provide
           a descriptive name that is still unique, combine the descriptive name
           with :func:`dask.base.tokenize` of the ``array_like``. See
           :ref:`graphs` for more.

    pure : bool, optional
        Indicates whether calling the resulting ``Delayed`` object is a pure
        operation. If True, arguments to the call are hashed to produce
        deterministic keys. If not provided, the default is to check the global
        ``delayed_pure`` setting, and fallback to ``False`` if unset.
    nout : int, optional
        The number of outputs returned from calling the resulting ``Delayed``
        object. If provided, the ``Delayed`` output of the call can be iterated
        into ``nout`` objects, allowing for unpacking of results. By default
        iteration over ``Delayed`` objects will error. Note, that ``nout=1``
        expects ``obj`` to return a tuple of length 1, and consequently for
        ``nout=0``, ``obj`` should return an empty tuple.
    traverse : bool, optional
        By default dask traverses builtin python collections looking for dask
        objects passed to ``delayed``. For large collections this can be
        expensive. If ``obj`` doesn't contain any dask objects, set
        ``traverse=False`` to avoid doing this traversal.

    Examples
    --------
    Apply to functions to delay execution:

    >>> from dask import delayed
    >>> def inc(x):
    ...     return x + 1

    >>> inc(10)
    11

    >>> x = delayed(inc, pure=True)(10)
    >>> type(x) == Delayed
    True
    >>> x.compute()
    11

    Can be used as a decorator:

    >>> @delayed(pure=True)
    ... def add(a, b):
    ...     return a + b
    >>> add(1, 2).compute()
    3

    ``delayed`` also accepts an optional keyword ``pure``. If False, then
    subsequent calls will always produce a different ``Delayed``. This is
    useful for non-pure functions (such as ``time`` or ``random``).

    >>> from random import random
    >>> out1 = delayed(random, pure=False)()
    >>> out2 = delayed(random, pure=False)()
    >>> out1.key == out2.key
    False

    If you know a function is pure (output only depends on the input, with no
    global state), then you can set ``pure=True``. This will attempt to apply a
    consistent name to the output, but will fallback on the same behavior of
    ``pure=False`` if this fails.

    >>> @delayed(pure=True)
    ... def add(a, b):
    ...     return a + b
    >>> out1 = add(1, 2)
    >>> out2 = add(1, 2)
    >>> out1.key == out2.key
    True

    Instead of setting ``pure`` as a property of the callable, you can also set
    it contextually using the ``delayed_pure`` setting. Note that this
    influences the *call* and not the *creation* of the callable:

    >>> @delayed
    ... def mul(a, b):
    ...     return a * b
    >>> import dask
    >>> with dask.config.set(delayed_pure=True):
    ...     print(mul(1, 2).key == mul(1, 2).key)
    True
    >>> with dask.config.set(delayed_pure=False):
    ...     print(mul(1, 2).key == mul(1, 2).key)
    False

    The key name of the result of calling a delayed object is determined by
    hashing the arguments by default. To explicitly set the name, you can use
    the ``dask_key_name`` keyword when calling the function:

    >>> add(1, 2)   # doctest: +SKIP
    Delayed('add-3dce7c56edd1ac2614add714086e950f')
    >>> add(1, 2, dask_key_name='three')
    Delayed('three')

    Note that objects with the same key name are assumed to have the same
    result. If you set the names explicitly you should make sure your key names
    are different for different results.

    >>> add(1, 2, dask_key_name='three')
    Delayed('three')
    >>> add(2, 1, dask_key_name='three')
    Delayed('three')
    >>> add(2, 2, dask_key_name='four')
    Delayed('four')

    ``delayed`` can also be applied to objects to make operations on them lazy:

    >>> a = delayed([1, 2, 3])
    >>> isinstance(a, Delayed)
    True
    >>> a.compute()
    [1, 2, 3]

    The key name of a delayed object is hashed by default if ``pure=True`` or
    is generated randomly if ``pure=False`` (default).  To explicitly set the
    name, you can use the ``name`` keyword. To ensure that the key is unique
    you should include the tokenized value as well, or otherwise ensure that
    it's unique:

    >>> from dask.base import tokenize
    >>> data = [1, 2, 3]
    >>> a = delayed(data, name='mylist-' + tokenize(data))
    >>> a  # doctest: +SKIP
    Delayed('mylist-55af65871cb378a4fa6de1660c3e8fb7')

    Delayed results act as a proxy to the underlying object. Many operators
    are supported:

    >>> (a + [1, 2]).compute()
    [1, 2, 3, 1, 2]
    >>> a[1].compute()
    2

    Method and attribute access also works:

    >>> a.count(2).compute()
    1

    Note that if a method doesn't exist, no error will be thrown until runtime:

    >>> res = a.not_a_real_method() # doctest: +SKIP
    >>> res.compute()  # doctest: +SKIP
    AttributeError("'list' object has no attribute 'not_a_real_method'")

    "Magic" methods (e.g. operators and attribute access) are assumed to be
    pure, meaning that subsequent calls must return the same results. This
    behavior is not overridable through the ``delayed`` call, but can be
    modified using other ways as described below.

    To invoke an impure attribute or operator, you'd need to use it in a
    delayed function with ``pure=False``:

    >>> class Incrementer:
    ...     def __init__(self):
    ...         self._n = 0
    ...     @property
    ...     def n(self):
    ...         self._n += 1
    ...         return self._n
    ...
    >>> x = delayed(Incrementer())
    >>> x.n.key == x.n.key
    True
    >>> get_n = delayed(lambda x: x.n, pure=False)
    >>> get_n(x).key == get_n(x).key
    False

    In contrast, methods are assumed to be impure by default, meaning that
    subsequent calls may return different results. To assume purity, set
    ``pure=True``. This allows sharing of any intermediate values.

    >>> a.count(2, pure=True).key == a.count(2, pure=True).key
    True

    As with function calls, method calls also respect the global
    ``delayed_pure`` setting and support the ``dask_key_name`` keyword:

    >>> a.count(2, dask_key_name="count_2")
    Delayed('count_2')
    >>> import dask
    >>> with dask.config.set(delayed_pure=True):
    ...     print(a.count(2).key == a.count(2).key)
    True
    """
    if isinstance(obj, Delayed):
        return obj

    if is_dask_collection(obj) or traverse:
        task, collections = unpack_collections(obj)
    else:
        task = quote(obj)
        collections = set()

    if not (nout is None or (type(nout) is int and nout >= 0)):
        raise ValueError("nout must be None or a non-negative integer, got %s" % nout)
    if task is obj:
        if not name:
            try:
                prefix = obj.__name__
            except AttributeError:
                prefix = type(obj).__name__
            token = tokenize(obj, nout, pure=pure)
            name = f"{prefix}-{token}"
        return DelayedLeaf(obj, name, pure=pure, nout=nout)
    else:
        if not name:
            name = f"{type(obj).__name__}-{tokenize(task, pure=pure)}"
        layer = {name: task}
        graph = HighLevelGraph.from_collections(name, layer, dependencies=collections)
        return Delayed(name, graph, nout)


def _swap(method, self, other):
    return method(other, self)


def right(method):
    """Wrapper to create 'right' version of operator given left version"""
    return partial(_swap, method)


def _fuse_delayed(dsk, keys, **kwargs):
    dependencies = dsk.get_all_dependencies()
    dsk = ensure_dict(dsk)

    dsk, _ = fuse(
        dsk,
        keys,
        dependencies=dependencies,
    )
    return dsk


def optimize(dsk, keys, **kwargs):
    if not isinstance(keys, (list, set)):
        keys = [keys]

    if config.get("optimization.fuse.delayed"):
        dsk = _fuse_delayed(dsk, keys, **kwargs)

    if not isinstance(dsk, HighLevelGraph):
        dsk = HighLevelGraph.from_collections(id(dsk), dsk, dependencies=())
    dsk = dsk.cull(set(flatten(keys)))
    return dsk


class Delayed(DaskMethodsMixin, OperatorMethodMixin):
    """Represents a value to be computed by dask.

    Equivalent to the output from a single key in a dask graph.
    """

    __slots__ = ("_key", "_dask", "_length", "_layer")

    def __init__(self, key, dsk, length=None, layer=None):
        self._key = key
        self._dask = dsk
        self._length = length

        # NOTE: Layer is used by `to_delayed` in other collections, but not in normal Delayed use
        self._layer = layer or key
        if isinstance(dsk, HighLevelGraph) and self._layer not in dsk.layers:
            raise ValueError(
                f"Layer {self._layer} not in the HighLevelGraph's layers: {list(dsk.layers)}"
            )

    @property
    def key(self):
        return self._key

    @property
    def dask(self):
        return self._dask

    def __dask_graph__(self) -> Graph:
        return self.dask

    def __dask_keys__(self) -> NestedKeys:
        return [self.key]

    def __dask_layers__(self) -> Sequence[str]:
        return (self._layer,)

    def __dask_tokenize__(self):
        return self.key

    __dask_scheduler__ = staticmethod(DEFAULT_GET)
    __dask_optimize__ = globalmethod(optimize, key="delayed_optimize")

    def __dask_postcompute__(self):
        return single_key, ()

    def __dask_postpersist__(self):
        return self._rebuild, ()

    def _rebuild(self, dsk, *, rename=None):
        key = replace_name_in_key(self.key, rename) if rename else self.key
        if isinstance(dsk, HighLevelGraph) and len(dsk.layers) == 1:
            # FIXME Delayed is currently the only collection type that supports both high- and low-level graphs.
            # The HLG output of `optimize` will have a layer name that doesn't match `key`.
            # Remove this when Delayed is HLG-only (because `optimize` will only be passed HLGs, so it won't have
            # to generate random layer names).
            layer = next(iter(dsk.layers))
        else:
            layer = None
        return Delayed(key, dsk, self._length, layer=layer)

    def __repr__(self):
        return f"Delayed({repr(self.key)})"

    def __hash__(self):
        return hash(self.key)

    def __dir__(self):
        return dir(type(self))

    def __getattr__(self, attr):
        if attr.startswith("_"):
            raise AttributeError(f"Attribute {attr} not found")

        if attr == "visualise":
            # added to warn users in case of spelling error
            # for more details: https://github.com/dask/dask/issues/5721
            warnings.warn(
                "dask.delayed objects have no `visualise` method. "
                "Perhaps you meant `visualize`?"
            )

        return DelayedAttr(self, attr)

    def __setattr__(self, attr, val):
        try:
            object.__setattr__(self, attr, val)
        except AttributeError:
            # attr is neither in type(self).__slots__ nor in the __slots__ of any of its
            # parent classes, and all the parent classes define __slots__ too.
            # This last bit needs to be unit tested: if any of the parent classes omit
            # the __slots__ declaration, self will gain a __dict__ and this branch will
            # become unreachable.
            raise TypeError("Delayed objects are immutable")

    def __setitem__(self, index, val):
        raise TypeError("Delayed objects are immutable")

    def __iter__(self):
        if self._length is None:
            raise TypeError("Delayed objects of unspecified length are not iterable")
        for i in range(self._length):
            yield self[i]

    def __len__(self):
        if self._length is None:
            raise TypeError("Delayed objects of unspecified length have no len()")
        return self._length

    def __call__(self, *args, pure=None, dask_key_name=None, **kwargs):
        func = delayed(apply, pure=pure)
        if dask_key_name is not None:
            return func(self, args, kwargs, dask_key_name=dask_key_name)
        return func(self, args, kwargs)

    def __bool__(self):
        raise TypeError("Truth of Delayed objects is not supported")

    __nonzero__ = __bool__

    def __get__(self, instance, cls):
        if instance is None:
            return self
        return types.MethodType(self, instance)

    @classmethod
    def _get_binary_operator(cls, op, inv=False):
        method = delayed(right(op) if inv else op, pure=True)
        return lambda *args, **kwargs: method(*args, **kwargs)

    _get_unary_operator = _get_binary_operator


def call_function(func, func_token, args, kwargs, pure=None, nout=None):
    dask_key_name = kwargs.pop("dask_key_name", None)
    pure = kwargs.pop("pure", pure)

    if dask_key_name is None:
        name = "{}-{}".format(
            funcname(func),
            tokenize(func_token, *args, pure=pure, **kwargs),
        )
    else:
        name = dask_key_name

    args2, collections = unzip(map(unpack_collections, args), 2)
    collections = list(concat(collections))

    if kwargs:
        dask_kwargs, collections2 = unpack_collections(kwargs)
        collections.extend(collections2)
        task = (apply, func, list(args2), dask_kwargs)
    else:
        task = (func,) + args2

    graph = HighLevelGraph.from_collections(
        name, {name: task}, dependencies=collections
    )
    nout = nout if nout is not None else None
    return Delayed(name, graph, length=nout)


class DelayedLeaf(Delayed):
    __slots__ = ("_obj", "_pure", "_nout")

    def __init__(self, obj, key, pure=None, nout=None):
        super().__init__(key, None)
        self._obj = obj
        self._pure = pure
        self._nout = nout

    @property
    def dask(self):
        return HighLevelGraph.from_collections(
            self._key, {self._key: self._obj}, dependencies=()
        )

    def __call__(self, *args, **kwargs):
        return call_function(
            self._obj, self._key, args, kwargs, pure=self._pure, nout=self._nout
        )

    @property
    def __name__(self):
        return self._obj.__name__

    @property
    def __doc__(self):
        return self._obj.__doc__

    @property
    def __wrapped__(self):
        return self._obj


class DelayedAttr(Delayed):
    __slots__ = ("_obj", "_attr")

    def __init__(self, obj, attr):
        key = "getattr-%s" % tokenize(obj, attr, pure=True)
        super().__init__(key, None)
        self._obj = obj
        self._attr = attr

    def __getattr__(self, attr):
        # Calling np.dtype(dask.delayed(...)) used to result in a segfault, as
        # numpy recursively tries to get `dtype` from the object. This is
        # likely a bug in numpy. For now, we can do a dumb for if
        # `x.dtype().dtype()` is called (which shouldn't ever show up in real
        # code). See https://github.com/dask/dask/pull/4374#issuecomment-454381465
        if attr == "dtype" and self._attr == "dtype":
            raise AttributeError("Attribute dtype not found")
        return super().__getattr__(attr)

    @property
    def dask(self):
        layer = {self._key: (getattr, self._obj._key, self._attr)}
        return HighLevelGraph.from_collections(
            self._key, layer, dependencies=[self._obj]
        )

    def __call__(self, *args, **kwargs):
        return call_function(
            methodcaller(self._attr), self._attr, (self._obj,) + args, kwargs
        )


for op in [
    operator.abs,
    operator.neg,
    operator.pos,
    operator.invert,
    operator.add,
    operator.sub,
    operator.mul,
    operator.floordiv,
    operator.truediv,
    operator.mod,
    operator.pow,
    operator.and_,
    operator.or_,
    operator.xor,
    operator.lshift,
    operator.rshift,
    operator.eq,
    operator.ge,
    operator.gt,
    operator.ne,
    operator.le,
    operator.lt,
    operator.getitem,
]:
    Delayed._bind_operator(op)


try:
    Delayed._bind_operator(operator.matmul)
except AttributeError:
    pass


def single_key(seq):
    """Pick out the only element of this list, a list of keys"""
    return seq[0]
