from __future__ import annotations

""" Task specification for dask

This module contains the task specification for dask. It is used to represent
runnable (task) and non-runnable (data) nodes in a dask graph.

Simple examples of how to express tasks in dask
-----------------------------------------------

.. code-block:: python

    func("a", "b") ~ Task("key", func, "a", "b")

    [func("a"), func("b")] ~ [Task("key-1", func, "a"), Task("key-2", func, "b")]

    {"a": func("b")} ~ {"a": Task("a", func, "b")}

    "literal-string" ~ DataNode("key", "literal-string")


Keys, Aliases and TaskRefs
-------------------------

Keys are used to identify tasks in a dask graph. Every `GraphNode` instance has a
key attribute that _should_ reference the key in the dask graph.

.. code-block:: python

    {"key": Task("key", func, "a")}

Referencing other tasks is possible by using either one of `Alias` or a
`TaskRef`.

.. code-block:: python

    # TaskRef can be used to provide the name of the reference explicitly
    t = Task("key", func, TaskRef("key-1"))

    # If a task is still in scope, the method `ref` can be used for convenience
    t2 = Task("key2", func2, t.ref())


Executing a task
----------------

A task can be executed by calling it with a dictionary of values. The values
should contain the dependencies of the task.

.. code-block:: python

    t = Task("key", add, TaskRef("a"), TaskRef("b"))
    assert t.dependencies == {"a", "b"}
    t({"a": 1, "b": 2}) == 3

"""
import itertools
import sys
from collections import defaultdict
from collections.abc import Callable, Container, Iterable, Mapping, MutableMapping
from functools import lru_cache, partial
from typing import Any, TypeVar, cast

from dask.sizeof import sizeof
from dask.typing import Key as KeyType
from dask.utils import funcname, is_namedtuple_instance

_T = TypeVar("_T")
_T_GraphNode = TypeVar("_T_GraphNode", bound="GraphNode")
_T_Iterable = TypeVar("_T_Iterable", bound=Iterable)


# Ported from more-itertools
# https://github.com/more-itertools/more-itertools/blob/c8153e2801ade2527f3a6c8b623afae93f5a1ce1/more_itertools/recipes.py#L944-L973
def _batched(iterable, n, *, strict=False):
    """Batch data into tuples of length *n*. If the number of items in
    *iterable* is not divisible by *n*:
    * The last batch will be shorter if *strict* is ``False``.
    * :exc:`ValueError` will be raised if *strict* is ``True``.

    >>> list(batched('ABCDEFG', 3))
    [('A', 'B', 'C'), ('D', 'E', 'F'), ('G',)]

    On Python 3.13 and above, this is an alias for :func:`itertools.batched`.
    """
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := tuple(itertools.islice(it, n)):
        if strict and len(batch) != n:
            raise ValueError("batched(): incomplete batch")
        yield batch


if sys.hexversion >= 0x30D00A2:

    def batched(iterable, n, *, strict=False):
        return itertools.batched(iterable, n, strict=strict)

else:
    batched = _batched

    batched.__doc__ = _batched.__doc__
# End port


def identity(*args):
    return args


def _identity_cast(*args, typ):
    return typ(args)


_anom_count = itertools.count()


def parse_input(obj: Any) -> object:
    """Tokenize user input into GraphNode objects

    Note: This is similar to `convert_legacy_task` but does not
    - compare any values to a global set of known keys to infer references/futures
    - parse tuples and interprets them as runnable tasks
    - Deal with SubgraphCallables

    Parameters
    ----------
    obj : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    if isinstance(obj, GraphNode):
        return obj

    if isinstance(obj, TaskRef):
        return Alias(obj.key)

    if isinstance(obj, dict):
        parsed_dict = {k: parse_input(v) for k, v in obj.items()}
        if any(isinstance(v, GraphNode) for v in parsed_dict.values()):
            return Dict(parsed_dict)

    if isinstance(obj, (list, set, tuple)):
        parsed_collection = tuple(parse_input(o) for o in obj)
        if any(isinstance(o, GraphNode) for o in parsed_collection):
            if isinstance(obj, list):
                return List(*parsed_collection)
            if isinstance(obj, set):
                return Set(*parsed_collection)
            if isinstance(obj, tuple):
                if is_namedtuple_instance(obj):
                    return _wrap_namedtuple_task(None, obj, parse_input)
                return Tuple(*parsed_collection)

    return obj


def _wrap_namedtuple_task(k, obj, parser):
    if hasattr(obj, "__getnewargs_ex__"):
        new_args, kwargs = obj.__getnewargs_ex__()
        kwargs = {k: parser(v) for k, v in kwargs.items()}
    elif hasattr(obj, "__getnewargs__"):
        new_args = obj.__getnewargs__()
        kwargs = {}

    args_converted = parse_input(type(new_args)(map(parser, new_args)))

    return Task(
        k, partial(_instantiate_named_tuple, type(obj)), args_converted, Dict(kwargs)
    )


def _instantiate_named_tuple(typ, args, kwargs):
    return typ(*args, **kwargs)


class _MultiContainer(Container):
    container: tuple
    __slots__ = ("container",)

    def __init__(self, *container):
        self.container = container

    def __contains__(self, o: object) -> bool:
        for c in self.container:
            if o in c:
                return True
        return False


SubgraphType = None


def _execute_subgraph(inner_dsk, outkey, inkeys, *dependencies):
    final = {}
    final.update(inner_dsk)
    for k, v in zip(inkeys, dependencies):
        final[k] = DataNode(None, v)
    res = execute_graph(final, keys=[outkey])
    return res[outkey]


def convert_legacy_task(
    key: KeyType | None,
    task: _T,
    all_keys: Container,
) -> GraphNode | _T:
    if isinstance(task, GraphNode):
        return task

    if type(task) is tuple and task and callable(task[0]):
        func, args = task[0], task[1:]
        new_args = []
        new: object
        for a in args:
            if isinstance(a, dict):
                new = Dict(a)
            else:
                new = convert_legacy_task(None, a, all_keys)
            new_args.append(new)
        return Task(key, func, *new_args)
    try:
        if isinstance(task, (int, float, str, tuple)):
            if task in all_keys:
                if key is None:
                    return Alias(task)
                else:
                    return Alias(key, target=task)
    except TypeError:
        # Unhashable
        pass

    if isinstance(task, (list, tuple, set, frozenset)):
        if is_namedtuple_instance(task):
            return _wrap_namedtuple_task(
                key,
                task,
                partial(
                    convert_legacy_task,
                    None,
                    all_keys=all_keys,
                ),
            )
        else:
            parsed_args = tuple(convert_legacy_task(None, t, all_keys) for t in task)
            if any(isinstance(a, GraphNode) for a in parsed_args):
                return Task(key, _identity_cast, *parsed_args, typ=type(task))
            else:
                return cast(_T, type(task)(parsed_args))
    elif isinstance(task, TaskRef):
        if key is None:
            return Alias(task.key)
        else:
            return Alias(key, target=task.key)
    else:
        return task


def convert_legacy_graph(
    dsk: Mapping,
    all_keys: Container | None = None,
):
    if all_keys is None:
        all_keys = set(dsk)
    new_dsk = {}
    for k, arg in dsk.items():
        t = convert_legacy_task(k, arg, all_keys)
        if isinstance(t, Alias) and t.target == k:
            continue
        elif not isinstance(t, GraphNode):
            t = DataNode(k, t)
        new_dsk[k] = t
    return new_dsk


def resolve_aliases(dsk: dict, keys: set, dependents: dict) -> dict:
    """Remove trivial sequential alias chains

    Example:

        dsk = {'x': 1, 'y': Alias('x'), 'z': Alias('y')}

        resolve_aliases(dsk, {'z'}, {'x': {'y'}, 'y': {'z'}}) == {'z': 1}

    """
    if not keys:
        raise ValueError("No keys provided")
    dsk = dict(dsk)
    work = list(keys)
    seen = set()
    while work:
        k = work.pop()
        if k in seen or k not in dsk:
            continue
        seen.add(k)
        t = dsk[k]
        if isinstance(t, Alias):
            target_key = t.target
            # Rules for when we allow to collapse an alias
            # 1. The target key is not in the keys set. The keys set is what the
            #    user is requesting and by collapsing we'd no longer be able to
            #    return that result.
            # 2. The target key is in fact part of dsk. If it isnt' this could
            #    point to a persisted dependency and we cannot collapse it.
            # 3. The target key has only one dependent which is the key we're
            #    currently looking at. This means that there is a one to one
            #    relation between this and the target key in which case we can
            #    collapse them.
            #    Note: If target was an alias as well, we could continue with
            #    more advanced optimizations but this isn't implemented, yet
            if (
                target_key not in keys
                and target_key in dsk
                # Note: whenever we're performing a collapse, we're not updating
                # the dependents. The length == 1 should still be sufficient for
                # chains of these aliases
                and len(dependents[target_key]) == 1
            ):
                tnew = dsk.pop(target_key).copy()

                dsk[k] = tnew
                tnew.key = k
                if isinstance(tnew, Alias):
                    work.append(k)
                    seen.discard(k)
                else:
                    work.extend(tnew.dependencies)

        work.extend(t.dependencies)
    return dsk


class TaskRef:
    val: KeyType
    __slots__ = ("key",)

    def __init__(self, key: KeyType):
        self.key = key

    def __str__(self):
        return str(self.key)

    def __repr__(self):
        return f"{type(self).__name__}({self.key!r})"

    def __hash__(self) -> int:
        return hash(self.key)

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, TaskRef):
            return False
        return self.key == value.key


class GraphNode:
    key: KeyType
    _dependencies: frozenset

    __slots__ = tuple(__annotations__)

    def ref(self):
        return Alias(self.key)

    def copy(self):
        raise NotImplementedError

    @property
    def data_producer(self) -> bool:
        return False

    @property
    def dependencies(self) -> frozenset:
        return self._dependencies

    @property
    def block_fusion(self) -> bool:
        return False

    def _verify_values(self, values: tuple | dict) -> None:
        if not self.dependencies:
            return
        if missing := set(self.dependencies) - set(values):
            raise RuntimeError(f"Not enough arguments provided: missing keys {missing}")

    def __call__(self, values) -> Any:
        raise NotImplementedError("Not implemented")

    def __eq__(self, value: object) -> bool:
        if type(value) is not type(self):
            return False

        from dask.tokenize import tokenize

        return tokenize(self) == tokenize(value)

    @property
    def is_coro(self) -> bool:
        return False

    def __sizeof__(self) -> int:
        all_slots = self.get_all_slots()
        return sum(sizeof(getattr(self, sl)) for sl in all_slots) + sys.getsizeof(
            type(self)
        )

    def substitute(
        self, subs: dict[KeyType, KeyType | GraphNode], key: KeyType | None = None
    ) -> GraphNode:
        """Substitute a dependency with a new value. The new value either has to
        be a new valid key or a GraphNode to replace the dependency entirely.

        The GraphNode will not be mutated but instead a shallow copy will be
        returned. The substitution will be performed eagerly.

        Parameters
        ----------
        subs : dict[KeyType, KeyType  |  GraphNode]
            The mapping describing the substitutions to be made.
        key : KeyType | None, optional
            The key of the new GraphNode object. If None provided, the key of
            the old one will be reused.
        """
        raise NotImplementedError

    @staticmethod
    def fuse(*tasks: GraphNode, key: KeyType | None = None) -> GraphNode:
        """Fuse a set of tasks into a single task.

        The tasks are fused into a single task that will execute the tasks in a
        subgraph. The internal tasks are no longer accessible from the outside.

        All provided tasks must form a valid subgraph that will reduce to a
        single key. If multiple outputs are possible with the provided tasks, an
        exception will be raised.

        The tasks will not be rewritten but instead a new Task will be created
        that will merely reference the old task objects. This way, Task objects
        may be reused in multiple fused tasks.

        Parameters
        ----------
        key : KeyType | None, optional
            The key of the new Task object. If None provided, the key of the
            final task will be used.

        See also
        --------
        GraphNode.substitute : Easer substitution of dependencies
        """
        if any(t.key is None for t in tasks):
            raise ValueError("Cannot fuse tasks with missing keys")
        if len(tasks) == 1:
            return tasks[0].substitute({}, key=key)
        all_keys = set()
        all_deps: set[KeyType] = set()
        for t in tasks:
            all_deps.update(t.dependencies)
            all_keys.add(t.key)
        external_deps = tuple(sorted(all_deps - all_keys, key=hash))
        leafs = all_keys - all_deps
        if len(leafs) > 1:
            raise ValueError(f"Cannot fuse tasks with multiple outputs {leafs}")

        outkey = leafs.pop()
        return Task(
            key or outkey,
            _execute_subgraph,
            {t.key: t for t in tasks},
            outkey,
            external_deps,
            *(TaskRef(k) for k in external_deps),
            _data_producer=any(t.data_producer for t in tasks),
        )

    @classmethod
    @lru_cache
    def get_all_slots(cls):
        slots = list()
        for c in cls.mro():
            slots.extend(getattr(c, "__slots__", ()))
        # Interestingly, sorting this causes the nested containers to pickle
        # more efficiently
        return sorted(set(slots))


_no_deps: frozenset = frozenset()


class Alias(GraphNode):
    target: KeyType
    __slots__ = tuple(__annotations__)

    def __init__(
        self, key: KeyType | TaskRef, target: Alias | TaskRef | KeyType | None = None
    ):
        if isinstance(key, TaskRef):
            key = key.key
        self.key = key
        if target is None:
            target = key
        if isinstance(target, Alias):
            target = target.target
        if isinstance(target, TaskRef):
            target = target.key
        self.target = target
        self._dependencies = frozenset((self.target,))

    def __reduce__(self):
        return Alias, (self.key, self.target)

    def copy(self):
        return Alias(self.key, self.target)

    def substitute(
        self, subs: dict[KeyType, KeyType | GraphNode], key: KeyType | None = None
    ) -> GraphNode:
        if self.key in subs or self.target in subs:
            sub_key = subs.get(self.key, self.key)
            val = subs.get(self.target, self.target)
            if sub_key == self.key and val == self.target:
                return self
            if isinstance(val, GraphNode):
                return val.substitute({}, key=key)
            if key is None and isinstance(sub_key, GraphNode):
                raise RuntimeError(
                    f"Invalid substitution encountered {self.key!r} -> {sub_key}"
                )
            return Alias(key or sub_key, val)  # type: ignore
        return self

    def __dask_tokenize__(self):
        return (type(self).__name__, self.key, self.target)

    def __call__(self, values=()):
        self._verify_values(values)
        return values[self.target]

    def __repr__(self):
        if self.key != self.target:
            return f"Alias({self.key!r}->{self.target!r})"
        else:
            return f"Alias({self.key!r})"

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, Alias):
            return False
        if self.key != value.key:
            return False
        if self.target != value.target:
            return False
        return True


class DataNode(GraphNode):
    value: Any
    typ: type
    __slots__ = tuple(__annotations__)

    def __init__(self, key: Any, value: Any):
        if key is None:
            key = (type(value).__name__, next(_anom_count))
        self.key = key
        self.value = value
        self.typ = type(value)
        self._dependencies = _no_deps

    @property
    def data_producer(self) -> bool:
        return True

    def copy(self):
        return DataNode(self.key, self.value)

    def __call__(self, values=()):
        return self.value

    def __repr__(self):
        return f"DataNode({self.value!r})"

    def __reduce__(self):
        return (DataNode, (self.key, self.value))

    def __dask_tokenize__(self):
        from dask.base import tokenize

        return (type(self).__name__, tokenize(self.value))

    def substitute(
        self, subs: dict[KeyType, KeyType | GraphNode], key: KeyType | None = None
    ) -> DataNode:
        if key is not None and key != self.key:
            return DataNode(key, self.value)
        return self

    def __iter__(self):
        return iter(self.value)


def _call_recursively(obj, values):
    if isinstance(obj, GraphNode):
        return obj(values)
    elif isinstance(obj, dict):
        return {k: _call_recursively(v, values) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple, frozenset, set)):
        return type(obj)(_call_recursively(v, values) for v in obj)
    return obj


def _get_dependencies(obj: object) -> set | frozenset:
    if isinstance(obj, GraphNode):
        return obj.dependencies
    elif isinstance(obj, dict):
        if not obj:
            return _no_deps
        return set().union(*map(_get_dependencies, obj.values()))
    elif isinstance(obj, (list, tuple, frozenset, set)):
        if not obj:
            return _no_deps
        return set().union(*map(_get_dependencies, obj))
    return _no_deps


class Task(GraphNode):
    func: Callable
    args: tuple
    kwargs: dict
    _data_producer: bool
    _token: str | None
    _is_coro: bool | None
    _repr: str | None

    __slots__ = tuple(__annotations__)

    def __init__(
        self,
        key: Any,
        func: Callable,
        /,
        *args: Any,
        _data_producer: bool = False,
        _dependencies: set | frozenset | None = None,
        **kwargs: Any,
    ):
        self.key = key
        self.func = func
        if isinstance(func, Task):
            raise TypeError("Cannot nest tasks")
        self.args = tuple(
            Alias(obj.key) if isinstance(obj, TaskRef) else obj for obj in args
        )
        self.kwargs = {
            k: Alias(v.key) if isinstance(v, TaskRef) else v for k, v in kwargs.items()
        }
        if _dependencies is None:
            _dependencies = set()
            for a in itertools.chain(self.args, self.kwargs.values()):
                if isinstance(a, GraphNode):
                    _dependencies.update(a.dependencies)
        if _dependencies:
            self._dependencies = frozenset(_dependencies)
        else:
            self._dependencies = _no_deps
        self._is_coro = None
        self._token = None
        self._repr = None
        self._data_producer = _data_producer

    @property
    def data_producer(self) -> bool:
        return self._data_producer

    def copy(self):
        return type(self)(
            self.key,
            self.func,
            *self.args,
            _dependencies=self._dependencies,
            **self.kwargs,
        )

    def __hash__(self):
        return hash(self._get_token())

    def _get_token(self) -> str:
        if self._token:
            return self._token
        from dask.base import tokenize

        self._token = tokenize(
            (
                type(self).__name__,
                self.func,
                self.args,
                self.kwargs,
            )
        )
        return self._token

    def __dask_tokenize__(self):
        return self._get_token()

    def __repr__(self) -> str:
        # When `Task` is deserialized the constructor will not run and
        # `self._repr` is thus undefined.
        if not hasattr(self, "_repr") or not self._repr:
            head = funcname(self.func)
            tail = ")"
            label_size = 40
            args = self.args
            kwargs = self.kwargs
            if args or kwargs:
                label_size2 = int(
                    (label_size - len(head) - len(tail) - len(str(self.key)))
                    // (len(args) + len(kwargs))
                )
            if args:
                if label_size2 > 5:
                    args_repr = ", ".join(repr(t) for t in args)
                else:
                    args_repr = "..."
            else:
                args_repr = ""
            if kwargs:
                if label_size2 > 5:
                    kwargs_repr = ", " + ", ".join(
                        f"{k}={repr(v)}" for k, v in sorted(kwargs.items())
                    )
                else:
                    kwargs_repr = ", ..."
            else:
                kwargs_repr = ""
            self._repr = f"<Task {self.key!r} {head}({args_repr}{kwargs_repr}{tail}>"
        return self._repr

    def __call__(self, values=()):
        self._verify_values(values)
        new_argspec = tuple(
            a({k: values[k] for k in a.dependencies}) if isinstance(a, GraphNode) else a
            for a in self.args
        )
        if self.kwargs:
            kwargs = {
                k: (
                    kw({k: values[k] for k in kw.dependencies})
                    if isinstance(kw, GraphNode)
                    else kw
                )
                for k, kw in self.kwargs.items()
            }
            return self.func(*new_argspec, **kwargs)
        return self.func(*new_argspec)

    def __setstate__(self, state):
        slots = self.__class__.get_all_slots()
        for sl, val in zip(slots, state):
            setattr(self, sl, val)

    def __getstate__(self):
        slots = self.__class__.get_all_slots()
        return tuple(getattr(self, sl) for sl in slots)

    @property
    def is_coro(self):
        if self._is_coro is None:
            # Note: Can't use cached_property on objects without __dict__
            try:
                from distributed.utils import iscoroutinefunction

                self._is_coro = iscoroutinefunction(self.func)
            except Exception:
                self._is_coro = False
        return self._is_coro

    def substitute(
        self, subs: dict[KeyType, KeyType | GraphNode], key: KeyType | None = None
    ) -> Task:
        subs_filtered = {
            k: v for k, v in subs.items() if k in self.dependencies and k != v
        }
        if subs_filtered:
            new_args = tuple(
                a.substitute(subs_filtered) if isinstance(a, GraphNode) else a
                for a in self.args
            )
            new_kwargs = {
                k: v.substitute(subs_filtered) if isinstance(v, GraphNode) else v
                for k, v in self.kwargs.items()
            }
            return type(self)(
                key or self.key,
                self.func,
                *new_args,
                _data_producer=self.data_producer,
                **new_kwargs,
            )
        elif key is None or key == self.key:
            return self
        else:
            # Rename
            return type(self)(
                key,
                self.func,
                *self.args,
                _data_producer=self.data_producer,
                **self.kwargs,
            )


class NestedContainer(Task):
    constructor: Callable
    klass: type
    __slots__ = tuple(__annotations__)

    def __init__(
        self,
        /,
        *args: Any,
        _dependencies: set | frozenset | None = None,
        **kwargs: Any,
    ):
        if len(args) == 1 and isinstance(args[0], self.klass):
            args = args[0]  # type: ignore
        super().__init__(
            None,
            self.to_container,
            *args,
            constructor=self.constructor,
            _dependencies=_dependencies,
            **kwargs,
        )

    def __getstate__(self):
        state = super().__getstate__()
        state = list(state)
        slots = self.__class__.get_all_slots()
        ix = slots.index("kwargs")
        # The constructor as a kwarg is redundant since this is encoded in the
        # class itself. Serializing the builtin types is not trivial
        # This saves about 15% of overhead
        state[ix] = state[ix].copy()
        state[ix].pop("constructor", None)
        return state

    def __setstate__(self, state):
        super().__setstate__(state)
        self.kwargs["constructor"] = self.__class__.constructor
        return self

    def __repr__(self):
        return f"{type(self).__name__}({self.args})"

    def substitute(
        self, subs: dict[KeyType, KeyType | GraphNode], key: KeyType | None = None
    ) -> NestedContainer:
        subs_filtered = {
            k: v for k, v in subs.items() if k in self.dependencies and k != v
        }
        if not subs_filtered:
            return self
        return type(self)(
            *(
                a.substitute(subs_filtered) if isinstance(a, GraphNode) else a
                for a in self.args
            )
        )

    def __dask_tokenize__(self):
        from dask.tokenize import tokenize

        return (
            type(self).__name__,
            self.klass,
            sorted(tokenize(a) for a in self.args),
        )

        return super().__dask_tokenize__()

    @staticmethod
    def to_container(*args, constructor):
        return constructor(args)


class List(NestedContainer):
    constructor = klass = list


class Tuple(NestedContainer):
    constructor = klass = tuple


class Set(NestedContainer):
    constructor = klass = set


class Dict(NestedContainer):
    klass = dict

    def __init__(
        self, /, *args: Any, _dependencies: set | frozenset | None = None, **kwargs: Any
    ):
        if args:
            if len(args) > 1:
                raise ValueError("Dict can only take one positional argument")
            kwargs = args[0]
        super().__init__(
            *(itertools.chain(*kwargs.items())), _dependencies=_dependencies
        )

    def __repr__(self):
        values = ", ".join(f"{k}: {v}" for k, v in batched(self.args, 2, strict=True))
        return f"Dict({values})"

    def substitute(
        self, subs: dict[KeyType, KeyType | GraphNode], key: KeyType | None = None
    ) -> Dict:
        subs_filtered = {
            k: v for k, v in subs.items() if k in self.dependencies and k != v
        }
        if not subs_filtered:
            return self

        new = {}
        for k, v in batched(self.args, 2, strict=True):
            v = v.substitute(subs_filtered) if isinstance(v, GraphNode) else v
            new[k] = v
        return type(self)(new)

    @staticmethod
    def constructor(args):
        return dict(batched(args, 2, strict=True))


class DependenciesMapping(MutableMapping):
    def __init__(self, dsk):
        self.dsk = dsk
        self._removed = set()
        # Set a copy of dsk to avoid dct resizing
        self._cache = dsk.copy()
        self._cache.clear()

    def __getitem__(self, key):
        if (val := self._cache.get(key)) is not None:
            return val
        else:
            v = self.dsk[key]
            try:
                deps = v.dependencies
            except AttributeError:
                from dask.core import get_dependencies

                deps = get_dependencies(self.dsk, task=v)

            if self._removed:
                # deps is a frozenset but for good measure, let's not use -= since
                # that _may_ perform an inplace mutation
                deps = deps - self._removed
            self._cache[key] = deps
            return deps

    def __iter__(self):
        return iter(self.dsk)

    def __delitem__(self, key: Any) -> None:
        self._cache.clear()
        self._removed.add(key)

    def __setitem__(self, key: Any, value: Any) -> None:
        raise NotImplementedError

    def __len__(self) -> int:
        return len(self.dsk)


class _DevNullMapping(MutableMapping):
    def __getitem__(self, key):
        raise KeyError(key)

    def __setitem__(self, key, value):
        pass

    def __delitem__(self, key):
        pass

    def __len__(self):
        return 0

    def __iter__(self):
        return iter(())


def execute_graph(
    dsk: Iterable[GraphNode] | Mapping[KeyType, GraphNode],
    cache: MutableMapping[KeyType, object] | None = None,
    keys: Container[KeyType] | None = None,
) -> MutableMapping[KeyType, object]:
    """Execute a given graph.

    The graph is exceuted in topological order as defined by dask.order until
    all leaf nodes, i.e. nodes without any dependents, are reached. The returned
    dictionary contains the results of the leaf nodes.

    If keys are required that are not part of the graph, they can be provided in the `cache` argument.

    If `keys` is provided, the result will contain only values that are part of the `keys` set.

    """
    if isinstance(dsk, (list, tuple, set, frozenset)):
        dsk = {t.key: t for t in dsk}
    else:
        assert isinstance(dsk, dict)

    refcount: defaultdict[KeyType, int] = defaultdict(int)
    for vals in DependenciesMapping(dsk).values():
        for val in vals:
            refcount[val] += 1

    cache = cache or {}
    from dask.order import order

    priorities = order(dsk)

    for key, node in sorted(dsk.items(), key=lambda it: priorities[it[0]]):
        cache[key] = node(cache)
        for dep in node.dependencies:
            refcount[dep] -= 1
            if refcount[dep] == 0 and keys and dep not in keys:
                del cache[dep]

    return cache


def fuse_linear_task_spec(dsk, keys):
    """
    keys are the keys from the graph that are requested by a computation. We
    can't fuse those together.
    """
    from dask.core import reverse_dict
    from dask.optimization import default_fused_keys_renamer

    keys = set(keys)
    dependencies = DependenciesMapping(dsk)
    dependents = reverse_dict(dependencies)

    seen = set()
    result = {}

    for key in dsk:
        if key in seen:
            continue

        seen.add(key)

        deps = dependencies[key]
        dependents_key = dependents[key]

        if len(deps) != 1 and len(dependents_key) != 1 or dsk[key].block_fusion:
            result[key] = dsk[key]
            continue

        linear_chain = [dsk[key]]
        top_key = key

        # Walk towards the leafs as long as the nodes have a single dependency
        # and a single dependent, we can't fuse two nodes of an intermediate node
        # is the source for 2 dependents
        while len(deps) == 1:
            (new_key,) = deps
            if new_key in seen:
                break
            seen.add(new_key)
            if new_key not in dsk:
                # This can happen if a future is in the graph, the dependency mapping
                # adds the key that is referenced by the future as a dependency
                # see test_futures_to_delayed_array
                break
            if (
                len(dependents[new_key]) != 1
                or dsk[new_key].block_fusion
                or new_key in keys
            ):
                result[new_key] = dsk[new_key]
                break
            # backwards comp for new names, temporary until is_rootish is removed
            linear_chain.insert(0, dsk[new_key])
            deps = dependencies[new_key]

        # Walk the tree towards the root as long as the nodes have a single dependent
        # and a single dependency, we can't fuse two nodes if node has multiple
        # dependencies
        while len(dependents_key) == 1 and top_key not in keys:
            new_key = dependents_key.pop()
            if new_key in seen:
                break
            seen.add(new_key)
            if len(dependencies[new_key]) != 1 or dsk[new_key].block_fusion:
                # Exit if the dependent has multiple dependencies, triangle
                result[new_key] = dsk[new_key]
                break
            linear_chain.append(dsk[new_key])
            top_key = new_key
            dependents_key = dependents[new_key]

        if len(linear_chain) == 1:
            result[top_key] = linear_chain[0]
        else:
            # Renaming the keys is necessary to preserve the rootish detection for now
            renamed_key = default_fused_keys_renamer([tsk.key for tsk in linear_chain])
            result[renamed_key] = Task.fuse(*linear_chain, key=renamed_key)
            if renamed_key != top_key:
                # Having the same prefixes can result in the same key, i.e. getitem-hash -> getitem-hash
                result[top_key] = Alias(top_key, target=renamed_key)
    return result


def cull(
    dsk: dict[KeyType, GraphNode], keys: Iterable[KeyType]
) -> dict[KeyType, GraphNode]:
    work = set(keys)
    seen: set[KeyType] = set()
    dsk2 = {}
    wpop = work.pop
    wupdate = work.update
    sadd = seen.add
    while work:
        k = wpop()
        if k in seen or k not in dsk:
            continue
        sadd(k)
        dsk2[k] = v = dsk[k]
        wupdate(v.dependencies)
    return dsk2
