from __future__ import annotations

import io
import itertools
import math
import operator
import uuid
import warnings
from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from functools import partial, reduce, wraps
from random import Random
from urllib.request import urlopen

import tlz as toolz
from fsspec.core import open_files
from tlz import (
    accumulate,
    compose,
    count,
    curry,
    first,
    frequencies,
    groupby,
    join,
    merge,
    merge_with,
    partition_all,
    peek,
    pluck,
    reduceby,
    remove,
    second,
    take,
    topk,
    unique,
    valmap,
)

from dask import config
from dask.bag import chunk
from dask.bag.avro import to_avro
from dask.base import (
    DaskMethodsMixin,
    dont_optimize,
    named_schedulers,
    replace_name_in_key,
    tokenize,
)
from dask.blockwise import blockwise
from dask.context import globalmethod
from dask.core import flatten, get_dependencies, istask, quote, reverse_dict
from dask.delayed import Delayed, unpack_collections
from dask.highlevelgraph import HighLevelGraph
from dask.optimization import cull, fuse, inline
from dask.sizeof import sizeof
from dask.typing import Graph, NestedKeys, no_default
from dask.utils import (
    apply,
    digit,
    ensure_bytes,
    ensure_dict,
    ensure_unicode,
    funcname,
    get_default_shuffle_method,
    insert,
    iter_chunks,
    key_split,
    parse_bytes,
    system_encoding,
    takes_multiple_arguments,
)

DEFAULT_GET = named_schedulers.get("processes", named_schedulers["sync"])

no_result = type(
    "no_result", (object,), {"__slots__": (), "__reduce__": lambda self: "no_result"}
)


def lazify_task(task, start=True):
    """
    Given a task, remove unnecessary calls to ``list`` and ``reify``.

    This traverses tasks and small lists.  We choose not to traverse down lists
    of size >= 50 because it is unlikely that sequences this long contain other
    sequences in practice.

    Examples
    --------
    >>> def inc(x):
    ...     return x + 1
    >>> task = (sum, (list, (map, inc, [1, 2, 3])))
    >>> lazify_task(task)  # doctest: +ELLIPSIS
    (<built-in function sum>, (<class 'map'>, <function inc at ...>, [1, 2, 3]))
    """
    if type(task) is list and len(task) < 50:
        return [lazify_task(arg, False) for arg in task]
    if not istask(task):
        return task
    head, tail = task[0], task[1:]
    if not start and head in (list, reify):
        task = task[1]
        return lazify_task(*tail, start=False)
    else:
        return (head,) + tuple(lazify_task(arg, False) for arg in tail)


def lazify(dsk):
    """
    Remove unnecessary calls to ``list`` in tasks.

    See Also
    --------
    dask.bag.core.lazify_task
    """
    return valmap(lazify_task, dsk)


def inline_singleton_lists(dsk, keys, dependencies=None):
    """Inline lists that are only used once.

    >>> d = {'b': (list, 'a'),
    ...      'c': (sum, 'b', 1)}
    >>> inline_singleton_lists(d, 'c')
    {'c': (<built-in function sum>, (<class 'list'>, 'a'), 1)}

    Pairs nicely with lazify afterwards.
    """
    if dependencies is None:
        dependencies = {k: get_dependencies(dsk, task=v) for k, v in dsk.items()}
    dependents = reverse_dict(dependencies)

    inline_keys = {
        k
        for k, v in dsk.items()
        if istask(v) and v and v[0] is list and len(dependents[k]) == 1
    }
    inline_keys.difference_update(flatten(keys))
    dsk = inline(dsk, inline_keys, inline_constants=False)
    for k in inline_keys:
        del dsk[k]
    return dsk


def optimize(dsk, keys, fuse_keys=None, rename_fused_keys=None, **kwargs):
    """Optimize a dask from a dask Bag."""
    dsk = ensure_dict(dsk)
    dsk2, dependencies = cull(dsk, keys)
    kwargs = {}
    if rename_fused_keys is not None:
        kwargs["rename_keys"] = rename_fused_keys
    dsk3, dependencies = fuse(dsk2, keys + (fuse_keys or []), dependencies, **kwargs)
    dsk4 = inline_singleton_lists(dsk3, keys, dependencies)
    dsk5 = lazify(dsk4)
    return dsk5


def _to_textfiles_chunk(data, lazy_file, last_endline):
    with lazy_file as f:
        if isinstance(f, io.TextIOWrapper):
            endline = "\n"
            ensure = ensure_unicode
        else:
            endline = b"\n"
            ensure = ensure_bytes
        started = False
        for d in data:
            if started:
                f.write(endline)
            else:
                started = True
            f.write(ensure(d))
        if last_endline:
            f.write(endline)


def to_textfiles(
    b,
    path,
    name_function=None,
    compression="infer",
    encoding=system_encoding,
    compute=True,
    storage_options=None,
    last_endline=False,
    **kwargs,
):
    """Write dask Bag to disk, one filename per partition, one line per element.

    **Paths**: This will create one file for each partition in your bag. You
    can specify the filenames in a variety of ways.

    Use a globstring

    >>> b.to_textfiles('/path/to/data/*.json.gz')  # doctest: +SKIP

    The * will be replaced by the increasing sequence 1, 2, ...

    ::

        /path/to/data/0.json.gz
        /path/to/data/1.json.gz

    Use a globstring and a ``name_function=`` keyword argument.  The
    name_function function should expect an integer and produce a string.
    Strings produced by name_function must preserve the order of their
    respective partition indices.

    >>> from datetime import date, timedelta
    >>> def name(i):
    ...     return str(date(2015, 1, 1) + i * timedelta(days=1))

    >>> name(0)
    '2015-01-01'
    >>> name(15)
    '2015-01-16'

    >>> b.to_textfiles('/path/to/data/*.json.gz', name_function=name)  # doctest: +SKIP

    ::

        /path/to/data/2015-01-01.json.gz
        /path/to/data/2015-01-02.json.gz
        ...

    You can also provide an explicit list of paths.

    >>> paths = ['/path/to/data/alice.json.gz', '/path/to/data/bob.json.gz', ...]  # doctest: +SKIP
    >>> b.to_textfiles(paths) # doctest: +SKIP

    **Compression**: Filenames with extensions corresponding to known
    compression algorithms (gz, bz2) will be compressed accordingly.

    **Bag Contents**: The bag calling ``to_textfiles`` must be a bag of
    text strings. For example, a bag of dictionaries could be written to
    JSON text files by mapping ``json.dumps`` on to the bag first, and
    then calling ``to_textfiles`` :

    >>> b_dict.map(json.dumps).to_textfiles("/path/to/data/*.json")  # doctest: +SKIP

    **Last endline**: By default the last line does not end with a newline
    character. Pass ``last_endline=True`` to invert the default.
    """
    mode = "wb" if encoding is None else "wt"
    files = open_files(
        path,
        compression=compression,
        mode=mode,
        encoding=encoding,
        name_function=name_function,
        num=b.npartitions,
        **(storage_options or {}),
    )

    name = "to-textfiles-" + uuid.uuid4().hex
    dsk = {
        (name, i): (_to_textfiles_chunk, (b.name, i), f, last_endline)
        for i, f in enumerate(files)
    }
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[b])
    out = type(b)(graph, name, b.npartitions)

    if compute:
        out.compute(**kwargs)
        return [f.path for f in files]
    else:
        return out.to_delayed()


def finalize(results):
    if not results:
        return results
    if isinstance(results, Iterator):
        results = list(results)
    if isinstance(results[0], Iterable) and not isinstance(results[0], str):
        results = toolz.concat(results)
    if isinstance(results, Iterator):
        results = list(results)
    return results


def finalize_item(results):
    return results[0]


class StringAccessor:
    """String processing functions

    Examples
    --------

    >>> import dask.bag as db
    >>> b = db.from_sequence(['Alice Smith', 'Bob Jones', 'Charlie Smith'])
    >>> list(b.str.lower())
    ['alice smith', 'bob jones', 'charlie smith']

    >>> list(b.str.match('*Smith'))
    ['Alice Smith', 'Charlie Smith']

    >>> list(b.str.split(' '))
    [['Alice', 'Smith'], ['Bob', 'Jones'], ['Charlie', 'Smith']]
    """

    def __init__(self, bag):
        self._bag = bag

    def __dir__(self):
        return sorted(set(dir(type(self)) + dir(str)))

    def _strmap(self, key, *args, **kwargs):
        return self._bag.map(operator.methodcaller(key, *args, **kwargs))

    def __getattr__(self, key):
        try:
            return object.__getattribute__(self, key)
        except AttributeError:
            if key in dir(str):
                func = getattr(str, key)
                return robust_wraps(func)(partial(self._strmap, key))
            else:
                raise

    def match(self, pattern):
        """Filter strings by those that match a pattern.

        Examples
        --------

        >>> import dask.bag as db
        >>> b = db.from_sequence(['Alice Smith', 'Bob Jones', 'Charlie Smith'])
        >>> list(b.str.match('*Smith'))
        ['Alice Smith', 'Charlie Smith']

        See Also
        --------
        fnmatch.fnmatch
        """
        from fnmatch import fnmatch

        return self._bag.filter(partial(fnmatch, pat=pattern))


def robust_wraps(wrapper):
    """A weak version of wraps that only copies doc."""

    def _(wrapped):
        wrapped.__doc__ = wrapper.__doc__
        return wrapped

    return _


class Item(DaskMethodsMixin):
    def __init__(self, dsk, key, layer=None):
        self.dask = dsk
        self.key = key
        self.name = key

        # NOTE: Layer only used by `Item.from_delayed`, to handle Delayed objects created by other collections.
        # e.g.: Item.from_delayed(da.ones(1).to_delayed()[0])
        # See Delayed.__init__
        self._layer = layer or key
        if isinstance(dsk, HighLevelGraph) and self._layer not in dsk.layers:
            raise ValueError(
                f"Layer {self._layer} not in the HighLevelGraph's layers: {list(dsk.layers)}"
            )

    def __dask_graph__(self) -> Graph:
        return self.dask

    def __dask_keys__(self) -> NestedKeys:
        return [self.key]

    def __dask_layers__(self) -> Sequence[str]:
        return (self._layer,)

    def __dask_tokenize__(self):
        return self.key

    __dask_optimize__ = globalmethod(optimize, key="bag_optimize", falsey=dont_optimize)
    __dask_scheduler__ = staticmethod(DEFAULT_GET)

    def __dask_postcompute__(self):
        return finalize_item, ()

    def __dask_postpersist__(self):
        return self._rebuild, ()

    def _rebuild(self, dsk, *, rename=None):
        key = replace_name_in_key(self.key, rename) if rename else self.key
        return Item(dsk, key)

    @staticmethod
    def from_delayed(value):
        """Create bag item from a dask.delayed value.

        See ``dask.bag.from_delayed`` for details
        """
        from dask.delayed import Delayed, delayed

        if not isinstance(value, Delayed) and hasattr(value, "key"):
            value = delayed(value)
        assert isinstance(value, Delayed)
        return Item(value.dask, value.key, layer=value.__dask_layers__()[0])

    @property
    def _args(self):
        return (self.dask, self.key)

    def __getstate__(self):
        return self._args

    def __setstate__(self, state):
        self.dask, self.key = state

    def apply(self, func):
        name = "{}-{}".format(funcname(func), tokenize(self, func, "apply"))
        dsk = {name: (func, self.key)}
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return Item(graph, name)

    __int__ = __float__ = __complex__ = __bool__ = DaskMethodsMixin.compute

    def to_delayed(self, optimize_graph=True):
        """Convert into a ``dask.delayed`` object.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.
        """
        from dask.delayed import Delayed

        dsk = self.__dask_graph__()
        if optimize_graph:
            dsk = self.__dask_optimize__(dsk, self.__dask_keys__())
        return Delayed(self.key, dsk, layer=self._layer)


class Bag(DaskMethodsMixin):
    """Parallel collection of Python objects

    Examples
    --------
    Create Bag from sequence

    >>> import dask.bag as db
    >>> b = db.from_sequence(range(5))
    >>> list(b.filter(lambda x: x % 2 == 0).map(lambda x: x * 10))
    [0, 20, 40]

    Create Bag from filename or globstring of filenames

    >>> b = db.read_text('/path/to/mydata.*.json.gz').map(json.loads)  # doctest: +SKIP

    Create manually (expert use)

    >>> dsk = {('x', 0): (range, 5),
    ...        ('x', 1): (range, 5),
    ...        ('x', 2): (range, 5)}
    >>> b = db.Bag(dsk, 'x', npartitions=3)

    >>> sorted(b.map(lambda x: x * 10))
    [0, 0, 0, 10, 10, 10, 20, 20, 20, 30, 30, 30, 40, 40, 40]

    >>> int(b.fold(lambda x, y: x + y))
    30
    """

    def __init__(self, dsk: Graph, name: str, npartitions: int):
        if not isinstance(dsk, HighLevelGraph):
            dsk = HighLevelGraph.from_collections(name, dsk, dependencies=[])
        self.dask = dsk
        self.name = name
        self.npartitions = npartitions

    def __dask_graph__(self) -> Graph:
        return self.dask

    def __dask_keys__(self) -> NestedKeys:
        return [(self.name, i) for i in range(self.npartitions)]

    def __dask_layers__(self) -> Sequence[str]:
        return (self.name,)

    def __dask_tokenize__(self):
        return self.name

    __dask_optimize__ = globalmethod(optimize, key="bag_optimize", falsey=dont_optimize)
    __dask_scheduler__ = staticmethod(DEFAULT_GET)

    def __dask_postcompute__(self):
        return finalize, ()

    def __dask_postpersist__(self):
        return self._rebuild, ()

    def _rebuild(self, dsk, *, rename=None):
        name = self.name
        if rename:
            name = rename.get(name, name)
        return type(self)(dsk, name, self.npartitions)

    def __str__(self):
        return "dask.bag<%s, npartitions=%d>" % (key_split(self.name), self.npartitions)

    __repr__ = __str__

    str = property(fget=StringAccessor)

    def map(self, func, *args, **kwargs):
        """Apply a function elementwise across one or more bags.

        Note that all ``Bag`` arguments must be partitioned identically.

        Parameters
        ----------
        func : callable
        *args, **kwargs : Bag, Item, or object
            Extra arguments and keyword arguments to pass to ``func`` *after*
            the calling bag instance. Non-Bag args/kwargs are broadcasted
            across all calls to ``func``.

        Notes
        -----
        For calls with multiple `Bag` arguments, corresponding partitions
        should have the same length; if they do not, the call will error at
        compute time.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(range(5), npartitions=2)
        >>> b2 = db.from_sequence(range(5, 10), npartitions=2)

        Apply a function to all elements in a bag:

        >>> b.map(lambda x: x + 1).compute()
        [1, 2, 3, 4, 5]

        Apply a function with arguments from multiple bags:

        >>> from operator import add
        >>> b.map(add, b2).compute()
        [5, 7, 9, 11, 13]

        Non-bag arguments are broadcast across all calls to the mapped
        function:

        >>> b.map(add, 1).compute()
        [1, 2, 3, 4, 5]

        Keyword arguments are also supported, and have the same semantics as
        regular arguments:

        >>> def myadd(x, y=0):
        ...     return x + y
        >>> b.map(myadd, y=b2).compute()
        [5, 7, 9, 11, 13]
        >>> b.map(myadd, y=1).compute()
        [1, 2, 3, 4, 5]

        Both arguments and keyword arguments can also be instances of
        ``dask.bag.Item``. Here we'll add the max value in the bag to each
        element:

        >>> b.map(myadd, b.max()).compute()
        [4, 5, 6, 7, 8]
        """
        return bag_map(func, self, *args, **kwargs)

    def starmap(self, func, **kwargs):
        """Apply a function using argument tuples from the given bag.

        This is similar to ``itertools.starmap``, except it also accepts
        keyword arguments. In pseudocode, this is could be written as:

        >>> def starmap(func, bag, **kwargs):
        ...     return (func(*args, **kwargs) for args in bag)

        Parameters
        ----------
        func : callable
        **kwargs : Item, Delayed, or object, optional
            Extra keyword arguments to pass to ``func``. These can either be
            normal objects, ``dask.bag.Item``, or ``dask.delayed.Delayed``.

        Examples
        --------
        >>> import dask.bag as db
        >>> data = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10)]
        >>> b = db.from_sequence(data, npartitions=2)

        Apply a function to each argument tuple:

        >>> from operator import add
        >>> b.starmap(add).compute()
        [3, 7, 11, 15, 19]

        Apply a function to each argument tuple, with additional keyword
        arguments:

        >>> def myadd(x, y, z=0):
        ...     return x + y + z
        >>> b.starmap(myadd, z=10).compute()
        [13, 17, 21, 25, 29]

        Keyword arguments can also be instances of ``dask.bag.Item`` or
        ``dask.delayed.Delayed``:

        >>> max_second = b.pluck(1).max()
        >>> max_second.compute()
        10
        >>> b.starmap(myadd, z=max_second).compute()
        [13, 17, 21, 25, 29]
        """
        name = "{}-{}".format(funcname(func), tokenize(self, func, "starmap", **kwargs))
        dependencies = [self]
        if kwargs:
            kwargs, collections = unpack_scalar_dask_kwargs(kwargs)
            dependencies.extend(collections)

        dsk = {
            (name, i): (reify, (starmap_chunk, func, (self.name, i), kwargs))
            for i in range(self.npartitions)
        }
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=dependencies)
        return type(self)(graph, name, self.npartitions)

    @property
    def _args(self):
        return (self.dask, self.name, self.npartitions)

    def __getstate__(self):
        return self._args

    def __setstate__(self, state):
        self.dask, self.name, self.npartitions = state

    def filter(self, predicate):
        """Filter elements in collection by a predicate function.

        >>> def iseven(x):
        ...     return x % 2 == 0

        >>> import dask.bag as db
        >>> b = db.from_sequence(range(5))
        >>> list(b.filter(iseven))
        [0, 2, 4]
        """
        name = f"filter-{funcname(predicate)}-{tokenize(self, predicate)}"
        dsk = {
            (name, i): (reify, (filter, predicate, (self.name, i)))
            for i in range(self.npartitions)
        }
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return type(self)(graph, name, self.npartitions)

    def random_sample(self, prob, random_state=None):
        """Return elements from bag with probability of ``prob``.

        Parameters
        ----------
        prob : float
            A float between 0 and 1, representing the probability that each
            element will be returned.
        random_state : int or random.Random, optional
            If an integer, will be used to seed a new ``random.Random`` object.
            If provided, results in deterministic sampling.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(range(10))
        >>> b.random_sample(0.5, 43).compute()
        [0, 1, 3, 4, 7, 9]
        >>> b.random_sample(0.5, 43).compute()
        [0, 1, 3, 4, 7, 9]
        """
        if not 0 <= prob <= 1:
            raise ValueError("prob must be a number in the interval [0, 1]")
        if not isinstance(random_state, Random):
            random_state = Random(random_state)

        name = "random-sample-%s" % tokenize(self, prob, random_state.getstate())
        state_data = random_state_data_python(self.npartitions, random_state)
        dsk = {
            (name, i): (reify, (random_sample, (self.name, i), state, prob))
            for i, state in zip(range(self.npartitions), state_data)
        }
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return type(self)(graph, name, self.npartitions)

    def remove(self, predicate):
        """Remove elements in collection that match predicate.

        >>> def iseven(x):
        ...     return x % 2 == 0

        >>> import dask.bag as db
        >>> b = db.from_sequence(range(5))
        >>> list(b.remove(iseven))
        [1, 3]
        """
        name = f"remove-{funcname(predicate)}-{tokenize(self, predicate)}"
        dsk = {
            (name, i): (reify, (remove, predicate, (self.name, i)))
            for i in range(self.npartitions)
        }
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return type(self)(graph, name, self.npartitions)

    def map_partitions(self, func, *args, **kwargs):
        """Apply a function to every partition across one or more bags.

        Note that all ``Bag`` arguments must be partitioned identically.

        Parameters
        ----------
        func : callable
            The function to be called on every partition.
            This function should expect an ``Iterator`` or ``Iterable`` for
            every partition and should return an ``Iterator`` or ``Iterable``
            in return.
        *args, **kwargs : Bag, Item, Delayed, or object
            Arguments and keyword arguments to pass to ``func``.
            Partitions from this bag will be the first argument, and these will
            be passed *after*.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(range(1, 101), npartitions=10)
        >>> def div(nums, den=1):
        ...     return [num / den for num in nums]

        Using a python object:

        >>> hi = b.max().compute()
        >>> hi
        100
        >>> b.map_partitions(div, den=hi).take(5)
        (0.01, 0.02, 0.03, 0.04, 0.05)

        Using an ``Item``:

        >>> b.map_partitions(div, den=b.max()).take(5)
        (0.01, 0.02, 0.03, 0.04, 0.05)

        Note that while both versions give the same output, the second forms a
        single graph, and then computes everything at once, and in some cases
        may be more efficient.
        """
        return map_partitions(func, self, *args, **kwargs)

    def pluck(self, key, default=no_default):
        """Select item from all tuples/dicts in collection.

        >>> import dask.bag as db
        >>> b = db.from_sequence([{'name': 'Alice', 'credits': [1, 2, 3]},
        ...                       {'name': 'Bob',   'credits': [10, 20]}])
        >>> list(b.pluck('name'))
        ['Alice', 'Bob']
        >>> list(b.pluck('credits').pluck(0))
        [1, 10]
        """
        name = "pluck-" + tokenize(self, key, default)
        key = quote(key)
        if default is no_default:
            dsk = {
                (name, i): (list, (pluck, key, (self.name, i)))
                for i in range(self.npartitions)
            }
        else:
            dsk = {
                (name, i): (list, (pluck, key, (self.name, i), default))
                for i in range(self.npartitions)
            }
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return type(self)(graph, name, self.npartitions)

    def unzip(self, n):
        """Transform a bag of tuples to ``n`` bags of their elements.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence([(i, i + 1, i + 2) for i in range(10)])
        >>> first, second, third = b.unzip(3)
        >>> isinstance(first, db.Bag)
        True
        >>> first.compute()
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        Note that this is equivalent to:

        >>> first, second, third = (b.pluck(i) for i in range(3))
        """
        return tuple(self.pluck(i) for i in range(n))

    @wraps(to_textfiles)
    def to_textfiles(
        self,
        path,
        name_function=None,
        compression="infer",
        encoding=system_encoding,
        compute=True,
        storage_options=None,
        last_endline=False,
        **kwargs,
    ):
        return to_textfiles(
            self,
            path,
            name_function,
            compression,
            encoding,
            compute,
            storage_options=storage_options,
            last_endline=last_endline,
            **kwargs,
        )

    @wraps(to_avro)
    def to_avro(
        self,
        filename,
        schema,
        name_function=None,
        storage_options=None,
        codec="null",
        sync_interval=16000,
        metadata=None,
        compute=True,
        **kwargs,
    ):
        return to_avro(
            self,
            filename,
            schema,
            name_function,
            storage_options,
            codec,
            sync_interval,
            metadata,
            compute,
            **kwargs,
        )

    def fold(
        self, binop, combine=None, initial=no_default, split_every=None, out_type=Item
    ):
        """Parallelizable reduction

        Fold is like the builtin function ``reduce`` except that it works in
        parallel.  Fold takes two binary operator functions, one to reduce each
        partition of our dataset and another to combine results between
        partitions

        1.  ``binop``: Binary operator to reduce within each partition
        2.  ``combine``:  Binary operator to combine results from binop

        Sequentially this would look like the following:

        >>> intermediates = [reduce(binop, part) for part in partitions]  # doctest: +SKIP
        >>> final = reduce(combine, intermediates)  # doctest: +SKIP

        If only one function is given then it is used for both functions
        ``binop`` and ``combine`` as in the following example to compute the
        sum:

        >>> def add(x, y):
        ...     return x + y

        >>> import dask.bag as db
        >>> b = db.from_sequence(range(5))
        >>> b.fold(add).compute()
        10

        In full form we provide both binary operators as well as their default
        arguments

        >>> b.fold(binop=add, combine=add, initial=0).compute()
        10

        More complex binary operators are also doable

        >>> def add_to_set(acc, x):
        ...     ''' Add new element x to set acc '''
        ...     return acc | set([x])
        >>> b.fold(add_to_set, set.union, initial=set()).compute()
        {0, 1, 2, 3, 4}

        See Also
        --------

        Bag.foldby
        """
        combine = combine or binop
        if initial is not no_default:
            return self.reduction(
                curry(_reduce, binop, initial=initial),
                curry(_reduce, combine),
                split_every=split_every,
                out_type=out_type,
            )
        else:
            from tlz.curried import reduce

            return self.reduction(
                reduce(binop),
                reduce(combine),
                split_every=split_every,
                out_type=out_type,
            )

    def frequencies(self, split_every=None, sort=False):
        """Count number of occurrences of each distinct element.

        >>> import dask.bag as db
        >>> b = db.from_sequence(['Alice', 'Bob', 'Alice'])
        >>> dict(b.frequencies())       # doctest: +SKIP
        {'Alice': 2, 'Bob', 1}
        """
        result = self.reduction(
            frequencies,
            merge_frequencies,
            out_type=Bag,
            split_every=split_every,
            name="frequencies",
        ).map_partitions(dictitems)
        if sort:
            result = result.map_partitions(sorted, key=second, reverse=True)
        return result

    def topk(self, k, key=None, split_every=None):
        """K largest elements in collection

        Optionally ordered by some key function

        >>> import dask.bag as db
        >>> b = db.from_sequence([10, 3, 5, 7, 11, 4])
        >>> list(b.topk(2))
        [11, 10]

        >>> list(b.topk(2, lambda x: -x))
        [3, 4]
        """
        if key:
            if callable(key) and takes_multiple_arguments(key):
                key = partial(apply, key)
            func = partial(topk, k, key=key)
        else:
            func = partial(topk, k)
        return self.reduction(
            func,
            compose(func, toolz.concat),
            out_type=Bag,
            split_every=split_every,
            name="topk",
        )

    def distinct(self, key=None):
        """Distinct elements of collection

        Unordered without repeats.

        Parameters
        ----------
        key: {callable,str}
            Defines uniqueness of items in bag by calling ``key`` on each item.
            If a string is passed ``key`` is considered to be ``lambda x: x[key]``.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(['Alice', 'Bob', 'Alice'])
        >>> sorted(b.distinct())
        ['Alice', 'Bob']
        >>> b = db.from_sequence([{'name': 'Alice'}, {'name': 'Bob'}, {'name': 'Alice'}])
        >>> b.distinct(key=lambda x: x['name']).compute()
        [{'name': 'Alice'}, {'name': 'Bob'}]
        >>> b.distinct(key='name').compute()
        [{'name': 'Alice'}, {'name': 'Bob'}]
        """
        func = chunk_distinct if key is None else partial(chunk_distinct, key=key)
        agg = merge_distinct if key is None else partial(merge_distinct, key=key)
        return self.reduction(func, agg, out_type=Bag, name="distinct")

    def reduction(
        self, perpartition, aggregate, split_every=None, out_type=Item, name=None
    ):
        """Reduce collection with reduction operators.

        Parameters
        ----------
        perpartition: function
            reduction to apply to each partition
        aggregate: function
            reduction to apply to the results of all partitions
        split_every: int (optional)
            Group partitions into groups of this size while performing reduction
            Defaults to 8
        out_type: {Bag, Item}
            The out type of the result, Item if a single element, Bag if a list
            of elements.  Defaults to Item.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(range(10))
        >>> b.reduction(sum, sum).compute()
        45
        """
        if split_every is None:
            split_every = 8
        if split_every is False:
            split_every = self.npartitions

        token = tokenize(self, perpartition, aggregate, split_every)
        a = f"{name or funcname(perpartition)}-part-{token}"
        is_last = self.npartitions == 1
        dsk = {
            (a, i): (empty_safe_apply, perpartition, (self.name, i), is_last)
            for i in range(self.npartitions)
        }
        k = self.npartitions
        b = a
        fmt = f"{name or funcname(aggregate)}-aggregate-{token}"
        depth = 0

        while k > split_every:
            c = fmt + str(depth)
            for i, inds in enumerate(partition_all(split_every, range(k))):
                dsk[(c, i)] = (
                    empty_safe_aggregate,
                    aggregate,
                    [(b, j) for j in inds],
                    False,
                )

            k = i + 1
            b = c
            depth += 1

        dsk[(fmt, 0)] = (
            empty_safe_aggregate,
            aggregate,
            [(b, j) for j in range(k)],
            True,
        )

        graph = HighLevelGraph.from_collections(fmt, dsk, dependencies=[self])
        if out_type is Item:
            dsk[fmt] = dsk.pop((fmt, 0))
            return Item(graph, fmt)
        else:
            return Bag(graph, fmt, 1)

    def sum(self, split_every=None):
        """Sum all elements"""
        return self.reduction(sum, sum, split_every=split_every)

    def max(self, split_every=None):
        """Maximum element"""
        return self.reduction(max, max, split_every=split_every)

    def min(self, split_every=None):
        """Minimum element"""
        return self.reduction(min, min, split_every=split_every)

    def any(self, split_every=None):
        """Are any of the elements truthy?

        Examples
        --------
        >>> import dask.bag as db
        >>> bool_bag = db.from_sequence([True, True, False])
        >>> bool_bag.any().compute()
        True
        """
        return self.reduction(any, any, split_every=split_every)

    def all(self, split_every=None):
        """Are all elements truthy?

        Examples
        --------
        >>> import dask.bag as db
        >>> bool_bag = db.from_sequence([True, True, False])
        >>> bool_bag.all().compute()
        False
        """
        return self.reduction(all, all, split_every=split_every)

    def count(self, split_every=None):
        """Count the number of elements.

        Examples
        --------
        >>> import dask.bag as db
        >>> numbers = db.from_sequence([1, 2, 3])
        >>> numbers.count().compute()
        3
        """
        return self.reduction(count, sum, split_every=split_every)

    def mean(self):
        """Arithmetic mean"""

        def mean_chunk(seq):
            total, n = 0.0, 0
            for x in seq:
                total += x
                n += 1
            return total, n

        def mean_aggregate(x):
            totals, counts = list(zip(*x))
            return 1.0 * sum(totals) / sum(counts)

        return self.reduction(mean_chunk, mean_aggregate, split_every=False)

    def var(self, ddof=0):
        """Variance"""
        return self.reduction(
            chunk.var_chunk, partial(chunk.var_aggregate, ddof=ddof), split_every=False
        )

    def std(self, ddof=0):
        """Standard deviation"""
        return self.var(ddof=ddof).apply(math.sqrt)

    def join(self, other, on_self, on_other=None):
        """Joins collection with another collection.

        Other collection must be one of the following:

        1.  An iterable.  We recommend tuples over lists for internal
            performance reasons.
        2.  A delayed object, pointing to a tuple.  This is recommended if the
            other collection is sizable and you're using the distributed
            scheduler.  Dask is able to pass around data wrapped in delayed
            objects with greater sophistication.
        3.  A Bag with a single partition

        You might also consider Dask Dataframe, whose join operations are much
        more heavily optimized.

        Parameters
        ----------
        other: Iterable, Delayed, Bag
            Other collection on which to join
        on_self: callable
            Function to call on elements in this collection to determine a
            match
        on_other: callable (defaults to on_self)
            Function to call on elements in the other collection to determine a
            match

        Examples
        --------
        >>> import dask.bag as db
        >>> people = db.from_sequence(['Alice', 'Bob', 'Charlie'])
        >>> fruit = ['Apple', 'Apricot', 'Banana']
        >>> list(people.join(fruit, lambda x: x[0]))
        [('Apple', 'Alice'), ('Apricot', 'Alice'), ('Banana', 'Bob')]
        """
        name = "join-" + tokenize(self, other, on_self, on_other)
        dsk = {}
        if isinstance(other, Bag):
            if other.npartitions == 1:
                dsk.update(other.dask)
                other = other.__dask_keys__()[0]
                dsk["join-%s-other" % name] = (list, other)
            else:
                msg = (
                    "Multi-bag joins are not implemented. "
                    "We recommend Dask dataframe if appropriate"
                )
                raise NotImplementedError(msg)
        elif isinstance(other, Delayed):
            dsk.update(other.dask)
            other = other._key
        elif isinstance(other, Iterable):
            other = other
        else:
            msg = (
                "Joined argument must be single-partition Bag, "
                " delayed object, or Iterable, got %s" % type(other).__name
            )
            raise TypeError(msg)

        if on_other is None:
            on_other = on_self

        for i in range(self.npartitions):
            dsk[(name, i)] = (list, (join, on_other, other, on_self, (self.name, i)))

        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return type(self)(graph, name, self.npartitions)

    def product(self, other):
        """Cartesian product between two bags."""
        assert isinstance(other, Bag)
        name = "product-" + tokenize(self, other)
        n, m = self.npartitions, other.npartitions
        dsk = {
            (name, i * m + j): (
                list,
                (itertools.product, (self.name, i), (other.name, j)),
            )
            for i in range(n)
            for j in range(m)
        }
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self, other])
        return type(self)(graph, name, n * m)

    def foldby(
        self,
        key,
        binop,
        initial=no_default,
        combine=None,
        combine_initial=no_default,
        split_every=None,
    ):
        """Combined reduction and groupby.

        Foldby provides a combined groupby and reduce for efficient parallel
        split-apply-combine tasks.

        The computation

        >>> b.foldby(key, binop, init)                        # doctest: +SKIP

        is equivalent to the following:

        >>> def reduction(group):                               # doctest: +SKIP
        ...     return reduce(binop, group, init)               # doctest: +SKIP

        >>> b.groupby(key).map(lambda (k, v): (k, reduction(v)))# doctest: +SKIP

        But uses minimal communication and so is *much* faster.

        >>> import dask.bag as db
        >>> b = db.from_sequence(range(10))
        >>> iseven = lambda x: x % 2 == 0
        >>> add = lambda x, y: x + y
        >>> dict(b.foldby(iseven, add))
        {True: 20, False: 25}

        **Key Function**

        The key function determines how to group the elements in your bag.
        In the common case where your bag holds dictionaries then the key
        function often gets out one of those elements.

        >>> def key(x):
        ...     return x['name']

        This case is so common that it is special cased, and if you provide a
        key that is not a callable function then dask.bag will turn it into one
        automatically.  The following are equivalent:

        >>> b.foldby(lambda x: x['name'], ...)  # doctest: +SKIP
        >>> b.foldby('name', ...)  # doctest: +SKIP

        **Binops**

        It can be tricky to construct the right binary operators to perform
        analytic queries.  The ``foldby`` method accepts two binary operators,
        ``binop`` and ``combine``.  Binary operators two inputs and output must
        have the same type.

        Binop takes a running total and a new element and produces a new total:

        >>> def binop(total, x):
        ...     return total + x['amount']

        Combine takes two totals and combines them:

        >>> def combine(total1, total2):
        ...     return total1 + total2

        Each of these binary operators may have a default first value for
        total, before any other value is seen.  For addition binary operators
        like above this is often ``0`` or the identity element for your
        operation.

        **split_every**

        Group partitions into groups of this size while performing reduction.
        Defaults to 8.

        >>> b.foldby('name', binop, 0, combine, 0)  # doctest: +SKIP

        Examples
        --------

        We can compute the maximum of some ``(key, value)`` pairs, grouped
        by the ``key``. (You might be better off converting the ``Bag`` to
        a ``dask.dataframe`` and using its groupby).

        >>> import random
        >>> import dask.bag as db

        >>> tokens = list('abcdefg')
        >>> values = range(10000)
        >>> a = [(random.choice(tokens), random.choice(values))
        ...       for _ in range(100)]
        >>> a[:2]  # doctest: +SKIP
        [('g', 676), ('a', 871)]

        >>> a = db.from_sequence(a)

        >>> def binop(t, x):
        ...     return max((t, x), key=lambda x: x[1])

        >>> a.foldby(lambda x: x[0], binop).compute()  # doctest: +SKIP
        [('g', ('g', 984)),
         ('a', ('a', 871)),
         ('b', ('b', 999)),
         ('c', ('c', 765)),
         ('f', ('f', 955)),
         ('e', ('e', 991)),
         ('d', ('d', 854))]

        See Also
        --------

        toolz.reduceby
        pyspark.combineByKey
        """
        if split_every is None:
            split_every = 8
        if split_every is False:
            split_every = self.npartitions

        token = tokenize(self, key, binop, initial, combine, combine_initial)
        a = "foldby-a-" + token
        if combine is None:
            combine = binop
        if initial is not no_default:
            dsk = {
                (a, i): (reduceby, key, binop, (self.name, i), initial)
                for i in range(self.npartitions)
            }
        else:
            dsk = {
                (a, i): (reduceby, key, binop, (self.name, i))
                for i in range(self.npartitions)
            }

        combine2 = partial(chunk.foldby_combine2, combine)
        depth = 0
        k = self.npartitions
        b = a
        while k > split_every:
            c = b + str(depth)
            if combine_initial is not no_default:
                for i, inds in enumerate(partition_all(split_every, range(k))):
                    dsk[(c, i)] = (
                        reduceby,
                        0,
                        combine2,
                        (toolz.concat, (map, dictitems, [(b, j) for j in inds])),
                        combine_initial,
                    )
            else:
                for i, inds in enumerate(partition_all(split_every, range(k))):
                    dsk[(c, i)] = (
                        merge_with,
                        (partial, reduce, combine),
                        [(b, j) for j in inds],
                    )

            k = i + 1
            b = c
            depth += 1

        e = "foldby-b-" + token
        if combine_initial is not no_default:
            dsk[(e, 0)] = (
                dictitems,
                (
                    reduceby,
                    0,
                    combine2,
                    (toolz.concat, (map, dictitems, [(b, j) for j in range(k)])),
                    combine_initial,
                ),
            )
        else:
            dsk[(e, 0)] = (
                dictitems,
                (merge_with, (partial, reduce, combine), [(b, j) for j in range(k)]),
            )

        graph = HighLevelGraph.from_collections(e, dsk, dependencies=[self])
        return type(self)(graph, e, 1)

    def take(self, k, npartitions=1, compute=True, warn=True):
        """Take the first k elements.

        Parameters
        ----------
        k : int
            The number of elements to return
        npartitions : int, optional
            Elements are only taken from the first ``npartitions``, with a
            default of 1. If there are fewer than ``k`` rows in the first
            ``npartitions`` a warning will be raised and any found rows
            returned. Pass -1 to use all partitions.
        compute : bool, optional
            Whether to compute the result, default is True.
        warn : bool, optional
            Whether to warn if the number of elements returned is less than
            requested, default is True.

        >>> import dask.bag as db
        >>> b = db.from_sequence(range(1_000))
        >>> b.take(3)
        (0, 1, 2)
        """

        if npartitions <= -1:
            npartitions = self.npartitions
        if npartitions > self.npartitions:
            raise ValueError(
                "only {} partitions, take "
                "received {}".format(self.npartitions, npartitions)
            )

        token = tokenize(self, k, npartitions)
        name = "take-" + token

        if npartitions > 1:
            name_p = "take-partial-" + token

            dsk = {}
            for i in range(npartitions):
                dsk[(name_p, i)] = (list, (take, k, (self.name, i)))

            concat = (toolz.concat, ([(name_p, i) for i in range(npartitions)]))
            dsk[(name, 0)] = (safe_take, k, concat, warn)
        else:
            dsk = {(name, 0): (safe_take, k, (self.name, 0), warn)}

        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        b = Bag(graph, name, 1)

        if compute:
            return tuple(b.compute())
        else:
            return b

    def flatten(self):
        """Concatenate nested lists into one long list.

        >>> import dask.bag as db
        >>> b = db.from_sequence([[1], [2, 3]])
        >>> list(b)
        [[1], [2, 3]]

        >>> list(b.flatten())
        [1, 2, 3]
        """
        name = "flatten-" + tokenize(self)
        dsk = {
            (name, i): (list, (toolz.concat, (self.name, i)))
            for i in range(self.npartitions)
        }
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return type(self)(graph, name, self.npartitions)

    def __iter__(self):
        return iter(self.compute())

    def groupby(
        self,
        grouper,
        method=None,
        npartitions=None,
        blocksize=2**20,
        max_branch=None,
        shuffle=None,
    ):
        """Group collection by key function

        This requires a full dataset read, serialization and shuffle.
        This is expensive.  If possible you should use ``foldby``.

        Parameters
        ----------
        grouper: function
            Function on which to group elements
        shuffle: str
            Either 'disk' for an on-disk shuffle or 'tasks' to use the task
            scheduling framework.  Use 'disk' if you are on a single machine
            and 'tasks' if you are on a distributed cluster.
        npartitions: int
            If using the disk-based shuffle, the number of output partitions
        blocksize: int
            If using the disk-based shuffle, the size of shuffle blocks (bytes)
        max_branch: int
            If using the task-based shuffle, the amount of splitting each
            partition undergoes.  Increase this for fewer copies but more
            scheduler overhead.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(range(10))
        >>> iseven = lambda x: x % 2 == 0
        >>> dict(b.groupby(iseven))             # doctest: +SKIP
        {True: [0, 2, 4, 6, 8], False: [1, 3, 5, 7, 9]}

        See Also
        --------
        Bag.foldby
        """
        if method is not None:
            raise Exception("The method= keyword has been moved to shuffle=")
        if shuffle is None:
            try:
                shuffle = get_default_shuffle_method()
            except ImportError:
                shuffle = "tasks"
            else:
                if shuffle == "p2p":
                    # Not implemented for Bags
                    shuffle = "tasks"
        if shuffle == "disk":
            return groupby_disk(
                self, grouper, npartitions=npartitions, blocksize=blocksize
            )
        elif shuffle == "tasks":
            return groupby_tasks(self, grouper, max_branch=max_branch)
        else:
            msg = "Shuffle must be 'disk' or 'tasks'"
            raise NotImplementedError(msg)

    def to_dataframe(self, meta=None, columns=None, optimize_graph=True):
        """Create Dask Dataframe from a Dask Bag.

        Bag should contain tuples, dict records, or scalars.

        Index will not be particularly meaningful.  Use ``reindex`` afterwards
        if necessary.

        Parameters
        ----------
        meta : pd.DataFrame, dict, iterable, optional
            An empty ``pd.DataFrame`` that matches the dtypes and column names
            of the output. This metadata is necessary for many algorithms in
            dask dataframe to work.  For ease of use, some alternative inputs
            are also available. Instead of a ``DataFrame``, a ``dict`` of
            ``{name: dtype}`` or iterable of ``(name, dtype)`` can be provided.
            If not provided or a list, a single element from the first
            partition will be computed, triggering a potentially expensive call
            to ``compute``. This may lead to unexpected results, so providing
            ``meta`` is recommended. For more information, see
            ``dask.dataframe.utils.make_meta``.
        columns : sequence, optional
            Column names to use. If the passed data do not have names
            associated with them, this argument provides names for the columns.
            Otherwise this argument indicates the order of the columns in the
            result (any names not found in the data will become all-NA
            columns).  Note that if ``meta`` is provided, column names will be
            taken from there and this parameter is invalid.
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            :class:`dask.dataframe.DataFrame`.


        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence([{'name': 'Alice',   'balance': 100},
        ...                       {'name': 'Bob',     'balance': 200},
        ...                       {'name': 'Charlie', 'balance': 300}],
        ...                      npartitions=2)
        >>> df = b.to_dataframe()

        >>> df.compute()
              name  balance
        0    Alice      100
        1      Bob      200
        0  Charlie      300
        """
        import pandas as pd

        import dask.dataframe as dd

        if meta is None:
            head = self.take(1, warn=False)
            if len(head) == 0:
                raise ValueError(
                    "`dask.bag.Bag.to_dataframe` failed to "
                    "properly infer metadata, please pass in "
                    "metadata via the `meta` keyword"
                )
            meta = pd.DataFrame(list(head), columns=columns)
        elif columns is not None:
            raise ValueError("Can't specify both `meta` and `columns`")
        else:
            meta = dd.utils.make_meta(meta, parent_meta=pd.DataFrame())
        # Serializing the columns and dtypes is much smaller than serializing
        # the empty frame
        cols = list(meta.columns)
        dtypes = meta.dtypes.to_dict()

        dfs = self.map_partitions(to_dataframe, cols, dtypes)
        if optimize_graph:
            dsk = self.__dask_optimize__(dfs.dask, dfs.__dask_keys__())
        else:
            dsk = dfs.dask

        divisions = [None] * (self.npartitions + 1)
        if not dd._dask_expr_enabled():
            return dd.DataFrame(dsk, dfs.name, meta, divisions)
        else:
            from dask_expr import from_legacy_dataframe

            from dask.dataframe.core import DataFrame

            df = DataFrame(dsk, dfs.name, meta, divisions)
            return from_legacy_dataframe(df)

    def to_delayed(self, optimize_graph=True):
        """Convert into a list of ``dask.delayed`` objects, one per partition.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.

        See Also
        --------
        dask.bag.from_delayed
        """
        from dask.delayed import Delayed

        keys = self.__dask_keys__()
        dsk = self.__dask_graph__()
        layer = self.name
        if optimize_graph:
            dsk = self.__dask_optimize__(dsk, keys)
            layer = "delayed-" + layer
            dsk = HighLevelGraph.from_collections(layer, dsk, dependencies=())
        return [Delayed(k, dsk, layer=layer) for k in keys]

    def repartition(self, npartitions=None, partition_size=None):
        """Repartition Bag across new divisions.

        Parameters
        ----------
        npartitions : int, optional
            Number of partitions of output.
        partition_size : int or string, optional
            Max number of bytes of memory for each partition. Use numbers or
            strings like 5MB.

            .. warning::

               This keyword argument triggers computation to determine
               the memory size of each partition, which may be expensive.

        Notes
        -----
        Exactly one of ``npartitions`` or ``partition_size`` should be specified.
        A ``ValueError`` will be raised when that is not the case.

        Examples
        --------
        >>> b.repartition(5)  # set to have 5 partitions  # doctest: +SKIP
        """
        if sum([partition_size is not None, npartitions is not None]) != 1:
            raise ValueError(
                "Please provide exactly one ``npartitions`` or ``partition_size`` keyword arguments"
            )
        if npartitions is not None:
            return repartition_npartitions(self, npartitions)
        elif partition_size is not None:
            return repartition_size(self, partition_size)

    def accumulate(self, binop, initial=no_default):
        """Repeatedly apply binary function to a sequence, accumulating results.

        This assumes that the bag is ordered.  While this is typically the case
        not all Dask.bag functions preserve this property.

        Examples
        --------
        >>> import dask.bag as db
        >>> from operator import add
        >>> b = db.from_sequence([1, 2, 3, 4, 5], npartitions=2)
        >>> b.accumulate(add).compute()
        [1, 3, 6, 10, 15]

        Accumulate also takes an optional argument that will be used as the
        first value.

        >>> b.accumulate(add, initial=-1).compute()
        [-1, 0, 2, 5, 9, 14]
        """
        token = tokenize(self, binop, initial)
        binop_name = funcname(binop)
        a = f"{binop_name}-part-{token}"
        b = f"{binop_name}-first-{token}"
        c = f"{binop_name}-second-{token}"
        dsk = {
            (a, 0): (accumulate_part, binop, (self.name, 0), initial, True),
            (b, 0): (first, (a, 0)),
            (c, 0): (second, (a, 0)),
        }
        for i in range(1, self.npartitions):
            dsk[(a, i)] = (accumulate_part, binop, (self.name, i), (c, i - 1))
            dsk[(b, i)] = (first, (a, i))
            dsk[(c, i)] = (second, (a, i))
        graph = HighLevelGraph.from_collections(b, dsk, dependencies=[self])
        return Bag(graph, b, self.npartitions)


def accumulate_part(binop, seq, initial, is_first=False):
    if initial is no_default:
        res = list(accumulate(binop, seq))
    else:
        res = list(accumulate(binop, seq, initial=initial))
    if is_first:
        return res, res[-1] if res else [], initial
    return res[1:], res[-1]


def partition(grouper, sequence, npartitions, p, nelements=2**20):
    """Partition a bag along a grouper, store partitions on disk."""
    for block in partition_all(nelements, sequence):
        d = groupby(grouper, block)
        d2 = defaultdict(list)
        for k, v in d.items():
            d2[abs(int(tokenize(k), 16)) % npartitions].extend(v)
        p.append(d2, fsync=True)
    return p


def collect(grouper, group, p, barrier_token):
    """Collect partitions from disk and yield k,v group pairs."""
    d = groupby(grouper, p.get(group, lock=False))
    return list(d.items())


def from_sequence(seq, partition_size=None, npartitions=None):
    """Create a dask Bag from Python sequence.

    This sequence should be relatively small in memory.  Dask Bag works
    best when it handles loading your data itself.  Commonly we load a
    sequence of filenames into a Bag and then use ``.map`` to open them.

    Parameters
    ----------
    seq: Iterable
        A sequence of elements to put into the dask
    partition_size: int (optional)
        The length of each partition
    npartitions: int (optional)
        The number of desired partitions

    It is best to provide either ``partition_size`` or ``npartitions``
    (though not both.)

    Examples
    --------
    >>> import dask.bag as db
    >>> b = db.from_sequence(['Alice', 'Bob', 'Chuck'], partition_size=2)

    See Also
    --------
    read_text: Create bag from text files
    """
    seq = list(seq)
    if npartitions and not partition_size:
        if len(seq) <= 100:
            partition_size = int(math.ceil(len(seq) / npartitions))
        else:
            partition_size = max(1, int(math.floor(len(seq) / npartitions)))
    if npartitions is None and partition_size is None:
        if len(seq) <= 100:
            partition_size = 1
        else:
            partition_size = max(1, math.ceil(math.sqrt(len(seq)) / math.sqrt(100)))

    parts = list(partition_all(partition_size, seq))
    name = "from_sequence-" + tokenize(seq, partition_size)
    if len(parts) > 0:
        d = {(name, i): list(part) for i, part in enumerate(parts)}
    else:
        d = {(name, 0): []}

    return Bag(d, name, len(d))


def from_url(urls):
    """Create a dask Bag from a url.

    Examples
    --------
    >>> a = from_url('http://raw.githubusercontent.com/dask/dask/main/README.rst')
    >>> a.npartitions
    1

    >>> a.take(8)  # doctest: +SKIP
    (b'Dask\\n',
     b'====\\n',
     b'\\n',
     b'|Build Status| |Coverage| |Doc Status| |Discourse| |Version Status| |NumFOCUS|\\n',
     b'\\n',
     b'Dask is a flexible parallel computing library for analytics.  See\\n',
     b'documentation_ for more information.\\n',
     b'\\n')

    >>> b = from_url(['http://github.com', 'http://google.com'])
    >>> b.npartitions
    2
    """
    if isinstance(urls, str):
        urls = [urls]
    name = "from_url-" + uuid.uuid4().hex
    dsk = {}
    for i, u in enumerate(urls):
        dsk[(name, i)] = (list, (urlopen, u))
    return Bag(dsk, name, len(urls))


def dictitems(d):
    """A pickleable version of dict.items

    >>> dictitems({'x': 1})
    [('x', 1)]
    """
    return list(d.items())


def concat(bags):
    """Concatenate many bags together, unioning all elements.

    >>> import dask.bag as db
    >>> a = db.from_sequence([1, 2, 3])
    >>> b = db.from_sequence([4, 5, 6])
    >>> c = db.concat([a, b])

    >>> list(c)
    [1, 2, 3, 4, 5, 6]
    """
    name = "concat-" + tokenize(*bags)
    counter = itertools.count(0)
    dsk = {(name, next(counter)): key for bag in bags for key in bag.__dask_keys__()}
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=bags)
    return Bag(graph, name, len(dsk))


def reify(seq):
    if isinstance(seq, Iterator):
        seq = list(seq)
    if len(seq) and isinstance(seq[0], Iterator):
        seq = list(map(list, seq))
    return seq


def from_delayed(values):
    """Create bag from many dask Delayed objects.

    These objects will become the partitions of the resulting Bag.  They should
    evaluate to a ``list`` or some other concrete sequence.

    Parameters
    ----------
    values: list of delayed values
        An iterable of dask Delayed objects.  Each evaluating to a list.

    Returns
    -------
    Bag

    Examples
    --------
    >>> x, y, z = [delayed(load_sequence_from_file)(fn)
    ...             for fn in filenames] # doctest: +SKIP
    >>> b = from_delayed([x, y, z])  # doctest: +SKIP

    See also
    --------
    dask.delayed
    """
    from dask.delayed import Delayed, delayed

    if isinstance(values, Delayed):
        values = [values]
    values = [
        delayed(v) if not isinstance(v, Delayed) and hasattr(v, "key") else v
        for v in values
    ]

    name = "bag-from-delayed-" + tokenize(*values)
    names = [(name, i) for i in range(len(values))]
    values2 = [(reify, v.key) for v in values]
    dsk = dict(zip(names, values2))

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=values)
    return Bag(graph, name, len(values))


def chunk_distinct(seq, key=None):
    if key is not None and not callable(key):
        key = partial(chunk.getitem, key=key)
    return list(unique(seq, key=key))


def merge_distinct(seqs, key=None):
    return chunk_distinct(toolz.concat(seqs), key=key)


def merge_frequencies(seqs):
    if isinstance(seqs, Iterable):
        seqs = list(seqs)
    if not seqs:
        return {}
    first, rest = seqs[0], seqs[1:]
    if not rest:
        return first
    out = defaultdict(int)
    out.update(first)
    for d in rest:
        for k, v in d.items():
            out[k] += v
    return out


def bag_range(n, npartitions):
    """Numbers from zero to n

    Examples
    --------

    >>> import dask.bag as db
    >>> b = db.range(5, npartitions=2)
    >>> list(b)
    [0, 1, 2, 3, 4]
    """
    size = n // npartitions
    name = "range-%d-npartitions-%d" % (n, npartitions)
    ijs = list(enumerate(take(npartitions, range(0, n, size))))
    dsk = {(name, i): (reify, (range, j, min(j + size, n))) for i, j in ijs}

    if n % npartitions != 0:
        i, j = ijs[-1]
        dsk[(name, i)] = (reify, (range, j, n))

    return Bag(dsk, name, npartitions)


def bag_zip(*bags):
    """Partition-wise bag zip

    All passed bags must have the same number of partitions.

    NOTE: corresponding partitions should have the same length; if they do not,
    the "extra" elements from the longer partition(s) will be dropped.  If you
    have this case chances are that what you really need is a data alignment
    mechanism like pandas's, and not a missing value filler like zip_longest.

    Examples
    --------

    Correct usage:

    >>> import dask.bag as db
    >>> evens = db.from_sequence(range(0, 10, 2), partition_size=4)
    >>> odds = db.from_sequence(range(1, 10, 2), partition_size=4)
    >>> pairs = db.zip(evens, odds)
    >>> list(pairs)
    [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]

    Incorrect usage:

    >>> numbers = db.range(31, npartitions=1)
    >>> fizz = numbers.filter(lambda n: n % 3 == 0)
    >>> buzz = numbers.filter(lambda n: n % 5 == 0)
    >>> fizzbuzz = db.zip(fizz, buzz)
    >>> list(fizzbuzz)
    [(0, 0), (3, 5), (6, 10), (9, 15), (12, 20), (15, 25), (18, 30)]

    When what you really wanted was more along the lines of the following:

    >>> list(fizzbuzz) # doctest: +SKIP
    (0, 0), (3, None), (None, 5), (6, None), (9, None), (None, 10),
    (12, None), (15, 15), (18, None), (None, 20),
    (21, None), (24, None), (None, 25), (27, None), (30, 30)
    """
    npartitions = bags[0].npartitions
    assert all(bag.npartitions == npartitions for bag in bags)
    # TODO: do more checks

    name = "zip-" + tokenize(*bags)
    dsk = {
        (name, i): (reify, (zip,) + tuple((bag.name, i) for bag in bags))
        for i in range(npartitions)
    }
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=bags)
    return Bag(graph, name, npartitions)


def map_chunk(f, iters, iter_kwarg_keys=None, kwargs=None):
    """Map ``f`` across one or more iterables, maybe with keyword arguments.

    Low-level function used in ``bag_map``, not user facing.

    Arguments
    ---------
    f : callable
    iters : List[Iterable]
    iter_kwarg_keys : List[str] or None
        Keyword names to use for pair with the tail end of ``iters``, allowing
        keyword arguments to be passed in from iterators.
    kwargs : dict or None
        Additional constant keyword arguments to use on every call to ``f``.
    """
    if kwargs:
        f = partial(f, **kwargs)
    iters = [iter(a) for a in iters]
    return _MapChunk(f, iters, kwarg_keys=iter_kwarg_keys)


class _MapChunk(Iterator):
    def __init__(self, f, iters, kwarg_keys=None):
        self.f = f
        self.iters = iters
        self.kwarg_keys = kwarg_keys or ()
        self.nkws = len(self.kwarg_keys)

    def __next__(self):
        try:
            vals = [next(i) for i in self.iters]
        except StopIteration:
            self.check_all_iterators_consumed()
            raise

        if self.nkws:
            args = vals[: -self.nkws]
            kwargs = dict(zip(self.kwarg_keys, vals[-self.nkws :]))
            return self.f(*args, **kwargs)
        return self.f(*vals)

    def check_all_iterators_consumed(self):
        if len(self.iters) > 1:
            for i in self.iters:
                if isinstance(i, itertools.repeat):
                    continue
                try:
                    next(i)
                except StopIteration:
                    pass
                else:
                    msg = (
                        "map called with multiple bags that aren't identically "
                        "partitioned. Please ensure that all bag arguments "
                        "have the same partition lengths"
                    )
                    raise ValueError(msg)


def starmap_chunk(f, x, kwargs):
    if kwargs:
        f = partial(f, **kwargs)
    return itertools.starmap(f, x)


def unpack_scalar_dask_kwargs(kwargs):
    """Extracts dask values from kwargs.

    Currently only ``dask.bag.Item`` and ``dask.delayed.Delayed`` are
    supported.  Returns a merged dask graph and a task resulting in a keyword
    dict.
    """
    kwargs2 = {}
    dependencies = []
    for k, v in kwargs.items():
        vv, collections = unpack_collections(v)
        if not collections:
            kwargs2[k] = v
        else:
            kwargs2[k] = vv
            dependencies.extend(collections)
    if dependencies:
        kwargs2 = (dict, (zip, list(kwargs2), list(kwargs2.values())))
    return kwargs2, dependencies


def bag_map(func, *args, **kwargs):
    """Apply a function elementwise across one or more bags.

    Note that all ``Bag`` arguments must be partitioned identically.

    Parameters
    ----------
    func : callable
    *args, **kwargs : Bag, Item, Delayed, or object
        Arguments and keyword arguments to pass to ``func``. Non-Bag args/kwargs
        are broadcasted across all calls to ``func``.

    Notes
    -----
    For calls with multiple `Bag` arguments, corresponding partitions should
    have the same length; if they do not, the call will error at compute time.

    Examples
    --------
    >>> import dask.bag as db
    >>> b = db.from_sequence(range(5), npartitions=2)
    >>> b2 = db.from_sequence(range(5, 10), npartitions=2)

    Apply a function to all elements in a bag:

    >>> db.map(lambda x: x + 1, b).compute()
    [1, 2, 3, 4, 5]

    Apply a function with arguments from multiple bags:

    >>> from operator import add
    >>> db.map(add, b, b2).compute()
    [5, 7, 9, 11, 13]

    Non-bag arguments are broadcast across all calls to the mapped function:

    >>> db.map(add, b, 1).compute()
    [1, 2, 3, 4, 5]

    Keyword arguments are also supported, and have the same semantics as
    regular arguments:

    >>> def myadd(x, y=0):
    ...     return x + y
    >>> db.map(myadd, b, y=b2).compute()
    [5, 7, 9, 11, 13]
    >>> db.map(myadd, b, y=1).compute()
    [1, 2, 3, 4, 5]

    Both arguments and keyword arguments can also be instances of
    ``dask.bag.Item`` or ``dask.delayed.Delayed``. Here we'll add the max value
    in the bag to each element:

    >>> db.map(myadd, b, b.max()).compute()
    [4, 5, 6, 7, 8]
    """
    name = "{}-{}".format(funcname(func), tokenize(func, "map", *args, **kwargs))
    dependencies = []

    bags = []
    args2 = []
    for a in args:
        if isinstance(a, Bag):
            bags.append(a)
            args2.append(a)
        elif isinstance(a, (Item, Delayed)):
            dependencies.append(a)
            args2.append((itertools.repeat, a.key))
        else:
            args2.append((itertools.repeat, a))

    bag_kwargs = {}
    other_kwargs = {}
    for k, v in kwargs.items():
        if isinstance(v, Bag):
            bag_kwargs[k] = v
            bags.append(v)
        else:
            other_kwargs[k] = v

    other_kwargs, collections = unpack_scalar_dask_kwargs(other_kwargs)
    dependencies.extend(collections)

    if not bags:
        raise ValueError("At least one argument must be a Bag.")

    npartitions = {b.npartitions for b in bags}
    if len(npartitions) > 1:
        raise ValueError("All bags must have the same number of partitions.")
    npartitions = npartitions.pop()

    def build_iters(n):
        args = [(a.name, n) if isinstance(a, Bag) else a for a in args2]
        if bag_kwargs:
            args.extend((b.name, n) for b in bag_kwargs.values())
        return args

    if bag_kwargs:
        iter_kwarg_keys = list(bag_kwargs)
    else:
        iter_kwarg_keys = None

    dsk = {
        (name, n): (
            reify,
            (map_chunk, func, build_iters(n), iter_kwarg_keys, other_kwargs),
        )
        for n in range(npartitions)
    }

    # If all bags are the same type, use that type, otherwise fallback to Bag
    return_type = set(map(type, bags))
    return_type = return_type.pop() if len(return_type) == 1 else Bag

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=bags + dependencies)

    return return_type(graph, name, npartitions)


def map_partitions(func, *args, **kwargs):
    """Apply a function to every partition across one or more bags.

    Note that all ``Bag`` arguments must be partitioned identically.

    Parameters
    ----------
    func : callable
    *args, **kwargs : Bag, Item, Delayed, or object
        Arguments and keyword arguments to pass to ``func``.

    Examples
    --------
    >>> import dask.bag as db
    >>> b = db.from_sequence(range(1, 101), npartitions=10)
    >>> def div(nums, den=1):
    ...     return [num / den for num in nums]

    Using a python object:

    >>> hi = b.max().compute()
    >>> hi
    100
    >>> b.map_partitions(div, den=hi).take(5)
    (0.01, 0.02, 0.03, 0.04, 0.05)

    Using an ``Item``:

    >>> b.map_partitions(div, den=b.max()).take(5)
    (0.01, 0.02, 0.03, 0.04, 0.05)

    Note that while both versions give the same output, the second forms a
    single graph, and then computes everything at once, and in some cases
    may be more efficient.
    """
    name = kwargs.pop("token", None) or funcname(func)
    name = "{}-{}".format(name, tokenize(func, "map-partitions", *args, **kwargs))
    bags = []
    args2 = []
    dependencies = []
    for a in args:
        if isinstance(a, Bag):
            bags.append(a)
            args2.append(a)
        elif isinstance(a, (Item, Delayed)):
            args2.append(a.key)
            dependencies.append(a)
        else:
            args2.append(a)

    bag_kwargs = {}
    other_kwargs = {}
    for k, v in kwargs.items():
        if isinstance(v, Bag):
            bag_kwargs[k] = v
            bags.append(v)
        else:
            other_kwargs[k] = v

    other_kwargs, collections = unpack_scalar_dask_kwargs(other_kwargs)
    dependencies.extend(collections)

    if not bags:
        raise ValueError("At least one argument must be a Bag.")

    npartitions = {b.npartitions for b in bags}
    if len(npartitions) > 1:
        raise ValueError("All bags must have the same number of partitions.")
    npartitions = npartitions.pop()

    def build_args(n):
        return [(a.name, n) if isinstance(a, Bag) else a for a in args2]

    def build_bag_kwargs(n):
        if not bag_kwargs:
            return {}
        return (
            dict,
            (zip, list(bag_kwargs), [(b.name, n) for b in bag_kwargs.values()]),
        )

    if bag_kwargs:
        # Avoid using `blockwise` when a key-word
        # argument is being used to refer to a collection.
        dsk = {
            (name, n): (
                apply,
                func,
                build_args(n),
                (merge, build_bag_kwargs(n), other_kwargs),
            )
            for n in range(npartitions)
        }
    else:
        pairs = []
        numblocks = {}
        for arg in args2:
            # TODO: Allow interaction with Frame/Array
            # collections with proper partitioning
            if isinstance(arg, Bag):
                pairs.extend([arg.name, "i"])
                numblocks[arg.name] = (arg.npartitions,)
            else:
                pairs.extend([arg, None])
        if other_kwargs and isinstance(other_kwargs, tuple):
            # `other_kwargs` is a nested subgraph,
            # designed to generate the kwargs lazily.
            # We need to convert this to a dictionary
            # before passing to `blockwise`
            other_kwargs = other_kwargs[0](other_kwargs[1][0](*other_kwargs[1][1:]))
        dsk = blockwise(
            func,
            name,
            "i",
            *pairs,
            numblocks=numblocks,
            concatenate=True,
            dependencies=dependencies,
            **other_kwargs,
        )

    # If all bags are the same type, use that type, otherwise fallback to Bag
    return_type = set(map(type, bags))
    return_type = return_type.pop() if len(return_type) == 1 else Bag

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=bags + dependencies)

    return return_type(graph, name, npartitions)


def _reduce(binop, sequence, initial=no_default):
    if initial is not no_default:
        return reduce(binop, sequence, initial)
    else:
        return reduce(binop, sequence)


def make_group(k, stage):
    def h(x):
        return x[0] // k**stage % k

    return h


def groupby_tasks(b, grouper, hash=lambda x: int(tokenize(x), 16), max_branch=32):
    max_branch = max_branch or 32
    n = b.npartitions

    stages = int(math.ceil(math.log(n) / math.log(max_branch))) or 1
    if stages > 1:
        k = int(math.ceil(n ** (1 / stages)))
    else:
        k = n

    groups = []
    splits = []
    joins = []

    inputs = [tuple(digit(i, j, k) for j in range(stages)) for i in range(k**stages)]

    b2 = b.map(partial(chunk.groupby_tasks_group_hash, hash=hash, grouper=grouper))

    token = tokenize(b, grouper, hash, max_branch)

    shuffle_join_name = "shuffle-join-" + token
    shuffle_group_name = "shuffle-group-" + token
    shuffle_split_name = "shuffle-split-" + token

    start = {}

    for idx, inp in enumerate(inputs):
        group = {}
        split = {}
        if idx < b.npartitions:
            start[(shuffle_join_name, 0, inp)] = (b2.name, idx)
        else:
            start[(shuffle_join_name, 0, inp)] = []

        for stage in range(1, stages + 1):
            _key_tuple = (shuffle_group_name, stage, inp)
            group[_key_tuple] = (
                groupby,
                (make_group, k, stage - 1),
                (shuffle_join_name, stage - 1, inp),
            )

            for i in range(k):
                split[(shuffle_split_name, stage, i, inp)] = (
                    dict.get,
                    _key_tuple,
                    i,
                    {},
                )

        groups.append(group)
        splits.append(split)

    for stage in range(1, stages + 1):
        join = {
            (shuffle_join_name, stage, inp): (
                list,
                (
                    toolz.concat,
                    [
                        (
                            shuffle_split_name,
                            stage,
                            inp[stage - 1],
                            insert(inp, stage - 1, j),
                        )
                        for j in range(k)
                    ],
                ),
            )
            for inp in inputs
        }

        joins.append(join)

    name = "shuffle-" + token

    end = {
        (name, i): (list, (dict.items, (groupby, grouper, (pluck, 1, j))))
        for i, j in enumerate(join)
    }

    groups.extend(splits)
    groups.extend(joins)

    dsk = merge(start, end, *(groups))
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[b2])
    return type(b)(graph, name, len(inputs))


def groupby_disk(b, grouper, npartitions=None, blocksize=2**20):
    if npartitions is None:
        npartitions = b.npartitions
    token = tokenize(b, grouper, npartitions, blocksize)

    import partd

    p = ("partd-" + token,)
    dirname = config.get("temporary_directory", None)
    if dirname:
        file = (apply, partd.File, (), {"dir": dirname})
    else:
        file = (partd.File,)
    try:
        dsk1 = {p: (partd.Python, (partd.Snappy, file))}
    except AttributeError:
        dsk1 = {p: (partd.Python, file)}

    # Partition data on disk
    name = f"groupby-part-{funcname(grouper)}-{token}"
    dsk2 = {
        (name, i): (partition, grouper, (b.name, i), npartitions, p, blocksize)
        for i in range(b.npartitions)
    }

    # Barrier
    barrier_token = "groupby-barrier-" + token

    dsk3 = {barrier_token: (chunk.barrier,) + tuple(dsk2)}

    # Collect groups
    name = "groupby-collect-" + token
    dsk4 = {
        (name, i): (collect, grouper, i, p, barrier_token) for i in range(npartitions)
    }

    dsk = merge(dsk1, dsk2, dsk3, dsk4)
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[b])
    return type(b)(graph, name, npartitions)


def empty_safe_apply(func, part, is_last):
    if isinstance(part, Iterator):
        try:
            _, part = peek(part)
        except StopIteration:
            if not is_last:
                return no_result
        return func(part)
    elif not is_last and len(part) == 0:
        return no_result
    else:
        return func(part)


def empty_safe_aggregate(func, parts, is_last):
    parts2 = (p for p in parts if p is not no_result)
    return empty_safe_apply(func, parts2, is_last)


def safe_take(n, b, warn=True):
    r = list(take(n, b))
    if len(r) != n and warn:
        warnings.warn(
            f"Insufficient elements for `take`. {n} elements requested, only {len(r)} "
            "elements available. Try passing larger `npartitions` to `take`."
        )
    return r


def random_sample(x, state_data, prob):
    """Filter elements of `x` by a probability `prob`.

    Parameters
    ----------
    x : iterable
    state_data : tuple
        A tuple that can be passed to ``random.Random.setstate``.
    prob : float
        A float between 0 and 1, representing the probability that each
        element will be yielded.
    """
    random_state = Random()
    random_state.setstate(state_data)
    for i in x:
        if random_state.random() < prob:
            yield i


def random_state_data_python(
    n: int, random_state: int | Random | None = None
) -> list[tuple[int, tuple[int, ...], None]]:
    """Return a list of tuples that can be passed to
    ``random.Random.setstate``.

    Parameters
    ----------
    n : int
        Number of tuples to return.
    random_state : int or ``random.Random``, optional
        If an int, is used to seed a new ``random.Random``.

    See Also
    --------
    dask.utils.random_state_data
    """
    maxuint32 = 1 << 32

    try:
        import numpy as np

        if isinstance(random_state, Random):
            random_state = random_state.randint(0, maxuint32)
        np_rng = np.random.default_rng(random_state)

        random_data = np_rng.bytes(624 * n * 4)  # `n * 624` 32-bit integers
        arr = np.frombuffer(random_data, dtype=np.uint32).reshape((n, -1))
        return [(3, tuple(row) + (624,), None) for row in arr.tolist()]

    except ImportError:
        # Pure python (much slower)
        if not isinstance(random_state, Random):
            random_state = Random(random_state)

        return [
            (
                3,
                tuple(random_state.randint(0, maxuint32) for _ in range(624)) + (624,),
                None,
            )
            for _ in range(n)
        ]


def split(seq, n):
    """Split apart a sequence into n equal pieces.

    >>> split(range(10), 3)
    [[0, 1, 2], [3, 4, 5], [6, 7, 8, 9]]
    """
    if not isinstance(seq, (list, tuple)):
        seq = list(seq)

    part = len(seq) / n
    L = [seq[int(part * i) : int(part * (i + 1))] for i in range(n - 1)]
    L.append(seq[int(part * (n - 1)) :])
    return L


def to_dataframe(seq, columns, dtypes):
    import pandas as pd
    from packaging.version import Version

    seq = reify(seq)
    # pd.DataFrame expects lists, only copy if necessary
    if not isinstance(seq, list):
        seq = list(seq)

    kwargs = {} if Version(pd.__version__).major >= 3 else {"copy": False}
    res = pd.DataFrame(seq, columns=list(columns))
    return res.astype(dtypes, **kwargs)


def repartition_npartitions(bag, npartitions):
    """Changes the number of partitions of the bag.

    This can be used to reduce or increase the number of partitions
    of the bag.
    """
    if npartitions == bag.npartitions:
        return bag

    new_name = "repartition-%d-%s" % (npartitions, tokenize(bag, npartitions))
    if bag.npartitions > npartitions:
        ratio = bag.npartitions / npartitions
        new_partitions_boundaries = [
            int(old_partition_index * ratio)
            for old_partition_index in range(npartitions + 1)
        ]
        return _repartition_from_boundaries(bag, new_partitions_boundaries, new_name)
    else:  # npartitions > bag.npartitions
        div, mod = divmod(npartitions, bag.npartitions)
        nsplits = [div] * bag.npartitions
        nsplits[-1] += mod
        return _split_partitions(bag, nsplits, new_name)


def total_mem_usage(partition):
    from copy import deepcopy

    # if repartition is called multiple times prior to calling compute(), the partitions
    # will be an Iterable. Copy the object to avoid consuming the iterable.
    if isinstance(partition, Iterable):
        partition = reify(deepcopy(partition))
    return sizeof(partition)


def repartition_size(bag, size):
    """
    Repartition bag so that new partitions have approximately `size` memory usage each
    """
    if isinstance(size, str):
        size = parse_bytes(size)
    size = int(size)
    mem_usages = bag.map_partitions(total_mem_usage).compute()

    # 1. split each partition that is larger than partition size
    nsplits = [1 + mem_usage // size for mem_usage in mem_usages]
    if any(nsplit > 1 for nsplit in nsplits):
        split_name = f"repartition-split-{tokenize(bag, size)}"
        bag = _split_partitions(bag, nsplits, split_name)
        # update mem_usages to account for the split partitions
        split_mem_usages = []
        for n, usage in zip(nsplits, mem_usages):
            split_mem_usages.extend([usage / n] * n)
        mem_usages = split_mem_usages

    # 2. now that all partitions are less than size, concat them up to size
    assert all(mem_usage <= size for mem_usage in mem_usages)
    new_npartitions = list(map(len, iter_chunks(mem_usages, size)))
    new_partitions_boundaries = accumulate(operator.add, new_npartitions)
    new_name = f"repartition-{tokenize(bag, size)}"
    return _repartition_from_boundaries(bag, new_partitions_boundaries, new_name)


def _split_partitions(bag, nsplits, new_name):
    """Split a Dask bag into new partitions

    Parameters
    ----------
    bag: Dask bag
    nsplits: List[int]
        Number of target bags for each partition
        The length of nsplits should be the same as bag.npartitions
    new_name: str

    See Also
    --------
    repartition_npartitions
    repartition_size
    """
    if len(nsplits) != bag.npartitions:
        raise ValueError(f"nsplits should have len={bag.npartitions}")
    dsk = {}
    split_name = f"split-{tokenize(bag, nsplits)}"
    j = 0
    for i, k in enumerate(nsplits):
        if k == 1:
            dsk[new_name, j] = (bag.name, i)
            j += 1
        else:
            dsk[split_name, i] = (split, (bag.name, i), k)
            for jj in range(k):
                dsk[new_name, j] = (operator.getitem, (split_name, i), jj)
                j += 1

    graph = HighLevelGraph.from_collections(new_name, dsk, dependencies=[bag])
    return Bag(graph, name=new_name, npartitions=sum(nsplits))


def _repartition_from_boundaries(bag, new_partitions_boundaries, new_name):
    if not isinstance(new_partitions_boundaries, list):
        new_partitions_boundaries = list(new_partitions_boundaries)
    if new_partitions_boundaries[0] > 0:
        new_partitions_boundaries.insert(0, 0)
    if new_partitions_boundaries[-1] < bag.npartitions:
        new_partitions_boundaries.append(bag.npartitions)
    num_new_partitions = len(new_partitions_boundaries) - 1
    dsk = {}
    for new_partition_index in range(num_new_partitions):
        value = (
            list,
            (
                toolz.concat,
                [
                    (bag.name, old_partition_index)
                    for old_partition_index in range(
                        new_partitions_boundaries[new_partition_index],
                        new_partitions_boundaries[new_partition_index + 1],
                    )
                ],
            ),
        )
        dsk[new_name, new_partition_index] = value
    graph = HighLevelGraph.from_collections(new_name, dsk, dependencies=[bag])
    return Bag(graph, name=new_name, npartitions=num_new_partitions)
