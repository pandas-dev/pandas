from __future__ import annotations

import dataclasses
import inspect
import uuid
import warnings
from collections import OrderedDict
from collections.abc import Hashable, Iterator, Mapping
from concurrent.futures import Executor
from contextlib import contextmanager, suppress
from contextvars import ContextVar
from functools import partial
from numbers import Integral, Number
from operator import getitem
from typing import Any, Literal, TypeVar

from tlz import groupby, merge

from dask import config, local
from dask._compatibility import EMSCRIPTEN
from dask._task_spec import DataNode, Dict, List, Task, TaskRef
from dask.core import flatten
from dask.core import get as simple_get
from dask.system import CPU_COUNT
from dask.typing import Key, SchedulerGetCallable
from dask.utils import ensure_dict, is_namedtuple_instance, key_split, shorten_traceback

_DistributedClient = None
_get_distributed_client = None
_DISTRIBUTED_AVAILABLE = None


def _distributed_available() -> bool:
    # Lazy import in get_scheduler can be expensive
    global _DistributedClient, _get_distributed_client, _DISTRIBUTED_AVAILABLE
    if _DISTRIBUTED_AVAILABLE is not None:
        return _DISTRIBUTED_AVAILABLE  # type: ignore[unreachable]
    try:
        from distributed import Client as _DistributedClient
        from distributed.worker import get_client as _get_distributed_client

        _DISTRIBUTED_AVAILABLE = True
    except ImportError:
        _DISTRIBUTED_AVAILABLE = False
    return _DISTRIBUTED_AVAILABLE


__all__ = (
    "DaskMethodsMixin",
    "annotate",
    "get_annotations",
    "is_dask_collection",
    "compute",
    "persist",
    "optimize",
    "visualize",
    "tokenize",
    "normalize_token",
    "get_collection_names",
    "get_name_from_key",
    "replace_name_in_key",
    "clone_key",
)

# Backwards compat
from dask.tokenize import TokenizationError, normalize_token, tokenize  # noqa: F401

_annotations: ContextVar[dict[str, Any] | None] = ContextVar(
    "annotations", default=None
)


def get_annotations() -> dict[str, Any]:
    """Get current annotations.

    Returns
    -------
    Dict of all current annotations

    See Also
    --------
    annotate
    """
    return _annotations.get() or {}


@contextmanager
def annotate(**annotations: Any) -> Iterator[None]:
    """Context Manager for setting HighLevelGraph Layer annotations.

    Annotations are metadata or soft constraints associated with
    tasks that dask schedulers may choose to respect: They signal intent
    without enforcing hard constraints. As such, they are
    primarily designed for use with the distributed scheduler.

    Almost any object can serve as an annotation, but small Python objects
    are preferred, while large objects such as NumPy arrays are discouraged.

    Callables supplied as an annotation should take a single *key* argument and
    produce the appropriate annotation. Individual task keys in the annotated collection
    are supplied to the callable.

    Parameters
    ----------
    **annotations : key-value pairs

    Examples
    --------

    All tasks within array A should have priority 100 and be retried 3 times
    on failure.

    >>> import dask
    >>> import dask.array as da
    >>> with dask.annotate(priority=100, retries=3):
    ...     A = da.ones((10000, 10000))

    Prioritise tasks within Array A on flattened block ID.

    >>> nblocks = (10, 10)
    >>> with dask.annotate(priority=lambda k: k[1]*nblocks[1] + k[2]):
    ...     A = da.ones((1000, 1000), chunks=(100, 100))

    Annotations may be nested.

    >>> with dask.annotate(priority=1):
    ...     with dask.annotate(retries=3):
    ...         A = da.ones((1000, 1000))
    ...     B = A + 1

    See Also
    --------
    get_annotations
    """

    # Sanity check annotations used in place of
    # legacy distributed Client.{submit, persist, compute} keywords
    if "workers" in annotations:
        if isinstance(annotations["workers"], (list, set, tuple)):
            annotations["workers"] = list(annotations["workers"])
        elif isinstance(annotations["workers"], str):
            annotations["workers"] = [annotations["workers"]]
        elif callable(annotations["workers"]):
            pass
        else:
            raise TypeError(
                "'workers' annotation must be a sequence of str, a str or a callable, but got %s."
                % annotations["workers"]
            )

    if (
        "priority" in annotations
        and not isinstance(annotations["priority"], Number)
        and not callable(annotations["priority"])
    ):
        raise TypeError(
            "'priority' annotation must be a Number or a callable, but got %s"
            % annotations["priority"]
        )

    if (
        "retries" in annotations
        and not isinstance(annotations["retries"], Number)
        and not callable(annotations["retries"])
    ):
        raise TypeError(
            "'retries' annotation must be a Number or a callable, but got %s"
            % annotations["retries"]
        )

    if (
        "resources" in annotations
        and not isinstance(annotations["resources"], dict)
        and not callable(annotations["resources"])
    ):
        raise TypeError(
            "'resources' annotation must be a dict, but got %s"
            % annotations["resources"]
        )

    if (
        "allow_other_workers" in annotations
        and not isinstance(annotations["allow_other_workers"], bool)
        and not callable(annotations["allow_other_workers"])
    ):
        raise TypeError(
            "'allow_other_workers' annotations must be a bool or a callable, but got %s"
            % annotations["allow_other_workers"]
        )
    ctx_annot = _annotations.get()
    if ctx_annot is None:
        ctx_annot = {}
    token = _annotations.set(merge(ctx_annot, annotations))
    try:
        yield
    finally:
        _annotations.reset(token)


def is_dask_collection(x) -> bool:
    """Returns ``True`` if ``x`` is a dask collection.

    Parameters
    ----------
    x : Any
        Object to test.

    Returns
    -------
    result : bool
        ``True`` if `x` is a Dask collection.

    Notes
    -----
    The DaskCollection typing.Protocol implementation defines a Dask
    collection as a class that returns a Mapping from the
    ``__dask_graph__`` method. This helper function existed before the
    implementation of the protocol.

    """
    if (
        isinstance(x, type)
        or not hasattr(x, "__dask_graph__")
        or not callable(x.__dask_graph__)
    ):
        return False

    pkg_name = getattr(type(x), "__module__", "")
    if pkg_name.split(".")[0] in ("dask_cudf",):
        # Temporary hack to avoid graph materialization. Note that this won't work with
        # dask_expr.array objects wrapped by xarray or pint. By the time dask_expr.array
        # is published, we hope to be able to rewrite this method completely.
        # Read: https://github.com/dask/dask/pull/10676
        return True
    elif pkg_name.startswith("dask.dataframe.dask_expr"):
        return True
    elif pkg_name.startswith("dask.array._array_expr"):
        return True

    # xarray, pint, and possibly other wrappers always define a __dask_graph__ method,
    # but it may return None if they wrap around a non-dask object.
    # In all known dask collections other than dask-expr,
    # calling __dask_graph__ is cheap.
    return x.__dask_graph__() is not None


class DaskMethodsMixin:
    """A mixin adding standard dask collection methods"""

    __slots__ = ("__weakref__",)

    def visualize(self, filename="mydask", format=None, optimize_graph=False, **kwargs):
        """Render the computation of this object's task graph using graphviz.

        Requires ``graphviz`` to be installed.

        Parameters
        ----------
        filename : str or None, optional
            The name of the file to write to disk. If the provided `filename`
            doesn't include an extension, '.png' will be used by default.
            If `filename` is None, no file will be written, and we communicate
            with dot using only pipes.
        format : {'png', 'pdf', 'dot', 'svg', 'jpeg', 'jpg'}, optional
            Format in which to write output file.  Default is 'png'.
        optimize_graph : bool, optional
            If True, the graph is optimized before rendering.  Otherwise,
            the graph is displayed as is. Default is False.
        color: {None, 'order'}, optional
            Options to color nodes.  Provide ``cmap=`` keyword for additional
            colormap
        **kwargs
           Additional keyword arguments to forward to ``to_graphviz``.

        Examples
        --------
        >>> x.visualize(filename='dask.pdf')  # doctest: +SKIP
        >>> x.visualize(filename='dask.pdf', color='order')  # doctest: +SKIP

        Returns
        -------
        result : IPython.display.Image, IPython.display.SVG, or None
            See dask.dot.dot_graph for more information.

        See Also
        --------
        dask.visualize
        dask.dot.dot_graph

        Notes
        -----
        For more information on optimization see here:

        https://docs.dask.org/en/latest/optimize.html
        """
        return visualize(
            self,
            filename=filename,
            format=format,
            optimize_graph=optimize_graph,
            **kwargs,
        )

    def persist(self, **kwargs):
        """Persist this dask collection into memory

        This turns a lazy Dask collection into a Dask collection with the same
        metadata, but now with the results fully computed or actively computing
        in the background.

        The action of function differs significantly depending on the active
        task scheduler.  If the task scheduler supports asynchronous computing,
        such as is the case of the dask.distributed scheduler, then persist
        will return *immediately* and the return value's task graph will
        contain Dask Future objects.  However if the task scheduler only
        supports blocking computation then the call to persist will *block*
        and the return value's task graph will contain concrete Python results.

        This function is particularly useful when using distributed systems,
        because the results will be kept in distributed memory, rather than
        returned to the local process as with compute.

        Parameters
        ----------
        scheduler : string, optional
            Which scheduler to use like "threads", "synchronous" or "processes".
            If not provided, the default is to check the global settings first,
            and then fall back to the collection defaults.
        optimize_graph : bool, optional
            If True [default], the graph is optimized before computation.
            Otherwise the graph is run as is. This can be useful for debugging.
        **kwargs
            Extra keywords to forward to the scheduler function.

        Returns
        -------
        New dask collections backed by in-memory data

        See Also
        --------
        dask.persist
        """
        (result,) = persist(self, traverse=False, **kwargs)
        return result

    def compute(self, **kwargs):
        """Compute this dask collection

        This turns a lazy Dask collection into its in-memory equivalent.
        For example a Dask array turns into a NumPy array and a Dask dataframe
        turns into a Pandas dataframe.  The entire dataset must fit into memory
        before calling this operation.

        Parameters
        ----------
        scheduler : string, optional
            Which scheduler to use like "threads", "synchronous" or "processes".
            If not provided, the default is to check the global settings first,
            and then fall back to the collection defaults.
        optimize_graph : bool, optional
            If True [default], the graph is optimized before computation.
            Otherwise the graph is run as is. This can be useful for debugging.
        kwargs
            Extra keywords to forward to the scheduler function.

        See Also
        --------
        dask.compute
        """
        (result,) = compute(self, traverse=False, **kwargs)
        return result

    def __await__(self):
        try:
            from distributed import futures_of, wait
        except ImportError as e:
            raise ImportError(
                "Using async/await with dask requires the `distributed` package"
            ) from e

        async def f():
            if futures_of(self):
                await wait(self)
            return self

        return f().__await__()


def compute_as_if_collection(cls, dsk, keys, scheduler=None, get=None, **kwargs):
    """Compute a graph as if it were of type cls.

    Allows for applying the same optimizations and default scheduler."""
    schedule = get_scheduler(scheduler=scheduler, cls=cls, get=get)
    dsk2 = optimization_function(cls)(dsk, keys, **kwargs)
    return schedule(dsk2, keys, **kwargs)


def dont_optimize(dsk, keys, **kwargs):
    return dsk


def optimization_function(x):
    return getattr(x, "__dask_optimize__", dont_optimize)


def collections_to_dsk(collections, optimize_graph=True, optimizations=(), **kwargs):
    """
    Convert many collections into a single dask graph, after optimization
    """
    from dask.highlevelgraph import HighLevelGraph

    optimizations = tuple(optimizations) + tuple(config.get("optimizations", ()))

    if optimize_graph:
        groups = groupby(optimization_function, collections)

        graphs = []
        for opt, val in groups.items():
            dsk, keys = _extract_graph_and_keys(val)
            dsk = opt(dsk, keys, **kwargs)

            for opt_inner in optimizations:
                dsk = opt_inner(dsk, keys, **kwargs)

            graphs.append(dsk)

        # Merge all graphs
        if any(isinstance(graph, HighLevelGraph) for graph in graphs):
            dsk = HighLevelGraph.merge(*graphs)
        else:
            dsk = merge(*map(ensure_dict, graphs))
    else:
        dsk, _ = _extract_graph_and_keys(collections)

    return dsk


def _extract_graph_and_keys(vals):
    """Given a list of dask vals, return a single graph and a list of keys such
    that ``get(dsk, keys)`` is equivalent to ``[v.compute() for v in vals]``."""
    from dask.highlevelgraph import HighLevelGraph

    graphs, keys = [], []
    for v in vals:
        graphs.append(v.__dask_graph__())
        keys.append(v.__dask_keys__())

    if any(isinstance(graph, HighLevelGraph) for graph in graphs):
        graph = HighLevelGraph.merge(*graphs)
    else:
        graph = merge(*map(ensure_dict, graphs))

    return graph, keys


def unpack_collections(*args, traverse=True):
    """Extract collections in preparation for compute/persist/etc...

    Intended use is to find all collections in a set of (possibly nested)
    python objects, do something to them (compute, etc...), then repackage them
    in equivalent python objects.

    Parameters
    ----------
    *args
        Any number of objects. If it is a dask collection, it's extracted and
        added to the list of collections returned. By default, python builtin
        collections are also traversed to look for dask collections (for more
        information see the ``traverse`` keyword).
    traverse : bool, optional
        If True (default), builtin python collections are traversed looking for
        any dask collections they might contain.

    Returns
    -------
    collections : list
        A list of all dask collections contained in ``args``
    repack : callable
        A function to call on the transformed collections to repackage them as
        they were in the original ``args``.
    """

    collections = []
    repack_dsk = {}

    collections_token = uuid.uuid4().hex

    def _unpack(expr):
        if is_dask_collection(expr):
            tok = tokenize(expr)
            if tok not in repack_dsk:
                repack_dsk[tok] = Task(
                    tok, getitem, TaskRef(collections_token), len(collections)
                )
                collections.append(expr)
            return TaskRef(tok)

        tok = uuid.uuid4().hex
        tsk: DataNode | Task
        if not traverse:
            tsk = DataNode(None, expr)
        else:
            # Treat iterators like lists
            typ = list if isinstance(expr, Iterator) else type(expr)
            if typ in (list, tuple, set):
                tsk = Task(tok, typ, List(*[_unpack(i) for i in expr]))
            elif typ in (dict, OrderedDict):
                tsk = Task(
                    tok, typ, Dict({_unpack(k): _unpack(v) for k, v in expr.items()})
                )
            elif dataclasses.is_dataclass(expr) and not isinstance(expr, type):
                tsk = Task(
                    tok,
                    typ,
                    *[_unpack(getattr(expr, f.name)) for f in dataclasses.fields(expr)],
                )
            elif is_namedtuple_instance(expr):
                tsk = Task(tok, typ, *[_unpack(i) for i in expr])
            else:
                return expr

        repack_dsk[tok] = tsk
        return TaskRef(tok)

    out = uuid.uuid4().hex
    repack_dsk[out] = Task(out, tuple, List(*[_unpack(i) for i in args]))

    def repack(results):
        dsk = repack_dsk.copy()
        dsk[collections_token] = DataNode(collections_token, results)
        return simple_get(dsk, out)

    # The original `collections` is kept alive by the closure
    # This causes the collection to be only freed by the garbage collector
    collections2 = list(collections)
    collections.clear()
    return collections2, repack


def optimize(*args, traverse=True, **kwargs):
    """Optimize several dask collections at once.

    Returns equivalent dask collections that all share the same merged and
    optimized underlying graph. This can be useful if converting multiple
    collections to delayed objects, or to manually apply the optimizations at
    strategic points.

    Note that in most cases you shouldn't need to call this method directly.

    Parameters
    ----------
    *args : objects
        Any number of objects. If a dask object, its graph is optimized and
        merged with all those of all other dask objects before returning an
        equivalent dask collection. Non-dask arguments are passed through
        unchanged.
    traverse : bool, optional
        By default dask traverses builtin python collections looking for dask
        objects passed to ``optimize``. For large collections this can be
        expensive. If none of the arguments contain any dask objects, set
        ``traverse=False`` to avoid doing this traversal.
    optimizations : list of callables, optional
        Additional optimization passes to perform.
    **kwargs
        Extra keyword arguments to forward to the optimization passes.

    Examples
    --------
    >>> import dask
    >>> import dask.array as da
    >>> a = da.arange(10, chunks=2).sum()
    >>> b = da.arange(10, chunks=2).mean()
    >>> a2, b2 = dask.optimize(a, b)

    >>> a2.compute() == a.compute()
    np.True_
    >>> b2.compute() == b.compute()
    np.True_
    """
    collections, repack = unpack_collections(*args, traverse=traverse)
    if not collections:
        return args

    dsk = collections_to_dsk(collections, **kwargs)

    postpersists = []
    for a in collections:
        r, s = a.__dask_postpersist__()
        postpersists.append(r(dsk, *s))

    return repack(postpersists)


def compute(
    *args, traverse=True, optimize_graph=True, scheduler=None, get=None, **kwargs
):
    """Compute several dask collections at once.

    Parameters
    ----------
    args : object
        Any number of objects. If it is a dask object, it's computed and the
        result is returned. By default, python builtin collections are also
        traversed to look for dask objects (for more information see the
        ``traverse`` keyword). Non-dask arguments are passed through unchanged.
    traverse : bool, optional
        By default dask traverses builtin python collections looking for dask
        objects passed to ``compute``. For large collections this can be
        expensive. If none of the arguments contain any dask objects, set
        ``traverse=False`` to avoid doing this traversal.
    scheduler : string, optional
        Which scheduler to use like "threads", "synchronous" or "processes".
        If not provided, the default is to check the global settings first,
        and then fall back to the collection defaults.
    optimize_graph : bool, optional
        If True [default], the optimizations for each collection are applied
        before computation. Otherwise the graph is run as is. This can be
        useful for debugging.
    get : ``None``
        Should be left to ``None`` The get= keyword has been removed.
    kwargs
        Extra keywords to forward to the scheduler function.

    Examples
    --------
    >>> import dask
    >>> import dask.array as da
    >>> a = da.arange(10, chunks=2).sum()
    >>> b = da.arange(10, chunks=2).mean()
    >>> dask.compute(a, b)
    (np.int64(45), np.float64(4.5))

    By default, dask objects inside python collections will also be computed:

    >>> dask.compute({'a': a, 'b': b, 'c': 1})
    ({'a': np.int64(45), 'b': np.float64(4.5), 'c': 1},)
    """

    collections, repack = unpack_collections(*args, traverse=traverse)
    if not collections:
        return args

    schedule = get_scheduler(
        scheduler=scheduler,
        collections=collections,
        get=get,
    )

    dsk = collections_to_dsk(collections, optimize_graph, **kwargs)
    keys, postcomputes = [], []
    for x in collections:
        keys.append(x.__dask_keys__())
        postcomputes.append(x.__dask_postcompute__())

    with shorten_traceback():
        results = schedule(dsk, keys, **kwargs)

    return repack([f(r, *a) for r, (f, a) in zip(results, postcomputes)])


def visualize(
    *args,
    filename="mydask",
    traverse=True,
    optimize_graph=False,
    maxval=None,
    engine: Literal["cytoscape", "ipycytoscape", "graphviz"] | None = None,
    **kwargs,
):
    """
    Visualize several dask graphs simultaneously.

    Requires ``graphviz`` to be installed. All options that are not the dask
    graph(s) should be passed as keyword arguments.

    Parameters
    ----------
    args : object
        Any number of objects. If it is a dask collection (for example, a
        dask DataFrame, Array, Bag, or Delayed), its associated graph
        will be included in the output of visualize. By default, python builtin
        collections are also traversed to look for dask objects (for more
        information see the ``traverse`` keyword). Arguments lacking an
        associated graph will be ignored.
    filename : str or None, optional
        The name of the file to write to disk. If the provided `filename`
        doesn't include an extension, '.png' will be used by default.
        If `filename` is None, no file will be written, and we communicate
        with dot using only pipes.
    format : {'png', 'pdf', 'dot', 'svg', 'jpeg', 'jpg'}, optional
        Format in which to write output file.  Default is 'png'.
    traverse : bool, optional
        By default, dask traverses builtin python collections looking for dask
        objects passed to ``visualize``. For large collections this can be
        expensive. If none of the arguments contain any dask objects, set
        ``traverse=False`` to avoid doing this traversal.
    optimize_graph : bool, optional
        If True, the graph is optimized before rendering.  Otherwise,
        the graph is displayed as is. Default is False.
    color : {None, 'order', 'ages', 'freed', 'memoryincreases', 'memorydecreases', 'memorypressure'}, optional
        Options to color nodes. colormap:

        - None, the default, no colors.
        - 'order', colors the nodes' border based on the order they appear in the graph.
        - 'ages', how long the data of a node is held.
        - 'freed', the number of dependencies released after running a node.
        - 'memoryincreases', how many more outputs are held after the lifetime of a node.
          Large values may indicate nodes that should have run later.
        - 'memorydecreases', how many fewer outputs are held after the lifetime of a node.
          Large values may indicate nodes that should have run sooner.
        - 'memorypressure', the number of data held when the node is run (circle), or
          the data is released (rectangle).
    maxval : {int, float}, optional
        Maximum value for colormap to normalize form 0 to 1.0. Default is ``None``
        will make it the max number of values
    collapse_outputs : bool, optional
        Whether to collapse output boxes, which often have empty labels.
        Default is False.
    verbose : bool, optional
        Whether to label output and input boxes even if the data aren't chunked.
        Beware: these labels can get very long. Default is False.
    engine : {"graphviz", "ipycytoscape", "cytoscape"}, optional.
        The visualization engine to use. If not provided, this checks the dask config
        value "visualization.engine". If that is not set, it tries to import ``graphviz``
        and ``ipycytoscape``, using the first one to succeed.
    **kwargs
       Additional keyword arguments to forward to the visualization engine.

    Examples
    --------
    >>> x.visualize(filename='dask.pdf')  # doctest: +SKIP
    >>> x.visualize(filename='dask.pdf', color='order')  # doctest: +SKIP

    Returns
    -------
    result : IPython.display.Image, IPython.display.SVG, or None
        See dask.dot.dot_graph for more information.

    See Also
    --------
    dask.dot.dot_graph

    Notes
    -----
    For more information on optimization see here:

    https://docs.dask.org/en/latest/optimize.html
    """
    args, _ = unpack_collections(*args, traverse=traverse)

    dsk = dict(collections_to_dsk(args, optimize_graph=optimize_graph))

    return visualize_dsk(
        dsk=dsk,
        filename=filename,
        traverse=traverse,
        optimize_graph=optimize_graph,
        maxval=maxval,
        engine=engine,
        **kwargs,
    )


def visualize_dsk(
    dsk,
    filename="mydask",
    traverse=True,
    optimize_graph=False,
    maxval=None,
    o=None,
    engine: Literal["cytoscape", "ipycytoscape", "graphviz"] | None = None,
    limit=None,
    **kwargs,
):
    color = kwargs.get("color")
    from dask.order import diagnostics, order

    if color in {
        "order",
        "order-age",
        "order-freed",
        "order-memoryincreases",
        "order-memorydecreases",
        "order-memorypressure",
        "age",
        "freed",
        "memoryincreases",
        "memorydecreases",
        "memorypressure",
        "critical",
        "cpath",
    }:
        import matplotlib.pyplot as plt

        if o is None:
            o_stats = order(dsk, return_stats=True)
            o = {k: v.priority for k, v in o_stats.items()}
        elif isinstance(next(iter(o.values())), int):
            o_stats = order(dsk, return_stats=True)
        else:
            o_stats = o
            o = {k: v.priority for k, v in o.items()}

        try:
            cmap = kwargs.pop("cmap")
        except KeyError:
            cmap = plt.cm.plasma
        if isinstance(cmap, str):
            import matplotlib.pyplot as plt

            cmap = getattr(plt.cm, cmap)

        def label(x):
            return str(values[x])

        data_values = None
        if color != "order":
            info = diagnostics(dsk, o)[0]
            if color.endswith("age"):
                values = {key: val.age for key, val in info.items()}
            elif color.endswith("freed"):
                values = {key: val.num_dependencies_freed for key, val in info.items()}
            elif color.endswith("memorypressure"):
                values = {key: val.num_data_when_run for key, val in info.items()}
                data_values = {
                    key: val.num_data_when_released for key, val in info.items()
                }
            elif color.endswith("memoryincreases"):
                values = {
                    key: max(0, val.num_data_when_released - val.num_data_when_run)
                    for key, val in info.items()
                }
            elif color.endswith("memorydecreases"):
                values = {
                    key: max(0, val.num_data_when_run - val.num_data_when_released)
                    for key, val in info.items()
                }
            elif color.split("-")[-1] in {"critical", "cpath"}:
                values = {key: val.critical_path for key, val in o_stats.items()}
            else:
                raise NotImplementedError(color)

            if color.startswith("order-"):

                def label(x):
                    return str(o[x]) + "-" + str(values[x])

        else:
            values = o
        if maxval is None:
            maxval = max(1, max(values.values()))
        colors = {
            k: _colorize(tuple(map(int, cmap(v / maxval, bytes=True))))
            for k, v in values.items()
        }
        if data_values is None:
            data_colors = colors
        else:
            data_colors = {
                k: _colorize(tuple(map(int, cmap(v / maxval, bytes=True))))
                for k, v in values.items()
            }

        kwargs["function_attributes"] = {
            k: {"color": v, "label": label(k)} for k, v in colors.items()
        }
        kwargs["data_attributes"] = {k: {"color": v} for k, v in data_colors.items()}
    elif color:
        raise NotImplementedError("Unknown value color=%s" % color)

    # Determine which engine to dispatch to, first checking the kwarg, then config,
    # then whichever of graphviz or ipycytoscape are installed, in that order.
    engine = engine or config.get("visualization.engine", None)

    if not engine:
        try:
            import graphviz  # noqa: F401

            engine = "graphviz"
        except ImportError:
            try:
                import ipycytoscape  # noqa: F401

                engine = "cytoscape"
            except ImportError:
                pass
    if engine == "graphviz":
        from dask.dot import dot_graph

        return dot_graph(dsk, filename=filename, **kwargs)
    elif engine in ("cytoscape", "ipycytoscape"):
        from dask.dot import cytoscape_graph

        return cytoscape_graph(dsk, filename=filename, **kwargs)
    elif engine is None:
        raise RuntimeError(
            "No visualization engine detected, please install graphviz or ipycytoscape"
        )
    else:
        raise ValueError(f"Visualization engine {engine} not recognized")


def persist(*args, traverse=True, optimize_graph=True, scheduler=None, **kwargs):
    """Persist multiple Dask collections into memory

    This turns lazy Dask collections into Dask collections with the same
    metadata, but now with their results fully computed or actively computing
    in the background.

    For example a lazy dask.array built up from many lazy calls will now be a
    dask.array of the same shape, dtype, chunks, etc., but now with all of
    those previously lazy tasks either computed in memory as many small :class:`numpy.array`
    (in the single-machine case) or asynchronously running in the
    background on a cluster (in the distributed case).

    This function operates differently if a ``dask.distributed.Client`` exists
    and is connected to a distributed scheduler.  In this case this function
    will return as soon as the task graph has been submitted to the cluster,
    but before the computations have completed.  Computations will continue
    asynchronously in the background.  When using this function with the single
    machine scheduler it blocks until the computations have finished.

    When using Dask on a single machine you should ensure that the dataset fits
    entirely within memory.

    Examples
    --------
    >>> df = dd.read_csv('/path/to/*.csv')  # doctest: +SKIP
    >>> df = df[df.name == 'Alice']  # doctest: +SKIP
    >>> df['in-debt'] = df.balance < 0  # doctest: +SKIP
    >>> df = df.persist()  # triggers computation  # doctest: +SKIP

    >>> df.value().min()  # future computations are now fast  # doctest: +SKIP
    -10
    >>> df.value().max()  # doctest: +SKIP
    100

    >>> from dask import persist  # use persist function on multiple collections
    >>> a, b = persist(a, b)  # doctest: +SKIP

    Parameters
    ----------
    *args: Dask collections
    scheduler : string, optional
        Which scheduler to use like "threads", "synchronous" or "processes".
        If not provided, the default is to check the global settings first,
        and then fall back to the collection defaults.
    traverse : bool, optional
        By default dask traverses builtin python collections looking for dask
        objects passed to ``persist``. For large collections this can be
        expensive. If none of the arguments contain any dask objects, set
        ``traverse=False`` to avoid doing this traversal.
    optimize_graph : bool, optional
        If True [default], the graph is optimized before computation.
        Otherwise the graph is run as is. This can be useful for debugging.
    **kwargs
        Extra keywords to forward to the scheduler function.

    Returns
    -------
    New dask collections backed by in-memory data
    """
    collections, repack = unpack_collections(*args, traverse=traverse)
    if not collections:
        return args

    schedule = get_scheduler(scheduler=scheduler, collections=collections)

    if inspect.ismethod(schedule):
        try:
            from distributed.client import default_client
        except ImportError:
            pass
        else:
            try:
                client = default_client()
            except ValueError:
                pass
            else:
                if client.get == schedule:
                    results = client.persist(
                        collections, optimize_graph=optimize_graph, **kwargs
                    )
                    return repack(results)

    dsk = collections_to_dsk(collections, optimize_graph, **kwargs)
    keys, postpersists = [], []
    for a in collections:
        a_keys = list(flatten(a.__dask_keys__()))
        rebuild, state = a.__dask_postpersist__()
        keys.extend(a_keys)
        postpersists.append((rebuild, a_keys, state))

    with shorten_traceback():
        results = schedule(dsk, keys, **kwargs)

    d = dict(zip(keys, results))
    results2 = [r({k: d[k] for k in ks}, *s) for r, ks, s in postpersists]
    return repack(results2)


def _colorize(t):
    """Convert (r, g, b) triple to "#RRGGBB" string

    For use with ``visualize(color=...)``

    Examples
    --------
    >>> _colorize((255, 255, 255))
    '#FFFFFF'
    >>> _colorize((0, 32, 128))
    '#002080'
    """
    t = t[:3]
    i = sum(v * 256 ** (len(t) - i - 1) for i, v in enumerate(t))
    h = hex(int(i))[2:].upper()
    h = "0" * (6 - len(h)) + h
    return "#" + h


named_schedulers: dict[str, SchedulerGetCallable] = {
    "sync": local.get_sync,
    "synchronous": local.get_sync,
    "single-threaded": local.get_sync,
}

if not EMSCRIPTEN:
    from dask import threaded

    named_schedulers.update(
        {
            "threads": threaded.get,
            "threading": threaded.get,
        }
    )

    from dask import multiprocessing as dask_multiprocessing

    named_schedulers.update(
        {
            "processes": dask_multiprocessing.get,
            "multiprocessing": dask_multiprocessing.get,
        }
    )


get_err_msg = """
The get= keyword has been removed.

Please use the scheduler= keyword instead with the name of
the desired scheduler like 'threads' or 'processes'

    x.compute(scheduler='single-threaded')
    x.compute(scheduler='threads')
    x.compute(scheduler='processes')

or with a function that takes the graph and keys

    x.compute(scheduler=my_scheduler_function)

or with a Dask client

    x.compute(scheduler=client)
""".strip()


def _ensure_not_async(client):
    if client.asynchronous:
        if fallback := config.get("admin.async-client-fallback", None):
            warnings.warn(
                "Distributed Client detected but Client instance is "
                f"asynchronous. Falling back to `{fallback}` scheduler. "
                "To use an asynchronous Client, please use "
                "``Client.compute`` and ``Client.gather`` "
                "instead of the top level ``dask.compute``",
                UserWarning,
            )
            return get_scheduler(scheduler=fallback)
        else:
            raise RuntimeError(
                "Attempting to use an asynchronous "
                "Client in a synchronous context of `dask.compute`"
            )
    return client.get


def get_scheduler(get=None, scheduler=None, collections=None, cls=None):
    """Get scheduler function

    There are various ways to specify the scheduler to use:

    1.  Passing in scheduler= parameters
    2.  Passing these into global configuration
    3.  Using a dask.distributed default Client
    4.  Using defaults of a dask collection

    This function centralizes the logic to determine the right scheduler to use
    from those many options
    """
    if get:
        raise TypeError(get_err_msg)

    if scheduler is not None:
        if callable(scheduler):
            return scheduler
        elif "Client" in type(scheduler).__name__ and hasattr(scheduler, "get"):
            return _ensure_not_async(scheduler)
        elif isinstance(scheduler, str):
            scheduler = scheduler.lower()

            client_available = False
            if _distributed_available():
                assert _DistributedClient is not None
                with suppress(ValueError):
                    _DistributedClient.current(allow_global=True)
                    client_available = True
            if scheduler in named_schedulers:
                if client_available:
                    warnings.warn(
                        "Running on a single-machine scheduler when a distributed client "
                        "is active might lead to unexpected results."
                    )
                return named_schedulers[scheduler]
            elif scheduler in ("dask.distributed", "distributed"):
                if not client_available:
                    raise RuntimeError(
                        f"Requested {scheduler} scheduler but no Client active."
                    )
                assert _get_distributed_client is not None
                client = _get_distributed_client()
                return _ensure_not_async(client)
            else:
                raise ValueError(
                    "Expected one of [distributed, %s]"
                    % ", ".join(sorted(named_schedulers))
                )
        elif isinstance(scheduler, Executor):
            # Get `num_workers` from `Executor`'s `_max_workers` attribute.
            # If undefined, fallback to `config` or worst case CPU_COUNT.
            num_workers = getattr(scheduler, "_max_workers", None)
            if num_workers is None:
                num_workers = config.get("num_workers", CPU_COUNT)
            assert isinstance(num_workers, Integral) and num_workers > 0
            return partial(local.get_async, scheduler.submit, num_workers)
        else:
            raise ValueError("Unexpected scheduler: %s" % repr(scheduler))
        # else:  # try to connect to remote scheduler with this name
        #     return get_client(scheduler).get

    if config.get("scheduler", None):
        return get_scheduler(scheduler=config.get("scheduler", None))

    if config.get("get", None):
        raise ValueError(get_err_msg)

    try:
        from distributed import get_client

        return _ensure_not_async(get_client())
    except (ImportError, ValueError):
        pass

    if cls is not None:
        return cls.__dask_scheduler__

    if collections:
        collections = [c for c in collections if c is not None]
    if collections:
        get = collections[0].__dask_scheduler__
        if not all(c.__dask_scheduler__ == get for c in collections):
            raise ValueError(
                "Compute called on multiple collections with "
                "differing default schedulers. Please specify a "
                "scheduler=` parameter explicitly in compute or "
                "globally with `dask.config.set`."
            )
        return get

    return None


def wait(x, timeout=None, return_when="ALL_COMPLETED"):
    """Wait until computation has finished

    This is a compatibility alias for ``dask.distributed.wait``.
    If it is applied onto Dask collections without Dask Futures or if Dask
    distributed is not installed then it is a no-op
    """
    try:
        from distributed import wait

        return wait(x, timeout=timeout, return_when=return_when)
    except (ImportError, ValueError):
        return x


def get_collection_names(collection) -> set[str]:
    """Infer the collection names from the dask keys, under the assumption that all keys
    are either tuples with matching first element, and that element is a string, or
    there is exactly one key and it is a string.

    Examples
    --------
    >>> a.__dask_keys__()  # doctest: +SKIP
    ["foo", "bar"]
    >>> get_collection_names(a)  # doctest: +SKIP
    {"foo", "bar"}
    >>> b.__dask_keys__()  # doctest: +SKIP
    [[("foo-123", 0, 0), ("foo-123", 0, 1)], [("foo-123", 1, 0), ("foo-123", 1, 1)]]
    >>> get_collection_names(b)  # doctest: +SKIP
    {"foo-123"}
    """
    if not is_dask_collection(collection):
        raise TypeError(f"Expected Dask collection; got {type(collection)}")
    return {get_name_from_key(k) for k in flatten(collection.__dask_keys__())}


def get_name_from_key(key: Key) -> str:
    """Given a dask collection's key, extract the collection name.

    Parameters
    ----------
    key: string or tuple
        Dask collection's key, which must be either a single string or a tuple whose
        first element is a string (commonly referred to as a collection's 'name'),

    Examples
    --------
    >>> get_name_from_key("foo")
    'foo'
    >>> get_name_from_key(("foo-123", 1, 2))
    'foo-123'
    """
    if isinstance(key, tuple) and key and isinstance(key[0], str):
        return key[0]
    if isinstance(key, str):
        return key
    raise TypeError(f"Expected str or a tuple starting with str; got {key!r}")


KeyOrStrT = TypeVar("KeyOrStrT", Key, str)


def replace_name_in_key(key: KeyOrStrT, rename: Mapping[str, str]) -> KeyOrStrT:
    """Given a dask collection's key, replace the collection name with a new one.

    Parameters
    ----------
    key: string or tuple
        Dask collection's key, which must be either a single string or a tuple whose
        first element is a string (commonly referred to as a collection's 'name'),
    rename:
        Mapping of zero or more names from : to. Extraneous names will be ignored.
        Names not found in this mapping won't be replaced.

    Examples
    --------
    >>> replace_name_in_key("foo", {})
    'foo'
    >>> replace_name_in_key("foo", {"foo": "bar"})
    'bar'
    >>> replace_name_in_key(("foo-123", 1, 2), {"foo-123": "bar-456"})
    ('bar-456', 1, 2)
    """
    if isinstance(key, tuple) and key and isinstance(key[0], str):
        return (rename.get(key[0], key[0]),) + key[1:]
    if isinstance(key, str):
        return rename.get(key, key)
    raise TypeError(f"Expected str or a tuple starting with str; got {key!r}")


def clone_key(key: KeyOrStrT, seed: Hashable) -> KeyOrStrT:
    """Clone a key from a Dask collection, producing a new key with the same prefix and
    indices and a token which is a deterministic function of the previous key and seed.

    Examples
    --------
    >>> clone_key("x", 123)  # doctest: +SKIP
    'x-c4fb64ccca807af85082413d7ef01721'
    >>> clone_key("inc-cbb1eca3bafafbb3e8b2419c4eebb387", 123)  # doctest: +SKIP
    'inc-bc629c23014a4472e18b575fdaf29ee7'
    >>> clone_key(("sum-cbb1eca3bafafbb3e8b2419c4eebb387", 4, 3), 123)  # doctest: +SKIP
    ('sum-c053f3774e09bd0f7de6044dbc40e71d', 4, 3)
    """
    if isinstance(key, tuple) and key and isinstance(key[0], str):
        return (clone_key(key[0], seed),) + key[1:]
    if isinstance(key, str):
        prefix = key_split(key)
        return prefix + "-" + tokenize(key, seed)
    raise TypeError(f"Expected str or a tuple starting with str; got {key!r}")
