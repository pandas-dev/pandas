from __future__ import annotations

import abc
from collections.abc import Callable, Hashable, Mapping, Sequence
from enum import Enum
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Protocol,
    TypeVar,
    Union,
    runtime_checkable,
)

if TYPE_CHECKING:
    # IPython import is relatively slow. Avoid if not necessary
    # TODO import from typing (requires Python >=3.10)
    from typing import TypeAlias

    from IPython.display import DisplayObject

CollType = TypeVar("CollType", bound="DaskCollection")
CollType_co = TypeVar("CollType_co", bound="DaskCollection", covariant=True)
PostComputeCallable = Callable


Key: TypeAlias = Union[str, bytes, int, float, tuple["Key", ...]]
# FIXME: This type is a little misleading. Low level graphs are often
# MutableMappings but HLGs are not
Graph: TypeAlias = Mapping[Key, Any]
# Potentially nested list of Dask keys
NestedKeys: TypeAlias = list[Union[Key, "NestedKeys"]]


class SchedulerGetCallable(Protocol):
    """Protocol defining the signature of a ``__dask_scheduler__`` callable."""

    def __call__(
        self,
        dsk: Graph,
        keys: Sequence[Key] | Key,
        **kwargs: Any,
    ) -> Any:
        """Method called as the default scheduler for a collection.

        Parameters
        ----------
        dsk :
            The task graph.
        keys :
            Key(s) corresponding to the desired data.
        **kwargs :
            Additional arguments.

        Returns
        -------
        Any
            Result(s) associated with `keys`

        """
        raise NotImplementedError("Inheriting class must implement this method.")


class PostPersistCallable(Protocol[CollType_co]):
    """Protocol defining the signature of a ``__dask_postpersist__`` callable."""

    def __call__(
        self,
        dsk: Graph,
        *args: Any,
        rename: Mapping[str, str] | None = None,
    ) -> CollType_co:
        """Method called to rebuild a persisted collection.

        Parameters
        ----------
        dsk: Mapping
            A mapping which contains at least the output keys returned
            by __dask_keys__().
        *args : Any
            Additional optional arguments If no extra arguments are
            necessary, it must be an empty tuple.
        rename : Mapping[str, str], optional
            If defined, it indicates that output keys may be changing
            too; e.g. if the previous output of :meth:`__dask_keys__`
            was ``[("a", 0), ("a", 1)]``, after calling
            ``rebuild(dsk, *extra_args, rename={"a": "b"})``
            it must become ``[("b", 0), ("b", 1)]``.
            The ``rename`` mapping may not contain the collection
            name(s); in such case the associated keys do not change.
            It may contain replacements for unexpected names, which
            must be ignored.

        Returns
        -------
        Collection
            An equivalent Dask collection with the same keys as
            computed through a different graph.

        """
        raise NotImplementedError("Inheriting class must implement this method.")


@runtime_checkable
class DaskCollection(Protocol):
    """Protocol defining the interface of a Dask collection."""

    @abc.abstractmethod
    def __dask_graph__(self) -> Graph:
        """The Dask task graph.

        The core Dask collections (Array, DataFrame, Bag, and Delayed)
        use a :py:class:`~dask.highlevelgraph.HighLevelGraph` to
        represent the collection task graph. It is also possible to
        represent the task graph as a low level graph using a Python
        dictionary.

        Returns
        -------
        Mapping
            The Dask task graph. If the instance returns a
            :py:class:`dask.highlevelgraph.HighLevelGraph` then the
            :py:func:`__dask_layers__` method must be implemented, as
            defined by the :py:class:`~dask.typing.HLGDaskCollection`
            protocol.

        """
        raise NotImplementedError("Inheriting class must implement this method.")

    @abc.abstractmethod
    def __dask_keys__(self) -> NestedKeys:
        """The output keys of the task graph.

        Note that there are additional constraints on keys for a Dask
        collection than those described in the :doc:`task graph
        specification documentation <spec>`. These additional
        constraints are described below.

        All keys must either be non-empty strings or tuples where the first element is a
        non-empty string, followed by zero or more arbitrary str, bytes, int, float, or
        tuples thereof. The non-empty string is commonly known as the *collection name*.
        All collections embedded in the dask package have exactly one name, but this is
        not a requirement.

        These are all valid outputs:

        - ``[]``
        - ``["x", "y"]``
        - ``[[("y", "a", 0), ("y", "a", 1)], [("y", "b", 0), ("y", "b", 1)]``

        Returns
        -------
        list
            A possibly nested list of keys that represent the outputs
            of the graph. After computation, the results will be
            returned in the same layout, with the keys replaced with
            their corresponding outputs.

        """
        raise NotImplementedError("Inheriting class must implement this method.")

    @abc.abstractmethod
    def __dask_postcompute__(self) -> tuple[PostComputeCallable, tuple]:
        """Finalizer function and optional arguments to construct final result.

        Upon computation each key in the collection will have an in
        memory result, the postcompute function combines each key's
        result into a final in memory representation. For example,
        dask.array.Array concatenates the arrays at each chunk into a
        final in-memory array.

        Returns
        -------
        PostComputeCallable
            Callable that receives the sequence of the results of each
            final key along with optional arguments. An example signature
            would be ``finalize(results: Sequence[Any], *args)``.
        tuple[Any, ...]
            Optional arguments passed to the function following the
            key results (the `*args` part of the
            ``PostComputeCallable``. If no additional arguments are to
            be passed then this must be an empty tuple.

        """
        raise NotImplementedError("Inheriting class must implement this method.")

    @abc.abstractmethod
    def __dask_postpersist__(self) -> tuple[PostPersistCallable, tuple]:
        """Rebuilder function and optional arguments to construct a persisted collection.

        See also the documentation for :py:class:`dask.typing.PostPersistCallable`.

        Returns
        -------
        PostPersistCallable
            Callable that rebuilds the collection. The signature
            should be
            ``rebuild(dsk: Mapping, *args: Any, rename: Mapping[str, str] | None)``
            (as defined by the
            :py:class:`~dask.typing.PostPersistCallable` protocol).
            The callable should return an equivalent Dask collection
            with the same keys as `self`, but with results that are
            computed through a different graph. In the case of
            :py:func:`dask.persist`, the new graph will have just the
            output keys and the values already computed.
        tuple[Any, ...]
            Optional arguments passed to the rebuild callable. If no
            additional arguments are to be passed then this must be an
            empty tuple.

        """
        raise NotImplementedError("Inheriting class must implement this method.")

    @abc.abstractmethod
    def __dask_tokenize__(self) -> Hashable:
        """Value that must fully represent the object."""
        raise NotImplementedError("Inheriting class must implement this method.")

    __dask_optimize__: Any
    """Given a graph and keys, return a new optimized graph.

    This method can be either a ``staticmethod`` or a ``classmethod``,
    but not an ``instancemethod``. For example implementations see the
    definitions of ``__dask_optimize__`` in the core Dask collections:
    ``dask.array.Array``, ``dask.dataframe.DataFrame``, etc.

    Note that graphs and keys are merged before calling
    ``__dask_optimize__``; as such, the graph and keys passed to
    this method may represent more than one collection sharing the
    same optimize method.

    Parameters
    ----------
    dsk : Graph
        The merged graphs from all collections sharing the same
        ``__dask_optimize__`` method.
    keys : Sequence[Key]
        A list of the outputs from ``__dask_keys__`` from all
        collections sharing the same ``__dask_optimize__`` method.
    **kwargs : Any
        Extra keyword arguments forwarded from the call to
        ``compute`` or ``persist``. Can be used or ignored as
        needed.

    Returns
    -------
    MutableMapping
        The optimized Dask graph.

    """

    # FIXME: It is currently not possible to define a staticmethod from a callback protocol
    # Also, the definition in `is_dask_collection` cannot be satisfied by a
    # protocol / typing check
    # staticmethod[SchedulerGetCallable]
    __dask_scheduler__: staticmethod
    """The default scheduler ``get`` to use for this object.

    Usually attached to the class as a staticmethod, e.g.:

    >>> import dask.threaded
    >>> class MyCollection:
    ...     # Use the threaded scheduler by default
    ...     __dask_scheduler__ = staticmethod(dask.threaded.get)

    """

    @abc.abstractmethod
    def compute(self, **kwargs: Any) -> Any:
        """Compute this dask collection.

        This turns a lazy Dask collection into its in-memory
        equivalent. For example a Dask array turns into a NumPy array
        and a Dask dataframe turns into a Pandas dataframe. The entire
        dataset must fit into memory before calling this operation.

        Parameters
        ----------
        scheduler : string, optional
            Which scheduler to use like "threads", "synchronous" or
            "processes". If not provided, the default is to check the
            global settings first, and then fall back to the
            collection defaults.
        optimize_graph : bool, optional
            If True [default], the graph is optimized before
            computation. Otherwise the graph is run as is. This can be
            useful for debugging.
        kwargs :
            Extra keywords to forward to the scheduler function.

        Returns
        -------
        The collection's computed result.

        See Also
        --------
        dask.compute

        """
        raise NotImplementedError("Inheriting class must implement this method.")

    @abc.abstractmethod
    def persist(self: CollType, **kwargs: Any) -> CollType:
        """Persist this dask collection into memory

        This turns a lazy Dask collection into a Dask collection with
        the same metadata, but now with the results fully computed or
        actively computing in the background.

        The action of function differs significantly depending on the
        active task scheduler. If the task scheduler supports
        asynchronous computing, such as is the case of the
        dask.distributed scheduler, then persist will return
        *immediately* and the return value's task graph will contain
        Dask Future objects. However if the task scheduler only
        supports blocking computation then the call to persist will
        *block* and the return value's task graph will contain
        concrete Python results.

        This function is particularly useful when using distributed
        systems, because the results will be kept in distributed
        memory, rather than returned to the local process as with
        compute.

        Parameters
        ----------
        scheduler : string, optional
            Which scheduler to use like "threads", "synchronous" or
            "processes". If not provided, the default is to check the
            global settings first, and then fall back to the
            collection defaults.
        optimize_graph : bool, optional
            If True [default], the graph is optimized before
            computation. Otherwise the graph is run as is. This can be
            useful for debugging.
        **kwargs
            Extra keywords to forward to the scheduler function.

        Returns
        -------
        New dask collections backed by in-memory data

        See Also
        --------
        dask.persist

        """
        raise NotImplementedError("Inheriting class must implement this method.")

    @abc.abstractmethod
    def visualize(
        self,
        filename: str = "mydask",
        format: str | None = None,
        optimize_graph: bool = False,
        **kwargs: Any,
    ) -> DisplayObject | None:
        """Render the computation of this object's task graph using graphviz.

        Requires ``graphviz`` to be installed.

        Parameters
        ----------
        filename : str or None, optional
            The name of the file to write to disk. If the provided
            `filename` doesn't include an extension, '.png' will be
            used by default. If `filename` is None, no file will be
            written, and we communicate with dot using only pipes.
        format : {'png', 'pdf', 'dot', 'svg', 'jpeg', 'jpg'}, optional
            Format in which to write output file. Default is 'png'.
        optimize_graph : bool, optional
            If True, the graph is optimized before rendering.
            Otherwise, the graph is displayed as is. Default is False.
        color: {None, 'order'}, optional
            Options to color nodes. Provide ``cmap=`` keyword for
            additional colormap
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
        raise NotImplementedError("Inheriting class must implement this method.")


@runtime_checkable
class HLGDaskCollection(DaskCollection, Protocol):
    """Protocol defining a Dask collection that uses HighLevelGraphs.

    This protocol is nearly identical to
    :py:class:`~dask.typing.DaskCollection`, with the addition of the
    ``__dask_layers__`` method (required for collections backed by
    high level graphs).

    """

    @abc.abstractmethod
    def __dask_layers__(self) -> Sequence[str]:
        """Names of the HighLevelGraph layers."""
        raise NotImplementedError("Inheriting class must implement this method.")


class _NoDefault(Enum):
    """typing-aware constant to detect when the user omits a parameter and you can't use
    None.

    Copied from pandas._libs.lib._NoDefault.

    Usage
    -----
    from dask.typing import NoDefault, no_default

    def f(x: int | None | NoDefault = no_default) -> int:
        if x is no_default:
            ...
    """

    no_default = "NO_DEFAULT"

    def __repr__(self) -> str:
        return "<no_default>"


no_default = _NoDefault.no_default
NoDefault = Literal[_NoDefault.no_default]
