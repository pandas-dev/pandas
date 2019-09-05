"""
Metadata propagation through pandas operations.

This module contains the infrastructure for propagating ``NDFrame._metadata``
through operations. We perform an operation (say :meth:`pandas.Series.copy`) that
returns an ``NDFrame`` and would like to propagate the metadata (say ``Series.name``)
from ``self`` to the new ``NDFrame``.

.. note::

   Currently, pandas doesn't provide a clean, documented API on

   * which methods call finalize
   * the types passed to finalize for each method

   This is a known limitation we would like to address in the future.
"""
from collections import defaultdict
from functools import wraps
from typing import TYPE_CHECKING, Any, Callable, Union

from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries

if TYPE_CHECKING:
    from pandas.core.generic import NDFrame

dispatch = defaultdict(dict)
dispatch_method_type = Union[Callable[..., "NDFrame"], str]


def key_of(method):
    if isinstance(method, str):
        # TODO: figure out if this is OK. May be necessary when we have
        #   things like pd.merge and DataFrame.merge that hit the same finalize.
        return method
    elif method:
        return method.__module__, method.__name__


class PandasMetadata:
    """
    Dispatch metadata finalization for pandas metadata.

    Users should instantiate a single `PandasMetadata` instance
    for their piece of metadata and register finalizers for various
    pandas methods using :meth:`PandsaMetadata.register`.

    Parameters
    ----------
    name : str
        The name of the attribute being finalized.

    Examples
    --------
    >>> maxmeta = PandasMetadata("attr")

    Register a finalizer for a given pandas method:

    >>> @maxmeta.register(pd.concat)
    ... def _(new, concatenator):
    ...     new.attr = max(x.attr_meta for x in concatenator.objs)

    >>> pd.DataFrame._metadata = ['attr']
    >>> x = pd.DataFrame({"x"}); x.attr = 1
    >>> y = pd.DataFrame({"y"}); y.attr = 2
    >>> pd.concat([x, y]).attr
    2
    """

    def __init__(self, name: str):
        self.name = name

    def register(self, pandas_method: dispatch_method_type):
        """
        A decorator to register a finalizer for a specific pandas method.

        Parameters
        ----------
        pandas_method : callable or str
            A pandas method, like :meth:`pandas.concat`, that this finalizer
            should be used for. The function being decorated will be called
            with the relevant arguments (typically the output and the source NDFrame).
            When `NDFrame.__finalize__` is called as a result of `pandas_method`,
            the registered finalizer will be called.
        """

        def decorate(func):
            # TODO: warn of collisions?
            dispatch[key_of(pandas_method)][self.name] = func

            @wraps(func)
            def wrapper(*args, **kwargs):
                return func(*args, **kwargs)

            return wrapper

        return decorate


def default_finalizer(new: "NDFrame", other: Any, *, name: str):
    """
    The default finalizer when this method, attribute hasn't been overridden.

    This copies the ``_metadata`` attribute from ``other`` to ``self``, modifying
    ``self`` inplace.

    Parameters
    ----------
    new : NDFrame
        The newly created NDFrame being finalized.
    other : NDFrame
        The source NDFrame attributes will be extracted from.
    """
    object.__setattr__(new, name, getattr(other, name, None))


# ----------------------------------------------------------------------------
# Pandas Internals.


def ndframe_finalize(new: "NDFrame", other: Any, method: dispatch_method_type):
    """
    Finalize a new NDFrame.

    The finalizer is looked up from finalizers registered with PandasMetadata.
    `new` is modified inplace, and nothing is returned.

    Parameters
    ----------
    new : NDFrame
    other : NDFrame
        Or a list of them? TBD
    method : callable or str
    """
    # To avoid one isinstance per _metadata name, we check up front.
    # Most of the time `other` is an ndframe, but in some cases (e.g. concat)
    # it's `_Concatenator` object
    other_is_ndframe = isinstance(other, (ABCSeries, ABCDataFrame))

    for name in new._metadata:
        finalizer = dispatch.get(key_of(method), {}).get(name)

        if finalizer:
            finalizer(new, other)
        elif other_is_ndframe:
            default_finalizer(new, other, name=name)
