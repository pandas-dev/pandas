"""
Metadata propagation through pandas operations.

This module contains the infrastructure for propagating ``NDFrame._metadata``
through operations. We perform an operation (say :meth:`pandas.Series.copy`) that
returns an ``NDFrame`` and would like to propagate the metadata (say ``Series.name``)
from ``self`` to the new ``NDFrame``.
"""
from typing import Any, Dict

from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries

from pandas._typing import FrameOrSeries


class PandasMetadataType(type):
    """
    Metaclass controlling creation of metadata finalizers.

    This ensures we have one finalizer instance per name, and
    provides a place to look up finalizer per name.
    """

    # TODO(Py35): Replace metaclass with __subclass_init__

    _instances = {}  # type: Dict[str, "PandasMetadata"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        name = args[0]
        if name in self._instances:
            return self._instances[name]
        else:
            new = super().__call__(name, *args, **kwds)
            self._instances[name] = new
            return new


class PandasMetadata(metaclass=PandasMetadataType):
    """
    Dispatch metadata finalization for pandas metadata.

    Parameters
    ----------
    name : str
        The name of the attribute being finalized.

    Examples
    --------
    If you want the default resolution (copy from a source NDFrame
    to a new NDFrame), you can just create an instance

    >>> mymeta = PandasMetadata("mymeta")

    If you need custom metadata resolution, you'll need to subclass.

    >>> class IncrementMetadata:
    ...     def finalize_default(self, new, other):
    ...         setattr(new, self.attr, getattr(other, self.name, -1) + 1)

    >>> increment_metadata = IncrementMetadata("attr")
    """

    def __init__(self, name: str):
        self.name = name

    def __repr__(self):
        return "PandasMetadata(name='{}')".format(self.name)

    def finalize_default(self, new: "FrameOrSeries", other: Any):
        # TODO: tighten the type on Any. Union[NDFrame, _Concatenator, ...]
        """
        The default finalizer when this method, attribute hasn't been overridden.

        This copies the ``_metadata`` attribute from ``other`` to ``self``, modifying
        ``self`` inplace.

        Parameters
        ----------
        new : NDFrame
            The newly created NDFrame being finalized.
        other : Any
            The source object attributes will be extracted from.
        """
        # TODO: check perf on this isinstance.
        if isinstance(other, (ABCSeries, ABCDataFrame)):
            object.__setattr__(new, self.name, getattr(other, self.name, None))

    def finalize_copy(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.finalize_default(new, other)

    def finalize_concat(self, new: "FrameOrSeries", other: Any):
        # TODO: type other as _Concatenator
        self.finalize_default(new, other)


class NameMetadata(PandasMetadata):
    """Finalization for Series.name"""


# TODO: having to create this here feels weird.
name_metadata = NameMetadata("name")

# For backwards compat. Do we care about this?
# We can pretty easily deprecate, require subclasses to make their
# own instance.
default_finalizer = PandasMetadata("pandas")
