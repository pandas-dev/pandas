"""
Metadata propagation through pandas operations.

This module contains the infrastructure for propagating ``NDFrame._metadata``
through operations. We perform an operation (say :meth:`pandas.Series.copy`) that
returns an ``NDFrame`` and would like to propagate the metadata (say ``Series.name``)
from ``self`` to the new ``NDFrame``.
"""
from typing import TYPE_CHECKING, Any, Dict, Union

from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries

if TYPE_CHECKING:
    from pandas._typing import (
        FrameOrSeries,
        Series,
        Concatenator,
        OrderedMerge,
        MergeOperation,
        AsOfMerge,
    )


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

    def __call__(self, name, *args: Any, **kwds: Any) -> Any:  # type: ignore
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
    ...     def default(self, new, other):
    ...         setattr(new, self.attr, getattr(other, self.name, -1) + 1)

    >>> increment_metadata = IncrementMetadata("attr")
    """

    def __init__(self, name: str):
        self.name = name

    def __repr__(self):
        return "PandasMetadata(name='{}')".format(self.name)

    def default(
        self, new: "FrameOrSeries", other: Union["FrameOrSeries", "Concatenator"]
    ):
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

    def sort_index(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        self.default(new, other)

    def sort_values(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        self.default(new, other)

    def copy(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def rename(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    # TODO: these are private methods.
    def getitem_multilevel(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        self.default(new, other)

    def simple_replace(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def array_wrap(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def align_frame(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        # TODO: Check the types here. We end up calling __finalize__ twice.
        #   left.__finalize__(self)
        #   right.__finalize__(other)
        #  I think the type on right is Frame->Frame.
        return self.default(new, other)

    def align_series(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        # TODO: Check the types here. We end up calling __finalize__ twice.
        #   left.__finalize__(self)
        #   right.__finalize__(other)
        #  I think the type on right is Series->Series.
        return self.default(new, other)

    def consolidate(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def convert(self, new: "FrameOrSeries", other: "FrameOrSeries"):  # XXX: private
        return self.default(new, other)

    def get_bool_data(
        self, new: "FrameOrSeries", other: "FrameOrSeries"
    ):  # XXX: private
        return self.default(new, other)

    def get_numeric_data(
        self, new: "FrameOrSeries", other: "FrameOrSeries"
    ):  # XXX: private
        return self.default(new, other)

    def reindex_with_indexers(
        self, new: "FrameOrSeries", other: "FrameOrSeries"
    ):  # XXX: private
        return self.default(new, other)

    def slice(self, new: "FrameOrSeries", other: "FrameOrSeries"):  # XXX: private
        return self.default(new, other)

    def to_dict_of_blocks(
        self, new: "FrameOrSeries", other: "FrameOrSeries"
    ):  # XXX: private
        return self.default(new, other)

    def where(self, new: "FrameOrSeries", other: "FrameOrSeries"):  # XXX: private
        return self.default(new, other)

    def astype(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def fillna(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def infer_objects(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def interpolate(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def rank(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def replace(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def shift(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def slice_shift(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def swapaxes(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def take(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def transpose(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def tshift(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def tz_convert(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def tz_localize(self, new: "FrameOrSeries", other: "FrameOrSeries"):
        return self.default(new, other)

    def cum_func(self, new: "FrameOrSeries", other: "FrameOrSeries"):  # XXX: private
        return self.default(new, other)

    # ops
    # TODO: constrcut_result, wrapper
    # def construct_result(self):
    def flex_wrapper(self, new: "FrameOrSeries", other: Any):
        return self.default(new, other)

    # ----------------------------------------------------------------------------------
    # Series
    # TODO: How should we namespace these?
    def series_getitem(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_get_values(self, new: "Series", other: "Series"):  # XXX: private
        return self.default(new, other)

    def series_get_values_tuple(self, new: "Series", other: "Series"):  # XXX: private
        return self.default(new, other)

    def series_apply(self, new: "Series", other: "Series"):
        # TODO: We don't call __finalize__ when expanding to a DataFrame. May need to
        # update the type
        return self.default(new, other)

    def series_argsort(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_count(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_diff(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_dot(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_isin(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_map(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_repeat(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_reset_index(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_round(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_sort_index(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_sort_values(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_swaplevel(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_take(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_to_period(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_to_sparse(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_to_timestamp(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_view(self, new: "Series", other: "Series"):
        return self.default(new, other)

    def series_duplicated(self, new: "Series", other: "Series"):
        return self.default(new, other)

    # Other finalizers
    def concat(self, new: "FrameOrSeries", other: "Concatenator"):
        self.default(new, other)

    def merge(
        self, new: "FrameOrSeries", other: "MergeOperation"
    ):  # todo: specialize? merge_asof?
        self.default(new, other)

    def ordered_merge(
        self, new: "FrameOrSeries", other: "OrderedMerge"
    ):  # todo: specialize? merge_asof?
        self.default(new, other)

    def asof_merge(
        self, new: "FrameOrSeries", other: "AsOfMerge"
    ):  # todo: specialize? merge_asof?
        self.default(new, other)

    # TODO: Sparse? Just ignoring for now I think.


class NameMetadata(PandasMetadata):
    """Finalization for Series.name"""


# TODO: having to create this here feels weird.
name_metadata = NameMetadata("name")

# For backwards compat. Do we care about this?
# We can pretty easily deprecate, require subclasses to make their
# own instance.
default_finalizer = PandasMetadata("pandas")
