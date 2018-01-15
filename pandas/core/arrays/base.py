"""An interface for extending pandas with custom arrays."""
import abc

import numpy as np

from pandas.compat import add_metaclass


_not_implemented_message = "{} does not implement {}."


@add_metaclass(abc.ABCMeta)
class ExtensionArray(object):
    """Abstract base class for custom array types

    pandas will recognize instances of this class as proper arrays
    with a custom type and will not attempt to coerce them to objects.

    Subclasses are expected to implement the following methods.
    """
    # ------------------------------------------------------------------------
    # Must be a Sequence
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def __getitem__(self, item):
        """Select a subset of self

        Notes
        -----
        As a sequence, __getitem__ should expect integer or slice ``key``.

        For slice ``key``, you should return an instance of yourself, even
        if the slice is length 0 or 1.

        For scalar ``key``, you may return a scalar suitable for your type.
        The scalar need not be an instance or subclass of your array type.
        """
        # type (Any) -> Any

    def __setitem__(self, key, value):
        # type: (Any, Any) -> None
        raise NotImplementedError(_not_implemented_message.format(
            type(self), '__setitem__')
        )

    @abc.abstractmethod
    def __iter__(self):
        # type: () -> Iterator
        pass

    @abc.abstractmethod
    def __len__(self):
        # type: () -> int
        pass

    # ------------------------------------------------------------------------
    # Required attributes
    # ------------------------------------------------------------------------
    @property
    def base(self):
        """The base array I am a view of. None by default."""

    @property
    @abc.abstractmethod
    def dtype(self):
        """An instance of 'ExtensionDtype'."""
        # type: () -> ExtensionDtype
        pass

    @property
    def shape(self):
        # type: () -> Tuple[int, ...]
        return (len(self),)

    @property
    def ndim(self):
        # type: () -> int
        """Extension Arrays are only allowed to be 1-dimensional."""
        return 1

    @property
    @abc.abstractmethod
    def nbytes(self):
        """The number of bytes needed to store this object in memory."""
        # type: () -> int
        pass

    # ------------------------------------------------------------------------
    # Additional Methods
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def isna(self):
        """Boolean NumPy array indicating if each value is missing."""
        # type: () -> np.ndarray
        pass

    # ------------------------------------------------------------------------
    # Indexing methods
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def take(self, indexer, allow_fill=True, fill_value=None):
        # type: (Sequence, bool, Optional[Any]) -> ExtensionArray
        """For slicing"""

    def take_nd(self, indexer, allow_fill=True, fill_value=None):
        """For slicing"""
        # TODO: this isn't really nescessary for 1-D
        return self.take(indexer, allow_fill=allow_fill,
                         fill_value=fill_value)

    @abc.abstractmethod
    def copy(self, deep=False):
        # type: (bool) -> ExtensionArray
        """Return a copy of the array."""

    # ------------------------------------------------------------------------
    # Block-related methods
    # ------------------------------------------------------------------------
    @property
    def _fill_value(self):
        """The missing value for this type, e.g. np.nan"""
        # type: () -> Any
        return None

    @abc.abstractmethod
    def _formatting_values(self):
        # type: () -> np.ndarray
        # At the moment, this has to be an array since we use result.dtype
        """An array of values to be printed in, e.g. the Series repr"""

    @classmethod
    @abc.abstractmethod
    def _concat_same_type(cls, to_concat):
        # type: (Sequence[ExtensionArray]) -> ExtensionArray
        """Concatenate multiple array

        Parameters
        ----------
        to_concat : sequence of this type

        Returns
        -------
        ExtensionArray
        """

    @abc.abstractmethod
    def get_values(self):
        # type: () -> np.ndarray
        """Get the underlying values backing your data
        """
        pass

    def _can_hold_na(self):
        """Whether your array can hold missing values. True by default.

        Notes
        -----
        Setting this to false will optimize some operations like fillna.
        """
        # type: () -> bool
        return True

    @property
    def is_sparse(self):
        """Whether your array is sparse. True by default."""
        # type: () -> bool
        return False

    def _slice(self, slicer):
        # type: (Union[tuple, Sequence, int]) -> 'ExtensionArray'
        """Return a new array sliced by `slicer`.

        Parameters
        ----------
        slicer : slice or np.ndarray
            If an array, it should just be a boolean mask

        Returns
        -------
        array : ExtensionArray
            Should return an ExtensionArray, even if ``self[slicer]``
            would return a scalar.
        """
        return type(self)(self[slicer])

    def value_counts(self, dropna=True):
        """Optional method for computing the histogram of the counts.

        Parameters
        ----------
        dropna : bool, default True
            whether to exclude missing values from the computation

        Returns
        -------
        counts : Series
        """
        from pandas.core.algorithms import value_counts
        mask = ~np.asarray(self.isna())
        values = self[mask]  # XXX: this imposes boolean indexing
        return value_counts(np.asarray(values), dropna=dropna)
