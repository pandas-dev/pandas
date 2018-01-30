"""An interface for extending pandas with custom arrays."""
import abc

from pandas.compat import add_metaclass


_not_implemented_message = "{} does not implement {}."


@add_metaclass(abc.ABCMeta)
class ExtensionArray(object):
    """Abstract base class for custom 1-D array types.

    Notes
    -----
    pandas will recognize instances of this class as proper arrays
    with a custom type and will not attempt to coerce them to objects.

    ExtensionArrays are limited to 1 dimension.

    They may be backed by none, one, or many NumPy ararys. For example,
    ``pandas.Categorical`` is an extension array backed by two arrays,
    one for codes and one for categories. An array of IPv6 address may
    be backed by a NumPy structured array with two fields, one for the
    lower 64 bits and one for the upper 64 bits. Or they may be backed
    by some other storage type, like Python lists. Pandas makes no
    assumptions on how the data are stored, just that it can be converted
    to a NumPy array.

    Extension arrays should be able to be constructed with instances of
    the class, i.e. ``ExtensionArray(extension_array)`` should return
    an instance, not error.

    Additionally, certain methods and interfaces are required for proper
    this array to be properly stored inside a ``DataFrame`` or ``Series``.
    """
    # ------------------------------------------------------------------------
    # Must be a Sequence
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def __getitem__(self, item):
        # type (Any) -> Any
        """Select a subset of self.

        Parameters
        ----------
        item : int, slice, or ndarray
            * int: The position in 'self' to get.

            * slice: A slice object, where 'start', 'stop', and 'step' are
              integers or None

            * ndarray: A 1-d boolean NumPy ndarray the same length as 'self'

        Returns
        -------
        item : scalar or ExtensionArray

        Notes
        -----
        For scalar ``item``, return a scalar value suitable for the array's
        type. This should be an instance of ``self.dtype.type``.

        For slice ``key``, return an instance of ``ExtensionArray``, even
        if the slice is length 0 or 1.

        For a boolean mask, return an instance of ``ExtensionArray``, filtered
        to the values where ``item`` is True.
        """

    def __setitem__(self, key, value):
        # type: (Any, Any) -> None
        raise NotImplementedError(_not_implemented_message.format(
            type(self), '__setitem__')
        )

    @abc.abstractmethod
    def __len__(self):
        """Length of this array

        Returns
        -------
        length : int
        """
        # type: () -> int
        pass

    # ------------------------------------------------------------------------
    # Required attributes
    # ------------------------------------------------------------------------
    @property
    @abc.abstractmethod
    def dtype(self):
        # type: () -> ExtensionDtype
        """An instance of 'ExtensionDtype'."""

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
        # type: () -> int
        """The number of bytes needed to store this object in memory.

        If this is expensive to compute, return an approximate lower bound
        on the number of bytes needed.
        """

    # ------------------------------------------------------------------------
    # Additional Methods
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def isna(self):
        # type: () -> np.ndarray
        """Boolean NumPy array indicating if each value is missing."""

    # ------------------------------------------------------------------------
    # Indexing methods
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def take(self, indexer, allow_fill=True, fill_value=None):
        # type: (Sequence[int], bool, Optional[Any]) -> ExtensionArray
        """Take elements from an array.

        Parameters
        ----------
        indexer : sequence of integers
            indices to be taken. -1 is used to indicate values
            that are missing.
        allow_fill : bool, default True
            If False, indexer is assumed to contain no -1 values so no filling
            will be done. This short-circuits computation of a mask. Result is
            undefined if allow_fill == False and -1 is present in indexer.
        fill_value : any, default None
            Fill value to replace -1 values with. By default, this uses
            the missing value sentinel for this type, ``self._fill_value``.

        Notes
        -----
        This should follow pandas' semantics where -1 indicates missing values.
        Positions where indexer is ``-1`` should be filled with the missing
        value for this type.

        This is called by ``Series.__getitem__``, ``.loc``, ``iloc``, when the
        indexer is a sequence of values.

        Examples
        --------
        Suppose the extension array somehow backed by a NumPy structured array
        and that the underlying structured array is stored as ``self.data``.
        Then ``take`` may be written as

        .. code-block:: python

           def take(self, indexer, allow_fill=True, fill_value=None):
               mask = indexer == -1
               result = self.data.take(indexer)
               result[mask] = self._fill_value
               return type(self)(result)
        """

    @abc.abstractmethod
    def copy(self, deep=False):
        # type: (bool) -> ExtensionArray
        """Return a copy of the array.

        Parameters
        ----------
        deep : bool, default False
            Also copy the underlying data backing this array.

        Returns
        -------
        ExtensionArray
        """

    # ------------------------------------------------------------------------
    # Block-related methods
    # ------------------------------------------------------------------------
    @property
    def _fill_value(self):
        # type: () -> Any
        """The missing value for this type, e.g. np.nan"""
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

    def _can_hold_na(self):
        # type: () -> bool
        """Whether your array can hold missing values. True by default.

        Notes
        -----
        Setting this to false will optimize some operations like fillna.
        """
        return True
