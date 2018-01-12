"""Extend pandas with custom array types.
"""
import abc
import typing as T

import numpy as np


class ExtensionDtype(metaclass=abc.ABCMeta):
    """A custom data type for your array.
    """
    @property
    @abc.abstractmethod
    def type(self):
        # type: () -> T.Any
        """Typically a metaclass inheriting from 'type' with no methods."""

    @property
    @abc.abstractmethod
    def base(self):
        # type: () -> np.dtype
        # TODO: what do we need from this?
        pass

    @property
    def kind(self):
        # type: () -> str
        """A character code (one of 'biufcmMOSUV'), default 'O'

        See Also
        --------
        numpy.dtype.kind
        """
        return 'O'

    @property
    @abc.abstractmethod
    def name(self):
        # type: () -> str
        """An string identifying the data type.

        Will be used in, e.g. ``Series.dtype``
        """

    @property
    def names(self):
        # type: () -> T.Optional[T.List[str]]
        """Ordered list of field names, or None if there are no fields"""
        return None

    @classmethod
    def construct_from_string(cls, string):
        # type: (str) -> ExtensionDtype
        """Attempt to construct this type from a string.

        Parameters
        ----------
        string : str

        Returns
        -------
        self : instance of 'cls'

        Raises
        ------
        TypeError

        Notes
        -----
        The default implementation checks if 'string' matches your
        type's name. If so, it calls your class with no arguments.
        """
        if string == cls.name:
            return cls()
        else:
            raise TypeError("Cannot construct a '{}' from "
                            "'{}'".format(cls, string))

    @classmethod
    def is_dtype(cls, dtype):
        # type: (T.Union[str, type]) -> bool
        """Check if we match 'dtype'

        Parameters
        ----------
        dtype : str or dtype

        Returns
        -------
        is_dtype : bool

        Notes
        -----
        The default implementation is True if

        1. 'dtype' is a string that returns true for
           ``cls.construct_from_string``
        2. 'dtype' is ``cls`` or a subclass of ``cls``.
        """
        if isinstance(dtype, str):
            try:
                return isinstance(cls.construct_from_string(dtype), cls)
            except TypeError:
                return False
        else:
            return issubclass(dtype, cls)


class ExtensionArray(metaclass=abc.ABCMeta):
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
        pass

    @abc.abstractmethod
    def __iter__(self):
        pass

    @abc.abstractmethod
    def __len__(self):
        pass

    # ------------------------------------------------------------------------
    # Required attributes
    # ------------------------------------------------------------------------
    @property
    @abc.abstractmethod
    def dtype(self):
        # type: () -> ExtensionDtype
        pass

    @property
    def shape(self):
        # type: () -> T.Tuple[int, ...]
        return (len(self),)

    @property
    def ndim(self):
        # type: () -> int
        """Extension Arrays are only allowed to be 1-dimensional"""
        return 1

    @property
    @abc.abstractmethod
    def nbytes(self):
        # type: () -> int
        # TODO: default impl?
        pass

    # ------------------------------------------------------------------------
    # Additional Methods
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def isna(self):
        # type: () -> T.Sequence[bool]
        # TODO: narrow this type?
        pass

    # ------------------------------------------------------------------------
    # Indexing methods
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def take(self, indexer, allow_fill=True, fill_value=None):
        # type: (T.Sequence, bool, T.Optional[T.Any]) -> ExtensionArray
        """For slicing"""

    @abc.abstractmethod
    def take_nd(self, indexer, allow_fill=True, fill_value=None):
        """For slicing"""
        # TODO: this isn't nescesary if we only allow 1D (though maybe
        # impelment it).

    @abc.abstractmethod
    def copy(self, deep=False):
        # type: (bool) -> ExtensionArray
        """For slicing"""

    # ------------------------------------------------------------------------
    # Block-related methods
    # ------------------------------------------------------------------------
    @property
    def fill_value(self):
        # type: () -> T.Any
        # TODO
        return None

    @abc.abstractmethod
    def formatting_values(self):
        # type: () -> np.ndarray
        # TODO: must this be an array? Can it be any sequence?
        """An array of values to be printed in, e.g. the Series repr"""

    @classmethod
    @abc.abstractmethod
    def concat_same_type(cls, to_concat):
        # type: (T.Sequence[ExtensionArray]) -> ExtensionArray
        """Concatenate multiple array

        Parameters
        ----------
        to_concat : sequence of this type

        Returns
        -------
        cls
        """

    @abc.abstractmethod
    def get_values(self):
        # type: () -> ExtensionArray
        # TODO: What is the required return value? Sequence? ndarray?, ...?
        """Get the underlying values backing your data
        """
        pass

    @abc.abstractmethod
    def to_dense(self):
        # type: () -> ExtensionArray
        # TODO: this shouldn't be abstract.
        pass

    @property
    @abc.abstractmethod
    def can_hold_na(self):
        # type: () -> bool
        pass

    @property
    def is_sparse(self):
        # type: () -> bool
        return False

    def slice(self, slicer):
        # TODO: is this right?
        # In general, no. Probably just remove it?
        return self.get_values()[slicer]
