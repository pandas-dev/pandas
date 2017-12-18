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
    def type(self) -> type:
        """Typically a metaclass inheriting from 'type' with no methods."""

    @property
    @abc.abstractmethod
    def base(self) -> np.dtype:
        # TODO
        pass

    @property
    def kind(self) -> str:
        """A character code (one of 'biufcmMOSUV'), default 'O'

        See Also
        --------
        numpy.dtype.kind
        """
        return 'O'

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """An string identifying the data type.

        Will be used in, e.g. ``Series.dtype``
        """

    @property
    def names(self) -> T.Optional[T.List[str]]:
        """Ordered list of field names, or None if there are no fields"""
        return None

    @classmethod
    def construct_from_string(cls, string: str) -> 'ExtensionDtype':
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
    def is_dtype(cls, dtype: T.Union[str, np.dtype]) -> bool:
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
        pass

    @property
    @abc.abstractmethod
    def shape(self) -> T.Tuple[int, ...]:
        pass

    @property
    @abc.abstractmethod
    def ndim(self) -> int:
        pass

    @property
    @abc.abstractmethod
    def nbytes(self) -> int:
        # TODO: default impl?
        pass

    # ------------------------------------------------------------------------
    # Additional Methods
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def isna(self) -> T.Sequence[bool]:
        # TODO: narrow this type?
        pass

    # ------------------------------------------------------------------------
    # Indexing methods
    # ------------------------------------------------------------------------
    @abc.abstractmethod
    def take(self, indexer, allow_fill=True, fill_value=None):
        """For slicing"""

    @abc.abstractmethod
    def take_nd(self, indexer, allow_fill=True, fill_value=None):
        """For slicing"""

    @abc.abstractmethod
    def copy(self, deep=False):
        """For slicing"""

    # ------------------------------------------------------------------------
    # Block-related methods
    # ------------------------------------------------------------------------
    @property
    def fill_value(self):
        "TODO"
        return None

    @abc.abstractmethod
    def formatting_values(self) -> np.ndarray:
        """An array of values to be printed in, e.g. the Series repr"""

    @classmethod
    @abc.abstractmethod
    def concat_same_type(to_concat):
        """Concatenate multiple array

        Parameters
        ----------
        to_concat : sequence of this type

        Returns
        -------
        cls
        """

    @abc.abstractmethod
    def get_values(self) -> T.Sequence:
        """Get the underlying values backing your data
        """
        # TODO: This is probably the trickiest method to make work in general.
        # The question is how ndarray-like does this have to look?
        pass

    @abc.abstractmethod
    def to_dense(self):
        pass

    @property
    @abc.abstractmethod
    def can_hold_na(self) -> bool:
        pass

    @property
    def is_sparse(self):
        return False

    def slice(self, slicer):
        return self.get_values()[slicer]
