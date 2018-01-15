"""Extend pandas with custom array types"""
import abc

from pandas.compat import add_metaclass


@add_metaclass(abc.ABCMeta)
class ExtensionDtype(object):
    """A custom data type for your array.
    """
    @property
    def type(self):
        """Typically a metaclass inheriting from 'type' with no methods."""
        return type(self.name, (), {})

    @property
    def kind(self):
        """A character code (one of 'biufcmMOSUV'), default 'O'

        See Also
        --------
        numpy.dtype.kind
        """
        return 'O'

    @property
    @abc.abstractmethod
    def name(self):
        """An string identifying the data type.

        Will be used in, e.g. ``Series.dtype``
        """

    @property
    def names(self):
        """Ordered list of field names, or None if there are no fields"""
        return None

    @classmethod
    def construct_from_string(cls, string):
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
