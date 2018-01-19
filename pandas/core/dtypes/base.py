"""Extend pandas with custom array types"""
import abc

from pandas.compat import add_metaclass


@add_metaclass(abc.ABCMeta)
class ExtensionDtype(object):
    """A custom data type for your array.
    """
    @property
    @abc.abstractmethod
    def type(self):
        # type: () -> type
        """The scalar type for your array, e.g. ``int``

        It's expected ``ExtensionArray[item]`` returns an instance
        of ``ExtensionDtype.type`` for scalar ``item``.
        """

    @property
    def kind(self):
        # type () -> str
        """A character code (one of 'biufcmMOSUV'), default 'O'

        This should match the NumPy dtype used when your array is
        converted to an ndarray, which is probably 'O' for object if
        your extension type cannot be represented as a built-in NumPy
        type.

        See Also
        --------
        numpy.dtype.kind
        """
        return 'O'

    @property
    @abc.abstractmethod
    def name(self):
        # type: () -> str
        """A string identifying the data type.

        Will be used for display in, e.g. ``Series.dtype``
        """

    @property
    def names(self):
        # type: () -> Optional[List[str]]
        """Ordered list of field names, or None if there are no fields."""
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
            # XXX: Better to mandate a ``.from_empty`` classmethod
            # rather than imposing this on the constructor?
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
