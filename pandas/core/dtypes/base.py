"""Extend pandas with custom array types"""
from typing import List, Optional, Type

import numpy as np

from pandas.errors import AbstractMethodError

from pandas.core.dtypes.generic import ABCDataFrame, ABCIndexClass, ABCSeries


class _DtypeOpsMixin(object):
    # Not all of pandas' extension dtypes are compatibile with
    # the new ExtensionArray interface. This means PandasExtensionDtype
    # can't subclass ExtensionDtype yet, as is_extension_array_dtype would
    # incorrectly say that these types are extension types.
    #
    # In the interim, we put methods that are shared between the two base
    # classes ExtensionDtype and PandasExtensionDtype here. Both those base
    # classes will inherit from this Mixin. Once everything is compatible, this
    # class's methods can be moved to ExtensionDtype and removed.

    # na_value is the default NA value to use for this type. This is used in
    # e.g. ExtensionArray.take. This should be the user-facing "boxed" version
    # of the NA value, not the physical NA vaalue for storage.
    # e.g. for JSONArray, this is an empty dictionary.
    na_value = np.nan
    _metadata = ()

    def __eq__(self, other):
        """Check whether 'other' is equal to self.

        By default, 'other' is considered equal if either

        * it's a string matching 'self.name'.
        * it's an instance of this type and all of the
          the attributes in ``self._metadata`` are equal between
          `self` and `other`.

        Parameters
        ----------
        other : Any

        Returns
        -------
        bool
        """
        if isinstance(other, str):
            try:
                other = self.construct_from_string(other)
            except TypeError:
                return False
        if isinstance(other, type(self)):
            return all(
                getattr(self, attr) == getattr(other, attr)
                for attr in self._metadata
            )
        return False

    def __hash__(self):
        return hash(tuple(getattr(self, attr) for attr in self._metadata))

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def names(self):
        # type: () -> Optional[List[str]]
        """Ordered list of field names, or None if there are no fields.

        This is for compatibility with NumPy arrays, and may be removed in the
        future.
        """
        return None

    @classmethod
    def is_dtype(cls, dtype):
        """Check if we match 'dtype'.

        Parameters
        ----------
        dtype : object
            The object to check.

        Returns
        -------
        is_dtype : bool

        Notes
        -----
        The default implementation is True if

        1. ``cls.construct_from_string(dtype)`` is an instance
           of ``cls``.
        2. ``dtype`` is an object and is an instance of ``cls``
        3. ``dtype`` has a ``dtype`` attribute, and any of the above
           conditions is true for ``dtype.dtype``.
        """
        dtype = getattr(dtype, 'dtype', dtype)

        if isinstance(dtype, (ABCSeries, ABCIndexClass,
                              ABCDataFrame, np.dtype)):
            # https://github.com/pandas-dev/pandas/issues/22960
            # avoid passing data to `construct_from_string`. This could
            # cause a FutureWarning from numpy about failing elementwise
            # comparison from, e.g., comparing DataFrame == 'category'.
            return False
        elif dtype is None:
            return False
        elif isinstance(dtype, cls):
            return True
        try:
            return cls.construct_from_string(dtype) is not None
        except TypeError:
            return False

    @property
    def _is_numeric(self):
        # type: () -> bool
        """
        Whether columns with this dtype should be considered numeric.

        By default ExtensionDtypes are assumed to be non-numeric.
        They'll be excluded from operations that exclude non-numeric
        columns, like (groupby) reductions, plotting, etc.
        """
        return False

    @property
    def _is_boolean(self):
        # type: () -> bool
        """
        Whether this dtype should be considered boolean.

        By default, ExtensionDtypes are assumed to be non-numeric.
        Setting this to True will affect the behavior of several places,
        e.g.

        * is_bool
        * boolean indexing

        Returns
        -------
        bool
        """
        return False


class ExtensionDtype(_DtypeOpsMixin):
    """
    A custom data type, to be paired with an ExtensionArray.

    .. versionadded:: 0.23.0

    See Also
    --------
    extensions.register_extension_dtype
    extensions.ExtensionArray

    Notes
    -----
    The interface includes the following abstract methods that must
    be implemented by subclasses:

    * type
    * name
    * construct_from_string

    The following attributes influence the behavior of the dtype in
    pandas operations

    * _is_numeric
    * _is_boolean

    Optionally one can override construct_array_type for construction
    with the name of this dtype via the Registry. See
    :meth:`extensions.register_extension_dtype`.

    * construct_array_type

    The `na_value` class attribute can be used to set the default NA value
    for this type. :attr:`numpy.nan` is used by default.

    ExtensionDtypes are required to be hashable. The base class provides
    a default implementation, which relies on the ``_metadata`` class
    attribute. ``_metadata`` should be a tuple containing the strings
    that define your data type. For example, with ``PeriodDtype`` that's
    the ``freq`` attribute.

    **If you have a parametrized dtype you should set the ``_metadata``
    class property**.

    Ideally, the attributes in ``_metadata`` will match the
    parameters to your ``ExtensionDtype.__init__`` (if any). If any of
    the attributes in ``_metadata`` don't implement the standard
    ``__eq__`` or ``__hash__``, the default implementations here will not
    work.

    .. versionchanged:: 0.24.0

       Added ``_metadata``, ``__hash__``, and changed the default definition
       of ``__eq__``.

    This class does not inherit from 'abc.ABCMeta' for performance reasons.
    Methods and properties required by the interface raise
    ``pandas.errors.AbstractMethodError`` and no ``register`` method is
    provided for registering virtual subclasses.
    """

    def __str__(self):
        return self.name

    @property
    def type(self):
        # type: () -> Type
        """
        The scalar type for the array, e.g. ``int``

        It's expected ``ExtensionArray[item]`` returns an instance
        of ``ExtensionDtype.type`` for scalar ``item``, assuming
        that value is valid (not NA). NA values do not need to be
        instances of `type`.
        """
        raise AbstractMethodError(self)

    @property
    def kind(self):
        # type () -> str
        """
        A character code (one of 'biufcmMOSUV'), default 'O'

        This should match the NumPy dtype used when the array is
        converted to an ndarray, which is probably 'O' for object if
        the extension type cannot be represented as a built-in NumPy
        type.

        See Also
        --------
        numpy.dtype.kind
        """
        return 'O'

    @property
    def name(self):
        # type: () -> str
        """
        A string identifying the data type.

        Will be used for display in, e.g. ``Series.dtype``
        """
        raise AbstractMethodError(self)

    @classmethod
    def construct_array_type(cls):
        """
        Return the array type associated with this dtype

        Returns
        -------
        type
        """
        raise NotImplementedError

    @classmethod
    def construct_from_string(cls, string):
        """
        Attempt to construct this type from a string.

        Parameters
        ----------
        string : str

        Returns
        -------
        self : instance of 'cls'

        Raises
        ------
        TypeError
            If a class cannot be constructed from this 'string'.

        Examples
        --------
        If the extension dtype can be constructed without any arguments,
        the following may be an adequate implementation.

        >>> @classmethod
        ... def construct_from_string(cls, string)
        ...     if string == cls.name:
        ...         return cls()
        ...     else:
        ...         raise TypeError("Cannot construct a '{}' from "
        ...                         "'{}'".format(cls, string))
        """
        raise AbstractMethodError(cls)
