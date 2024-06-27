"""Implementation of enumerated types.

This module provides the `Enum` class, which can be used to construct
enumerated types.  Those types are defined by providing an *exhaustive
set or list* of possible, named values for a variable of that type.
Enumerated variables of the same type are usually compared between them
for equality and sometimes for order, but are not usually operated upon.

Enumerated values have an associated *name* and *concrete value*.  Every
name is unique and so are concrete values.  An enumerated variable
always takes the concrete value, not its name.  Usually, the concrete
value is not used directly, and frequently it is entirely irrelevant.
For the same reason, an enumerated variable is not usually compared with
concrete values out of its enumerated type.  For that kind of use,
standard variables and constants are more adequate.

"""


__docformat__ = 'reStructuredText'
"""The format of documentation strings in this module."""


class Enum:
    """Enumerated type.

    Each instance of this class represents an enumerated type. The
    values of the type must be declared
    *exhaustively* and named with
    *strings*, and they might be given explicit
    concrete values, though this is not compulsory. Once the type is
    defined, it can not be modified.

    There are three ways of defining an enumerated type. Each one
    of them corresponds to the type of the only argument in the
    constructor of Enum:

    - *Sequence of names*: each enumerated
      value is named using a string, and its order is determined by
      its position in the sequence; the concrete value is assigned
      automatically::

          >>> boolEnum = Enum(['True', 'False'])

    - *Mapping of names*: each enumerated
      value is named by a string and given an explicit concrete value.
      All of the concrete values must be different, or a
      ValueError will be raised::

          >>> priority = Enum({'red': 20, 'orange': 10, 'green': 0})
          >>> colors = Enum({'red': 1, 'blue': 1})
          Traceback (most recent call last):
          ...
          ValueError: enumerated values contain duplicate concrete values: 1

    - *Enumerated type*: in that case, a copy
      of the original enumerated type is created. Both enumerated
      types are considered equal::

          >>> prio2 = Enum(priority)
          >>> priority == prio2
          True

    Please note that names starting with _ are
    not allowed, since they are reserved for internal usage::

        >>> prio2 = Enum(['_xx'])
        Traceback (most recent call last):
        ...
        ValueError: name of enumerated value can not start with ``_``: '_xx'

    The concrete value of an enumerated value is obtained by
    getting its name as an attribute of the Enum
    instance (see __getattr__()) or as an item (see
    __getitem__()). This allows comparisons between
    enumerated values and assigning them to ordinary Python
    variables::

        >>> redv = priority.red
        >>> redv == priority['red']
        True
        >>> redv > priority.green
        True
        >>> priority.red == priority.orange
        False

    The name of the enumerated value corresponding to a concrete
    value can also be obtained by using the
    __call__() method of the enumerated type. In this
    way you get the symbolic name to use it later with
    __getitem__()::

        >>> priority(redv)
        'red'
        >>> priority.red == priority[priority(priority.red)]
        True

    (If you ask, the __getitem__() method is
    not used for this purpose to avoid ambiguity in the case of using
    strings as concrete values.)

    """

    def __init__(self, enum):
        mydict = self.__dict__

        mydict['_names'] = {}
        mydict['_values'] = {}

        if isinstance(enum, list) or isinstance(enum, tuple):
            for (value, name) in enumerate(enum):  # values become 0, 1, 2...
                self._check_and_set_pair(name, value)
        elif isinstance(enum, dict):
            for (name, value) in enum.items():
                self._check_and_set_pair(name, value)
        elif isinstance(enum, Enum):
            for (name, value) in enum._names.items():
                self._check_and_set_pair(name, value)
        else:
            raise TypeError("""\
enumerations can only be created from \
sequences, mappings and other enumerations""")

    def _check_and_set_pair(self, name, value):
        """Check validity of enumerated value and insert it into type."""

        names = self._names
        values = self._values

        if not isinstance(name, str):
            raise TypeError(
                f"name of enumerated value is not a string: {name!r}")
        if name.startswith('_'):
            raise ValueError(
                "name of enumerated value can not start with ``_``: %r"
                % name)
        # This check is only necessary with a sequence base object.
        if name in names:
            raise ValueError(
                "enumerated values contain duplicate names: %r" % name)
        # This check is only necessary with a mapping base object.
        if value in values:
            raise ValueError(
                "enumerated values contain duplicate concrete values: %r"
                % value)

        names[name] = value
        values[value] = name
        self.__dict__[name] = value

    def __getitem__(self, name):
        """Get the concrete value of the enumerated value with that name.

        The name of the enumerated value must be a string. If there is no value
        with that name in the enumeration, a KeyError is raised.

        Examples
        --------

        Let ``enum`` be an enumerated type defined as:

        >>> enum = Enum({'T0': 0, 'T1': 2, 'T2': 5})

        then:

        >>> enum['T1']
        2
        >>> enum['foo']
        Traceback (most recent call last):
          ...
        KeyError: "no enumerated value with that name: 'foo'"

        """

        try:
            return self._names[name]
        except KeyError:
            raise KeyError(f"no enumerated value with that name: {name!r}")

    def __setitem__(self, name, value):
        """This operation is forbidden."""
        raise IndexError("operation not allowed")

    def __delitem__(self, name):
        """This operation is forbidden."""
        raise IndexError("operation not allowed")

    def __getattr__(self, name):
        """Get the concrete value of the enumerated value with that name.

        The name of the enumerated value must be a string. If there is no value
        with that name in the enumeration, an AttributeError is raised.

        Examples
        --------
        Let ``enum`` be an enumerated type defined as:

        >>> enum = Enum({'T0': 0, 'T1': 2, 'T2': 5})

        then:

        >>> enum.T1
        2
        >>> enum.foo
        Traceback (most recent call last):
          ...
        AttributeError: no enumerated value with that name: 'foo'

        """

        try:
            return self[name]
        except KeyError as ke:
            raise AttributeError(*ke.args)

    def __setattr__(self, name, value):
        """This operation is forbidden."""
        raise AttributeError("operation not allowed")

    def __delattr__(self, name):
        """This operation is forbidden."""
        raise AttributeError("operation not allowed")

    def __contains__(self, name):
        """Is there an enumerated value with that name in the type?

        If the enumerated type has an enumerated value with that name, True is
        returned.  Otherwise, False is returned. The name must be a string.

        This method does *not* check for concrete values matching a value in an
        enumerated type. For that, please use the :meth:`Enum.__call__` method.

        Examples
        --------
        Let ``enum`` be an enumerated type defined as:

        >>> enum = Enum({'T0': 0, 'T1': 2, 'T2': 5})

        then:

        >>> 'T1' in enum
        True
        >>> 'foo' in enum
        False
        >>> 0 in enum
        Traceback (most recent call last):
          ...
        TypeError: name of enumerated value is not a string: 0
        >>> enum.T1 in enum  # Be careful with this!
        Traceback (most recent call last):
          ...
        TypeError: name of enumerated value is not a string: 2

        """

        if not isinstance(name, str):
            raise TypeError(
                f"name of enumerated value is not a string: {name!r}")
        return name in self._names

    def __call__(self, value, *default):
        """Get the name of the enumerated value with that concrete value.

        If there is no value with that concrete value in the enumeration and a
        second argument is given as a default, this is returned. Else, a
        ValueError is raised.

        This method can be used for checking that a concrete value belongs to
        the set of concrete values in an enumerated type.

        Examples
        --------
        Let ``enum`` be an enumerated type defined as:

        >>> enum = Enum({'T0': 0, 'T1': 2, 'T2': 5})

        then:

        >>> enum(5)
        'T2'
        >>> enum(42, None) is None
        True
        >>> enum(42)
        Traceback (most recent call last):
          ...
        ValueError: no enumerated value with that concrete value: 42

        """

        try:
            return self._values[value]
        except KeyError:
            if len(default) > 0:
                return default[0]
            raise ValueError(
                f"no enumerated value with that concrete value: {value!r}")

    def __len__(self):
        """Return the number of enumerated values in the enumerated type.

        Examples
        --------
        >>> len(Enum(['e%d' % i for i in range(10)]))
        10

        """

        return len(self._names)

    def __iter__(self):
        """Iterate over the enumerated values.

        Enumerated values are returned as (name, value) pairs *in no particular
        order*.

        Examples
        --------
        >>> enumvals = {'red': 4, 'green': 2, 'blue': 1}
        >>> enum = Enum(enumvals)
        >>> enumdict = dict([(name, value) for (name, value) in enum])
        >>> enumvals == enumdict
        True

        """

        yield from self._names.items()

    def __eq__(self, other):
        """Is the other enumerated type equivalent to this one?

        Two enumerated types are equivalent if they have exactly the same
        enumerated values (i.e. with the same names and concrete values).

        Examples
        --------

        Let ``enum*`` be enumerated types defined as:

        >>> enum1 = Enum({'T0': 0, 'T1': 2})
        >>> enum2 = Enum(enum1)
        >>> enum3 = Enum({'T1': 2, 'T0': 0})
        >>> enum4 = Enum({'T0': 0, 'T1': 2, 'T2': 5})
        >>> enum5 = Enum({'T0': 0})
        >>> enum6 = Enum({'T0': 10, 'T1': 20})

        then:

        >>> enum1 == enum1
        True
        >>> enum1 == enum2 == enum3
        True
        >>> enum1 == enum4
        False
        >>> enum5 == enum1
        False
        >>> enum1 == enum6
        False

        Comparing enumerated types with other kinds of objects produces
        a false result:

        >>> enum1 == {'T0': 0, 'T1': 2}
        False
        >>> enum1 == ['T0', 'T1']
        False
        >>> enum1 == 2
        False

        """

        if not isinstance(other, Enum):
            return False
        return self._names == other._names

    def __ne__(self, other):
        """Is the `other` enumerated type different from this one?

        Two enumerated types are different if they don't have exactly
        the same enumerated values (i.e. with the same names and
        concrete values).

        Examples
        --------

        Let ``enum*`` be enumerated types defined as:

        >>> enum1 = Enum({'T0': 0, 'T1': 2})
        >>> enum2 = Enum(enum1)
        >>> enum3 = Enum({'T1': 2, 'T0': 0})
        >>> enum4 = Enum({'T0': 0, 'T1': 2, 'T2': 5})
        >>> enum5 = Enum({'T0': 0})
        >>> enum6 = Enum({'T0': 10, 'T1': 20})

        then:

        >>> enum1 != enum1
        False
        >>> enum1 != enum2 != enum3
        False
        >>> enum1 != enum4
        True
        >>> enum5 != enum1
        True
        >>> enum1 != enum6
        True

        """

        return not self.__eq__(other)

    # XXX: API incompatible change for PyTables 3 line
    # Overriding __eq__ blocks inheritance of __hash__ in 3.x
    # def __hash__(self):
    #    return hash((self.__class__, tuple(self._names.items())))
    def __repr__(self):
        """Return the canonical string representation of the enumeration. The
        output of this method can be evaluated to give a new enumeration object
        that will compare equal to this one.

        Examples
        --------
        >>> repr(Enum({'name': 10}))
        "Enum({'name': 10})"

        """

        return 'Enum(%s)' % self._names


def _test():
    import doctest
    return doctest.testmod()


if __name__ == '__main__':
    _test()
