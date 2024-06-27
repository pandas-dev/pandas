"""Miscellaneous mappings used to avoid circular imports.

Variables:

`class_name_dict`
    Node class name to class object mapping.
`class_id_dict`
    Class identifier to class object mapping.

Misc variables:

`__docformat__`
    The format of documentation strings in this module.

"""


# Important: no modules from PyTables should be imported here
# (but standard modules are OK), since the main reason for this module
# is avoiding circular imports!

__docformat__ = 'reStructuredText'
"""The format of documentation strings in this module."""

class_name_dict = {}
"""Node class name to class object mapping.

This dictionary maps class names (e.g. ``'Group'``) to actual class
objects (e.g. `Group`).  Classes are registered here when they are
defined, and they are not expected to be unregistered (by now), but they
can be replaced when the module that defines them is reloaded.

.. versionchanged:: 3.0
   The *classNameDict* dictionary has been renamed into *class_name_dict*.

"""

class_id_dict = {}
"""Class identifier to class object mapping.

This dictionary maps class identifiers (e.g. ``'GROUP'``) to actual
class objects (e.g. `Group`).  Classes defining a new ``_c_classid``
attribute are registered here when they are defined, and they are not
expected to be unregistered (by now), but they can be replaced when the
module that defines them is reloaded.

.. versionchanged:: 3.0
   The *classIdDict* dictionary has been renamed into *class_id_dict*.

"""

# Deprecated API
classNameDict = class_name_dict
classIdDict = class_id_dict


def get_class_by_name(classname):
    """Get the node class matching the `classname`.

    If the name is not registered, a ``TypeError`` is raised.  The empty
    string and ``None`` are also accepted, and mean the ``Node`` class.

    .. versionadded:: 3.0

    """

    # The empty string is accepted for compatibility
    # with old default arguments.
    if classname is None or classname == '':
        classname = 'Node'

    # Get the class object corresponding to `classname`.
    if classname not in class_name_dict:
        raise TypeError("there is no registered node class named ``%s``"
                        % (classname,))

    return class_name_dict[classname]
