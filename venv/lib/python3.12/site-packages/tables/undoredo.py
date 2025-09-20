"""Support for undoing and redoing actions.

Functions:

* undo(file, operation, *args)
* redo(file, operation, *args)
* move_to_shadow(file, path)
* move_from_shadow(file, path)
* attr_to_shadow(file, path, name)
* attr_from_shadow(file, path, name)

Misc variables:

`__docformat__`
    The format of documentation strings in this module.

"""

from typing import Literal, TYPE_CHECKING

from .path import split_path

if TYPE_CHECKING:
    from .file import File

__docformat__ = "reStructuredText"
"""The format of documentation strings in this module."""


def undo(
    file_: "File",
    operation: Literal["ADDATTR", "CREATE", "DELATTR", "MOVE", "REMOVE"],
    *args: str,
) -> None:
    """Undo action."""
    if operation == "CREATE":
        undo_create(file_, args[0])
    elif operation == "REMOVE":
        undo_remove(file_, args[0])
    elif operation == "MOVE":
        undo_move(file_, args[0], args[1])
    elif operation == "ADDATTR":
        undo_add_attr(file_, args[0], args[1])
    elif operation == "DELATTR":
        undo_del_attr(file_, args[0], args[1])
    else:
        raise NotImplementedError(
            "the requested unknown operation %r can "
            "not be undone; please report this to the "
            "authors" % operation
        )


def redo(
    file_: "File",
    operation: Literal["ADDATTR", "CREATE", "DELATTR", "MOVE", "REMOVE"],
    *args: str,
) -> None:
    """Re-do action."""
    if operation == "CREATE":
        redo_create(file_, args[0])
    elif operation == "REMOVE":
        redo_remove(file_, args[0])
    elif operation == "MOVE":
        redo_move(file_, args[0], args[1])
    elif operation == "ADDATTR":
        redo_add_attr(file_, args[0], args[1])
    elif operation == "DELATTR":
        redo_del_attr(file_, args[0], args[1])
    else:
        raise NotImplementedError(
            "the requested unknown operation %r can "
            "not be redone; please report this to the "
            "authors" % operation
        )


def move_to_shadow(file_: "File", path: str) -> None:
    """Move a node to the set of shadowed ones."""
    node = file_._get_node(path)

    (shparent, shname) = file_._shadow_name()
    node._g_move(shparent, shname)


def move_from_shadow(file_: "File", path: str) -> None:
    """Move a node fro the set of shadowe dones back to foreground."""
    (shparent, shname) = file_._shadow_name()
    node = shparent._f_get_child(shname)

    (pname, name) = split_path(path)
    parent = file_._get_node(pname)
    node._g_move(parent, name)


def undo_create(file_: "File", path: str) -> None:
    """Undo create node."""
    move_to_shadow(file_, path)


def redo_create(file_: "File", path: str) -> None:
    """Re-do create node."""
    move_from_shadow(file_, path)


def undo_remove(file_: "File", path: str) -> None:
    """Undo remove node."""
    move_from_shadow(file_, path)


def redo_remove(file_: "File", path: str) -> None:
    """Re-do remove node."""
    move_to_shadow(file_, path)


def undo_move(file_: "File", origpath: str, destpath: str) -> None:
    """Undo move node."""
    (origpname, origname) = split_path(origpath)

    node = file_._get_node(destpath)
    origparent = file_._get_node(origpname)
    node._g_move(origparent, origname)


def redo_move(file_: "File", origpath: str, destpath: str) -> None:
    """Re-do move node."""
    (destpname, destname) = split_path(destpath)

    node = file_._get_node(origpath)
    destparent = file_._get_node(destpname)
    node._g_move(destparent, destname)


def attr_to_shadow(file_: "File", path: str, name: str) -> None:
    """Move an attribute to the shadowed attribute set."""
    node = file_._get_node(path)
    attrs = node._v_attrs
    value = getattr(attrs, name)

    (shparent, shname) = file_._shadow_name()
    shattrs = shparent._v_attrs

    # Set the attribute only if it has not been kept in the shadow.
    # This avoids re-pickling complex attributes on REDO.
    if shname not in shattrs:
        shattrs._g__setattr(shname, value)

    attrs._g__delattr(name)


def attr_from_shadow(file_: "File", path: str, name: str) -> None:
    """Move an attribute from shadowed attribute set to foreground."""
    (shparent, shname) = file_._shadow_name()
    shattrs = shparent._v_attrs
    value = getattr(shattrs, shname)

    node = file_._get_node(path)
    node._v_attrs._g__setattr(name, value)

    # Keeping the attribute in the shadow allows reusing it on Undo/Redo.
    # shattrs._g__delattr(shname)


def undo_add_attr(file_: "File", path: str, name: str) -> None:
    """Undo add attribute."""
    attr_to_shadow(file_, path, name)


def redo_add_attr(file_: "File", path: str, name: str) -> None:
    """Re-do add attribute."""
    attr_from_shadow(file_, path, name)


def undo_del_attr(file_: "File", path: str, name: str) -> None:
    """Undo delete attribute."""
    attr_from_shadow(file_, path, name)


def redo_del_attr(file_: "File", path: str, name: str) -> None:
    """Re-do delete attribute."""
    attr_to_shadow(file_, path, name)
