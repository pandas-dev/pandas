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

from .path import split_path


__docformat__ = 'reStructuredText'
"""The format of documentation strings in this module."""


def undo(file_, operation, *args):
    if operation == 'CREATE':
        undo_create(file_, args[0])
    elif operation == 'REMOVE':
        undo_remove(file_, args[0])
    elif operation == 'MOVE':
        undo_move(file_, args[0], args[1])
    elif operation == 'ADDATTR':
        undo_add_attr(file_, args[0], args[1])
    elif operation == 'DELATTR':
        undo_del_attr(file_, args[0], args[1])
    else:
        raise NotImplementedError("the requested unknown operation %r can "
                                  "not be undone; please report this to the "
                                  "authors" % operation)


def redo(file_, operation, *args):
    if operation == 'CREATE':
        redo_create(file_, args[0])
    elif operation == 'REMOVE':
        redo_remove(file_, args[0])
    elif operation == 'MOVE':
        redo_move(file_, args[0], args[1])
    elif operation == 'ADDATTR':
        redo_add_attr(file_, args[0], args[1])
    elif operation == 'DELATTR':
        redo_del_attr(file_, args[0], args[1])
    else:
        raise NotImplementedError("the requested unknown operation %r can "
                                  "not be redone; please report this to the "
                                  "authors" % operation)


def move_to_shadow(file_, path):
    node = file_._get_node(path)

    (shparent, shname) = file_._shadow_name()
    node._g_move(shparent, shname)


def move_from_shadow(file_, path):
    (shparent, shname) = file_._shadow_name()
    node = shparent._f_get_child(shname)

    (pname, name) = split_path(path)
    parent = file_._get_node(pname)
    node._g_move(parent, name)


def undo_create(file_, path):
    move_to_shadow(file_, path)


def redo_create(file_, path):
    move_from_shadow(file_, path)


def undo_remove(file_, path):
    move_from_shadow(file_, path)


def redo_remove(file_, path):
    move_to_shadow(file_, path)


def undo_move(file_, origpath, destpath):
    (origpname, origname) = split_path(origpath)

    node = file_._get_node(destpath)
    origparent = file_._get_node(origpname)
    node._g_move(origparent, origname)


def redo_move(file_, origpath, destpath):
    (destpname, destname) = split_path(destpath)

    node = file_._get_node(origpath)
    destparent = file_._get_node(destpname)
    node._g_move(destparent, destname)


def attr_to_shadow(file_, path, name):
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


def attr_from_shadow(file_, path, name):
    (shparent, shname) = file_._shadow_name()
    shattrs = shparent._v_attrs
    value = getattr(shattrs, shname)

    node = file_._get_node(path)
    node._v_attrs._g__setattr(name, value)

    # Keeping the attribute in the shadow allows reusing it on Undo/Redo.
    # shattrs._g__delattr(shname)


def undo_add_attr(file_, path, name):
    attr_to_shadow(file_, path, name)


def redo_add_attr(file_, path, name):
    attr_from_shadow(file_, path, name)


def undo_del_attr(file_, path, name):
    attr_from_shadow(file_, path, name)


def redo_del_attr(file_, path, name):
    attr_to_shadow(file_, path, name)
