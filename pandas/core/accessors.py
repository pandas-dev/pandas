#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pandas.core.base import PandasObject
from pandas.core.common import AbstractMethodError


class PandasDelegate(PandasObject):
    """ an abstract base class for delegating methods/properties

    Usage: To make a custom accessor, subclass `PandasDelegate`, overriding
    the methods below.  Then decorate this subclass with
    `accessors.wrap_delegate_names` describing the methods and properties
    that should be delegated.

    Examples can be found in:

    pandas.core.accessors.CategoricalAccessor
    pandas.core.indexes.accessors (complicated example)
    pandas.core.indexes.category.CategoricalIndex
    pandas.core.strings.StringMethods
    pandas.tests.test_accessors

    """

    def __init__(self, values):
        """
        The subclassed constructor will generally only be called by
        _make_accessor.  See _make_accessor.__doc__.
        """
        self.values = values

    @classmethod
    def _make_accessor(cls, data):  # pragma: no cover
        """
        _make_accessor should implement any necessary validation on the
        data argument to ensure that the properties/methods being
        accessed will be available.

        _make_accessor should return cls(data).  If necessary, the arguments
        to the constructor can be expanded.  In this case, __init__ will
        need to be overrided as well.

        Parameters
        ----------
        data : the underlying object being accessed, usually Series or Index

        Returns
        -------
        Delegate : instance of PandasDelegate or subclass

        """
        raise AbstractMethodError(
            'It is up to subclasses to implement '
            '_make_accessor.  This does input validation on the object to '
            'which the accessor is being pinned.  '
            'It should return an instance of `cls`.')
        # return cls(data)

    def _delegate_property_get(self, name, *args, **kwargs):
        raise TypeError("You cannot access the "
                        "property {name}".format(name=name))

    def _delegate_property_set(self, name, value, *args, **kwargs):
        """
        Overriding _delegate_property_set is discouraged.  It is generally
        better to directly interact with the underlying data than to
        alter it via the accessor.

        An example that ignores this advice can be found in
        tests.test_accessors.TestVectorizedAccessor
        """
        raise TypeError("The property {name} cannot be set".format(name=name))

    def _delegate_method(self, name, *args, **kwargs):
        raise TypeError("You cannot call method {name}".format(name=name))


class AccessorProperty(object):
    """Descriptor for implementing accessor properties like Series.str
    """

    def __init__(self, accessor_cls, construct_accessor=None):
        self.accessor_cls = accessor_cls

        if construct_accessor is None:
            # accessor_cls._make_accessor must be a classmethod
            construct_accessor = accessor_cls._make_accessor

        self.construct_accessor = construct_accessor
        self.__doc__ = accessor_cls.__doc__

    def __get__(self, instance, owner=None):
        if instance is None:
            # this ensures that Series.str.<method> is well defined
            return self.accessor_cls
        return self.construct_accessor(instance)

    def __set__(self, instance, value):
        raise AttributeError("can't set attribute")

    def __delete__(self, instance):
        raise AttributeError("can't delete attribute")


class Delegator(object):
    """ Delegator class contains methods that are used by PandasDelegate
    and Accesor subclasses, but that so not ultimately belong in
    the namespaces of user-facing classes.

    Many of these methods *could* be module-level functions, but are
    retained as staticmethods for organization purposes.
    """

    @staticmethod
    def create_delegator_property(name, delegate):
        # Note: we really only need the `delegate` here for the docstring

        def _getter(self):
            return self._delegate_property_get(name)

        def _setter(self, new_values):
            return self._delegate_property_set(name, new_values)
            # TODO: not hit in tests; not sure this is something we
            # really want anyway

        _getter.__name__ = name
        _setter.__name__ = name
        _doc = getattr(delegate, name).__doc__
        return property(fget=_getter, fset=_setter, doc=_doc)

    @staticmethod
    def create_delegator_method(name, delegate):
        # Note: we really only need the `delegate` here for the docstring

        def func(self, *args, **kwargs):
            return self._delegate_method(name, *args, **kwargs)

        func.__name__ = name
        func.__doc__ = getattr(delegate, name).__doc__
        return func

    @staticmethod
    def delegate_names(delegate, accessors, typ, overwrite=False):
        """
        delegate_names decorates class definitions, e.g:

        @delegate_names(Categorical, ["categories", "ordered"], "property")
        class CategoricalAccessor(PandasDelegate):

            @classmethod
            def _make_accessor(cls, data):
                [...]


        The motivation is that we would like to keep as much of a class's
        internals inside the class definition.  For things that we cannot
        keep directly in the class definition, a decorator is more directly
        tied to the definition than a method call outside the definition.

        """
        # Note: we really only need the `delegate` here for the docstring

        def add_delegate_accessors(cls):
            """
            add accessors to cls from the delegate class

            Parameters
            ----------
            cls : the class to add the methods/properties to
            delegate : the class to get methods/properties & doc-strings
            acccessors : string list of accessors to add
            typ : 'property' or 'method'
            overwrite : boolean, default False
                overwrite the method/property in the target class if it exists
            """
            for name in accessors:
                if typ == "property":
                    func = Delegator.create_delegator_property(name, delegate)
                else:
                    func = Delegator.create_delegator_method(name, delegate)

                # don't overwrite existing methods/properties unless
                # specifically told to do so
                if overwrite or not hasattr(cls, name):
                    setattr(cls, name, func)

            return cls

        return add_delegate_accessors


wrap_delegate_names = Delegator.delegate_names
# TODO: the `delegate` arg to `wrap_delegate_names` is really only relevant
# for a docstring.  It'd be nice if we didn't require it and could duck-type
# instead.

# TODO: There are 2-3 implementations of `_delegate_method`
# and `_delegate_property` that are common enough that we should consider
# making them the defaults.  First, if the series being accessed has `name`
# method/property:
#
# def _delegate_method(self, name, *args, **kwargs):
#    result = getattr(self.values, name)(*args, **kwargs)
#    return result
#
# def _delegate_property_get(self, name):
#    result = getattr(self.values, name)
#    return result
#
#
# Alternately if the series being accessed does not have this attribute,
# but is a series of objects that do have the attribute:
#
# def _delegate_method(self, name, *args, **kwargs):
#    meth = lambda x: getattr(x, name)(*args, **kwargs)
#    return self.values.apply(meth)
#
# def _delegate_property_get(self, name):
#    prop = lambda x: getattr(x, name)
#    return self.values.apply(prop)
#
#
# `apply` would need to be changed to `map` if self.values is an Index.
#
# The third thing to consider moving into the general case is
# core.strings.StringMethods._wrap_result, which handles a bunch of cases
# for how to wrap delegated outputs.
