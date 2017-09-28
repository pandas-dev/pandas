# -*- coding: utf-8 -*-
"""

accessor.py contains base classes for implementing accessor properties
that can be mixed into or pinned onto other pandas classes.

"""
from pandas.core.common import AbstractMethodError


class DirNamesMixin(object):
    _accessors = frozenset([])
    _deprecations = frozenset([])

    def _dir_deletions(self):
        """ delete unwanted __dir__ for this object """
        return self._accessors | self._deprecations

    def _dir_additions(self):
        """ add addtional __dir__ for this object """
        rv = set()
        for accessor in self._accessors:
            try:
                getattr(self, accessor)
                rv.add(accessor)
            except AttributeError:
                pass
        return rv

    def __dir__(self):
        """
        Provide method name lookup and completion
        Only provide 'public' methods
        """
        rv = set(dir(type(self)))
        rv = (rv - self._dir_deletions()) | self._dir_additions()
        return sorted(rv)


class AccessorProperty(object):
    """Descriptor for implementing accessor properties like Series.str
    """

    def __init__(self, accessor_cls, construct_accessor=None):
        self.accessor_cls = accessor_cls
        self.construct_accessor = (construct_accessor or
                                   accessor_cls._make_accessor)
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


class PandasDelegate(object):
    """ an abstract base class for delegating methods/properties """

    @classmethod
    def _make_accessor(cls, data):
        raise AbstractMethodError("_make_accessor should be implemented"
                                  "by subclass and return an instance"
                                  "of `cls`.")

    def _delegate_property_get(self, name, *args, **kwargs):
        raise TypeError("You cannot access the "
                        "property {name}".format(name=name))

    def _delegate_property_set(self, name, value, *args, **kwargs):
        raise TypeError("The property {name} cannot be set".format(name=name))

    def _delegate_method(self, name, *args, **kwargs):
        raise TypeError("You cannot call method {name}".format(name=name))

    @classmethod
    def _add_delegate_accessors(cls, delegate, accessors, typ,
                                overwrite=False):
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

            if typ == 'property':
                f = Delegator.create_delegator_property(name, delegate)
            else:
                f = Delegator.create_delegator_method(name, delegate)

            # don't overwrite existing methods/properties
            if overwrite or not hasattr(cls, name):
                setattr(cls, name, f)


class Delegator(object):
    """ Delegator class contains methods that are used by PandasDelegate
    and Accessor subclasses, but that so not ultimately belong in
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
