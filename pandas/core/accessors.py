#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

An example/recipe for creating a custom accessor.


The primary use case for accessors is when a Series contains instances
of a particular class and we want to access properties/methods of these
instances in Series form.

Suppose we have a custom State class representing US states:

class State(object):
    def __repr__(self):
        return repr(self.name)

    def __init__(self, name):
        self.name = name
        self._abbrev_dict = {'California': 'CA', 'Alabama': 'AL'}

    @property
    def abbrev(self):
        return self._abbrev_dict[self.name]

    @abbrev.setter
    def abbrev(self, value):
        self._abbrev_dict[self.name] = value

    def fips(self):
        return {'California': 6, 'Alabama': 1}[self.name]


We can construct a series of these objects:

>>> ser = pd.Series([State('Alabama'), State('California')])
>>> ser
0       'Alabama'
1    'California'
dtype: object

We would like direct access to the `abbrev` property and `fips` method.
One option is to access these manually with `apply`:

>>> ser.apply(lambda x: x.fips())
0    1
1    6
dtype: int64

But doing that repeatedly gets old in a hurry, so we decide to make a
custom accessor.  This entails subclassing `PandasDelegate` to specify
what should be accessed and how.

There are four methods that *may* be defined in this subclass, one of which
*must* be defined.  The mandatory method is a classmethod called
`_make_accessor`.  `_make_accessor` is responsible doing any validation on
inputs for the accessor.  In this case, the inputs must be a Series
containing State objects.


class StateDelegate(PandasDelegate):

    def __init__(self, values):
        self.values = values

    @classmethod
    def _make_accessor(cls, data):
        if not isinstance(data, pd.Series):
            raise ValueError('Input must be a Series of States')
        elif not data.apply(lambda x: isinstance(x, State)).all():
            raise ValueError('All entries must be State objects')
        return StateDelegate(data)


With `_make_accessor` defined, we have enough to create the accessor, but
not enough to actually do anything useful with it.  In order to access
*methods* of State objects, we implement `_delegate_method`.
`_delegate_method` calls the underlying method for each object in the
series and wraps these in a new Series.  The simplest version looks like:

    def _delegate_method(self, name, *args, **kwargs):
        state_method = lambda x: getattr(x, name)(*args, **kwargs)
        return self.values.apply(state_method)

Similarly in order to access *properties* of State objects, we need to
implement `_delegate_property_get`:

    def _delegate_property_get(self, name):
        state_property = lambda x: getattr(x, name)
        return self.values.apply(state_property)


On ocassion, we may want to be able to *set* property being accessed.
This is discouraged, but allowed (as long as the class being accessed
allows the property to be set).  Doing so requires implementing
`_delegate_property_set`:

    def _delegate_property_set(self, name, new_values):
        for (obj, val) in zip(self.values, new_values):
            setattr(obj, name, val)


With these implemented, `StateDelegate` knows how to handle methods and
properties.  We just need to tell it what names and properties it is
supposed to handle.  This is done by decorating the `StateDelegate`
class with `pd.accessors.wrap_delegate_names`.  We apply the decorator
once with a list of all the methods the accessor should recognize and
once with a list of all the properties the accessor should recognize.


@wrap_delegate_names(delegate=State,
                     accessors=["fips"],
                     typ="method")
@wrap_delegate_names(delegate=State,
                     accessors=["abbrev"],
                     typ="property")
class StateDelegate(PandasDelegate):
    [...]


We can now pin the `state` accessor to the pd.Series class (we could
alternatively pin it to the pd.Index class with a slightly different
implementation above):

pd.Series.state = accessors.AccessorProperty(StateDelegate)


>>> ser = pd.Series([State('Alabama'), State('California')])
>>> isinstance(ser.state, StateDelegate)
True

>>> ser.state.abbrev
0    AL
1    CA
dtype: object

>>> ser.state.fips()
0    1
1    6

>>> ser.state.abbrev = ['Foo', 'Bar']
>>> ser.state.abbrev
0    Foo
1    Bar
dtype: object



"""
from pandas.core.base import PandasObject
from pandas.core import common as com


class PandasDelegate(PandasObject):
    """ an abstract base class for delegating methods/properties

    Usage: To make a custom accessor, start by subclassing `Delegate`.
    See example in the module-level docstring.

    """

    def __init__(self, values):
        self.values = values
    #    #self._freeze()

    @classmethod
    def _make_accessor(cls, data):  # pragma: no cover
        raise NotImplementedError(
            'It is up to subclasses to implement '
            '_make_accessor.  This does input validation on the object to '
            'which the accessor is being pinned.  '
            'It should return an instance of `cls`.')

    def _delegate_property_get(self, name, *args, **kwargs):
        raise TypeError("You cannot access the "
                        "property {name}".format(name=name))

    def _delegate_property_set(self, name, value, *args, **kwargs):
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

        if callable(name):
            # A function/method was passed directly instead of a name
            # This may also render the `delegate` arg unnecessary.
            func.__name__ = name.__name__  # TODO: is this generally valid?
            func.__doc__ = name.__doc__
        else:
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


        This replaces the older usage in which following a class definition
        we would use `Foo._add_delegate_accessors(...)`.  The motivation
        is that we would like to keep as much of a class's internals inside
        the class definition.  For things that we cannot keep directly
        in the class definition, a decorator is more directly tied to
        the definition than a method call outside the definition.

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

                # Allow for a callable to be passed instead of a name.
                title = com._get_callable_name(name)
                title = title or name
                # don't overwrite existing methods/properties unless
                # specifically told to do so
                if overwrite or not hasattr(cls, title):
                    setattr(cls, title, func)

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
