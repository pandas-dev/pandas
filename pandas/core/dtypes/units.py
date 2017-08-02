

import pint

import numpy as np

from pandas import compat
from .dtypes import ExtensionDtype

unit_registry = pint.UnitRegistry()

class DimensionedFloatDtypeType(type):
    """
    The type of UnitDtype.
    """
    pass

class DimensionedFloatDtype(ExtensionDtype):

    """
    A dtype for holding float64 nubers with units

    THIS IS NOT A REAL NUMPY DTYPE
    """
    type = DimensionedFloatDtypeType
    _metadata = ['unit']
    _cache = {}
    kind = "f"
    str="f8"
    base = np.dtype('f8')

    def __new__(cls, unit=None):
        """ Create a new unit if needed, otherwise return from the cache

        Parameters
        ----------
        unit : string unit that this represents.
        """
        if unit is None:
            # we are called as an empty constructor
            # generally for pickle compat
            return object.__new__(cls)

        # Assume unit is a string.
        unit_object = getattr(unit_registry, unit) #Raises TypeError if Unit is not a string.

        # set/retrieve from cache
        try:
            return cls._cache[unit]
        except KeyError:
            u = object.__new__(cls)
            u.unit = unit_object
            cls._cache[unit] = u
            return u
    def __hash__(self):
        return hash(str(self))

    def __unicode__(self):
        return "dimensionedFloat[{unit}]".format(unit=str(self.unit))

    def __setstate__(self, state):
        # Use the same unit but from our registry, not the pickled unit
        # Mixing units from different registries causes errors.
        self.unit = getattr(unit_registry, str(state["unit"]))


    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == str(self)
        if not isinstance(other, DimensionedFloatDtype):
            return NotImplemented
        return self.unit==other.unit

    @classmethod
    def construct_from_string(cls, string):
        """ attempt to construct this type from a string, raise a TypeError if
        it's not possible """
        try:
            typename, unit = string.split("[")
            if unit[-1]=="]" and typename == "dimensionedFloat":
                return cls(unit[:-1])
        except:
            pass
        raise TypeError("cannot construct a DimensionedFloatDtype")
