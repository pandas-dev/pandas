"""
This module contains an interface that external libraries can use to define
their own dtypes compatible with pandas (but NOT NUMPY).
"""

import numpy as np

from .dtypes import ExtensionDtype

class NumpyDtypeWithMetadataType(type): # Do we need this?
    """
    The type of NumpyDtypeWithMetadata
    """
    pass

class NumpyDtypeWithMetadata(ExtensionDtype):

    """
    An ExtentionDtype for data where the data
    can be stored in a numpy dtype, but the dtype itself
    contains meta-data and may redefine arithmetic operations.

    To properly implement caching behaviour,
    you might have to implement a __new__ method.
    """
    type = NumpyDtypeWithMetadataType
    # What attributes should be stored during pickling?
    # If this is provided, you usually do not have to
    # override __getstate__
    _metadata = []

    def base(self):
        """
        In what numpy-compatible dtype the actual data is stored.

        Example: np.dtype('f8')
        """
        raise NotImplementedError("'base' must be implemented by subclass "
                                  "(probably as class-level variable)")


    @classmethod
    def construct_from_string(cls, string):
        """ attempt to construct this type from a string, raise a TypeError if
        it's not possible """
        raise NotImplementedError("'construct_from_string' must be implemented by subclass.")

    def operation_typecompatible(self, operation_name, other_dtype, is_left=True):
        """
        Is the desired operation possible between this dtype and other_dtype?

        Parameters
        ----------
        opertation_name: The name of the desired operation, e.g. '__eq__'
        other_dtype: The dtype of the other operand
        is_left: If this dtype is on the left-hand side of the binary operation.

        Returns
        -------
        Boolean or NotImplemented
        """
        return False

    def get_operation_wrapper(self):
        """
        This is called by `pandas.ops._Op.get_op` to get an object
        responsible for type-coercion (which should have the same interface as _Op)

        Returns
        -------
        A class implementing the same interface as pandas.ops._Op or None
        It should return None, if the default _Op class should be used.
        """
        return None

    def to_dtype(self, data):
        """
        Convert arbitrary data to this dtype.

        Override this, if you need any additional conversions.

        Parameters
        ----------
        data: array-like

        Returns
        -------
        An numpy array with the same dtype as self.base
        """
        return np.asarray(data, dtype = self.base)

class AlwaysSame(NumpyDtypeWithMetadata):
    """
    This is an example how a library could implement a
    subclass of NumpyDtypeWithMetadata, but is it (except for testing)
    not useful for anything else.
    """
    _metadata = [ "_target_value", "base"]
    def __new__(cls, target_value=None):
        if target_value is None:
            #We are unpickling
            return object.__new__(cls)
        try:
            return cls._cache[target_value]
        except KeyError:
            d = object.__new__(cls)
            d._target_value = target_value
            # In this case, we set the base numpy dtype upon object construction.
            d.base = np.dtype(type(target_value)) #Raises, if target_value is not a simple number
            cls._cache[target_value] = d
            return d

    def __hash__(self):
        return hash(self._target_value)

    def __unicode__(self):
        return "always[{}]".format(repr(self._target_value))

    def __setstate__(self, state):
        try:
            self._target_value = state["_target_value"]
        except KeyError:
            print("state", state)
            raise
        self.base = np.dtype(type(self._target_value))

    def __eq__(self, other):
        if not isinstance(other, AlwaysSame):
            return NotImplemented
        return self._target_value == other._target_value

    def to_dtype(self, data):
        """
        Fill the array with the target value.
        """
        # Since performance is irrelevant for this Test-dtype, we
        # do not try to modify data in-place
        data = np.ones(np.shape(data), dtype=self.base)
        data = data*self._target_value
        return data

    def get_operation_wrapper(self):
        """
        This is called by `pandas.ops._Op.get_op` to get an object
        responsible for type-coercion (which should have the same interface as _Op)

        Returns
        -------
        A class implementing the same interface as pandas.ops._Op or None
        It should return None, if the default _Op class should be used.
        """
        class AlwaysSameOp():
            dtype = None
            fill_value = self._target_value
            def __init__(self, left, right, name, na_op):
                self.left = left
                self.right = right

                self.name = name
                self.na_op = na_op

                # Here, a conversion of left and right to lvalues and rvalues could take place.
                # lvalues must be a type that has the desired operator defined.
                self.lvalues = left
                self.rvalues = right
                return None
            def wrap_results(self, results):
                print("l,r ", type(self.left), type(self.right))
                #  Comparison operators return dtype bool.
                if self.name in ["__eq__", "__lt__", "__gt__", "__ge__", "__le__", "__ne__"]:
                    return results
                #  All other operators return dtype AlwaysSame
                if isinstance(self.left.dtype, AlwaysSame):
                    target_dtype = self.left.dtype
                else:
                    assert isinstance(self.right.dtype, AlwaysSame)
                    target_dtype = self.right.dtype
                return target_dtype.to_dtype(results)
        return AlwaysSameOp

    def operation_typecompatible(self, name, other, is_left=True):
        if isinstance(other, AlwaysSame):
            if other._target_value != self._target_value:
                if type(other) != AlwaysSame:
                    return NotImplemented #Allow
                return False
        return True
