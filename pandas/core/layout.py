"""
Implements ArrayLayout copy factory to change memory layout
of `numpy.ndarrays`.
Depending on the use case, operations on DataFrames can be much
faster if the appropriate memory layout is set and preserved.

The implementation allows for changing the desired layout. Changes apply when
copies or new objects are created, as for example, when slicing or aggregating
via groupby ...

This implementation tries to solve the issue raised on GitHub
https://github.com/pandas-dev/pandas/issues/26502

"""
import numpy as np


_DEFAULT_MEMORY_LAYOUT = 'F'


class ArrayLayout(object):
    """
    Array layout management for numpy.ndarrays.

    Singleton implementation.

    Example:
    >>> from pandas.core.layout import array_layout
    >>> array_layout.order = 'K'  #
    >>> # K ... keep array layout from input
    >>> # C ... set to c-contiguous / column order
    >>> # F ... set to f-contiguous / row order
    >>> array = array_layout.apply(array)
    >>> array = array_layout.apply(array, 'C')
    >>> array = array_layout.copy(array)
    >>> array = array_layout.apply_on_transpose(array)

    """

    _order = _DEFAULT_MEMORY_LAYOUT
    _instance = None

    @property
    def order(self):
        """
        Return memory layout ordering.

        :return: `str`
        """
        if self.__class__._order is None:
            raise AssertionError("Array layout order not set.")
        return self.__class__._order

    @order.setter
    def order(self, order):
        """
        Set memory layout order.
        Allowed values are 'C', 'F', and 'K'. Raises AssertionError
        when trying to set other values.

        :param order: `str`
        :return: `None`
        """
        assert order in ['C', 'F', 'K'], "Only 'C', 'F' and 'K' supported."
        self.__class__._order = order

    def __new__(cls):
        """
        Create only one instance throughout the lifetime of this process.

        :return: `ArrayLayout` instance as singleton
        """
        return cls._instance if cls._instance is not None \
            else super(ArrayLayout, cls).__new__(ArrayLayout)

    def __init__(self):
        """
        Kind of a singleton implementation for the memory layout order in use.
        Works together with the __new__ definition.

        Returns initialized singleton instance of ArrayLayout

        :return: `ArrayLayout` instance, singleton, the one and only.
        """
        if self._instance is None:
            self._instance = self

    @staticmethod
    def get_from(array):
        """
        Get memory layout from array

        Possible values:
           'C' ... only C-contiguous or column order
           'F' ... only F-contiguous or row order
           'O' ... other: both, C- and F-contiguous or both
           not C- or F-contiguous (as on empty arrays).

        :param array: `numpy.ndarray`
        :return: `str`
        """
        if array.flags.c_contiguous == array.flags.f_contiguous:
            return 'O'
        return {True: 'C', False: 'F'}[array.flags.c_contiguous]

    def apply(self, array, order=None):
        """
        Apply the order set or the order given as input on the array
        given as input.

        Possible values:
           'C' ... apply C-contiguous layout or column order
           'F' ... apply F-contiguous layout or row order
           'K' ... keep the given layout

        :param array: `numpy.ndarray`
        :param order: `str`
        :return: `np.ndarray`
        """
        order = self.__class__._order if order is None else order

        if order == 'K':
            return array

        array_order = ArrayLayout.get_from(array)
        if array_order == order:
            return array

        return np.reshape(np.ravel(array), array.shape, order=order)

    def copy(self, array, order=None):
        """
        Return a copy of the input array with the memory layout set.
        Layout set:
           'C' ... return C-contiguous copy
           'F' ... return F-contiguous copy
           'K' ... return copy with same layout as
           given by the input array.

        :param array: `np.ndarray`
        :return: `np.ndarray`
        """
        order = order if order is not None else self.__class__._order
        return array.copy(order=self.get_from(array)) if order == 'K' \
            else array.copy(order=order)

    def copy_transposed(self, array):
        """
        Return a copy of the input array in order that its transpose
        has the memory layout set.

        Note: numpy simply changes the memory layout from row to column
        order instead of reshuffling the data in memory.

        Layout set:
           'C' ... return F-contiguous copy
           'F' ... return C-contiguous copy
           'K' ... return copy with oposite (C versus F) layout as
           given by the input array.

        :param array: `np.ndarray`
        :return: `np.ndarray`

        :param array:
        :return:
        """
        if self.__class__._order == 'K':
            return array.copy(
                order={'C': 'C', 'F': 'F', 'O': None}[self.get_from(array)])
        else:
            return array.copy(
                order={'C': 'F', 'F': 'C'}[self.__class__._order])

    def __str__(self):
        return str(self.__class__._order)


array_layout = ArrayLayout()  # Singleton
