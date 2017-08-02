import numpy as np

from pandas.core.base import (PandasObject)
from pandas.util._decorators import cache_readonly
from pandas import compat
from pandas.core.common import is_null_slice


class Dimensional(PandasObject):
    """
    """

    __array_priority__ = 10
    _typ = 'dimensional'

    def __init__(self, values, dtype):
        # TODO: Sanitize
        self.values = values
        self.dtype = dtype

    @property
    def _constructor(self):
        return Dimensional

    def copy(self):
        """ Copy constructor. """
        return self._constructor(self.values.copy(), self.dtype)

    def astype(self, dtype, copy=True):
        """
        Coerce this type to another dtype
        """
        return np.array(self, dtype=dtype, copy=copy)

    @cache_readonly
    def ndim(self):
        """Number of dimensions """
        return self.values.ndim

    @cache_readonly
    def size(self):
        """ return the len of myself """
        return len(self)

    @property
    def base(self):
        """ compat, we are always our own object """
        return None

    # for Series/ndarray like compat
    @property
    def shape(self):
        """ Shape of the Categorical.

        For internal compatibility with numpy arrays.

        Returns
        -------
        shape : tuple
        """
        return tuple([len(self.values)])

    def __array__(self, dtype=None):
        """
        The numpy array interface.

        Returns
        -------
        values : numpy array
            A numpy array of either the specified dtype or,
            if dtype==None (default), the same dtype as
            categorical.categories.dtype
        """
        if dtype:
            return np.asarray(self.values, dtype)
        return self.values

    @property
    def T(self):
        return self

    def isna(self):
        raise NotImplementedError
    isnull = isna

    def notna(self):
        """
        Inverse of isna

        Both missing values (-1 in .codes) and NA as a category are detected as
        null.

        Returns
        -------
        a boolean array of whether my values are not null

        See also
        --------
        notna : top-level notna
        notnull : alias of notna
        Categorical.isna : boolean inverse of Categorical.notna

        """
        return ~self.isna()
    notnull = notna

    def put(self, *args, **kwargs):
        """
        Replace specific elements in the Categorical with given values.
        """
        raise NotImplementedError(("'put' is not yet implemented "
                                   "for Categorical"))

    def dropna(self):
        raise NotImplementedError

    def get_values(self):
        """ Return the values.

        For internal compatibility with pandas formatting.

        Returns
        -------
        values : numpy array
            A numpy array of the same dtype as categorical.categories.dtype or
            Index if datetime / periods
        """
        return np.array(self)

    def ravel(self, order='C'):
        """ Return a flattened (numpy) array.

        For internal compatibility with numpy arrays.

        Returns
        -------
        raveled : numpy array
        """
        return np.array(self)

    def view(self):
        """Return a view of myself.

        For internal compatibility with numpy arrays.

        Returns
        -------
        view : Categorical
           Returns `self`!
        """
        return self

    def to_dense(self):
        """Return my 'dense' representation

        For internal compatibility with numpy arrays.

        Returns
        -------
        dense : array
        """
        return np.asarray(self)

    def fillna(self, value=None, method=None, limit=None):
        """ Fill NA/NaN values using the specified method.

        Parameters
        ----------
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        value : scalar
            Value to use to fill holes (e.g. 0)
        limit : int, default None
            (Not implemented yet for Categorical!)
            If method is specified, this is the maximum number of consecutive
            NaN values to forward/backward fill. In other words, if there is
            a gap with more than this number of consecutive NaNs, it will only
            be partially filled. If method is not specified, this is the
            maximum number of entries along the entire axis where NaNs will be
            filled.

        Returns
        -------
        filled : Categorical with NA/NaN filled
        """
        raise NotImplementedError

    def _slice(self, slicer):
        """ Return a slice of myself.

        For internal compatibility with numpy arrays.
        """

        # only allow 1 dimensional slicing, but can
        # in a 2-d case be passd (slice(None),....)
        if isinstance(slicer, tuple) and len(slicer) == 2:
            if not is_null_slice(slicer[0]):
                raise AssertionError("invalid slicing for a 1-ndim "
                                     "categorical")
            slicer = slicer[1]

        return self._constructor(self.values[slicer], self.dtype)

    def __len__(self):
        """The length of this Categorical."""
        return len(self.values)

    def __iter__(self):
        """Returns an Iterator over the values of this Categorical."""
        return iter(self.get_values())

    def _tidy_repr(self, max_vals=10, footer=True):
        """ a short repr displaying only max_vals and an optional (but default
        footer)
        """
        num = max_vals // 2
        head = self[:num]._get_repr(length=False, footer=False)
        tail = self[-(max_vals - num):]._get_repr(length=False, footer=False)

        result = '%s, ..., %s' % (head[:-1], tail[1:])
        if footer:
            result = '%s\n%s' % (result, self._repr_footer())

        return compat.text_type(result)

    def _repr_footer(self):
        return 'Length: %d' % (len(self))

    def _get_repr(self, length=True, na_rep='NaN', footer=True):
        return "Dimensional {}".format(self.__array__())
        # TODO: Implement properly

    def __unicode__(self):
        """ Unicode representation. """
        # TODO: implement
        return self._tidy_repr()

    def __getitem__(self, key):
        """ Return an item. """
        return Dimensional(values=self.values[key], dtype=self.dtype)

    def __setitem__(self, key, value):
        raise NotImplementedError
