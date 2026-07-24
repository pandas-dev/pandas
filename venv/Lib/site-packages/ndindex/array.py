import warnings

from .ndindex import NDIndex
from .shapetools import asshape

class ArrayIndex(NDIndex):
    """
    Superclass for array indices

    This class should not be instantiated directly. Rather, use one of its
    subclasses, :class:`~.IntegerArray` or :class:`~.BooleanArray`.

    To subclass this, define the `dtype` attribute, as well as all the usual
    ndindex methods.
    """
    __slots__ = ()

    # Subclasses should redefine this
    dtype = None

    def _typecheck(self, idx, shape=None, _copy=True):
        try:
            from numpy import ndarray, asarray, integer, bool_, empty, intp
        except ImportError: # pragma: no cover
            raise ImportError("NumPy must be installed to create array indices")
        try:
            from numpy import VisibleDeprecationWarning
        except ImportError: # pragma: no cover
            from numpy.exceptions import VisibleDeprecationWarning

        if self.dtype is None:
            raise TypeError("Do not instantiate the superclass ArrayIndex directly")

        if shape is not None:
            if idx != []:
                raise ValueError("The shape argument is only allowed for empty arrays (idx=[])")
            shape = asshape(shape)
            if 0 not in shape:
                raise ValueError("The shape argument must be an empty shape")
            idx = empty(shape, dtype=self.dtype)

        if isinstance(idx, (list, ndarray, bool, integer, int, bool_)):
            # Ignore deprecation warnings for things like [1, []]. These will be
            # filtered out anyway since they produce object arrays.
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',
                                        category=VisibleDeprecationWarning,
                                        message='Creating an ndarray from ragged nested sequences')
                a = asarray(idx)
                if a is idx and _copy:
                    a = a.copy()
                if isinstance(idx, list) and 0 in a.shape:
                    if not _copy:
                        raise ValueError("_copy=False is not allowed with list input")
                    a = a.astype(self.dtype)
            if self.dtype == intp and issubclass(a.dtype.type, integer):
                if a.dtype != self.dtype:
                    if not _copy:
                        raise ValueError("If _copy=False, the input array dtype must already be intp")
                    a = a.astype(self.dtype)
            if a.dtype != self.dtype:
                raise TypeError(f"The input array to {self.__class__.__name__} must have dtype {self.dtype.__name__}, not {a.dtype}")
            a.flags.writeable = False
            return (a,)
        raise TypeError(f"{self.__class__.__name__} must be created with an array with dtype {self.dtype.__name__}")

    # These will allow array == ArrayIndex to give True or False instead of
    # returning an array.
    __array_ufunc__ = None
    def __array_function__(self, func, types, args, kwargs):
        return NotImplemented

    def __array__(self, **kwargs):
        raise TypeError(f"Cannot convert {self.__class__.__name__} to an array. Use .array instead.")


    @property
    def raw(self):
        return self.args[0]

    @property
    def array(self):
        """
        Return the NumPy array of self.

        This is the same as `self.args[0]`.

        >>> from ndindex import IntegerArray, BooleanArray
        >>> IntegerArray([0, 1]).array
        array([0, 1])
        >>> BooleanArray([False, True]).array
        array([False, True])

        """
        return self.args[0]

    @property
    def shape(self):
        """
        Return the shape of the array of self.

        This is the same as `self.array.shape`. Note that this is **not** the
        same as the shape of an array that is indexed by `self`. Use
        :meth:`~.NDIndex.newshape` to get that.

        >>> from ndindex import IntegerArray, BooleanArray
        >>> IntegerArray([[0], [1]]).shape
        (2, 1)
        >>> BooleanArray([[False], [True]]).shape
        (2, 1)

        """
        return self.array.shape

    @property
    def ndim(self):
        """
        Return the number of dimensions of the array of self.

        This is the same as `self.array.ndim`. Note that this is **not** the
        same as the number of dimensions of an array that is indexed by
        `self`. Use `len` on :meth:`~.NDIndex.newshape` to get that.

        >>> from ndindex import IntegerArray, BooleanArray
        >>> IntegerArray([[0], [1]]).ndim
        2
        >>> BooleanArray([[False], [True]]).ndim
        2

        """
        return self.array.ndim

    @property
    def size(self):
        """
        Return the number of elements of the array of self.

        This is the same as `self.array.size`. Note that this is **not** the
        same as the number of elements of an array that is indexed by `self`.
        Use `np.prod` on :meth:`~.NDIndex.newshape` to get that.

        >>> from ndindex import IntegerArray, BooleanArray
        >>> IntegerArray([[0], [1]]).size
        2
        >>> BooleanArray([[False], [True]]).size
        2

        """
        return self.array.size

    # The repr form recreates the object. The str form gives the truncated
    # array string and is explicitly non-valid Python (doesn't have commas).
    def __repr__(self):
        if 0 not in self.shape:
            arg = repr(self.array.tolist())
        else:
            arg = f"[], shape={self.shape}"
        return f"{self.__class__.__name__}({arg})"

    def __str__(self):
        from numpy import array2string

        return (self.__class__.__name__
                + "("
                + array2string(self.array).replace('\n', '')
                + ")")

    def __hash__(self):
        return hash(self.array.tobytes())

    def isvalid(self, shape, _axis=0):
        shape = asshape(shape)
        try:
            # The logic is in _raise_indexerror because the error message uses
            # the additional information that is computed when checking if the
            # array is valid.
            self._raise_indexerror(shape, _axis)
        except IndexError:
            return False
        return True
