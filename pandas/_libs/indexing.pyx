cdef class _NDFrameIndexerBase:
    """
    A base class for _NDFrameIndexer for fast instantiation and attribute
    access.
    """
    cdef public object obj, name, _ndim

    def __init__(self, name, obj):
        self.obj = obj
        self.name = name
        self._ndim = None

    @property
    def ndim(self) -> int:
        # Delay `ndim` instantiation until required as reading it
        # from `obj` isn't entirely cheap.
        ndim = self._ndim
        if ndim is None:
            ndim = self._ndim = self.obj.ndim
            if ndim > 2:
                raise ValueError(
                    "NDFrameIndexer does not support NDFrame objects with ndim > 2"
                )
        return ndim
