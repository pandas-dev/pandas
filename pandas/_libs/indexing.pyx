# -*- coding: utf-8 -*-
import cython


@cython.cclass
class _NDFrameIndexerBase:
    """
    A base class for _NDFrameIndexer for fast instantiation and attribute
    access.
    """
    obj = cython.declare(object, visibility='public')
    name = cython.declare(object, visibility='public')
    _ndim = cython.declare(object, visibility='public')

    def __init__(self, name, obj):
        self.obj = obj
        self.name = name
        self._ndim = None

    @property
    def ndim(self):
        # Delay `ndim` instantiation until required as reading it
        # from `obj` isn't entirely cheap.
        ndim = self._ndim
        if ndim is None:
            ndim = self._ndim = self.obj.ndim
        return ndim
