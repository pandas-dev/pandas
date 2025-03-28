from __future__ import annotations

from functools import cached_property

from dask._expr import Expr
from dask.array.core import T_IntOrNaN
from dask.utils import cached_cumsum


class ArrayExpr(Expr):
    _cached_keys = None

    @cached_property
    def shape(self) -> tuple[T_IntOrNaN, ...]:
        return tuple(cached_cumsum(c, initial_zero=True)[-1] for c in self.chunks)

    @cached_property
    def ndim(self):
        return len(self.shape)

    @cached_property
    def chunksize(self) -> tuple[T_IntOrNaN, ...]:
        return tuple(max(c) for c in self.chunks)

    @cached_property
    def dtype(self):
        if isinstance(self._meta, tuple):
            dtype = self._meta[0].dtype
        else:
            dtype = self._meta.dtype
        return dtype

    @cached_property
    def chunks(self):
        raise NotImplementedError("Subclass must implement 'chunks'")

    @cached_property
    def numblocks(self):
        return tuple(map(len, self.chunks))

    @property
    def name(self):
        return self._name

    def __dask_keys__(self):
        out = self.lower_completely()
        if self._cached_keys is not None:
            return self._cached_keys

        name, chunks, numblocks = out.name, out.chunks, out.numblocks

        def keys(*args):
            if not chunks:
                return [(name,)]
            ind = len(args)
            if ind + 1 == len(numblocks):
                result = [(name,) + args + (i,) for i in range(numblocks[ind])]
            else:
                result = [keys(*(args + (i,))) for i in range(numblocks[ind])]
            return result

        self._cached_keys = result = keys()
        return result

    def __hash__(self):
        return hash(self._name)

    def optimize(self):
        return self.simplify().lower_completely()

    def __getattr__(self, key):
        try:
            return object.__getattribute__(self, key)
        except AttributeError as err:
            if key.startswith("_meta"):
                # Avoid a recursive loop if/when `self._meta*`
                # produces an `AttributeError`
                raise RuntimeError(
                    f"Failed to generate metadata for {self}. "
                    "This operation may not be supported by the current backend."
                )

            # Allow operands to be accessed as attributes
            # as long as the keys are not already reserved
            # by existing methods/properties
            _parameters = type(self)._parameters
            if key in _parameters:
                idx = _parameters.index(key)
                return self.operands[idx]

            raise AttributeError(
                f"{err}\n\n"
                "This often means that you are attempting to use an unsupported "
                f"API function.."
            )
