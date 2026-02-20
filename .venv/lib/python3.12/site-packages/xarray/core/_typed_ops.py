"""Mixin classes with arithmetic operators."""

# This file was generated using xarray.util.generate_ops. Do not edit manually.

from __future__ import annotations

import operator
from collections.abc import Callable
from typing import TYPE_CHECKING, Any, overload

from xarray.computation import ops
from xarray.core import nputils
from xarray.core.types import (
    DaCompatible,
    DsCompatible,
    DtCompatible,
    Self,
    T_Xarray,
    VarCompatible,
)

if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.datatree import DataTree
    from xarray.core.types import T_DataArray as T_DA


class DataTreeOpsMixin:
    __slots__ = ()

    def _binary_op(
        self, other: DtCompatible, f: Callable, reflexive: bool = False
    ) -> Self:
        raise NotImplementedError

    def __add__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.add)

    def __sub__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.sub)

    def __mul__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.mul)

    def __pow__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.pow)

    def __truediv__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.truediv)

    def __floordiv__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.floordiv)

    def __mod__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.mod)

    def __and__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.and_)

    def __xor__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.xor)

    def __or__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.or_)

    def __lshift__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.lshift)

    def __rshift__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.rshift)

    def __lt__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.lt)

    def __le__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.le)

    def __gt__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.gt)

    def __ge__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.ge)

    def __eq__(self, other: DtCompatible) -> Self:  # type:ignore[override]
        return self._binary_op(other, nputils.array_eq)

    def __ne__(self, other: DtCompatible) -> Self:  # type:ignore[override]
        return self._binary_op(other, nputils.array_ne)

    # When __eq__ is defined but __hash__ is not, then an object is unhashable,
    # and it should be declared as follows:
    __hash__: None  # type:ignore[assignment]

    def __radd__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.add, reflexive=True)

    def __rsub__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.sub, reflexive=True)

    def __rmul__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.mul, reflexive=True)

    def __rpow__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.pow, reflexive=True)

    def __rtruediv__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.truediv, reflexive=True)

    def __rfloordiv__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.floordiv, reflexive=True)

    def __rmod__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.mod, reflexive=True)

    def __rand__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.and_, reflexive=True)

    def __rxor__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.xor, reflexive=True)

    def __ror__(self, other: DtCompatible) -> Self:
        return self._binary_op(other, operator.or_, reflexive=True)

    def _unary_op(self, f: Callable, *args: Any, **kwargs: Any) -> Self:
        raise NotImplementedError

    def __neg__(self) -> Self:
        return self._unary_op(operator.neg)

    def __pos__(self) -> Self:
        return self._unary_op(operator.pos)

    def __abs__(self) -> Self:
        return self._unary_op(operator.abs)

    def __invert__(self) -> Self:
        return self._unary_op(operator.invert)

    def round(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.round_, *args, **kwargs)

    def argsort(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.argsort, *args, **kwargs)

    def conj(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.conj, *args, **kwargs)

    def conjugate(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.conjugate, *args, **kwargs)

    __add__.__doc__ = operator.add.__doc__
    __sub__.__doc__ = operator.sub.__doc__
    __mul__.__doc__ = operator.mul.__doc__
    __pow__.__doc__ = operator.pow.__doc__
    __truediv__.__doc__ = operator.truediv.__doc__
    __floordiv__.__doc__ = operator.floordiv.__doc__
    __mod__.__doc__ = operator.mod.__doc__
    __and__.__doc__ = operator.and_.__doc__
    __xor__.__doc__ = operator.xor.__doc__
    __or__.__doc__ = operator.or_.__doc__
    __lshift__.__doc__ = operator.lshift.__doc__
    __rshift__.__doc__ = operator.rshift.__doc__
    __lt__.__doc__ = operator.lt.__doc__
    __le__.__doc__ = operator.le.__doc__
    __gt__.__doc__ = operator.gt.__doc__
    __ge__.__doc__ = operator.ge.__doc__
    __eq__.__doc__ = nputils.array_eq.__doc__
    __ne__.__doc__ = nputils.array_ne.__doc__
    __radd__.__doc__ = operator.add.__doc__
    __rsub__.__doc__ = operator.sub.__doc__
    __rmul__.__doc__ = operator.mul.__doc__
    __rpow__.__doc__ = operator.pow.__doc__
    __rtruediv__.__doc__ = operator.truediv.__doc__
    __rfloordiv__.__doc__ = operator.floordiv.__doc__
    __rmod__.__doc__ = operator.mod.__doc__
    __rand__.__doc__ = operator.and_.__doc__
    __rxor__.__doc__ = operator.xor.__doc__
    __ror__.__doc__ = operator.or_.__doc__
    __neg__.__doc__ = operator.neg.__doc__
    __pos__.__doc__ = operator.pos.__doc__
    __abs__.__doc__ = operator.abs.__doc__
    __invert__.__doc__ = operator.invert.__doc__
    round.__doc__ = ops.round_.__doc__
    argsort.__doc__ = ops.argsort.__doc__
    conj.__doc__ = ops.conj.__doc__
    conjugate.__doc__ = ops.conjugate.__doc__


class DatasetOpsMixin:
    __slots__ = ()

    def _binary_op(
        self, other: DsCompatible, f: Callable, reflexive: bool = False
    ) -> Self:
        raise NotImplementedError

    @overload
    def __add__(self, other: DataTree) -> DataTree: ...

    @overload
    def __add__(self, other: DsCompatible) -> Self: ...

    def __add__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.add)

    @overload
    def __sub__(self, other: DataTree) -> DataTree: ...

    @overload
    def __sub__(self, other: DsCompatible) -> Self: ...

    def __sub__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.sub)

    @overload
    def __mul__(self, other: DataTree) -> DataTree: ...

    @overload
    def __mul__(self, other: DsCompatible) -> Self: ...

    def __mul__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.mul)

    @overload
    def __pow__(self, other: DataTree) -> DataTree: ...

    @overload
    def __pow__(self, other: DsCompatible) -> Self: ...

    def __pow__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.pow)

    @overload
    def __truediv__(self, other: DataTree) -> DataTree: ...

    @overload
    def __truediv__(self, other: DsCompatible) -> Self: ...

    def __truediv__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.truediv)

    @overload
    def __floordiv__(self, other: DataTree) -> DataTree: ...

    @overload
    def __floordiv__(self, other: DsCompatible) -> Self: ...

    def __floordiv__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.floordiv)

    @overload
    def __mod__(self, other: DataTree) -> DataTree: ...

    @overload
    def __mod__(self, other: DsCompatible) -> Self: ...

    def __mod__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.mod)

    @overload
    def __and__(self, other: DataTree) -> DataTree: ...

    @overload
    def __and__(self, other: DsCompatible) -> Self: ...

    def __and__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.and_)

    @overload
    def __xor__(self, other: DataTree) -> DataTree: ...

    @overload
    def __xor__(self, other: DsCompatible) -> Self: ...

    def __xor__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.xor)

    @overload
    def __or__(self, other: DataTree) -> DataTree: ...

    @overload
    def __or__(self, other: DsCompatible) -> Self: ...

    def __or__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.or_)

    @overload
    def __lshift__(self, other: DataTree) -> DataTree: ...

    @overload
    def __lshift__(self, other: DsCompatible) -> Self: ...

    def __lshift__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.lshift)

    @overload
    def __rshift__(self, other: DataTree) -> DataTree: ...

    @overload
    def __rshift__(self, other: DsCompatible) -> Self: ...

    def __rshift__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.rshift)

    @overload
    def __lt__(self, other: DataTree) -> DataTree: ...

    @overload
    def __lt__(self, other: DsCompatible) -> Self: ...

    def __lt__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.lt)

    @overload
    def __le__(self, other: DataTree) -> DataTree: ...

    @overload
    def __le__(self, other: DsCompatible) -> Self: ...

    def __le__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.le)

    @overload
    def __gt__(self, other: DataTree) -> DataTree: ...

    @overload
    def __gt__(self, other: DsCompatible) -> Self: ...

    def __gt__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.gt)

    @overload
    def __ge__(self, other: DataTree) -> DataTree: ...

    @overload
    def __ge__(self, other: DsCompatible) -> Self: ...

    def __ge__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, operator.ge)

    @overload  # type:ignore[override]
    def __eq__(self, other: DataTree) -> DataTree: ...

    @overload
    def __eq__(self, other: DsCompatible) -> Self: ...

    def __eq__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, nputils.array_eq)

    @overload  # type:ignore[override]
    def __ne__(self, other: DataTree) -> DataTree: ...

    @overload
    def __ne__(self, other: DsCompatible) -> Self: ...

    def __ne__(self, other: DsCompatible) -> Self | DataTree:
        return self._binary_op(other, nputils.array_ne)

    # When __eq__ is defined but __hash__ is not, then an object is unhashable,
    # and it should be declared as follows:
    __hash__: None  # type:ignore[assignment]

    def __radd__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.add, reflexive=True)

    def __rsub__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.sub, reflexive=True)

    def __rmul__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.mul, reflexive=True)

    def __rpow__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.pow, reflexive=True)

    def __rtruediv__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.truediv, reflexive=True)

    def __rfloordiv__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.floordiv, reflexive=True)

    def __rmod__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.mod, reflexive=True)

    def __rand__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.and_, reflexive=True)

    def __rxor__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.xor, reflexive=True)

    def __ror__(self, other: DsCompatible) -> Self:
        return self._binary_op(other, operator.or_, reflexive=True)

    def _inplace_binary_op(self, other: DsCompatible, f: Callable) -> Self:
        raise NotImplementedError

    def __iadd__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.iadd)

    def __isub__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.isub)

    def __imul__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.imul)

    def __ipow__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ipow)

    def __itruediv__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.itruediv)

    def __ifloordiv__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ifloordiv)

    def __imod__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.imod)

    def __iand__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.iand)

    def __ixor__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ixor)

    def __ior__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ior)

    def __ilshift__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ilshift)

    def __irshift__(self, other: DsCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.irshift)

    def _unary_op(self, f: Callable, *args: Any, **kwargs: Any) -> Self:
        raise NotImplementedError

    def __neg__(self) -> Self:
        return self._unary_op(operator.neg)

    def __pos__(self) -> Self:
        return self._unary_op(operator.pos)

    def __abs__(self) -> Self:
        return self._unary_op(operator.abs)

    def __invert__(self) -> Self:
        return self._unary_op(operator.invert)

    def round(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.round_, *args, **kwargs)

    def argsort(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.argsort, *args, **kwargs)

    def conj(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.conj, *args, **kwargs)

    def conjugate(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.conjugate, *args, **kwargs)

    __add__.__doc__ = operator.add.__doc__
    __sub__.__doc__ = operator.sub.__doc__
    __mul__.__doc__ = operator.mul.__doc__
    __pow__.__doc__ = operator.pow.__doc__
    __truediv__.__doc__ = operator.truediv.__doc__
    __floordiv__.__doc__ = operator.floordiv.__doc__
    __mod__.__doc__ = operator.mod.__doc__
    __and__.__doc__ = operator.and_.__doc__
    __xor__.__doc__ = operator.xor.__doc__
    __or__.__doc__ = operator.or_.__doc__
    __lshift__.__doc__ = operator.lshift.__doc__
    __rshift__.__doc__ = operator.rshift.__doc__
    __lt__.__doc__ = operator.lt.__doc__
    __le__.__doc__ = operator.le.__doc__
    __gt__.__doc__ = operator.gt.__doc__
    __ge__.__doc__ = operator.ge.__doc__
    __eq__.__doc__ = nputils.array_eq.__doc__
    __ne__.__doc__ = nputils.array_ne.__doc__
    __radd__.__doc__ = operator.add.__doc__
    __rsub__.__doc__ = operator.sub.__doc__
    __rmul__.__doc__ = operator.mul.__doc__
    __rpow__.__doc__ = operator.pow.__doc__
    __rtruediv__.__doc__ = operator.truediv.__doc__
    __rfloordiv__.__doc__ = operator.floordiv.__doc__
    __rmod__.__doc__ = operator.mod.__doc__
    __rand__.__doc__ = operator.and_.__doc__
    __rxor__.__doc__ = operator.xor.__doc__
    __ror__.__doc__ = operator.or_.__doc__
    __iadd__.__doc__ = operator.iadd.__doc__
    __isub__.__doc__ = operator.isub.__doc__
    __imul__.__doc__ = operator.imul.__doc__
    __ipow__.__doc__ = operator.ipow.__doc__
    __itruediv__.__doc__ = operator.itruediv.__doc__
    __ifloordiv__.__doc__ = operator.ifloordiv.__doc__
    __imod__.__doc__ = operator.imod.__doc__
    __iand__.__doc__ = operator.iand.__doc__
    __ixor__.__doc__ = operator.ixor.__doc__
    __ior__.__doc__ = operator.ior.__doc__
    __ilshift__.__doc__ = operator.ilshift.__doc__
    __irshift__.__doc__ = operator.irshift.__doc__
    __neg__.__doc__ = operator.neg.__doc__
    __pos__.__doc__ = operator.pos.__doc__
    __abs__.__doc__ = operator.abs.__doc__
    __invert__.__doc__ = operator.invert.__doc__
    round.__doc__ = ops.round_.__doc__
    argsort.__doc__ = ops.argsort.__doc__
    conj.__doc__ = ops.conj.__doc__
    conjugate.__doc__ = ops.conjugate.__doc__


class DataArrayOpsMixin:
    __slots__ = ()

    def _binary_op(
        self, other: DaCompatible, f: Callable, reflexive: bool = False
    ) -> Self:
        raise NotImplementedError

    @overload
    def __add__(self, other: Dataset) -> Dataset: ...

    @overload
    def __add__(self, other: DataTree) -> DataTree: ...

    @overload
    def __add__(self, other: DaCompatible) -> Self: ...

    def __add__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.add)

    @overload
    def __sub__(self, other: Dataset) -> Dataset: ...

    @overload
    def __sub__(self, other: DataTree) -> DataTree: ...

    @overload
    def __sub__(self, other: DaCompatible) -> Self: ...

    def __sub__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.sub)

    @overload
    def __mul__(self, other: Dataset) -> Dataset: ...

    @overload
    def __mul__(self, other: DataTree) -> DataTree: ...

    @overload
    def __mul__(self, other: DaCompatible) -> Self: ...

    def __mul__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.mul)

    @overload
    def __pow__(self, other: Dataset) -> Dataset: ...

    @overload
    def __pow__(self, other: DataTree) -> DataTree: ...

    @overload
    def __pow__(self, other: DaCompatible) -> Self: ...

    def __pow__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.pow)

    @overload
    def __truediv__(self, other: Dataset) -> Dataset: ...

    @overload
    def __truediv__(self, other: DataTree) -> DataTree: ...

    @overload
    def __truediv__(self, other: DaCompatible) -> Self: ...

    def __truediv__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.truediv)

    @overload
    def __floordiv__(self, other: Dataset) -> Dataset: ...

    @overload
    def __floordiv__(self, other: DataTree) -> DataTree: ...

    @overload
    def __floordiv__(self, other: DaCompatible) -> Self: ...

    def __floordiv__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.floordiv)

    @overload
    def __mod__(self, other: Dataset) -> Dataset: ...

    @overload
    def __mod__(self, other: DataTree) -> DataTree: ...

    @overload
    def __mod__(self, other: DaCompatible) -> Self: ...

    def __mod__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.mod)

    @overload
    def __and__(self, other: Dataset) -> Dataset: ...

    @overload
    def __and__(self, other: DataTree) -> DataTree: ...

    @overload
    def __and__(self, other: DaCompatible) -> Self: ...

    def __and__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.and_)

    @overload
    def __xor__(self, other: Dataset) -> Dataset: ...

    @overload
    def __xor__(self, other: DataTree) -> DataTree: ...

    @overload
    def __xor__(self, other: DaCompatible) -> Self: ...

    def __xor__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.xor)

    @overload
    def __or__(self, other: Dataset) -> Dataset: ...

    @overload
    def __or__(self, other: DataTree) -> DataTree: ...

    @overload
    def __or__(self, other: DaCompatible) -> Self: ...

    def __or__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.or_)

    @overload
    def __lshift__(self, other: Dataset) -> Dataset: ...

    @overload
    def __lshift__(self, other: DataTree) -> DataTree: ...

    @overload
    def __lshift__(self, other: DaCompatible) -> Self: ...

    def __lshift__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.lshift)

    @overload
    def __rshift__(self, other: Dataset) -> Dataset: ...

    @overload
    def __rshift__(self, other: DataTree) -> DataTree: ...

    @overload
    def __rshift__(self, other: DaCompatible) -> Self: ...

    def __rshift__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.rshift)

    @overload
    def __lt__(self, other: Dataset) -> Dataset: ...

    @overload
    def __lt__(self, other: DataTree) -> DataTree: ...

    @overload
    def __lt__(self, other: DaCompatible) -> Self: ...

    def __lt__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.lt)

    @overload
    def __le__(self, other: Dataset) -> Dataset: ...

    @overload
    def __le__(self, other: DataTree) -> DataTree: ...

    @overload
    def __le__(self, other: DaCompatible) -> Self: ...

    def __le__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.le)

    @overload
    def __gt__(self, other: Dataset) -> Dataset: ...

    @overload
    def __gt__(self, other: DataTree) -> DataTree: ...

    @overload
    def __gt__(self, other: DaCompatible) -> Self: ...

    def __gt__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.gt)

    @overload
    def __ge__(self, other: Dataset) -> Dataset: ...

    @overload
    def __ge__(self, other: DataTree) -> DataTree: ...

    @overload
    def __ge__(self, other: DaCompatible) -> Self: ...

    def __ge__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, operator.ge)

    @overload  # type:ignore[override]
    def __eq__(self, other: Dataset) -> Dataset: ...

    @overload
    def __eq__(self, other: DataTree) -> DataTree: ...

    @overload
    def __eq__(self, other: DaCompatible) -> Self: ...

    def __eq__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, nputils.array_eq)

    @overload  # type:ignore[override]
    def __ne__(self, other: Dataset) -> Dataset: ...

    @overload
    def __ne__(self, other: DataTree) -> DataTree: ...

    @overload
    def __ne__(self, other: DaCompatible) -> Self: ...

    def __ne__(self, other: DaCompatible) -> Self | Dataset | DataTree:
        return self._binary_op(other, nputils.array_ne)

    # When __eq__ is defined but __hash__ is not, then an object is unhashable,
    # and it should be declared as follows:
    __hash__: None  # type:ignore[assignment]

    def __radd__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.add, reflexive=True)

    def __rsub__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.sub, reflexive=True)

    def __rmul__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.mul, reflexive=True)

    def __rpow__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.pow, reflexive=True)

    def __rtruediv__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.truediv, reflexive=True)

    def __rfloordiv__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.floordiv, reflexive=True)

    def __rmod__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.mod, reflexive=True)

    def __rand__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.and_, reflexive=True)

    def __rxor__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.xor, reflexive=True)

    def __ror__(self, other: DaCompatible) -> Self:
        return self._binary_op(other, operator.or_, reflexive=True)

    def _inplace_binary_op(self, other: DaCompatible, f: Callable) -> Self:
        raise NotImplementedError

    def __iadd__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.iadd)

    def __isub__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.isub)

    def __imul__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.imul)

    def __ipow__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ipow)

    def __itruediv__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.itruediv)

    def __ifloordiv__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ifloordiv)

    def __imod__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.imod)

    def __iand__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.iand)

    def __ixor__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ixor)

    def __ior__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ior)

    def __ilshift__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ilshift)

    def __irshift__(self, other: DaCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.irshift)

    def _unary_op(self, f: Callable, *args: Any, **kwargs: Any) -> Self:
        raise NotImplementedError

    def __neg__(self) -> Self:
        return self._unary_op(operator.neg)

    def __pos__(self) -> Self:
        return self._unary_op(operator.pos)

    def __abs__(self) -> Self:
        return self._unary_op(operator.abs)

    def __invert__(self) -> Self:
        return self._unary_op(operator.invert)

    def round(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.round_, *args, **kwargs)

    def argsort(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.argsort, *args, **kwargs)

    def conj(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.conj, *args, **kwargs)

    def conjugate(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.conjugate, *args, **kwargs)

    __add__.__doc__ = operator.add.__doc__
    __sub__.__doc__ = operator.sub.__doc__
    __mul__.__doc__ = operator.mul.__doc__
    __pow__.__doc__ = operator.pow.__doc__
    __truediv__.__doc__ = operator.truediv.__doc__
    __floordiv__.__doc__ = operator.floordiv.__doc__
    __mod__.__doc__ = operator.mod.__doc__
    __and__.__doc__ = operator.and_.__doc__
    __xor__.__doc__ = operator.xor.__doc__
    __or__.__doc__ = operator.or_.__doc__
    __lshift__.__doc__ = operator.lshift.__doc__
    __rshift__.__doc__ = operator.rshift.__doc__
    __lt__.__doc__ = operator.lt.__doc__
    __le__.__doc__ = operator.le.__doc__
    __gt__.__doc__ = operator.gt.__doc__
    __ge__.__doc__ = operator.ge.__doc__
    __eq__.__doc__ = nputils.array_eq.__doc__
    __ne__.__doc__ = nputils.array_ne.__doc__
    __radd__.__doc__ = operator.add.__doc__
    __rsub__.__doc__ = operator.sub.__doc__
    __rmul__.__doc__ = operator.mul.__doc__
    __rpow__.__doc__ = operator.pow.__doc__
    __rtruediv__.__doc__ = operator.truediv.__doc__
    __rfloordiv__.__doc__ = operator.floordiv.__doc__
    __rmod__.__doc__ = operator.mod.__doc__
    __rand__.__doc__ = operator.and_.__doc__
    __rxor__.__doc__ = operator.xor.__doc__
    __ror__.__doc__ = operator.or_.__doc__
    __iadd__.__doc__ = operator.iadd.__doc__
    __isub__.__doc__ = operator.isub.__doc__
    __imul__.__doc__ = operator.imul.__doc__
    __ipow__.__doc__ = operator.ipow.__doc__
    __itruediv__.__doc__ = operator.itruediv.__doc__
    __ifloordiv__.__doc__ = operator.ifloordiv.__doc__
    __imod__.__doc__ = operator.imod.__doc__
    __iand__.__doc__ = operator.iand.__doc__
    __ixor__.__doc__ = operator.ixor.__doc__
    __ior__.__doc__ = operator.ior.__doc__
    __ilshift__.__doc__ = operator.ilshift.__doc__
    __irshift__.__doc__ = operator.irshift.__doc__
    __neg__.__doc__ = operator.neg.__doc__
    __pos__.__doc__ = operator.pos.__doc__
    __abs__.__doc__ = operator.abs.__doc__
    __invert__.__doc__ = operator.invert.__doc__
    round.__doc__ = ops.round_.__doc__
    argsort.__doc__ = ops.argsort.__doc__
    conj.__doc__ = ops.conj.__doc__
    conjugate.__doc__ = ops.conjugate.__doc__


class VariableOpsMixin:
    __slots__ = ()

    def _binary_op(
        self, other: VarCompatible, f: Callable, reflexive: bool = False
    ) -> Self:
        raise NotImplementedError

    @overload
    def __add__(self, other: T_DA) -> T_DA: ...

    @overload
    def __add__(self, other: Dataset) -> Dataset: ...

    @overload
    def __add__(self, other: DataTree) -> DataTree: ...

    @overload
    def __add__(self, other: VarCompatible) -> Self: ...

    def __add__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.add)

    @overload
    def __sub__(self, other: T_DA) -> T_DA: ...

    @overload
    def __sub__(self, other: Dataset) -> Dataset: ...

    @overload
    def __sub__(self, other: DataTree) -> DataTree: ...

    @overload
    def __sub__(self, other: VarCompatible) -> Self: ...

    def __sub__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.sub)

    @overload
    def __mul__(self, other: T_DA) -> T_DA: ...

    @overload
    def __mul__(self, other: Dataset) -> Dataset: ...

    @overload
    def __mul__(self, other: DataTree) -> DataTree: ...

    @overload
    def __mul__(self, other: VarCompatible) -> Self: ...

    def __mul__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.mul)

    @overload
    def __pow__(self, other: T_DA) -> T_DA: ...

    @overload
    def __pow__(self, other: Dataset) -> Dataset: ...

    @overload
    def __pow__(self, other: DataTree) -> DataTree: ...

    @overload
    def __pow__(self, other: VarCompatible) -> Self: ...

    def __pow__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.pow)

    @overload
    def __truediv__(self, other: T_DA) -> T_DA: ...

    @overload
    def __truediv__(self, other: Dataset) -> Dataset: ...

    @overload
    def __truediv__(self, other: DataTree) -> DataTree: ...

    @overload
    def __truediv__(self, other: VarCompatible) -> Self: ...

    def __truediv__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.truediv)

    @overload
    def __floordiv__(self, other: T_DA) -> T_DA: ...

    @overload
    def __floordiv__(self, other: Dataset) -> Dataset: ...

    @overload
    def __floordiv__(self, other: DataTree) -> DataTree: ...

    @overload
    def __floordiv__(self, other: VarCompatible) -> Self: ...

    def __floordiv__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.floordiv)

    @overload
    def __mod__(self, other: T_DA) -> T_DA: ...

    @overload
    def __mod__(self, other: Dataset) -> Dataset: ...

    @overload
    def __mod__(self, other: DataTree) -> DataTree: ...

    @overload
    def __mod__(self, other: VarCompatible) -> Self: ...

    def __mod__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.mod)

    @overload
    def __and__(self, other: T_DA) -> T_DA: ...

    @overload
    def __and__(self, other: Dataset) -> Dataset: ...

    @overload
    def __and__(self, other: DataTree) -> DataTree: ...

    @overload
    def __and__(self, other: VarCompatible) -> Self: ...

    def __and__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.and_)

    @overload
    def __xor__(self, other: T_DA) -> T_DA: ...

    @overload
    def __xor__(self, other: Dataset) -> Dataset: ...

    @overload
    def __xor__(self, other: DataTree) -> DataTree: ...

    @overload
    def __xor__(self, other: VarCompatible) -> Self: ...

    def __xor__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.xor)

    @overload
    def __or__(self, other: T_DA) -> T_DA: ...

    @overload
    def __or__(self, other: Dataset) -> Dataset: ...

    @overload
    def __or__(self, other: DataTree) -> DataTree: ...

    @overload
    def __or__(self, other: VarCompatible) -> Self: ...

    def __or__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.or_)

    @overload
    def __lshift__(self, other: T_DA) -> T_DA: ...

    @overload
    def __lshift__(self, other: Dataset) -> Dataset: ...

    @overload
    def __lshift__(self, other: DataTree) -> DataTree: ...

    @overload
    def __lshift__(self, other: VarCompatible) -> Self: ...

    def __lshift__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.lshift)

    @overload
    def __rshift__(self, other: T_DA) -> T_DA: ...

    @overload
    def __rshift__(self, other: Dataset) -> Dataset: ...

    @overload
    def __rshift__(self, other: DataTree) -> DataTree: ...

    @overload
    def __rshift__(self, other: VarCompatible) -> Self: ...

    def __rshift__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.rshift)

    @overload
    def __lt__(self, other: T_DA) -> T_DA: ...

    @overload
    def __lt__(self, other: Dataset) -> Dataset: ...

    @overload
    def __lt__(self, other: DataTree) -> DataTree: ...

    @overload
    def __lt__(self, other: VarCompatible) -> Self: ...

    def __lt__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.lt)

    @overload
    def __le__(self, other: T_DA) -> T_DA: ...

    @overload
    def __le__(self, other: Dataset) -> Dataset: ...

    @overload
    def __le__(self, other: DataTree) -> DataTree: ...

    @overload
    def __le__(self, other: VarCompatible) -> Self: ...

    def __le__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.le)

    @overload
    def __gt__(self, other: T_DA) -> T_DA: ...

    @overload
    def __gt__(self, other: Dataset) -> Dataset: ...

    @overload
    def __gt__(self, other: DataTree) -> DataTree: ...

    @overload
    def __gt__(self, other: VarCompatible) -> Self: ...

    def __gt__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.gt)

    @overload
    def __ge__(self, other: T_DA) -> T_DA: ...

    @overload
    def __ge__(self, other: Dataset) -> Dataset: ...

    @overload
    def __ge__(self, other: DataTree) -> DataTree: ...

    @overload
    def __ge__(self, other: VarCompatible) -> Self: ...

    def __ge__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, operator.ge)

    @overload  # type:ignore[override]
    def __eq__(self, other: T_DA) -> T_DA: ...

    @overload
    def __eq__(self, other: Dataset) -> Dataset: ...

    @overload
    def __eq__(self, other: DataTree) -> DataTree: ...

    @overload
    def __eq__(self, other: VarCompatible) -> Self: ...

    def __eq__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, nputils.array_eq)

    @overload  # type:ignore[override]
    def __ne__(self, other: T_DA) -> T_DA: ...

    @overload
    def __ne__(self, other: Dataset) -> Dataset: ...

    @overload
    def __ne__(self, other: DataTree) -> DataTree: ...

    @overload
    def __ne__(self, other: VarCompatible) -> Self: ...

    def __ne__(self, other: VarCompatible) -> Self | T_DA | Dataset | DataTree:
        return self._binary_op(other, nputils.array_ne)

    # When __eq__ is defined but __hash__ is not, then an object is unhashable,
    # and it should be declared as follows:
    __hash__: None  # type:ignore[assignment]

    def __radd__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.add, reflexive=True)

    def __rsub__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.sub, reflexive=True)

    def __rmul__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.mul, reflexive=True)

    def __rpow__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.pow, reflexive=True)

    def __rtruediv__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.truediv, reflexive=True)

    def __rfloordiv__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.floordiv, reflexive=True)

    def __rmod__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.mod, reflexive=True)

    def __rand__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.and_, reflexive=True)

    def __rxor__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.xor, reflexive=True)

    def __ror__(self, other: VarCompatible) -> Self:
        return self._binary_op(other, operator.or_, reflexive=True)

    def _inplace_binary_op(self, other: VarCompatible, f: Callable) -> Self:
        raise NotImplementedError

    def __iadd__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.iadd)

    def __isub__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.isub)

    def __imul__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.imul)

    def __ipow__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ipow)

    def __itruediv__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.itruediv)

    def __ifloordiv__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ifloordiv)

    def __imod__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.imod)

    def __iand__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.iand)

    def __ixor__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ixor)

    def __ior__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ior)

    def __ilshift__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.ilshift)

    def __irshift__(self, other: VarCompatible) -> Self:  # type:ignore[misc]
        return self._inplace_binary_op(other, operator.irshift)

    def _unary_op(self, f: Callable, *args: Any, **kwargs: Any) -> Self:
        raise NotImplementedError

    def __neg__(self) -> Self:
        return self._unary_op(operator.neg)

    def __pos__(self) -> Self:
        return self._unary_op(operator.pos)

    def __abs__(self) -> Self:
        return self._unary_op(operator.abs)

    def __invert__(self) -> Self:
        return self._unary_op(operator.invert)

    def round(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.round_, *args, **kwargs)

    def argsort(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.argsort, *args, **kwargs)

    def conj(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.conj, *args, **kwargs)

    def conjugate(self, *args: Any, **kwargs: Any) -> Self:
        return self._unary_op(ops.conjugate, *args, **kwargs)

    __add__.__doc__ = operator.add.__doc__
    __sub__.__doc__ = operator.sub.__doc__
    __mul__.__doc__ = operator.mul.__doc__
    __pow__.__doc__ = operator.pow.__doc__
    __truediv__.__doc__ = operator.truediv.__doc__
    __floordiv__.__doc__ = operator.floordiv.__doc__
    __mod__.__doc__ = operator.mod.__doc__
    __and__.__doc__ = operator.and_.__doc__
    __xor__.__doc__ = operator.xor.__doc__
    __or__.__doc__ = operator.or_.__doc__
    __lshift__.__doc__ = operator.lshift.__doc__
    __rshift__.__doc__ = operator.rshift.__doc__
    __lt__.__doc__ = operator.lt.__doc__
    __le__.__doc__ = operator.le.__doc__
    __gt__.__doc__ = operator.gt.__doc__
    __ge__.__doc__ = operator.ge.__doc__
    __eq__.__doc__ = nputils.array_eq.__doc__
    __ne__.__doc__ = nputils.array_ne.__doc__
    __radd__.__doc__ = operator.add.__doc__
    __rsub__.__doc__ = operator.sub.__doc__
    __rmul__.__doc__ = operator.mul.__doc__
    __rpow__.__doc__ = operator.pow.__doc__
    __rtruediv__.__doc__ = operator.truediv.__doc__
    __rfloordiv__.__doc__ = operator.floordiv.__doc__
    __rmod__.__doc__ = operator.mod.__doc__
    __rand__.__doc__ = operator.and_.__doc__
    __rxor__.__doc__ = operator.xor.__doc__
    __ror__.__doc__ = operator.or_.__doc__
    __iadd__.__doc__ = operator.iadd.__doc__
    __isub__.__doc__ = operator.isub.__doc__
    __imul__.__doc__ = operator.imul.__doc__
    __ipow__.__doc__ = operator.ipow.__doc__
    __itruediv__.__doc__ = operator.itruediv.__doc__
    __ifloordiv__.__doc__ = operator.ifloordiv.__doc__
    __imod__.__doc__ = operator.imod.__doc__
    __iand__.__doc__ = operator.iand.__doc__
    __ixor__.__doc__ = operator.ixor.__doc__
    __ior__.__doc__ = operator.ior.__doc__
    __ilshift__.__doc__ = operator.ilshift.__doc__
    __irshift__.__doc__ = operator.irshift.__doc__
    __neg__.__doc__ = operator.neg.__doc__
    __pos__.__doc__ = operator.pos.__doc__
    __abs__.__doc__ = operator.abs.__doc__
    __invert__.__doc__ = operator.invert.__doc__
    round.__doc__ = ops.round_.__doc__
    argsort.__doc__ = ops.argsort.__doc__
    conj.__doc__ = ops.conj.__doc__
    conjugate.__doc__ = ops.conjugate.__doc__


class DatasetGroupByOpsMixin:
    __slots__ = ()

    def _binary_op(
        self, other: Dataset | DataArray, f: Callable, reflexive: bool = False
    ) -> Dataset:
        raise NotImplementedError

    def __add__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.add)

    def __sub__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.sub)

    def __mul__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.mul)

    def __pow__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.pow)

    def __truediv__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.truediv)

    def __floordiv__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.floordiv)

    def __mod__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.mod)

    def __and__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.and_)

    def __xor__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.xor)

    def __or__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.or_)

    def __lshift__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.lshift)

    def __rshift__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.rshift)

    def __lt__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.lt)

    def __le__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.le)

    def __gt__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.gt)

    def __ge__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.ge)

    def __eq__(self, other: Dataset | DataArray) -> Dataset:  # type:ignore[override]
        return self._binary_op(other, nputils.array_eq)

    def __ne__(self, other: Dataset | DataArray) -> Dataset:  # type:ignore[override]
        return self._binary_op(other, nputils.array_ne)

    # When __eq__ is defined but __hash__ is not, then an object is unhashable,
    # and it should be declared as follows:
    __hash__: None  # type:ignore[assignment]

    def __radd__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.add, reflexive=True)

    def __rsub__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.sub, reflexive=True)

    def __rmul__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.mul, reflexive=True)

    def __rpow__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.pow, reflexive=True)

    def __rtruediv__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.truediv, reflexive=True)

    def __rfloordiv__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.floordiv, reflexive=True)

    def __rmod__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.mod, reflexive=True)

    def __rand__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.and_, reflexive=True)

    def __rxor__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.xor, reflexive=True)

    def __ror__(self, other: Dataset | DataArray) -> Dataset:
        return self._binary_op(other, operator.or_, reflexive=True)

    __add__.__doc__ = operator.add.__doc__
    __sub__.__doc__ = operator.sub.__doc__
    __mul__.__doc__ = operator.mul.__doc__
    __pow__.__doc__ = operator.pow.__doc__
    __truediv__.__doc__ = operator.truediv.__doc__
    __floordiv__.__doc__ = operator.floordiv.__doc__
    __mod__.__doc__ = operator.mod.__doc__
    __and__.__doc__ = operator.and_.__doc__
    __xor__.__doc__ = operator.xor.__doc__
    __or__.__doc__ = operator.or_.__doc__
    __lshift__.__doc__ = operator.lshift.__doc__
    __rshift__.__doc__ = operator.rshift.__doc__
    __lt__.__doc__ = operator.lt.__doc__
    __le__.__doc__ = operator.le.__doc__
    __gt__.__doc__ = operator.gt.__doc__
    __ge__.__doc__ = operator.ge.__doc__
    __eq__.__doc__ = nputils.array_eq.__doc__
    __ne__.__doc__ = nputils.array_ne.__doc__
    __radd__.__doc__ = operator.add.__doc__
    __rsub__.__doc__ = operator.sub.__doc__
    __rmul__.__doc__ = operator.mul.__doc__
    __rpow__.__doc__ = operator.pow.__doc__
    __rtruediv__.__doc__ = operator.truediv.__doc__
    __rfloordiv__.__doc__ = operator.floordiv.__doc__
    __rmod__.__doc__ = operator.mod.__doc__
    __rand__.__doc__ = operator.and_.__doc__
    __rxor__.__doc__ = operator.xor.__doc__
    __ror__.__doc__ = operator.or_.__doc__


class DataArrayGroupByOpsMixin:
    __slots__ = ()

    def _binary_op(
        self, other: T_Xarray, f: Callable, reflexive: bool = False
    ) -> T_Xarray:
        raise NotImplementedError

    def __add__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.add)

    def __sub__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.sub)

    def __mul__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.mul)

    def __pow__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.pow)

    def __truediv__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.truediv)

    def __floordiv__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.floordiv)

    def __mod__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.mod)

    def __and__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.and_)

    def __xor__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.xor)

    def __or__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.or_)

    def __lshift__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.lshift)

    def __rshift__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.rshift)

    def __lt__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.lt)

    def __le__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.le)

    def __gt__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.gt)

    def __ge__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.ge)

    def __eq__(self, other: T_Xarray) -> T_Xarray:  # type:ignore[override]
        return self._binary_op(other, nputils.array_eq)

    def __ne__(self, other: T_Xarray) -> T_Xarray:  # type:ignore[override]
        return self._binary_op(other, nputils.array_ne)

    # When __eq__ is defined but __hash__ is not, then an object is unhashable,
    # and it should be declared as follows:
    __hash__: None  # type:ignore[assignment]

    def __radd__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.add, reflexive=True)

    def __rsub__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.sub, reflexive=True)

    def __rmul__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.mul, reflexive=True)

    def __rpow__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.pow, reflexive=True)

    def __rtruediv__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.truediv, reflexive=True)

    def __rfloordiv__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.floordiv, reflexive=True)

    def __rmod__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.mod, reflexive=True)

    def __rand__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.and_, reflexive=True)

    def __rxor__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.xor, reflexive=True)

    def __ror__(self, other: T_Xarray) -> T_Xarray:
        return self._binary_op(other, operator.or_, reflexive=True)

    __add__.__doc__ = operator.add.__doc__
    __sub__.__doc__ = operator.sub.__doc__
    __mul__.__doc__ = operator.mul.__doc__
    __pow__.__doc__ = operator.pow.__doc__
    __truediv__.__doc__ = operator.truediv.__doc__
    __floordiv__.__doc__ = operator.floordiv.__doc__
    __mod__.__doc__ = operator.mod.__doc__
    __and__.__doc__ = operator.and_.__doc__
    __xor__.__doc__ = operator.xor.__doc__
    __or__.__doc__ = operator.or_.__doc__
    __lshift__.__doc__ = operator.lshift.__doc__
    __rshift__.__doc__ = operator.rshift.__doc__
    __lt__.__doc__ = operator.lt.__doc__
    __le__.__doc__ = operator.le.__doc__
    __gt__.__doc__ = operator.gt.__doc__
    __ge__.__doc__ = operator.ge.__doc__
    __eq__.__doc__ = nputils.array_eq.__doc__
    __ne__.__doc__ = nputils.array_ne.__doc__
    __radd__.__doc__ = operator.add.__doc__
    __rsub__.__doc__ = operator.sub.__doc__
    __rmul__.__doc__ = operator.mul.__doc__
    __rpow__.__doc__ = operator.pow.__doc__
    __rtruediv__.__doc__ = operator.truediv.__doc__
    __rfloordiv__.__doc__ = operator.floordiv.__doc__
    __rmod__.__doc__ = operator.mod.__doc__
    __rand__.__doc__ = operator.and_.__doc__
    __rxor__.__doc__ = operator.xor.__doc__
    __ror__.__doc__ = operator.or_.__doc__
