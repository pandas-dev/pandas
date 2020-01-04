"""
Functions to generate methods and pin them to the appropriate classes.
"""
import operator

from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries, ABCSparseArray

from pandas.core.ops.roperator import (
    radd,
    rand_,
    rdivmod,
    rfloordiv,
    rmod,
    rmul,
    ror_,
    rpow,
    rsub,
    rtruediv,
    rxor,
)


def _get_method_wrappers(cls):
    """
    Find the appropriate operation-wrappers to use when defining flex/special
    arithmetic, boolean, and comparison operations with the given class.

    Parameters
    ----------
    cls : class

    Returns
    -------
    arith_flex : function or None
    comp_flex : function or None
    arith_special : function
    comp_special : function
    bool_special : function

    Notes
    -----
    None is only returned for SparseArray
    """
    # TODO: make these non-runtime imports once the relevant functions
    #  are no longer in __init__
    from pandas.core.ops import (
        _arith_method_FRAME,
        _arith_method_SERIES,
        _bool_method_SERIES,
        _comp_method_FRAME,
        _comp_method_SERIES,
        _flex_comp_method_FRAME,
        _flex_method_SERIES,
    )

    if issubclass(cls, ABCSeries):
        # Just Series
        arith_flex = _flex_method_SERIES
        comp_flex = _flex_method_SERIES
        arith_special = _arith_method_SERIES
        comp_special = _comp_method_SERIES
        bool_special = _bool_method_SERIES
    elif issubclass(cls, ABCDataFrame):
        arith_flex = _arith_method_FRAME
        comp_flex = _flex_comp_method_FRAME
        arith_special = _arith_method_FRAME
        comp_special = _comp_method_FRAME
        bool_special = _arith_method_FRAME
    return arith_flex, comp_flex, arith_special, comp_special, bool_special


def add_special_arithmetic_methods(cls):
    """
    Adds the full suite of special arithmetic methods (``__add__``,
    ``__sub__``, etc.) to the class.

    Parameters
    ----------
    cls : class
        special methods will be defined and pinned to this class
    """
    _, _, arith_method, comp_method, bool_method = _get_method_wrappers(cls)
    new_methods = _create_methods(
        cls, arith_method, comp_method, bool_method, special=True
    )
    # inplace operators (I feel like these should get passed an `inplace=True`
    # or just be removed

    def _wrap_inplace_method(method):
        """
        return an inplace wrapper for this method
        """

        def f(self, other):
            result = method(self, other)

            # this makes sure that we are aligned like the input
            # we are updating inplace so we want to ignore is_copy
            self._update_inplace(
                result.reindex_like(self, copy=False)._data, verify_is_copy=False
            )

            return self

        name = method.__name__.strip("__")
        f.__name__ = f"__i{name}__"
        return f

    new_methods.update(
        dict(
            __iadd__=_wrap_inplace_method(new_methods["__add__"]),
            __isub__=_wrap_inplace_method(new_methods["__sub__"]),
            __imul__=_wrap_inplace_method(new_methods["__mul__"]),
            __itruediv__=_wrap_inplace_method(new_methods["__truediv__"]),
            __ifloordiv__=_wrap_inplace_method(new_methods["__floordiv__"]),
            __imod__=_wrap_inplace_method(new_methods["__mod__"]),
            __ipow__=_wrap_inplace_method(new_methods["__pow__"]),
        )
    )

    new_methods.update(
        dict(
            __iand__=_wrap_inplace_method(new_methods["__and__"]),
            __ior__=_wrap_inplace_method(new_methods["__or__"]),
            __ixor__=_wrap_inplace_method(new_methods["__xor__"]),
        )
    )

    _add_methods(cls, new_methods=new_methods)


def add_flex_arithmetic_methods(cls):
    """
    Adds the full suite of flex arithmetic methods (``pow``, ``mul``, ``add``)
    to the class.

    Parameters
    ----------
    cls : class
        flex methods will be defined and pinned to this class
    """
    flex_arith_method, flex_comp_method, _, _, _ = _get_method_wrappers(cls)
    new_methods = _create_methods(
        cls, flex_arith_method, flex_comp_method, bool_method=None, special=False
    )
    new_methods.update(
        dict(
            multiply=new_methods["mul"],
            subtract=new_methods["sub"],
            divide=new_methods["div"],
        )
    )
    # opt out of bool flex methods for now
    assert not any(kname in new_methods for kname in ("ror_", "rxor", "rand_"))

    _add_methods(cls, new_methods=new_methods)


def _create_methods(cls, arith_method, comp_method, bool_method, special):
    # creates actual methods based upon arithmetic, comp and bool method
    # constructors.

    have_divmod = issubclass(cls, ABCSeries)
    # divmod is available for Series

    new_methods = dict(
        add=arith_method(cls, operator.add, special),
        radd=arith_method(cls, radd, special),
        sub=arith_method(cls, operator.sub, special),
        mul=arith_method(cls, operator.mul, special),
        truediv=arith_method(cls, operator.truediv, special),
        floordiv=arith_method(cls, operator.floordiv, special),
        # Causes a floating point exception in the tests when numexpr enabled,
        # so for now no speedup
        mod=arith_method(cls, operator.mod, special),
        pow=arith_method(cls, operator.pow, special),
        # not entirely sure why this is necessary, but previously was included
        # so it's here to maintain compatibility
        rmul=arith_method(cls, rmul, special),
        rsub=arith_method(cls, rsub, special),
        rtruediv=arith_method(cls, rtruediv, special),
        rfloordiv=arith_method(cls, rfloordiv, special),
        rpow=arith_method(cls, rpow, special),
        rmod=arith_method(cls, rmod, special),
    )
    new_methods["div"] = new_methods["truediv"]
    new_methods["rdiv"] = new_methods["rtruediv"]
    if have_divmod:
        # divmod doesn't have an op that is supported by numexpr
        new_methods["divmod"] = arith_method(cls, divmod, special)
        new_methods["rdivmod"] = arith_method(cls, rdivmod, special)

    new_methods.update(
        dict(
            eq=comp_method(cls, operator.eq, special),
            ne=comp_method(cls, operator.ne, special),
            lt=comp_method(cls, operator.lt, special),
            gt=comp_method(cls, operator.gt, special),
            le=comp_method(cls, operator.le, special),
            ge=comp_method(cls, operator.ge, special),
        )
    )

    if bool_method:
        new_methods.update(
            dict(
                and_=bool_method(cls, operator.and_, special),
                or_=bool_method(cls, operator.or_, special),
                # For some reason ``^`` wasn't used in original.
                xor=bool_method(cls, operator.xor, special),
                rand_=bool_method(cls, rand_, special),
                ror_=bool_method(cls, ror_, special),
                rxor=bool_method(cls, rxor, special),
            )
        )

    if special:
        dunderize = lambda x: f"__{x.strip('_')}__"
    else:
        dunderize = lambda x: x
    new_methods = {dunderize(k): v for k, v in new_methods.items()}
    return new_methods


def _add_methods(cls, new_methods):
    for name, method in new_methods.items():
        # For most methods, if we find that the class already has a method
        # of the same name, it is OK to over-write it.  The exception is
        # inplace methods (__iadd__, __isub__, ...) for SparseArray, which
        # retain the np.ndarray versions.
        force = not (issubclass(cls, ABCSparseArray) and name.startswith("__i"))
        if force or name not in cls.__dict__:
            setattr(cls, name, method)
