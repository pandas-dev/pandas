"""
Expressions
-----------

Offer fast expression evaluation through numexpr

"""
import operator
from typing import (
    List,
    Optional,
    Set,
)
import warnings

import numpy as np

from pandas._config import get_option

from pandas._typing import FuncType

from pandas.core.dtypes.generic import ABCDataFrame

from pandas.core.computation.check import NUMEXPR_INSTALLED
from pandas.core.ops import roperator

if NUMEXPR_INSTALLED:
    import numexpr as ne

_TEST_MODE: Optional[bool] = None
_TEST_RESULT: List[bool] = []
USE_NUMEXPR = NUMEXPR_INSTALLED
_evaluate: Optional[FuncType] = None
_where: Optional[FuncType] = None

# the set of dtypes that we will allow pass to numexpr
_ALLOWED_DTYPES = {
    "evaluate": {"int64", "int32", "float64", "float32", "bool"},
    "where": {"int64", "float64", "bool"},
}

# the minimum prod shape that we will use numexpr
_MIN_ELEMENTS = 10000


def set_use_numexpr(v=True):
    # set/unset to use numexpr
    global USE_NUMEXPR
    if NUMEXPR_INSTALLED:
        USE_NUMEXPR = v

    # choose what we are going to do
    global _evaluate, _where

    _evaluate = _evaluate_numexpr if USE_NUMEXPR else _evaluate_standard
    _where = _where_numexpr if USE_NUMEXPR else _where_standard


def set_numexpr_threads(n=None):
    # if we are using numexpr, set the threads to n
    # otherwise reset
    if NUMEXPR_INSTALLED and USE_NUMEXPR:
        if n is None:
            n = ne.detect_number_of_cores()
        ne.set_num_threads(n)


def _evaluate_standard(op, op_str, a, b):
    """
    Standard evaluation.
    """
    if _TEST_MODE:
        _store_test_result(False)
    return op(a, b)


def _can_use_numexpr(op, op_str, a, b, dtype_check):
    """ return a boolean if we WILL be using numexpr """
    if op_str is not None:

        # required min elements (otherwise we are adding overhead)
        if a.size > _MIN_ELEMENTS:
            # check for dtype compatibility
            dtypes: Set[str] = set()
            for o in [a, b]:
                # Series implements dtypes, check for dimension count as well
                if hasattr(o, "dtypes") and o.ndim > 1:
                    s = o.dtypes.value_counts()
                    if len(s) > 1:
                        return False
                    dtypes |= set(s.index.astype(str))
                # ndarray and Series Case
                elif hasattr(o, "dtype"):
                    dtypes |= {o.dtype.name}

            # allowed are a superset
            if not len(dtypes) or _ALLOWED_DTYPES[dtype_check] >= dtypes:
                return True

    return False


def _evaluate_numexpr(op, op_str, a, b):
    result = None

    if _can_use_numexpr(op, op_str, a, b, "evaluate"):
        is_reversed = op.__name__.strip("_").startswith("r")
        if is_reversed:
            # we were originally called by a reversed op method
            a, b = b, a

        a_value = a
        b_value = b

        result = ne.evaluate(
            f"a_value {op_str} b_value",
            local_dict={"a_value": a_value, "b_value": b_value},
            casting="safe",
        )

    if _TEST_MODE:
        _store_test_result(result is not None)

    if result is None:
        result = _evaluate_standard(op, op_str, a, b)

    return result


_op_str_mapping = {
    operator.add: "+",
    roperator.radd: "+",
    operator.mul: "*",
    roperator.rmul: "*",
    operator.sub: "-",
    roperator.rsub: "-",
    operator.truediv: "/",
    roperator.rtruediv: "/",
    operator.floordiv: "//",
    roperator.rfloordiv: "//",
    # we require Python semantics for mod of negative for backwards compatibility
    # see https://github.com/pydata/numexpr/issues/365
    # so sticking with unaccelerated for now
    operator.mod: None,
    roperator.rmod: "%",
    operator.pow: "**",
    roperator.rpow: "**",
    operator.eq: "==",
    operator.ne: "!=",
    operator.le: "<=",
    operator.lt: "<",
    operator.ge: ">=",
    operator.gt: ">",
    operator.and_: "&",
    roperator.rand_: "&",
    operator.or_: "|",
    roperator.ror_: "|",
    operator.xor: "^",
    roperator.rxor: "^",
    divmod: None,
    roperator.rdivmod: None,
}


def _where_standard(cond, a, b):
    # Caller is responsible for extracting ndarray if necessary
    return np.where(cond, a, b)


def _where_numexpr(cond, a, b):
    # Caller is responsible for extracting ndarray if necessary
    result = None

    if _can_use_numexpr(None, "where", a, b, "where"):

        result = ne.evaluate(
            "where(cond_value, a_value, b_value)",
            local_dict={"cond_value": cond, "a_value": a, "b_value": b},
            casting="safe",
        )

    if result is None:
        result = _where_standard(cond, a, b)

    return result


# turn myself on
set_use_numexpr(get_option("compute.use_numexpr"))


def _has_bool_dtype(x):
    if isinstance(x, ABCDataFrame):
        return "bool" in x.dtypes
    try:
        return x.dtype == bool
    except AttributeError:
        return isinstance(x, (bool, np.bool_))


def _bool_arith_check(
    op_str, a, b, not_allowed=frozenset(("/", "//", "**")), unsupported=None
):
    if unsupported is None:
        unsupported = {"+": "|", "*": "&", "-": "^"}

    if _has_bool_dtype(a) and _has_bool_dtype(b):
        if op_str in unsupported:
            warnings.warn(
                f"evaluating in Python space because the {repr(op_str)} "
                "operator is not supported by numexpr for "
                f"the bool dtype, use {repr(unsupported[op_str])} instead"
            )
            return False

        if op_str in not_allowed:
            raise NotImplementedError(
                f"operator {repr(op_str)} not implemented for bool dtypes"
            )
    return True


def evaluate(op, a, b, use_numexpr: bool = True):
    """
    Evaluate and return the expression of the op on a and b.

    Parameters
    ----------
    op : the actual operand
    a : left operand
    b : right operand
    use_numexpr : bool, default True
        Whether to try to use numexpr.
    """
    op_str = _op_str_mapping[op]
    if op_str is not None:
        use_numexpr = use_numexpr and _bool_arith_check(op_str, a, b)
        if use_numexpr:
            # error: "None" not callable
            return _evaluate(op, op_str, a, b)  # type: ignore[misc]
    return _evaluate_standard(op, op_str, a, b)


def where(cond, a, b, use_numexpr=True):
    """
    Evaluate the where condition cond on a and b.

    Parameters
    ----------
    cond : np.ndarray[bool]
    a : return if cond is True
    b : return if cond is False
    use_numexpr : bool, default True
        Whether to try to use numexpr.
    """
    assert _where is not None
    return _where(cond, a, b) if use_numexpr else _where_standard(cond, a, b)


def set_test_mode(v: bool = True) -> None:
    """
    Keeps track of whether numexpr was used.

    Stores an additional ``True`` for every successful use of evaluate with
    numexpr since the last ``get_test_result``.
    """
    global _TEST_MODE, _TEST_RESULT
    _TEST_MODE = v
    _TEST_RESULT = []


def _store_test_result(used_numexpr: bool) -> None:
    global _TEST_RESULT
    if used_numexpr:
        _TEST_RESULT.append(used_numexpr)


def get_test_result() -> List[bool]:
    """
    Get test result and reset test_results.
    """
    global _TEST_RESULT
    res = _TEST_RESULT
    _TEST_RESULT = []
    return res
