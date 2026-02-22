import warnings
import numpy as np
from numpy import nanmean, nansum

__all__ = [
    "median",
    "nanmedian",
    "nansum",
    "nanmean",
    "nanvar",
    "nanstd",
    "nanmin",
    "nanmax",
    "nanargmin",
    "nanargmax",
    "ss",
    "anynan",
    "allnan",
]


def nanargmin(a, axis=None):
    "Slow nanargmin function used for unaccelerated dtypes."
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanargmin(a, axis=axis)


def nanargmax(a, axis=None):
    "Slow nanargmax function used for unaccelerated dtypes."
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanargmax(a, axis=axis)


def nanvar(a, axis=None, ddof=0):
    "Slow nanvar function used for unaccelerated dtypes."
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanvar(a, axis=axis, ddof=ddof)


def nanstd(a, axis=None, ddof=0):
    "Slow nanstd function used for unaccelerated dtypes."
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanstd(a, axis=axis, ddof=ddof)


def nanmin(a, axis=None):
    "Slow nanmin function used for unaccelerated dtypes."
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmin(a, axis=axis)


def nanmax(a, axis=None):
    "Slow nanmax function used for unaccelerated dtypes."
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmax(a, axis=axis)


def median(a, axis=None):
    "Slow median function used for unaccelerated dtypes."
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.median(a, axis=axis)


def nanmedian(a, axis=None):
    "Slow nanmedian function used for unaccelerated dtypes."
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.nanmedian(a, axis=axis)


def ss(a, axis=None):
    "Slow sum of squares used for unaccelerated dtypes."
    a = np.asarray(a)
    y = np.multiply(a, a).sum(axis)
    return y


def anynan(a, axis=None):
    "Slow check for Nans used for unaccelerated dtypes."
    return np.isnan(a).any(axis)


def allnan(a, axis=None):
    "Slow check for all Nans used for unaccelerated dtypes."
    return np.isnan(a).all(axis)
