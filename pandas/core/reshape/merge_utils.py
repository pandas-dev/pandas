from __future__ import annotations

import numpy as np
import numpy.typing as npt

import pandas._libs.hashtable as libhashtable

import pandas.core.algorithms as algos

_factorizers = {
    np.int64: libhashtable.Int64Factorizer,
    np.longlong: libhashtable.Int64Factorizer,
    np.int32: libhashtable.Int32Factorizer,
    np.int16: libhashtable.Int16Factorizer,
    np.int8: libhashtable.Int8Factorizer,
    np.uint64: libhashtable.UInt64Factorizer,
    np.uint32: libhashtable.UInt32Factorizer,
    np.uint16: libhashtable.UInt16Factorizer,
    np.uint8: libhashtable.UInt8Factorizer,
    np.bool_: libhashtable.UInt8Factorizer,
    np.float64: libhashtable.Float64Factorizer,
    np.float32: libhashtable.Float32Factorizer,
    np.complex64: libhashtable.Complex64Factorizer,
    np.complex128: libhashtable.Complex128Factorizer,
    np.object_: libhashtable.ObjectFactorizer,
}

# See https://github.com/pandas-dev/pandas/issues/52451
if np.intc is not np.int32:
    _factorizers[np.intc] = libhashtable.Int64Factorizer


def _sort_labels(
    uniques: np.ndarray, left: npt.NDArray[np.intp], right: npt.NDArray[np.intp]
) -> tuple[npt.NDArray[np.intp], npt.NDArray[np.intp]]:
    llength = len(left)
    labels = np.concatenate([left, right])

    _, new_labels = algos.safe_sort(uniques, labels, use_na_sentinel=True)
    new_left, new_right = new_labels[:llength], new_labels[llength:]

    return new_left, new_right


def factorize_arrays(
    lk: np.ndarray,
    rk: np.ndarray,
    sort: bool = False,
    lk_mask: np.ndarray | None = None,
    rk_mask: np.ndarray | None = None,
):
    factorizer = _factorizers[lk.dtype.type](max(len(lk), len(rk)))
    llab = factorizer.factorize(lk, mask=lk_mask)
    rlab = factorizer.factorize(rk, mask=rk_mask)
    count = factorizer.get_count()

    if sort:
        uniques = factorizer.uniques.to_array()
        llab, rlab = _sort_labels(uniques, llab, rlab)

    lmask = llab == -1
    lany = lmask.any()
    rmask = rlab == -1
    rany = rmask.any()

    if lany or rany:
        if lany:
            np.putmask(llab, lmask, count)
        if rany:
            np.putmask(rlab, rmask, count)
        count += 1
    return llab, rlab, count
