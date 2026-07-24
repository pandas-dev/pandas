import numpy as np
from .._lib._bunch import _make_tuple_bunch

__all__ = ['_find_repeats']

# This is not a namedtuple for backwards compatibility. See PR #12983
TheilslopesResult = _make_tuple_bunch('TheilslopesResult',
                                      ['slope', 'intercept',
                                       'low_slope', 'high_slope'])
SiegelslopesResult = _make_tuple_bunch('SiegelslopesResult',
                                       ['slope', 'intercept'])


def _n_samples_optional_x(kwargs):
    return 2 if kwargs.get('x', None) is not None else 1


def _find_repeats(arr):
    # This function assumes it may clobber its input.
    if len(arr) == 0:
        return np.array(0, np.float64), np.array(0, np.intp)

    # XXX This cast was previously needed for the Fortran implementation,
    # should we ditch it?
    arr = np.asarray(arr, np.float64).ravel()
    arr.sort()

    # Taken from NumPy 1.9's np.unique.
    change = np.concatenate(([True], arr[1:] != arr[:-1]))
    unique = arr[change]
    change_idx = np.concatenate(np.nonzero(change) + ([arr.size],))
    freq = np.diff(change_idx)
    atleast2 = freq > 1
    return unique[atleast2], freq[atleast2]
