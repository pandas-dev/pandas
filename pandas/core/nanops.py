import numpy as np

from pandas.core.common import isnull, notnull

def nansum(values, axis=0, skipna=True, copy=True):
    if values.dtype == np.object_:
        the_sum = values.sum(axis)
    else:
        mask = notnull(values)

        if skipna and not issubclass(values.dtype.type, np.integer):
            if copy:
                values = values.copy()
            np.putmask(values, -mask, 0)

        the_sum = values.sum(axis)
        the_count = mask.sum(axis)

        ct_mask = the_count == 0
        if ct_mask.any():
            the_sum[ct_mask] = np.nan

    return the_sum

def nanmean(values, axis=0, skipna=True, copy=True):
    if values.dtype == np.object_:
        the_mean = values.sum(axis) / float(values.shape[axis])
    else:
        mask = notnull(values)

        if skipna and not issubclass(values.dtype.type, np.integer):
            if copy:
                values = values.copy()
            np.putmask(values, -mask, 0)

        the_sum = values.sum(axis)
        the_count = mask.sum(axis)
        the_mean = the_sum / the_count.astype('f8')

        ct_mask = the_count == 0
        if ct_mask.any():
            the_mean[ct_mask] = np.nan

    return the_mean

def nanvar(values, axis=0, skipna=True, copy=True):
    mask = isnull(values)
    count = (values.shape[axis] - mask.sum(axis)).astype(float)

    if skipna:
        if copy:
            values = values.copy()
        np.putmask(values, mask, 0)

    X = values.sum(axis)
    XX = (values ** 2).sum(axis)
    return (XX - X ** 2 / count) / (count - 1)

def nanskew(values, axis=0, skipna=True, copy=True):
    if not isinstance(values.dtype.type, np.floating):
        values = values.astype('f8')

    mask = isnull(values)
    count = (values.shape[axis] - mask.sum(axis)).astype(float)

    if skipna:
        if copy:
            values = values.copy()
        np.putmask(values, mask, 0)

    A = values.sum(axis) / count
    B = (values ** 2).sum(axis) / count - A ** 2
    C = (values ** 3).sum(axis) / count - A ** 3 - 3 * A * B

    # floating point error
    B = np.where(np.abs(B) < 1e-14, 0, B)
    C = np.where(np.abs(C) < 1e-14, 0, C)

    result = ((np.sqrt((count ** 2 - count)) * C) /
              ((count - 2) * np.sqrt(B) ** 3))

    result = np.where(B == 0, 0, result)

    return result

def nanmin(values, axis=0, skipna=True, copy=True):
    if skipna and not issubclass(values.dtype.type, np.integer):
        if copy:
            values = values.copy()
        np.putmask(values, isnull(values), np.inf)
    return values.min(axis)

def nanmax(values, axis=0, skipna=True, copy=True):
    if skipna and not issubclass(values.dtype.type, np.integer):
        if copy:
            values = values.copy()
        np.putmask(values, isnull(values), -np.inf)

    return values.max(axis)

