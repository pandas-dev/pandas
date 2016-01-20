"""
Functions for converting object to other types
"""

import numpy as np

import pandas as pd
from pandas.core.common import (_possibly_cast_to_datetime, is_object_dtype,
                                isnull)
import pandas.lib as lib


# TODO: Remove in 0.18 or 2017, which ever is sooner
def _possibly_convert_objects(values, convert_dates=True, convert_numeric=True,
                              convert_timedeltas=True, copy=True):
    """ if we have an object dtype, try to coerce dates and/or numbers """

    # if we have passed in a list or scalar
    if isinstance(values, (list, tuple)):
        values = np.array(values, dtype=np.object_)
    if not hasattr(values, 'dtype'):
        values = np.array([values], dtype=np.object_)

    # convert dates
    if convert_dates and values.dtype == np.object_:

        # we take an aggressive stance and convert to datetime64[ns]
        if convert_dates == 'coerce':
            new_values = _possibly_cast_to_datetime(values, 'M8[ns]',
                                                    errors='coerce')

            # if we are all nans then leave me alone
            if not isnull(new_values).all():
                values = new_values

        else:
            values = lib.maybe_convert_objects(values,
                                               convert_datetime=convert_dates)

    # convert timedeltas
    if convert_timedeltas and values.dtype == np.object_:

        if convert_timedeltas == 'coerce':
            from pandas.tseries.timedeltas import to_timedelta
            new_values = to_timedelta(values, coerce=True)

            # if we are all nans then leave me alone
            if not isnull(new_values).all():
                values = new_values

        else:
            values = lib.maybe_convert_objects(
                values, convert_timedelta=convert_timedeltas)

    # convert to numeric
    if values.dtype == np.object_:
        if convert_numeric:
            try:
                new_values = lib.maybe_convert_numeric(values, set(),
                                                       coerce_numeric=True)

                # if we are all nans then leave me alone
                if not isnull(new_values).all():
                    values = new_values

            except:
                pass
        else:
            # soft-conversion
            values = lib.maybe_convert_objects(values)

    values = values.copy() if copy else values

    return values


def _soft_convert_objects(values, datetime=True, numeric=True, timedelta=True,
                          coerce=False, copy=True):
    """ if we have an object dtype, try to coerce dates and/or numbers """

    conversion_count = sum((datetime, numeric, timedelta))
    if conversion_count == 0:
        raise ValueError('At least one of datetime, numeric or timedelta must '
                         'be True.')
    elif conversion_count > 1 and coerce:
        raise ValueError("Only one of 'datetime', 'numeric' or "
                         "'timedelta' can be True when when coerce=True.")

    if isinstance(values, (list, tuple)):
        # List or scalar
        values = np.array(values, dtype=np.object_)
    elif not hasattr(values, 'dtype'):
        values = np.array([values], dtype=np.object_)
    elif not is_object_dtype(values.dtype):
        # If not object, do not attempt conversion
        values = values.copy() if copy else values
        return values

    # If 1 flag is coerce, ensure 2 others are False
    if coerce:
        # Immediate return if coerce
        if datetime:
            return pd.to_datetime(values, errors='coerce', box=False)
        elif timedelta:
            return pd.to_timedelta(values, errors='coerce', box=False)
        elif numeric:
            return pd.to_numeric(values, errors='coerce')

    # Soft conversions
    if datetime:
        values = lib.maybe_convert_objects(values, convert_datetime=datetime)

    if timedelta and is_object_dtype(values.dtype):
        # Object check to ensure only run if previous did not convert
        values = lib.maybe_convert_objects(values, convert_timedelta=timedelta)

    if numeric and is_object_dtype(values.dtype):
        try:
            converted = lib.maybe_convert_numeric(values, set(),
                                                  coerce_numeric=True)
            # If all NaNs, then do not-alter
            values = converted if not isnull(converted).all() else values
            values = values.copy() if copy else values
        except:
            pass

    return values
