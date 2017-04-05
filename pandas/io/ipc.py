""" ipc format compat """

from pandas.types.generic import ABCIndexClass, ABCSeries, ABCDataFrame
from pandas.compat import BytesIO, string_types
from pandas._libs.lib import is_string_array, is_unicode_array
from pandas.types.common import is_object_dtype


def _try_import():
    # since pandas
    # we need to import on first use

    try:
        import pyarrow
    except ImportError:

        # give a nice error message
        raise ImportError("the pyarrow is not installed\n"
                          "you can install via conda\n"
                          "conda install pyarrow -c conda-forge")

    return pyarrow


def to_ipc(obj, engine='infer'):
    """
    Write a DataFrame to the ipc format

    Parameters
    ----------
    obj : Index, Series, DataFrame
    engine : string, optional
        string to indicate the engine {'infer', 'pickle', 'pyarrow'}
        'infer' will pick an engine based upon performance considerations

    Returns
    -------
    dict-of-metadata and bytes

    """
    if engine == 'pickle':
        return _to_pickle(obj)
    elif engine == 'pyarrow':
        try:
            return _to_pyarrow(obj)
        except:  # pragma
            pass

    if isinstance(obj, (ABCIndexClass, ABCSeries)):
        return _to_pickle(obj)
    elif isinstance(obj, ABCDataFrame):

        # decide quickly if we can serialize using
        # pyarrow or pickle

        # smallish, just pickle
        if len(obj) <= 100000:
            return _to_pickle(obj)

        # check our object columns
        for c, col in obj.iteritems():
            if not is_object_dtype(col):
                continue

            # if we discover we have actual python objects
            # embedded with strings/unicode, then pickle
            values = col.values
            if isinstance(values[0], string_types):
                if not is_string_array(values):
                    return _to_pickle(obj)
            else:
                if not is_unicode_array(values):
                    return _to_pickle(obj)

        return _to_pyarrow(obj)

    raise ValueError("ipc only supports IO with Index,"
                     "Series, DataFrames, a {} was "
                     "passed".format(type(obj)))


def _to_pyarrow(df):
    """ helper routine to return via pyarrow """
    pyarrow = _try_import()
    return pyarrow.write_ipc(df)


def _to_pickle(obj):
    """ helper routine to return a pickle of an object """
    from pandas import to_pickle
    db = BytesIO()
    to_pickle(obj, db)
    return db.getvalue()


def read_ipc(db, engine='infer'):
    """
    Load a pyarrow ipc format object from the file dict-of-bytes

    .. versionadded 0.20.0

    Parameters
    ----------
    dict-of-meta-and-bytes : a dictionary of meta data & bytes
    engine : string, optional
        string to indicate the engine {'infer', 'pickle', 'pyarrow'}
        'infer' will pick an engine based upon performance considerations

    Returns
    -------
    DataFrame

    """
    if engine == 'pickle':
        return _read_pickle(db)
    try:
        return _read_pyarrow(db)
    except:  # pragma
        return _read_pickle(db)


def _read_pyarrow(db):
    """ helper to return via pyarrow """
    pyarrow = _try_import()
    return pyarrow.read_ipc(db)


def _read_pickle(db):
    """ helper to return via pickle """
    from pandas import read_pickle

    db = BytesIO(db)
    return read_pickle(db)
