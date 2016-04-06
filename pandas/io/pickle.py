""" pickle compat """

import numpy as np
from numpy.lib.format import read_array, write_array
from pandas.compat import BytesIO, cPickle as pkl, pickle_compat as pc, PY3
import pandas.core.common as com


def to_pickle(obj, path):
    """
    Pickle (serialize) object to input file path

    Parameters
    ----------
    obj : any object
    path : string
        File path
    """
    with open(path, 'wb') as f:
        pkl.dump(obj, f, protocol=pkl.HIGHEST_PROTOCOL)


def read_pickle(path):
    """
    Load pickled pandas object (or any other pickled object) from the specified
    file path

    Warning: Loading pickled data received from untrusted sources can be
    unsafe. See: http://docs.python.org/2.7/library/pickle.html

    Parameters
    ----------
    path : string
        File path

    Returns
    -------
    unpickled : type of object stored in file
    """

    def try_read(path, encoding=None):
        # try with cPickle
        # try with current pickle, if we have a Type Error then
        # try with the compat pickle to handle subclass changes
        # pass encoding only if its not None as py2 doesn't handle
        # the param

        # cpickle
        # GH 6899
        try:
            with open(path, 'rb') as fh:
                return pkl.load(fh)
        except Exception:
            # reg/patched pickle
            try:
                with open(path, 'rb') as fh:
                    return pc.load(fh, encoding=encoding, compat=False)

            # compat pickle
            except:
                with open(path, 'rb') as fh:
                    return pc.load(fh, encoding=encoding, compat=True)

    try:
        return try_read(path)
    except:
        if PY3:
            return try_read(path, encoding='latin1')
        raise

# compat with sparse pickle / unpickle


def _pickle_array(arr):
    arr = arr.view(np.ndarray)

    buf = BytesIO()
    write_array(buf, arr)

    return buf.getvalue()


def _unpickle_array(bytes):
    arr = read_array(BytesIO(bytes))

    # All datetimes should be stored as M8[ns].  When unpickling with
    # numpy1.6, it will read these as M8[us].  So this ensures all
    # datetime64 types are read as MS[ns]
    if com.is_datetime64_dtype(arr):
        arr = arr.view(com._NS_DTYPE)

    return arr
