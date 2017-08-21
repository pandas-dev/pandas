""" pickle compat """

import numpy as np
from numpy.lib.format import read_array, write_array
from pandas.compat import BytesIO, cPickle as pkl, pickle_compat as pc, PY3
from pandas.core.dtypes.common import is_datetime64_dtype, _NS_DTYPE
from pandas.io.common import _get_handle, _infer_compression, _stringify_path


def to_pickle(obj, path, compression='infer', protocol=pkl.HIGHEST_PROTOCOL):
    """
    Pickle (serialize) object to input file path

    Parameters
    ----------
    obj : any object
    path : string
        File path
    compression : {'infer', 'gzip', 'bz2', 'xz', None}, default 'infer'
        a string representing the compression to use in the output file

        .. versionadded:: 0.20.0
    protocol : int
        Int which indicates which protocol should be used by the pickler,
        default HIGHEST_PROTOCOL (see [1], paragraph 12.1.2). The possible
        values for this parameter depend on the version of Python. For Python
        2.x, possible values are 0, 1, 2. For Python>=3.0, 3 is a valid value.
        For Python >= 3.4, 4 is a valid value. A negative value for the
        protocol parameter is equivalent to setting its value to
        HIGHEST_PROTOCOL.

        .. [1] https://docs.python.org/3/library/pickle.html
        .. versionadded:: 0.21.0


    """
    path = _stringify_path(path)
    inferred_compression = _infer_compression(path, compression)
    f, fh = _get_handle(path, 'wb',
                        compression=inferred_compression,
                        is_text=False)
    if protocol < 0:
        protocol = pkl.HIGHEST_PROTOCOL
    try:
        pkl.dump(obj, f, protocol=protocol)
    finally:
        for _f in fh:
            _f.close()


def read_pickle(path, compression='infer'):
    """
    Load pickled pandas object (or any other pickled object) from the specified
    file path

    Warning: Loading pickled data received from untrusted sources can be
    unsafe. See: http://docs.python.org/2.7/library/pickle.html

    Parameters
    ----------
    path : string
        File path
    compression : {'infer', 'gzip', 'bz2', 'xz', 'zip', None}, default 'infer'
        For on-the-fly decompression of on-disk data. If 'infer', then use
        gzip, bz2, xz or zip if path ends in '.gz', '.bz2', '.xz',
        or '.zip' respectively, and no decompression otherwise.
        Set to None for no decompression.

        .. versionadded:: 0.20.0

    Returns
    -------
    unpickled : type of object stored in file
    """
    path = _stringify_path(path)
    inferred_compression = _infer_compression(path, compression)

    def read_wrapper(func):
        # wrapper file handle open/close operation
        f, fh = _get_handle(path, 'rb',
                            compression=inferred_compression,
                            is_text=False)
        try:
            return func(f)
        finally:
            for _f in fh:
                _f.close()

    def try_read(path, encoding=None):
        # try with cPickle
        # try with current pickle, if we have a Type Error then
        # try with the compat pickle to handle subclass changes
        # pass encoding only if its not None as py2 doesn't handle
        # the param

        # cpickle
        # GH 6899
        try:
            return read_wrapper(lambda f: pkl.load(f))
        except Exception:
            # reg/patched pickle
            try:
                return read_wrapper(
                    lambda f: pc.load(f, encoding=encoding, compat=False))
            # compat pickle
            except:
                return read_wrapper(
                    lambda f: pc.load(f, encoding=encoding, compat=True))
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
    if is_datetime64_dtype(arr):
        arr = arr.view(_NS_DTYPE)

    return arr
