from pandas.compat import cPickle as pkl, pickle_compat as pc, PY3
from pandas.io.common import _get_handle, get_compression_type

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
    compression: {'gzip', 'bz2', 'infer', None}, default 'infer'
        Compression type,  ('infer' looks for the file extensions .gz and .bz2, using gzip and bz2 to decompress
        respectively).

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
        _compression = get_compression_type(path, compression)
        try:
            with _get_handle(path, 'rb', encoding, _compression) as fh:
                return pkl.load(fh)
        except (Exception) as e:

            # reg/patched pickle
            try:
                with _get_handle(path, 'rb', encoding, _compression) as fh:
                    return pc.load(fh, encoding=encoding, compat=False)

            # compat pickle
            except:
                with _get_handle(path, 'rb', encoding, _compression) as fh:
                    return pc.load(fh, encoding=encoding, compat=True)

    try:
        return try_read(path)
    except:
        if PY3:
            return try_read(path, encoding='latin1')
        raise
