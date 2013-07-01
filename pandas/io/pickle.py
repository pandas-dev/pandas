import cPickle as pkl


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

    Warning: Loading pickled data received from untrusted sources can be unsafe.
    See: http://docs.python.org/2.7/library/pickle.html

    Parameters
    ----------
    path : string
        File path

    Returns
    -------
    unpickled : type of object stored in file
    """
    try:
        with open(path, 'rb') as fh:
            return pkl.load(fh)
    except:
        from pandas.util.py3compat import PY3
        if PY3:
            with open(path, 'rb') as fh:
                return pkl.load(fh, encoding='latin1')
        raise
