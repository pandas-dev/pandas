from pandas.compat import cPickle as pkl, pickle_compat as pc, PY3


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
        # try with current pickle, if we have a Type Error then
        # try with the compat pickle to handle subclass changes
        # pass encoding only if its not None as py2 doesn't handle
        # the param
        try:
            with open(path, 'rb') as fh:
                return pc.load(fh, encoding=encoding, compat=False)
        except:
            with open(path, 'rb') as fh:
                return pc.load(fh, encoding=encoding, compat=True)

    try:
        return try_read(path)
    except:
        if PY3:
            return try_read(path, encoding='latin1')
        raise
