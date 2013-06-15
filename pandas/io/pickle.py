# XXX: HACK for NumPy 1.5.1 to suppress warnings
try:
    import cPickle as pickle
except ImportError:  # pragma: no cover
    import pickle

def to_pickle(obj, path):
    """
    Pickle (serialize) object to input file path

    Parameters
    ----------
    obj : any object
    path : string
        File path
    """
    f = open(path, 'wb')
    try:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    finally:
        f.close()

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
        with open(path,'rb') as fh:
            return pickle.load(fh)
    except:
        from pandas.util import py3compat
        if not py3compat.PY3:
            raise
        with open(path,'rb') as fh:
            return pickle.load(fh, encoding='latin1')