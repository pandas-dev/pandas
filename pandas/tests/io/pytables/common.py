import os
import tempfile
from contextlib import contextmanager
from distutils.version import LooseVersion

import pytest
import pandas.util.testing as tm
from pandas.io.pytables import HDFStore

tables = pytest.importorskip('tables')
_default_compressor = ('blosc' if LooseVersion(tables.__version__) >=
                       LooseVersion('2.2') else 'zlib')


def get_file_name():
    """ Create a file name for test data files """
    return 'tmp.__%s__.h5' % tm.rands(10)


def create_tempfile(file_name=None):
    """ Create an unopened named temporary file """
    if file_name is None:
        file_name = get_file_name()
    return os.path.join(tempfile.gettempdir(), file_name)


@contextmanager
def ensure_clean_path(path=None):
    """
    return essentially a named temporary file that is not opened
    and deleted on existing; if path is a list, then create and
    return list of filenames
    """
    try:
        if isinstance(path, list):
            filenames = [create_tempfile(p) for p in path]
            yield filenames
        else:
            filenames = [create_tempfile(path)]
            yield filenames[0]
    finally:
        for f in filenames:
            safe_remove(f)


@contextmanager
def ensure_clean_store(path=None, mode='a', complevel=None, complib=None,
                       fletcher32=False):

    try:
        if path is None:
            path = create_tempfile(path)

        store = HDFStore(path, mode=mode, complevel=complevel,
                         complib=complib, fletcher32=False)
        yield store
    finally:
        safe_close(store)
        if mode == 'w' or mode == 'a':
            safe_remove(path)


def safe_remove(path):
    if path is not None:
        try:
            os.remove(path)
        except:
            pass


def safe_close(store):
    try:
        if store is not None:
            store.close()
    except:
        pass


def _maybe_remove(store, key):
    """For tests using tables, try removing the table to be sure there is
    no content from previous tests using the same table name."""
    try:
        store.remove(key)
    except:
        pass


def _check_roundtrip(obj, comparator, compression=False, **kwargs):
    options = {}
    if compression:
        options['complib'] = _default_compressor

    with ensure_clean_store('w', **options) as store:
        store['obj'] = obj
        retrieved = store['obj']
        comparator(retrieved, obj, **kwargs)


def _check_double_roundtrip(obj, comparator, compression=False, **kwargs):
    options = {}
    if compression:
        options['complib'] = compression or _default_compressor

    with ensure_clean_store('w', **options) as store:
        store['obj'] = obj
        retrieved = store['obj']
        comparator(retrieved, obj, **kwargs)
        store['obj'] = retrieved
        again = store['obj']
        comparator(again, obj, **kwargs)


def _check_roundtrip_table(obj, comparator, compression=False):
    options = {}
    if compression:
        options['complib'] = _default_compressor

    with ensure_clean_store('w', **options) as store:
        store.put('obj', obj, format='table')
        retrieved = store['obj']

        comparator(retrieved, obj)
