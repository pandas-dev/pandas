from contextlib import contextmanager
import os
import tempfile

import pandas.util.testing as tm

from pandas.io.pytables import HDFStore


class Base(object):

    @classmethod
    def setup_class(cls):

        # Pytables 3.0.0 deprecates lots of things
        tm.reset_testing_mode()

    @classmethod
    def teardown_class(cls):

        # Pytables 3.0.0 deprecates lots of things
        tm.set_testing_mode()

    def setup_method(self, method):
        self.path = 'tmp.__%s__.h5' % tm.rands(10)

    def teardown_method(self, method):
        pass


def safe_close(store):
    try:
        if store is not None:
            store.close()
    except IOError:
        pass


@contextmanager
def ensure_clean_store(path, mode='a', complevel=None, complib=None,
                       fletcher32=False):

    store = None
    try:

        # put in the temporary path if we don't have one already
        if not len(os.path.dirname(path)):
            path = create_tempfile(path)

        store = HDFStore(path, mode=mode, complevel=complevel,
                         complib=complib, fletcher32=False)
        yield store
    finally:
        safe_close(store)
        if mode == 'w' or mode == 'a':
            safe_remove(path)


@contextmanager
def ensure_clean_path(path):
    """
    return essentially a named temporary file that is not opened
    and deleted on existing; if path is a list, then create and
    return list of filenames
    """
    filenames = []
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


def safe_remove(path):
    if path is not None:
        try:
            os.remove(path)
        except OSError:
            pass


def create_tempfile(path):
    """ create an unopened named temporary file """
    return os.path.join(tempfile.gettempdir(), path)


def maybe_remove(store, key):
    """For tests using tables, try removing the table to be sure there is
    no content from previous tests using the same table name."""
    try:
        store.remove(key)
    except (ValueError, KeyError):
        pass
