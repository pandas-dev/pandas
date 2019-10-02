from contextlib import contextmanager
import os
import pytest
import numpy as np
from distutils.version import LooseVersion
import tempfile

from pandas.io.pytables import HDFStore


# TODO:
# remove when gh-24839 is fixed; this affects numpy 1.16
# and pytables 3.4.4
tables = pytest.importorskip("tables")
xfail_non_writeable = pytest.mark.xfail(
    LooseVersion(np.__version__) >= LooseVersion("1.16")
    and LooseVersion(tables.__version__) < LooseVersion("3.5.1"),
    reason=(
        "gh-25511, gh-24839. pytables needs a "
        "release beyong 3.4.4 to support numpy 1.16x"
    ),
)

# set these parameters so we don't have file sharing
tables.parameters.MAX_NUMEXPR_THREADS = 1
tables.parameters.MAX_BLOSC_THREADS = 1
tables.parameters.MAX_THREADS = 1


def safe_remove(path):
    if path is not None:
        try:
            os.remove(path)
        except OSError:
            pass


def safe_close(store):
    try:
        if store is not None:
            store.close()
    except IOError:
        pass


def create_tempfile(path):
    """ create an unopened named temporary file """
    return os.path.join(tempfile.gettempdir(), path)


# contextmanager to ensure the file cleanup
@contextmanager
def ensure_clean_store(path, mode="a", complevel=None, complib=None, fletcher32=False):

    try:

        # put in the temporary path if we don't have one already
        if not len(os.path.dirname(path)):
            path = create_tempfile(path)

        store = HDFStore(
            path, mode=mode, complevel=complevel, complib=complib, fletcher32=False
        )
        yield store
    finally:
        safe_close(store)
        if mode == "w" or mode == "a":
            safe_remove(path)


@contextmanager
def ensure_clean_path(path):
    """
    return essentially a named temporary file that is not opened
    and deleted on exiting; if path is a list, then create and
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


def _maybe_remove(store, key):
    """For tests using tables, try removing the table to be sure there is
    no content from previous tests using the same table name."""
    try:
        store.remove(key)
    except (ValueError, KeyError):
        pass
