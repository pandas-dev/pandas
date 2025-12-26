import pytest

tables = pytest.importorskip("tables")
# set these parameters so we don't have file sharing
tables.parameters.MAX_NUMEXPR_THREADS = 1
tables.parameters.MAX_BLOSC_THREADS = 1
tables.parameters.MAX_THREADS = 1


def safe_close(store):
    try:
        if store is not None:
            store.close()
    except OSError:
        pass


def _maybe_remove(store, key):
    """
    For tests using tables, try removing the table to be sure there is
    no content from previous tests using the same table name.
    """
    try:
        store.remove(key)
    except (ValueError, KeyError):
        pass
