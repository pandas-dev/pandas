import logging
import warnings

from importlib import metadata as importlib_metadata


_already_initialized = False
logger = logging.getLogger(__name__)


def init_all():
    """Execute all `numba_extensions` entry points with the name `init`

    If extensions have already been initialized, this function does nothing.
    """
    global _already_initialized
    if _already_initialized:
        return

    # Must put this here to avoid extensions re-triggering initialization
    _already_initialized = True

    def load_ep(entry_point):
        """Loads a given entry point. Warns and logs on failure.
        """
        logger.debug('Loading extension: %s', entry_point)
        try:
            func = entry_point.load()
            func()
        except Exception as e:
            msg = (f"Numba extension module '{entry_point.module}' "
                   f"failed to load due to '{type(e).__name__}({str(e)})'.")
            warnings.warn(msg, stacklevel=3)
            logger.debug('Extension loading failed for: %s', entry_point)

    eps = importlib_metadata.entry_points()
    # Split, Python 3.10+ and importlib_metadata 3.6+ have the "selectable"
    # interface, versions prior to that do not. See "compatibility note" in:
    # https://docs.python.org/3.10/library/importlib.metadata.html#entry-points
    if hasattr(eps, 'select'):
        for entry_point in eps.select(group="numba_extensions", name="init"):
            load_ep(entry_point)
    else:
        for entry_point in eps.get("numba_extensions", ()):
            if entry_point.name == "init":
                load_ep(entry_point)
