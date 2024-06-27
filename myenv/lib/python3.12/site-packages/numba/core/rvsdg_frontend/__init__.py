from numba.core.utils import PYVERSION

def _ensure_rvsdg_supported():
    """Check that rvsdg_frontend is usable
    """
    # Only support Python 3.11.
    if PYVERSION != (3, 11):
        raise ImportError("rvsdg_frontend is only supported on python 3.11")
    # Require that numba_rvsdg to be installed.
    try:
        import numba_rvsdg
    except ImportError:
        raise ImportError("Failed to import numba_rvsdg package, "
                          "which is required for the rvsdg_frontend.")


_ensure_rvsdg_supported()
