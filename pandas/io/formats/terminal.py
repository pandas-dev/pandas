"""
Terminal utilities.
"""

__all__ = ['is_terminal']


def is_terminal():
    """
    Detect if Python is running in a terminal.

    Returns True if Python is running in a terminal or False if not.
    """
    try:
        ip = get_ipython()
    except NameError:  # assume standard Python interpreter in a terminal
        return True
    else:
        if hasattr(ip, 'kernel'):  # IPython as a Jupyter kernel
            return False
        else:  # IPython in a terminal
            return True
