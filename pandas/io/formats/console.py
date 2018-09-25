"""
Internal module for console introspection
"""

import sys
import locale
from pandas.io.formats.terminal import get_terminal_size

# -----------------------------------------------------------------------------
# Global formatting options
_initial_defencoding = None


def detect_console_encoding():
    """
    Try to find the most capable encoding supported by the console.
    slightly modified from the way IPython handles the same issue.
    """
    global _initial_defencoding

    encoding = None
    try:
        encoding = sys.stdout.encoding or sys.stdin.encoding
    except (AttributeError, IOError):
        pass

    # try again for something better
    if not encoding or 'ascii' in encoding.lower():
        try:
            encoding = locale.getpreferredencoding()
        except Exception:
            pass

    # when all else fails. this will usually be "ascii"
    if not encoding or 'ascii' in encoding.lower():
        encoding = sys.getdefaultencoding()

    # GH3360, save the reported defencoding at import time
    # MPL backends may change it. Make available for debugging.
    if not _initial_defencoding:
        _initial_defencoding = sys.getdefaultencoding()

    return encoding


def get_console_size():
    """Return console size as tuple = (width, height).

    Returns (None,None) in non-interactive session.
    """
    from pandas import get_option

    display_width = get_option('display.width')
    # deprecated.
    display_height = get_option('display.max_rows')

    # Consider
    # interactive shell terminal, can detect term size
    # interactive non-shell terminal (ipnb/ipqtconsole), cannot detect term
    # size non-interactive script, should disregard term size

    # in addition
    # width,height have default values, but setting to 'None' signals
    # should use Auto-Detection, But only in interactive shell-terminal.
    # Simple. yeah.

    if in_interactive_session():
        if in_ipython_frontend():
            # sane defaults for interactive non-shell terminal
            # match default for width,height in config_init
            from pandas.core.config import get_default_val
            terminal_width = get_default_val('display.width')
            terminal_height = get_default_val('display.max_rows')
        else:
            # pure terminal
            terminal_width, terminal_height = get_terminal_size()
    else:
        terminal_width, terminal_height = None, None

    # Note if the User sets width/Height to None (auto-detection)
    # and we're in a script (non-inter), this will return (None,None)
    # caller needs to deal.
    return (display_width or terminal_width, display_height or terminal_height)


# ----------------------------------------------------------------------
# Detect our environment

def in_interactive_session():
    """ check if we're running in an interactive shell

    returns True if running under python/ipython interactive shell
    """
    from pandas import get_option

    def check_main():
        import __main__ as main
        return (not hasattr(main, '__file__') or
                get_option('mode.sim_interactive'))

    try:
        return __IPYTHON__ or check_main()  # noqa
    except:
        return check_main()


def in_qtconsole():
    """
    check if we're inside an IPython qtconsole

    .. deprecated:: 0.14.1
       This is no longer needed, or working, in IPython 3 and above.
    """
    try:
        ip = get_ipython()  # noqa
        front_end = (
            ip.config.get('KernelApp', {}).get('parent_appname', "") or
            ip.config.get('IPKernelApp', {}).get('parent_appname', ""))
        if 'qtconsole' in front_end.lower():
            return True
    except:
        return False
    return False


def in_ipnb():
    """
    check if we're inside an IPython Notebook

    .. deprecated:: 0.14.1
       This is no longer needed, or working, in IPython 3 and above.
    """
    try:
        ip = get_ipython()  # noqa
        front_end = (
            ip.config.get('KernelApp', {}).get('parent_appname', "") or
            ip.config.get('IPKernelApp', {}).get('parent_appname', ""))
        if 'notebook' in front_end.lower():
            return True
    except:
        return False
    return False


def in_ipython_frontend():
    """
    check if we're inside an an IPython zmq frontend
    """
    try:
        ip = get_ipython()  # noqa
        return 'zmq' in str(type(ip)).lower()
    except:
        pass

    return False
