"""
Internal module for console introspection
"""

from __future__ import annotations

from shutil import get_terminal_size


def get_console_size() -> tuple[int | None, int | None]:
    """
    Return console size as tuple = (width, height).

    In a non-interactive session an unset width is the terminal width (80 when
    not a tty) and an unset height is None.
    """
    from pandas._config.config import _global_config as config

    display_width = config["display"]["width"]
    display_height = config["display"]["max_rows"]

    # Consider
    # interactive shell terminal, can detect term size
    # interactive non-shell terminal (ipnb/ipqtconsole), cannot detect term
    # size non-interactive script, uses term width but not height

    # in addition
    # setting width/height to None signals auto-detection

    if in_interactive_session():
        if in_ipython_frontend():
            # sane defaults for interactive non-shell terminal, which cannot
            # auto-detect a width (GH#21337)
            from pandas._config.config import get_default_val

            terminal_width = 80
            terminal_height = get_default_val("display.max_rows")
        else:
            # pure terminal
            terminal_width, terminal_height = get_terminal_size()
    else:
        # fit script output to the terminal width, GH#21337
        terminal_width, terminal_height = get_terminal_size()[0], None

    # Note if the User sets height to None (auto-detection)
    # and we're in a script (non-inter), height will be None
    # caller needs to deal.
    return display_width or terminal_width, display_height or terminal_height


# ----------------------------------------------------------------------
# Detect our environment


def in_interactive_session() -> bool:
    """
    Check if we're running in an interactive shell.

    Returns
    -------
    bool
        True if running under python/ipython interactive shell.
    """
    from pandas._config.config import _global_config as config

    def check_main() -> bool:
        try:
            import __main__ as main
        except ModuleNotFoundError:
            return config["mode"]["sim_interactive"]
        return not hasattr(main, "__file__") or config["mode"]["sim_interactive"]

    try:
        # error: Name '__IPYTHON__' is not defined
        return __IPYTHON__ or check_main()  # type: ignore[name-defined]
    except NameError:
        return check_main()


def in_ipython_frontend() -> bool:
    """
    Check if we're inside an IPython zmq frontend.

    Returns
    -------
    bool
    """
    try:
        # error: Name 'get_ipython' is not defined
        ip = get_ipython()  # type: ignore[name-defined]
        return "zmq" in str(type(ip)).lower()
    except NameError:
        pass

    return False
