# -*- coding: utf-8 -*-

"""General constants used to spawn a PTY."""


class Backend:
    """Available PTY backends."""
    ConPTY = 0
    WinPTY = 1


class Encoding:
    """Available byte encodings to communicate with a PTY."""
    UTF8 = 'utf-8'
    UTF16 = 'utf-16'


class MouseMode:
    """Mouse capture settings for the winpty backend."""

    # QuickEdit mode is initially disabled, and the agent does not send mouse
    # mode sequences to the terminal.  If it receives mouse input, though, it
    # still writes MOUSE_EVENT_RECORD values into CONIN.
    WINPTY_MOUSE_MODE_NONE = 0

    # QuickEdit mode is initially enabled.  As CONIN enters or leaves mouse
    # input mode (i.e. where ENABLE_MOUSE_INPUT is on and
    # ENABLE_QUICK_EDIT_MODE is off), the agent enables or disables mouse
    # input on the terminal.
    WINPTY_MOUSE_MODE_AUTO = 1

    # QuickEdit mode is initially disabled, and the agent enables the
    # terminal's mouse input mode.  It does not disable terminal
    # mouse mode (until exit).
    WINPTY_MOUSE_MODE_FORCE = 2


class AgentConfig:
    """General configuration settings for the winpty backend."""

    # Create a new screen buffer (connected to the "conerr" terminal pipe) and
    # pass it to child processes as the STDERR handle.  This flag also prevents
    # the agent from reopening CONOUT$ when it polls -- regardless of whether
    # the active screen buffer changes, winpty continues to monitor the
    # original primary screen buffer.
    WINPTY_FLAG_CONERR = 0x1

    # Don't output escape sequences.
    WINPTY_FLAG_PLAIN_OUTPUT = 0x2

    # Do output color escape sequences.  These escapes are output by default,
    # but are suppressed with WINPTY_FLAG_PLAIN_OUTPUT.
    # Use this flag to re-enable them.
    WINPTY_FLAG_COLOR_ESCAPES = 0x4
