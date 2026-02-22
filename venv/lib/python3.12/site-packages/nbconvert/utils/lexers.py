"""Deprecated as of 5.0; import from IPython.lib.lexers instead."""

from warnings import warn

warn("nbconvert.utils.lexers is deprecated as of 5.0. Use IPython.lib.lexers", stacklevel=2)

from IPython.lib.lexers import *  # noqa: F403, E402
