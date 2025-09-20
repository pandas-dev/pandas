# Licensed under a 3-clause BSD style license - see LICENSE.rst

from importlib_metadata import version as get_version

from asv import plugin_manager  # noqa: F401 Needed to load the plugins

__version__ = get_version("asv")

__all__ = ('__version__',)
