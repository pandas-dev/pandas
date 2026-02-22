from . import backend_inline, config  # noqa

__version__ = "0.2.1"

# we can't ''.join(...) otherwise finding the version number at build time requires
# import which introduces IPython and matplotlib at build time, and thus circular
# dependencies.
version_info = tuple(int(s) for s in __version__.split(".")[:3])
