from pkg_resources import require, DistributionNotFound

try:
    __version__ = require('scikits.timeseries')[0].version
except DistributionNotFound:
    # package hasn't actually been installed. Importing directly from source
    # folder. Explicitly import setup.py to extract version number
    # This should only happen for developers of the package
    import imp
    import os
    _parent_dir = os.path.split(os.path.dirname(__file__))[0]
    _setup_dir = os.path.split(_parent_dir)[0]
    _setup_py = os.path.join(_setup_dir, "setup.py")

    if os.path.exists(_setup_py):
        _setup = imp.load_source("setup", os.path.join(_setup_dir, "setup.py"))
        __version__ = _setup.version
    else:
        # package not installed through setup tools, just leave version undefined
        __version__ = None

# Make an alias (sometimes it's easier to drop the _).
version = __version__