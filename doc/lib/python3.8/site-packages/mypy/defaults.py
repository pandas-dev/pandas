import os

from typing_extensions import Final

PYTHON2_VERSION = (2, 7)  # type: Final
PYTHON3_VERSION = (3, 6)  # type: Final
PYTHON3_VERSION_MIN = (3, 4)  # type: Final
CACHE_DIR = '.mypy_cache'  # type: Final
CONFIG_FILE = ['mypy.ini', '.mypy.ini']  # type: Final
PYPROJECT_CONFIG_FILES = ['pyproject.toml', ]  # type: Final
SHARED_CONFIG_FILES = ['setup.cfg', ]  # type: Final
USER_CONFIG_FILES = ['~/.config/mypy/config', '~/.mypy.ini', ]  # type: Final
if os.environ.get('XDG_CONFIG_HOME'):
    USER_CONFIG_FILES.insert(0, os.path.join(os.environ['XDG_CONFIG_HOME'], 'mypy/config'))

CONFIG_FILES = (CONFIG_FILE + PYPROJECT_CONFIG_FILES + SHARED_CONFIG_FILES +
               USER_CONFIG_FILES)  # type: Final

# This must include all reporters defined in mypy.report. This is defined here
# to make reporter names available without importing mypy.report -- this speeds
# up startup.
REPORTER_NAMES = ['linecount',
                  'any-exprs',
                  'linecoverage',
                  'memory-xml',
                  'cobertura-xml',
                  'xml',
                  'xslt-html',
                  'xslt-txt',
                  'html',
                  'txt',
                  'lineprecision']  # type: Final

# Threshold after which we sometimes filter out most errors to avoid very
# verbose output
MANY_ERRORS_THRESHOLD = 200  # type: Final
