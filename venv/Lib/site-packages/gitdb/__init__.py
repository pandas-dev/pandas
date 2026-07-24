# Copyright (C) 2010, 2011 Sebastian Thiel (byronimo@gmail.com) and contributors
#
# This module is part of GitDB and is released under
# the New BSD License: https://opensource.org/license/bsd-3-clause/
"""Initialize the object database module"""

__author__ = "Sebastian Thiel"
__contact__ = "byronimo@gmail.com"
__homepage__ = "https://github.com/gitpython-developers/gitdb"
version_info = (4, 0, 12)
__version__ = '.'.join(str(i) for i in version_info)

# default imports
from gitdb.base import *
from gitdb.db import *
from gitdb.stream import *
