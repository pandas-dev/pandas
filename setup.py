#!/usr/bin/env python3

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) which have been permitted for use under the
BSD license. Parts are from lxml (https://github.com/lxml/lxml)
"""

import pathlib
from Cython.Build import cythonize


if __name__ == "__main__":
    # TODO: this creates files but never cleans them up. Also runs serially
    # should integrate into CMake process so it is included as part of make /
    # make clean
    for x in pathlib.Path(".").glob("**/*.pyx"):
        cythonize(str(x), language_level=3)
