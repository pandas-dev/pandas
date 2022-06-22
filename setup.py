#!/usr/bin/env python3

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) which have been permitted for use under the
BSD license. Parts are from lxml (https://github.com/lxml/lxml)
"""
import os
import sys

from setuptools import setup

# uncomment to enable pep517 after versioneer problem is fixed.
# https://github.com/python-versioneer/python-versioneer/issues/193
sys.path.insert(0, os.path.dirname(__file__))
import versioneer

cmdclass = versioneer.get_cmdclass()


if __name__ == "__main__":
    setup(
        version=versioneer.get_version(),
        cmdclass=cmdclass,
    )
