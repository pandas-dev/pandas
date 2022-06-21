#!/usr/bin/env python3

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) which have been permitted for use under the
BSD license. Parts are from lxml (https://github.com/lxml/lxml)
"""
from setuptools import setup

import versioneer


cmdclass = versioneer.get_cmdclass()


if __name__ == "__main__":
    setup(
        version=versioneer.get_version(),
        cmdclass=cmdclass,
    )
