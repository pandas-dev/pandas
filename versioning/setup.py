# @Author: richard
# @Date:   2018-12-04T17:54:43+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-04T19:14:15+00:00
import os
from setuptools import setup, find_packages
from vortexa_utils.versioning import __version__

namespace = 'vortexa_utils'

# Setup boilerplate below

package_root = os.path.abspath(os.path.dirname(__file__))

packages = [
    package for package in find_packages()
    if package.startswith(namespace)
]

setup(
    name="vortexa_utils_versioning",
    version=__version__,
    description="",
    long_description="",

    author="Richard Mathie",
    author_email="richard.mathie@vortexa.com",

    zip_safe=False,
    tests_require=['nose2'],
    test_suite='nose2.collector.collector',

    packages=packages,
)
