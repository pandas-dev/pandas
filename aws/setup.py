# @Author: richard
# @Date:   2018-12-04T17:54:43+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-04T20:16:54+00:00
import os
import io
from setuptools import setup, find_packages

namespace = 'vortexa_utils'
name = 'vortexa_utils_aws'
version = '1'
description = 'Vortexa AWS utils helper library',

dependencies = [
    'boto3',
    'pycryptodomex'
]

# Setup boilerplate below

package_root = os.path.abspath(os.path.dirname(__file__))

readme_filename = os.path.join(package_root, 'README.rst')
with io.open(readme_filename, encoding='utf-8') as readme_file:
    readme = readme_file.read()

packages = [
    package for package in find_packages()
    if package.startswith(namespace)
]

setup(
    name=name,
    version=version,
    description=description,
    long_description=readme,

    author='Richard Mathie',
    author_email='richard.mathie@vortexa.com',

    zip_safe=False,
    test_suite='nose2.collector.collector',
    tests_require=['nose2', 'pandas'],

    packages=packages,
    install_requires=dependencies,
    extras_require={
        'pandas': ['pandas']
    }
)
