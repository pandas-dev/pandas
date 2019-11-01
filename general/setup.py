import io
import os

from setuptools import setup, find_packages

namespace = 'vortexa_utils'
name = 'vortexa_utils_general'
version = '1.0.0'
description = 'Vortexa general utils helper library',

dependencies = [
    'gitpython',
    'logzero',
    'tenacity'
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
    author='Marcin Szymanski',
    author_email='marcin.szymanski@vortexa.com',
    zip_safe=False,
    packages=packages,
    install_requires=dependencies,
)
