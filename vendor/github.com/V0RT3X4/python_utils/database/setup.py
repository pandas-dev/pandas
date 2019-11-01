# @Author: richard
# @Date:   2018-12-04T17:54:43+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-04T19:14:15+00:00
import io
import os

from setuptools import find_packages, setup

namespace = "vortexa_utils"
description = ("Vortexa Database Engine Factory",)

dependencies = ["boto3", "SqlAlchemy", "psycopg2-binary", "requests"]

# Setup boilerplate below

package_root = os.path.abspath(os.path.dirname(__file__))

readme_filename = os.path.join(package_root, "README.rst")
with io.open(readme_filename, encoding="utf-8") as readme_file:
    readme = readme_file.read()

packages = [
    package for package in find_packages() if package.startswith(namespace)
]

setup(
    name="vortexa_utils_database",
    version="0.0.1",
    description=description,
    long_description=readme,
    author="Richard Mathie",
    author_email="richard.mathie@vortexa.com",
    zip_safe=False,
    tests_require=["nose2"],
    test_suite="nose2.collector.collector",
    packages=packages,
    install_requires=dependencies,
    extras_require={"query_cache": ["pandas", "pyarrow"]},
)
