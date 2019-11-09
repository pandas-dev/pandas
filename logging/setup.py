import io
import os

from setuptools import find_packages, setup

namespace = "vortexa_utils"
description = ("Vortexa Error Logging",)

# Setup boilerplate below

package_root = os.path.abspath(os.path.dirname(__file__))

readme_filename = os.path.join(package_root, "README.md")
with io.open(readme_filename, encoding="utf-8") as readme_file:
    readme = readme_file.read()

packages = [
    package for package in find_packages() if package.startswith(namespace)
]

requirements = [
    "logzero",
    "psutil"
]

setup(
    name="vortexa_utils_logging",
    version="0.0.1",
    description=description,
    long_description=readme,
    author="Tino von Stegmann",
    author_email="constantin.vonstegmann@vortexa.com",
    zip_safe=False,
    tests_require=["nose2"],
    install_requires=requirements,
    test_suite="nose2.collector.collector",
    packages=packages,
)
