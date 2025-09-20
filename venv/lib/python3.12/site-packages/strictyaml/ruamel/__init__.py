# coding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

if False:  # MYPY
    from typing import Dict, Any  # NOQA

_package_data = dict(
    full_package_name="strictyaml.ruamel",
    version_info=(0, 16, 13),
    __version__="0.16.13",
    author="Anthon van der Neut",
    author_email="a.van.der.neut@ruamel.eu",
    description="strictyaml.ruamel is a YAML parser/emitter that supports roundtrip preservation of comments, seq/map flow style, and map key order",  # NOQA
    entry_points=None,
    since=2014,
    extras_require={
        ':platform_python_implementation=="CPython" and python_version<="2.7"': [
            "ruamel.ordereddict"
        ],  # NOQA
        ':platform_python_implementation=="CPython" and python_version<"3.10"': [
            "strictyaml.ruamel.clib>=0.1.2"
        ],  # NOQA
        "jinja2": ["strictyaml.ruamel.jinja2>=0.2"],
        "docs": ["ryd"],
    },
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Programming Language :: Python :: Implementation :: Jython",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Text Processing :: Markup",
        "Typing :: Typed",
    ],
    keywords="yaml 1.2 parser round-trip preserve quotes order config",
    read_the_docs="yaml",
    supported=[(2, 7), (3, 5)],  # minimum
    tox=dict(
        env="*",  # remove 'pn', no longer test narrow Python 2.7 for unicode patterns and PyPy
        deps="ruamel.std.pathlib",
        fl8excl="_test/lib",
    ),
    universal=True,
    rtfd="yaml",
)  # type: Dict[Any, Any]


version_info = _package_data["version_info"]
__version__ = _package_data["__version__"]

try:
    from .cyaml import *  # NOQA

    __with_libyaml__ = True
except (ImportError, ValueError):  # for Jython
    __with_libyaml__ = False

from strictyaml.ruamel.main import *  # NOQA
