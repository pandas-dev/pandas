from . import (
    bdist as bdist,
    bdist_rpm as bdist_rpm,
    build as build,
    build_clib as build_clib,
    build_ext as build_ext,
    build_py as build_py,
    # build_scripts as build_scripts,
    # check as check,
    # clean as clean,
    install as install,
    install_data as install_data,
    # install_headers as install_headers,
    install_lib as install_lib,
    install_scripts as install_scripts,
    sdist as sdist,
)

# Commented out commands are not stubbed.
# (Many of these may be implementation details,
# but they can be added if people ask for them)
__all__ = [
    "build",
    "build_py",
    "build_ext",
    "build_clib",
    # "build_scripts",
    # "clean",
    "install",
    "install_lib",
    # "install_headers",
    "install_scripts",
    "install_data",
    "sdist",
    "bdist",
    # "bdist_dumb",
    "bdist_rpm",
    # "check",
]
