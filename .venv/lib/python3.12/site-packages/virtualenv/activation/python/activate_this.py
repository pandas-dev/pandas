"""Activate virtualenv for current interpreter:

import runpy runpy.run_path(this_file)

This can be used when you must use an existing Python interpreter, not the virtualenv bin/python.

"""  # noqa: D415

from __future__ import annotations

import os
import site
import sys

try:
    abs_file = os.path.abspath(__file__)
except NameError as exc:
    msg = "You must use import runpy; runpy.run_path(this_file)"
    raise AssertionError(msg) from exc

bin_dir = os.path.dirname(abs_file)
base = bin_dir[: -len(__BIN_NAME__) - 1]  # ty: ignore[unresolved-reference]

# prepend bin to PATH (this file is inside the bin directory)
os.environ["PATH"] = os.pathsep.join([bin_dir, *os.environ.get("PATH", "").split(os.pathsep)])
os.environ["VIRTUAL_ENV"] = base  # virtual env is right above bin directory
os.environ["VIRTUAL_ENV_PROMPT"] = __VIRTUAL_PROMPT__ or os.path.basename(base)  # ty: ignore[unresolved-reference]

# Set PKG_CONFIG_PATH to include the virtualenv's pkgconfig directory
pkg_config_path = os.path.join(base, "lib", "pkgconfig")
existing_pkg_config_path = os.environ.get("PKG_CONFIG_PATH", "")
if existing_pkg_config_path:
    os.environ["PKG_CONFIG_PATH"] = os.pathsep.join([pkg_config_path, existing_pkg_config_path])
else:
    os.environ["PKG_CONFIG_PATH"] = pkg_config_path

# add the virtual environments libraries to the host python import mechanism
prev_length = len(sys.path)
for lib in __LIB_FOLDERS__.split(os.pathsep):  # ty: ignore[unresolved-reference]
    path = os.path.realpath(os.path.join(bin_dir, lib))
    site.addsitedir(path.decode("utf-8") if __DECODE_PATH__ else path)  # ty: ignore[unresolved-reference,unresolved-attribute]
sys.path[:] = sys.path[prev_length:] + sys.path[0:prev_length]

sys.real_prefix = sys.prefix  # ty: ignore[unresolved-attribute]
sys.prefix = base
