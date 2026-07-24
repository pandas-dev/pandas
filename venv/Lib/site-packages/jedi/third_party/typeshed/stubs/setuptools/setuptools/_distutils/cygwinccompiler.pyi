from collections.abc import Callable

from .compilers.C import cygwin
from .compilers.C.cygwin import (
    CONFIG_H_NOTOK as CONFIG_H_NOTOK,
    CONFIG_H_OK as CONFIG_H_OK,
    CONFIG_H_UNCERTAIN as CONFIG_H_UNCERTAIN,
    check_config_h as check_config_h,
    get_msvcr as get_msvcr,
    is_cygwincc as is_cygwincc,
)
from .version import LooseVersion

__all__ = [
    "CONFIG_H_NOTOK",
    "CONFIG_H_OK",
    "CONFIG_H_UNCERTAIN",
    "CygwinCCompiler",
    "Mingw32CCompiler",
    "check_config_h",
    "get_msvcr",
    "is_cygwincc",
]

CygwinCCompiler = cygwin.Compiler
Mingw32CCompiler = cygwin.MinGW32Compiler

get_versions: Callable[[], tuple[LooseVersion | None, ...]] | None
