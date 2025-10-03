from __future__ import annotations

import plumbum


def __getattr__(name: str) -> plumbum.machines.LocalCommand:
    """The module-hack that allows us to use ``from plumbum.cmd import some_program``"""
    try:
        return plumbum.local[name]
    except plumbum.CommandNotFound:
        raise AttributeError(name) from None
