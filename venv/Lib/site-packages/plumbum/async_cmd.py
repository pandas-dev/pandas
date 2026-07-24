"""Async command module for convenient imports.

This module provides a convenient way to import async local commands directly::

    from plumbum.async_cmd import ls, grep, echo

    async def main():
        result = await ls("-la")
        print(result)

        # Pipeline support
        output = await (ls | grep["py"])()
        print(output)

This is the async equivalent of plumbum.cmd, using the same module-hack pattern
to allow direct imports of commands as async versions.
"""

from __future__ import annotations

import plumbum


def __getattr__(name: str) -> plumbum.commands.async_.AsyncLocalCommand:
    """Module-hack that allows ``from plumbum.async_cmd import some_program``

    This returns async local commands instead of sync commands.

    Args:
        name: The command name to look up

    Returns:
        AsyncLocalCommand instance

    Raises:
        AttributeError: If command is not found in PATH

    Example::

        from plumbum.async_cmd import ls, grep

        async def main():
            result = await ls("-la")
            output = await (ls | grep["test"])()
    """
    try:
        return plumbum.async_local[name]
    except plumbum.commands.CommandNotFound:
        raise AttributeError(name) from None
