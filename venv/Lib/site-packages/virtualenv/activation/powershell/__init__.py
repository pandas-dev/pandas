from __future__ import annotations

from typing import TYPE_CHECKING

from virtualenv.activation.via_template import ViaTemplateActivator

if TYPE_CHECKING:
    from collections.abc import Iterator


class PowerShellActivator(ViaTemplateActivator):
    def templates(self) -> Iterator[str]:
        yield "activate.ps1"

    @staticmethod
    def quote(string: str) -> str:
        """This should satisfy PowerShell quoting rules [1], unless the quoted string is passed directly to Windows native commands [2].

        [1]: https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_quoting_rules
        [2]:
        https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_parsing#passing-arguments-that-contain-quote-characters

        """
        string = string.replace("'", "''")
        return f"'{string}'"


__all__ = [
    "PowerShellActivator",
]
