from __future__ import annotations

from typing import TYPE_CHECKING

from virtualenv.activation.via_template import ViaTemplateActivator

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from virtualenv.create.creator import Creator


class NushellActivator(ViaTemplateActivator):
    def templates(self) -> Iterator[str]:
        yield "activate.nu"

    @staticmethod
    def quote(string: str) -> str:
        """Nushell supports raw strings like: r###'this is a string'###.

        https://github.com/nushell/nushell.github.io/blob/main/book/working_with_strings.md

        This method finds the maximum continuous sharps in the string and then quote it with an extra sharp.

        """
        max_sharps = 0
        current_sharps = 0
        for char in string:
            if char == "#":
                current_sharps += 1
                max_sharps = max(current_sharps, max_sharps)
            else:
                current_sharps = 0
        wrapping = "#" * (max_sharps + 1)
        return f"r{wrapping}'{string}'{wrapping}"

    def replacements(self, creator: Creator, dest_folder: Path) -> dict[str, str]:  # ruff:ignore[unused-method-argument]
        return {
            "__VIRTUAL_PROMPT__": "" if self.flag_prompt is None else self.flag_prompt,
            "__VIRTUAL_ENV__": str(creator.dest),
            "__VIRTUAL_NAME__": creator.env_name,
            "__BIN_NAME__": str(creator.bin_dir.relative_to(creator.dest)),
            "__TCL_LIBRARY__": getattr(creator.interpreter, "tcl_lib", None) or "",
            "__TK_LIBRARY__": getattr(creator.interpreter, "tk_lib", None) or "",
        }


__all__ = [
    "NushellActivator",
]
