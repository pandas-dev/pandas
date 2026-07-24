from __future__ import annotations

from rattler.rattler import PyPatchInstructions


class PatchInstructions:
    _patch_instructions: PyPatchInstructions

    @classmethod
    def _from_py_patch_instructions(cls, py_patch_instructions: PyPatchInstructions) -> PatchInstructions:
        """
        Construct Rattler PatchInstructions from FFI PyPatchInstructions object.
        """
        patch_instructions = cls.__new__(cls)
        patch_instructions._patch_instructions = py_patch_instructions
        return patch_instructions

    def __repr__(self) -> str:
        """
        Returns a representation of the PatchInstructions.
        """
        return f"{type(self).__name__}()"
