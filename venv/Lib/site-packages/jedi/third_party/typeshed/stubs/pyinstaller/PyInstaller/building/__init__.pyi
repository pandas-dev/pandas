from typing_extensions import TypeAlias

# PyiBlockCipher is deprecated and misleads users into thinking it adds any security. Runtime deprecation warning:
# DEPRECATION: Bytecode encryption will be removed in PyInstaller v6.
# Please remove cipher and block_cipher parameters from your spec file to avoid breakages on upgrade.
# For the rationale/alternatives see https://github.com/pyinstaller/pyinstaller/pull/6999
_PyiBlockCipher: TypeAlias = None  # noqa: Y047  # Used by other modules
