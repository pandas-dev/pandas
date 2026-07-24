from typing import Final

from .core import Shell as Shell, make_click_shell as make_click_shell
from .decorators import shell as shell

__all__ = ["make_click_shell", "shell", "Shell", "__version__"]
__version__: Final[str]
