from typing import Any

from .setopt import option_base

def shquote(arg): ...

class alias(option_base):
    description: str
    command_consumes_arguments: bool
    user_options: Any
    boolean_options: Any
    args: Any
    remove: Any
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...

def format_alias(name, aliases): ...
